/* Based heavily on public-domain code by Romain Dolbeau
 * Different handling of nonce+counter than original version using
 * separated 64-bit nonce and internal 64-bit counter, starting from zero
 * Public Domain */

#include <immintrin.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "aes256ctr.h"

typedef struct aes256ctr_struct {
  __m128i rkeys[16];
  __m128i n;
} aes256ctr_ctx;

static inline void aesni_encrypt4(uint8_t out[64], __m128i *n,
                                  const __m128i rkeys[16]) {
  __m128i f, f0, f1, f2, f3;
  const __m128i idx =
      _mm_set_epi8(8, 9, 10, 11, 12, 13, 14, 15, 7, 6, 5, 4, 3, 2, 1, 0);

  /* Load current counter value */
  f = _mm_load_si128(n);

  /* Increase counter in 4 consecutive blocks */
  f0 = _mm_shuffle_epi8(_mm_add_epi64(f, _mm_set_epi64x(0, 0)), idx);
  f1 = _mm_shuffle_epi8(_mm_add_epi64(f, _mm_set_epi64x(1, 0)), idx);
  f2 = _mm_shuffle_epi8(_mm_add_epi64(f, _mm_set_epi64x(2, 0)), idx);
  f3 = _mm_shuffle_epi8(_mm_add_epi64(f, _mm_set_epi64x(3, 0)), idx);

  /* Write counter for next iteration, increased by 4 */
  _mm_store_si128(n, _mm_add_epi64(f, _mm_set_epi64x(4, 0)));

  /* Actual AES encryption, 4x interleaved */
  f = _mm_load_si128(&rkeys[0]);
  f0 = _mm_xor_si128(f0, f);
  f1 = _mm_xor_si128(f1, f);
  f2 = _mm_xor_si128(f2, f);
  f3 = _mm_xor_si128(f3, f);

  for (int i = 1; i < 14; i++) {
    f = _mm_load_si128(&rkeys[i]);
    f0 = _mm_aesenc_si128(f0, f);
    f1 = _mm_aesenc_si128(f1, f);
    f2 = _mm_aesenc_si128(f2, f);
    f3 = _mm_aesenc_si128(f3, f);
  }

  f = _mm_load_si128(&rkeys[14]);
  f0 = _mm_aesenclast_si128(f0, f);
  f1 = _mm_aesenclast_si128(f1, f);
  f2 = _mm_aesenclast_si128(f2, f);
  f3 = _mm_aesenclast_si128(f3, f);

  /* Write results */
  _mm_storeu_si128((__m128i *)(out + 0), f0);
  _mm_storeu_si128((__m128i *)(out + 16), f1);
  _mm_storeu_si128((__m128i *)(out + 32), f2);
  _mm_storeu_si128((__m128i *)(out + 48), f3);
}

void *aes256ctr_init(const uint8_t key[32], uint8_t nonce[16]) {
  aes256ctr_ctx *state = (aes256ctr_ctx *)malloc(sizeof(aes256ctr_ctx));

  __m128i key0, key1, temp0, temp1, temp2, temp4;
  int idx = 0;

  key0 = _mm_loadu_si128((__m128i *)(key + 0));
  key1 = _mm_loadu_si128((__m128i *)(key + 16));
  state->n = _mm_loadl_epi64((__m128i *)nonce);

  state->rkeys[idx++] = key0;
  temp0 = key0;
  temp2 = key1;
  temp4 = _mm_setzero_si128();

#define BLOCK1(IMM)                                                            \
  temp1 = _mm_aeskeygenassist_si128(temp2, IMM);                               \
  state->rkeys[idx++] = temp2;                                                 \
  temp4 = (__m128i)_mm_shuffle_ps((__m128)temp4, (__m128)temp0, 0x10);         \
  temp0 = _mm_xor_si128(temp0, temp4);                                         \
  temp4 = (__m128i)_mm_shuffle_ps((__m128)temp4, (__m128)temp0, 0x8c);         \
  temp0 = _mm_xor_si128(temp0, temp4);                                         \
  temp1 = (__m128i)_mm_shuffle_ps((__m128)temp1, (__m128)temp1, 0xff);         \
  temp0 = _mm_xor_si128(temp0, temp1)

#define BLOCK2(IMM)                                                            \
  temp1 = _mm_aeskeygenassist_si128(temp0, IMM);                               \
  state->rkeys[idx++] = temp0;                                                 \
  temp4 = (__m128i)_mm_shuffle_ps((__m128)temp4, (__m128)temp2, 0x10);         \
  temp2 = _mm_xor_si128(temp2, temp4);                                         \
  temp4 = (__m128i)_mm_shuffle_ps((__m128)temp4, (__m128)temp2, 0x8c);         \
  temp2 = _mm_xor_si128(temp2, temp4);                                         \
  temp1 = (__m128i)_mm_shuffle_ps((__m128)temp1, (__m128)temp1, 0xaa);         \
  temp2 = _mm_xor_si128(temp2, temp1)

  BLOCK1(0x01);
  BLOCK2(0x01);

  BLOCK1(0x02);
  BLOCK2(0x02);

  BLOCK1(0x04);
  BLOCK2(0x04);

  BLOCK1(0x08);
  BLOCK2(0x08);

  BLOCK1(0x10);
  BLOCK2(0x10);

  BLOCK1(0x20);
  BLOCK2(0x20);

  BLOCK1(0x40);
  state->rkeys[idx++] = temp0;

  return state;
}

void aes256ctr_deinit(void *ctx) { free((aes256ctr_ctx *)ctx); }

void aes256ctr_squeezeblocks(void *ctx, uint8_t *out, size_t nblocks) {
  aes256ctr_ctx *state = (aes256ctr_ctx *)ctx;
  size_t i;
  for (i = 0; i < nblocks; i++) {
    aesni_encrypt4(out, &state->n, state->rkeys);
    out += 64;
  }
}

typedef struct {
  __m256i rkeys2[16];
  __m128i rkeys[16];
} aes256_ctx;

static __m128i TIMES_TWO[2];

/** @brief naive gf256 multiplication */
static uint8_t mul(uint8_t x, uint8_t y) {
  static const uint16_t B = 0x11b;
  uint16_t xx = x;
  uint16_t r = 0;
  for (uint8_t i = 0; i < 8; ++i)
    r ^= (-(uint16_t)((y >> i) & 1)) & (xx << i);
  for (uint8_t i = 15; i >= 8; --i)
    r ^= (-(uint16_t)((r >> i) & 1)) & (B << (i - 8));
  return r;
}

// here, we precompute the 2x16 bytes that correspond to the multiplication by
// B=0x1b we could actually hardcode this constant.
EXPORT void fixedkeyaes_prf_init() {
  uint8_t(*t)[16] = (uint8_t(*)[16])TIMES_TWO;
  for (uint64_t i = 0; i < 16; ++i) {
    t[0][i] = mul(2, i);
    t[1][i] = mul(2, i << 4);
  }
}

EXPORT void fixedkeyaes_prf(FIXEDKEYAES_CTX *ctx, void *out,
                                  void const *in, const uint64_t tweak0) {
  aes256_ctx *key = (aes256_ctx *)ctx;
  const __m128i TIMES_TWO_LO = _mm_load_si128(TIMES_TWO + 0);
  const __m128i TIMES_TWO_HI = _mm_load_si128(TIMES_TWO + 1);
  const __m128i MASK_LO = _mm_set1_epi8(0xf);
  const __m128i tweak = _mm_set_epi64x(tweak0, 0);
  const __m128i min = _mm_loadu_si128((__m128i *)in);
  const __m128i min_lo = _mm_and_si128(min, MASK_LO);
  const __m128i min_hi = _mm_and_si128(_mm_srli_epi16(min, 4), MASK_LO);
  const __m128i m = _mm_xor_si128(
      tweak, _mm_xor_si128(_mm_shuffle_epi8(TIMES_TWO_LO, min_lo),
                           _mm_shuffle_epi8(TIMES_TWO_HI,
                                            min_hi))); // TODO tweaks + 2*min

  /* Actual AES encryption, 4x interleaved */
  __m128i k = _mm_load_si128(key->rkeys + 0);
  __m128i m0 = _mm_xor_si128(k, m);

  for (int i = 1; i < 14; i++) {
    k = _mm_load_si128(key->rkeys + i);
    m0 = _mm_aesenc_si128(m0, k);
  }

  k = _mm_load_si128(key->rkeys + 14);
  m0 = _mm_aesenclast_si128(m0, k);

  /* xor it with m */
  m0 = _mm_xor_si128(m0, m);

  /* Write results */
  _mm_storeu_si128((__m128i *)out, m0);
}

EXPORT FIXEDKEYAES_CTX *new_fixedkeyaes_prf(const void *key) {
  aes256_ctx *state = aligned_alloc(16, sizeof(aes256_ctx));
  uint8_t nonce[16] = {0};
  aes256ctr_ctx *c = aes256ctr_init(key, nonce);
  memcpy(state->rkeys, c->rkeys, 256);
  uint8_t(*a)[16] = (uint8_t(*)[16])c->rkeys;
  uint8_t(*b)[2][16] = (uint8_t(*)[2][16])state->rkeys2;
  for (uint64_t i = 0; i < 16; ++i) {
    memcpy(b[i][0], a[i], 16);
    memcpy(b[i][1], a[i], 16);
  }
  aes256ctr_deinit(c);
  return (FIXEDKEYAES_CTX *)state;
}

EXPORT void delete_fixedkeyaes_prf(FIXEDKEYAES_CTX *ctx) {
  free((aes256_ctx *)ctx);
}
