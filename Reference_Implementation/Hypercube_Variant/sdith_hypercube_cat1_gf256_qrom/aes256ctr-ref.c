#include "aes256ctr.h"

#include <openssl/evp.h>
#include <openssl/ossl_typ.h>
#include <stdint.h>
#include <string.h>

void *aes256ctr_init(const uint8_t key[32], uint8_t nonce[16]) {
  EVP_CIPHER_CTX *ctx = EVP_CIPHER_CTX_new();
  EVP_CIPHER_CTX_init(ctx);
  EVP_EncryptInit_ex(ctx, EVP_aes_256_ctr(), 0, (unsigned char const *)key,
                     nonce);
  return ctx;
}

void aes256ctr_deinit(void *ctx) { EVP_CIPHER_CTX_free((EVP_CIPHER_CTX *)ctx); }

static inline void aes_encrypt4(EVP_CIPHER_CTX *ctx, unsigned char *out) {
  static const uint8_t ZERO16[16] = {0};
  int outl; // must be int since it is a return of ssl
  for (uint64_t i = 0; i < 4; ++i) {
    EVP_EncryptUpdate(ctx, out + i * 16, &outl, ZERO16, 16);
  }
}

void aes256ctr_squeezeblocks(void *ctx, uint8_t *out, size_t nblocks) {
  size_t i;
  for (i = 0; i < nblocks; i++) {
    aes_encrypt4((EVP_CIPHER_CTX *)ctx, out);
    out += 64;
  }
}

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

EXPORT void fixedkeyaes_prf_init() {}

// prf_tweak(in) = aes256_key(Y) XOR Y where Y=((GF256(2) * in) XOR tweak)
EXPORT void fixedkeyaes_prf(FIXEDKEYAES_CTX *ctx, void *out, void const *in,
                            const uint64_t tweak0) {
  EVP_CIPHER_CTX *evp_ctx = (EVP_CIPHER_CTX *)ctx;

  // Compute Y.
  uint8_t y[16] = {0};

  // Apply tweak.
  memcpy(y, &tweak0, sizeof(tweak0));
  memcpy(&y[8], &tweak0, sizeof(tweak0));

  // GF256(2) * in.
  uint8_t *in8 = (uint8_t *)in;
  for (uint64_t i = 0; i < 16; ++i) {
    y[i] ^= mul(in8[i], 2);
  }

  int len;
  EVP_EncryptUpdate(evp_ctx, out, &len, y, 16);

  uint8_t *out8 = (uint8_t *)out;
  for (uint64_t i = 0; i < 16; ++i) {
    out8[i] ^= y[i];
  }
}

EXPORT FIXEDKEYAES_CTX *new_fixedkeyaes_prf(const void *key) {
  EVP_CIPHER_CTX *evp_ctx = EVP_CIPHER_CTX_new();
  EVP_EncryptInit_ex(evp_ctx, EVP_aes_256_ecb(), NULL, key, NULL);
  return (FIXEDKEYAES_CTX *)evp_ctx;
}

EXPORT void delete_fixedkeyaes_prf(FIXEDKEYAES_CTX *ctx) {
  EVP_CIPHER_CTX_free((EVP_CIPHER_CTX *)ctx);
}
