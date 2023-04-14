#include "rng.h"

#include "param.h"

#include "aes256ctr.h"
#include "stdlib.h"

RNG_CTX *sdith_rng_create_rng_ctx(void const *const seed, void const *const nonce) {
  return (RNG_CTX *)aes256ctr_init(seed, (uint8_t*)nonce);
}

void sdith_rng_free_rng_ctx(RNG_CTX *ctx) { aes256ctr_deinit(ctx); }

void sdith_rng_next_bytes(RNG_CTX *ctx, void *out, int outLen) {
  uint8_t *out_ptr = (uint8_t *)out;
  aes256ctr_squeezeblocks(ctx, out_ptr, outLen >> 6);
  uint16_t rest = outLen & 0x3F;
  if (rest > 0) {
    uint8_t buf[64];
    aes256ctr_squeezeblocks(ctx, buf, 1);
    out_ptr += (outLen - rest);
    for (uint8_t i = 0; i < rest; i++)
      out_ptr[i] = buf[i];
  }
}

void sdith_rng_next_bytes_mod251(RNG_CTX *ctx, void *out, int outLen) {
  // Roughly sample 1.03x of original length, which is greater than 256/251.
  int len = outLen + (outLen >> 5);
  uint8_t *buf = (uint8_t *)malloc(len);
  sdith_rng_next_bytes(ctx, buf, len);
  int bytes_remaining = len;
  uint8_t *buf_ptr = buf;
  uint8_t *out_buf = (uint8_t *)out;
  int counter = 0;
  while (counter < outLen) {
    if (*buf_ptr < 251) {
      out_buf[counter++] = *buf_ptr;
    }
    buf_ptr++;
    bytes_remaining--;
    if (bytes_remaining == 0) {
      sdith_rng_next_bytes(ctx, buf, len);
      bytes_remaining = len;
      buf_ptr = buf;
    }
  }
  free(buf);
}
