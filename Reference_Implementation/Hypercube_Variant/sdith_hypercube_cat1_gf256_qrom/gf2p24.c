#include "gf2p24.h"

#include <stdlib.h>
#include "assertions.h"
#include "gf256.h"
#define nullptr 0x0

uint32_t (*gf2p24_mul_ct)(uint32_t x, uint32_t y) = nullptr;
uint32_t (*gf2p24_mul)(uint32_t x, uint32_t y) = nullptr;

EXPORT void gf2p24_init() {
  static int gf2p24_is_initialized = 0;
  if (gf2p24_is_initialized) return;
  gf256_init(0);
  //TODO currently, only the naive version of gf2p24_mul is available
  gf2p24_mul = gf2p24_mul_table;
#ifdef AVX2
  if (__builtin_cpu_supports("avx2")) {
    gf2p24_mul_ct = gf2p24_mul_pclmul_ct;
  } else {
    gf2p24_mul_ct = mul_gf2p24_naive;
  }
#else
  gf2p24_mul_ct = mul_gf2p24_naive;
#endif
  gf2p24_create_log_tables();
  gf2p24_is_initialized = 1;
}

uint32_t const* sdith_gf2p24_dexp_table = nullptr;
uint32_t const* sdith_gf2p24_dlog_table = nullptr;

/** @brief multiplication in gf2p24 (pclmul constant time) */
EXPORT uint32_t gf2p24_mul_pclmul_ct(uint32_t x, uint32_t y);

uint32_t gf2p24_dlog(uint32_t x) {
  ASSERT_DRAMATICALLY(sdith_gf2p24_dlog_table, "dlog tables not loaded!");
  return sdith_gf2p24_dlog_table[x];
}

uint32_t gf2p24_dexp(uint32_t x) {
  ASSERT_DRAMATICALLY(sdith_gf2p24_dexp_table, "dexp tables not loaded!");
  return sdith_gf2p24_dexp_table[x];
}

/** log(x^p) where x is in log-scale */
uint32_t gf2p24_log_pow(uint32_t logx, uint32_t p) {
  if (logx == 0xffffff) return 0xffffff;
  // we need 64-bits in this product to avoid rare overflows
  return ((uint64_t)p * logx) % 0xffffffUL;
}

uint32_t gf2p24_log_mul_log_log(uint32_t logx, uint32_t logy) {
  if (logx == 0xffffff || logy == 0xffffff) return 0xffffff;
  static const uint32_t order = (1ul << 24) - 1;
  uint32_t l = logx + logy;
  if (l >= order) l -= order;
  return l;
}

uint32_t gf2p24_mul_table(uint32_t x, uint32_t y) {
  ASSERT_DRAMATICALLY(sdith_gf2p24_dlog_table, "dlog tables not loaded!");
  static const uint32_t order = (1ul << 24) - 1;
  if (x == 0 || y == 0) return 0;
  uint32_t l = sdith_gf2p24_dlog_table[x] + sdith_gf2p24_dlog_table[y];
  if (l >= order) l -= order;
  return sdith_gf2p24_dexp_table[l];
}

void gf2p24_create_log_tables() {
  // don't re-create the tables if they are already initialized
  if (sdith_gf2p24_dlog_table) return;
  gf256_create_log_tables();
  uint32_t* dlog24_table = malloc(sizeof(uint32_t) * (1ul << 24));
  uint32_t* dexp24_table = malloc(sizeof(uint32_t) * (1ul << 24));
  // create the dlog table over gf2p24
  static const uint64_t order = (1UL << 24) - 1;
  uint32_t z24 = 1;  // GF2p24
  for (uint64_t i = 0, j = 0; i < 255; i++, j += 0x10101) {
    dexp24_table[j] = sdith_gf256_dexp_table[i];
  }
  for (uint64_t i = 1; i < 0x10101; ++i) {
    z24 = gf2p24_mul_ct(z24, SDITH_GEN_GF2P24);
    dexp24_table[i] = z24;
    uint8_t le0 = sdith_gf256_dlog_table[(uint8_t)(z24)];
    uint8_t le1 = sdith_gf256_dlog_table[(uint8_t)(z24 >> 8)];
    uint8_t le2 = sdith_gf256_dlog_table[(uint8_t)(z24 >> 16)];
    for (uint64_t j = 1; j < 255; ++j) {
      le0 = (le0 == 255) ? 255 : (le0 == 254) ? 0 : (le0 + 1);
      le1 = (le1 == 255) ? 255 : (le1 == 254) ? 0 : (le1 + 1);
      le2 = (le2 == 255) ? 255 : (le2 == 254) ? 0 : (le2 + 1);
      uint32_t v = (uint32_t)(sdith_gf256_dexp_table[le0]) | ((uint32_t)(sdith_gf256_dexp_table[le1]) << 8) |
                   ((uint32_t)(sdith_gf256_dexp_table[le2]) << 16);
      dexp24_table[i + j * 0x10101] = v;
    }
  }
  dexp24_table[order] = 0;
  for (uint64_t i = 0; i < (1ul << 24); ++i) {
    dlog24_table[dexp24_table[i]] = i;
  }
#ifndef NDEBUG
  for (uint64_t i = 0; i < (1ul << 8); ++i) {
    uint32_t el = dexp24_table[i * 0x10101];
    REQUIRE_DRAMATICALLY(el <= 255, "dexp[%ld]=%d is not in GF2p8\n", i * 0x10101, el);
  }
  for (uint64_t i = 0; i < (1ul << 8); ++i) {
    uint32_t le = dlog24_table[i];
    REQUIRE_DRAMATICALLY(le % 0x10101 == 0, "dlog[%ld]=%d\n", i, le);
  }
#endif  // NDEBUG
  sdith_gf2p24_dexp_table = dexp24_table;
  sdith_gf2p24_dlog_table = dlog24_table;
}
