#include "p251p3.h"

#include <stdlib.h>

#include "assertions.h"
#define nullptr 0x0;

uint32_t (*p251p3_mul_ct)(uint32_t x, uint32_t y) = nullptr;
uint32_t (*p251p3_mul)(uint32_t x, uint32_t y) = nullptr;

EXPORT void p251p3_mul_init() {
  p251p3_create_log_tables();
  p251p3_mul_ct = p251p3_mul_naive;
  p251p3_mul = p251p3_mul_table;
}

EXPORT uint32_t p251p3_mul_naive(uint32_t x, uint32_t y) {
  // modulus is x^3 - x - 3
  // hence, x^3 = x + 3
  uint32_t x0 = x & 0xFF;
  uint32_t x1 = (x >> 8) & 0xFF;
  uint32_t x2 = (x >> 16) & 0xFF;

  uint32_t y0 = y & 0xFF;
  uint32_t y1 = (y >> 8) & 0xFF;
  uint32_t y2 = (y >> 16) & 0xFF;

  uint32_t z0 = x0 * y0;
  uint32_t z1 = x0 * y1 + x1 * y0;
  uint32_t z2 = x0 * y2 + x2 * y0 + x1 * y1;
  uint32_t z3 = x1 * y2 + x2 * y1;
  uint32_t z4 = x2 * y2;

  z2 = z2 + z4;
  z1 = z1 + 3 * z4;
  z1 = z1 + z3;
  z0 = z0 + 3 * z3;

  return (z0 % 251) | ((z1 % 251) << 8) | ((z2 % 251) << 16);
}

uint32_t const* sdith_p251p3_dexp_table = nullptr;
uint32_t const* sdith_p251p3_dlog_table = nullptr;

EXPORT uint32_t p251p3_mul_table(uint32_t x, uint32_t y) {
  ASSERT_DRAMATICALLY(sdith_p251p3_dlog_table, "dlog tables not loaded!");
  if (x == 0 || y == 0) return 0;
  uint32_t l = sdith_p251p3_dlog_table[x] + sdith_p251p3_dlog_table[y];
  if (l >= SDITH_ORDER_P251P3) l -= SDITH_ORDER_P251P3;
  return sdith_p251p3_dexp_table[l];
}

/** log(xy) where x,y are in log-scale */
uint32_t p251p3_dlog_mul_l32_l32(uint32_t logx, uint32_t logy) {
  if (logx == SDITH_ORDER_P251P3 || logy == SDITH_ORDER_P251P3) return SDITH_ORDER_P251P3;
  uint32_t l = logx + logy;
  if (l >= SDITH_ORDER_P251P3) l -= SDITH_ORDER_P251P3;
  return l;
}


/** log(x^p) where x is in log-scale */
uint32_t p251p3_dlog_pow_l32(uint32_t logx, uint32_t p) {
  if (logx == SDITH_ORDER_P251P3) return SDITH_ORDER_P251P3;
  return ((uint64_t)p * (uint64_t)logx) % SDITH_ORDER_P251P3;
}

/** x^p where x is in log-scale */
uint32_t p251p3_pow_log(uint32_t logx, uint32_t p) {
  return sdith_p251p3_dexp_table[p251p3_dlog_pow_l32(logx, p)];
}

EXPORT uint32_t p251p3_pow(uint32_t x, uint32_t p) {
  if (p == 0) return 1;
  return p251p3_pow_log(p251p3_dlog(x), p);
}

/** @brief discrete log in P251^3. */
EXPORT uint32_t p251p3_dlog(uint32_t x) {
  ASSERT_DRAMATICALLY(sdith_p251p3_dlog_table, "dlog tables not loaded!");
  return sdith_p251p3_dlog_table[x];
}

EXPORT uint32_t p251p3_dexp(uint32_t x) {
  ASSERT_DRAMATICALLY(sdith_p251p3_dexp_table, "dlog tables not loaded!");
  return sdith_p251p3_dexp_table[x];
}

EXPORT uint32_t p251p3_add(uint32_t x, uint32_t y) {
  uint32_t x0 = x & 0xFF;
  uint32_t x1 = (x >> 8) & 0xFF;
  uint32_t x2 = (x >> 16) & 0xFF;

  uint32_t y0 = y & 0xFF;
  uint32_t y1 = (y >> 8) & 0xFF;
  uint32_t y2 = (y >> 16) & 0xFF;

  return ((x0 + y0) % 251) | ((x1 + y1) % 251) << 8 | ((x2 + y2) % 251) << 16;
}

EXPORT uint32_t p251p3_sub(uint32_t x, uint32_t y) {
  uint32_t x0 = x & 0xFF;
  uint32_t x1 = (x >> 8) & 0xFF;
  uint32_t x2 = (x >> 16) & 0xFF;

  uint32_t y0 = y & 0xFF;
  uint32_t y1 = (y >> 8) & 0xFF;
  uint32_t y2 = (y >> 16) & 0xFF;

  return ((251 + x0 - y0) % 251) | ((251 + x1 - y1) % 251) << 8 | ((251 + x2 - y2) % 251) << 16;
}

void p251p3_create_log_tables() {
  if (sdith_p251p3_dlog_table) return;
  p251_create_log_tables();
  // don't re-create the tables if they are already initialized
  uint32_t* dlog24_table = malloc(sizeof(uint32_t)*(1ul << 24));
  uint32_t* dexp24_table = malloc(sizeof(uint32_t)*(SDITH_ORDER_P251P3+1));
  // create the dlog table over p251^3
  static const uint64_t ppp_power = SDITH_ORDER_P251P3 / 250;
  for (uint64_t i = 0, j = 0; i < 250; i++, j += ppp_power) {
    dexp24_table[j] = sdith_p251_dexp_table[i];
  }
  uint32_t z24 = 1; // p251^3
  for (uint64_t i = 1; i < ppp_power; ++i) {
    z24 = p251p3_mul_naive(z24, SDITH_GEN_P251P3);
    dexp24_table[i] = z24;
    uint8_t le0 = sdith_p251_dlog_table[(uint8_t)(z24)];
    uint8_t le1 = sdith_p251_dlog_table[(uint8_t)(z24 >> 8)];
    uint8_t le2 = sdith_p251_dlog_table[(uint8_t)(z24 >> 16)];
    for (uint64_t j = 1; j < 250; ++j) {
      le0 = (le0 == 250) ? 250 : (le0 == 249) ? 0 : (le0 + 1);
      le1 = (le1 == 250) ? 250 : (le1 == 249) ? 0 : (le1 + 1);
      le2 = (le2 == 250) ? 250 : (le2 == 249) ? 0 : (le2 + 1);
      uint32_t v = (uint32_t)(sdith_p251_dexp_table[le0]) | ((uint32_t)(sdith_p251_dexp_table[le1]) << 8) | ((uint32_t)(sdith_p251_dexp_table[le2]) << 16);
      dexp24_table[i + j * ppp_power] = v;
    }
  }
  dexp24_table[SDITH_ORDER_P251P3] = 0;
  for (uint64_t i = 0; i <= SDITH_ORDER_P251P3; ++i) {
    dlog24_table[dexp24_table[i]] = i;
  }
  sdith_p251p3_dexp_table = dexp24_table;
  sdith_p251p3_dlog_table = dlog24_table;
}
