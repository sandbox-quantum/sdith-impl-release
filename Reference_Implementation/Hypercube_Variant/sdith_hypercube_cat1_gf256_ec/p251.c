#include "p251.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "assertions.h"

#define nullptr 0x0

void (*p251_vec_mat16cols_muladd)(void *vz, void const *vx, void const *my,
                                  uint64_t m) = nullptr;
void (*p251_vec_mat16cols_muladd_ct)(void *vz, void const *vx, void const *my,
                                     uint64_t m) = nullptr;
void (*p251_vec_mat128cols_muladd)(void *vz, void const *vx, void const *my,
                                   uint64_t m) = nullptr;
void (*p251_vec_mat128cols_muladd_ct)(void *vz, void const *vx, void const *my,
                                      uint64_t m) = nullptr;

/// Performs "z[] += x[] * y" bulk memory operation without reduction
EXPORT void p251_muladd_mem_lazy(void *vz, uint8_t y, const void *vx,
                                 int bytes) {
  uint32_t *const vzz = (uint32_t *)vz;
  const uint8_t *const vxx = (const uint8_t *)vx;
  for (int64_t i = 0; i < bytes; ++i) {
    vzz[i] = vzz[i] + (uint32_t)vxx[i] * (uint32_t)y;
  }
}

/// Performs "z[16] += vx[m] * my[m][16]"
/// parallelizaion over m
EXPORT void p251_vec_mat16cols_muladd_ref_ct(void *vz, void const *vx,
                                          void const *my, uint64_t m) {
  uint8_t const *const myy = (uint8_t const *)my;
  uint8_t const *const vxx = (uint8_t const *)vx;
  uint32_t scratch[16] = {0};
  for (uint64_t j = 0; j < m; ++j) {
    for (uint64_t i = 0; i < 16; ++i) {
      scratch[i] += (uint32_t)vxx[j] * (uint32_t)myy[16 * j + i];
    }
  }
  uint8_t *const vzz = (uint8_t *)vz;
  for (uint64_t i = 0; i < 16; i++) {
    vzz[i] = (scratch[i] + vzz[i]) % 251;
  }
}

/// Performs "z[n] += vx[m] * my[m][n]" where n is a multiple of 32
/// parallelization over n
EXPORT void p251_vec_mat128cols_muladd_ref_ct(void *vz, void const *vx,
                                           void const *my, uint64_t m) {
  const uint8_t (*const y)[128] = (uint8_t (*)[128]) my;
  const uint8_t *const x = (uint8_t *) vx;
  uint8_t *const z = (uint8_t*) vz;
  uint32_t vzz[128] = {0};
  for (uint64_t j = 0; j < 128; ++j) {
    vzz[j] = z[j];
  }
  for (uint64_t j = 0; j < m; ++j) {
    p251_muladd_mem_lazy(vzz, x[j], y[j], 128);
  }
  for (uint64_t j = 0; j < 128; ++j) {
    z[j] = vzz[j] % 251;
  }
}

EXPORT uint8_t p251_mul_naive(uint8_t x, uint8_t y) {
  uint16_t xx = x;
  uint16_t yy = y;
  return (xx * yy) % 251;
}

EXPORT uint8_t p251_mul_table(uint8_t x, uint8_t y) {
  ASSERT_DRAMATICALLY(sdith_p251_dlog_table, "dlog tables not loaded!");
  if (x == 0 || y == 0) return 0;
  uint16_t l = (uint16_t)(sdith_p251_dlog_table[x]) + (uint16_t)(sdith_p251_dlog_table[y]);
  if (l >= 250) l -= 250;
  return sdith_p251_dexp_table[l];
}

/** log(x^p) where x is in log-scale */
uint8_t p251_log_pow_log(uint8_t logx, uint8_t p) {
  if (logx == 250) return 250;
  return (p * logx) % 250;
}

/** x^p where x is in log-scale */
uint8_t p251_pow_log(uint8_t logx, uint8_t p) {
  ASSERT_DRAMATICALLY(sdith_p251_dlog_table, "dlog tables not loaded!");
  return sdith_p251_dexp_table[p251_log_pow_log(logx, p)];
}

EXPORT uint8_t p251_pow(uint8_t x, uint8_t p) {
  if (p == 0) return 1;
  return p251_pow_log(p251_dlog(x), p);
}

EXPORT uint8_t p251_dlog(uint8_t x) {
  ASSERT_DRAMATICALLY(sdith_p251_dlog_table, "dlog tables not loaded!");
  return sdith_p251_dlog_table[x];
}

EXPORT int p251_init_(__attribute__((unused)) int version) {
  p251_create_log_tables();
  p251_vec_mat16cols_muladd = p251_vec_mat16cols_muladd_ref_ct;
  p251_vec_mat16cols_muladd_ct = p251_vec_mat16cols_muladd_ref_ct;
  p251_vec_mat128cols_muladd = p251_vec_mat128cols_muladd_ref_ct;
  p251_vec_mat128cols_muladd_ct = p251_vec_mat128cols_muladd_ref_ct;
  return 0;
}

uint8_t const* sdith_p251_dexp_table = nullptr;
uint8_t const* sdith_p251_dlog_table = nullptr;

void p251_create_log_tables() {
  // don't re-create the tables if they are already initialized
  if (sdith_p251_dlog_table) return;
  uint8_t* dlog8_table = malloc(sizeof(uint8_t)*251);
  uint8_t* dexp8_table = malloc(sizeof(uint8_t)*251);
  // create the dlog table over p251
  dexp8_table[0] = 1;
  for (uint64_t i = 1; i < 250; ++i) {
    dexp8_table[i] = p251_mul_naive(dexp8_table[i-1], SDITH_GEN_P251);
  }
  dexp8_table[250] = 0;
  for (uint64_t i = 0; i < 251; ++i) {
    dlog8_table[dexp8_table[i]] = i;
  }
  sdith_p251_dexp_table = dexp8_table;
  sdith_p251_dlog_table = dlog8_table;
}
