#ifndef SIMPLE_EXAMPLE_SDITH_ARITHMETIC_P251P3_H_
#define SIMPLE_EXAMPLE_SDITH_ARITHMETIC_P251P3_H_

#include "p251.h"

#define SDITH_GEN_P251P3 0x4e2e54 // multiplicative generator of p251p3
#define SDITH_ORDER_P251P3 15813250U // order of p251p3

EXPORT_DECL uint32_t (*p251p3_mul_ct)(uint32_t x, uint32_t y);
EXPORT_DECL uint32_t (*p251p3_mul)(uint32_t x, uint32_t y);

EXPORT void p251p3_mul_init();

/** @brief multiplication in P251^3. */
EXPORT uint32_t p251p3_mul_naive(uint32_t x, uint32_t y);
/** @brief multiplication using LUT in P251^3. */
EXPORT uint32_t p251p3_mul_table(uint32_t x, uint32_t y);
/** @brief discrete log in P251^3. */
EXPORT uint32_t p251p3_dlog(uint32_t x);
/** @brief discrete exp in P251^3. */
EXPORT uint32_t p251p3_dexp(uint32_t x);
/** @brief log(x^p) where x is in log-scale */
EXPORT uint32_t p251p3_dlog_pow_l32(uint32_t logx, uint32_t p);
/** @brief log(xy) where x and y are in log-scale */
EXPORT uint32_t p251p3_dlog_mul_l32_l32(uint32_t logx,uint32_t logy);

/** @brief extended power in P251^3. */
EXPORT uint32_t p251p3_pow(uint32_t x, uint32_t p);
/** @brief addition in P251^3. */
EXPORT uint32_t p251p3_add(uint32_t x, uint32_t y);
/** @brief subtraction in P251^3. */
EXPORT uint32_t p251p3_sub(uint32_t x, uint32_t y);

EXPORT_DECL uint32_t const* sdith_p251p3_dexp_table;
EXPORT_DECL uint32_t const* sdith_p251p3_dlog_table;

void p251p3_create_log_tables();


#endif  // SIMPLE_EXAMPLE_SDITH_ARITHMETIC_P251P3_H_
