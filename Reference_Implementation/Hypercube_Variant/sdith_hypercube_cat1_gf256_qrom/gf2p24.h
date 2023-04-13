#ifndef SIMPLE_EXAMPLE_SDITH_ARITHMETIC_GF2P24_H_
#define SIMPLE_EXAMPLE_SDITH_ARITHMETIC_GF2P24_H_

#include "gf256.h"

#define SDITH_GEN_GF2P24 0x761852 // multiplicative generator of gf2p24

/**  @brief constant-time multiplication in gf256^3 (mod Y^3+GF256(2))  */
EXPORT_DECL uint32_t (*gf2p24_mul_ct)(uint32_t x, uint32_t y);
/**  @brief multiplication in gf256^3 (mod Y^3+GF256(2))  */
EXPORT_DECL uint32_t (*gf2p24_mul)(uint32_t x, uint32_t y);

// due to the huge space, gf2p24 log tables are allocated on the heap on demand
EXPORT_DECL uint32_t const* sdith_gf2p24_dexp_table;
EXPORT_DECL uint32_t const* sdith_gf2p24_dlog_table;

/** @brief initializes gf2p24 multiplications */
EXPORT void gf2p24_init();

EXPORT uint32_t gf2p24_mul_pclmul_ct(uint32_t x, uint32_t y);

/** @brief extended discrete log in GF2p24. log(0)=255 by convention, otherwise,
 * log in base r */
EXPORT uint32_t gf2p24_dlog(uint32_t x);
EXPORT uint32_t gf2p24_dexp(uint32_t x);

/** @brief log(x^p) where x is in log-scale in GF2p24. */
EXPORT uint32_t gf2p24_log_pow(uint32_t logx, uint32_t p);
EXPORT uint32_t gf2p24_log_mul_log_log(uint32_t logx, uint32_t logy);

/** @brief multiplication in GF2p24. (lookup table based) */
EXPORT uint32_t gf2p24_mul_table(uint32_t x, uint32_t y);

EXPORT void gf2p24_create_log_tables();


#endif  // SIMPLE_EXAMPLE_SDITH_ARITHMETIC_GF2P24_H_
