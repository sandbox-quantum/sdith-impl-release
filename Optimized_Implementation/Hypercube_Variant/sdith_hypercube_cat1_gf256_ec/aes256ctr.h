#ifndef AES256CTR_H
#define AES256CTR_H

#ifdef __cplusplus
#define EXPORT extern "C"
#define EXPORT_DECL extern "C"
#else
#define EXPORT
#define EXPORT_DECL extern
#endif

#include <stddef.h>
#include <stdint.h>

#define AES256CTR_BLOCKBYTES 64

/**
 * @brief initialize an aes-256-ctr PRNG context
 * @param key 256-bit key
 * @param nonce initial counter
 * @param state state to initialize (created on the stack, or with new/malloc)
 * @returns context of aes256ctr. must call aes256ctr_deinit() to destroy.
 */
EXPORT void *aes256ctr_init(const uint8_t key[32], uint8_t nonce[16]);

EXPORT void aes256ctr_deinit(void *ctx);

/**
 * @brief next PRNG bytes (encryptions of zero)
 * @param state the PRNG state initialized with aes256ctr_init
 * @param nblocks number of blocks of 64 bytes to generate
 */
EXPORT void aes256ctr_squeezeblocks(void *ctx, uint8_t *out, size_t nblocks);

typedef struct FIXEDKEYAES_CTX FIXEDKEYAES_CTX;

EXPORT void fixedkeyaes_prf_init();

/** @brief generates a fixed key aes prf context from a 256-bit key */
EXPORT FIXEDKEYAES_CTX *new_fixedkeyaes_prf(const void *key);

/**
 * @brief out = AESFIXED_KEY_PRF_tweak(in)
 * here: prf_tweak(in) = aes256_key(Y) XOR Y      where    Y=((GF256(2) * in)
 * XOR tweak)
 * @param key  key can be fixed, zero, or salted.
 * @param in  must be 16 bytes long.
 * @param out must be 16 bytes long.
 * @param tweak0 counter that must be unique per node. */
EXPORT void fixedkeyaes_prf(FIXEDKEYAES_CTX *key, void *out,
                                  void const *in, const uint64_t tweak0);

/** @brief deletes a fixed key aes prf context */
EXPORT void delete_fixedkeyaes_prf(FIXEDKEYAES_CTX *ctx);

#endif // AES256CTR_H
