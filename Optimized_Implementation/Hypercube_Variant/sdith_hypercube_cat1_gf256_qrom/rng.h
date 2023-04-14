#ifndef SIMPLE_EXAMPLE_RNG_H
#define SIMPLE_EXAMPLE_RNG_H

#ifdef __cplusplus
#include <cstdint>
#define EXPORT extern "C"
#else
#include <stdint.h>
#define EXPORT
#endif

/** @brief opaque structure that represents a RNG context */
typedef struct RNG_struct RNG_CTX;  // opaque structure
/** @brief creates a rng context out of a 32-bytes key/seed (256-bit) */
EXPORT RNG_CTX* sdith_rng_create_rng_ctx(void const* const seed, void const* const nonce);
/** @brief deletes a rng context instantiated with create_rng_ctx */
EXPORT void sdith_rng_free_rng_ctx(RNG_CTX* ctx);
/** @brief produces the next random bytes out of this context */
EXPORT void sdith_rng_next_bytes(RNG_CTX* ctx, void* out, int outLen);
/** @brief produces the next random bytes out of this context, but each byte is sampled within [0, 251) using rejection sampling. */
EXPORT void sdith_rng_next_bytes_mod251(RNG_CTX* ctx, void* out, int outLen);

/** @brief opaque structure that represents a HASH context */
typedef struct HASH_struct HASH_CTX;
HASH_CTX* sdith_hash_create_hash_ctx();
void sdith_hash_free_hash_ctx(HASH_CTX* ctx);
void sdith_hash_digest_update(HASH_CTX* ctx, void const* in, int inBytes);
void sdith_hash_finalize(HASH_CTX* ctx, void* dest, int destBytes);
void sdith_hash(void* dest, int destBytes, void const* data, int dataBytes);


/** @brief opaque structure that represents a TREE_PRG context */
typedef struct TREE_PRG_struct TREE_PRG_CTX;
/** @brief takes as input n/2 seeds and produces n seeds: out[i]=F_{first_tweak+i}(in[i/2]) */
EXPORT void sdith_tree_prg_seed_expand(TREE_PRG_CTX* key, void* out, const void *in, const uint64_t first_tweak, uint64_t n);
/** @brief takes as input n/2 seeds and produces n seeds: out[i:i+1]=PRNG(in[i//2]) */
EXPORT void sdith_tree_prg_seed_prng_expand(TREE_PRG_CTX* key, void* out, const void *in, uint64_t n);
/** @brief takes as input n seeds and produces n commitments: out[i]=F_{first_tweak+i}(in[i]) */
EXPORT void sdith_tree_prg_seed_commit(TREE_PRG_CTX* key, void* out, const void *in, const uint64_t first_tweak, uint64_t n);
/** @brief takes a tree_prg context */
EXPORT TREE_PRG_CTX* sdith_create_tree_prg_ctx(void const* const root_salt);
/** @brief free a tree_prg context */
EXPORT void sdith_free_tree_prg_ctx(TREE_PRG_CTX* key);

#ifdef __cplusplus
namespace rng {
    /** @brief opaque structure that represents a RNG context */
    using RNG_CTX = struct RNG_struct;  // opaque structure
    /** @brief creates a rng context out of a 32-bytes key/seed (256-bit) */
    inline RNG_CTX* create_rng_ctx(void const* const seed) { return sdith_rng_create_rng_ctx(seed); }
    /** @brief deletes a rng context instantiated with create_rng_ctx */
    inline void free_rng_ctx(RNG_CTX* ctx) { return sdith_rng_free_rng_ctx(ctx); }
    /** @brief produces the next random bytes out of this context */
    inline void next_bytes(RNG_CTX* ctx, void* out, int outLen) { return sdith_rng_next_bytes(ctx, out, outLen); }
    /** @brief produces the next random bytes out of this context, but each byte is sampled within [0, 251) using rejection sampling. */
    inline void next_bytes_mod251(RNG_CTX* ctx, void* out, int outLen) { return sdith_rng_next_bytes_mod251(ctx, out, outLen); }

    /** @brief opaque structure that represents a TREE_PRG context */
    using TREE_PRG_CTX = struct TREE_PRG_struct;
    /** @brief takes as input n/2 seeds and produces n seeds: out[i]=F_{first_tweak+i}(in[i/2]) */
    inline void tree_prg_seed_expand(TREE_PRG_CTX* key, void* out, const void *in, const uint64_t first_tweak, uint64_t n) {
      return sdith_tree_prg_seed_expand(key, out, in, first_tweak, n);
    }
    /** @brief takes as input n seeds and produces n commitments: out[i]=F_{first_tweak+i}(in[i]) */
    inline void tree_prg_seed_commit(TREE_PRG_CTX* key, void* out, const void *in, const uint64_t first_tweak, uint64_t n) {
      return sdith_tree_prg_seed_commit(key, out, in, first_tweak, n);
    }
    /** @brief takes a tree_prg context */
    inline TREE_PRG_CTX* create_tree_prg_ctx(void const* const root_salt) { return sdith_create_tree_prg_ctx(root_salt); }
    /** @brief free a tree_prg context */
    inline void free_tree_prg_ctx(TREE_PRG_CTX* key) { return sdith_free_tree_prg_ctx(key); }

}  // namespace rng
#endif // C++

#endif //SIMPLE_EXAMPLE_RNG_H
