#include <stdlib.h>
#include <string.h>

#include "assertions.h"
#include "param.h"
#include "types.h"

#include "aes256ctr.h"
#include "rng.h"

#ifndef AES_COMMITMENT
typedef struct treeprg_sha3_context_struct {
  salt_t salt;
} treeprg_sha3_context_t;
#endif

/** @brief takes as input n/2 seeds and produces n seeds: out[i:i+1]=PRNG(in[i//2]) */
EXPORT void sdith_tree_prg_seed_prng_expand(TREE_PRG_CTX* key, void* out, const void *in, uint64_t n) {
  treeprg_sha3_context_t *context = (treeprg_sha3_context_t *)key;
  const seed_t *in_seeds = (seed_t *)in;
  seed_t(*out_seeds)[2] = (seed_t(*)[2])out;

  for (uint64_t i = 0; i < n/2; ++i) {
    RNG_CTX* leaf_seed_rng_ctx = sdith_rng_create_rng_ctx(in_seeds[i], context->salt);
    sdith_rng_next_bytes(leaf_seed_rng_ctx, out_seeds[i][0], PARAM_seed_size);
    sdith_rng_next_bytes(leaf_seed_rng_ctx, out_seeds[i][1], PARAM_seed_size);
    sdith_rng_free_rng_ctx(leaf_seed_rng_ctx);
  }
}

/** @brief takes as input n/2 seeds and produces n seeds:
 * out[i]=F_{first_tweak+i}(in[i/2]) */
EXPORT void sdith_tree_prg_seed_expand(TREE_PRG_CTX *key, void *out,
                                       const void *in,
                                       const uint64_t first_tweak, uint64_t n) {
  ASSERT_DRAMATICALLY(n % 2 == 0, "n must be even for this function");
#ifdef AES_COMMITMENT
  // TODO: use avx2
  uint8_t(*sin)[16] = (uint8_t(*)[16])in;
  uint8_t(*sout)[16] = (uint8_t(*)[16])out;
  uint64_t tweak = first_tweak;
  for (uint64_t i = 0; i < n; i += 2) {
    fixedkeyaes_prf((FIXEDKEYAES_CTX *)key, sout, sin, tweak);
    fixedkeyaes_prf((FIXEDKEYAES_CTX *)key, sout + 1, sin, tweak + 1);
    tweak += 2;
    ++sin;
    sout += 2;
  }
#else
  static const uint64_t PP = (3 * (1 << PARAM_D));
  uint64_t cur_tweak = first_tweak % PP;
  ASSERT_DRAMATICALLY(n % 2 == 0, "only even first tweaks are supported");
  uint64_t base_tweak = first_tweak - cur_tweak;
  uint64_t tweak = cur_tweak / 2;
  uint64_t ns2 = n / 2;
  treeprg_sha3_context_t *context = (treeprg_sha3_context_t *)key;

  const seed_t *in_seeds = (seed_t *)in;
  seed_t(*out_seeds)[2] = (seed_t(*)[2])out;

  // salt + base_tweak
  uint8_t vec_to_hash[PARAM_salt_size + 4 + PARAM_seed_size]
      __attribute__((aligned(8)));
  uint64_t *vv = (uint64_t *)vec_to_hash;
  memcpy(vec_to_hash, context->salt, PARAM_salt_size);
  *vv += base_tweak;

  for (uint64_t i = 0; i < ns2; ++i) {
    // out[i, i+1] = sha3(base_salt||(tweak+i)||in[i])
    memcpy(vec_to_hash + PARAM_salt_size, &tweak, 4);
    memcpy(vec_to_hash + PARAM_salt_size + 4, in_seeds + i, PARAM_seed_size);
    sdith_hash(out_seeds[i], 2 * PARAM_seed_size, vec_to_hash,
               sizeof(vec_to_hash));
    ++tweak;
  }
#endif
}

/** @brief takes as input n seeds and produces n commitments:
 * out[i]=F_{first_tweak+i}(in[i]) */
EXPORT void sdith_tree_prg_seed_commit(TREE_PRG_CTX *key, void *out, const void *in,
                                 const uint64_t first_tweak, uint64_t n) {
#ifdef AES_COMMITMENT
  // TODO: use avx2
  uint8_t(*sin)[16] = (uint8_t(*)[16])in;
  uint8_t(*sout)[16] = (uint8_t(*)[16])out;
  uint64_t tweak = first_tweak;
  for (uint64_t i = 0; i < n; i++) {
    fixedkeyaes_prf((FIXEDKEYAES_CTX *)key, sout, sin, tweak);
    ++tweak;
    ++sin;
    ++sout;
  }
#else
  static const uint64_t PP = (3 * (1 << PARAM_D));
  uint64_t cur_tweak = first_tweak % PP;
  ASSERT_DRAMATICALLY(n % 2 == 0, "only even first tweaks are supported");
  uint64_t base_tweak = first_tweak - cur_tweak;
  uint64_t tweak = cur_tweak;
  treeprg_sha3_context_t *context = (treeprg_sha3_context_t *)key;

  const seed_t *in_seeds = (seed_t *)in;
  commit_t *out_commits = (commit_t *)out;

  // base_salt = salt + base_tweak
  uint8_t vec_to_hash[PARAM_salt_size + 4 + PARAM_seed_size]
      __attribute__((aligned(8)));
  uint64_t *vv = (uint64_t *)vec_to_hash;
  memcpy(vec_to_hash, context->salt, PARAM_salt_size);
  *vv += base_tweak;

  for (uint64_t i = 0; i < n; ++i) {
    // out[i, i+1] = sha3(base_salt||(tweak+i)||in[i])
    memcpy(vec_to_hash + PARAM_salt_size, &tweak, 4);
    memcpy(vec_to_hash + PARAM_salt_size + 4, in_seeds + i, PARAM_seed_size);
    sdith_hash(out_commits + i, 2 * PARAM_commit_size, vec_to_hash,
               sizeof(vec_to_hash));
  }
#endif
}

/** @brief takes a tree_prg context */
EXPORT TREE_PRG_CTX *sdith_create_tree_prg_ctx(void const *const root_salt) {
#ifdef AES_COMMITMENT
  return (TREE_PRG_CTX *)new_fixedkeyaes_prf(root_salt);
#else
  treeprg_sha3_context_t *context =
      (treeprg_sha3_context_t *)malloc(sizeof(treeprg_sha3_context_t));
  memcpy(context->salt, root_salt, PARAM_salt_size);
  return (TREE_PRG_CTX *)context;
#endif
}

/** @brief free a tree_prg context */
EXPORT void sdith_free_tree_prg_ctx(TREE_PRG_CTX *key) {
#ifdef AES_COMMITMENT
  delete_fixedkeyaes_prf((FIXEDKEYAES_CTX *)key);
#else
  free(key);
#endif
}
