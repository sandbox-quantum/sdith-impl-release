#include <openssl/evp.h>
#include <stdlib.h>

#include "param.h"

#ifdef SUPERCOP
#include <libkeccak.a.headers/KeccakHash.h>
#else
#include <libXKCP.a.headers/KeccakHash.h>
#endif

#include "rng.h"

HASH_CTX* sdith_hash_create_hash_ctx() {
#ifndef XKCP
  const EVP_MD* md = EVP_shake128();
  EVP_MD_CTX* ctx = EVP_MD_CTX_new();
  EVP_DigestInit_ex(ctx, md, NULL);
  return (HASH_CTX*)ctx;
#else
  Keccak_HashInstance* ctx = (Keccak_HashInstance*) malloc(sizeof(Keccak_HashInstance));
  Keccak_HashInitialize_SHAKE128(ctx);
  return (HASH_CTX*) ctx;
#endif
}

void sdith_hash_free_hash_ctx(HASH_CTX* ctx) {
#ifndef XKCP
  EVP_MD_CTX_free((EVP_MD_CTX*)ctx);
#else
  free(ctx);
#endif
}

void sdith_hash_digest_update(HASH_CTX* ctx, void const* in, int inBytes) {
#ifndef XKCP
  EVP_DigestUpdate((EVP_MD_CTX*)ctx, in, inBytes);
#else
  Keccak_HashUpdate((Keccak_HashInstance*)ctx, (uint8_t const*)in, inBytes << 3);
#endif
}

void sdith_hash_finalize(HASH_CTX* ctx, void* dest, int destBytes) {
#ifndef XKCP
  EVP_DigestFinalXOF((EVP_MD_CTX*)ctx, (unsigned char*)dest, destBytes);
#else
  Keccak_HashFinal((Keccak_HashInstance*)ctx, NULL);
  Keccak_HashSqueeze((Keccak_HashInstance*)ctx, (uint8_t*)dest, destBytes << 3);
#endif
}

void sdith_hash(void* dest, int destBytes, void const* data, int dataBytes) {
  HASH_CTX* ctx = sdith_hash_create_hash_ctx();
  sdith_hash_digest_update(ctx, data, dataBytes);
  sdith_hash_finalize(ctx, dest, destBytes);
  sdith_hash_free_hash_ctx(ctx);
}

