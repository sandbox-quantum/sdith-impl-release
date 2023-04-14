#include "api.h"

#include "gf256.h"
#include "gf2p24.h"
#include "gf2p32.h"
#include "p251.h"
#include "p251p3.h"
#include "p251p4.h"
#include "sdith.h"

#include <string.h>

void fixedkeyaes_prf_init();
int randombytes(unsigned char *x, unsigned long long xlen);

#ifdef SUPERCOP
#define crypto_sign_keypair CRYPTO_NAMESPACE(keypair)
#define crypto_sign_open CRYPTO_NAMESPACE(open)
#define crypto_sign CRYPTO_NAMESPACETOP
#endif

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk) {
  gf256_init(0);
  gf2p24_init();
  gf2p32_mul_init();
  p251_init_(0);
  p251p4_mul_init();
  p251p3_mul_init();
  fixedkeyaes_prf_init();

  sdith_compressed_pubkey_t c_pk;
  sdith_compressed_key_t c_sk;
  keygen(&c_pk, &c_sk);

  sdith_full_pubkey_t u_pk;
  sdith_full_key_t u_sk;
  //memset(&u_pk, 0, sizeof(sdith_full_pubkey_t));
  memset(&u_sk, 0, sizeof(sdith_full_key_t));
  uncompress_key(&c_pk, &c_sk, &u_pk, &u_sk);

  memcpy(pk, &c_pk, sizeof(sdith_compressed_pubkey_t));
  memcpy(sk, &u_sk, sizeof(sdith_full_key_t));

  return 0;
}

int crypto_sign(unsigned char *sm, unsigned long long *smlen,
                const unsigned char *m, unsigned long long mlen,
                const unsigned char *sk) {
  uint8_t signature[sig_size()];
  int sigBytes;

  sdith_full_key_t *u_sk = (sdith_full_key_t*)sk;
  sdith_full_pubkey_t u_pk;
  uncompress_pubkey(&u_sk->compressed_pubkey, &u_pk);
  sign(&u_pk, u_sk, m, mlen, signature, &sigBytes);

  memcpy(sm, signature, sigBytes);
  memcpy(sm + sigBytes, m, mlen);
  *smlen = sigBytes + mlen;

  return 0;
}

int crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                     const unsigned char *sm, unsigned long long smlen,
                     const unsigned char *pk) {
  sdith_full_pubkey_t u_pk;
  uncompress_pubkey((sdith_compressed_pubkey_t*) pk, &u_pk);

  int ret = verify(&u_pk, sm + sig_size(), smlen - sig_size(), sm);
  if (ret) {
    return -1;
  }

  *mlen = smlen - sig_size();
  memcpy(m, sm + sig_size(), *mlen);

  return 0;
}
