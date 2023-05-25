#include "sdith.h"
#include "timer.h"

#include <string.h>

int main() {
  // Initialize arithmetic libraries and prepopulate tables.
  field_init();

  sdith_bench_timer keygen_timer;
  sdith_bench_timer_init(&keygen_timer);

  sdith_compressed_pubkey_t pk;
  sdith_compressed_key_t sk;
  sdith_bench_timer_start(&keygen_timer);
  keygen(&pk, &sk);
  sdith_bench_timer_end(&keygen_timer);
  sdith_bench_timer_count(&keygen_timer);
  fprintf(stdout, "keygen (ms): %lf\n", sdith_bench_timer_get(&keygen_timer));

  sdith_bench_timer key_expand;
  sdith_bench_timer_init(&key_expand);

  sdith_full_pubkey_t u_pk;
  sdith_full_key_t u_sk;

  sdith_bench_timer_start(&key_expand);
  uncompress_key(&pk, &sk, &u_pk, &u_sk);
  sdith_bench_timer_end(&key_expand);
  sdith_bench_timer_count(&key_expand);
  fprintf(stdout, "key_expand (ms): %lf\n", sdith_bench_timer_get(&key_expand));

  char const msg[] = "Hello, creators of SDitHitH at SandboxAQ.";
  uint8_t signature[sig_size()];
  int sigBytes;

  sign(&u_pk, &u_sk, msg, strlen(msg), signature, &sigBytes);

  sdith_bench_timer verify_timer;
  sdith_bench_timer_init(&verify_timer);
  sdith_bench_timer_start(&verify_timer);
  int res = verify(&u_pk, msg, strlen(msg), signature);
  sdith_bench_timer_end(&verify_timer);
  sdith_bench_timer_count(&verify_timer);
  fprintf(stdout, "verify (ms): %lf\n", sdith_bench_timer_get(&verify_timer));
  fprintf(stdout, "verification result: %d\n", res);

  return 0;
}
