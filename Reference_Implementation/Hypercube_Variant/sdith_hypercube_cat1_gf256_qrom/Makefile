# This Makefile builds sdith_hypercube_<SEC_LEVEL>_<VARIANT>.
# SEC_LEVEL: CAT_1|CAT_3|CAT_5
# VARIANT: generic64|AVX2

.POSIX:

VARIANT = generic64
SEC_LEVEL = CAT_1

ifeq ($(VARIANT), AVX2)
	AVXFLAGS = -mavx2 -mpclmul -mgfni -mavx -maes
else
	AVXFLAGS =
endif
CC = gcc
CFLAGS = -W -Wall -O3 -fPIC -D${SEC_LEVEL}

XKCP_PATH = ../../../Library/XKCP
KAT_MAIN_PATH = ../../../KAT/generator

ifeq ($(VARIANT), AVX2)
	GF256_SRC = gf2p24.c gf2p32.c gf256.c gf256-avx2-polytable-ct.c gf256-avx2.c gf256-avx2-gfni.c gf256-avx-pclmul.c
	P251_SRC = p251.c p251p3.c p251p4.c p251-avx2-ct.c
	AES_SRC = aes256ctr-avx.c
else
	GF256_SRC = gf2p24.c gf2p32.c gf256.c
	P251_SRC = p251.c p251p3.c p251p4.c
	AES_SRC = aes256ctr-ref.c
endif
CRYPTO_SRC = hash-sha3.c rng.c treeprg.c
SDITH_SRC = sdith.c sign.c
KAT_MAIN_SRC = ${KAT_MAIN_PATH}/PQCgenKAT_sign.c ${KAT_MAIN_PATH}/rng.c

all: sign

libXKCP.a:
	-cd ${XKCP_PATH} && ${MAKE} ${VARIANT}/libXKCP.a

sign: libXKCP.a ${SRC} ${AES_SRC} ${GF256_SRC} ${P251_SRC} ${CRYPTO_SRC} ${SDITH_SRC} ${KAT_MAIN_SRC}
	-${CC} ${CFLAGS} ${AVXFLAGS} -o $@ ${SRC} ${AES_SRC} ${GF256_SRC} ${P251_SRC} ${CRYPTO_SRC} ${SDITH_SRC} ${KAT_MAIN_SRC} ${XKCP_PATH}/bin/${VARIANT}/libXKCP.a -I${XKCP_PATH}/bin/${VARIANT} -I. -lcrypto

clean:
	-rm -f *.o *.a sign
	-cd ${XKCP_PATH} && ${MAKE} clean