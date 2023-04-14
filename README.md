# Code release for Hypercube-SDitH

To build, go to a particular variant folder under Reference_Implementation and run `make -j`.
It will build a binary called `sign` which will generate the KAT.

This library depends on libXKCP which needs to be initialized in the submodule.

## Variants

This library contains two variants of the Hypercube-SDitH reference implementation.

- Variant 1 (`sdith_hypercube_cat1_gf256_ec`): [The Return of the SDitH](https://eprint.iacr.org/2022/1645). This variant matches the specs of the original Hypercube-SDitH paper.
- Variant 2 (`sdith_hypercube_cat1_gf256_qrom`): Security of Hypercube-SDitH against the Quantum Menace. This variant contains Variant 1 plus the 3-round IDS protocol as well as the Proof-of-Work mechanism introduced in the following paper.
