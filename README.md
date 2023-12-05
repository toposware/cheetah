# Cheetah 🐆

[![codecov](https://codecov.io/gh/toposware/cheetah/branch/main/graph/badge.svg?token=0NWRGPBE8Q)](https://codecov.io/gh/toposware/cheetah)
![example workflow](https://github.com/toposware/cheetah/actions/workflows/ci.yml/badge.svg)

## NOTICE

This repository is unmaintained and therefore publicly archived.
If you have any interest in the Cheetah curve or the work from its associated preprint [eprint.iacr.org/2022/277](https://eprint.iacr.org/2022/277)
on extension-based elliptic curves, feel free to contact the authors at "research [at] toposware [dot] com".

---

This crate provides an implementation of the Cheetah curve over a sextic extension of the prime field `Fp`,
with p = 2<sup>64</sup> - 2<sup>32</sup> + 1.

* This implementation can be made `no_std` by relying on the `alloc` crate instead.
* Arithmetic operations are all constant time unless "_vartime" is explicited mentioned

**WARNING:** This is an ongoing, prototype implementation subject to changes. In particular, it has not been audited and may contain bugs and security flaws. This implementation is NOT ready for production use.

## Features

* `std` (on by default): Enables use of the Rust standard library.
* `serialize` (on by default): Enables Serde serialization.

## Description

Cheetah is a STARK-friendly elliptic curve defined over a sextic extension of `Fp`, p = 2<sup>64</sup> - 2<sup>32</sup> + 1 defined by
`E: y^2 = x^3 + x + B` with
`B = u + 395`
where

* `u^6 - 7 = 0` is the polynomial defining `Fp6 / Fp`

Cheetah defines a subgroup G of prime order
<p align="center">
q = 55610362957290864006699123731285679659474893560816383126640993521607086746831
</p>
of 255-bits.

The extension `Fp6` has been specifically constructed with a sparse polynomial of the form `X^6 - A`, where `A` is a small quadratic and cubic non-residue. The current implementation may however not be fully optimal with respect to the number of multiplications in the base field.

The Cheetah curve has been generated with the Sagemath utility script `sextic_search.sage` available [here](https://github.com/Toposware/cheetah_evidence).

## Curve security

Elliptic curves based on extension fields may suffer from specific attacks that do not apply to common elliptic curves constructed over large prime fields and may outperform regular Pollard-Rho attacks, and hence require more scrutiny when evaluating their estimated security level. To verify the security level
of Cheetah against generic attacks as well as cover and decomposition attacks, please use the Sagemath utility script `verify.sage` available
[here](https://github.com/Toposware/cheetah_evidence).

## License

Licensed under either of

* Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
* MIT license ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>)

at your option.
