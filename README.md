# Cheetah 🐆

This crates provide an implementation of the Cheetah curve over the field extension $\mathbb{F}_{p^6}$, with p = 2<sup>62</sup> + 2<sup>56</sup> + 2<sup>55</sup> + 1.

* This implementation does not require the Rust standard library
* Arithmetic operations are all constant time unless "_vartime" is explicited mentioned

**WARNING:** This is an ongoing, prototype implementation subject to changes. In particular, it has not been audited and may contain bugs and security flaws. This implementation is NOT ready for production use.

## Features

* `serialize` (on by default): Enables Serde serialization

## Description

Cheetah is a STARK-friendly elliptic curve defined over a sextic extension of $\mathbb{F}_p$, p = 2<sup>62</sup> + 2<sup>56</sup> + 2<sup>55</sup> + 1 defined by
`E: y^2 = x^3 + x + B` with 
`B = (1200866201009650596 * u + 1935817186716799185) * v^2 + (3999205700308519553 * u + 3518137720867787056) * v + 2508413708960025374 * u + 1526905369741321712`
where
- `u^2 - 2u - 2 = 0` is the polynomial defining $\mathbb{F}_{p^2} / \mathbb{F}_p$
- `v^3 + v + 1 = 0` is the polynomial defining $\mathbb{F}_{p^6} / \mathbb{F}_{p^2}$.

Cheetah defines a subgroup G of prime order `q = 0x26337f752795f77cb6b6ebb9a18fecc9f2f264f035242b271e13aee130956aa5` of 254-bits.

Extension towering $\mathbb{F}_{p^2}$ and $\mathbb{F}_{p^6}$ have been specifically constructed with polynomials of small coefficients in order to reduce the cost of multiplication, squaring and inversion in the extension fields. The current implementation may however not be fully optimal with respect to the number of multiplications in the base field.

The Cheetah curve has been generated with the Sagemath utility script `find_curve_extension.sage` available [here](https://github.com/Nashtare/curve_experiments).


## Curve security

Elliptic curves based on extension fields may suffer from specific attacks that do not apply to common elliptic curves constructed over large prime fields and may outperform regular Pollard-Rho attacks, and hence require more scrutiny when evaluating their estimated security level.

## License

Licensed under either of

 * Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.
