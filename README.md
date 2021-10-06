# cheetah

This crates provide an ongoing implementation of the prime-order curve Cheetah over the field extension $\mathbb{F}_{p^6}$, with p = 2<sup>62</sup> - 111. 2<sup>39</sup> + 1.

* This implementation does not require the Rust standard library
* Arithmetic operations are all constant time unless "_vartime" is explicited mentioned

## Features

* None

## Planned Features

* `serde`

## Description

Cheetah is a STARK-friendly elliptic curve defined over a sextic extension of $\mathbb{F}_p$, p = 2<sup>62</sup> - 111. 2<sup>39</sup> + 1 defined by
`E: y^2 = x^3 + x + B` with 
`B = (4198454498232167043u + 1236584124511164798)v^2 + (109675342022917428u + 2744078705171299893)v + 4212198449936917436u + 2883908786857436727`
where
- `u^2 - u -1 = 0` is the polynomial defining $\mathbb{F}_{p^2} / \mathbb{F}_p$
- `v^3 - v - 2 = 0` is the polynomial defining $\mathbb{F}_{p^6} / \mathbb{F}_{p^2}$.

Cheetah defines a subgroup G of prime order `q` of 372-bits.

Extension towering $\mathbb{F}_{p^2}$ and $\mathbb{F}_{p^6}$ have been specifically constructed with polynomials of small coefficients in order to reduce the cost of multiplication, squaring and inversion in the extension fields. The current implementation may however not be fully optimal with respect to the number of multiplications in the base field.

The Cheetah curve has been generated with the Sagemath utility script `find_curve_extension.sage` available [here](https://github.com/Nashtare/curve_experiments).


## Curve security

Elliptic curves based on extension fields may suffer from specific attacks that do not apply to common elliptic curves constructed over large prime fields and may outperform regular Pollard-Rho attacks, and hence require more scrutiny when evaluating their estimated security level.

## License

Licensed under
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)
