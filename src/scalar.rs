// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This module provides an implementation of the scalar field Fq
//! of the STARK-friendly cheetah curve, with characteristic
//! q = 0x7af2599b3b3f22d0563fbf0f990a37b5327aa72330157722d443623eaed4accf.

use core::{
    borrow::Borrow,
    convert::{TryFrom, TryInto},
    fmt::{self, Debug, Display, Formatter},
    hash::{Hash, Hasher},
    iter::Sum,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::utils::{
    add64_with_carry, mul64_with_carry, shl64_by_u32_with_carry, square_assign_multi,
    sub64_with_carry,
};

use bitvec::{order::Lsb0, slice::BitSlice};
use group::ff::{Field, PrimeField};
use rand_core::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "serialize")]
use serde::de::Visitor;
#[cfg(feature = "serialize")]
use serde::{self, Deserialize, Deserializer, Serialize, Serializer};

// CONSTANTS
// ================================================================================================

// Field modulus = 55610362957290864006699123731285679659474893560816383126640993521607086746831
const M: Scalar = Scalar([
    0xd443623eaed4accf,
    0x327aa72330157722,
    0x563fbf0f990a37b5,
    0x7af2599b3b3f22d0,
]);

// 2^256 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R: Scalar = Scalar([
    0x57793b82a256a662,
    0x9b0ab1b99fd511ba,
    0x538081e0cdeb9095,
    0x0a1b4cc98981ba5f,
]);

// 2^512 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R2: Scalar = Scalar([
    0x5a93b1562a974d84,
    0x5a5314649f39eecd,
    0xa2780be8a30fdfc6,
    0x2a68a265fd96d1bb,
]);

// 2^768 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R3: Scalar = Scalar([
    0x2d80dd65de381f7a,
    0x3a30adaf7b70fb37,
    0x5154791dc29dc568,
    0x7219959569912364,
]);

// Multiplicative generator g of order q-1
// g = 6
//   = 0x3ca3ccb9390a5e3bf5030b44d3856381a2402a59befe6a5e0cd7650fce07e64c in Montgomery form
const GENERATOR: Scalar = Scalar([
    0x0cd7650fce07e64c,
    0xa2402a59befe6a5e,
    0xf5030b44d3856381,
    0x3ca3ccb9390a5e3b,
]);

// Two-adicity of the field: (q-1) % 2^1 = 0
const TWO_ADICITY: u32 = 1;

// 2^1 root of unity = 0x7af2599b3b3f22d0563fbf0f990a37b5327aa72330157722d443623eaed4acce
//                   = 0x70d70cd1b1bd687102bf3d2ecb1ea71f976ff569904065687cca26bc0c7e066d in Montgomery form
const TWO_ADIC_ROOT_OF_UNITY: Scalar = Scalar([
    0x7cca26bc0c7e066d,
    0x976ff56990406568,
    0x02bf3d2ecb1ea71f,
    0x70d70cd1b1bd6871,
]);

// -M^{-1} mod 2^64; this is used during element multiplication.
const U: u64 = 7208734935082542545;

// SCALAR FIELD ELEMENT
// ================================================================================================

/// Represents a scalar field element.
///
/// Internal values are stored in Montgomery representation.
/// The backing type is `[u64; 4]`.
#[derive(Copy, Clone, Eq)]
pub struct Scalar(pub(crate) [u64; 4]);

impl Default for Scalar {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

impl Debug for Scalar {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        let tmp = self.to_bytes();
        write!(f, "0x")?;
        for &b in tmp.iter().rev() {
            write!(f, "{:02x}", b)?;
        }
        Ok(())
    }
}

impl Display for Scalar {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl Hash for Scalar {
    fn hash<H: Hasher>(&self, hasher: &mut H) {
        self.0.hash(hasher);
    }
}

impl ConstantTimeEq for Scalar {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0[0].ct_eq(&other.0[0])
            & self.0[1].ct_eq(&other.0[1])
            & self.0[2].ct_eq(&other.0[2])
            & self.0[3].ct_eq(&other.0[3])
    }
}

impl PartialEq for Scalar {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl ConditionallySelectable for Scalar {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Scalar([
            u64::conditional_select(&a.0[0], &b.0[0], choice),
            u64::conditional_select(&a.0[1], &b.0[1], choice),
            u64::conditional_select(&a.0[2], &b.0[2], choice),
            u64::conditional_select(&a.0[3], &b.0[3], choice),
        ])
    }
}

impl zeroize::DefaultIsZeroes for Scalar {}

impl Scalar {
    /// Creates a new field element from a [u64; 4] value.
    /// The value is converted to Montgomery form by computing
    /// (a.R^0 * R^2) / R = a.R
    pub const fn new(value: [u64; 4]) -> Self {
        (&Scalar(value)).mul(&R2)
    }

    /// Returns zero, the additive identity.
    #[inline]
    pub const fn zero() -> Self {
        Scalar([0, 0, 0, 0])
    }

    /// Returns one, the multiplicative identity.
    #[inline]
    pub const fn one() -> Self {
        R
    }

    /// Checks whether `self` is zero or not
    pub fn is_zero(&self) -> Choice {
        self.ct_eq(&Scalar::zero())
    }

    /// Generates a random element
    pub fn random(mut rng: impl RngCore) -> Self {
        let mut buf = [0; 64];
        rng.fill_bytes(&mut buf);
        Self::from_bytes_wide(&buf)
    }

    #[inline(always)]
    #[allow(clippy::too_many_arguments)]
    const fn montgomery_reduce(
        r0: u64,
        r1: u64,
        r2: u64,
        r3: u64,
        r4: u64,
        r5: u64,
        r6: u64,
        r7: u64,
    ) -> Self {
        let k = r0.wrapping_mul(U);
        let (_, carry) = mul64_with_carry(r0, k, M.0[0], 0);
        let (r1, carry) = mul64_with_carry(r1, k, M.0[1], carry);
        let (r2, carry) = mul64_with_carry(r2, k, M.0[2], carry);
        let (r3, carry) = mul64_with_carry(r3, k, M.0[3], carry);
        let (r4, carry2) = add64_with_carry(r4, 0, carry);

        let k = r1.wrapping_mul(U);
        let (_, carry) = mul64_with_carry(r1, k, M.0[0], 0);
        let (r2, carry) = mul64_with_carry(r2, k, M.0[1], carry);
        let (r3, carry) = mul64_with_carry(r3, k, M.0[2], carry);
        let (r4, carry) = mul64_with_carry(r4, k, M.0[3], carry);
        let (r5, carry2) = add64_with_carry(r5, carry2, carry);

        let k = r2.wrapping_mul(U);
        let (_, carry) = mul64_with_carry(r2, k, M.0[0], 0);
        let (r3, carry) = mul64_with_carry(r3, k, M.0[1], carry);
        let (r4, carry) = mul64_with_carry(r4, k, M.0[2], carry);
        let (r5, carry) = mul64_with_carry(r5, k, M.0[3], carry);
        let (r6, carry2) = add64_with_carry(r6, carry2, carry);

        let k = r3.wrapping_mul(U);
        let (_, carry) = mul64_with_carry(r3, k, M.0[0], 0);
        let (r4, carry) = mul64_with_carry(r4, k, M.0[1], carry);
        let (r5, carry) = mul64_with_carry(r5, k, M.0[2], carry);
        let (r6, carry) = mul64_with_carry(r6, k, M.0[3], carry);
        let (r7, _) = add64_with_carry(r7, carry2, carry);

        // The result may be within M of the correct value,
        // hence substracting the modulus
        (&Scalar([r4, r5, r6, r7])).sub(&M)
    }

    /// Computes the summation of two scalar elements
    #[inline]
    pub const fn add(&self, rhs: &Self) -> Self {
        let (d0, carry) = add64_with_carry(self.0[0], rhs.0[0], 0);
        let (d1, carry) = add64_with_carry(self.0[1], rhs.0[1], carry);
        let (d2, carry) = add64_with_carry(self.0[2], rhs.0[2], carry);
        let (d3, _) = add64_with_carry(self.0[3], rhs.0[3], carry);

        // Attempt to subtract the modulus, to ensure the value
        // is smaller than the modulus.
        (&Scalar([d0, d1, d2, d3])).sub(&M)
    }

    /// Computes the difference of two scalar elements
    #[inline]
    pub const fn sub(&self, rhs: &Self) -> Self {
        let (d0, borrow) = sub64_with_carry(self.0[0], rhs.0[0], 0);
        let (d1, borrow) = sub64_with_carry(self.0[1], rhs.0[1], borrow);
        let (d2, borrow) = sub64_with_carry(self.0[2], rhs.0[2], borrow);
        let (d3, borrow) = sub64_with_carry(self.0[3], rhs.0[3], borrow);

        // If underflow occurred on the final limb,
        // borrow = 0xfff...fff, otherwise borrow = 0x000...000.
        let (d0, carry) = add64_with_carry(d0, M.0[0] & borrow, 0);
        let (d1, carry) = add64_with_carry(d1, M.0[1] & borrow, carry);
        let (d2, carry) = add64_with_carry(d2, M.0[2] & borrow, carry);
        let (d3, _) = add64_with_carry(d3, M.0[3] & borrow, carry);

        Scalar([d0, d1, d2, d3])
    }

    /// Computes the negation of a scalar element
    #[inline]
    pub const fn neg(&self) -> Self {
        // Subtract `self` from `M` to negate. Ignore the final
        // borrow because it cannot underflow; self is guaranteed to
        // be in the field.
        let (d0, borrow) = sub64_with_carry(M.0[0], self.0[0], 0);
        let (d1, borrow) = sub64_with_carry(M.0[1], self.0[1], borrow);
        let (d2, borrow) = sub64_with_carry(M.0[2], self.0[2], borrow);
        let (d3, _) = sub64_with_carry(M.0[3], self.0[3], borrow);

        // `tmp` could be `M` if `self` was zero. Create a mask that is
        // zero if `self` was zero, and `u64::max_value()` if self was nonzero.
        let mask = (((self.0[0] | self.0[1] | self.0[2] | self.0[3]) == 0) as u64).wrapping_sub(1);

        Scalar([d0 & mask, d1 & mask, d2 & mask, d3 & mask])
    }

    /// Computes the multiplication of two scalar elements
    #[inline]
    pub const fn mul(&self, rhs: &Self) -> Self {
        // Schoolbook multiplication

        let (r0, carry) = mul64_with_carry(0, self.0[0], rhs.0[0], 0);
        let (r1, carry) = mul64_with_carry(0, self.0[0], rhs.0[1], carry);
        let (r2, carry) = mul64_with_carry(0, self.0[0], rhs.0[2], carry);
        let (r3, r4) = mul64_with_carry(0, self.0[0], rhs.0[3], carry);

        let (r1, carry) = mul64_with_carry(r1, self.0[1], rhs.0[0], 0);
        let (r2, carry) = mul64_with_carry(r2, self.0[1], rhs.0[1], carry);
        let (r3, carry) = mul64_with_carry(r3, self.0[1], rhs.0[2], carry);
        let (r4, r5) = mul64_with_carry(r4, self.0[1], rhs.0[3], carry);

        let (r2, carry) = mul64_with_carry(r2, self.0[2], rhs.0[0], 0);
        let (r3, carry) = mul64_with_carry(r3, self.0[2], rhs.0[1], carry);
        let (r4, carry) = mul64_with_carry(r4, self.0[2], rhs.0[2], carry);
        let (r5, r6) = mul64_with_carry(r5, self.0[2], rhs.0[3], carry);

        let (r3, carry) = mul64_with_carry(r3, self.0[3], rhs.0[0], 0);
        let (r4, carry) = mul64_with_carry(r4, self.0[3], rhs.0[1], carry);
        let (r5, carry) = mul64_with_carry(r5, self.0[3], rhs.0[2], carry);
        let (r6, r7) = mul64_with_carry(r6, self.0[3], rhs.0[3], carry);

        Scalar::montgomery_reduce(r0, r1, r2, r3, r4, r5, r6, r7)
    }

    /// Computes the square of a scalar element
    #[inline]
    pub const fn square(&self) -> Self {
        let (r1, carry) = mul64_with_carry(0, self.0[0], self.0[1], 0);
        let (r2, carry) = mul64_with_carry(0, self.0[0], self.0[2], carry);
        let (r3, r4) = mul64_with_carry(0, self.0[0], self.0[3], carry);

        let (r3, carry) = mul64_with_carry(r3, self.0[1], self.0[2], 0);
        let (r4, r5) = mul64_with_carry(r4, self.0[1], self.0[3], carry);

        let (r5, r6) = mul64_with_carry(r5, self.0[2], self.0[3], 0);

        let r7 = r6 >> 63;
        let r6 = (r6 << 1) | (r5 >> 63);
        let r5 = (r5 << 1) | (r4 >> 63);
        let r4 = (r4 << 1) | (r3 >> 63);
        let r3 = (r3 << 1) | (r2 >> 63);
        let r2 = (r2 << 1) | (r1 >> 63);
        let r1 = r1 << 1;

        let (r0, carry) = mul64_with_carry(0, self.0[0], self.0[0], 0);
        let (r1, carry) = add64_with_carry(0, r1, carry);
        let (r2, carry) = mul64_with_carry(r2, self.0[1], self.0[1], carry);
        let (r3, carry) = add64_with_carry(0, r3, carry);
        let (r4, carry) = mul64_with_carry(r4, self.0[2], self.0[2], carry);
        let (r5, carry) = add64_with_carry(0, r5, carry);
        let (r6, carry) = mul64_with_carry(r6, self.0[3], self.0[3], carry);
        let (r7, _) = add64_with_carry(0, r7, carry);

        Scalar::montgomery_reduce(r0, r1, r2, r3, r4, r5, r6, r7)
    }

    /// Computes the square root of this element, if it exists.
    pub fn sqrt(&self) -> CtOption<Self> {
        // Tonelli-Shank's algorithm for q mod 4 = 3
        // See https://eprint.iacr.org/2012/685.pdf

        // Compute s^((q+1)/4)
        let s = self.exp_vartime(&[
            0xb510d88fabb52b34,
            0x4c9ea9c8cc055dc8,
            0x158fefc3e6428ded,
            0x1ebc9666cecfc8b4,
        ]);

        CtOption::new(s, (s * s).ct_eq(self))
    }

    /// Computes the double of a scalar element
    #[inline]
    pub const fn double(&self) -> Self {
        let (d0, carry) = shl64_by_u32_with_carry(self.0[0], 1, 0);
        let (d1, carry) = shl64_by_u32_with_carry(self.0[1], 1, carry);
        let (d2, carry) = shl64_by_u32_with_carry(self.0[2], 1, carry);
        let (d3, _carry) = shl64_by_u32_with_carry(self.0[3], 1, carry);

        // Attempt to subtract the modulus, to ensure the value
        // is smaller than the modulus.
        (&Scalar([d0, d1, d2, d3])).sub(&M)
    }

    /// Attempts to convert a little-endian byte representation of
    /// a scalar into a `Scalar`, failing if the input is not canonical.
    pub fn from_bytes(bytes: &[u8; 32]) -> CtOption<Scalar> {
        let mut tmp = Scalar([0, 0, 0, 0]);

        tmp.0[0] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[0..8]).unwrap());
        tmp.0[1] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[8..16]).unwrap());
        tmp.0[2] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[16..24]).unwrap());
        tmp.0[3] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[24..32]).unwrap());

        // Try to subtract the modulus
        let (_, borrow) = sub64_with_carry(tmp.0[0], M.0[0], 0);
        let (_, borrow) = sub64_with_carry(tmp.0[1], M.0[1], borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[2], M.0[2], borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[3], M.0[3], borrow);

        // If the element is smaller than M then the
        // subtraction will underflow, producing a borrow value
        // of 0xffff...ffff. Otherwise, it'll be zero.
        let is_some = (borrow as u8) & 1;

        // Convert to Montgomery form by computing
        // (a.R^0 * R^2) / R = a.R
        tmp *= &R2;

        CtOption::new(tmp, Choice::from(is_some))
    }

    /// Converts a 256-bit little endian integer into
    /// a `Scalar`. The element does ***NOT*** have to
    /// be canonical.
    /// This is to be used when providing an unconstrained
    /// slice of bytes, for instance a hash digest.
    pub fn from_bytes_non_canonical(bytes: &[u8; 32]) -> Self {
        let mut tmp = Scalar([0, 0, 0, 0]);

        tmp.0[0] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[0..8]).unwrap());
        tmp.0[1] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[8..16]).unwrap());
        tmp.0[2] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[16..24]).unwrap());
        tmp.0[3] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[24..32]).unwrap());

        // Convert to Montgomery form by computing
        // (a.R^0 * R^2) / R = a.R
        tmp * R2
    }

    /// Converts a little-endian bit sequence into a Scalar element
    pub fn from_bits(bit_slice: &BitSlice<Lsb0, u8>) -> Scalar {
        assert_eq!(bit_slice.len(), 256);

        let mut result = Scalar::zero();
        for i in (0..256).rev() {
            result = result.double();
            let tmp = Scalar::conditional_select(
                &Scalar::zero(),
                &Scalar::one(),
                Choice::from(bit_slice[i] as u8),
            );
            result += tmp;
        }

        result
    }

    /// Converts a little-endian bit sequence into a Scalar element
    ///
    /// **This operation is variable time with respect
    /// to the binary slice.** If the slice is fixed,
    /// this operation is effectively constant time.
    pub fn from_bits_vartime(bit_slice: &BitSlice<Lsb0, u8>) -> Scalar {
        assert_eq!(bit_slice.len(), 256);

        let mut result = Scalar::zero();
        for i in (0..256).rev() {
            result = result.double();
            if bit_slice[i] {
                result += Scalar::one();
            }
        }

        result
    }

    /// Outputs the internal representation as 4 64-bit limbs after Montgomery reduction
    pub const fn output_reduced_limbs(&self) -> [u64; 4] {
        Scalar::montgomery_reduce(self.0[0], self.0[1], self.0[2], self.0[3], 0, 0, 0, 0).0
    }

    /// Outputs the internal representation as 4 64-bit limbs without Montgomery reduction
    /// This is intended for uses like re-interpreting the type containing the internal value.
    pub const fn output_unreduced_limbs(&self) -> [u64; 4] {
        self.0
    }

    /// Converts a `Scalar` element into a byte representation in
    /// little-endian byte order.
    pub fn to_bytes(&self) -> [u8; 32] {
        // Turn into canonical form by computing
        // (a.R) / R = a
        let tmp = Scalar::montgomery_reduce(self.0[0], self.0[1], self.0[2], self.0[3], 0, 0, 0, 0);

        let mut res = [0; 32];
        res[0..8].copy_from_slice(&tmp.0[0].to_le_bytes());
        res[8..16].copy_from_slice(&tmp.0[1].to_le_bytes());
        res[16..24].copy_from_slice(&tmp.0[2].to_le_bytes());
        res[24..32].copy_from_slice(&tmp.0[3].to_le_bytes());

        res
    }

    /// Converts a 512-bit little endian integer into
    /// a `Scalar` by reducing by the modulus.
    pub fn from_bytes_wide(bytes: &[u8; 64]) -> Self {
        Scalar::from_u512([
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[0..8]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[8..16]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[16..24]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[24..32]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[32..40]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[40..48]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[48..56]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[56..64]).unwrap()),
        ])
    }

    fn from_u512(limbs: [u64; 8]) -> Self {
        // We reduce an arbitrary 512-bit number by decomposing it into two 256-bit digits
        // with the higher bits multiplied by 2^256. Thus, we perform two reductions
        //
        // 1. the lower bits are multiplied by R^2, as normal
        // 2. the upper bits are multiplied by R^2 * 2^256 = R^3
        //
        // and computing their sum in the field.
        // The reduction works so long as the product is less than R=2^256 multiplied by
        // the modulus. This holds because for any `c` smaller than the modulus, we have
        // that (2^256 - 1)*c is an acceptable product for the reduction. Therefore, the
        // reduction always works so long as `c` is in the field; in this case it is either the
        // constant `R2` or `R3`.
        let d0 = Scalar([limbs[0], limbs[1], limbs[2], limbs[3]]);
        let d1 = Scalar([limbs[4], limbs[5], limbs[6], limbs[7]]);

        // Convert to Montgomery form
        d0 * R2 + d1 * R3
    }

    /// Converts a `Scalar` element given as byte representation into a radix-16
    /// representation, where each resulting coefficient is in [-8; 8).
    ///
    /// The resulting decomposition `[a_0, ..., a_63]` is such that
    /// `sum(a_j * 2^(j * 4)) == a`.
    pub(crate) fn bytes_to_radix_16(bytes: &[u8; 32]) -> [i8; 64] {
        let mut result = [0i8; 64];

        // Convert from bytes to radix-16
        for i in 0..32 {
            result[2 * i] = (bytes[i] & 0xf) as i8;
            result[2 * i + 1] = ((bytes[i] >> 4) & 0xf) as i8;
        }

        // Shift every coefficients from [0; 16) to [-8; 8)
        for i in 0..63 {
            let carry = (result[i] + 8) >> 4;
            result[i] -= carry << 4;
            result[i + 1] += carry;
        }

        result
    }

    /// Converts a `Scalar` element given as byte representation into a w-NAF
    /// representation, where each resulting coefficient is odd and in (-2^(w-1); 2^(w-1)).
    /// In addition, the leading coefficient is non-zero, and there cannot be
    /// more than one non-zero coefficient in any w consecutive set of coefficients.
    ///
    /// **This operation is variable time with respect to the scalar.**
    /// If the scalar is fixed, this operation is effectively constant time.
    pub(crate) fn bytes_to_wnaf_vartime(bytes: &[u8; 32], w: usize) -> [i8; 256] {
        // Taken from https://github.com/dalek-cryptography/curve25519-dalek/blob/main/src/scalar.rs
        // from an adaptation of Algorithm 3.35 in Guide to Elliptic Curve Cryptography by
        // Hankerson, Menezes and Vanstone.

        debug_assert!(w >= 2);
        debug_assert!(w <= 8);

        let mut naf = [0i8; 256];

        let mut x_u64 = [0u64; 5];

        x_u64[0] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[0..8]).unwrap());
        x_u64[1] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[8..16]).unwrap());
        x_u64[2] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[16..24]).unwrap());
        x_u64[3] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[24..32]).unwrap());

        let width = 1 << w;
        let window_mask = width - 1;

        let mut pos = 0;
        let mut carry = 0;
        while pos < 256 {
            // Construct a buffer of bits of the scalar, starting at bit `pos`
            let u64_idx = pos / 64;
            let bit_idx = pos % 64;
            let bit_buf = if bit_idx < 64 - w {
                // This window's bits are contained in a single u64
                x_u64[u64_idx] >> bit_idx
            } else {
                // Combine the current u64's bits with the bits from the next u64
                (x_u64[u64_idx] >> bit_idx) | (x_u64[1 + u64_idx] << (64 - bit_idx))
            };

            // Add the carry into the current window
            let window = carry + (bit_buf & window_mask);

            if window & 1 == 0 {
                // If the window value is even, preserve the carry and continue.
                // Why is the carry preserved?
                // If carry == 0 and window & 1 == 0, then the next carry should be 0
                // If carry == 1 and window & 1 == 0, then bit_buf & 1 == 1 so the next carry should be 1
                pos += 1;
                continue;
            }

            if window < width / 2 {
                carry = 0;
                naf[pos] = window as i8;
            } else {
                carry = 1;
                naf[pos] = (window as i8).wrapping_sub(width as i8);
            }

            pos += w;
        }

        naf
    }

    /// Returns whether or not this element is strictly lexicographically
    /// larger than its negation.
    pub fn lexicographically_largest(&self) -> Choice {
        // This can be determined by checking to see if the element is
        // larger than (M - 1) // 2. If we subtract by ((M - 1) // 2) + 1
        // and there is no underflow, then the element must be larger than
        // (M - 1) // 2.

        // First, because self is in Montgomery form we need to reduce it
        let tmp = Scalar::montgomery_reduce(self.0[0], self.0[1], self.0[2], self.0[3], 0, 0, 0, 0);
        let (_, borrow) = sub64_with_carry(tmp.0[0], 0x6a21b11f576a5668, 0);
        let (_, borrow) = sub64_with_carry(tmp.0[1], 0x993d5391980abb91, borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[2], 0x2b1fdf87cc851bda, borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[3], 0x3d792ccd9d9f9168, borrow);

        // If the element was smaller, the subtraction will underflow
        // producing a borrow value of 0xffff...ffff, otherwise it will
        // be zero. We create a Choice representing true if there was
        // overflow (and so this element is not lexicographically larger
        // than its negation) and then negate it.

        !Choice::from((borrow as u8) & 1)
    }

    /// Exponentiates `self` by `power`, where `power` is a
    /// little-endian order integer exponent.
    pub fn exp(self, power: &[u64; 4]) -> Self {
        let mut res = Self::one();
        for e in power.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();
                let mut tmp = res;
                tmp *= self;
                res.conditional_assign(&tmp, (((*e >> i) & 1) as u8).into());
            }
        }
        res
    }

    /// Exponentiates `self` by `power`, where `power` is a
    /// little-endian order integer exponent.
    ///
    /// **This operation is variable time with respect
    /// to the exponent.** If the exponent is fixed,
    /// this operation is effectively constant time.
    pub fn exp_vartime(&self, power: &[u64; 4]) -> Self {
        let mut res = Self::one();
        for e in power.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();

                if ((*e >> i) & 1) == 1 {
                    res.mul_assign(self);
                }
            }
        }
        res
    }

    /// Computes the multiplicative inverse of this element,
    /// failing if the element is zero.
    pub fn invert(&self) -> CtOption<Self> {
        // found using https://github.com/kwantam/addchain for M - 2
        let mut t0 = self.square(); //    1: 2
        let t13 = t0.square(); //    2: 4
        let t4 = t13 * self; //    3: 5
        let t12 = t4 * t0; //    4: 7
        let t14 = t4 * t13; //    5: 9
        let t9 = t12 * t13; //    6: 11
        let t1 = t14 * t13; //    7: 13
        let t6 = t9 * t13; //    8: 15
        let t7 = t1 * t13; //    9: 17
        let t2 = t6 * t13; //   10: 19
        let t3 = t7 * t13; //   11: 21
        let t15 = t2 * t13; //   12: 23
        let t10 = t3 * t13; //   13: 25
        let t5 = t15 * t13; //   14: 27
        let t11 = t10 * t13; //   15: 29
        t0 = t6.square(); //   16: 30
        let t13 = t5 * t13; //   17: 31
        square_assign_multi(&mut t0, 5); //   22: 960
        t0 *= t15; //   23: 983
        square_assign_multi(&mut t0, 4); //   27: 15728
        t0 *= t14; //   28: 15737
        square_assign_multi(&mut t0, 6); //   34: 1007168
        t0 *= t9; //   35: 1007179
        square_assign_multi(&mut t0, 7); //   42: 128918912
        t0 *= t10; //   43: 128918937
        square_assign_multi(&mut t0, 4); //   47: 2062702992
        t0 *= t9; //   48: 2062703003
        square_assign_multi(&mut t0, 7); //   55: 264025984384
        t0 *= t11; //   56: 264025984413
        square_assign_multi(&mut t0, 5); //   61: 8448831501216
        t0 *= t2; //   62: 8448831501235
        square_assign_multi(&mut t0, 4); //   66: 135181304019760
        t0 *= t6; //   67: 135181304019775
        square_assign_multi(&mut t0, 7); //   74: 17303206914531200
        t0 *= t7; //   75: 17303206914531217
        square_assign_multi(&mut t0, 5); //   80: 553702621264998944
        t0 *= t1; //   81: 553702621264998957
        square_assign_multi(&mut t0, 10); //   91: 566991484175358931968
        t0 *= t3; //   92: 566991484175358931989
        square_assign_multi(&mut t0, 5); //   97: 18143727493611485823648
        t0 *= t7; //   98: 18143727493611485823665
        square_assign_multi(&mut t0, 5); //  103: 580599279795567546357280
        t0 *= t13; //  104: 580599279795567546357311
        square_assign_multi(&mut t0, 4); //  108: 9289588476729080741716976
        t0 *= t9; //  109: 9289588476729080741716987
        square_assign_multi(&mut t0, 4); //  113: 148633415627665291867471792
        t0 *= t6; //  114: 148633415627665291867471807
        square_assign_multi(&mut t0, 9); //  123: 76100308801364629436145565184
        t0 *= t13; //  124: 76100308801364629436145565215
        square_assign_multi(&mut t0, 7); //  131: 9740839526574672567826632347520
        t0 *= t10; //  132: 9740839526574672567826632347545
        square_assign_multi(&mut t0, 7); //  139: 1246827459401558088681808940485760
        t0 *= t4; //  140: 1246827459401558088681808940485765
        square_assign_multi(&mut t0, 8); //  148: 319187829606798870702543088764355840
        t0 *= t5; //  149: 319187829606798870702543088764355867
        square_assign_multi(&mut t0, 5); //  154: 10214010547417563862481378840459387744
        t0 *= t5; //  155: 10214010547417563862481378840459387771
        square_assign_multi(&mut t0, 4); //  159: 163424168758681021799702061447350204336
        t0 *= t4; //  160: 163424168758681021799702061447350204341
        square_assign_multi(&mut t0, 7); //  167: 20918293601111170790361863865260826155648
        t0 *= t10; //  168: 20918293601111170790361863865260826155673
        square_assign_multi(&mut t0, 6); //  174: 1338770790471114930583159287376692873963072
        t0 *= t6; //  175: 1338770790471114930583159287376692873963087
        square_assign_multi(&mut t0, 6); //  181: 85681330590151355557322194392108343933637568
        t0 *= t3; //  182: 85681330590151355557322194392108343933637589
        square_assign_multi(&mut t0, 5); //  187: 2741802578884843377834310220547467005876402848
        t0 *= t12; //  188: 2741802578884843377834310220547467005876402855
        square_assign_multi(&mut t0, 7); //  195: 350950730097259952362791708230075776752179565440
        t0 *= t7; //  196: 350950730097259952362791708230075776752179565457
        square_assign_multi(&mut t0, 5); //  201: 11230423363112318475609334663362424856069746094624
        t0 *= t2; //  202: 11230423363112318475609334663362424856069746094643
        square_assign_multi(&mut t0, 12); //  214: 45999814095308056476095834781132492210461680003657728
        t0 *= t3; //  215: 45999814095308056476095834781132492210461680003657749
        square_assign_multi(&mut t0, 6); //  221: 2943988102099715614470133425992479501469547520234095936
        t0 *= t11; //  222: 2943988102099715614470133425992479501469547520234095965
        square_assign_multi(&mut t0, 5); //  227: 94207619267190899663044269631759344047025520647491070880
        t0 *= t10; //  228: 94207619267190899663044269631759344047025520647491070905
        square_assign_multi(&mut t0, 7); //  235: 12058575266200435156869666512865196038019266642878857075840
        t0 *= t9; //  236: 12058575266200435156869666512865196038019266642878857075851
        square_assign_multi(&mut t0, 4); //  240: 192937204259206962509914664205843136608308266286061713213616
        t0 *= t4; //  241: 192937204259206962509914664205843136608308266286061713213621
        square_assign_multi(&mut t0, 4); //  245: 3086995268147311400158634627293490185732932260576987411417936
        t0 *= self; //  246: 3086995268147311400158634627293490185732932260576987411417937
        square_assign_multi(&mut t0, 9); //  255: 1580541577291423436881220929174266975095261317415417554645983744
        t0 *= t5; //  256: 1580541577291423436881220929174266975095261317415417554645983771
        square_assign_multi(&mut t0, 8); //  264: 404618643786604399841592557868612345624386897258346893989371845376
        t0 *= t7; //  265: 404618643786604399841592557868612345624386897258346893989371845393
        square_assign_multi(&mut t0, 4); //  269: 6473898300585670397465480925897797529990190356133550303829949526288
        t0 *= t6; //  270: 6473898300585670397465480925897797529990190356133550303829949526303
        square_assign_multi(&mut t0, 6); //  276: 414329491237482905437790779257459041919372182792547219445116769683392
        t0 *= t3; //  277: 414329491237482905437790779257459041919372182792547219445116769683413
        square_assign_multi(&mut t0, 5); //  282: 13258543719599452974009304936238689341419909849361511022243736629869216
        t0 *= t5; //  283: 13258543719599452974009304936238689341419909849361511022243736629869243
        square_assign_multi(&mut t0, 4); //  287: 212136699513591247584148878979819029462718557589784176355899786077907888
        t0 *= t4; //  288: 212136699513591247584148878979819029462718557589784176355899786077907893
        square_assign_multi(&mut t0, 7); //  295: 27153497537739679690771056509416835771227975371492374573555172617972210304
        t0 *= t3; //  296: 27153497537739679690771056509416835771227975371492374573555172617972210325
        square_assign_multi(&mut t0, 5); //  301: 868911921207669750104673808301338744679295211887755986353765523775110730400
        t0 *= t2; //  302: 868911921207669750104673808301338744679295211887755986353765523775110730419
        square_assign_multi(&mut t0, 6); //  308: 55610362957290864006699123731285679659474893560816383126640993521607086746816
        t0 *= t1; //  309: 55610362957290864006699123731285679659474893560816383126640993521607086746829

        CtOption::new(t0, !self.ct_eq(&Self::zero()))
    }

    /// Constructs a `Scalar` element without checking that it is
    /// canonical.
    pub const fn from_raw_unchecked(v: [u64; 4]) -> Self {
        Scalar(v)
    }

    /// Outputs a `Scalar` element of multiplicative order equals to 2^n
    pub fn get_root_of_unity(n: u32) -> Self {
        assert!(n != 0, "cannot get root of unity for n = 0");
        assert!(n <= TWO_ADICITY, "order cannot exceed 2^{}", TWO_ADICITY);
        let power = 1u64 << (TWO_ADICITY - n);

        TWO_ADIC_ROOT_OF_UNITY.exp(&[power, 0, 0, 0])
    }
}

// OVERLOADED OPERATORS
// ================================================================================================

impl<'a> Neg for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn neg(self) -> Scalar {
        self.neg()
    }
}

impl Neg for Scalar {
    type Output = Scalar;

    #[inline]
    fn neg(self) -> Scalar {
        -&self
    }
}

impl<'a, 'b> Sub<&'b Scalar> for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn sub(self, rhs: &'b Scalar) -> Scalar {
        self.sub(rhs)
    }
}

impl<'a, 'b> Add<&'b Scalar> for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn add(self, rhs: &'b Scalar) -> Scalar {
        self.add(rhs)
    }
}

impl<T> Sum<T> for Scalar
where
    T: Borrow<Scalar>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Self::zero(), |acc, item| acc + item.borrow())
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn mul(self, rhs: &'b Scalar) -> Scalar {
        self.mul(rhs)
    }
}

impl_binops_additive!(Scalar, Scalar);
impl_binops_multiplicative!(Scalar, Scalar);

impl Div for Scalar {
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Scalar {
        self.mul(&rhs.invert().unwrap_or(Self::zero()))
    }
}

impl DivAssign for Scalar {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs
    }
}

// TYPE CONVERSIONS
// ================================================================================================

impl From<u128> for Scalar {
    /// Converts a 128-bit value into a field element. If the value is greater than or equal to
    /// the field modulus, modular reduction is silently performed.
    fn from(value: u128) -> Self {
        let value_high: u64 = (value >> 64).try_into().unwrap();
        let value_low: u64 = (value & (u64::MAX as u128)).try_into().unwrap();
        Scalar::new([value_low, value_high, 0, 0])
    }
}

impl From<u64> for Scalar {
    /// Converts a 64-bit value into a field element.
    fn from(value: u64) -> Self {
        Scalar([value, 0, 0, 0]) * R2
    }
}

impl From<u32> for Scalar {
    /// Converts a 32-bit value into a field element.
    fn from(value: u32) -> Self {
        Scalar([value as u64, 0, 0, 0]) * R2
    }
}

impl From<u16> for Scalar {
    /// Converts a 16-bit value into a field element.
    fn from(value: u16) -> Self {
        Scalar([value as u64, 0, 0, 0]) * R2
    }
}

impl From<u8> for Scalar {
    /// Converts an 8-bit value into a field element.
    fn from(value: u8) -> Self {
        Scalar([value as u64, 0, 0, 0]) * R2
    }
}

// FIELD TRAITS IMPLEMENTATION
// ================================================================================================

impl Field for Scalar {
    fn random(mut rng: impl RngCore) -> Self {
        Self::random(&mut rng)
    }

    fn zero() -> Self {
        Self::zero()
    }

    fn one() -> Self {
        Self::one()
    }

    fn is_zero(&self) -> Choice {
        self.ct_eq(&Self::zero())
    }

    #[must_use]
    fn square(&self) -> Self {
        self.square()
    }

    #[must_use]
    fn double(&self) -> Self {
        self.double()
    }

    fn invert(&self) -> CtOption<Self> {
        self.invert()
    }

    fn sqrt(&self) -> CtOption<Self> {
        self.sqrt()
    }
}

impl PrimeField for Scalar {
    type Repr = [u8; 32];

    fn from_repr(r: Self::Repr) -> CtOption<Self> {
        Self::from_bytes(&r)
    }

    fn to_repr(&self) -> Self::Repr {
        self.to_bytes()
    }

    fn is_odd(&self) -> Choice {
        (self.to_bytes()[0] & 1).ct_eq(&1)
    }

    const NUM_BITS: u32 = 255;
    const CAPACITY: u32 = Self::NUM_BITS - 1;

    fn multiplicative_generator() -> Self {
        GENERATOR
    }

    const S: u32 = TWO_ADICITY;

    fn root_of_unity() -> Self {
        TWO_ADIC_ROOT_OF_UNITY
    }
}

// SERDE SERIALIZATION
// ================================================================================================

#[cfg(feature = "serialize")]
impl Serialize for Scalar {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        use serde::ser::SerializeTuple;
        let mut tup = serializer.serialize_tuple(32)?;
        for byte in self.to_bytes().iter() {
            tup.serialize_element(byte)?;
        }
        tup.end()
    }
}

#[cfg(feature = "serialize")]
impl<'de> Deserialize<'de> for Scalar {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct ScalarVisitor;

        impl<'de> Visitor<'de> for ScalarVisitor {
            type Value = Scalar;

            fn expecting(&self, formatter: &mut Formatter) -> fmt::Result {
                formatter.write_str("a valid field element")
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<Scalar, A::Error>
            where
                A: serde::de::SeqAccess<'de>,
            {
                let mut bytes = [0u8; 32];
                for (i, byte) in bytes.iter_mut().enumerate() {
                    *byte = seq
                        .next_element()?
                        .ok_or_else(|| serde::de::Error::invalid_length(i, &"expected 32 bytes"))?;
                }
                let elem = Scalar::from_bytes(&bytes);
                if bool::from(elem.is_none()) {
                    Err(serde::de::Error::custom("decompression failed"))
                } else {
                    Ok(elem.unwrap())
                }
            }
        }

        deserializer.deserialize_tuple(32, ScalarVisitor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bitvec::view::AsBits;
    use rand_core::OsRng;

    const LARGEST: Scalar = Scalar([
        0xd443623eaed4acce,
        0x327aa72330157722,
        0x563fbf0f990a37b5,
        0x7af2599b3b3f22d0,
    ]);

    // DISPLAY
    // ================================================================================================

    #[test]
    fn test_debug() {
        assert_eq!(
            format!("{:?}", Scalar::zero()),
            "0x0000000000000000000000000000000000000000000000000000000000000000"
        );
        assert_eq!(
            format!("{:?}", Scalar::one()),
            "0x0000000000000000000000000000000000000000000000000000000000000001"
        );
        assert_eq!(
            format!("{:?}", R2),
            "0x0a1b4cc98981ba5f538081e0cdeb90959b0ab1b99fd511ba57793b82a256a662"
        );
    }

    #[test]
    fn test_output_reduced_limbs() {
        assert_eq!(
            format!("{:?}", Scalar::zero().output_reduced_limbs()),
            "[0, 0, 0, 0]"
        );
        assert_eq!(
            format!("{:?}", Scalar::one().output_reduced_limbs()),
            "[1, 0, 0, 0]"
        );
        assert_eq!(
            format!("{:?}", R2.output_reduced_limbs()),
            "[6303134585737094754, 11171937236454543802, 6016951904694407317, 728260193229584991]",
        );
    }

    // BASIC ALGEBRA
    // ================================================================================================

    #[test]
    fn test_equality() {
        assert_eq!(Scalar::default(), Scalar::zero());
        assert_eq!(Scalar::zero(), Scalar::zero());
        assert_eq!(Scalar::one(), Scalar::one());

        assert!(bool::from(Scalar::default().is_zero()));
        assert!(!bool::from(Scalar::zero().ct_eq(&Scalar::one())));

        assert!(Scalar::zero() != Scalar::one());
        assert!(Scalar::one() != R2);
    }

    #[test]
    fn test_addition() {
        let mut tmp = LARGEST;
        tmp += &LARGEST;

        assert_eq!(
            tmp,
            Scalar([
                0xd443623eaed4accd,
                0x327aa72330157722,
                0x563fbf0f990a37b5,
                0x7af2599b3b3f22d0,
            ])
        );

        assert_eq!(tmp, LARGEST.double());

        let mut tmp = LARGEST;
        tmp += &Scalar([1, 0, 0, 0]);

        assert_eq!(tmp, Scalar::zero());
    }

    #[test]
    fn test_subtraction() {
        let mut tmp = LARGEST;
        tmp -= &LARGEST;

        assert_eq!(tmp, Scalar::zero());

        let mut tmp = Scalar::zero();
        tmp -= &LARGEST;

        let mut tmp2 = M;
        tmp2 -= &LARGEST;

        assert_eq!(tmp, tmp2);
    }

    #[test]
    fn test_negation() {
        let mut rng = OsRng;

        let tmp = -&LARGEST;

        assert_eq!(tmp, Scalar([1, 0, 0, 0]));

        let tmp = -&Scalar::zero();
        assert_eq!(tmp, Scalar::zero());
        let tmp = -&Scalar([1, 0, 0, 0]);
        assert_eq!(tmp, LARGEST);

        for _ in 0..100 {
            let a = Scalar::random(&mut rng);
            let b = -a;

            assert_eq!(a + b, Scalar::zero());
        }
    }

    #[test]
    fn test_multiplication() {
        let mut cur = LARGEST;

        for _ in 0..100 {
            let mut tmp = cur;
            tmp *= &cur;

            let mut tmp2 = Scalar::zero();
            for b in cur
                .to_bytes()
                .iter()
                .rev()
                .flat_map(|byte| (0..8).rev().map(move |i| ((byte >> i) & 1u8) == 1u8))
            {
                let tmp3 = tmp2;
                tmp2.add_assign(&tmp3);

                if b {
                    tmp2.add_assign(&cur);
                }
            }

            assert_eq!(tmp, tmp2);

            cur.add_assign(&LARGEST);
        }
    }

    #[test]
    fn test_inversion() {
        assert!(bool::from(Scalar::zero().invert().is_none()));
        assert_eq!(Scalar::one().invert().unwrap(), Scalar::one());
        assert_eq!((-&Scalar::one()).invert().unwrap(), -&Scalar::one());

        let mut tmp = Scalar::random(&mut OsRng);

        for _ in 0..100 {
            let mut tmp2 = tmp.invert().unwrap();
            tmp2.mul_assign(&tmp);

            assert_eq!(tmp2, Scalar::one());

            tmp.add_assign(&Scalar::random(&mut OsRng));
        }
    }

    #[test]
    fn test_squaring() {
        let mut cur = LARGEST;

        for _ in 0..100 {
            let mut tmp = cur;
            let pow2 = tmp.exp(&[2, 0, 0, 0]);
            tmp = tmp.square();

            let mut tmp2 = Scalar::zero();
            for b in cur
                .to_bytes()
                .iter()
                .rev()
                .flat_map(|byte| (0..8).rev().map(move |i| ((byte >> i) & 1u8) == 1u8))
            {
                let tmp3 = tmp2;
                tmp2.add_assign(&tmp3);

                if b {
                    tmp2.add_assign(&cur);
                }
            }

            assert_eq!(tmp, tmp2);
            assert_eq!(tmp, pow2);

            cur.add_assign(&LARGEST);
        }
    }

    #[test]
    fn test_sqrt() {
        for _ in 0..100 {
            let a = Scalar::random(&mut OsRng).square();
            let b = a.sqrt().unwrap();
            assert_eq!(a, b.square());
        }

        assert_eq!(Scalar::zero().sqrt().unwrap(), Scalar::zero());
        assert_eq!(Scalar::one().sqrt().unwrap(), Scalar::one());

        // 3 is not a quadratic residue in Scalar
        assert!(bool::from(Scalar::new([3, 0, 0, 0]).sqrt().is_none()));
    }

    #[test]
    fn test_invert_is_pow() {
        let q_minus_2 = [
            0xd443623eaed4accd,
            0x327aa72330157722,
            0x563fbf0f990a37b5,
            0x7af2599b3b3f22d0,
        ];

        let mut r1 = R;
        let mut r2 = R;
        let mut r3 = R;

        for _ in 0..100 {
            r1 = r1.invert().unwrap();
            r2 = r2.exp(&q_minus_2);
            r3 = r3.exp_vartime(&q_minus_2);

            assert_eq!(r1, r2);
            assert_eq!(r2, r3);
            // Add R so we check something different next time around
            r1.add_assign(&R);
            r2 = r1;
            r3 = r1;
        }
    }

    // ROOTS OF UNITY
    // ================================================================================================

    #[test]
    fn test_get_root_of_unity() {
        let root_1 = Scalar::get_root_of_unity(1);
        assert_eq!(TWO_ADIC_ROOT_OF_UNITY, root_1);
        assert_eq!(Scalar::one(), root_1.exp(&[2, 0, 0, 0]));
    }

    #[test]
    fn test_lexicographically_largest() {
        // a = 18150892113463577006064251079316678276376639745798754653559980246979387702623
        let a = Scalar([
            0x2385774fb320cf85,
            0xd6831bd1db2dbef1,
            0xfca8a81ed3272e6d,
            0x4f36676f7b6b7531,
        ]);

        // b = 37459470843827287000634872651969001383098253815017628473081013274627699044208
        let b = Scalar([
            0xb0bdeaeefbb3dd4a,
            0x5bf78b5154e7b831,
            0x599716f0c5e30947,
            0x2bbbf22bbfd3ad9e,
        ]);

        assert_eq!(a.square(), b.square());
        assert!(!bool::from(a.lexicographically_largest()));
        assert!(bool::from(b.lexicographically_largest()));
    }

    #[test]
    fn test_zeroize() {
        use zeroize::Zeroize;

        let mut a = Scalar::one();
        a.zeroize();
        assert_eq!(a, Scalar::zero());
    }

    // INITIALIZATION
    // ================================================================================================

    #[test]
    fn test_from_int() {
        let n = 42u8;
        let element = Scalar::from(n);

        assert_eq!(element, Scalar::from(n as u16));
        assert_eq!(element, Scalar::from(n as u32));
        assert_eq!(element, Scalar::from(n as u64));
        assert_eq!(element, Scalar::from(n as u128));

        assert_eq!(element, Scalar::new([n as u64, 0, 0, 0]));
    }

    #[test]
    fn test_from_raw_unchecked() {
        let mut element = Scalar::from_raw_unchecked([
            0x57793b82a256a662,
            0x9b0ab1b99fd511ba,
            0x538081e0cdeb9095,
            0x0a1b4cc98981ba5f,
        ]);

        let element_normalized = Scalar::new([
            0x57793b82a256a662,
            0x9b0ab1b99fd511ba,
            0x538081e0cdeb9095,
            0x0a1b4cc98981ba5f,
        ]);

        assert_eq!(element, Scalar::one());
        element *= &R2;

        assert!(element != Scalar::one());
        assert_eq!(element, element_normalized);
    }

    // FIELD TRAIT
    // ================================================================================================

    #[test]
    fn test_field_trait_methods() {
        assert_eq!(<Scalar as Field>::zero(), Scalar::new([0, 0, 0, 0]));
        assert_eq!(<Scalar as Field>::one(), Scalar::new([1, 0, 0, 0]));

        assert_eq!(
            bool::from(<Scalar as Field>::zero().is_zero()),
            bool::from(Scalar::new([0, 0, 0, 0]).is_zero())
        );
        assert_eq!(
            bool::from(<Scalar as Field>::one().is_zero()),
            bool::from(Scalar::new([1, 0, 0, 0]).is_zero())
        );

        let mut rng = OsRng;
        let e = Scalar::random(&mut rng).square();

        assert_eq!(<Scalar as Field>::square(&e), e.square());
        assert_eq!(<Scalar as Field>::double(&e), e.double());

        assert_eq!(<Scalar as Field>::invert(&e).unwrap(), e.invert().unwrap());
        assert!(bool::from(
            <Scalar as Field>::invert(&Scalar::zero()).is_none()
        ));

        assert_eq!(<Scalar as Field>::sqrt(&e).unwrap(), e.sqrt().unwrap());
        assert!(bool::from(
            <Scalar as Field>::sqrt(&Scalar::new([3, 0, 0, 0])).is_none()
        ));
    }

    #[test]
    fn test_primefield_trait_methods() {
        assert_eq!(
            Scalar::from_repr([
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ])
            .unwrap(),
            Scalar::zero()
        );
        assert_eq!(
            Scalar::from_repr([
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ])
            .unwrap(),
            Scalar::one()
        );

        let mut rng = OsRng;
        let e = Scalar::random(&mut rng).square();

        assert_eq!(Scalar::from_repr(Scalar::to_repr(&e)).unwrap(), e);

        assert_eq!(
            Scalar::get_root_of_unity(<Scalar as PrimeField>::S),
            <Scalar as PrimeField>::root_of_unity()
        )
    }

    // SERIALIZATION / DESERIALIZATION
    // ================================================================================================

    #[test]
    fn test_to_radix16() {
        let mut rng = OsRng;

        for _ in 0..100 {
            let a = Scalar::random(&mut rng);
            let digits = Scalar::bytes_to_radix_16(&a.to_bytes());

            let radix = Scalar::from(16u64);
            let mut term = Scalar::one();
            let mut a_bis = Scalar::zero();
            for &digit in digits.iter() {
                if digit < 0 {
                    a_bis += -Scalar::from((-(digit as i64)) as u64) * term;
                } else {
                    a_bis += Scalar::from(digit as u64) * term;
                };
                term *= radix;
            }

            assert_eq!(a_bis, a);
        }
    }

    #[test]
    fn test_to_naf() {
        let mut rng = OsRng;
        for _ in 0..100 {
            let a = Scalar::random(&mut rng);

            for w in [2, 3, 4, 5, 6, 7, 8] {
                let digits = Scalar::bytes_to_wnaf_vartime(&a.to_bytes(), w);

                let mut b = Scalar::zero();
                for &digit in digits.iter().rev() {
                    if digit < 0 {
                        b = b.double() - Scalar::from((-digit as i64) as u64);
                    } else {
                        b = b.double() + Scalar::from(digit as u64);
                    };
                }

                assert_eq!(a, b);
            }
        }
    }

    #[test]
    fn test_to_bytes() {
        assert_eq!(
            Scalar::zero().to_bytes(),
            [
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ]
        );

        assert_eq!(
            Scalar::one().to_bytes(),
            [
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ]
        );

        assert_eq!(
            R2.to_bytes(),
            [
                98, 166, 86, 162, 130, 59, 121, 87, 186, 17, 213, 159, 185, 177, 10, 155, 149, 144,
                235, 205, 224, 129, 128, 83, 95, 186, 129, 137, 201, 76, 27, 10
            ]
        );

        assert_eq!(
            (-&Scalar::one()).to_bytes(),
            [
                206, 172, 212, 174, 62, 98, 67, 212, 34, 119, 21, 48, 35, 167, 122, 50, 181, 55,
                10, 153, 15, 191, 63, 86, 208, 34, 63, 59, 155, 89, 242, 122
            ]
        );
    }

    #[test]
    fn test_from_bytes_non_canonical() {
        let mut rng = OsRng;
        for _ in 0..100 {
            let a = Scalar::random(&mut rng);
            let bytes = a.to_bytes();
            assert_eq!(a, Scalar::from_bytes_non_canonical(&bytes));
        }

        assert_eq!(
            Scalar::from_bytes_non_canonical(&[
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ]),
            Scalar::zero()
        );

        assert_eq!(
            Scalar::from_bytes_non_canonical(&[
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ]),
            Scalar::one()
        );

        assert_eq!(
            Scalar::from_bytes_non_canonical(&[
                98, 166, 86, 162, 130, 59, 121, 87, 186, 17, 213, 159, 185, 177, 10, 155, 149, 144,
                235, 205, 224, 129, 128, 83, 95, 186, 129, 137, 201, 76, 27, 10
            ]),
            R2
        );

        assert_eq!(
            Scalar::from_bytes_non_canonical(&[
                206, 172, 212, 174, 62, 98, 67, 212, 34, 119, 21, 48, 35, 167, 122, 50, 181, 55,
                10, 153, 15, 191, 63, 86, 208, 34, 63, 59, 155, 89, 242, 122
            ]),
            -Scalar::one(),
        );

        // M will yield zero
        assert_eq!(
            Scalar::from_bytes_non_canonical(&[
                207, 172, 212, 174, 62, 98, 67, 212, 34, 119, 21, 48, 35, 167, 122, 50, 181, 55,
                10, 153, 15, 191, 63, 86, 208, 34, 63, 59, 155, 89, 242, 122
            ]),
            Scalar::zero(),
        );

        // M+1 will yield one
        assert_eq!(
            Scalar::from_bytes_non_canonical(&[
                208, 172, 212, 174, 62, 98, 67, 212, 34, 119, 21, 48, 35, 167, 122, 50, 181, 55,
                10, 153, 15, 191, 63, 86, 208, 34, 63, 59, 155, 89, 242, 122
            ]),
            Scalar::one(),
        );
    }

    #[test]
    fn test_from_bytes() {
        let mut rng = OsRng;
        for _ in 0..100 {
            let a = Scalar::random(&mut rng);
            let bytes = a.to_bytes();
            assert_eq!(a, Scalar::from_bytes(&bytes).unwrap());
        }

        assert_eq!(
            Scalar::from_bytes(&[
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ])
            .unwrap(),
            Scalar::zero()
        );

        assert_eq!(
            Scalar::from_bytes(&[
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ])
            .unwrap(),
            Scalar::one()
        );

        assert_eq!(
            Scalar::from_bytes(&[
                98, 166, 86, 162, 130, 59, 121, 87, 186, 17, 213, 159, 185, 177, 10, 155, 149, 144,
                235, 205, 224, 129, 128, 83, 95, 186, 129, 137, 201, 76, 27, 10
            ])
            .unwrap(),
            R2
        );

        // -1 should work
        assert_eq!(
            Scalar::from_bytes(&[
                206, 172, 212, 174, 62, 98, 67, 212, 34, 119, 21, 48, 35, 167, 122, 50, 181, 55,
                10, 153, 15, 191, 63, 86, 208, 34, 63, 59, 155, 89, 242, 122
            ])
            .unwrap(),
            -Scalar::one(),
        );

        // M is invalid
        assert!(bool::from(
            Scalar::from_bytes(&[
                207, 172, 212, 174, 62, 98, 67, 212, 34, 119, 21, 48, 35, 167, 122, 50, 181, 55,
                10, 153, 15, 191, 63, 86, 208, 34, 63, 59, 155, 89, 242, 122
            ])
            .is_none()
        ));

        // Anything larger than M is invalid
        assert!(bool::from(
            Scalar::from_bytes(&[
                208, 172, 212, 174, 62, 98, 67, 212, 34, 119, 21, 48, 35, 167, 122, 50, 181, 55,
                10, 153, 15, 191, 63, 86, 208, 34, 63, 59, 155, 89, 242, 122
            ])
            .is_none()
        ));
        assert!(bool::from(
            Scalar::from_bytes(&[
                206, 173, 212, 174, 62, 98, 67, 212, 34, 119, 21, 48, 35, 167, 122, 50, 181, 55,
                10, 153, 15, 191, 63, 86, 208, 34, 63, 59, 155, 89, 242, 122
            ])
            .is_none()
        ));
        assert!(bool::from(
            Scalar::from_bytes(&[
                0, 0, 0, 174, 62, 98, 67, 212, 34, 119, 21, 48, 35, 167, 122, 50, 181, 55, 10, 153,
                15, 191, 63, 86, 208, 34, 63, 59, 155, 89, 242, 255
            ])
            .is_none()
        ));
    }

    #[test]
    fn test_from_u512_zero() {
        assert_eq!(
            Scalar::zero(),
            Scalar::from_u512([M.0[0], M.0[1], M.0[2], M.0[3], 0, 0, 0, 0])
        );
    }

    #[test]
    fn test_from_u512_r() {
        assert_eq!(R, Scalar::from_u512([1, 0, 0, 0, 0, 0, 0, 0]));
    }

    #[test]
    fn test_from_u512_r2() {
        assert_eq!(R2, Scalar::from_u512([0, 0, 0, 0, 1, 0, 0, 0]));
    }

    #[test]
    fn test_from_u512_max() {
        let max_u64 = 0xffff_ffff_ffff_ffff;
        assert_eq!(
            R3 - R,
            Scalar::from_u512([
                max_u64, max_u64, max_u64, max_u64, max_u64, max_u64, max_u64, max_u64
            ])
        );
    }

    #[test]
    fn test_from_bytes_wide_r2() {
        assert_eq!(
            R2,
            Scalar::from_bytes_wide(&[
                98, 166, 86, 162, 130, 59, 121, 87, 186, 17, 213, 159, 185, 177, 10, 155, 149, 144,
                235, 205, 224, 129, 128, 83, 95, 186, 129, 137, 201, 76, 27, 10, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
        );
    }

    #[test]
    fn test_from_bytes_wide_negative_one() {
        assert_eq!(
            -&Scalar::one(),
            Scalar::from_bytes_wide(&[
                206, 172, 212, 174, 62, 98, 67, 212, 34, 119, 21, 48, 35, 167, 122, 50, 181, 55,
                10, 153, 15, 191, 63, 86, 208, 34, 63, 59, 155, 89, 242, 122, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
        );
    }

    #[test]
    fn test_from_bytes_wide_maximum() {
        assert_eq!(
            Scalar([
                0xd607a1e33be17918,
                0x9f25fbf5db9be97c,
                0xfdd3f73cf4b234d2,
                0x67fe48cbe00f6904,
            ]),
            Scalar::from_bytes_wide(&[0xff; 64])
        );
    }

    #[test]
    fn test_from_bits() {
        let bytes = Scalar::zero().to_bytes();
        assert_eq!(Scalar::from_bits(bytes.as_bits::<Lsb0>()), Scalar::zero());
        assert_eq!(
            Scalar::from_bits_vartime(bytes.as_bits::<Lsb0>()),
            Scalar::zero()
        );

        let bytes = Scalar::one().to_bytes();
        assert_eq!(Scalar::from_bits(bytes.as_bits::<Lsb0>()), Scalar::one());
        assert_eq!(
            Scalar::from_bits_vartime(bytes.as_bits::<Lsb0>()),
            Scalar::one()
        );

        let bytes = R2.to_bytes();
        assert_eq!(Scalar::from_bits(bytes.as_bits::<Lsb0>()), R2);
        assert_eq!(Scalar::from_bits_vartime(bytes.as_bits::<Lsb0>()), R2);

        // -1 should work
        let bytes = (-Scalar::one()).to_bytes();
        assert_eq!(Scalar::from_bits(bytes.as_bits::<Lsb0>()), -Scalar::one());
        assert_eq!(
            Scalar::from_bits_vartime(bytes.as_bits::<Lsb0>()),
            -Scalar::one()
        );

        // Modulus results in Scalar::zero()
        let bytes = [
            207, 172, 212, 174, 62, 98, 67, 212, 34, 119, 21, 48, 35, 167, 122, 50, 181, 55, 10,
            153, 15, 191, 63, 86, 208, 34, 63, 59, 155, 89, 242, 122,
        ];
        assert_eq!(Scalar::from_bits(bytes.as_bits::<Lsb0>()), Scalar::zero());
        assert_eq!(
            Scalar::from_bits_vartime(bytes.as_bits::<Lsb0>()),
            Scalar::zero()
        );
    }

    #[test]
    #[cfg(feature = "serialize")]
    fn test_serde() {
        let mut rng = OsRng;
        let element = Scalar::random(&mut rng);
        let encoded = bincode::serialize(&element).unwrap();
        let parsed: Scalar = bincode::deserialize(&encoded).unwrap();
        assert_eq!(parsed, element);

        // Check that the encoding is 32 bytes exactly
        assert_eq!(encoded.len(), 32);

        // Check that the encoding itself matches the usual one
        assert_eq!(element, bincode::deserialize(&element.to_bytes()).unwrap());

        // Check that invalid encodings fail
        let element = Scalar::random(&mut rng);
        let mut encoded = bincode::serialize(&element).unwrap();
        encoded[31] = 127;
        assert!(bincode::deserialize::<Scalar>(&encoded).is_err());

        let encoded = bincode::serialize(&element).unwrap();
        assert!(bincode::deserialize::<Scalar>(&encoded[0..31]).is_err());
    }
}
