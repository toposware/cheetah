// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This module provides an implementation of the scalar field Fq
//! of the STARK-friendly cheetah curve, with characteristic
//! q = 0x26337f752795f77cb6b6ebb9a18fecc9f2f264f035242b271e13aee130956aa5.

use core::{
    convert::{TryFrom, TryInto},
    fmt::{self, Debug, Display, Formatter},
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

// Field modulus = 17278877126736494933592566161653303514319447234579276854188469089485337225893
const M: Scalar = Scalar([
    0x1e13aee130956aa5,
    0xf2f264f035242b27,
    0xb6b6ebb9a18fecc9,
    0x26337f752795f77c,
]);

// 2^256 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R: Scalar = Scalar([
    0x4b89e6b8dc7f8022,
    0x4e51a25ec126fd15,
    0xb7b679a636a07344,
    0x1acb0341127c3313,
]);

// 2^512 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R2: Scalar = Scalar([
    0xd081a0324c5e8462,
    0x800b297bdb4d56da,
    0x9d4305b1df84c934,
    0x1216a239925b192a,
]);

// 2^768 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R3: Scalar = Scalar([
    0xbbca5047a6bb9424,
    0x2f7707c1f6730698,
    0x1dfca6ccdfbca422,
    0x20d20982d647c982,
]);

// Multiplicative generator g of order q-1
// g = 2
//   = 0xf62870cfd626eaab8b60792cbb0f9bea9b0dfcd4d29cf0379001e908869959f in Montgomery form
const GENERATOR: Scalar = Scalar([
    0x79001e908869959f,
    0xa9b0dfcd4d29cf03,
    0xb8b60792cbb0f9be,
    0x0f62870cfd626eaa,
]);

// Two-adicity of the field: (q-1) % 2^2 = 0
const TWO_ADICITY: u32 = 2;

// 2^2 root of unity = 0x2d13d2b2496ec9fc67e21dae0a414bb8abd9c41c985b3108f825c4185efb292
//                   = 0x1e31bde7bb98bb8d8df21c9eb075903b797d57622ea1c8fd2f0afc7b183ab212 in Montgomery form
const TWO_ADIC_ROOT_OF_UNITY: Scalar = Scalar([
    0x2f0afc7b183ab212,
    0x797d57622ea1c8fd,
    0x8df21c9eb075903b,
    0x1e31bde7bb98bb8d,
]);

// T = 2^[(M-5)/8]
//   = 0x11b12125017f856e781c64ef6075ec07341a645735cf3c0b4748a94fd552dc09
//   = 0x9b51ee0c08b802c13e2a0972e04eb0a0f0ae80fbd3fc81de0c93d473638518b in Montgomery form
const T: Scalar = Scalar([
    0xe0c93d473638518b,
    0x0f0ae80fbd3fc81d,
    0x13e2a0972e04eb0a,
    0x09b51ee0c08b802c,
]);

// -M^{-1} mod 2^64; this is used during element multiplication.
const U: u64 = 11594536110104744659;

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
        // Atkin's algorithm for q mod 8 = 5
        // See https://eprint.iacr.org/2012/685.pdf, algorithm 3

        // a1 = self^((q - 5) // 8)
        //   = self^0x4c66feea4f2beef96d6dd773431fd993e5e4c9e06a48564e3c275dc2612ad54
        let a1 = self.exp_vartime(&[
            0xe3c275dc2612ad54,
            0x3e5e4c9e06a48564,
            0x96d6dd773431fd99,
            0x04c66feea4f2beef,
        ]);

        let a0 = (a1.square() * self).square();
        let a0_is_not_minus_one = !a0.ct_eq(&Scalar::one().neg());

        let b = T * a1;
        let ab = self * b;
        let i = ab.double() * b - Scalar::one();
        let x = ab * i;

        CtOption::new(x, (x * x).ct_eq(self) & a0_is_not_minus_one)
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
    /// representation, where each resulting coefficient is in [-8;8).
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

    /// Returns whether or not this element is strictly lexicographically
    /// larger than its negation.
    pub fn lexicographically_largest(&self) -> Choice {
        // This can be determined by checking to see if the element is
        // larger than (M - 1) // 2. If we subtract by ((M - 1) // 2) + 1
        // and there is no underflow, then the element must be larger than
        // (M - 1) // 2.

        // First, because self is in Montgomery form we need to reduce it
        let tmp = Scalar::montgomery_reduce(self.0[0], self.0[1], self.0[2], self.0[3], 0, 0, 0, 0);

        let (_, borrow) = sub64_with_carry(tmp.0[0], 0x8f09d770984ab552, 0);
        let (_, borrow) = sub64_with_carry(tmp.0[1], 0xf97932781a921593, borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[2], 0x5b5b75dcd0c7f664, borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[3], 0x1319bfba93cafbbe, borrow);

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
        let t3 = self.square(); //             1: 2
        let t1 = t3 * self; //                 2: 3
        let mut t0 = t3.square(); //           3: 4
        let t2 = t1 * t3; //                   4: 5
        let t6 = t0 * t1; //                   5: 7
        let t4 = t2 * t0; //                   6: 9
        let t15 = t6 * t0; //                  7: 11
        let t9 = t4 * t0; //                   8: 13
        let t7 = t15 * t0; //                  9: 15
        let t5 = t7 * t0; //                  10: 19
        let t3 = t5 * t3; //                  11: 21
        let t11 = t5 * t0; //                 12: 23
        let t10 = t3 * t0; //                 13: 25
        let t13 = t11 * t0; //                14: 27
        let t14 = t10 * t0; //                15: 29
        let t12 = t13 * t0; //                16: 31
        t0 = t5.square(); //                  17: 38
        square_assign_multi(&mut t0, 7); //   24: 4864
        t0 *= t10; //                         25: 4889
        square_assign_multi(&mut t0, 5); //   30: 156448
        t0 *= t11; //                         31: 156471
        square_assign_multi(&mut t0, 4); //   35: 2503536
        t0 *= t7; //                          36: 2503551
        square_assign_multi(&mut t0, 6); //   42: 160227264
        t0 *= t14; //                         43: 160227293
        square_assign_multi(&mut t0, 5); //   48: 5127273376
        t0 *= t4; //                          49: 5127273385
        square_assign_multi(&mut t0, 6); //   55: 328145496640
        t0 *= t7; //                          56: 328145496655
        square_assign_multi(&mut t0, 7); //   63: 42002623571840
        t0 *= t3; //                          64: 42002623571861
        square_assign_multi(&mut t0, 4); //   68: 672041977149776
        t0 *= t7; //                          69: 672041977149791
        square_assign_multi(&mut t0, 6); //   75: 43010686537586624
        t0 *= t14; //                         76: 43010686537586653
        square_assign_multi(&mut t0, 4); //   80: 688170984601386448
        t0 *= t7; //                          81: 688170984601386463
        square_assign_multi(&mut t0, 6); //   87: 44042943014488733632
        t0 *= t15; //                         88: 44042943014488733643
        square_assign_multi(&mut t0, 5); //   93: 1409374176463639476576
        t0 *= t9; //                          94: 1409374176463639476589
        square_assign_multi(&mut t0, 6); //  100: 90199947293672926501696
        t0 *= t13; //                        101: 90199947293672926501723
        square_assign_multi(&mut t0, 6); //  107: 5772796626795067296110272
        t0 *= t14; //                        108: 5772796626795067296110301
        square_assign_multi(&mut t0, 6); //  114: 369458984114884306951059264
        t0 *= t14; //                        115: 369458984114884306951059293
        square_assign_multi(&mut t0, 2); //  117: 1477835936459537227804237172
        t0 *= t1; //                         118: 1477835936459537227804237175
        square_assign_multi(&mut t0, 6); //  124: 94581499933410382579471179200
        t0 *= t9; //                         125: 94581499933410382579471179213
        square_assign_multi(&mut t0, 6); //  131: 6053215995738264485086155469632
        t0 *= t1; //                         132: 6053215995738264485086155469635
        square_assign_multi(&mut t0, 8); //  140: 1549623294908995708182055800226560
        t0 *= t12; //                        141: 1549623294908995708182055800226591
        square_assign_multi(&mut t0, 5); //  146: 49587945437087862661825785607250912
        t0 *= t13; //                        147: 49587945437087862661825785607250939
        square_assign_multi(&mut t0, 7); //  154: 6347257015947246420713700557728120192
        t0 *= t10; //                        155: 6347257015947246420713700557728120217
        square_assign_multi(&mut t0, 7); //  162: 812448898041247541851353671389199387776
        t0 *= t12; //                        163: 812448898041247541851353671389199387807
        square_assign_multi(&mut t0, 7); //  170: 103993458949279685356973269937817521639296
        t0 *= t11; //                        171: 103993458949279685356973269937817521639319
        square_assign_multi(&mut t0, 4); //  175: 1663895343188474965711572319005080346229104
        t0 *= t4; //                         176: 1663895343188474965711572319005080346229113
        square_assign_multi(&mut t0, 7); //  183: 212978603928124795611081256832650284317326464
        t0 *= t10; //                        184: 212978603928124795611081256832650284317326489
        square_assign_multi(&mut t0, 6); //  190: 13630630651399986919109200437289618196308895296
        t0 *= t7; //                         191: 13630630651399986919109200437289618196308895311
        square_assign_multi(&mut t0, 10); // 201: 13957765787033586605167821247784569033020308798464
        t0 *= t9; //                         202: 13957765787033586605167821247784569033020308798477
        square_assign_multi(&mut t0, 5); //  207: 446648505185074771365370279929106209056649881551264
        t0 *= t4; //                         208: 446648505185074771365370279929106209056649881551273
        square_assign_multi(&mut t0, 3); //  211: 3573188041480598170922962239432849672453199052410184
        t0 *= self; //                       212: 3573188041480598170922962239432849672453199052410185
        square_assign_multi(&mut t0, 9); //  221: 1829472277238066263512556666589619032296037914834014720
        t0 *= t3; //                         222: 1829472277238066263512556666589619032296037914834014741
        square_assign_multi(&mut t0, 4); //  226: 29271556435809060216200906665433904516736606637344235856
        t0 *= t4; //                         227: 29271556435809060216200906665433904516736606637344235865
        square_assign_multi(&mut t0, 5); //  232: 936689805945889926918429013293884944535571412395015547680
        t0 *= t6; //                         233: 936689805945889926918429013293884944535571412395015547687
        square_assign_multi(&mut t0, 7); //  240: 119896295161073910645558913701617272900553140786561990103936
        t0 *= t7; //                         241: 119896295161073910645558913701617272900553140786561990103951
        square_assign_multi(&mut t0, 9); //  250: 61386903122469842250526163815228043725083208082719738933222912
        t0 *= t5; //                         251: 61386903122469842250526163815228043725083208082719738933222931
        square_assign_multi(&mut t0, 5); //  256: 1964380899919034952016837242087297399202662658647031645863133792
        t0 *= t3; //                         257: 1964380899919034952016837242087297399202662658647031645863133813
        square_assign_multi(&mut t0, 2); //  259: 7857523599676139808067348968349189596810650634588126583452535252
        t0 *= t1; //                         260: 7857523599676139808067348968349189596810650634588126583452535255
        square_assign_multi(&mut t0, 4); //  264: 125720377594818236929077583493587033548970410153410025335240564080
        t0 *= t6; //                         265: 125720377594818236929077583493587033548970410153410025335240564087
        square_assign_multi(&mut t0, 9); //  274: 64368833328546937307687722748716561177072849998545932971643168812544
        t0 *= t5; //                         275: 64368833328546937307687722748716561177072849998545932971643168812563
        square_assign_multi(&mut t0, 8); //  283: 16478421332108015950768057023671439661330649599627758840740651216016128
        t0 *= t4; //                         284: 16478421332108015950768057023671439661330649599627758840740651216016137
        square_assign_multi(&mut t0, 6); //  290: 1054618965254913020849155649514972138325161574376176565807401677825032768
        t0 *= t3; //                         291: 1054618965254913020849155649514972138325161574376176565807401677825032789
        square_assign_multi(&mut t0, 5); //  296: 33747806888157216667172980784479108426405170380037650105836853690401049248
        t0 *= t3; //                         297: 33747806888157216667172980784479108426405170380037650105836853690401049269
        square_assign_multi(&mut t0, 4); //  301: 539964910210515466674767692551665734822482726080602401693389659046416788304
        t0 *= t2; //                         302: 539964910210515466674767692551665734822482726080602401693389659046416788309
        square_assign_multi(&mut t0, 5); //  307: 17278877126736494933592566161653303514319447234579276854188469089485337225888
        t0 *= t1; //                         308: 17278877126736494933592566161653303514319447234579276854188469089485337225891 = M - 2

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
        let mut buf = [0; 64];
        rng.fill_bytes(&mut buf);
        Self::from_bytes_wide(&buf)
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

    const NUM_BITS: u32 = 254;
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
        0x1e13aee130956aa4,
        0xf2f264f035242b27,
        0xb6b6ebb9a18fecc9,
        0x26337f752795f77c,
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
            "0x1acb0341127c3313b7b679a636a073444e51a25ec126fd154b89e6b8dc7f8022"
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
            "[5443135306301669410, 5643470335923125525, 13237901909490168644, 1930640443276276499]"
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
                0x1e13aee130956aa3,
                0xf2f264f035242b27,
                0xb6b6ebb9a18fecc9,
                0x26337f752795f77c,
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
            0x1e13aee130956aa3,
            0xf2f264f035242b27,
            0xb6b6ebb9a18fecc9,
            0x26337f752795f77c,
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
        let root_2 = Scalar::get_root_of_unity(2);
        assert_eq!(TWO_ADIC_ROOT_OF_UNITY, root_2);
        assert_eq!(Scalar::one(), root_2.exp(&[4, 0, 0, 0]));

        let root_1 = Scalar::get_root_of_unity(1);
        let expected = root_2.exp(&[2, 0, 0, 0]);
        assert_eq!(expected, root_1);
        assert_eq!(Scalar::one(), root_1.exp(&[2, 0, 0, 0]));
    }

    #[test]
    fn test_lexicographically_largest() {
        // a = 1475249844745184162945519477422622175802656699045099721608794143422599449058
        let a = Scalar([
            0x42e5c81ee0c2f208,
            0x7357405276648944,
            0xe256db77e0fb5f75,
            0x1ef35384edf87a42,
        ]);

        // b = 15803627281991310770647046684230681338516790535534177132579674946062737776835
        let b = Scalar([
            0xdb2de6c24fd2789d,
            0x7f9b249dbebfa1e2,
            0xd4601041c0948d54,
            0x07402bf0399d7d39,
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
            5443135306301669410,
            5643470335923125525,
            13237901909490168644,
            1930640443276276499,
        ]);

        let element_normalized = Scalar::new([
            5443135306301669410,
            5643470335923125525,
            13237901909490168644,
            1930640443276276499,
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
                34, 128, 127, 220, 184, 230, 137, 75, 21, 253, 38, 193, 94, 162, 81, 78, 68, 115,
                160, 54, 166, 121, 182, 183, 19, 51, 124, 18, 65, 3, 203, 26
            ]
        );

        assert_eq!(
            (-&Scalar::one()).to_bytes(),
            [
                164, 106, 149, 48, 225, 174, 19, 30, 39, 43, 36, 53, 240, 100, 242, 242, 201, 236,
                143, 161, 185, 235, 182, 182, 124, 247, 149, 39, 117, 127, 51, 38
            ]
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
                34, 128, 127, 220, 184, 230, 137, 75, 21, 253, 38, 193, 94, 162, 81, 78, 68, 115,
                160, 54, 166, 121, 182, 183, 19, 51, 124, 18, 65, 3, 203, 26
            ])
            .unwrap(),
            R2
        );

        // -1 should work
        assert_eq!(
            Scalar::from_bytes(&[
                164, 106, 149, 48, 225, 174, 19, 30, 39, 43, 36, 53, 240, 100, 242, 242, 201, 236,
                143, 161, 185, 235, 182, 182, 124, 247, 149, 39, 117, 127, 51, 38
            ])
            .unwrap(),
            -Scalar::one(),
        );

        // M is invalid
        assert!(bool::from(
            Scalar::from_bytes(&[
                165, 106, 149, 48, 225, 174, 19, 30, 39, 43, 36, 53, 240, 100, 242, 242, 201, 236,
                143, 161, 185, 235, 182, 182, 124, 247, 149, 39, 117, 127, 51, 38
            ])
            .is_none()
        ));

        // Anything larger than M is invalid
        assert!(bool::from(
            Scalar::from_bytes(&[
                166, 106, 149, 48, 225, 174, 19, 30, 39, 43, 36, 53, 240, 100, 242, 242, 201, 236,
                143, 161, 185, 235, 182, 182, 124, 247, 149, 39, 117, 127, 51, 38
            ])
            .is_none()
        ));
        assert!(bool::from(
            Scalar::from_bytes(&[
                164, 255, 149, 48, 225, 174, 19, 30, 39, 43, 36, 53, 240, 100, 242, 242, 201, 236,
                143, 161, 185, 235, 182, 182, 124, 247, 149, 39, 117, 127, 51, 38
            ])
            .is_none()
        ));
        assert!(bool::from(
            Scalar::from_bytes(&[
                0, 0, 0, 48, 225, 174, 19, 30, 39, 43, 36, 53, 240, 100, 242, 242, 201, 236, 143,
                161, 185, 235, 182, 182, 124, 247, 149, 39, 117, 127, 51, 39
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
                34, 128, 127, 220, 184, 230, 137, 75, 21, 253, 38, 193, 94, 162, 81, 78, 68, 115,
                160, 54, 166, 121, 182, 183, 19, 51, 124, 18, 65, 3, 203, 26, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
        );
    }

    #[test]
    fn test_from_bytes_wide_negative_one() {
        assert_eq!(
            -&Scalar::one(),
            Scalar::from_bytes_wide(&[
                164, 106, 149, 48, 225, 174, 19, 30, 39, 43, 36, 53, 240, 100, 242, 242, 201, 236,
                143, 161, 185, 235, 182, 182, 124, 247, 149, 39, 117, 127, 51, 38, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
        );
    }

    #[test]
    fn test_from_bytes_wide_maximum() {
        assert_eq!(
            Scalar([
                0x7040698eca3c1402,
                0xe1256563354c0983,
                0x66462d26a91c30dd,
                0x06070641c3cb966e,
            ]),
            Scalar::from_bytes_wide(&[0xff; 64])
        );
    }

    #[test]
    fn test_from_bits() {
        let bytes = Scalar::zero().to_bytes();
        assert_eq!(Scalar::from_bits(bytes.as_bits::<Lsb0>()), Scalar::zero());

        let bytes = Scalar::one().to_bytes();
        assert_eq!(Scalar::from_bits(bytes.as_bits::<Lsb0>()), Scalar::one());

        let bytes = R2.to_bytes();
        assert_eq!(Scalar::from_bits(bytes.as_bits::<Lsb0>()), R2);

        // -1 should work
        let bytes = (-Scalar::one()).to_bytes();
        assert_eq!(Scalar::from_bits(bytes.as_bits::<Lsb0>()), -Scalar::one());

        // Modulus results in Scalar::zero()
        let bytes = [
            165, 106, 149, 48, 225, 174, 19, 30, 39, 43, 36, 53, 240, 100, 242, 242, 201, 236, 143,
            161, 185, 235, 182, 182, 124, 247, 149, 39, 117, 127, 51, 38,
        ];
        assert_eq!(Scalar::from_bits(bytes.as_bits::<Lsb0>()), Scalar::zero());
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
