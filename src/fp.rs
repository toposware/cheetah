// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This module provides an implementation of the finite prime
//! field Fp of characteristic p = 2^62 + 2^56 + 2^55 + 1.

use core::{
    convert::TryFrom,
    fmt::{self, Debug, Display, Formatter},
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::utils::{add64_with_carry, mul64_with_carry, square_assign_multi, sub64_with_carry};

use group::ff::{Field, PrimeField};
use rand_core::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "serialize")]
use serde::de::Visitor;
#[cfg(feature = "serialize")]
use serde::{self, Deserialize, Deserializer, Serialize, Serializer};

// CONSTANTS
// ================================================================================================

// ******************************** //
// ********* FP CONSTANTS ********* //
// ******************************** //

// Field modulus = 2^62 + 2^56 + 2^55 + 1
const M: Fp = Fp(4719772409484279809);

// 2^64 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R: Fp = Fp(4287426845256712189);

// 2^128 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R2: Fp = Fp(3635333122111952146);

// 2^192 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R3: Fp = Fp(4670764763728573444);

// Multiplicative generator g of order p-1
// g = 3
//   = 3422735716801576949 in Montgomery form
const GENERATOR: Fp = Fp(3422735716801576949);

// Two-adicity of the field: (p-1) % 2^55 = 0
pub(crate) const TWO_ADICITY: u32 = 55;

// 2^55 root of unity = 90479342105353296
//                    = 1519868260574363836 in Montgomery form
const TWO_ADIC_ROOT_OF_UNITY: Fp = Fp(0x1517a82160ed00bc);

// -M^{-1} mod 2^64; this is used during element multiplication.
const U: u64 = 4719772409484279807;

// FIELD ELEMENT
// ================================================================================================

/// Represents a base field element.
///
/// Internal values are stored in Montgomery representation.
/// The backing type is `u64`.
#[derive(Copy, Clone, Eq, Default)]
pub struct Fp(pub(crate) u64);

impl Debug for Fp {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        let tmp = self.output_reduced_limbs();
        write!(f, "{:?}", tmp)
    }
}

impl Display for Fp {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl ConstantTimeEq for Fp {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0.ct_eq(&other.0)
    }
}

impl PartialEq for Fp {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl ConditionallySelectable for Fp {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Fp(u64::conditional_select(&a.0, &b.0, choice))
    }
}

impl zeroize::DefaultIsZeroes for Fp {}

impl Fp {
    /// Creates a new field element from a `u64` value.
    /// The value is converted to Montgomery form by computing
    /// (a.R^0 * R^2) / R = a.R
    pub const fn new(value: u64) -> Self {
        (&Fp(value)).mul(&R2)
    }

    /// Returns zero, the additive identity.
    #[inline]
    pub const fn zero() -> Self {
        Fp(0)
    }

    /// Returns one, the multiplicative identity.
    #[inline]
    pub const fn one() -> Self {
        R
    }

    /// Checks whether `self` is zero or not
    pub fn is_zero(&self) -> Choice {
        self.ct_eq(&Fp::zero())
    }

    #[inline(always)]
    pub(crate) const fn montgomery_reduce(r0: u64, r1: u64) -> Self {
        let k = r0.wrapping_mul(U);
        let (_, carry) = mul64_with_carry(r0, k, M.0, 0);
        let (r1, _) = add64_with_carry(r1, 0, carry);

        // The result may be within M of the correct value,
        // hence substracting the modulus
        (&Fp(r1)).sub(&M)
    }

    /// Computes the summation of two field elements
    #[inline]
    pub const fn add(&self, rhs: &Self) -> Self {
        let (d0, _) = add64_with_carry(self.0, rhs.0, 0);

        // The result may be within M of the correct value,
        // hence substracting the modulus
        (&Fp(d0)).sub(&M)
    }

    /// Computes the difference of two field elements
    #[inline]
    pub const fn sub(&self, rhs: &Self) -> Self {
        let (d0, borrow) = sub64_with_carry(self.0, rhs.0, 0);

        // If underflow occurred,
        // borrow = 0xfff...fff, otherwise borrow = 0x000...000.
        let (d0, _) = add64_with_carry(d0, M.0 & borrow, 0);

        Fp(d0)
    }

    /// Computes the negation of a field element
    #[inline]
    pub const fn neg(&self) -> Self {
        // Subtract `self` from `M` to negate. Ignore the borrow
        // because it cannot underflow; self is guaranteed to
        // be in the field.
        let (d0, _) = sub64_with_carry(M.0, self.0, 0);

        // `tmp` could be `M` if `self` was zero. Create a mask that is
        // zero if `self` was zero, and `u64::max_value()` if self was nonzero.
        let mask = ((self.0 == 0) as u64).wrapping_sub(1);

        Fp(d0 & mask)
    }

    /// Computes the multiplication of two field elements
    #[inline]
    pub const fn mul(&self, rhs: &Self) -> Self {
        let (r0, r1) = mul64_with_carry(0, self.0, rhs.0, 0);

        Fp::montgomery_reduce(r0, r1)
    }

    /// Computes the square of a field element
    #[inline]
    pub const fn square(&self) -> Self {
        self.mul(self)
    }

    /// Computes the square root of this element, if it exists.
    pub fn sqrt(&self) -> CtOption<Self> {
        // Tonelli-Shank's algorithm for q mod 16 = 1
        // See https://eprint.iacr.org/2020/1497.pdf, page 3 for a
        // constant time specification of the algorithm.

        // Compute the progenitor y of self
        // y = self^((t - 1) // 2)
        //   = self^0x41
        let y = self.exp_vartime(0x41);

        let mut s = self * y;
        let mut t = s * y;

        let mut z = TWO_ADIC_ROOT_OF_UNITY;

        for k in (2..=TWO_ADICITY).rev() {
            let mut b = t;

            square_assign_multi(&mut b, (k - 2) as usize);

            let new_s = s * z;
            s = Fp::conditional_select(&new_s, &s, b.ct_eq(&Fp::one()));
            z = z.square();
            let new_t = t * z;
            t = Fp::conditional_select(&new_t, &t, b.ct_eq(&Fp::one()));
        }

        CtOption::new(s, (s * s).ct_eq(self))
    }

    /// Computes the double of a field element
    #[inline]
    pub const fn double(&self) -> Self {
        // Left-shifting the current value cannot overflow, provided that
        // the element has been constructed via a safe method, or that it
        // has been checked prior using an unchecked one.

        // The doubled value may be within M of the correct value,
        // hence substracting the modulus
        (&Fp(self.0 << 1)).sub(&M)
    }

    /// Outputs the internal representation as a 64-bit limb after Montgomery reduction
    pub const fn output_reduced_limbs(&self) -> u64 {
        Fp::montgomery_reduce(self.0, 0).0
    }

    /// Outputs the internal representation as a 64-bit limb without Montgomery reduction
    /// This is intended for uses like re-interpreting the type containing the internal value.
    pub const fn output_unreduced_limbs(&self) -> u64 {
        self.0
    }

    /// Converts an `Fp` element into a byte representation in
    /// little-endian byte order.
    pub fn to_bytes(&self) -> [u8; 8] {
        // Turn into canonical form by computing
        // (a.R) / R = a
        let tmp = Fp::montgomery_reduce(self.0, 0);

        tmp.0.to_le_bytes()
    }

    /// Attempts to convert a little-endian byte representation of
    /// a scalar into a `Fp` element, failing if the input is not canonical.
    pub fn from_bytes(bytes: &[u8; 8]) -> CtOption<Self> {
        let mut tmp = Fp(u64::from_le_bytes(*bytes));

        // Try to subtract the modulus M
        let (_, borrow) = sub64_with_carry(tmp.0, M.0, 0);

        // If the element is smaller than M then the
        // subtraction will underflow, producing a borrow value
        // of 0xffff...ffff. Otherwise, it'll be zero.
        let is_some = (borrow as u8) & 1;

        // Convert to Montgomery form by computing
        // (a.R^0 * R^2) / R = a.R
        tmp *= &R2;

        CtOption::new(tmp, Choice::from(is_some))
    }

    /// Converts a 128-bit little endian integer into
    /// a `Fp` by reducing by the modulus.
    pub fn from_bytes_wide(bytes: &[u8; 16]) -> Self {
        Fp::from_u128([
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[0..8]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[8..16]).unwrap()),
        ])
    }

    fn from_u128(limbs: [u64; 2]) -> Self {
        // We reduce an arbitrary 128-bit number by decomposing it into two 64-bit digits
        // with the higher bits multiplied by 2^64. Thus, we perform two reductions
        //
        // 1. the lower bits are multiplied by R^2, as normal
        // 2. the upper bits are multiplied by R^2 * 2^64 = R^3
        //
        // and computing their sum in the field.
        // The reduction works so long as the product is less than R=2^64 multiplied by
        // the modulus. This holds because for any `c` smaller than the modulus, we have
        // that (2^64 - 1)*c is an acceptable product for the reduction. Therefore, the
        // reduction always works so long as `c` is in the field; in this case it is either the
        // constant `R2` or `R3`.
        let d0 = Fp(limbs[0]);
        let d1 = Fp(limbs[1]);

        // Convert to Montgomery form
        d0 * R2 + d1 * R3
    }

    /// Returns whether or not this element is strictly lexicographically
    /// larger than its negation.
    pub fn lexicographically_largest(&self) -> Choice {
        // This can be determined by checking to see if the element is
        // larger than (p - 1) // 2. If we subtract by ((p - 1) // 2) + 1
        // and there is no underflow, then the element must be larger than
        // (p - 1) // 2.

        // First, because self is in Montgomery form we need to reduce it
        let tmp = Fp::montgomery_reduce(self.0, 0);

        // (p-1) // 2 + 1 = 0x20c0000000000001
        let (_, borrow) = sub64_with_carry(tmp.0, 0x20c0000000000001, 0);

        // If the element was smaller, the subtraction will underflow
        // producing a borrow value of 0xffff...ffff, otherwise it will
        // be zero. We create a Choice representing true if there was
        // overflow (and so this element is not lexicographically larger
        // than its negation) and then negate it.

        !Choice::from((borrow as u8) & 1)
    }

    /// Exponentiates `self` by `power`, where `power` is a
    /// little-endian order integer exponent.
    pub fn exp(self, power: u64) -> Self {
        let mut res = Self::one();
        for i in (0..64).rev() {
            res = res.square();
            let mut tmp = res;
            tmp *= self;
            res.conditional_assign(&tmp, (((power >> i) & 1) as u8).into());
        }
        res
    }

    /// Exponentiates `self` by `power`, where `power` is a
    /// little-endian order integer exponent.
    ///
    /// **This operation is variable time with respect
    /// to the exponent.** If the exponent is fixed,
    /// this operation is effectively constant time.
    pub fn exp_vartime(&self, power: u64) -> Self {
        let mut res = Self::one();
        for i in (0..64).rev() {
            res = res.square();

            if ((power >> i) & 1) == 1 {
                res.mul_assign(self);
            }
        }
        res
    }

    /// Computes the multiplicative inverse of this element,
    /// failing if the element is zero.
    pub fn invert(&self) -> CtOption<Self> {
        // found using https://github.com/kwantam/addchain for M - 2
        let mut t2 = self.square(); //          1: 2
        let mut t0 = t2 * self; //              2: 3
        let t4 = t0.square(); //                3: 6
        let t3 = t4.square(); //                4: 12
        let mut t1 = t3.square(); //            5: 24
        t1 = t1.square(); //                    6: 48
        let mut t3 = t1 * t3; //                7: 60
        t1 = t3.square(); //                    8: 120
        t1 *= t4; //                            9: 126
        square_assign_multi(&mut t1, 5); //    14: 4032
        t3 = t1 * t3; //                       15: 4092
        t1 = t3.square(); //                   16: 8184
        t1 *= t4; //                           17: 8190
        square_assign_multi(&mut t1, 11); //   28: 16773120
        t1 *= t3; //                           29: 16777212
        t1 *= t2; //                           30: 16777214
        t0 = t1 * t0; //                       31: 16777217
        t2 = t0 * t1; //                       32: 33554431
        t0 = t2 * t0; //                       33: 50331648
        t3 = t0.square(); //                   34: 100663296
        t1 = t3.square(); //                   35: 201326592
        t1 *= t0; //                           36: 251658240
        square_assign_multi(&mut t1, 3); //    39: 2013265920
        t1 *= t3; //                           40: 2113929216
        t1 *= t2; //                           41: 2147483647
        t0 = t1 * t0; //                       42: 2197815295
        square_assign_multi(&mut t0, 31); //   73: 4719772407336796160
        t0 *= t1; //                           74: 4719772409484279807 = M - 2

        CtOption::new(t0, !self.ct_eq(&Self::zero()))
    }

    /// Constructs an element of `Fp` without checking that it is
    /// canonical.
    pub const fn from_raw_unchecked(v: u64) -> Self {
        Fp(v)
    }

    /// Outputs a `Fp` element of multiplicative order equals to 2^n
    pub fn get_root_of_unity(n: u32) -> Self {
        assert!(n != 0, "cannot get root of unity for n = 0");
        assert!(n <= TWO_ADICITY, "order cannot exceed 2^{}", TWO_ADICITY);
        let power = 1u64 << (TWO_ADICITY - n);

        TWO_ADIC_ROOT_OF_UNITY.exp(power)
    }

    /// Outputs a `Fp` element of multiplicative order equals to 2^n
    pub fn get_root_of_unity_vartime(n: u32) -> Self {
        assert!(n != 0, "cannot get root of unity for n = 0");
        assert!(n <= TWO_ADICITY, "order cannot exceed 2^{}", TWO_ADICITY);
        let power = 1u64 << (TWO_ADICITY - n);

        TWO_ADIC_ROOT_OF_UNITY.exp_vartime(power)
    }
}

// OVERLOADED OPERATORS
// ================================================================================================

impl<'a> Neg for &'a Fp {
    type Output = Fp;

    #[inline]
    fn neg(self) -> Fp {
        self.neg()
    }
}

impl Neg for Fp {
    type Output = Fp;

    #[inline]
    fn neg(self) -> Fp {
        -&self
    }
}

impl<'a, 'b> Sub<&'b Fp> for &'a Fp {
    type Output = Fp;

    #[inline]
    fn sub(self, rhs: &'b Fp) -> Fp {
        self.sub(rhs)
    }
}

impl<'a, 'b> Add<&'b Fp> for &'a Fp {
    type Output = Fp;

    #[inline]
    fn add(self, rhs: &'b Fp) -> Fp {
        self.add(rhs)
    }
}

impl<'a, 'b> Mul<&'b Fp> for &'a Fp {
    type Output = Fp;

    #[inline]
    fn mul(self, rhs: &'b Fp) -> Fp {
        self.mul(rhs)
    }
}

impl_binops_additive!(Fp, Fp);
impl_binops_multiplicative!(Fp, Fp);

impl Div for Fp {
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Fp {
        self.mul(&rhs.invert().unwrap_or(Self::zero()))
    }
}

impl DivAssign for Fp {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs
    }
}

// TYPE CONVERSIONS
// ================================================================================================

impl From<u64> for Fp {
    /// Converts a 64-bit value into a field element. If the value is greater than or equal to
    /// the field modulus, modular reduction is silently performed.
    fn from(value: u64) -> Self {
        Fp::new(value)
    }
}

impl From<u32> for Fp {
    /// Converts a 32-bit value into a field element.
    fn from(value: u32) -> Self {
        Fp::new(value as u64)
    }
}

impl From<u16> for Fp {
    /// Converts a 16-bit value into a field element.
    fn from(value: u16) -> Self {
        Fp::new(value as u64)
    }
}

impl From<u8> for Fp {
    /// Converts an 8-bit value into a field element.
    fn from(value: u8) -> Self {
        Fp::new(value as u64)
    }
}

// FIELD TRAITS IMPLEMENTATION
// ================================================================================================

impl Field for Fp {
    fn random(mut rng: impl RngCore) -> Self {
        Fp::new(rng.next_u64())
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

impl PrimeField for Fp {
    type Repr = [u8; 8];

    fn from_repr(r: Self::Repr) -> CtOption<Self> {
        Self::from_bytes(&r)
    }

    fn to_repr(&self) -> Self::Repr {
        self.to_bytes()
    }

    fn is_odd(&self) -> Choice {
        (self.to_bytes()[0] & 1).ct_eq(&1)
    }

    const NUM_BITS: u32 = 63;
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
impl Serialize for Fp {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        use serde::ser::SerializeTuple;
        let mut tup = serializer.serialize_tuple(8)?;
        for byte in self.to_bytes().iter() {
            tup.serialize_element(byte)?;
        }
        tup.end()
    }
}

#[cfg(feature = "serialize")]
impl<'de> Deserialize<'de> for Fp {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct FpVisitor;

        impl<'de> Visitor<'de> for FpVisitor {
            type Value = Fp;

            fn expecting(&self, formatter: &mut Formatter) -> fmt::Result {
                formatter.write_str("a valid field element")
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<Fp, A::Error>
            where
                A: serde::de::SeqAccess<'de>,
            {
                let mut bytes = [0u8; 8];
                for (i, byte) in bytes.iter_mut().enumerate() {
                    *byte = seq
                        .next_element()?
                        .ok_or_else(|| serde::de::Error::invalid_length(i, &"expected 8 bytes"))?;
                }
                let elem = Fp::from_bytes(&bytes);
                if bool::from(elem.is_none()) {
                    Err(serde::de::Error::custom("decompression failed"))
                } else {
                    Ok(elem.unwrap())
                }
            }
        }

        deserializer.deserialize_tuple(8, FpVisitor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::OsRng;

    const LARGEST: Fp = Fp(4719772409484279808);
    const TWO_POW_55: u64 = 36028797018963968;
    const TWO_POW_54: u64 = 18014398509481984;

    // DISPLAY
    // ================================================================================================

    #[test]
    fn test_debug() {
        assert_eq!(format!("{:?}", Fp::zero()), "0");
        assert_eq!(format!("{:?}", Fp::one()), "1");
        assert_eq!(format!("{:?}", R2), "4287426845256712189");
    }

    #[test]
    fn test_output_reduced_limbs() {
        assert_eq!(format!("{:?}", Fp::zero().output_reduced_limbs()), "0");
        assert_eq!(format!("{:?}", Fp::one().output_reduced_limbs()), "1");
        assert_eq!(
            format!("{:?}", R2.output_reduced_limbs()),
            "4287426845256712189"
        );

        let r = Fp::from_raw_unchecked(4287426845256712189);
        assert!(format!("{:?}", r.0) != format!("{:?}", R.output_reduced_limbs()));

        assert_eq!(
            format!("{:?}", r.output_reduced_limbs()),
            format!("{:?}", R.output_reduced_limbs())
        );
    }

    // BASIC ALGEBRA
    // ================================================================================================

    #[test]
    fn test_conditional_selection() {
        let a = Fp::from_raw_unchecked(1);
        let b = Fp::from_raw_unchecked(2);

        assert_eq!(
            ConditionallySelectable::conditional_select(&a, &b, Choice::from(0u8)),
            a
        );
        assert_eq!(
            ConditionallySelectable::conditional_select(&a, &b, Choice::from(1u8)),
            b
        );
    }

    #[test]
    fn test_equality() {
        assert_eq!(Fp::default(), Fp::zero());
        assert_eq!(Fp::zero(), Fp::zero());
        assert_eq!(Fp::one(), Fp::one());

        assert!(bool::from(Fp::default().is_zero()));
        assert!(bool::from(Fp::zero().ct_eq(&Fp::zero())));
        assert!(!bool::from(Fp::zero().ct_eq(&Fp::one())));

        assert!(Fp::zero() != Fp::one());
        assert!(Fp::one() != R2);
        assert!(Fp::one() == R);

        assert!(!Fp::zero().eq(&Fp::one()));
        assert!(!Fp::one().eq(&R2));
        assert!(Fp::one().eq(&R));

        assert_eq!(Fp::zero(), Fp::new(0));
        assert_eq!(Fp::one(), Fp::new(1));
    }

    #[test]
    fn test_addition() {
        let mut tmp = LARGEST;
        tmp += &LARGEST;

        assert_eq!(tmp, Fp(4719772409484279807));

        assert_eq!(tmp, LARGEST.double());

        let mut tmp = LARGEST;
        tmp += &Fp(1);

        assert_eq!(tmp, Fp::zero());
    }

    #[test]
    fn test_subtraction() {
        let mut tmp = LARGEST;
        tmp -= &LARGEST;

        assert_eq!(tmp, Fp::zero());

        let mut tmp = Fp::zero();
        tmp -= &LARGEST;

        let mut tmp2 = M;
        tmp2 -= &LARGEST;

        assert_eq!(tmp, tmp2);
    }

    #[test]
    fn test_negation() {
        let mut rng = OsRng;

        let tmp = -&LARGEST;

        assert_eq!(tmp, Fp(1));

        let tmp = -&Fp::zero();
        assert_eq!(tmp, Fp::zero());
        let tmp = -&Fp(1);
        assert_eq!(tmp, LARGEST);

        for _ in 0..100 {
            let a = Fp::random(&mut rng);
            let b = -a;

            assert_eq!(a + b, Fp::zero());
        }
    }

    #[test]
    fn test_multiplication() {
        let mut cur = LARGEST;

        for _ in 0..100 {
            let mut tmp = cur;
            tmp *= &cur;

            let mut tmp2 = Fp::zero();
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
    fn test_division() {
        let mut rng = OsRng;
        let scale = Fp::random(&mut rng);

        for _ in 0..100 {
            let a = Fp::random(&mut rng);
            let a_scaled = a * scale;

            assert_eq!(a, a_scaled.div(scale));
        }

        assert_eq!(Fp::one().div(Fp::zero()), Fp::zero());
    }

    #[test]

    fn test_inversion() {
        assert!(bool::from(Fp::zero().invert().is_none()));
        assert_eq!(Fp::one().invert().unwrap(), Fp::one());
        assert_eq!((-&Fp::one()).invert().unwrap(), -&Fp::one());

        let mut tmp = Fp::random(&mut OsRng);

        for _ in 0..100 {
            let mut tmp2 = tmp.invert().unwrap();
            tmp2.mul_assign(&tmp);

            assert_eq!(tmp2, Fp::one());

            tmp.add_assign(&Fp::random(&mut OsRng));
        }
    }

    #[test]
    fn test_squaring() {
        let mut cur = LARGEST;

        for _ in 0..100 {
            let mut tmp = cur;
            let pow2 = tmp.exp(2);
            tmp = tmp.square();

            let mut tmp2 = Fp::zero();
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
            let a = Fp::random(&mut OsRng).square();
            let b = a.sqrt().unwrap();
            assert_eq!(a, b.square());
        }

        assert_eq!(Fp::zero().sqrt().unwrap(), Fp::zero());
        assert_eq!(Fp::one().sqrt().unwrap(), Fp::one());

        // 3 is not a quadratic residue in Fp
        assert!(bool::from(Fp::new(3).sqrt().is_none()));
    }

    #[test]
    fn test_invert_is_pow() {
        let p_minus_2 = 4719772409484279807;

        let mut r1 = R;
        let mut r2 = R;
        let mut r3 = R;

        for _ in 0..100 {
            r1 = r1.invert().unwrap();
            r2 = r2.exp(p_minus_2);
            r3 = r3.exp_vartime(p_minus_2);

            assert_eq!(r1, r2);
            assert_eq!(r2, r3);

            // Double so we check a different element each time
            r1 = r1.double();
            r2 = r1;
            r3 = r1;
        }
    }

    // ROOTS OF UNITY
    // ================================================================================================

    #[test]
    fn test_get_root_of_unity() {
        let root_55 = Fp::get_root_of_unity(55);
        let root_55_vartime = Fp::get_root_of_unity_vartime(55);
        assert_eq!(TWO_ADIC_ROOT_OF_UNITY, root_55);
        assert_eq!(TWO_ADIC_ROOT_OF_UNITY, root_55_vartime);
        assert_eq!(Fp::one(), root_55.exp(TWO_POW_55));

        let root_54 = Fp::get_root_of_unity(54);
        let root_54_vartime = Fp::get_root_of_unity_vartime(54);
        let expected = root_55.exp(2);
        assert_eq!(expected, root_54);
        assert_eq!(expected, root_54_vartime);
        assert_eq!(Fp::one(), root_54.exp(TWO_POW_54));
    }

    #[test]
    #[should_panic]
    fn test_get_root_of_unity_zero() {
        let _ = Fp::get_root_of_unity(0);
    }

    #[test]
    #[should_panic]
    fn test_get_root_of_unity_too_large() {
        let _ = Fp::get_root_of_unity(TWO_ADICITY + 1);
    }

    #[test]
    fn test_lexicographically_largest() {
        let a = Fp::new(1565036539327771067);

        let b = Fp::new(3154735870156508742);

        assert_eq!(a.square(), b.square());
        assert!(!bool::from(a.lexicographically_largest()));
        assert!(bool::from(b.lexicographically_largest()));
    }

    #[test]
    fn test_zeroize() {
        use zeroize::Zeroize;

        let mut a = Fp::one();
        a.zeroize();
        assert_eq!(a, Fp::zero());
    }

    // INITIALIZATION
    // ================================================================================================

    #[test]
    fn test_from_int() {
        let n = 42u8;
        let element = Fp::from(n);

        assert_eq!(element, Fp::from(n as u16));
        assert_eq!(element, Fp::from(n as u32));
        assert_eq!(element, Fp::from(n as u64));

        assert_eq!(element, Fp::new(n as u64));
    }

    #[test]
    fn test_from_raw_unchecked() {
        let mut element = Fp::from_raw_unchecked(4287426845256712189);

        let element_normalized = Fp::new(4287426845256712189);

        assert_eq!(element, Fp::one());
        element *= &R2;

        assert!(element != Fp::one());
        assert_eq!(element, element_normalized);
    }

    // FIELD TRAIT
    // ================================================================================================

    #[test]
    fn test_field_trait_methods() {
        assert_eq!(<Fp as Field>::zero(), Fp::new(0));
        assert_eq!(<Fp as Field>::one(), Fp::new(1));

        assert_eq!(
            bool::from(<Fp as Field>::zero().is_zero()),
            bool::from(Fp::new(0).is_zero())
        );
        assert_eq!(
            bool::from(<Fp as Field>::one().is_zero()),
            bool::from(Fp::new(1).is_zero())
        );

        let mut rng = OsRng;
        let e = Fp::random(&mut rng).square();

        assert_eq!(<Fp as Field>::square(&e), e.square());
        assert_eq!(<Fp as Field>::double(&e), e.double());

        assert_eq!(<Fp as Field>::invert(&e).unwrap(), e.invert().unwrap());
        assert!(bool::from(<Fp as Field>::invert(&Fp::zero()).is_none()));

        assert_eq!(<Fp as Field>::sqrt(&e).unwrap(), e.sqrt().unwrap());
        assert!(bool::from(<Fp as Field>::sqrt(&Fp::new(3)).is_none()));
    }

    #[test]
    fn test_primefield_trait_methods() {
        assert_eq!(Fp::from_repr([0, 0, 0, 0, 0, 0, 0, 0]).unwrap(), Fp::zero());
        assert_eq!(Fp::from_repr([1, 0, 0, 0, 0, 0, 0, 0]).unwrap(), Fp::one());

        let mut rng = OsRng;
        let e = Fp::random(&mut rng).square();

        assert_eq!(Fp::from_repr(Fp::to_repr(&e)).unwrap(), e);

        assert_eq!(
            Fp::get_root_of_unity(<Fp as PrimeField>::S),
            <Fp as PrimeField>::root_of_unity()
        )
    }

    // SERIALIZATION / DESERIALIZATION
    // ================================================================================================

    #[test]
    fn test_to_bytes() {
        assert_eq!(Fp::zero().to_bytes(), [0, 0, 0, 0, 0, 0, 0, 0]);

        assert_eq!(Fp::one().to_bytes(), [1, 0, 0, 0, 0, 0, 0, 0]);

        assert_eq!(R2.to_bytes(), [253, 255, 255, 255, 255, 255, 127, 59]);

        assert_eq!((-&Fp::one()).to_bytes(), [0, 0, 0, 0, 0, 0, 128, 65]);
    }

    #[test]
    fn test_from_bytes() {
        let mut rng = OsRng;
        for _ in 0..100 {
            let a = Fp::random(&mut rng);
            let bytes = a.to_bytes();
            assert_eq!(a, Fp::from_bytes(&bytes).unwrap());
        }

        assert_eq!(
            Fp::from_bytes(&[0, 0, 0, 0, 0, 0, 0, 0]).unwrap(),
            Fp::zero()
        );

        assert_eq!(
            Fp::from_bytes(&[1, 0, 0, 0, 0, 0, 0, 0]).unwrap(),
            Fp::one()
        );

        assert_eq!(
            Fp::from_bytes(&[253, 255, 255, 255, 255, 255, 127, 59]).unwrap(),
            R2
        );

        // -1 should work
        assert_eq!(
            Fp::from_bytes(&[0, 0, 0, 0, 0, 0, 128, 65]).unwrap(),
            -Fp::one()
        );

        // M is invalid
        assert!(bool::from(
            Fp::from_bytes(&[1, 0, 0, 0, 0, 0, 128, 65]).is_none()
        ));

        // Anything larger than M is invalid
        assert!(bool::from(
            Fp::from_bytes(&[2, 0, 0, 0, 0, 0, 128, 65]).is_none()
        ));

        assert!(bool::from(
            Fp::from_bytes(&[0, 0, 0, 0, 255, 255, 255, 255]).is_none()
        ));
    }

    #[test]
    fn test_from_u128_max() {
        let max_u64 = 0xffff_ffff_ffff_ffff;
        assert_eq!(R3 - R, Fp::from_u128([max_u64, max_u64]));
    }

    #[test]
    fn test_from_bytes_wide_r2() {
        assert_eq!(
            R2,
            Fp::from_bytes_wide(&[253, 255, 255, 255, 255, 255, 127, 59, 0, 0, 0, 0, 0, 0, 0, 0])
        );
    }

    #[test]
    fn test_from_bytes_wide_negative_one() {
        assert_eq!(
            -&Fp::one(),
            Fp::from_bytes_wide(&[0, 0, 0, 0, 0, 0, 128, 65, 0, 0, 0, 0, 0, 0, 0, 0])
        );
    }

    #[test]
    fn test_from_bytes_wide_maximum() {
        assert_eq!(Fp(0x551e3ce4b73f407), Fp::from_bytes_wide(&[0xff; 16]));
    }

    #[test]
    #[cfg(feature = "serialize")]
    fn test_serde() {
        let mut rng = OsRng;
        let element = Fp::random(&mut rng);
        let encoded = bincode::serialize(&element).unwrap();
        let parsed: Fp = bincode::deserialize(&encoded).unwrap();
        assert_eq!(parsed, element);

        // Check that the encoding is 8 bytes exactly
        assert_eq!(encoded.len(), 8);

        // Check that the encoding itself matches the usual one
        assert_eq!(element, bincode::deserialize(&element.to_bytes()).unwrap());

        // Check that invalid encodings fail
        let element = Fp::random(&mut rng);
        let mut encoded = bincode::serialize(&element).unwrap();
        encoded[7] = 127;
        assert!(bincode::deserialize::<Fp>(&encoded).is_err());

        let encoded = bincode::serialize(&element).unwrap();
        assert!(bincode::deserialize::<Fp>(&encoded[0..7]).is_err());
    }
}
