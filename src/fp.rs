// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This module provides an implementation of the finite prime
//! field Fp of characteristic p = 2^64 - 2^32 + 1.

use core::{
    borrow::Borrow,
    fmt::{self, Debug, Display, Formatter},
    hash::{Hash, Hasher},
    iter::Sum,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::utils::{shl64_by_u32_with_carry, square_assign_multi, sub64_with_carry};

use group::ff::{Field, PrimeField};
use rand_core::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "serialize")]
use serde::de::Visitor;
#[cfg(feature = "serialize")]
use serde::{self, Deserialize, Deserializer, Serialize, Serializer};

// CONSTANTS
// ================================================================================================

// Field modulus = 2^64 - 2^32 + 1
const M: Fp = Fp(0xffffffff00000001);

// Multiplicative generator g of order p-1
// g = 7
pub(crate) const GENERATOR: Fp = Fp(7);

// Epsilon = 2^32 - 1;
const E: u64 = 0xffffffff;

// Two-adicity of the field: (p-1) % 2^32 = 0
pub(crate) const TWO_ADICITY: u32 = 32;

// 2^32 root of unity = 1753635133440165772
const TWO_ADIC_ROOT_OF_UNITY: Fp = Fp(1753635133440165772);

// FIELD ELEMENT
// ================================================================================================

/// Represents a base field element.
///
/// The backing type is `u64`.
#[derive(Copy, Clone, Eq, Default)]
pub struct Fp(pub(crate) u64);

impl Debug for Fp {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        let tmp = self.output_internal();
        write!(f, "{:?}", tmp)
    }
}

impl Display for Fp {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl Hash for Fp {
    fn hash<H: Hasher>(&self, hasher: &mut H) {
        self.0.hash(hasher);
    }
}

impl ConstantTimeEq for Fp {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.make_canonical().0.ct_eq(&other.make_canonical().0)
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
        Self(u64::conditional_select(&a.0, &b.0, choice))
    }
}

impl zeroize::DefaultIsZeroes for Fp {}

impl Fp {
    /// Creates a new field element from a `u64` value,
    /// reduced by M if necessary.
    pub const fn new(value: u64) -> Self {
        Self(value % M.0)
    }

    /// Returns zero, the additive identity.
    #[inline]
    pub const fn zero() -> Self {
        Self(0)
    }

    /// Returns one, the multiplicative identity.
    #[inline]
    pub const fn one() -> Self {
        Self(1)
    }

    /// Checks whether `self` is zero or not
    pub fn is_zero(&self) -> Choice {
        self.ct_eq(&Fp::zero())
    }

    #[inline(always)]
    /// Makes the element canonical by reducing by the modulus if needed
    pub const fn make_canonical(&self) -> Self {
        Self(self.0 % M.0)
    }

    /// Generates a random canonical element
    pub fn random(mut rng: impl RngCore) -> Self {
        Self::new(rng.next_u64())
    }

    /// Computes the summation of two field elements
    #[inline]
    pub const fn add(&self, rhs: &Self) -> Self {
        let (d0, is_overflow) = self.0.overflowing_add(rhs.0);
        let (d0, is_overflow) = d0.overflowing_add(E * (is_overflow as u64));

        Self(d0 + E * (is_overflow as u64))
    }

    /// Computes the double of a field element
    #[inline]
    pub const fn double(&self) -> Self {
        let (d0, is_overflow) = shl64_by_u32_with_carry(self.0, 1, 0);
        let (d0, is_overflow) = d0.overflowing_add(E * (is_overflow as u64));

        Self(d0 + E * (is_overflow as u64))
    }

    /// Computes the triple of a field element
    #[inline]
    pub const fn triple(&self) -> Self {
        let t = self.0 as u128;
        let t = t + (t << 1);

        Self(reduce_u96(t))
    }

    /// Computes the difference of two field elements
    #[inline]
    pub const fn sub(&self, rhs: &Self) -> Self {
        let (d0, is_overflow) = self.0.overflowing_sub(rhs.0);
        let (d0, is_overflow) = d0.overflowing_sub(E * (is_overflow as u64));

        Self(d0 - E * (is_overflow as u64))
    }

    /// Computes the negation of a field element
    #[inline]
    pub const fn neg(&self) -> Self {
        (&Self::zero()).sub(self)
    }

    /// Computes the multiplication of two field elements
    #[inline]
    pub const fn mul(&self, rhs: &Self) -> Self {
        let r0 = (self.0 as u128) * (rhs.0 as u128);

        Self(reduce_u128(r0))
    }

    /// Computes the multiplication of a field element with a u32 value
    #[inline]
    pub const fn mul_by_u32(&self, rhs: u32) -> Self {
        let r0 = (self.0 as u128) * (rhs as u128);

        Self(reduce_u96(r0))
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
        //   = self^0x7fffffff
        let y = self.exp_vartime(0x7fffffff);

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

    /// Outputs the internal representation as
    /// a 64-bit limb after canonical reduction.
    pub const fn output_internal(&self) -> u64 {
        self.make_canonical().0
    }

    /// Converts an `Fp` element into a byte representation in
    /// little-endian byte order.
    pub const fn to_bytes(&self) -> [u8; 8] {
        // Turn into canonical form by removing modulus
        // if self is greater.
        let tmp = self.make_canonical();

        tmp.0.to_le_bytes()
    }

    /// Attempts to convert a little-endian byte representation of
    /// a scalar into a `Fp` element, failing if the input is not canonical.
    pub fn from_bytes(bytes: &[u8; 8]) -> CtOption<Self> {
        let mut tmp = Self(u64::from_le_bytes(*bytes));

        // Try to subtract the modulus M
        let (_, borrow) = sub64_with_carry(tmp.0, M.0, 0);

        // If the element is smaller than M then the
        // subtraction will underflow, producing a borrow value
        // of 0xffff...ffff. Otherwise, it'll be zero.
        let is_some = (borrow as u8) & 1;

        // Convert to canonical form
        tmp = tmp.make_canonical();

        CtOption::new(tmp, Choice::from(is_some))
    }

    /// Converts a 128-bit little endian integer into
    /// a `Fp` element by reducing by the modulus.
    ///
    /// The result is always returned in canonical form,
    /// reduced by p if necessary.
    pub fn from_bytes_wide(bytes: [u8; 16]) -> Self {
        Self(reduce_u128(u128::from_le_bytes(bytes))).make_canonical()
    }

    /// Returns whether or not this element is strictly lexicographically
    /// larger than its negation.
    pub fn lexicographically_largest(&self) -> Choice {
        // This can be determined by checking to see if the element is
        // larger than (p - 1) // 2. If we subtract by ((p - 1) // 2) + 1
        // and there is no underflow, then the element must be larger than
        // (p - 1) // 2.

        // First, because self may not be canonical, we need to make_canonical it
        let tmp = self.make_canonical();

        // (p-1) // 2 + 1 = 0x7fffffff80000001
        let (_, borrow) = sub64_with_carry(tmp.0, 0x7fffffff80000001, 0);

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
        let mut t2 = self.square(); //        1: 2
        let mut t3 = t2 * self; //            2: 3
        let mut t0 = t3.square(); //          3: 6
        t0 = t0.square(); //                  4: 12
        t2 *= t0; //                          5: 14
        t3 *= t0; //                          6: 15
        t0 = t3.square(); //                  7: 30
        square_assign_multi(&mut t0, 3); //  10: 240
        t2 *= t0; //                         11: 254
        t3 *= t0; //                         12: 255
        t0 = t3.square(); //                 13: 510
        square_assign_multi(&mut t0, 7); //  20: 65280
        t2 *= t0; //                         21: 65534
        t0 *= t3; //                         22: 65535
        square_assign_multi(&mut t0, 16); // 38: 4294901760
        t2 *= t0; //                         39: 4294967294
        t0 = t2.square(); //                 40: 8589934588
        square_assign_multi(&mut t0, 31); // 71: 18446744065119617024
        t0 *= t2; //                         72: 18446744069414584318
        t0 *= self; //                       73: 18446744069414584319 = M - 2

        CtOption::new(t0, !self.ct_eq(&Self::zero()))
    }

    /// Constructs an element of `Fp` without checking that it is
    /// canonical.
    pub const fn from_raw_unchecked(v: u64) -> Self {
        Self(v)
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

impl<T> Sum<T> for Fp
where
    T: Borrow<Fp>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Self::zero(), |acc, item| acc + item.borrow())
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

    const NUM_BITS: u32 = 64;
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

/// Reduces a 128-bit value by M such that the output fits in a u64.
#[inline(always)]
pub(crate) const fn reduce_u128(x: u128) -> u64 {
    // See https://github.com/mir-protocol/plonky2/blob/main/plonky2.pdf
    // for a more complete description of the reduction.

    // Decompose x = a + b.2^32 + c.2^64 + d.2^96 with a,b,c and d u32 values
    let ab = x as u64;
    let cd = (x >> 64) as u64;
    let c = (cd as u32) as u64;
    let d = cd >> 32;

    // r0 = ab - d
    let (r0, is_overflow) = ab.overflowing_sub(d);
    // d > ab may happen, hence handling potential overflow
    let r0 = r0.wrapping_sub(E * (is_overflow as u64));

    // r1 = c * 2^32 - c
    // this cannot underflow
    let r1 = (c << 32) - c;

    // result = r0 + r1
    let (result, is_overflow) = r0.overflowing_add(r1);
    // handle potential overflow
    result.wrapping_add(E * (is_overflow as u64))
}

/// Reduces a 96-bit value (stored as u128) by M such that the output fits in a u64.
///
/// This is similar to reduce_u128() but is aimed to be used when we are guaranteed that
/// the value to be reduced is fitting in 96 bits.
#[inline(always)]
pub(crate) const fn reduce_u96(x: u128) -> u64 {
    // Decompose x = r0 + c.2^64 with r0 a u64 value and c a u32 value
    let c = ((x >> 64) as u32) as u64;

    let r0 = x as u64;

    // r1 = c * 2^32 - c
    // this cannot underflow
    let r1 = (c << 32) - c;

    // result = r0 + r1
    let (result, is_overflow) = r0.overflowing_add(r1);
    // handle potential overflow
    result.wrapping_add(E * (is_overflow as u64))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::OsRng;

    const LARGEST: Fp = Fp(18446744069414584320);
    const TWO_POW_32: u64 = 4294967296;
    const TWO_POW_31: u64 = 2147483648;

    // DISPLAY
    // ================================================================================================

    #[test]
    fn test_debug() {
        assert_eq!(format!("{:?}", Fp::zero()), "0");
        assert_eq!(format!("{:?}", Fp::one()), "1");
    }

    #[test]
    fn test_output_internal() {
        assert_eq!(format!("{:?}", Fp::zero().output_internal()), "0");
        assert_eq!(format!("{:?}", Fp::one().output_internal()), "1");
        assert_eq!(format!("{:?}", M.output_internal()), "0");
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

        assert!(!Fp::zero().eq(&Fp::one()));

        assert_eq!(Fp::zero(), Fp::new(0));
        assert_eq!(Fp::one(), Fp::new(1));
    }

    #[test]
    fn test_addition() {
        let mut tmp = LARGEST;
        tmp += &LARGEST;

        assert_eq!(tmp, Fp(18446744069414584319));

        assert_eq!(tmp, LARGEST.double());

        let mut tmp = LARGEST;
        tmp += &Fp(1);

        assert_eq!(tmp, Fp::zero());
    }

    #[test]
    fn test_double_and_triple() {
        let mut rng = OsRng;

        for _ in 0..100 {
            let e = Fp::random(&mut rng);
            assert_eq!(e + e, e.double());
            assert_eq!(e.triple(), e.double() + e);
        }
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

        // 7 is not a quadratic residue in Fp
        assert!(bool::from(Fp::new(7).sqrt().is_none()));
    }

    #[test]
    fn test_invert_is_pow() {
        let p_minus_2 = 18446744069414584319;

        let mut r1 = Fp::random(&mut OsRng);
        let mut r2 = r1;
        let mut r3 = r2;

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
        let root_32 = Fp::get_root_of_unity(32);
        let root_32_vartime = Fp::get_root_of_unity_vartime(32);
        assert_eq!(TWO_ADIC_ROOT_OF_UNITY, root_32);
        assert_eq!(TWO_ADIC_ROOT_OF_UNITY, root_32_vartime);
        assert_eq!(Fp::one(), root_32.exp(TWO_POW_32));

        let root_31 = Fp::get_root_of_unity(31);
        let root_31_vartime = Fp::get_root_of_unity_vartime(31);
        let expected = root_32.exp(2);
        assert_eq!(expected, root_31);
        assert_eq!(expected, root_31_vartime);
        assert_eq!(Fp::one(), root_31.exp(TWO_POW_31));
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
        let a = Fp::new(2293556613039705979);

        let b = Fp::new(16153187456374878342);

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
        let mut element = Fp::from_raw_unchecked(M.0 + 42);

        let element_normalized = Fp::new(42);

        assert_eq!(element, element_normalized);
        assert!(element.0 != element_normalized.0);

        element = element.make_canonical();
        assert_eq!(element, element_normalized);
        assert!(element.0 == element_normalized.0);
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
        assert!(bool::from(<Fp as Field>::sqrt(&Fp::new(7)).is_none()));
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

        assert_eq!((-&Fp::one()).to_bytes(), [0, 0, 0, 0, 255, 255, 255, 255]);
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

        // -1 should work
        assert_eq!(
            Fp::from_bytes(&[0, 0, 0, 0, 255, 255, 255, 255]).unwrap(),
            -Fp::one()
        );

        // M is invalid
        assert!(bool::from(
            Fp::from_bytes(&[1, 0, 0, 0, 255, 255, 255, 255]).is_none()
        ));

        // Anything larger than M is invalid
        assert!(bool::from(
            Fp::from_bytes(&[42, 0, 0, 0, 255, 255, 255, 255]).is_none()
        ));
    }

    #[test]
    fn test_from_bytes_wide_negative_one() {
        assert_eq!(
            -&Fp::one(),
            Fp::from_bytes_wide([0, 0, 0, 0, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0])
        );
    }

    #[test]
    fn test_from_bytes_wide_maximum() {
        assert_eq!(Fp(0xfffffffe00000000), Fp::from_bytes_wide([0xff; 16]));
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
        let encoded = [0xff; 8];
        assert!(bincode::deserialize::<Fp>(&encoded).is_err());

        assert!(bincode::deserialize::<Fp>(&encoded[0..7]).is_err());
    }
}
