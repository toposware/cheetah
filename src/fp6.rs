// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This module implements arithmetic over the extension field Fp6,
//! defined with irreducible polynomial u^6 - 7.

use core::{
    borrow::Borrow,
    fmt::{self, Formatter},
    hash::{Hash, Hasher},
    iter::Sum,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use group::ff::Field;
use rand_core::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "serialize")]
use serde::de::Visitor;
#[cfg(feature = "serialize")]
use serde::{self, Deserialize, Deserializer, Serialize, Serializer};

use crate::fp::reduce_u96;
use crate::fp::Fp;
use crate::utils::square_assign_multi;

use crate::fp::TWO_ADICITY;

const BETA: u128 = crate::fp::GENERATOR.0 as u128;

// 2^33 root of unity = 135436726719221772*u^3
const TWO_ADIC_ROOT_OF_UNITY_P6: Fp6 = Fp6 {
    c0: Fp::zero(),
    c1: Fp::zero(),
    c2: Fp::zero(),
    c3: Fp(135436726719221772),
    c4: Fp::zero(),
    c5: Fp::zero(),
};

const TWO_ADICITY_P6: u32 = TWO_ADICITY + 1;

#[derive(Copy, Clone)]
/// An element of the extension GF(p^6).
///
/// It represents the field extension element
/// c5.u^5 + c4.u^4 + c3.u^3 + c2.u^2 + c1.u + c0
/// where u is a root of the polynomial defining
/// the sextic extension.
pub struct Fp6 {
    /// First coefficient, lowest degree
    pub(crate) c0: Fp,
    /// Second coefficient
    pub(crate) c1: Fp,
    /// Third coefficient
    pub(crate) c2: Fp,
    /// Fourth coefficient
    pub(crate) c3: Fp,
    /// Fifth coefficient
    pub(crate) c4: Fp,
    /// Sixth coefficient, highest degree
    pub(crate) c5: Fp,
}

impl fmt::Debug for Fp6 {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(
            f,
            "{:?} + {:?}*u + {:?}*u^2 + {:?}*u^3 + {:?}*u^4 + {:?}*u^5",
            self.c0, self.c1, self.c2, self.c3, self.c4, self.c5
        )
    }
}

impl Default for Fp6 {
    fn default() -> Self {
        Self::zero()
    }
}

impl Hash for Fp6 {
    fn hash<H: Hasher>(&self, hasher: &mut H) {
        self.c0.hash(hasher);
        self.c1.hash(hasher);
        self.c2.hash(hasher);
        self.c3.hash(hasher);
        self.c4.hash(hasher);
        self.c5.hash(hasher);
    }
}

impl zeroize::DefaultIsZeroes for Fp6 {}

impl From<Fp> for Fp6 {
    fn from(f: Fp) -> Self {
        Self {
            c0: f,
            c1: Fp::zero(),
            c2: Fp::zero(),
            c3: Fp::zero(),
            c4: Fp::zero(),
            c5: Fp::zero(),
        }
    }
}

impl From<[Fp; 6]> for Fp6 {
    fn from(f: [Fp; 6]) -> Self {
        Self {
            c0: f[0],
            c1: f[1],
            c2: f[2],
            c3: f[3],
            c4: f[4],
            c5: f[5],
        }
    }
}

impl From<&Fp6> for [Fp; 6] {
    fn from(f: &Fp6) -> [Fp; 6] {
        [f.c0, f.c1, f.c2, f.c3, f.c4, f.c5]
    }
}

impl From<Fp6> for [Fp; 6] {
    fn from(f: Fp6) -> [Fp; 6] {
        Self::from(&f)
    }
}

impl From<u64> for Fp6 {
    /// Converts a 64-bit value into a field element. If the value is greater than or equal to
    /// the field modulus, modular reduction is silently performed.
    fn from(value: u64) -> Self {
        Self::from(Fp::new(value))
    }
}

impl From<u32> for Fp6 {
    /// Converts a 32-bit value into a field element.
    fn from(value: u32) -> Self {
        Self::from(Fp::new(value as u64))
    }
}

impl From<u16> for Fp6 {
    /// Converts a 16-bit value into a field element.
    fn from(value: u16) -> Self {
        Self::from(Fp::new(value as u64))
    }
}

impl From<u8> for Fp6 {
    /// Converts an 8-bit value into a field element.
    fn from(value: u8) -> Self {
        Self::from(Fp::new(value as u64))
    }
}

impl ConstantTimeEq for Fp6 {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.c0.ct_eq(&other.c0)
            & self.c1.ct_eq(&other.c1)
            & self.c2.ct_eq(&other.c2)
            & self.c3.ct_eq(&other.c3)
            & self.c4.ct_eq(&other.c4)
            & self.c5.ct_eq(&other.c5)
    }
}

impl Eq for Fp6 {}
impl PartialEq for Fp6 {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl ConditionallySelectable for Fp6 {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self {
            c0: Fp::conditional_select(&a.c0, &b.c0, choice),
            c1: Fp::conditional_select(&a.c1, &b.c1, choice),
            c2: Fp::conditional_select(&a.c2, &b.c2, choice),
            c3: Fp::conditional_select(&a.c3, &b.c3, choice),
            c4: Fp::conditional_select(&a.c4, &b.c4, choice),
            c5: Fp::conditional_select(&a.c5, &b.c5, choice),
        }
    }
}

impl<'a> Neg for &'a Fp6 {
    type Output = Fp6;

    #[inline]
    fn neg(self) -> Fp6 {
        self.neg()
    }
}

impl Neg for Fp6 {
    type Output = Fp6;

    #[inline]
    fn neg(self) -> Fp6 {
        -&self
    }
}

impl<'a, 'b> Sub<&'b Fp6> for &'a Fp6 {
    type Output = Fp6;

    #[inline]
    fn sub(self, rhs: &'b Fp6) -> Fp6 {
        self.sub(rhs)
    }
}

impl<'a, 'b> Add<&'b Fp6> for &'a Fp6 {
    type Output = Fp6;

    #[inline]
    fn add(self, rhs: &'b Fp6) -> Fp6 {
        self.add(rhs)
    }
}

impl<T> Sum<T> for Fp6
where
    T: Borrow<Fp6>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Self::zero(), |acc, item| acc + item.borrow())
    }
}

impl<'a, 'b> Mul<&'b Fp6> for &'a Fp6 {
    type Output = Fp6;

    #[inline]
    fn mul(self, rhs: &'b Fp6) -> Fp6 {
        self.mul(rhs)
    }
}

impl_binops_additive!(Fp6, Fp6);
impl_binops_multiplicative!(Fp6, Fp6);

impl Fp6 {
    /// Creates a new field element from a [u64; 6] value.
    /// The value is converted to canonical form by reducing
    /// each coordinate if necessary.
    pub const fn new(value: [u64; 6]) -> Self {
        Self {
            c0: Fp::new(value[0]),
            c1: Fp::new(value[1]),
            c2: Fp::new(value[2]),
            c3: Fp::new(value[3]),
            c4: Fp::new(value[4]),
            c5: Fp::new(value[5]),
        }
    }

    #[inline]
    /// The additive identity
    pub const fn zero() -> Self {
        Self {
            c0: Fp::zero(),
            c1: Fp::zero(),
            c2: Fp::zero(),
            c3: Fp::zero(),
            c4: Fp::zero(),
            c5: Fp::zero(),
        }
    }

    #[inline]
    /// The multiplicative identity
    pub const fn one() -> Self {
        Self {
            c0: Fp::one(),
            c1: Fp::zero(),
            c2: Fp::zero(),
            c3: Fp::zero(),
            c4: Fp::zero(),
            c5: Fp::zero(),
        }
    }

    /// Checks whether this element is zero or not
    pub fn is_zero(&self) -> Choice {
        self.c0.is_zero()
            & self.c1.is_zero()
            & self.c2.is_zero()
            & self.c3.is_zero()
            & self.c4.is_zero()
            & self.c5.is_zero()
    }

    #[inline(always)]
    /// Makes the element canonical by reducing each coordinate by the modulus if needed
    pub const fn make_canonical(&self) -> Self {
        Self {
            c0: self.c0.make_canonical(),
            c1: self.c1.make_canonical(),
            c2: self.c2.make_canonical(),
            c3: self.c3.make_canonical(),
            c4: self.c4.make_canonical(),
            c5: self.c5.make_canonical(),
        }
    }

    /// Generates a random canonical element
    pub fn random(mut rng: impl RngCore) -> Self {
        Self {
            c0: Fp::random(&mut rng),
            c1: Fp::random(&mut rng),
            c2: Fp::random(&mut rng),
            c3: Fp::random(&mut rng),
            c4: Fp::random(&mut rng),
            c5: Fp::random(&mut rng),
        }
    }

    /// Returns whether or not this element is strictly lexicographically
    /// larger than its negation.
    #[inline]
    pub fn lexicographically_largest(&self) -> Choice {
        // If this element's c5 coefficient is lexicographically largest
        // then it is lexicographically largest.
        // In the event the c5 coefficient is zero and the c4 coefficient is
        // lexicographically largest, then this element is lexicographically
        // largest.
        // Otherwise, in the event both the c5 and c4 coefficients are zero
        // and the c3 coefficient is lexicographically largest, then this
        // element is lexicographically largest.
        // Otherwise, in the event both the c5, c4 and c3 coefficients are
        // zero and the c2 coefficient is lexicographically largest, then
        // this element is lexicographically largest.
        // Otherwise, in the event both the c5, c4, c3 and c2 coefficients
        // are zero and the c1 coefficient is lexicographically largest,
        // then this element is lexicographically largest.
        // Otherwise, in the event both the c5, c4, c3, c2 and c1 coefficients
        // are zero and the c0 coefficient is lexicographically largest,
        // then this element is lexicographically largest.

        self.c5.lexicographically_largest()
            | (self.c5.is_zero() & self.c4.lexicographically_largest())
            | (self.c4.is_zero() & self.c3.lexicographically_largest())
            | (self.c3.is_zero() & self.c2.lexicographically_largest())
            | (self.c2.is_zero() & self.c1.lexicographically_largest())
            | (self.c1.is_zero() & self.c0.lexicographically_largest())
    }

    #[inline]
    /// Computes the multiplication of an Fp6 element with an Fp element
    pub const fn mul_by_fp(&self, other: &Fp) -> Fp6 {
        let c0 = (&self.c0).mul(other);
        let c1 = (&self.c1).mul(other);
        let c2 = (&self.c2).mul(other);
        let c3 = (&self.c3).mul(other);
        let c4 = (&self.c4).mul(other);
        let c5 = (&self.c5).mul(other);

        Self {
            c0,
            c1,
            c2,
            c3,
            c4,
            c5,
        }
    }

    #[inline]
    /// Computes the multiplication of two Fp6 elements
    pub const fn mul(&self, other: &Fp6) -> Fp6 {
        // All helper values computed below are seen as u128 after modular reduction.
        // This allows for a faster computation of the result coordinates,
        // by computing all operations in /ZZ, and then finally converting the
        // values back to Fp through a modular reduction.
        //
        // The reduction uses `reduce_u96()` as all final values are less than 96 bits.
        //
        // To speed-up sextic extension multiplication, we do not use the direct sextic
        // representation and schoolbook method. Instead, we switch to a towered
        // representation, as a quadratic extension over a cubic one. The conversion
        // between these representations is costless.
        // Indeed, an element c = c_0 + c_1.Z + + c_2.Z^2 + c_3.Z^3 + c_4.Z^4 + c_5.Z^5 in
        // GF(p^6), defined as GF(p)[Z]/[Z^6 - δ] (with δ = 7 here) can be represented as
        // the element a = (a_00 + a_01.X + a_02.X^2) + (a_10 + a_11.X + a_12.X^2).Y where
        // the cubic extension GF(p^3) is defined as Fp[X]/(X3 − β). β = δ^2 is a cubic
        // non-residue in GF(p), and the sextic extension is defined as GF(p^3)[Y]/(Y^2 − γ),
        // with γ = δ a quadratic non-residue. The correspondence is then:
        //      a_00 = c_0
        //      a_01 = c_2
        //      a_02 = c_4
        //      a_10 = c_1
        //      a_11 = c_3
        //      a_12 = c_5
        //
        // We then use the Karatsuba method for arithmetic on the cubic extension of GF(p),
        // and the Schoolbook method for the final quadratic extension.
        //
        // Multiplication of GF(p^3) elements by γ consists of a single multiplication
        // by β and a permutation of coordinates.
        // Indeed, for a = a_0 + a_1.X + a_2.X^2, we have a*γ = a_2.β + a_0.X + a_1.X^2.
        //
        // With a = a_0 + a_1.Y and b = b_0 + b_1.Y, we compute the product c = ab with:
        // c_0 = a_0*b_0 + a_1*b_1.γ;
        // c_1 = a_0*b_1 + a_1*b_0;

        // Precomputations for a_0 * b_0
        let t00 = (&self.c0).mul(&other.c0).0 as u128;
        let t01 = (&self.c2).mul(&other.c2).0 as u128;
        let t02 = (&self.c4).mul(&other.c4).0 as u128;

        let s012 = (&self.c2).add(&self.c4);
        let tmp = (&other.c2).add(&other.c4);
        let s012 = (&s012).mul(&tmp).0 as u128;

        let s001 = (&self.c0).add(&self.c2);
        let tmp = (&other.c0).add(&other.c2);
        let s001 = (&s001).mul(&tmp).0 as u128;

        let s002 = (&self.c0).add(&self.c4);
        let tmp = (&other.c0).add(&other.c4);
        let s002 = (&s002).mul(&tmp).0 as u128;

        // Precomputations for a_1 * b_1
        let t10 = (&self.c1).mul(&other.c1).0 as u128;
        let t11 = (&self.c3).mul(&other.c3).0 as u128;
        let t12 = (&self.c5).mul(&other.c5).0 as u128;

        let s112 = (&self.c3).add(&self.c5);
        let tmp = (&other.c3).add(&other.c5);
        let s112 = (&s112).mul(&tmp).0 as u128;

        let s101 = (&self.c1).add(&self.c3);
        let tmp = (&other.c1).add(&other.c3);
        let s101 = (&s101).mul(&tmp).0 as u128;

        let s102 = (&self.c1).add(&self.c5);
        let tmp = (&other.c1).add(&other.c5);
        let s102 = (&s102).mul(&tmp).0 as u128;

        // Precomputations for a_0 * b_1
        let t20 = (&self.c0).mul(&other.c1).0 as u128;
        let t21 = (&self.c2).mul(&other.c3).0 as u128;
        let t22 = (&self.c4).mul(&other.c5).0 as u128;

        let s212 = (&self.c2).add(&self.c4);
        let tmp = (&other.c3).add(&other.c5);
        let s212 = (&s212).mul(&tmp).0 as u128;

        let s201 = (&self.c0).add(&self.c2);
        let tmp = (&other.c1).add(&other.c3);
        let s201 = (&s201).mul(&tmp).0 as u128;

        let s202 = (&self.c0).add(&self.c4);
        let tmp = (&other.c1).add(&other.c5);
        let s202 = (&s202).mul(&tmp).0 as u128;

        // Precomputations for a_1 * b_0
        let t30 = (&self.c1).mul(&other.c0).0 as u128;
        let t31 = (&self.c3).mul(&other.c2).0 as u128;
        let t32 = (&self.c5).mul(&other.c4).0 as u128;

        let s312 = (&self.c3).add(&self.c5);
        let tmp = (&other.c2).add(&other.c4);
        let s312 = (&s312).mul(&tmp).0 as u128;

        let s301 = (&self.c1).add(&self.c3);
        let tmp = (&other.c0).add(&other.c2);
        let s301 = (&s301).mul(&tmp).0 as u128;

        let s302 = (&self.c1).add(&self.c5);
        let tmp = (&other.c0).add(&other.c4);
        let s302 = (&s302).mul(&tmp).0 as u128;

        // c_0 = a_0*b_0 + a_1*b_1.γ;
        // c_1 = a_0*b_1 + a_1*b_0;

        // Compute a_0 * b_0
        let d00 = t01 + t02;
        let d00 = s012 + 0x1fffffffe00000002 - d00;
        let d00 = d00 * BETA;
        let d00 = d00 + t00;

        let d01 = t02 * BETA;
        let tmp = t00 + t01;
        let d01 = d01 + 0x1fffffffe00000002 - tmp;
        let d01 = d01 + s001;

        let d02 = t00 + t02;
        let d02 = t01 + 0x1fffffffe00000002 - d02;
        let d02 = d02 + s002;

        // Compute a_1 * b_1
        let d10 = t11 + t12;
        let d10 = s112 + 0x1fffffffe00000002 - d10;
        let d10 = d10 * BETA;
        let d10 = d10 + t10;

        let d11 = t12 * BETA;
        let tmp = t10 + t11;
        let d11 = d11 + 0x1fffffffe00000002 - tmp;
        let d11 = d11 + s101;

        let d12 = t10 + t12;
        let d12 = t11 + 0x1fffffffe00000002 - d12;
        let d12 = d12 + s102;

        // Compute a_0 * b_1
        let d20 = t21 + t22;
        let d20 = s212 + 0x1fffffffe00000002 - d20;
        let d20 = d20 * BETA;
        let d20 = d20 + t20;

        let d21 = t22 * BETA;
        let tmp = t20 + t21;
        let d21 = d21 + 0x1fffffffe00000002 - tmp;
        let d21 = d21 + s201;

        let d22 = t20 + t22;
        let d22 = t21 + 0x1fffffffe00000002 - d22;
        let d22 = d22 + s202;

        // Compute a_1 * b_0
        let d30 = t31 + t32;
        let d30 = s312 + 0x1fffffffe00000002 - d30;
        let d30 = d30 * BETA;
        let d30 = d30 + t30;

        let d31 = t32 * BETA;
        let tmp = t30 + t31;
        let d31 = d31 + 0x1fffffffe00000002 - tmp;
        let d31 = d31 + s301;

        let d32 = t30 + t32;
        let d32 = t31 + 0x1fffffffe00000002 - d32;
        let d32 = d32 + s302;

        // Compute the final coordinates, reduced by the modulus
        let c0 = Fp(reduce_u96(d00 + d12 * BETA));
        let c2 = Fp(reduce_u96(d01 + d10));
        let c4 = Fp(reduce_u96(d02 + d11));
        let c1 = Fp(reduce_u96(d20 + d30));
        let c3 = Fp(reduce_u96(d21 + d31));
        let c5 = Fp(reduce_u96(d22 + d32));

        Self {
            c0,
            c1,
            c2,
            c3,
            c4,
            c5,
        }
    }

    /// Computes the square of an Fp6 element
    #[inline]
    pub const fn square(&self) -> Self {
        // All helper values computed below are seen as u128 after modular reduction.
        // This allows for a faster computation of the result coordinates,
        // by computing all operations in /ZZ, and then finally converting the
        // values back to Fp through a modular reduction.
        //
        // The reduction uses `reduce_u96()` as all final values are less than 96 bits.

        let aa = (&self.c0).square().0 as u128;
        let ab = (&self.c0).mul(&self.c1).0 as u128;
        let ac = (&self.c0).mul(&self.c2).0 as u128;
        let ad = (&self.c0).mul(&self.c3).0 as u128;
        let ae = (&self.c0).mul(&self.c4).0 as u128;
        let af = (&self.c0).mul(&self.c5).0 as u128;

        let bb = (&self.c1).square().0 as u128;
        let bc = (&self.c1).mul(&self.c2).0 as u128;
        let bd = (&self.c1).mul(&self.c3).0 as u128;
        let be = (&self.c1).mul(&self.c4).0 as u128;
        let bf = (&self.c1).mul(&self.c5).0 as u128;

        let cc = (&self.c2).square().0 as u128;
        let cd = (&self.c2).mul(&self.c3).0 as u128;
        let ce = (&self.c2).mul(&self.c4).0 as u128;
        let cf = (&self.c2).mul(&self.c5).0 as u128;

        let dd = (&self.c3).square().0 as u128;
        let de = (&self.c3).mul(&self.c4).0 as u128;
        let df = (&self.c3).mul(&self.c5).0 as u128;

        let ee = (&self.c4).square().0 as u128;
        let ef = (&self.c4).mul(&self.c5).0 as u128;

        let ff = (&self.c5).square().0 as u128;

        let c0 = bf + ce;
        let c0 = c0 << 1;
        let c0 = c0 + dd;
        let c0 = c0 * BETA;
        let c0 = c0 + aa;
        let c0 = Fp(reduce_u96(c0));

        let c1 = cf + de;
        let c1 = c1 * BETA;
        let c1 = c1 + ab;
        let c1 = c1 << 1;
        let c1 = Fp(reduce_u96(c1));

        let c2 = df << 1;
        let c2 = c2 + ee;
        let c2 = c2 * BETA;
        let t2 = ac << 1;
        let c2 = c2 + t2;
        let c2 = c2 + bb;
        let c2 = Fp(reduce_u96(c2));

        let c3 = ef * BETA;
        let c3 = c3 + ad;
        let c3 = c3 + bc;
        let c3 = c3 << 1;
        let c3 = Fp(reduce_u96(c3));

        let t4 = ff * BETA;
        let c4 = ae + bd;
        let c4 = c4 << 1;
        let c4 = c4 + cc;
        let c4 = c4 + t4;
        let c4 = Fp(reduce_u96(c4));

        let c5 = af + be;
        let c5 = c5 + cd;
        let c5 = c5 << 1;
        let c5 = Fp(reduce_u96(c5));

        Self {
            c0,
            c1,
            c2,
            c3,
            c4,
            c5,
        }
    }

    /// Computes the square root of this element, if it exists.
    pub fn sqrt(&self) -> CtOption<Self> {
        // Tonelli-Shank's algorithm for q mod 16 = 1
        // See https://eprint.iacr.org/2020/1497.pdf, page 3 for a
        // constant time specification of the algorithm.

        // Compute the progenitor y of self
        // y = self^((t - 1) // 2)
        //   = self^0x3ffffffe800000053ffffff3800000167fffffe0800000233fffffe0800000167ffffff3800000053ffffffe
        let y = self.exp_vartime(&[
            0x800000053ffffffe,
            0x800000167ffffff3,
            0x800000233fffffe0,
            0x800000167fffffe0,
            0x800000053ffffff3,
            0x000000003ffffffe,
        ]);

        let mut s = self * y;
        let mut t = s * y;

        let mut z = TWO_ADIC_ROOT_OF_UNITY_P6;

        for k in (2..=TWO_ADICITY_P6).rev() {
            let mut b = t;

            square_assign_multi(&mut b, (k - 2) as usize);

            let new_s = s * z;
            s = Self::conditional_select(&new_s, &s, b.ct_eq(&Self::one()));
            z = z.square();
            let new_t = t * z;
            t = Self::conditional_select(&new_t, &t, b.ct_eq(&Self::one()));
        }

        CtOption::new(s, (s * s).ct_eq(self))
    }

    /// Computes the double of a field element
    #[inline]
    pub const fn double(&self) -> Self {
        Self {
            c0: (&self.c0).double(),
            c1: (&self.c1).double(),
            c2: (&self.c2).double(),
            c3: (&self.c3).double(),
            c4: (&self.c4).double(),
            c5: (&self.c5).double(),
        }
    }

    /// Computes the triple of a field element
    #[inline]
    pub const fn triple(&self) -> Self {
        Self {
            c0: (&self.c0).triple(),
            c1: (&self.c1).triple(),
            c2: (&self.c2).triple(),
            c3: (&self.c3).triple(),
            c4: (&self.c4).triple(),
            c5: (&self.c5).triple(),
        }
    }

    /// Computes the summation of two field elements
    #[inline]
    pub const fn add(&self, rhs: &Self) -> Self {
        Self {
            c0: (&self.c0).add(&rhs.c0),
            c1: (&self.c1).add(&rhs.c1),
            c2: (&self.c2).add(&rhs.c2),
            c3: (&self.c3).add(&rhs.c3),
            c4: (&self.c4).add(&rhs.c4),
            c5: (&self.c5).add(&rhs.c5),
        }
    }

    /// Computes the difference of two field elements
    #[inline]
    pub const fn sub(&self, rhs: &Self) -> Self {
        Self {
            c0: (&self.c0).sub(&rhs.c0),
            c1: (&self.c1).sub(&rhs.c1),
            c2: (&self.c2).sub(&rhs.c2),
            c3: (&self.c3).sub(&rhs.c3),
            c4: (&self.c4).sub(&rhs.c4),
            c5: (&self.c5).sub(&rhs.c5),
        }
    }

    /// Computes the negation of a field element
    #[inline]
    pub const fn neg(&self) -> Self {
        Self {
            c0: (&self.c0).neg(),
            c1: (&self.c1).neg(),
            c2: (&self.c2).neg(),
            c3: (&self.c3).neg(),
            c4: (&self.c4).neg(),
            c5: (&self.c5).neg(),
        }
    }

    /// Computes the multiplicative inverse of this field
    /// element, returning None in the case that this element
    /// is zero.
    #[inline]
    pub fn invert(&self) -> CtOption<Self> {
        // Adapted from "A Fast Algorithm for Computing Multiplicative
        // Inverses in GF(2m) Using Normal Bases" from Itoh and Tsujii.

        let t0 = self.frobenius();
        let t1 = self.frobenius_triple();
        let t2 = (&t0).mul(&t1);
        let t3 = t2.frobenius();
        let t4 = t2.frobenius_double();
        let t5 = (&t4).mul(&t3);
        let t5 = (&t5).mul(&t0);
        let inv = (&t5).mul(self).c0;

        inv.invert().map(|t| t5.mul_by_fp(&t))
    }

    /// Exponentiates `self` by `power`, where `power` is a
    /// little-endian order integer exponent.
    pub fn exp(self, by: &[u64; 6]) -> Self {
        let mut res = Self::one();
        for e in by.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();
                let mut tmp = res;
                tmp *= self;
                res.conditional_assign(&tmp, (((e >> i) & 1) as u8).into());
            }
        }
        res
    }

    /// Although this is labeled "vartime", it is only
    /// variable time with respect to the exponent.
    pub fn exp_vartime(&self, by: &[u64; 6]) -> Self {
        let mut res = Self::one();
        for e in by.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();

                if ((*e >> i) & 1) == 1 {
                    res *= self;
                }
            }
        }
        res
    }

    /// Computes the Frobenius endomorphism
    #[inline]
    pub const fn frobenius(&self) -> Self {
        Self {
            c0: self.c0,
            c1: (&Fp(18446744065119617026)).mul(&self.c1),
            c2: (&Fp(18446744065119617025)).mul(&self.c2),
            c3: (&Fp(18446744069414584320)).mul(&self.c3),
            c4: (&Fp(4294967295)).mul(&self.c4),
            c5: (&Fp(4294967296)).mul(&self.c5),
        }
    }

    /// Computes the Frobenius endomorphism twice
    #[inline]
    pub const fn frobenius_double(&self) -> Self {
        Self {
            c0: self.c0,
            c1: (&Fp(18446744065119617025)).mul(&self.c1),
            c2: (&Fp(4294967295)).mul(&self.c2),
            c3: self.c3,
            c4: (&Fp(18446744065119617025)).mul(&self.c4),
            c5: (&Fp(4294967295)).mul(&self.c5),
        }
    }

    /// Computes the Frobenius endomorphism thrice
    #[inline]
    pub const fn frobenius_triple(&self) -> Self {
        Self {
            c0: self.c0,
            c1: (&self.c1).neg(),
            c2: self.c2,
            c3: (&self.c3).neg(),
            c4: self.c4,
            c5: (&self.c5).neg(),
        }
    }

    /// Outputs the internal representation as 6 64-bit limbs after canonical reduction
    pub const fn output_internal(&self) -> [u64; 6] {
        [
            self.c0.output_internal(),
            self.c1.output_internal(),
            self.c2.output_internal(),
            self.c3.output_internal(),
            self.c4.output_internal(),
            self.c5.output_internal(),
        ]
    }

    /// Converts an `Fp6` element into a byte representation in
    /// little-endian byte order.
    pub fn to_bytes(&self) -> [u8; 48] {
        let mut bytes = [0u8; 48];

        bytes[0..8].copy_from_slice(&self.c0.to_bytes());
        bytes[8..16].copy_from_slice(&self.c1.to_bytes());
        bytes[16..24].copy_from_slice(&self.c2.to_bytes());
        bytes[24..32].copy_from_slice(&self.c3.to_bytes());
        bytes[32..40].copy_from_slice(&self.c4.to_bytes());
        bytes[40..48].copy_from_slice(&self.c5.to_bytes());

        bytes
    }

    /// Attempts to convert a little-endian byte representation of
    /// a scalar into a `Fp6` element, failing if the input is not canonical.
    pub fn from_bytes(bytes: &[u8; 48]) -> CtOption<Self> {
        let mut array = [0u8; 8];

        array.copy_from_slice(&bytes[0..8]);
        let c0 = Fp::from_bytes(&array);

        array.copy_from_slice(&bytes[8..16]);
        let c1 = Fp::from_bytes(&array);

        array.copy_from_slice(&bytes[16..24]);
        let c2 = Fp::from_bytes(&array);

        array.copy_from_slice(&bytes[24..32]);
        let c3 = Fp::from_bytes(&array);

        array.copy_from_slice(&bytes[32..40]);
        let c4 = Fp::from_bytes(&array);

        array.copy_from_slice(&bytes[40..48]);
        let c5 = Fp::from_bytes(&array);

        let is_some =
            c0.is_some() & c1.is_some() & c2.is_some() & c3.is_some() & c4.is_some() & c5.is_some();

        CtOption::new(
            Self {
                c0: c0.unwrap_or(Fp::zero()),
                c1: c1.unwrap_or(Fp::zero()),
                c2: c2.unwrap_or(Fp::zero()),
                c3: c3.unwrap_or(Fp::zero()),
                c4: c4.unwrap_or(Fp::zero()),
                c5: c5.unwrap_or(Fp::zero()),
            },
            is_some,
        )
    }

    /// Constructs an element of `Fp6` without checking that it is
    /// canonical.
    pub const fn from_raw_unchecked(value: [u64; 6]) -> Self {
        Self {
            c0: Fp::from_raw_unchecked(value[0]),
            c1: Fp::from_raw_unchecked(value[1]),
            c2: Fp::from_raw_unchecked(value[2]),
            c3: Fp::from_raw_unchecked(value[3]),
            c4: Fp::from_raw_unchecked(value[4]),
            c5: Fp::from_raw_unchecked(value[5]),
        }
    }
}

// FIELD TRAITS IMPLEMENTATION
// ================================================================================================

impl Field for Fp6 {
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

// SERDE SERIALIZATION
// ================================================================================================

#[cfg(feature = "serialize")]
impl Serialize for Fp6 {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        use serde::ser::SerializeTuple;
        let mut tup = serializer.serialize_tuple(48)?;
        for byte in self.to_bytes().iter() {
            tup.serialize_element(byte)?;
        }
        tup.end()
    }
}

#[cfg(feature = "serialize")]
impl<'de> Deserialize<'de> for Fp6 {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct Fp6Visitor;

        impl<'de> Visitor<'de> for Fp6Visitor {
            type Value = Fp6;

            fn expecting(&self, formatter: &mut Formatter) -> fmt::Result {
                formatter.write_str("a valid field element")
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<Fp6, A::Error>
            where
                A: serde::de::SeqAccess<'de>,
            {
                let mut bytes = [0u8; 48];
                for (i, byte) in bytes.iter_mut().enumerate() {
                    *byte = seq
                        .next_element()?
                        .ok_or_else(|| serde::de::Error::invalid_length(i, &"expected 48 bytes"))?;
                }
                let elem = Fp6::from_bytes(&bytes);
                if bool::from(elem.is_none()) {
                    Err(serde::de::Error::custom("decompression failed"))
                } else {
                    Ok(elem.unwrap())
                }
            }
        }

        deserializer.deserialize_tuple(48, Fp6Visitor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::OsRng;

    // DISPLAY
    // ================================================================================================

    #[test]
    fn test_debug() {
        assert_eq!(
            format!("{:?}", Fp6::zero()),
            "0 + 0*u + 0*u^2 + 0*u^3 + 0*u^4 + 0*u^5"
        );
        assert_eq!(
            format!("{:?}", Fp6::one()),
            "1 + 0*u + 0*u^2 + 0*u^3 + 0*u^4 + 0*u^5"
        );
        assert_eq!(
            format!("{:?}", Fp6::new([1, 2, 3, 4, 5, 6])),
            "1 + 2*u + 3*u^2 + 4*u^3 + 5*u^4 + 6*u^5"
        );

        let a = Fp6::one().neg();
        assert_eq!(
            format!("{:?}", a),
            "18446744069414584320 + 0*u + 0*u^2 + 0*u^3 + 0*u^4 + 0*u^5"
        );
    }

    #[test]
    fn test_output_reduced_limbs() {
        assert_eq!(
            format!("{:?}", Fp6::zero().output_internal()),
            "[0, 0, 0, 0, 0, 0]"
        );
        assert_eq!(
            format!("{:?}", Fp6::one().output_internal()),
            "[1, 0, 0, 0, 0, 0]"
        );
        let a = Fp6::one().neg();
        assert_eq!(
            format!("{:?}", a.output_internal()),
            "[18446744069414584320, 0, 0, 0, 0, 0]"
        );
    }

    // BASIC ALGEBRA
    // ================================================================================================

    #[test]
    fn test_conditional_selection() {
        let a = Fp6 {
            c0: Fp::one(),
            c1: Fp::new(2),
            c2: Fp::new(3),
            c3: Fp::new(4),
            c4: Fp::new(5),
            c5: Fp::new(6),
        };
        let b = Fp6 {
            c0: Fp::new(7),
            c1: Fp::new(8),
            c2: Fp::new(9),
            c3: Fp::new(10),
            c4: Fp::new(11),
            c5: Fp::new(12),
        };

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
        fn is_equal(a: &Fp6, b: &Fp6) -> bool {
            let eq = a == b;
            let ct_eq = a.ct_eq(b);

            assert_eq!(eq, bool::from(ct_eq));

            eq
        }

        assert!(is_equal(
            &Fp6 {
                c0: Fp::one(),
                c1: Fp::new(2),
                c2: Fp::new(3),
                c3: Fp::new(4),
                c4: Fp::new(5),
                c5: Fp::new(6),
            },
            &Fp6 {
                c0: Fp::one(),
                c1: Fp::new(2),
                c2: Fp::new(3),
                c3: Fp::new(4),
                c4: Fp::new(5),
                c5: Fp::new(6),
            }
        ));

        assert!(!is_equal(
            &Fp6 {
                c0: Fp::new(2),
                c1: Fp::new(2),
                c2: Fp::new(3),
                c3: Fp::new(4),
                c4: Fp::new(5),
                c5: Fp::new(6),
            },
            &Fp6 {
                c0: Fp::one(),
                c1: Fp::new(2),
                c2: Fp::new(3),
                c3: Fp::new(4),
                c4: Fp::new(5),
                c5: Fp::new(6),
            }
        ));

        assert!(!is_equal(
            &Fp6 {
                c0: Fp::one(),
                c1: Fp::new(3),
                c2: Fp::new(3),
                c3: Fp::new(4),
                c4: Fp::new(5),
                c5: Fp::new(6),
            },
            &Fp6 {
                c0: Fp::one(),
                c1: Fp::new(2),
                c2: Fp::new(3),
                c3: Fp::new(4),
                c4: Fp::new(5),
                c5: Fp::new(6),
            }
        ));

        assert!(!is_equal(
            &Fp6 {
                c0: Fp::one(),
                c1: Fp::new(2),
                c2: Fp::new(4),
                c3: Fp::new(4),
                c4: Fp::new(5),
                c5: Fp::new(6),
            },
            &Fp6 {
                c0: Fp::one(),
                c1: Fp::new(2),
                c2: Fp::new(3),
                c3: Fp::new(4),
                c4: Fp::new(5),
                c5: Fp::new(6),
            }
        ));

        assert!(!is_equal(
            &Fp6 {
                c0: Fp::one(),
                c1: Fp::new(2),
                c2: Fp::new(3),
                c3: Fp::new(5),
                c4: Fp::new(5),
                c5: Fp::new(6),
            },
            &Fp6 {
                c0: Fp::one(),
                c1: Fp::new(2),
                c2: Fp::new(3),
                c3: Fp::new(4),
                c4: Fp::new(5),
                c5: Fp::new(6),
            }
        ));

        assert!(!is_equal(
            &Fp6 {
                c0: Fp::one(),
                c1: Fp::new(2),
                c2: Fp::new(3),
                c3: Fp::new(4),
                c4: Fp::new(6),
                c5: Fp::new(6),
            },
            &Fp6 {
                c0: Fp::one(),
                c1: Fp::new(2),
                c2: Fp::new(3),
                c3: Fp::new(4),
                c4: Fp::new(5),
                c5: Fp::new(6),
            }
        ));

        assert!(!is_equal(
            &Fp6 {
                c0: Fp::one(),
                c1: Fp::new(2),
                c2: Fp::new(3),
                c3: Fp::new(4),
                c4: Fp::new(5),
                c5: Fp::new(7),
            },
            &Fp6 {
                c0: Fp::one(),
                c1: Fp::new(2),
                c2: Fp::new(3),
                c3: Fp::new(4),
                c4: Fp::new(5),
                c5: Fp::new(6),
            }
        ));

        assert!(bool::from(Fp6::default().is_zero()));
        assert!(!bool::from(Fp6::zero().ct_eq(&Fp6::one())));

        assert_eq!(Fp6::zero(), Fp6::new([0, 0, 0, 0, 0, 0]));
        assert_eq!(Fp6::one(), Fp6::new([1, 0, 0, 0, 0, 0]));
    }

    #[test]
    fn test_squaring() {
        let a = Fp6 {
            c0: Fp::new(7311251207099623982),
            c1: Fp::new(11101156575191577169),
            c2: Fp::new(3471450285948349450),
            c3: Fp::new(18081717484378824777),
            c4: Fp::new(17271827452459074620),
            c5: Fp::new(13245881484175979175),
        };
        let b = Fp6 {
            c0: Fp::new(18442855563257030637),
            c1: Fp::new(416354386888665317),
            c2: Fp::new(12488382742185552506),
            c3: Fp::new(6173432951622501799),
            c4: Fp::new(9922025629450424125),
            c5: Fp::new(9218000164529386038),
        };

        assert_eq!(a.square(), b);
    }

    #[test]
    fn test_sqrt() {
        for _ in 0..100 {
            let a = Fp6::random(&mut OsRng).square();
            let b = a.sqrt().unwrap();
            assert_eq!(a, b.square());
        }

        assert_eq!(Fp6::zero().sqrt().unwrap(), Fp6::zero());
        assert_eq!(Fp6::one().sqrt().unwrap(), Fp6::one());

        // u + 2
        // is not a quadratic residue in Fp6
        assert!(bool::from(
            Fp6 {
                c0: Fp::new(2),
                c1: Fp::one(),
                c2: Fp::zero(),
                c3: Fp::zero(),
                c4: Fp::zero(),
                c5: Fp::zero(),
            }
            .sqrt()
            .is_none()
        ));
    }

    #[test]
    fn test_multiplication() {
        let a = Fp6 {
            c0: Fp::one(),
            c1: Fp::new(2),
            c2: Fp::new(3),
            c3: Fp::new(4),
            c4: Fp::new(5),
            c5: Fp::new(6),
        };
        let b = Fp6::one();
        let c = Fp6 {
            c0: Fp::one(),
            c1: Fp::new(2),
            c2: Fp::new(3),
            c3: Fp::new(4),
            c4: Fp::new(5),
            c5: Fp::new(6),
        };

        assert_eq!(a * b, c);

        let a = Fp6 {
            c0: Fp::new(7311251207099623982),
            c1: Fp::new(11101156575191577169),
            c2: Fp::new(3471450285948349450),
            c3: Fp::new(18081717484378824777),
            c4: Fp::new(17271827452459074620),
            c5: Fp::new(13245881484175979175),
        };
        let b = Fp6 {
            c0: Fp::new(10376570592177434053),
            c1: Fp::new(2266730037751779421),
            c2: Fp::new(17604547680137376634),
            c3: Fp::new(8522386576290610241),
            c4: Fp::new(4622828580943454349),
            c5: Fp::new(3310945140500303981),
        };
        let c = Fp6 {
            c0: Fp::new(16699204345089849089),
            c1: Fp::new(12290068205868725576),
            c2: Fp::new(15503892176250492698),
            c3: Fp::new(16153850030910369456),
            c4: Fp::new(4357822170351527454),
            c5: Fp::new(15609745236686299732),
        };

        assert_eq!(a * b, c);
    }

    #[test]
    fn test_addition() {
        let a = Fp6 {
            c0: Fp::one(),
            c1: Fp::new(2),
            c2: Fp::new(3),
            c3: Fp::new(4),
            c4: Fp::new(5),
            c5: Fp::new(6),
        };
        let b = Fp6 {
            c0: Fp::new(6),
            c1: Fp::new(5),
            c2: Fp::new(4),
            c3: Fp::new(3),
            c4: Fp::new(2),
            c5: Fp::one(),
        };
        let c = Fp6 {
            c0: Fp::new(7),
            c1: Fp::new(7),
            c2: Fp::new(7),
            c3: Fp::new(7),
            c4: Fp::new(7),
            c5: Fp::new(7),
        };

        assert_eq!(a + b, c);
    }

    #[test]
    fn test_subtraction() {
        let a = Fp6 {
            c0: Fp::new(6),
            c1: Fp::new(5),
            c2: Fp::new(4),
            c3: Fp::new(3),
            c4: Fp::new(2),
            c5: Fp::one(),
        };
        let b = Fp6 {
            c0: Fp::new(3),
            c1: Fp::new(3),
            c2: Fp::new(2),
            c3: Fp::new(2),
            c4: Fp::one(),
            c5: Fp::one(),
        };
        let c = Fp6 {
            c0: Fp::new(3),
            c1: Fp::new(2),
            c2: Fp::new(2),
            c3: Fp::one(),
            c4: Fp::one(),
            c5: Fp::zero(),
        };

        assert_eq!(a - b, c);
    }

    #[test]
    fn test_negation() {
        let mut rng = OsRng;

        let a = Fp6 {
            c0: Fp::one(),
            c1: Fp::new(2),
            c2: Fp::new(3),
            c3: Fp::new(4),
            c4: Fp::new(5),
            c5: Fp::new(6),
        };
        let b = Fp6 {
            c0: -Fp::one(),
            c1: -Fp::new(2),
            c2: -Fp::new(3),
            c3: -Fp::new(4),
            c4: -Fp::new(5),
            c5: -Fp::new(6),
        };

        assert_eq!(-a, b);

        for _ in 0..100 {
            let a = Fp6::random(&mut rng);
            let b = -a;

            assert_eq!(a + b, Fp6::zero());
        }
    }

    #[test]
    fn test_inversion() {
        let a = Fp6 {
            c0: Fp::new(5147025334779536289),
            c1: Fp::new(3732316503485825522),
            c2: Fp::new(9593647204665260945),
            c3: Fp::new(1415511037148678607),
            c4: Fp::new(15711614927439708015),
            c5: Fp::new(1050791595125110226),
        };

        let b = Fp6 {
            c0: Fp::new(6272010161093965568),
            c1: Fp::new(18017385922404754760),
            c2: Fp::new(8403606426326570881),
            c3: Fp::new(12472658335656979833),
            c4: Fp::new(14564178836325739191),
            c5: Fp::new(7388120411391452745),
        };

        assert_eq!(a.invert().unwrap(), b);

        assert_eq!(Fp6::one().invert().unwrap(), Fp6::one());

        assert!(bool::from(Fp6::zero().invert().is_none()));
    }

    #[test]
    fn test_invert_is_pow() {
        let mut rng = OsRng;

        let p6_minus_2 = [
            0xfffffff9ffffffff,
            0xffffffce00000014,
            0xffffff8200000059,
            0xffffff820000008c,
            0xffffffce00000059,
            0xfffffffa00000014,
        ];

        let mut r1 = Fp6::random(&mut rng);
        let mut r2 = r1;
        let mut r3 = r2;

        for _ in 0..100 {
            r1 = r1.invert().unwrap();
            r2 = r2.exp(&p6_minus_2);
            r3 = r3.exp_vartime(&p6_minus_2);

            assert_eq!(r1, r2);
            assert_eq!(r2, r3);

            // Call random() so we check a different element each time
            r1 = Fp6::random(&mut rng);
            r2 = r1;
            r3 = r1;
        }
    }

    #[test]
    fn test_frobenius() {
        let mut rng = OsRng;

        for _ in 0..100 {
            let a = Fp6::random(&mut rng);

            assert_eq!(a.frobenius(), a.exp(&[0xffffffff00000001, 0, 0, 0, 0, 0]));
            assert_eq!(a.frobenius_double(), a.frobenius().frobenius());
            assert_eq!(a.frobenius_triple(), a.frobenius_double().frobenius());
        }
    }

    // ROOTS OF UNITY
    // ================================================================================================

    #[test]
    fn test_get_root_of_unity() {
        let two_pow_33 = 1 << TWO_ADICITY_P6 as u64;
        assert_eq!(
            Fp6::one(),
            TWO_ADIC_ROOT_OF_UNITY_P6.exp(&[two_pow_33, 0, 0, 0, 0, 0])
        );
        assert_ne!(
            Fp6::one(),
            TWO_ADIC_ROOT_OF_UNITY_P6.exp(&[two_pow_33 - 1, 0, 0, 0, 0, 0])
        );
    }

    #[test]
    fn test_lexicographic_largest() {
        assert!(!bool::from(Fp6::zero().lexicographically_largest()));
        assert!(!bool::from(Fp6::one().lexicographically_largest()));
        // a = 13413783073807968269*u^5 + 15191285573968924193*u^4 + 6249462481928133327*u^3
        //      + 12517261035513395024*u^2 + 17772767032074153840*u + 16347749849157988583
        let a = Fp6 {
            c0: Fp::new(16347749849157988583),
            c1: Fp::new(17772767032074153840),
            c2: Fp::new(12517261035513395024),
            c3: Fp::new(6249462481928133327),
            c4: Fp::new(15191285573968924193),
            c5: Fp::new(13413783073807968269),
        };

        // b = -a = 5032960995606616052*u^5 + 3255458495445660128*u^4 + 12197281587486450994*u^3
        //      + 5929483033901189297*u^2 + 673977037340430481*u + 2098994220256595738
        let b = Fp6 {
            c0: Fp::new(2098994220256595738),
            c1: Fp::new(673977037340430481),
            c2: Fp::new(5929483033901189297),
            c3: Fp::new(12197281587486450994),
            c4: Fp::new(3255458495445660128),
            c5: Fp::new(5032960995606616052),
        };

        assert_eq!(a.square(), b.square());
        assert!(bool::from(a.lexicographically_largest()));
        assert!(!bool::from(b.lexicographically_largest()));

        assert!(bool::from(
            Fp6 {
                c0: Fp::new(16164148524715685436),
                c1: Fp::new(13417117976106719244),
                c2: Fp::zero(),
                c3: Fp::zero(),
                c4: Fp::zero(),
                c5: Fp::zero(),
            }
            .lexicographically_largest()
        ));
    }

    #[test]
    fn test_zeroize() {
        use zeroize::Zeroize;

        let mut a = Fp6::one();
        a.zeroize();
        assert!(bool::from(a.is_zero()));
    }

    #[test]
    fn test_from_fp() {
        let mut rng = OsRng;
        let v = rng.next_u64();
        let e = Fp::new(v);

        let e_fp6 = Fp6::new([v, 0, 0, 0, 0, 0]);
        let array: [Fp; 6] = [
            e,
            Fp::zero(),
            Fp::zero(),
            Fp::zero(),
            Fp::zero(),
            Fp::zero(),
        ];

        assert_eq!(e_fp6, e.into());
        assert_eq!(Fp6::from(array), e_fp6);
    }

    // FIELD TRAIT
    // ================================================================================================

    #[test]
    fn test_field_trait_methods() {
        assert_eq!(<Fp6 as Field>::zero(), Fp6::new([0, 0, 0, 0, 0, 0]));
        assert_eq!(<Fp6 as Field>::one(), Fp6::new([1, 0, 0, 0, 0, 0]));

        assert_eq!(
            bool::from(<Fp6 as Field>::zero().is_zero()),
            bool::from(Fp6::new([0, 0, 0, 0, 0, 0]).is_zero())
        );
        assert_eq!(
            bool::from(<Fp6 as Field>::one().is_zero()),
            bool::from(Fp6::new([1, 0, 0, 0, 0, 0]).is_zero())
        );

        let mut rng = OsRng;
        let e = Fp6::random(&mut rng).square();

        assert_eq!(<Fp6 as Field>::square(&e), e.square());
        assert_eq!(<Fp6 as Field>::double(&e), e.double());

        assert_eq!(<Fp6 as Field>::invert(&e).unwrap(), e.invert().unwrap());
        assert!(bool::from(<Fp6 as Field>::invert(&Fp6::zero()).is_none()));

        assert_eq!(<Fp6 as Field>::sqrt(&e).unwrap(), e.sqrt().unwrap());
        assert!(bool::from(
            <Fp6 as Field>::sqrt(&Fp6 {
                c0: Fp::new(7727692126874014667),
                c1: Fp::new(13893222280919250633),
                c2: Fp::new(10956350845007515521),
                c3: Fp::new(281659085176887693),
                c4: Fp::new(6673350866139507059),
                c5: Fp::new(5956049640144253975),
            })
            .is_none()
        ));
    }

    // SERIALIZATION / DESERIALIZATION
    // ================================================================================================

    #[test]
    fn test_to_bytes() {
        assert_eq!(
            Fp6::zero().to_bytes(),
            [
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ]
        );

        assert_eq!(
            Fp6::one().to_bytes(),
            [
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ]
        );

        assert_eq!(
            (-&Fp6::one()).to_bytes(),
            [
                0, 0, 0, 0, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ]
        );
    }

    #[test]
    fn test_from_bytes() {
        let mut rng = OsRng;
        for _ in 0..100 {
            let a = Fp6::random(&mut rng);
            let bytes = a.to_bytes();
            assert_eq!(a, Fp6::from_bytes(&bytes).unwrap());
        }

        assert_eq!(
            Fp6::from_bytes(&[
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
            .unwrap(),
            Fp6::zero()
        );

        assert_eq!(
            Fp6::from_bytes(&[
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
            .unwrap(),
            Fp6::one()
        );

        // -1 should work
        assert_eq!(
            Fp6::from_bytes(&[
                0, 0, 0, 0, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
            .unwrap(),
            -Fp6::one()
        );

        // Anything equal or larger than M in one of the members is invalid
        assert!(bool::from(
            Fp6::from_bytes(&[
                1, 0, 0, 0, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
            .is_none()
        ));

        assert!(bool::from(
            Fp6::from_bytes(&[
                1, 0, 0, 0, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
            .is_none()
        ));

        assert!(bool::from(
            Fp6::from_bytes(&[
                0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
            .is_none()
        ));

        assert!(bool::from(
            Fp6::from_bytes(&[
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 255, 255, 255, 255, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
            .is_none()
        ));

        assert!(bool::from(
            Fp6::from_bytes(&[
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
            .is_none()
        ));

        assert!(bool::from(
            Fp6::from_bytes(&[
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 1, 0, 0, 0, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0
            ])
            .is_none()
        ));

        assert!(bool::from(
            Fp6::from_bytes(&[
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 255, 255, 255, 255
            ])
            .is_none()
        ));
    }

    #[test]
    #[cfg(feature = "serialize")]
    fn test_serde() {
        let mut rng = OsRng;
        let element = Fp6::random(&mut rng);
        let encoded = bincode::serialize(&element).unwrap();
        let parsed: Fp6 = bincode::deserialize(&encoded).unwrap();
        assert_eq!(parsed, element);

        // Check that the encoding is 48 bytes exactly
        assert_eq!(encoded.len(), 48);

        // Check that the encoding itself matches the usual one
        assert_eq!(element, bincode::deserialize(&element.to_bytes()).unwrap());

        // Check that invalid encodings fail
        let wrong_encoding = [255; 48];
        assert!(bincode::deserialize::<Fp6>(&wrong_encoding).is_err());

        let element = Fp6::random(&mut rng);
        let encoded = bincode::serialize(&element).unwrap();
        assert!(bincode::deserialize::<Fp6>(&encoded[0..47]).is_err());
    }
}
