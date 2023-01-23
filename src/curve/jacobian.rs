// Copyright (c) 2021-2023 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This module provides an implementation of the cheetah curve in
//! jacobian coordinates.

use core::{
    borrow::Borrow,
    cmp::Ordering,
    fmt,
    hash::{Hash, Hasher},
    iter::Sum,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use super::B;

use crate::{AffinePoint, CompressedPoint, UncompressedPoint};
use crate::{Fp, Fp6, Scalar};

use crate::{MINUS_SHIFT_POINT_ARRAY, SHIFT_POINT_MODIFIED_JACOBIAN};

use crate::constants::ODD_MULTIPLES_BASEPOINT;
use crate::LookupTable;
use crate::NafLookupTable;

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;
use group::{Curve, Group};
use rand_core::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "serialize")]
use serde::{self, Deserialize, Deserializer, Serialize, Serializer};

impl_binops_additive!(JacobianPoint, AffinePoint);
impl_binops_additive_specify_output!(AffinePoint, JacobianPoint, JacobianPoint);

impl<'a, 'b> Add<&'b AffinePoint> for &'a JacobianPoint {
    type Output = JacobianPoint;

    #[inline]
    fn add(self, rhs: &'b AffinePoint) -> JacobianPoint {
        self.add_mixed(rhs)
    }
}

impl<'a, 'b> Sub<&'b AffinePoint> for &'a JacobianPoint {
    type Output = JacobianPoint;

    #[inline]
    fn sub(self, rhs: &'b AffinePoint) -> JacobianPoint {
        self + (-rhs)
    }
}

impl<T> Sum<T> for JacobianPoint
where
    T: Borrow<JacobianPoint>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Self::identity(), |acc, item| acc + item.borrow())
    }
}

/// A jacobian point
#[derive(Copy, Clone, Debug)]
pub struct JacobianPoint {
    pub(crate) x: Fp6,
    pub(crate) y: Fp6,
    pub(crate) z: Fp6,
}

impl Default for JacobianPoint {
    fn default() -> JacobianPoint {
        JacobianPoint::identity()
    }
}

impl zeroize::DefaultIsZeroes for JacobianPoint {}

impl fmt::Display for JacobianPoint {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl Hash for JacobianPoint {
    fn hash<H: Hasher>(&self, hasher: &mut H) {
        AffinePoint::from(self).hash(hasher);
    }
}

impl<'a> From<&'a AffinePoint> for JacobianPoint {
    fn from(p: &'a AffinePoint) -> JacobianPoint {
        JacobianPoint {
            x: p.x,
            y: p.y,
            z: Fp6::conditional_select(&Fp6::one(), &Fp6::zero(), p.infinity),
        }
    }
}

impl From<AffinePoint> for JacobianPoint {
    fn from(p: AffinePoint) -> JacobianPoint {
        JacobianPoint::from(&p)
    }
}

impl<'a> From<&'a ModifiedJacobianPoint> for JacobianPoint {
    fn from(p: &'a ModifiedJacobianPoint) -> JacobianPoint {
        JacobianPoint {
            x: p.x,
            y: p.y,
            z: p.z,
        }
    }
}

impl From<ModifiedJacobianPoint> for JacobianPoint {
    fn from(p: ModifiedJacobianPoint) -> JacobianPoint {
        JacobianPoint::from(&p)
    }
}

impl ConstantTimeEq for JacobianPoint {
    fn ct_eq(&self, other: &Self) -> Choice {
        // Is (xz^2, yz^3, z) equal to (x'z'^2, y'z'^3, z') when converted to affine?

        let z1_squared = self.z.square();
        let z2_squared = other.z.square();

        let x1 = self.x * z2_squared;
        let x2 = other.x * z1_squared;

        let y1 = self.y * z2_squared * other.z;
        let y2 = other.y * z1_squared * self.z;

        let self_is_zero = self.z.is_zero();
        let other_is_zero = other.z.is_zero();

        (self_is_zero & other_is_zero) // Both point at infinity
            | ((!self_is_zero) & (!other_is_zero) & x1.ct_eq(&x2) & y1.ct_eq(&y2))
        // Neither point at infinity, coordinates are the same
    }
}

impl ConditionallySelectable for JacobianPoint {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        JacobianPoint {
            x: Fp6::conditional_select(&a.x, &b.x, choice),
            y: Fp6::conditional_select(&a.y, &b.y, choice),
            z: Fp6::conditional_select(&a.z, &b.z, choice),
        }
    }
}

impl Eq for JacobianPoint {}
impl PartialEq for JacobianPoint {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl<'a> Neg for &'a JacobianPoint {
    type Output = JacobianPoint;

    #[inline]
    fn neg(self) -> JacobianPoint {
        self.neg()
    }
}

impl Neg for JacobianPoint {
    type Output = JacobianPoint;

    #[inline]
    fn neg(self) -> JacobianPoint {
        -&self
    }
}

impl<'a, 'b> Add<&'b JacobianPoint> for &'a JacobianPoint {
    type Output = JacobianPoint;

    #[inline]
    fn add(self, rhs: &'b JacobianPoint) -> JacobianPoint {
        self.add(rhs)
    }
}

impl<'a, 'b> Sub<&'b JacobianPoint> for &'a JacobianPoint {
    type Output = JacobianPoint;

    #[inline]
    fn sub(self, rhs: &'b JacobianPoint) -> JacobianPoint {
        self + (-rhs)
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a JacobianPoint {
    type Output = JacobianPoint;

    fn mul(self, other: &'b Scalar) -> Self::Output {
        self.multiply(&other.to_bytes())
    }
}

impl_binops_additive!(JacobianPoint, JacobianPoint);
impl_binops_multiplicative!(JacobianPoint, Scalar);

impl JacobianPoint {
    /// Returns the identity of the group: the point at infinity.
    pub const fn identity() -> JacobianPoint {
        JacobianPoint {
            x: Fp6::one(),
            y: Fp6::one(),
            z: Fp6::zero(),
        }
    }

    /// Returns the x coordinate of this JacobianPoint
    pub const fn get_x(&self) -> Fp6 {
        self.x
    }

    /// Returns the y coordinate of this JacobianPoint
    pub const fn get_y(&self) -> Fp6 {
        self.y
    }

    /// Returns the z coordinate of this JacobianPoint
    pub const fn get_z(&self) -> Fp6 {
        self.z
    }

    /// Computes a random `JacobianPoint` element
    pub fn random(mut rng: impl RngCore) -> Self {
        loop {
            let x = Fp6::random(&mut rng);
            let flip_sign = rng.next_u32() % 2 != 0;

            // Obtain the corresponding y-coordinate given x as y = sqrt(x^3 + x + B)
            let p = ((x.square() * x) + x + B).sqrt().map(|y| AffinePoint {
                x,
                y: if flip_sign { -y } else { y },
                infinity: 0.into(),
            });

            if p.is_some().into() {
                let p = JacobianPoint::from(p.unwrap()).clear_cofactor();

                if bool::from(!p.is_identity()) {
                    return p;
                }
            }
        }
    }

    /// Returns a fixed generator of the curve in jacobian coordinates
    /// The point has been generated from the Simplified Shallue-van
    /// de Woestijne-Ulas method for hashing to elliptic curves in
    /// Short Weierstrass form, applied on the integer value of the
    /// binary encoding of the string "Cheetah".
    ///
    /// See https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve
    /// for more.
    pub const fn generator() -> JacobianPoint {
        JacobianPoint {
            x: Fp6 {
                c0: Fp(0x263a588f4b0118a1),
                c1: Fp(0x7757a0bcb26a142d),
                c2: Fp(0x9215adfc1e925890),
                c3: Fp(0x430aad2ce14759a4),
                c4: Fp(0x534ece54de4b2c8),
                c5: Fp(0xb39050f01f7b1f33),
            },
            y: Fp6 {
                c0: Fp(0xd57f0d0d47482534),
                c1: Fp(0x26821d894fa8ea0f),
                c2: Fp(0xc77f564783ef13a1),
                c3: Fp(0x949c360784284ec2),
                c4: Fp(0xb7040bd639ef3cc4),
                c5: Fp(0x8aa635f2719d255f),
            },
            z: Fp6::one(),
        }
    }

    /// Outputs a compress byte representation of this `JacobianPoint` element
    pub fn to_compressed(&self) -> CompressedPoint {
        AffinePoint::from(self).to_compressed()
    }

    /// Outputs an uncompressed byte representation of this `JacobianPoint` element
    /// It is twice larger than when calling `JacobianPoint::to_uncompress()`
    pub fn to_uncompressed(&self) -> UncompressedPoint {
        AffinePoint::from(self).to_uncompressed()
    }

    /// Attempts to deserialize an uncompressed element.
    pub fn from_uncompressed(uncompressed_point: &UncompressedPoint) -> CtOption<Self> {
        AffinePoint::from_uncompressed(uncompressed_point)
            .and_then(|p| CtOption::new(JacobianPoint::from(p), 1.into()))
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is on the curve and not checking if it is in the correct subgroup.
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_uncompressed()` instead.
    pub fn from_uncompressed_unchecked(uncompressed_point: &UncompressedPoint) -> CtOption<Self> {
        AffinePoint::from_uncompressed_unchecked(uncompressed_point)
            .and_then(|p| CtOption::new(JacobianPoint::from(p), 1.into()))
    }

    /// Attempts to deserialize a compressed element.
    pub fn from_compressed(compressed_point: &CompressedPoint) -> CtOption<Self> {
        AffinePoint::from_compressed(compressed_point).map(JacobianPoint::from)
    }

    #[allow(unused)]
    /// Constructs a `JacobianPoint` element without checking that it is a valid point.
    pub const fn from_raw_coordinates(elems: [Fp6; 3]) -> Self {
        JacobianPoint {
            x: elems[0],
            y: elems[1],
            z: elems[2],
        }
    }

    /// Computes `n` iterated doubling of this point.
    pub fn double_multi(&self, n: u32) -> JacobianPoint {
        assert!(n >= 1);
        let mut output = self.double();

        for _ in 1..n {
            output = output.double();
        }

        output
    }

    /// Computes `n` iterated doubling of this point.
    /// **This is dangerous to call unless you know that the point to be doubled
    /// is not the identity point, otherwise, API invariants may be broken.**
    /// Please consider using `double_multi()` instead.
    pub fn double_multi_unchecked(&self, n: u32) -> JacobianPoint {
        assert!(n >= 1);
        let mut output = self.double_unchecked();

        for _ in 1..n {
            output = output.double_unchecked();
        }

        output
    }

    /// Computes the doubling of this point.
    pub fn double(&self) -> JacobianPoint {
        // Use formula given in Handbook of Elliptic and Hyperelliptic Curve Cryptography, part 13.2

        let x_sq = self.x.square();
        let y_sq = self.y.square();
        let y_sq2 = y_sq.square();
        let z_sq = self.z.square();
        let z_sq2 = z_sq.square();

        let a = y_sq.mul_by_u32(4);
        let a = a * self.x;

        let b = x_sq.triple();
        let b = b + z_sq2;

        let x3 = b.square();
        let x3 = x3 - a.double();

        let y3 = a - x3;
        let y3 = y3 * b;
        let y3 = y3 - y_sq2.mul_by_u32(8);

        let z3 = self.y.double();
        let z3 = z3 * self.z;

        let tmp = JacobianPoint {
            x: x3,
            y: y3,
            z: z3,
        };

        // The calculation above fails for doubling the infinity point,
        // hence we do a final conditional selection.
        JacobianPoint::conditional_select(&tmp, &JacobianPoint::identity(), self.is_identity())
    }

    /// Computes the doubling of this point. This formulae is faster than
    /// `JacobianPoint::double` but is not working for the infinity point.
    /// **This is dangerous to call unless you know that the point to be doubled
    /// is not the identity point, otherwise, API invariants may be broken.**
    /// Please consider using `double()` instead.
    pub const fn double_unchecked(&self) -> JacobianPoint {
        // Use formula given in Handbook of Elliptic and Hyperelliptic Curve Cryptography, part 13.2

        let x_sq = (&self.x).square();
        let y_sq = (&self.y).square();
        let y_sq2 = (&y_sq).square();
        let z_sq = (&self.z).square();
        let z_sq2 = (&z_sq).square();

        let a = (&y_sq).mul_by_u32(4);
        let a = (&a).mul(&self.x);

        let b = (&x_sq).triple();
        let b = (&b).add(&z_sq2);

        let x3 = (&b).square();
        let x3 = (&x3).sub(&a.double());

        let y3 = (&a).sub(&x3);
        let y3 = (&y3).mul(&b);
        let y3 = (&y3).sub(&y_sq2.mul_by_u32(8));

        let z3 = (&self.y).double();
        let z3 = (&z3).mul(&self.z);

        JacobianPoint {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    /// Computes the negation of a point in jacobian coordinates
    #[inline]
    pub const fn neg(&self) -> Self {
        Self {
            x: self.x,
            y: (&self.y).neg(),
            z: self.z,
        }
    }

    /// Adds this point to another point.
    pub fn add(&self, rhs: &JacobianPoint) -> JacobianPoint {
        // Use formula given in Handbook of Elliptic and Hyperelliptic Curve Cryptography, part 13.2

        let z1_sq = self.z.square();
        let z2_sq = rhs.z.square();

        let a = self.x * z2_sq;

        let b = rhs.x * z1_sq;

        let c = rhs.z * z2_sq;
        let c = c * self.y;

        let d = self.z * z1_sq;
        let d = d * rhs.y;

        let e = b - a;

        let f = d - c;

        let e_sq = e.square();
        let e_sq_times_e = e * e_sq;
        let a_e_sq = a * e_sq;

        let x3 = f.square();
        let x3 = x3 - e_sq_times_e;
        let x3 = x3 - a_e_sq.double();

        let y3 = a_e_sq - x3;
        let y3 = y3 * f;
        let t3 = c * e_sq_times_e;
        let y3 = y3 - t3;

        let z3 = self.z * rhs.z;
        let z3 = z3 * e;

        let tmp = JacobianPoint {
            x: x3,
            y: y3,
            z: z3,
        };

        let tmp = JacobianPoint::conditional_select(&tmp, rhs, self.is_identity());

        let tmp = JacobianPoint::conditional_select(&tmp, self, rhs.is_identity());

        let tmp = JacobianPoint::conditional_select(
            &tmp,
            &JacobianPoint::identity(),
            rhs.ct_eq(&self.neg()),
        );

        JacobianPoint::conditional_select(&tmp, &self.double(), rhs.ct_eq(self))
    }

    /// Adds this point to another point in the affine model.
    pub fn add_mixed(&self, rhs: &AffinePoint) -> JacobianPoint {
        // Use formula given in Handbook of Elliptic and Hyperelliptic Curve Cryptography, part 13.2

        let z1_sq = self.z.square();

        let b = rhs.x * z1_sq;

        let d = self.z * z1_sq;
        let d = d * rhs.y;

        let e = b - self.x;

        let f = d - self.y;

        let e_sq = e.square();
        let e_sq_times_e = e * e_sq;
        let x_e_sq = self.x * e_sq;

        let x3 = f.square();
        let x3 = x3 - e_sq_times_e;
        let x3 = x3 - x_e_sq.double();

        let y3 = x_e_sq - x3;
        let y3 = y3 * f;
        let t3 = self.y * e_sq_times_e;
        let y3 = y3 - t3;

        let z3 = self.z * e;

        let tmp = JacobianPoint {
            x: x3,
            y: y3,
            z: z3,
        };

        let tmp = JacobianPoint::conditional_select(&tmp, self, rhs.is_identity());

        let tmp = JacobianPoint::conditional_select(
            &tmp,
            &JacobianPoint::identity(),
            rhs.ct_eq(&self.neg().into()),
        );

        JacobianPoint::conditional_select(&tmp, &self.double(), rhs.ct_eq(&self.into()))
    }

    /// Adds this point to another point in the affine model. This formulae is faster than
    /// `JacobianPoint::add_mixed` but is not complete.
    /// **This is dangerous to call unless you know that the points to be added are not
    /// identical, opposite of each other, and neither of them is the identity point; otherwise,
    /// API invariants may be broken.** Please consider using `add_mixed()` instead.
    pub fn add_mixed_unchecked(&self, rhs: &AffinePoint) -> JacobianPoint {
        // Use formula given in Handbook of Elliptic and Hyperelliptic Curve Cryptography, part 13.2

        let z1_sq = self.z.square();

        let b = rhs.x * z1_sq;

        let d = self.z * z1_sq;
        let d = d * rhs.y;

        let e = b - self.x;

        let f = d - self.y;

        let e_sq = e.square();
        let e_sq_times_e = e * e_sq;
        let x_e_sq = self.x * e_sq;

        let x3 = f.square();
        let x3 = x3 - e_sq_times_e;
        let x3 = x3 - x_e_sq.double();

        let y3 = x_e_sq - x3;
        let y3 = y3 * f;
        let t3 = self.y * e_sq_times_e;
        let y3 = y3 - t3;

        let z3 = self.z * e;

        let tmp = JacobianPoint {
            x: x3,
            y: y3,
            z: z3,
        };

        JacobianPoint::conditional_select(&tmp, self, rhs.is_identity())
    }

    /// Performs a jacobian scalar multiplication from `by`
    /// given as byte representation of a `Scalar` element
    pub fn multiply(&self, by: &[u8; 32]) -> JacobianPoint {
        let table = LookupTable::<16>::from(self);
        let digits = Scalar::bytes_to_radix_16(by);

        let mut acc = *SHIFT_POINT_MODIFIED_JACOBIAN;
        for i in (0..64).rev() {
            acc = acc.double_multi_unchecked(4);
            acc = acc.add_mixed_unchecked(&table.get_point(digits[i]));
        }

        acc.add_mixed_unchecked(&MINUS_SHIFT_POINT_ARRAY[256])
            .into()
    }

    /// Performs a jacobian scalar multiplication from `by`
    /// given as byte representation of a `Scalar` element.
    ///
    /// **This operation is variable time with respect
    /// to the scalar.** If the scalar is fixed,
    /// this operation is effectively constant time.
    pub fn multiply_vartime(&self, by: &[u8; 32]) -> JacobianPoint {
        let digits = Scalar::bytes_to_wnaf_vartime(by, 5);

        // We skip unset digits
        let mut i: usize = 255;
        for j in (0..256).rev() {
            if digits[i] != 0 {
                i = j;
                break;
            }
        }
        let table = NafLookupTable::<8>::from(self);
        let mut acc = *SHIFT_POINT_MODIFIED_JACOBIAN;

        for j in (0..i + 1).rev() {
            acc = acc.double_unchecked();

            match digits[j].cmp(&0) {
                Ordering::Greater => {
                    acc = acc.add_mixed_unchecked(&table.get_point(digits[j] as usize))
                }
                Ordering::Less => {
                    acc = acc.add_mixed_unchecked(&table.get_point(-digits[j] as usize).neg())
                }
                Ordering::Equal => (),
            };
        }

        acc.add_mixed_unchecked(&MINUS_SHIFT_POINT_ARRAY[i + 1])
            .into()
    }

    /// Performs the jacobian sum [`by_lhs` * `self` + `by_rhs` * `rhs`] with
    /// `by_lhs` and `by_rhs` given as byte representations of `Scalar` elements.
    pub fn multiply_double(
        &self,
        rhs: &JacobianPoint,
        by_lhs: &[u8; 32],
        by_rhs: &[u8; 32],
    ) -> JacobianPoint {
        let table_lhs = LookupTable::<16>::from(self);
        let table_rhs = LookupTable::<16>::from(rhs);
        let digits_lhs = Scalar::bytes_to_radix_16(by_lhs);
        let digits_rhs = Scalar::bytes_to_radix_16(by_rhs);

        let mut acc = *SHIFT_POINT_MODIFIED_JACOBIAN;
        for i in (0..64).rev() {
            acc = acc.double_multi_unchecked(4);
            acc = acc.add_mixed_unchecked(&table_lhs.get_point(digits_lhs[i]));
            acc = acc.add_mixed_unchecked(&table_rhs.get_point(digits_rhs[i]));
        }

        acc.add_mixed_unchecked(&MINUS_SHIFT_POINT_ARRAY[256])
            .into()
    }

    /// Performs the jacobian sum [`by_lhs` * `self` + `by_rhs` * `rhs`] with
    /// `by_lhs` and `by_rhs` given as byte representations of `Scalar` elements.
    ///
    /// **This operation is variable time with respect
    /// to the scalars.** If the scalars are fixed,
    /// this operation is effectively constant time.
    pub fn multiply_double_vartime(
        &self,
        rhs: &JacobianPoint,
        by_lhs: &[u8; 32],
        by_rhs: &[u8; 32],
    ) -> JacobianPoint {
        let by_lhs_digits = Scalar::bytes_to_wnaf_vartime(by_lhs, 5);
        let by_rhs_digits = Scalar::bytes_to_wnaf_vartime(by_rhs, 5);

        // We skip unset digits
        let mut i: usize = 255;
        for j in (0..256).rev() {
            if by_lhs_digits[i] != 0 || by_rhs_digits[i] != 0 {
                i = j;
                break;
            }
        }
        let table_self = NafLookupTable::<8>::from(self);
        let table_rhs = NafLookupTable::<8>::from(rhs);
        let mut acc = *SHIFT_POINT_MODIFIED_JACOBIAN;

        for j in (0..i + 1).rev() {
            acc = acc.double_unchecked();

            match by_lhs_digits[j].cmp(&0) {
                Ordering::Greater => {
                    acc = acc.add_mixed_unchecked(&table_self.get_point(by_lhs_digits[j] as usize))
                }
                Ordering::Less => {
                    acc = acc.add_mixed_unchecked(
                        &table_self.get_point(-by_lhs_digits[j] as usize).neg(),
                    )
                }
                Ordering::Equal => (),
            };

            match by_rhs_digits[j].cmp(&0) {
                Ordering::Greater => {
                    acc = acc.add_mixed_unchecked(&table_rhs.get_point(by_rhs_digits[j] as usize))
                }
                Ordering::Less => {
                    acc = acc
                        .add_mixed_unchecked(&table_rhs.get_point(-by_rhs_digits[j] as usize).neg())
                }
                Ordering::Equal => (),
            };
        }

        acc.add_mixed_unchecked(&MINUS_SHIFT_POINT_ARRAY[i + 1])
            .into()
    }

    // TODO: It would be nice to have a constant-time variant of the method below,
    // though this may be tricky for accessing the point in the table.

    /// Performs the affine sum [`by_self` * `self` + `by_basepoint` * `ODD_MULTIPLES_BASEPOINT`]
    /// with `by_self` and `by_basepoint` given as byte representations of `Scalar` elements.
    ///
    /// This operation is useful to speed-up verification of signature schemes
    /// such as ECDSA or Schnorr signatures.
    ///
    /// **This operation is variable time with respect
    /// to the scalars.** If the scalars are fixed,
    /// this operation is effectively constant time.
    pub fn multiply_double_with_basepoint_vartime(
        &self,
        by_self: &[u8; 32],
        by_basepoint: &[u8; 32],
    ) -> JacobianPoint {
        let by_self_digits = Scalar::bytes_to_wnaf_vartime(by_self, 5);
        let by_basepoint_digits = Scalar::bytes_to_wnaf_vartime(by_basepoint, 8);

        // We skip unset digits
        let mut i: usize = 255;
        for j in (0..256).rev() {
            if by_self_digits[i] != 0 || by_basepoint_digits[i] != 0 {
                i = j;
                break;
            }
        }
        let table_self = NafLookupTable::<8>::from(self);
        let table_basepoint = &ODD_MULTIPLES_BASEPOINT;
        let mut acc = *SHIFT_POINT_MODIFIED_JACOBIAN;

        for j in (0..i + 1).rev() {
            acc = acc.double_unchecked();

            match by_self_digits[j].cmp(&0) {
                Ordering::Greater => {
                    acc = acc.add_mixed_unchecked(&table_self.get_point(by_self_digits[j] as usize))
                }
                Ordering::Less => {
                    acc = acc.add_mixed_unchecked(
                        &table_self.get_point(-by_self_digits[j] as usize).neg(),
                    )
                }
                Ordering::Equal => (),
            };

            match by_basepoint_digits[j].cmp(&0) {
                Ordering::Greater => {
                    acc = acc.add_mixed_unchecked(
                        &table_basepoint.get_point(by_basepoint_digits[j] as usize),
                    )
                }
                Ordering::Less => {
                    acc = acc.add_mixed_unchecked(
                        &table_basepoint
                            .get_point(-by_basepoint_digits[j] as usize)
                            .neg(),
                    )
                }
                Ordering::Equal => (),
            };
        }

        acc.add_mixed_unchecked(&MINUS_SHIFT_POINT_ARRAY[i + 1])
            .into()
    }

    /// Performs the jacobian multiscalar multiplication ∑ s[i].p[i] with
    /// the s[i] given as byte representations of `Scalar` elements.
    // TODO: this is significantly slower than its vartime counterpart, because of the
    // construction of LookupTable<16> vs NafLookupTable<5>, which is about twice slower.
    // It would be nice to implement a constant-time variation of the latter, including the
    // wnaf conversion of the scalar.
    pub fn multiply_many(points: &[JacobianPoint], scalars: &[[u8; 32]]) -> JacobianPoint {
        let digits: Vec<[i8; 64]> = scalars.iter().map(Scalar::bytes_to_radix_16).collect();

        let tables: Vec<LookupTable<16>> = points.iter().map(LookupTable::<16>::from).collect();

        let mut acc = *SHIFT_POINT_MODIFIED_JACOBIAN;

        for i in (0..64).rev() {
            acc = acc.double_multi_unchecked(4);
            for (table, digit) in tables.iter().zip(&digits) {
                acc = acc.add_mixed_unchecked(&table.get_point(digit[i]));
            }
        }

        acc.add_mixed_unchecked(&MINUS_SHIFT_POINT_ARRAY[256])
            .into()
    }

    /// Performs the jacobian multiscalar multiplication ∑ s[i].p[i] with
    /// the s[i] given as byte representations of `Scalar` elements.
    ///
    /// **This operation is variable time with respect
    /// to the scalars.** If the scalars are fixed,
    /// this operation is effectively constant time.
    pub fn multiply_many_vartime(points: &[JacobianPoint], scalars: &[[u8; 32]]) -> JacobianPoint {
        let digits: Vec<[i8; 256]> = scalars
            .iter()
            .map(|s| Scalar::bytes_to_wnaf_vartime(s, 5))
            .collect();

        let tables: Vec<NafLookupTable<8>> = points.iter().map(NafLookupTable::<8>::from).collect();

        let mut acc = *SHIFT_POINT_MODIFIED_JACOBIAN;

        for i in (0..256).rev() {
            acc = acc.double_unchecked();
            for (table, digit) in tables.iter().zip(&digits) {
                match digit[i].cmp(&0) {
                    Ordering::Greater => {
                        acc = acc.add_mixed_unchecked(&table.get_point(digit[i] as usize))
                    }
                    Ordering::Less => {
                        acc = acc.add_mixed_unchecked(&table.get_point(-digit[i] as usize).neg())
                    }
                    Ordering::Equal => (),
                };
            }
        }

        acc.add_mixed_unchecked(&MINUS_SHIFT_POINT_ARRAY[256])
            .into()
    }

    /// Multiplies by the curve cofactor
    pub fn clear_cofactor(&self) -> JacobianPoint {
        // cofactor = 708537115134665106932687062569690615370
        //          = 0x2150b48e071ef610049bc3f5d54304e4a
        const COFACTOR_BYTES: [u8; 32] = [
            74, 78, 48, 84, 93, 63, 188, 73, 0, 97, 239, 113, 224, 72, 11, 21, 2, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];

        self.multiply_vartime(&COFACTOR_BYTES)
    }

    /// Returns true if this point is free of an $h$-torsion component.
    /// This should always return true unless an "unchecked" API was used.
    pub fn is_torsion_free(&self) -> Choice {
        const FQ_MODULUS_BYTES: [u8; 32] = [
            207, 172, 212, 174, 62, 98, 67, 212, 34, 119, 21, 48, 35, 167, 122, 50, 181, 55, 10,
            153, 15, 191, 63, 86, 208, 34, 63, 59, 155, 89, 242, 122,
        ];

        // Clear the r-torsion from the point and check if it is the identity
        self.multiply_vartime(&FQ_MODULUS_BYTES).is_identity()
    }

    /// Converts a batch of `JacobianPoint` elements into `AffinePoint` elements. This
    /// function will panic if `p.len() != q.len()`.
    pub fn batch_normalize(p: &[Self], q: &mut [AffinePoint]) {
        assert_eq!(p.len(), q.len());

        let mut acc = Fp6::one();
        for (p, q) in p.iter().zip(q.iter_mut()) {
            // We use the `x` field of `AffinePoint` to store the product
            // of previous z-coordinates seen.
            q.x = acc;

            // We will end up skipping all identities in p
            acc = Fp6::conditional_select(&(acc * p.z), &acc, p.is_identity());
        }

        // This is the inverse, as all z-coordinates are nonzero and the ones
        // that are not are skipped.
        acc = acc.invert().unwrap();

        for (p, q) in p.iter().rev().zip(q.iter_mut().rev()) {
            let skip = p.is_identity();

            // Compute tmp = 1/z`
            let tmp = q.x * acc;
            let tmp2 = tmp.square();

            // Cancel out z-coordinate in denominator of `acc`
            acc = Fp6::conditional_select(&(acc * p.z), &acc, skip);

            // Set the coordinates to the correct value
            q.x = p.x * tmp2;
            q.y = p.y * tmp2 * tmp;
            q.infinity = Choice::from(0);

            *q = AffinePoint::conditional_select(q, &AffinePoint::identity(), skip);
        }
    }

    /// Returns true if this element is the identity (the point at infinity).
    #[inline]
    pub fn is_identity(&self) -> Choice {
        self.z.is_zero()
    }

    /// Returns true if this point is on the curve. This should always return
    /// true unless an "unchecked" API was used.
    pub fn is_on_curve(&self) -> Choice {
        // Y^2 = X^3 + X Z^4 + b Z^6

        let z2 = self.z.square();
        let z4 = z2.square();

        (self.y.square()).ct_eq(&(self.x.square() * self.x + self.x * z4 + z2 * z4 * B))
            | (self.z.is_zero())
    }
}

/// A modified jacobian point
#[derive(Copy, Clone, Debug)]
pub struct ModifiedJacobianPoint {
    pub(crate) x: Fp6,
    pub(crate) y: Fp6,
    pub(crate) z: Fp6,
    pub(crate) z4: Fp6,
}

impl Default for ModifiedJacobianPoint {
    fn default() -> ModifiedJacobianPoint {
        ModifiedJacobianPoint::identity()
    }
}

impl zeroize::DefaultIsZeroes for ModifiedJacobianPoint {}

impl fmt::Display for ModifiedJacobianPoint {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl Hash for ModifiedJacobianPoint {
    fn hash<H: Hasher>(&self, hasher: &mut H) {
        AffinePoint::from(self).hash(hasher);
    }
}

impl<'a> From<&'a AffinePoint> for ModifiedJacobianPoint {
    fn from(p: &'a AffinePoint) -> ModifiedJacobianPoint {
        ModifiedJacobianPoint {
            x: p.x,
            y: p.y,
            z: Fp6::conditional_select(&Fp6::one(), &Fp6::zero(), p.infinity),
            z4: Fp6::conditional_select(&Fp6::one(), &Fp6::zero(), p.infinity),
        }
    }
}

impl From<AffinePoint> for ModifiedJacobianPoint {
    fn from(p: AffinePoint) -> ModifiedJacobianPoint {
        ModifiedJacobianPoint::from(&p)
    }
}

impl<'a> From<&'a JacobianPoint> for ModifiedJacobianPoint {
    fn from(p: &'a JacobianPoint) -> ModifiedJacobianPoint {
        ModifiedJacobianPoint {
            x: p.x,
            y: p.y,
            z: p.z,
            z4: p.z.square().square(),
        }
    }
}

impl From<JacobianPoint> for ModifiedJacobianPoint {
    fn from(p: JacobianPoint) -> ModifiedJacobianPoint {
        ModifiedJacobianPoint::from(&p)
    }
}

impl ConstantTimeEq for ModifiedJacobianPoint {
    fn ct_eq(&self, other: &Self) -> Choice {
        // Is (xz^2, yz^3, z) equal to (x'z'^2, y'z'^3, z') when converted to affine?

        let z1_squared = self.z.square();
        let z2_squared = other.z.square();

        let x1 = self.x * z2_squared;
        let x2 = other.x * z1_squared;

        let y1 = self.y * z2_squared * other.z;
        let y2 = other.y * z1_squared * self.z;

        let self_is_zero = self.z.is_zero();
        let other_is_zero = other.z.is_zero();

        (self_is_zero & other_is_zero) // Both point at infinity
            | ((!self_is_zero) & (!other_is_zero) & x1.ct_eq(&x2) & y1.ct_eq(&y2))
        // Neither point at infinity, coordinates are the same
    }
}

impl ConditionallySelectable for ModifiedJacobianPoint {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        ModifiedJacobianPoint {
            x: Fp6::conditional_select(&a.x, &b.x, choice),
            y: Fp6::conditional_select(&a.y, &b.y, choice),
            z: Fp6::conditional_select(&a.z, &b.z, choice),
            z4: Fp6::conditional_select(&a.z4, &b.z4, choice),
        }
    }
}

impl Eq for ModifiedJacobianPoint {}
impl PartialEq for ModifiedJacobianPoint {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl<'a> Neg for &'a ModifiedJacobianPoint {
    type Output = ModifiedJacobianPoint;

    #[inline]
    fn neg(self) -> ModifiedJacobianPoint {
        self.neg()
    }
}

impl Neg for ModifiedJacobianPoint {
    type Output = ModifiedJacobianPoint;

    #[inline]
    fn neg(self) -> ModifiedJacobianPoint {
        -&self
    }
}

impl ModifiedJacobianPoint {
    /// Returns the identity of the group: the point at infinity.
    pub const fn identity() -> ModifiedJacobianPoint {
        ModifiedJacobianPoint {
            x: Fp6::one(),
            y: Fp6::one(),
            z: Fp6::zero(),
            z4: Fp6::zero(),
        }
    }

    /// Computes `n` iterated doubling of this point.
    /// **This is dangerous to call unless you know that the point to be doubled
    /// is not the identity point, otherwise, API invariants may be broken.**
    /// Please consider using `double_multi()` instead.
    pub fn double_multi_unchecked(&self, n: u32) -> ModifiedJacobianPoint {
        assert!(n >= 1);
        let mut output = self.double_unchecked();

        for _ in 1..n {
            output = output.double_unchecked();
        }

        output
    }

    /// Computes the doubling of this point. This formulae is faster than
    /// `ModifiedJacobianPoint::double` but is not working for the infinity point.
    /// **This is dangerous to call unless you know that the point to be doubled
    /// is not the identity point, otherwise, API invariants may be broken.**
    /// Please consider using `double()` instead.
    pub const fn double_unchecked(&self) -> ModifiedJacobianPoint {
        // Use formula given in Handbook of Elliptic and Hyperelliptic Curve Cryptography, part 13.2

        let x_sq = (&self.x).square();
        let y_sq = (&self.y).square();
        let y_sq2 = (&y_sq).square();

        let a = (&y_sq).mul_by_u32(4);
        let a = (&a).mul(&self.x);

        let b = (&x_sq).triple();
        let b = (&b).add(&self.z4);

        let c = (&y_sq2).mul_by_u32(8);

        let x3 = (&b).square();
        let x3 = (&x3).sub(&a.double());

        let y3 = (&a).sub(&x3);
        let y3 = (&y3).mul(&b);
        let y3 = (&y3).sub(&c);

        let z3 = (&self.y).double();
        let z3 = (&z3).mul(&self.z);

        let z34 = (&c).double();
        let z34 = (&z34).mul(&self.z4);

        ModifiedJacobianPoint {
            x: x3,
            y: y3,
            z: z3,
            z4: z34,
        }
    }

    /// Computes the negation of a point in modified jacobian coordinates
    #[inline]
    pub const fn neg(&self) -> Self {
        Self {
            x: self.x,
            y: (&self.y).neg(),
            z: self.z,
            z4: self.z4,
        }
    }

    /// Adds this point to another point in the affine model. This formulae is faster than
    /// `ModifiedJacobianPoint::add_mixed` but is not complete.
    /// **This is dangerous to call unless you know that the points to be added are not
    /// identical, opposite of each other, and neither of them is the identity point; otherwise,
    /// API invariants may be broken.** Please consider using `add_mixed()` instead.
    pub fn add_mixed_unchecked(&self, rhs: &AffinePoint) -> ModifiedJacobianPoint {
        // Use formula given in Handbook of Elliptic and Hyperelliptic Curve Cryptography, part 13.2

        let z1_sq = self.z.square();

        let b = rhs.x * z1_sq;

        let d = self.z * z1_sq;
        let d = d * rhs.y;

        let e = b - self.x;

        let f = d - self.y;

        let e_sq = e.square();
        let e_sq_times_e = e * e_sq;
        let x_e_sq = self.x * e_sq;

        let x3 = f.square();
        let x3 = x3 - e_sq_times_e;
        let x3 = x3 - x_e_sq.double();

        let y3 = x_e_sq - x3;
        let y3 = y3 * f;
        let t3 = self.y * e_sq_times_e;
        let y3 = y3 - t3;

        let z3 = self.z * e;

        let z34 = z3.square();
        let z34 = z34.square();

        let tmp = ModifiedJacobianPoint {
            x: x3,
            y: y3,
            z: z3,
            z4: z34,
        };

        ModifiedJacobianPoint::conditional_select(&tmp, self, rhs.is_identity())
    }

    /// Returns true if this element is the identity (the point at infinity).
    #[inline]
    pub fn is_identity(&self) -> Choice {
        self.z.is_zero()
    }
}

// GROUP TRAITS IMPLEMENTATION
// ================================================================================================

impl Group for JacobianPoint {
    type Scalar = Scalar;

    fn random(mut rng: impl RngCore) -> Self {
        Self::random(&mut rng)
    }

    fn identity() -> Self {
        Self::identity()
    }

    fn generator() -> Self {
        Self::generator()
    }

    fn is_identity(&self) -> Choice {
        self.is_identity()
    }

    #[must_use]
    fn double(&self) -> Self {
        self.double()
    }
}

impl Curve for JacobianPoint {
    type AffineRepr = AffinePoint;

    fn batch_normalize(p: &[Self], q: &mut [Self::AffineRepr]) {
        Self::batch_normalize(p, q);
    }

    fn to_affine(&self) -> Self::AffineRepr {
        self.into()
    }
}

// SERDE SERIALIZATION
// ================================================================================================

#[cfg(feature = "serialize")]
impl Serialize for JacobianPoint {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        self.to_compressed().serialize(serializer)
    }
}

#[cfg(feature = "serialize")]
impl<'de> Deserialize<'de> for JacobianPoint {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let compressed_point = CompressedPoint::deserialize(deserializer)?;
        let p = JacobianPoint::from_compressed(&compressed_point);
        if bool::from(p.is_some()) {
            Ok(p.unwrap())
        } else {
            Err(serde::de::Error::custom("decompression failed"))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::OsRng;

    use crate::BASEPOINT_LOOKUP;

    #[test]
    fn test_is_on_curve() {
        assert!(bool::from(JacobianPoint::identity().is_on_curve()));
        assert!(bool::from(JacobianPoint::generator().is_on_curve()));

        let z = Fp6 {
            c0: Fp(0x79b230ab69fdb493),
            c1: Fp(0xfe0bff52056ea77b),
            c2: Fp(0x68ae55ee9963bb1d),
            c3: Fp(0xaed3160aa9927b96),
            c4: Fp(0x57f9b8ca1372d78e),
            c5: Fp(0xcfb89c494391a555),
        };

        let gen = AffinePoint::generator();
        let mut test = JacobianPoint {
            x: gen.x * z.square(),
            y: gen.y * z.square() * z,
            z,
        };

        assert!(bool::from(test.is_on_curve()));

        test.x = z;
        assert!(!bool::from(test.is_on_curve()));
    }

    #[test]
    #[allow(clippy::eq_op)]
    fn test_jacobian_point_equality() {
        let a = JacobianPoint::generator();
        let b = JacobianPoint::identity();
        let c = JacobianPoint::default();

        assert!(a == a);
        assert!(b == b);
        assert!(b == c);
        assert!(a != b);
        assert!(b != a);

        assert!(bool::from(b.is_identity()));
        assert!(!bool::from(a.ct_eq(&b)));

        let z = Fp6 {
            c0: Fp(0x79b230ab69fdb493),
            c1: Fp(0xfe0bff52056ea77b),
            c2: Fp(0x68ae55ee9963bb1d),
            c3: Fp(0xaed3160aa9927b96),
            c4: Fp(0x57f9b8ca1372d78e),
            c5: Fp(0xcfb89c494391a555),
        };

        let mut c = JacobianPoint {
            x: a.x * z.square(),
            y: a.y * z.square() * z,
            z,
        };
        assert!(bool::from(c.is_on_curve()));

        assert!(a == c);
        assert!(b != c);
        assert!(c == a);
        assert!(c != b);

        c.y = -c.y;
        assert!(bool::from(c.is_on_curve()));

        assert!(a != c);
        assert!(b != c);
        assert!(c != a);
        assert!(c != b);

        c.y = -c.y;
        c.x = z;
        assert!(!bool::from(c.is_on_curve()));
        assert!(a != b);
        assert!(a != c);
        assert!(b != c);
    }

    #[test]
    fn test_jacobian_to_affine() {
        let a = JacobianPoint::generator();
        let b = JacobianPoint::identity();

        assert!(bool::from(AffinePoint::from(a).is_on_curve()));
        assert!(!bool::from(AffinePoint::from(a).is_identity()));
        assert!(bool::from(AffinePoint::from(b).is_on_curve()));
        assert!(bool::from(AffinePoint::from(b).is_identity()));

        let z = Fp6 {
            c0: Fp(0x79b230ab69fdb493),
            c1: Fp(0xfe0bff52056ea77b),
            c2: Fp(0x68ae55ee9963bb1d),
            c3: Fp(0xaed3160aa9927b96),
            c4: Fp(0x57f9b8ca1372d78e),
            c5: Fp(0xcfb89c494391a555),
        };

        let c = JacobianPoint {
            x: a.x * z.square(),
            y: a.y * z.square() * z,
            z,
        };

        assert_eq!(AffinePoint::from(c), AffinePoint::generator());

        let mut rng = OsRng;
        for _ in 0..100 {
            let a = AffinePoint::random(&mut rng);
            let z = Fp6::random(&mut rng);

            let b = JacobianPoint {
                x: a.x * z.square(),
                y: a.y * z.square() * z,
                z,
            };

            assert_eq!(AffinePoint::from(b), a);
        }
    }

    #[test]
    fn test_affine_to_jacobian() {
        let a = AffinePoint::generator();
        let b = AffinePoint::identity();

        assert!(bool::from(JacobianPoint::from(a).is_on_curve()));
        assert!(!bool::from(JacobianPoint::from(a).is_identity()));
        assert!(bool::from(JacobianPoint::from(b).is_on_curve()));
        assert!(bool::from(JacobianPoint::from(b).is_identity()));
    }

    #[test]
    fn test_doubling() {
        {
            let tmp = JacobianPoint::identity().double();
            assert!(bool::from(tmp.is_identity()));
            assert!(bool::from(tmp.is_on_curve()));
        }
        {
            let tmp = JacobianPoint::generator().double();
            assert!(!bool::from(tmp.is_identity()));
            assert!(bool::from(tmp.is_on_curve()));
            assert_eq!(tmp, JacobianPoint::generator().double_unchecked());

            assert_eq!(
                AffinePoint::from(tmp),
                AffinePoint {
                    x: Fp6 {
                        c0: Fp(0x7f4c1bfc52278ad8),
                        c1: Fp(0xfa8e921f7580e371),
                        c2: Fp(0x97252bf35d1c7668),
                        c3: Fp(0xe6d0901604cae95a),
                        c4: Fp(0xae36bba2ad2ee0d7),
                        c5: Fp(0x194b4e35a2a9c77),
                    },
                    y: Fp6 {
                        c0: Fp(0x144045efbce03ef8),
                        c1: Fp(0x8e5fe3f66f8b370d),
                        c2: Fp(0x3d54df63b96bfd20),
                        c3: Fp(0x2418219e37948caa),
                        c4: Fp(0xd4c1a40432582552),
                        c5: Fp(0x367b029f5f146e3d),
                    },
                    infinity: Choice::from(0u8),
                }
            );
        }
    }

    #[test]
    fn test_jacobian_addition() {
        {
            let a = JacobianPoint::identity();
            let b = JacobianPoint::identity();
            let c = a + b;
            assert!(bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
        }
        {
            let a = JacobianPoint::identity();
            let mut b = JacobianPoint::generator();
            {
                let z = Fp6 {
                    c0: Fp(0xa3c62f9770336022),
                    c1: Fp(0x69a1531152c92fc1),
                    c2: Fp(0x1636d9b90656a08a),
                    c3: Fp(0x635289066593aaf6),
                    c4: Fp(0x1e2178e3c54f6682),
                    c5: Fp(0xe15458e6d847f393),
                };

                b = JacobianPoint {
                    x: b.x * z.square(),
                    y: b.y * z.square() * z,
                    z,
                };
            }
            let c = a + b;
            assert!(!bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
            assert!(c == JacobianPoint::generator());
        }
        {
            let a = JacobianPoint::identity();
            let mut b = JacobianPoint::generator();
            {
                let z = Fp6 {
                    c0: Fp(0xa3c62f9770336022),
                    c1: Fp(0x69a1531152c92fc1),
                    c2: Fp(0x1636d9b90656a08a),
                    c3: Fp(0x635289066593aaf6),
                    c4: Fp(0x1e2178e3c54f6682),
                    c5: Fp(0xe15458e6d847f393),
                };

                b = JacobianPoint {
                    x: b.x * z.square(),
                    y: b.y * z.square() * z,
                    z,
                };
            }
            let c = b + a;
            assert!(!bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
            assert!(c == JacobianPoint::generator());
        }
        {
            let a = JacobianPoint::generator().double_multi(2); // 4P
            let b = JacobianPoint::generator().double(); // 2P
            let c = a + b;

            let mut d = JacobianPoint::generator();
            for _ in 0..5 {
                d += JacobianPoint::generator();
            }
            assert!(!bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
            assert!(!bool::from(d.is_identity()));
            assert!(bool::from(d.is_on_curve()));
            assert_eq!(c, d);
        }
    }

    #[test]
    fn test_mixed_addition() {
        {
            let a = AffinePoint::identity();
            let b = JacobianPoint::identity();
            let c = a + b;
            assert!(bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
        }
        {
            let a = AffinePoint::identity();
            let mut b = JacobianPoint::generator();
            {
                let z = Fp6 {
                    c0: Fp(0xa3c62f9770336022),
                    c1: Fp(0x69a1531152c92fc1),
                    c2: Fp(0x1636d9b90656a08a),
                    c3: Fp(0x635289066593aaf6),
                    c4: Fp(0x1e2178e3c54f6682),
                    c5: Fp(0xe15458e6d847f393),
                };

                b = JacobianPoint {
                    x: b.x * z.square(),
                    y: b.y * z.square() * z,
                    z,
                };
            }
            let c = a + b;
            assert!(!bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
            assert!(c == JacobianPoint::generator());
        }
        {
            let a = AffinePoint::identity();
            let mut b = JacobianPoint::generator();
            {
                let z = Fp6 {
                    c0: Fp(0xa3c62f9770336022),
                    c1: Fp(0x69a1531152c92fc1),
                    c2: Fp(0x1636d9b90656a08a),
                    c3: Fp(0x635289066593aaf6),
                    c4: Fp(0x1e2178e3c54f6682),
                    c5: Fp(0xe15458e6d847f393),
                };

                b = JacobianPoint {
                    x: b.x * z.square(),
                    y: b.y * z.square() * z,
                    z,
                };
            }
            let c = b + a;
            assert!(!bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
            assert!(c == JacobianPoint::generator());
        }
        {
            let a = JacobianPoint::generator().double_multi(2); // 4P
            let b = JacobianPoint::generator().double(); // 2P
            let c = a + b;

            let mut d = JacobianPoint::generator();
            for _ in 0..5 {
                d += AffinePoint::generator();
            }
            assert!(!bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
            assert!(!bool::from(d.is_identity()));
            assert!(bool::from(d.is_on_curve()));
            assert_eq!(c, d);
        }
        {
            let mut rng = OsRng;
            for _ in 0..100 {
                let a = JacobianPoint::random(&mut rng);
                let b = AffinePoint::random(&mut rng);

                // If this fails, we are very unlucky!
                assert_eq!(a.add_mixed(&b), a.add_mixed_unchecked(&b));
            }

            let a = JacobianPoint::random(&mut rng);
            assert_ne!(a.add(&a), a.add_mixed_unchecked(&a.to_affine()));
        }
    }

    #[test]
    #[allow(clippy::eq_op)]
    fn test_jacobian_negation_and_subtraction() {
        let a = JacobianPoint::generator().double();
        assert_eq!(a + (-a), JacobianPoint::identity());
        assert_eq!(a + (-a), a - a);
    }

    #[test]
    fn test_jacobian_scalar_multiplication() {
        let mut rng = OsRng;
        let g = JacobianPoint::generator();

        for _ in 0..100 {
            let a = Scalar::random(&mut rng);
            let b = Scalar::random(&mut rng);

            let c = a * b;

            assert_eq!((g * a) * b, g * c);
            assert_eq!(&BASEPOINT_LOOKUP * c, g * c);
            assert_eq!(g * c, g.multiply_vartime(&c.to_bytes()));
        }
    }

    #[test]
    fn test_jacobian_double_scalar_multiplication() {
        let mut rng = OsRng;
        let g = JacobianPoint::generator() * Scalar::random(&mut rng);
        let h = JacobianPoint::generator() * Scalar::random(&mut rng);

        for _ in 0..100 {
            let a = Scalar::random(&mut rng);
            let b = Scalar::random(&mut rng);

            assert_eq!(
                g.multiply_double(&h, &a.to_bytes(), &b.to_bytes()),
                g * a + h * b
            );
            assert_eq!(
                g.multiply_double_vartime(&h, &a.to_bytes(), &b.to_bytes()),
                g * a + h * b
            );
            assert_eq!(
                JacobianPoint::multiply_many(&[g, h], &[a.to_bytes(), b.to_bytes()]),
                g * a + h * b
            );
            assert_eq!(
                JacobianPoint::multiply_many_vartime(&[g, h], &[a.to_bytes(), b.to_bytes()]),
                g * a + h * b
            );
        }
    }

    #[test]
    fn test_jacobian_double_scalar_multiplication_with_basepoint() {
        let mut rng = OsRng;
        let g = JacobianPoint::generator() * Scalar::random(&mut rng);
        let h = JacobianPoint::generator();

        for _ in 0..100 {
            let a = Scalar::random(&mut rng);
            let b = Scalar::random(&mut rng);

            assert_eq!(
                g.multiply_double_with_basepoint_vartime(&a.to_bytes(), &b.to_bytes()),
                g * a + h * b
            );
        }
    }

    #[test]
    fn test_clear_cofactor() {
        // the generator (and the identity) are always on the curve
        let generator = JacobianPoint::generator();
        assert!(bool::from(generator.clear_cofactor().is_on_curve()));
        let id = JacobianPoint::identity();
        assert!(bool::from(id.clear_cofactor().is_on_curve()));

        let point = JacobianPoint {
            x: Fp6 {
                c0: Fp(0x9bfcd3244afcb637),
                c1: Fp(0x39005e478830b187),
                c2: Fp(0x7046f1c03b42c6cc),
                c3: Fp(0xb5eeac99193711e5),
                c4: Fp(0x7fd272e724307b98),
                c5: Fp(0xcc371dd6dd5d8625),
            },
            y: Fp6 {
                c0: Fp(0x9d03fdc216dfaae8),
                c1: Fp(0xbf4ade2a7665d9b8),
                c2: Fp(0xf08b022d5b3262b7),
                c3: Fp(0x2eaf583a3cf15c6f),
                c4: Fp(0xa92531e4b1338285),
                c5: Fp(0x5b8157814141a7a7),
            },
            z: Fp6::one(),
        };

        assert!(bool::from(point.is_on_curve()));
        assert!(!bool::from(AffinePoint::from(point).is_torsion_free()));
        let cleared_point = point.clear_cofactor();
        assert!(bool::from(cleared_point.is_on_curve()));
        assert!(bool::from(
            AffinePoint::from(cleared_point).is_torsion_free()
        ));
    }

    #[test]
    fn test_is_torsion_free() {
        let a = AffinePoint {
            x: Fp6 {
                c0: Fp(0x9bfcd3244afcb637),
                c1: Fp(0x39005e478830b187),
                c2: Fp(0x7046f1c03b42c6cc),
                c3: Fp(0xb5eeac99193711e5),
                c4: Fp(0x7fd272e724307b98),
                c5: Fp(0xcc371dd6dd5d8625),
            },
            y: Fp6 {
                c0: Fp(0x9d03fdc216dfaae8),
                c1: Fp(0xbf4ade2a7665d9b8),
                c2: Fp(0xf08b022d5b3262b7),
                c3: Fp(0x2eaf583a3cf15c6f),
                c4: Fp(0xa92531e4b1338285),
                c5: Fp(0x5b8157814141a7a7),
            },
            infinity: Choice::from(0u8),
        };

        let a = JacobianPoint::from(&a);
        assert!(bool::from(a.is_on_curve()));
        assert!(!bool::from(a.is_torsion_free()));
        assert!(bool::from(JacobianPoint::identity().is_torsion_free()));
        assert!(bool::from(JacobianPoint::generator().is_torsion_free()));
    }

    #[test]
    fn test_batch_normalize() {
        let a = JacobianPoint::generator().double();
        let b = a.double();
        let c = b.double();

        for a_identity in (0..1).map(|n| n == 1) {
            for b_identity in (0..1).map(|n| n == 1) {
                for c_identity in (0..1).map(|n| n == 1) {
                    let mut v = [a, b, c];
                    if a_identity {
                        v[0] = JacobianPoint::identity()
                    }
                    if b_identity {
                        v[1] = JacobianPoint::identity()
                    }
                    if c_identity {
                        v[2] = JacobianPoint::identity()
                    }

                    let mut t = [
                        AffinePoint::identity(),
                        AffinePoint::identity(),
                        AffinePoint::identity(),
                    ];
                    let expected = [
                        AffinePoint::from(v[0]),
                        AffinePoint::from(v[1]),
                        AffinePoint::from(v[2]),
                    ];

                    JacobianPoint::batch_normalize(&v[..], &mut t[..]);

                    assert_eq!(&t[..], &expected[..]);
                }
            }
        }
    }

    #[test]
    fn test_zeroize() {
        use zeroize::Zeroize;

        let mut a = JacobianPoint::generator();
        a.zeroize();
        assert_eq!(a, JacobianPoint::identity());
    }

    // POINT COMPRESSION
    // ================================================================================================

    #[test]
    fn test_point_compressed() {
        let mut rng = OsRng;
        // Random points
        for _ in 0..100 {
            let point = JacobianPoint::random(&mut rng);
            let bytes = point.to_compressed();
            let point_decompressed = JacobianPoint::from_compressed(&bytes).unwrap();
            assert_eq!(point, point_decompressed);
        }

        // Identity point
        {
            let bytes = JacobianPoint::identity().to_compressed();
            let point_decompressed = JacobianPoint::from_compressed(&bytes).unwrap();
            assert!(bool::from(point_decompressed.is_identity()));

            assert_eq!(
                JacobianPoint::identity().to_compressed().to_bytes(),
                [
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 128
                ]
            );
        }

        // Invalid points
        {
            let point = AffinePoint {
                x: Fp6::zero(),
                y: Fp6::zero(),
                infinity: 0.into(),
            };

            let point = JacobianPoint::from(&point);
            let bytes = point.to_compressed();
            let point_decompressed = JacobianPoint::from_compressed(&bytes);
            assert!(bool::from(point_decompressed.is_none()));
        }
    }

    #[test]
    fn test_point_uncompressed() {
        let mut rng = OsRng;

        // Random points
        for _ in 0..100 {
            let point = JacobianPoint::random(&mut rng);
            let bytes = point.to_uncompressed();
            let point_decompressed = JacobianPoint::from_uncompressed(&bytes).unwrap();
            assert_eq!(point, point_decompressed);
        }

        // Identity point
        {
            let bytes = JacobianPoint::identity().to_uncompressed();
            let point_decompressed = JacobianPoint::from_uncompressed(&bytes).unwrap();
            assert!(bool::from(point_decompressed.is_identity()));

            assert_eq!(
                JacobianPoint::identity().to_uncompressed().to_bytes(),
                [
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 128
                ]
            );
        }

        // Invalid points
        {
            let point = AffinePoint {
                x: Fp6::zero(),
                y: Fp6::zero(),
                infinity: 0.into(),
            };

            let point = JacobianPoint::from(&point);
            let bytes = point.to_uncompressed();
            let point_decompressed = JacobianPoint::from_uncompressed(&bytes);
            assert!(bool::from(point_decompressed.is_none()));
        }
        {
            let bytes = UncompressedPoint([
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 255,
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            ]);

            let point_decompressed = JacobianPoint::from_uncompressed_unchecked(&bytes);
            assert!(bool::from(point_decompressed.is_none()));
        }
        {
            let bytes = UncompressedPoint([
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255,
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
            ]);

            let point_decompressed = JacobianPoint::from_uncompressed_unchecked(&bytes);
            assert!(bool::from(point_decompressed.is_none()));
        }
    }

    // SERDE SERIALIZATIOIN
    // ================================================================================================

    #[test]
    #[cfg(feature = "serialize")]
    fn test_serde_jacobian() {
        let mut rng = OsRng;
        let point = JacobianPoint::random(&mut rng);
        let encoded = bincode::serialize(&point).unwrap();
        let parsed: JacobianPoint = bincode::deserialize(&encoded).unwrap();
        assert_eq!(parsed, point);

        // Check that the encoding is 49 bytes exactly
        assert_eq!(encoded.len(), 49);

        // Check that the encoding itself matches the usual one
        assert_eq!(
            point,
            bincode::deserialize(&point.to_compressed().0).unwrap()
        );

        // Check that invalid encodings fail
        let encoded = bincode::serialize(&point).unwrap();
        assert!(bincode::deserialize::<JacobianPoint>(&encoded[0..47]).is_err());

        let encoded = [255; 49];
        assert!(bincode::deserialize::<JacobianPoint>(&encoded).is_err());
    }
}
