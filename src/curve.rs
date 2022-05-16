// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This module provides an implementation of the STARK-friendly
//! cheetah curve over the sextic extension of the prime field Fp
//! of characteristic p = 2^64 - 2^32 + 1.

use core::{
    borrow::Borrow,
    cmp::Ordering,
    fmt,
    hash::{Hash, Hasher},
    iter::Sum,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::fp::reduce_u96;
use crate::fp::GENERATOR;
use crate::{Fp, Fp6, Scalar};

use crate::{MINUS_SHIFT_POINT_POW_256, SHIFT_POINT};

use crate::constants::ODD_MULTIPLES_BASEPOINT;
use crate::LookupTable;
use crate::NafLookupTable;

use group::{Curve, Group};
use rand_core::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "serialize")]
use serde::de::Visitor;
#[cfg(feature = "serialize")]
use serde::{self, Deserialize, Deserializer, Serialize, Serializer};

impl_binops_additive!(ProjectivePoint, AffinePoint);
impl_binops_additive_specify_output!(AffinePoint, ProjectivePoint, ProjectivePoint);

//  A = 1
/// B = u + 395
pub const B: Fp6 = Fp6 {
    c0: Fp(395),
    c1: Fp::one(),
    c2: Fp::zero(),
    c3: Fp::zero(),
    c4: Fp::zero(),
    c5: Fp::zero(),
};

const B3: Fp6 = (&B).mul(&Fp6::new([3, 0, 0, 0, 0, 0]));

/// An affine point
#[derive(Copy, Clone, Debug)]
pub struct AffinePoint {
    pub(crate) x: Fp6,
    pub(crate) y: Fp6,
    pub(crate) infinity: Choice,
}

impl Default for AffinePoint {
    fn default() -> AffinePoint {
        AffinePoint::identity()
    }
}

impl zeroize::DefaultIsZeroes for AffinePoint {}

impl fmt::Display for AffinePoint {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl Hash for AffinePoint {
    fn hash<H: Hasher>(&self, hasher: &mut H) {
        self.x.hash(hasher);
        self.y.hash(hasher);
    }
}

impl<'a> From<&'a ProjectivePoint> for AffinePoint {
    fn from(p: &'a ProjectivePoint) -> AffinePoint {
        let zinv = p.z.invert().unwrap_or(Fp6::zero());
        let x = p.x * zinv;
        let y = p.y * zinv;

        let tmp = AffinePoint {
            x,
            y,
            infinity: Choice::from(0u8),
        };

        AffinePoint::conditional_select(&tmp, &AffinePoint::identity(), zinv.ct_eq(&Fp6::zero()))
    }
}

impl From<ProjectivePoint> for AffinePoint {
    fn from(p: ProjectivePoint) -> AffinePoint {
        AffinePoint::from(&p)
    }
}

impl ConstantTimeEq for AffinePoint {
    fn ct_eq(&self, other: &Self) -> Choice {
        // The only cases in which two points are equal are
        // 1. infinity is set on both
        // 2. infinity is not set on both, and their coordinates are equal

        (self.infinity & other.infinity)
            | ((!self.infinity)
                & (!other.infinity)
                & self.x.ct_eq(&other.x)
                & self.y.ct_eq(&other.y))
    }
}

impl ConditionallySelectable for AffinePoint {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        AffinePoint {
            x: Fp6::conditional_select(&a.x, &b.x, choice),
            y: Fp6::conditional_select(&a.y, &b.y, choice),
            infinity: Choice::conditional_select(&a.infinity, &b.infinity, choice),
        }
    }
}

impl Eq for AffinePoint {}
impl PartialEq for AffinePoint {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl<'a> Neg for &'a AffinePoint {
    type Output = AffinePoint;

    #[inline]
    fn neg(self) -> AffinePoint {
        AffinePoint {
            x: self.x,
            y: Fp6::conditional_select(&-self.y, &Fp6::one(), self.infinity),
            infinity: self.infinity,
        }
    }
}

impl Neg for AffinePoint {
    type Output = AffinePoint;

    #[inline]
    fn neg(self) -> AffinePoint {
        -&self
    }
}

impl<'a, 'b> Add<&'b ProjectivePoint> for &'a AffinePoint {
    type Output = ProjectivePoint;

    #[inline]
    fn add(self, rhs: &'b ProjectivePoint) -> ProjectivePoint {
        rhs.add_mixed(self)
    }
}

impl<'a, 'b> Add<&'b AffinePoint> for &'a ProjectivePoint {
    type Output = ProjectivePoint;

    #[inline]
    fn add(self, rhs: &'b AffinePoint) -> ProjectivePoint {
        self.add_mixed(rhs)
    }
}

impl<'a, 'b> Sub<&'b ProjectivePoint> for &'a AffinePoint {
    type Output = ProjectivePoint;

    #[inline]
    fn sub(self, rhs: &'b ProjectivePoint) -> ProjectivePoint {
        self + (-rhs)
    }
}

impl<'a, 'b> Sub<&'b AffinePoint> for &'a ProjectivePoint {
    type Output = ProjectivePoint;

    #[inline]
    fn sub(self, rhs: &'b AffinePoint) -> ProjectivePoint {
        self + (-rhs)
    }
}

impl<T> Sum<T> for ProjectivePoint
where
    T: Borrow<ProjectivePoint>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Self::identity(), |acc, item| acc + item.borrow())
    }
}

impl AffinePoint {
    /// Returns the identity point in affine coordinates
    pub fn identity() -> AffinePoint {
        AffinePoint {
            x: Fp6::zero(),
            y: Fp6::one(),
            infinity: Choice::from(1u8),
        }
    }

    /// Returns the x coordinate of this AffinePoint
    pub const fn get_x(&self) -> Fp6 {
        self.x
    }

    /// Returns the y coordinate of this AffinePoint
    pub const fn get_y(&self) -> Fp6 {
        self.y
    }

    /// Computes a random `AffinePoint` element
    pub fn random(mut rng: impl RngCore) -> Self {
        ProjectivePoint::random(&mut rng).into()
    }

    /// Returns a fixed generator of the curve in affine coordinates
    /// The point has been generated from the Simplified Shallue-van
    /// de Woestijne-Ulas method for hashing to elliptic curves in
    /// Short Weierstrass form, applied on the integer value of the
    /// binary encoding of the string "Cheetah".
    ///
    /// See https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve
    /// for more.
    pub fn generator() -> AffinePoint {
        AffinePoint {
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
            infinity: Choice::from(0u8),
        }
    }

    /// Outputs a compress byte representation of this `AffinePoint` element
    pub fn to_compressed(&self) -> CompressedPoint {
        CompressedPoint::from_affine(self)
    }

    /// Outputs an uncompressed byte representation of this `AffinePoint` element
    /// It is twice larger than when calling `AffinePoint::to_compressed()`
    pub fn to_uncompressed(&self) -> UncompressedPoint {
        UncompressedPoint::from_affine(self)
    }

    /// Attempts to deserialize an uncompressed element.
    pub fn from_uncompressed(uncompressed_point: &UncompressedPoint) -> CtOption<Self> {
        Self::from_uncompressed_unchecked(uncompressed_point)
            .and_then(|p| CtOption::new(p, p.is_on_curve() & p.is_torsion_free()))
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is on the curve and not checking if it is in the correct subgroup.
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_uncompressed()` instead.
    pub fn from_uncompressed_unchecked(uncompressed_point: &UncompressedPoint) -> CtOption<Self> {
        uncompressed_point.to_affine()
    }

    /// Attempts to deserialize a compressed element.
    pub fn from_compressed(compressed_point: &CompressedPoint) -> CtOption<Self> {
        // We already know the point is on the curve because this is established
        // by the y-coordinate recovery procedure in from_compressed_unchecked().

        Self::from_compressed_unchecked(compressed_point)
            .and_then(|p| CtOption::new(p, p.is_torsion_free()))
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is in the correct subgroup.
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_compressed()` instead.
    pub fn from_compressed_unchecked(compressed_point: &CompressedPoint) -> CtOption<Self> {
        compressed_point.to_affine()
    }

    #[allow(unused)]
    /// Constructs an `AffinePoint` element without checking that it is a valid point.
    /// Assumes the coordinates do not represent the infinity point.
    pub fn from_raw_coordinates(elems: [Fp6; 2]) -> Self {
        AffinePoint {
            x: elems[0],
            y: elems[1],
            infinity: Choice::from(0u8),
        }
    }

    /// Returns true if this element is the identity (the point at infinity).
    #[inline]
    pub const fn is_identity(&self) -> Choice {
        self.infinity
    }

    /// Returns true if this point is on the curve. This should always return
    /// true unless an "unchecked" API was used.
    pub fn is_on_curve(&self) -> Choice {
        // y^2 - x^3 - x ?= B

        (self.y.square() - self.x * (self.x.square() + Fp6::one())).ct_eq(&B) | self.infinity
    }

    /// Performs an affine scalar multiplication from `by`
    /// given as byte representation of a `Scalar` element.
    #[inline]
    pub fn multiply(&self, by: &[u8; 32]) -> AffinePoint {
        ProjectivePoint::multiply(&self.into(), by).into()
    }

    /// Performs an affine scalar multiplication from `by`
    /// given as byte representation of a `Scalar` element.
    ///
    /// **This operation is variable time with respect
    /// to the scalar.** If the scalar is fixed,
    /// this operation is effectively constant time.
    #[inline]
    pub fn multiply_vartime(&self, by: &[u8; 32]) -> AffinePoint {
        ProjectivePoint::multiply_vartime(&self.into(), by).into()
    }

    /// Performs the affine sum [`by_lhs` * `self` + `by_rhs` * `rhs`] with `by_lhs`
    /// and `by_rhs` given as byte representations of `Scalar` elements.
    #[inline]
    pub fn multiply_double(
        &self,
        rhs: &AffinePoint,
        by_lhs: &[u8; 32],
        by_rhs: &[u8; 32],
    ) -> AffinePoint {
        ProjectivePoint::multiply_double(&self.into(), &rhs.into(), by_lhs, by_rhs).into()
    }

    /// Performs the affine sum [`by_lhs` * `self` + `by_rhs` * `rhs`] with `by_lhs`
    /// and `by_rhs` given as byte representations of `Scalar` elements.
    ///
    /// **This operation is variable time with respect
    /// to the scalars.** If the scalars are fixed,
    /// this operation is effectively constant time.
    #[inline]
    pub fn multiply_double_vartime(
        &self,
        rhs: &AffinePoint,
        by_lhs: &[u8; 32],
        by_rhs: &[u8; 32],
    ) -> AffinePoint {
        ProjectivePoint::multiply_double_vartime(&self.into(), &rhs.into(), by_lhs, by_rhs).into()
    }

    /// Performs the affine sum [`by_lhs` * `self` + `by_rhs` * `rhs`] with `by_lhs`
    /// and `by_rhs` given as byte representations of `Scalar` elements.
    ///
    /// **This operation is variable time with respect
    /// to the scalars.** If the scalars are fixed,
    /// this operation is effectively constant time.
    #[inline]
    pub fn multiply_double_with_basepoint_vartime(
        &self,
        by_self: &[u8; 32],
        by_basepoint: &[u8; 32],
    ) -> AffinePoint {
        ProjectivePoint::multiply_double_with_basepoint_vartime(&self.into(), by_self, by_basepoint)
            .into()
    }

    /// Multiplies by the curve cofactor
    pub fn clear_cofactor(&self) -> AffinePoint {
        let point: ProjectivePoint = self.into();

        point.clear_cofactor().into()
    }

    /// Returns true if this point is free of an $h$-torsion component.
    /// This should always return true unless an "unchecked" API was used.
    pub fn is_torsion_free(&self) -> Choice {
        let point: ProjectivePoint = self.into();

        point.is_torsion_free()
    }
}

/// A projective point
#[derive(Copy, Clone, Debug)]
pub struct ProjectivePoint {
    pub(crate) x: Fp6,
    pub(crate) y: Fp6,
    pub(crate) z: Fp6,
}

impl Default for ProjectivePoint {
    fn default() -> ProjectivePoint {
        ProjectivePoint::identity()
    }
}

impl zeroize::DefaultIsZeroes for ProjectivePoint {}

impl fmt::Display for ProjectivePoint {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl Hash for ProjectivePoint {
    fn hash<H: Hasher>(&self, hasher: &mut H) {
        self.to_affine().hash(hasher);
    }
}

impl<'a> From<&'a AffinePoint> for ProjectivePoint {
    fn from(p: &'a AffinePoint) -> ProjectivePoint {
        ProjectivePoint {
            x: p.x,
            y: p.y,
            z: Fp6::conditional_select(&Fp6::one(), &Fp6::zero(), p.infinity),
        }
    }
}

impl From<AffinePoint> for ProjectivePoint {
    fn from(p: AffinePoint) -> ProjectivePoint {
        ProjectivePoint::from(&p)
    }
}

impl ConstantTimeEq for ProjectivePoint {
    fn ct_eq(&self, other: &Self) -> Choice {
        // Is (xz, yz, z) equal to (x'z', y'z', z') when converted to affine?

        let x1 = self.x * other.z;
        let x2 = other.x * self.z;

        let y1 = self.y * other.z;
        let y2 = other.y * self.z;

        let self_is_zero = self.z.is_zero();
        let other_is_zero = other.z.is_zero();

        (self_is_zero & other_is_zero) // Both point at infinity
            | ((!self_is_zero) & (!other_is_zero) & x1.ct_eq(&x2) & y1.ct_eq(&y2))
        // Neither point at infinity, coordinates are the same
    }
}

impl ConditionallySelectable for ProjectivePoint {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        ProjectivePoint {
            x: Fp6::conditional_select(&a.x, &b.x, choice),
            y: Fp6::conditional_select(&a.y, &b.y, choice),
            z: Fp6::conditional_select(&a.z, &b.z, choice),
        }
    }
}

impl Eq for ProjectivePoint {}
impl PartialEq for ProjectivePoint {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl<'a> Neg for &'a ProjectivePoint {
    type Output = ProjectivePoint;

    #[inline]
    fn neg(self) -> ProjectivePoint {
        self.neg()
    }
}

impl Neg for ProjectivePoint {
    type Output = ProjectivePoint;

    #[inline]
    fn neg(self) -> ProjectivePoint {
        -&self
    }
}

impl<'a, 'b> Add<&'b ProjectivePoint> for &'a ProjectivePoint {
    type Output = ProjectivePoint;

    #[inline]
    fn add(self, rhs: &'b ProjectivePoint) -> ProjectivePoint {
        self.add(rhs)
    }
}

impl<'a, 'b> Sub<&'b ProjectivePoint> for &'a ProjectivePoint {
    type Output = ProjectivePoint;

    #[inline]
    fn sub(self, rhs: &'b ProjectivePoint) -> ProjectivePoint {
        self + (-rhs)
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a ProjectivePoint {
    type Output = ProjectivePoint;

    fn mul(self, other: &'b Scalar) -> Self::Output {
        self.multiply(&other.to_bytes())
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a AffinePoint {
    type Output = ProjectivePoint;

    fn mul(self, other: &'b Scalar) -> Self::Output {
        ProjectivePoint::from(self).multiply(&other.to_bytes())
    }
}

impl_binops_additive!(ProjectivePoint, ProjectivePoint);
impl_binops_multiplicative!(ProjectivePoint, Scalar);
impl_binops_multiplicative_mixed!(AffinePoint, Scalar, ProjectivePoint);

#[inline(always)]
const fn mul_by_3b(a: &Fp6) -> Fp6 {
    let aa = (&a.c0).mul(&B3.c0).0 as u128;
    let ab = a.c0.0 as u128;
    let ab = ab + (ab << 1);

    let ba = (&a.c1).mul(&B3.c0).0 as u128;
    let bb = a.c1.0 as u128;
    let bb = bb + (bb << 1);

    let ca = (&a.c2).mul(&B3.c0).0 as u128;
    let cb = a.c2.0 as u128;
    let cb = cb + (cb << 1);

    let da = (&a.c3).mul(&B3.c0).0 as u128;
    let db = a.c3.0 as u128;
    let db = db + (db << 1);

    let ea = (&a.c4).mul(&B3.c0).0 as u128;
    let eb = a.c4.0 as u128;
    let eb = eb + (eb << 1);

    let fa = (&a.c5).mul(&B3.c0).0 as u128;
    let fb = a.c5.0 as u128;
    let fb = fb + (fb << 1);

    let c0 = fb * GENERATOR.0 as u128;
    let c0 = c0 + aa;
    let c0 = Fp(reduce_u96(c0));

    let c1 = ab + ba;
    let c1 = Fp(reduce_u96(c1));

    let c2 = bb + ca;
    let c2 = Fp(reduce_u96(c2));

    let c3 = cb + da;
    let c3 = Fp(reduce_u96(c3));

    let c4 = db + ea;
    let c4 = Fp(reduce_u96(c4));

    let c5 = eb + fa;
    let c5 = Fp(reduce_u96(c5));

    Fp6 {
        c0,
        c1,
        c2,
        c3,
        c4,
        c5,
    }
}

impl ProjectivePoint {
    /// Returns the identity of the group: the point at infinity.
    pub const fn identity() -> ProjectivePoint {
        ProjectivePoint {
            x: Fp6::zero(),
            y: Fp6::one(),
            z: Fp6::zero(),
        }
    }

    /// Returns the x coordinate of this ProjectivePoint
    pub const fn get_x(&self) -> Fp6 {
        self.x
    }

    /// Returns the y coordinate of this ProjectivePoint
    pub const fn get_y(&self) -> Fp6 {
        self.y
    }

    /// Returns the z coordinate of this ProjectivePoint
    pub const fn get_z(&self) -> Fp6 {
        self.z
    }

    /// Computes a random `ProjectivePoint` element
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
                let p = ProjectivePoint::from(p.unwrap()).clear_cofactor();

                if bool::from(!p.is_identity()) {
                    return p;
                }
            }
        }
    }

    /// Returns a fixed generator of the curve in projective coordinates
    /// The point has been generated from the Simplified Shallue-van
    /// de Woestijne-Ulas method for hashing to elliptic curves in
    /// Short Weierstrass form, applied on the integer value of the
    /// binary encoding of the string "Cheetah".
    ///
    /// See https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve
    /// for more.
    pub const fn generator() -> ProjectivePoint {
        ProjectivePoint {
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

    /// Outputs a compress byte representation of this `ProjectivePoint` element
    pub fn to_compressed(&self) -> CompressedPoint {
        AffinePoint::from(self).to_compressed()
    }

    /// Outputs an uncompressed byte representation of this `ProjectivePoint` element
    /// It is twice larger than when calling `ProjectivePoint::to_uncompress()`
    pub fn to_uncompressed(&self) -> UncompressedPoint {
        AffinePoint::from(self).to_uncompressed()
    }

    /// Attempts to deserialize an uncompressed element.
    pub fn from_uncompressed(uncompressed_point: &UncompressedPoint) -> CtOption<Self> {
        AffinePoint::from_uncompressed(uncompressed_point)
            .and_then(|p| CtOption::new(ProjectivePoint::from(p), 1.into()))
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is on the curve and not checking if it is in the correct subgroup.
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_uncompressed()` instead.
    pub fn from_uncompressed_unchecked(uncompressed_point: &UncompressedPoint) -> CtOption<Self> {
        AffinePoint::from_uncompressed_unchecked(uncompressed_point)
            .and_then(|p| CtOption::new(ProjectivePoint::from(p), 1.into()))
    }

    /// Attempts to deserialize a compressed element.
    pub fn from_compressed(compressed_point: &CompressedPoint) -> CtOption<Self> {
        AffinePoint::from_compressed(compressed_point).map(ProjectivePoint::from)
    }

    #[allow(unused)]
    /// Constructs a `ProjectivePoint` element without checking that it is a valid point.
    pub const fn from_raw_coordinates(elems: [Fp6; 3]) -> Self {
        ProjectivePoint {
            x: elems[0],
            y: elems[1],
            z: elems[2],
        }
    }

    /// Computes `n` iterated doubling of this point.
    pub fn double_multi(&self, n: u32) -> ProjectivePoint {
        assert!(n >= 1);
        let mut output = self.double();

        for _ in 1..n {
            output = output.double();
        }

        output
    }

    /// Computes the doubling of this point.
    pub fn double(&self) -> ProjectivePoint {
        // Use formula given in Handbook of Elliptic and Hyperelliptic Curve Cryptography, part 13.2

        let x2 = self.x.square();
        let z2 = self.z.square();

        let a = x2.double();
        let a = (&a).add(&x2);
        let a = (&a).add(&z2);

        let b = (&self.y).mul(&self.z);
        let t3 = (&self.y).mul(&b);
        let c = (&t3).mul(&self.x);

        let d = c.double();
        let c4 = d.double();
        let d = c4.double();

        let t = a.square();
        let d = (&t).sub(&d);

        let x3 = (&b).mul(&d);
        let x3 = x3.double();

        let y3 = (&c4).sub(&d);
        let y3 = (&a).mul(&y3);
        let t3 = t3.square();

        let t3 = t3.double();
        let t3 = t3.double();
        let t3 = t3.double();

        let y3 = (&y3).sub(&t3);
        let z3 = b.square();
        let z3 = (&z3).mul(&b);

        let z3 = z3.double();
        let z3 = z3.double();
        let z3 = z3.double();

        let tmp = ProjectivePoint {
            x: x3,
            y: y3,
            z: z3,
        };

        // The calculation above fails for doubling the infinity point,
        // hence we do a final conditional selection.
        ProjectivePoint::conditional_select(&tmp, &ProjectivePoint::identity(), self.is_identity())
    }

    /// Computes the negation of a point in projective coordinates
    #[inline]
    pub const fn neg(&self) -> Self {
        Self {
            x: self.x,
            y: (&self.y).neg(),
            z: self.z,
        }
    }

    /// Adds this point to another point.
    pub const fn add(&self, rhs: &ProjectivePoint) -> ProjectivePoint {
        // Algorithm 1, https://eprint.iacr.org/2015/1060.pdf

        let t0 = (&self.x).mul(&rhs.x);
        let t1 = (&self.y).mul(&rhs.y);
        let t2 = (&self.z).mul(&rhs.z);

        let t3 = (&self.x).add(&self.y);
        let t4 = (&rhs.x).add(&rhs.y);
        let t3 = (&t3).mul(&t4);

        let t4 = (&t0).add(&t1);
        let t3 = (&t3).sub(&t4);
        let t4 = (&self.x).add(&self.z);

        let t5 = (&rhs.x).add(&rhs.z);
        let t4 = (&t4).mul(&t5);
        let t5 = (&t0).add(&t2);

        let t4 = (&t4).sub(&t5);
        let t5 = (&self.y).add(&self.z);
        let x3 = (&rhs.y).add(&rhs.z);

        let t5 = (&t5).mul(&x3);
        let x3 = (&t1).add(&t2);
        let t5 = (&t5).sub(&x3);

        let x3 = mul_by_3b(&t2);
        let z3 = (&x3).add(&t4);

        let x3 = (&t1).sub(&z3);
        let z3 = (&t1).add(&z3);
        let y3 = (&x3).mul(&z3);

        let t1 = t0.double();
        let t1 = (&t1).add(&t0);

        let t4 = mul_by_3b(&t4);
        let t1 = (&t1).add(&t2);
        let t2 = (&t0).sub(&t2);

        let t4 = (&t4).add(&t2);
        let t0 = (&t1).mul(&t4);

        let y3 = (&y3).add(&t0);
        let t0 = (&t5).mul(&t4);
        let x3 = (&t3).mul(&x3);

        let x3 = (&x3).sub(&t0);
        let t0 = (&t3).mul(&t1);
        let z3 = (&t5).mul(&z3);

        let z3 = (&z3).add(&t0);

        ProjectivePoint {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    /// Adds this point to another point in the affine model.
    pub fn add_mixed(&self, rhs: &AffinePoint) -> ProjectivePoint {
        // Algorithm 2, https://eprint.iacr.org/2015/1060.pdf

        let t0 = (&self.x).mul(&rhs.x);
        let t1 = (&self.y).mul(&rhs.y);
        let t3 = (&rhs.x).add(&rhs.y);

        let t4 = (&self.x).add(&self.y);
        let t3 = (&t3).mul(&t4);
        let t4 = (&t0).add(&t1);

        let t3 = (&t3).sub(&t4);
        let t4 = (&rhs.x).mul(&self.z);
        let t4 = (&t4).add(&self.x);

        let t5 = (&rhs.y).mul(&self.z);
        let t5 = (&t5).add(&self.y);

        let x3 = mul_by_3b(&self.z);
        let z3 = (&x3).add(&t4);
        let x3 = (&t1).sub(&z3);

        let z3 = (&t1).add(&z3);
        let y3 = (&x3).mul(&z3);
        let t1 = t0.double();

        let t1 = (&t1).add(&t0);
        let t4 = mul_by_3b(&t4);

        let t1 = (&t1).add(&self.z);
        let t2 = (&t0).sub(&self.z);

        let t4 = (&t4).add(&t2);
        let t0 = (&t1).mul(&t4);
        let y3 = (&y3).add(&t0);

        let t0 = (&t5).mul(&t4);
        let x3 = (&t3).mul(&x3);
        let x3 = (&x3).sub(&t0);

        let t0 = (&t3).mul(&t1);
        let z3 = (&t5).mul(&z3);
        let z3 = (&z3).add(&t0);

        let tmp = ProjectivePoint {
            x: x3,
            y: y3,
            z: z3,
        };

        ProjectivePoint::conditional_select(&tmp, self, rhs.is_identity())
    }

    /// Adds this point to another point in the affine model. This formulae is faster than
    /// `ProjectivePoint::add_mixed` but is not complete.
    /// **This is dangerous to call unless you know that the points to be added are not
    /// identical, opposite of each other, and neither of them is the identity point; otherwise,
    /// API invariants may be broken.** Please consider using `add_mixed()` instead.
    pub fn add_mixed_unchecked(&self, rhs: &AffinePoint) -> ProjectivePoint {
        // Use formula given in Handbook of Elliptic and Hyperelliptic Curve Cryptography, part 13.2

        let t0 = (&rhs.y).mul(&self.z);
        let t0 = (&t0).sub(&self.y);

        let t1 = (&rhs.x).mul(&self.z);
        let t1 = (&t1).sub(&self.x);

        let t2 = t0.square();
        let t2 = (&t2).mul(&self.z);
        let t3 = t1.square();
        let t4 = (&t1).mul(&t3);
        let t2 = (&t2).sub(&t4);
        let t5 = (&t3).mul(&self.x);
        let t6 = t5.double();
        let t2 = (&t2).sub(&t6);

        let x3 = (&t1).mul(&t2);

        let y3 = (&t5).sub(&t2);
        let y3 = (&y3).mul(&t0);
        let t7 = (&t4).mul(&self.y);
        let y3 = (&y3).sub(&t7);

        let z3 = (&t4).mul(&self.z);

        let tmp = ProjectivePoint {
            x: x3,
            y: y3,
            z: z3,
        };

        ProjectivePoint::conditional_select(&tmp, self, rhs.is_identity())
    }

    /// Performs a projective scalar multiplication from `by`
    /// given as byte representation of a `Scalar` element
    pub fn multiply(&self, by: &[u8; 32]) -> ProjectivePoint {
        let table = LookupTable::<16>::from(self);
        let digits = Scalar::bytes_to_radix_16(by);

        let mut acc = SHIFT_POINT;
        for i in (0..64).rev() {
            acc = acc.double_multi(4);
            acc = acc.add_mixed_unchecked(&table.get_point(digits[i]));
        }

        acc.add_mixed_unchecked(&MINUS_SHIFT_POINT_POW_256)
    }

    /// Performs a projective scalar multiplication from `by`
    /// given as byte representation of a `Scalar` element.
    ///
    /// **This operation is variable time with respect
    /// to the scalar.** If the scalar is fixed,
    /// this operation is effectively constant time.
    pub fn multiply_vartime(&self, by: &[u8; 32]) -> ProjectivePoint {
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
        let mut acc = ProjectivePoint::identity();

        for j in (0..i + 1).rev() {
            acc = acc.double();

            match digits[j].cmp(&0) {
                Ordering::Greater => acc += &table.get_point(digits[j] as usize),
                Ordering::Less => acc -= &table.get_point(-digits[j] as usize),
                Ordering::Equal => (),
            };
        }

        acc
    }

    /// Performs the projective sum [`by_lhs` * `self` + `by_rhs` * `rhs`] with
    /// `by_lhs` and `by_rhs` given as byte representations of `Scalar` elements.
    pub fn multiply_double(
        &self,
        rhs: &ProjectivePoint,
        by_lhs: &[u8; 32],
        by_rhs: &[u8; 32],
    ) -> ProjectivePoint {
        let table_lhs = LookupTable::<16>::from(self);
        let table_rhs = LookupTable::<16>::from(rhs);
        let digits_lhs = Scalar::bytes_to_radix_16(by_lhs);
        let digits_rhs = Scalar::bytes_to_radix_16(by_rhs);

        let mut acc = SHIFT_POINT;
        for i in (0..64).rev() {
            acc = acc.double_multi(4);
            acc = acc.add_mixed_unchecked(&table_lhs.get_point(digits_lhs[i]));
            acc = acc.add_mixed_unchecked(&table_rhs.get_point(digits_rhs[i]));
        }

        acc.add_mixed_unchecked(&MINUS_SHIFT_POINT_POW_256)
    }

    /// Performs the projective sum [`by_lhs` * `self` + `by_rhs` * `rhs`] with
    /// `by_lhs` and `by_rhs` given as byte representations of `Scalar` elements.
    ///
    /// **This operation is variable time with respect
    /// to the scalars.** If the scalars are fixed,
    /// this operation is effectively constant time.
    pub fn multiply_double_vartime(
        &self,
        rhs: &ProjectivePoint,
        by_lhs: &[u8; 32],
        by_rhs: &[u8; 32],
    ) -> ProjectivePoint {
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
        let mut acc = ProjectivePoint::identity();

        for j in (0..i + 1).rev() {
            acc = acc.double();

            match by_lhs_digits[j].cmp(&0) {
                Ordering::Greater => acc += &table_self.get_point(by_lhs_digits[j] as usize),
                Ordering::Less => acc -= &table_self.get_point(-by_lhs_digits[j] as usize),
                Ordering::Equal => (),
            };

            match by_rhs_digits[j].cmp(&0) {
                Ordering::Greater => acc += &table_rhs.get_point(by_rhs_digits[j] as usize),
                Ordering::Less => acc -= &table_rhs.get_point(-by_rhs_digits[j] as usize),
                Ordering::Equal => (),
            };
        }

        acc
    }

    // TODO: It would be nice to have a constant-time variant of the method below,
    // though this may be tricky for accessing the point in the table.

    /// Performs the projective sum [`by_lhs` * `self` + `by_rhs` * `g`] with
    /// `by_lhs` and `by_rhs` given as byte representations of `Scalar` elements
    /// and `g` being the curve basepoint.
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
    ) -> ProjectivePoint {
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
        let mut acc = ProjectivePoint::identity();

        for j in (0..i + 1).rev() {
            acc = acc.double();

            match by_self_digits[j].cmp(&0) {
                Ordering::Greater => acc += &table_self.get_point(by_self_digits[j] as usize),
                Ordering::Less => acc -= &table_self.get_point(-by_self_digits[j] as usize),
                Ordering::Equal => (),
            };

            match by_basepoint_digits[j].cmp(&0) {
                Ordering::Greater => {
                    acc += &table_basepoint.get_point(by_basepoint_digits[j] as usize)
                }
                Ordering::Less => {
                    acc -= &table_basepoint.get_point(-by_basepoint_digits[j] as usize)
                }
                Ordering::Equal => (),
            };
        }

        acc
    }

    /// Multiplies by the curve cofactor
    pub fn clear_cofactor(&self) -> ProjectivePoint {
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

    /// Converts a batch of `ProjectivePoint` elements into `AffinePoint` elements. This
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

            // Cancel out z-coordinate in denominator of `acc`
            acc = Fp6::conditional_select(&(acc * p.z), &acc, skip);

            // Set the coordinates to the correct value
            q.x = p.x * tmp;
            q.y = p.y * tmp;
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
        // Y^2 Z = X^3 + X Z^2 + b Z^3

        (self.y.square() * self.z).ct_eq(
            &(self.x.square() * self.x + self.x * self.z.square() + self.z.square() * self.z * B),
        ) | (self.z.is_zero())
    }
}

/// A compressed point, storing the `x` coordinate of
/// a point, along an extra byte storing metadata
/// to be used for decompression.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct CompressedPoint(pub [u8; 49]);

impl Default for CompressedPoint {
    fn default() -> Self {
        AffinePoint::identity().to_compressed()
    }
}

impl ConstantTimeEq for CompressedPoint {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.as_bytes().ct_eq(other.as_bytes())
    }
}

impl CompressedPoint {
    /// Converts an `AffinePoint` to a `CompressedPoint`
    fn from_affine(point: &AffinePoint) -> Self {
        let mut bytes = [0u8; 49];

        // Strictly speaking, point.x is zero already when point.infinity is true, but
        // to guard against implementation mistakes we do not assume this.
        bytes[0..48].copy_from_slice(
            &Fp6::conditional_select(&point.x, &Fp6::zero(), point.infinity).to_bytes(),
        );

        // Is this point at infinity? If so, set the most significant bit of the last byte.
        bytes[48] |= u8::conditional_select(&0u8, &(1u8 << 7), point.infinity);

        // Is the y-coordinate the lexicographically largest of the two associated with the
        // x-coordinate? If so, set the second most significant bit of the last bytes so long
        // as this is not the point at infinity.
        bytes[48] |= u8::conditional_select(
            &0u8,
            &(1u8 << 6),
            (!point.infinity) & point.y.lexicographically_largest(),
        );

        Self(bytes)
    }

    /// Attempts to convert a `CompressedPoint` to an `AffinePoint`
    /// The resulting point is ensured to be on the curve, but is
    /// not necessarily on the prime order subgroup.
    fn to_affine(self) -> CtOption<AffinePoint> {
        // Obtain the two flags stored in the last byte of `self`
        let infinity_flag_set = Choice::from((self.0[48] >> 7) & 1);
        let sort_flag_set = Choice::from((self.0[48] >> 6) & 1);

        // Attempt to obtain the x-coordinate
        let x = {
            let mut tmp = [0; 48];
            tmp.copy_from_slice(&self.0[0..48]);

            Fp6::from_bytes(&tmp)
        };

        x.and_then(|x| {
            // If the infinity flag is set, return the value assuming
            // the x-coordinate is zero and the sort flag is not set.
            //
            // Otherwise, return a recovered point (assuming the correct
            // y-coordinate can be found) so long as the infinity flag
            // was not set.
            CtOption::new(
                AffinePoint::identity(),
                infinity_flag_set & // Infinity flag should be set
                    (!sort_flag_set) & // Sort flag should not be set
                    x.is_zero(), // The x-coordinate should be zero
            )
            .or_else(|| {
                // Recover a y-coordinate given x by y = sqrt(x^3 + x + B)
                ((x.square() * x) + x + B).sqrt().and_then(|y| {
                    // Switch to the correct y-coordinate if necessary.
                    let y = Fp6::conditional_select(
                        &y,
                        &-y,
                        y.lexicographically_largest() ^ sort_flag_set,
                    );

                    CtOption::new(
                        AffinePoint {
                            x,
                            y,
                            infinity: infinity_flag_set,
                        },
                        !infinity_flag_set, // Infinity flag should not be set
                    )
                })
            })
        })
    }

    /// Copies the bytes of this `CompressedPoint`.
    pub fn to_bytes(&self) -> [u8; 49] {
        self.0
    }

    /// Views this `CompressedPoint` as an array of bytes.
    pub fn as_bytes(&self) -> &[u8; 49] {
        &self.0
    }

    /// Interprets the provided bytes as a `CompressedPoint`.
    /// This does not check the validity of the input, and may result
    /// in failure when calling `to_affine()`.
    pub fn from_bytes(bytes: &[u8; 49]) -> Self {
        Self(*bytes)
    }
}

/// A uncompressed point, storing the `x` and `y` coordinates
/// of a point, along an extra byte storing metadata to be
/// used for decompression.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct UncompressedPoint(pub [u8; 97]);

impl Default for UncompressedPoint {
    fn default() -> Self {
        AffinePoint::identity().to_uncompressed()
    }
}

impl ConstantTimeEq for UncompressedPoint {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.as_bytes().ct_eq(other.as_bytes())
    }
}

impl UncompressedPoint {
    /// Converts an `AffinePoint` to an `UncompressedPoint`
    fn from_affine(point: &AffinePoint) -> Self {
        let mut bytes = [0; 97];

        bytes[0..48].copy_from_slice(
            &Fp6::conditional_select(&point.x, &Fp6::zero(), point.infinity).to_bytes()[..],
        );
        bytes[48..96].copy_from_slice(
            &Fp6::conditional_select(&point.y, &Fp6::zero(), point.infinity).to_bytes()[..],
        );

        // Is this point at infinity? If so, set the most significant bit of the last byte.
        bytes[96] |= u8::conditional_select(&0u8, &(1u8 << 7), point.infinity);

        Self(bytes)
    }

    /// Attempts to convert an `UncompressedPoint` to an `AffinePoint`
    /// The resulting point is ensured to be on the curve, but is
    /// not necessarily on the prime order subgroup.
    fn to_affine(self) -> CtOption<AffinePoint> {
        // Obtain the two flags
        let infinity_flag_set = Choice::from((self.0[96] >> 7) & 1);

        // Attempt to obtain the x-coordinate
        let x = {
            let mut tmp = [0; 48];
            tmp.copy_from_slice(&self.0[0..48]);

            Fp6::from_bytes(&tmp)
        };

        // Attempt to obtain the y-coordinate
        let y = {
            let mut tmp = [0; 48];
            tmp.copy_from_slice(&self.0[48..96]);

            Fp6::from_bytes(&tmp)
        };

        x.and_then(|x| {
            y.and_then(|y| {
                let p = AffinePoint::conditional_select(
                    &AffinePoint {
                        x,
                        y,
                        infinity: infinity_flag_set,
                    },
                    &AffinePoint::identity(),
                    infinity_flag_set,
                );

                CtOption::new(
                    p,
                    // If the infinity flag is set, the x and y coordinates should have been zero.
                    (!infinity_flag_set) | (infinity_flag_set & x.is_zero() & y.is_zero()),
                )
            })
        })
    }

    /// Copies the bytes of this `UncompressedPoint`.
    pub fn to_bytes(&self) -> [u8; 97] {
        self.0
    }

    /// Views this `UncompressedPoint` as an array of bytes.
    pub fn as_bytes(&self) -> &[u8; 97] {
        &self.0
    }

    /// Interprets the provided bytes as a `UncompressedPoint`.
    /// This does not check the validity of the input, and may result
    /// in failure when calling `to_affine()`.
    pub fn from_bytes(bytes: &[u8; 97]) -> Self {
        Self(*bytes)
    }
}

// GROUP TRAITS IMPLEMENTATION
// ================================================================================================

impl Group for ProjectivePoint {
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

impl Curve for ProjectivePoint {
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
impl Serialize for CompressedPoint {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        use serde::ser::SerializeTuple;
        let mut tup = serializer.serialize_tuple(49)?;
        for byte in self.0.iter() {
            tup.serialize_element(byte)?;
        }
        tup.end()
    }
}

#[cfg(feature = "serialize")]
impl Serialize for UncompressedPoint {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        use serde::ser::SerializeTuple;
        let mut tup = serializer.serialize_tuple(97)?;
        for byte in self.0.iter() {
            tup.serialize_element(byte)?;
        }
        tup.end()
    }
}

#[cfg(feature = "serialize")]
impl Serialize for AffinePoint {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        self.to_compressed().serialize(serializer)
    }
}

#[cfg(feature = "serialize")]
impl Serialize for ProjectivePoint {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        self.to_compressed().serialize(serializer)
    }
}

#[cfg(feature = "serialize")]
impl<'de> Deserialize<'de> for CompressedPoint {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct CompressedPointVisitor;

        impl<'de> Visitor<'de> for CompressedPointVisitor {
            type Value = CompressedPoint;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("49 bytes of data")
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<CompressedPoint, A::Error>
            where
                A: serde::de::SeqAccess<'de>,
            {
                let mut bytes = [0u8; 49];
                for (i, byte) in bytes.iter_mut().enumerate() {
                    *byte = seq
                        .next_element()?
                        .ok_or_else(|| serde::de::Error::invalid_length(i, &"expected 49 bytes"))?;
                }

                Ok(CompressedPoint(bytes))
            }
        }

        deserializer.deserialize_tuple(49, CompressedPointVisitor)
    }
}

#[cfg(feature = "serialize")]
impl<'de> Deserialize<'de> for UncompressedPoint {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct UncompressedPointVisitor;

        impl<'de> Visitor<'de> for UncompressedPointVisitor {
            type Value = UncompressedPoint;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("97 bytes of data")
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<UncompressedPoint, A::Error>
            where
                A: serde::de::SeqAccess<'de>,
            {
                let mut bytes = [0u8; 97];
                for (i, byte) in bytes.iter_mut().enumerate() {
                    *byte = seq
                        .next_element()?
                        .ok_or_else(|| serde::de::Error::invalid_length(i, &"expected 97 bytes"))?;
                }

                Ok(UncompressedPoint(bytes))
            }
        }

        deserializer.deserialize_tuple(97, UncompressedPointVisitor)
    }
}

#[cfg(feature = "serialize")]
impl<'de> Deserialize<'de> for AffinePoint {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let compressed_point = CompressedPoint::deserialize(deserializer)?;
        let p = AffinePoint::from_compressed(&compressed_point);
        if bool::from(p.is_some()) {
            Ok(p.unwrap())
        } else {
            Err(serde::de::Error::custom("decompression failed"))
        }
    }
}

#[cfg(feature = "serialize")]
impl<'de> Deserialize<'de> for ProjectivePoint {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let compressed_point = CompressedPoint::deserialize(deserializer)?;
        let p = ProjectivePoint::from_compressed(&compressed_point);
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

    use crate::{BasePointTable, BASEPOINT_TABLE};

    #[test]
    fn test_is_on_curve() {
        assert!(bool::from(AffinePoint::identity().is_on_curve()));
        assert!(bool::from(AffinePoint::generator().is_on_curve()));
        assert!(bool::from(ProjectivePoint::identity().is_on_curve()));
        assert!(bool::from(ProjectivePoint::generator().is_on_curve()));

        let z = Fp6 {
            c0: Fp(0x79b230ab69fdb493),
            c1: Fp(0xfe0bff52056ea77b),
            c2: Fp(0x68ae55ee9963bb1d),
            c3: Fp(0xaed3160aa9927b96),
            c4: Fp(0x57f9b8ca1372d78e),
            c5: Fp(0xcfb89c494391a555),
        };

        let gen = AffinePoint::generator();
        let mut test = ProjectivePoint {
            x: gen.x * z,
            y: gen.y * z,
            z,
        };

        assert!(bool::from(test.is_on_curve()));

        test.x = z;
        assert!(!bool::from(test.is_on_curve()));
    }

    #[test]
    #[allow(clippy::eq_op)]
    fn test_affine_point_equality() {
        let a = AffinePoint::generator();
        let b = AffinePoint::identity();
        let c = AffinePoint::default();

        assert!(a == a);
        assert!(b == b);
        assert!(b == c);
        assert!(a != b);
        assert!(b != a);

        assert!(bool::from(b.is_identity()));
        assert!(!bool::from(a.ct_eq(&b)));
    }

    #[test]
    #[allow(clippy::eq_op)]
    fn test_projective_point_equality() {
        let a = ProjectivePoint::generator();
        let b = ProjectivePoint::identity();
        let c = ProjectivePoint::default();

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

        let mut c = ProjectivePoint {
            x: a.x * z,
            y: a.y * z,
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
    fn test_projective_to_affine() {
        let a = ProjectivePoint::generator();
        let b = ProjectivePoint::identity();

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

        let c = ProjectivePoint {
            x: a.x * z,
            y: a.y * z,
            z,
        };

        assert_eq!(AffinePoint::from(c), AffinePoint::generator());
    }

    #[test]
    fn test_affine_to_projective() {
        let a = AffinePoint::generator();
        let b = AffinePoint::identity();

        assert!(bool::from(ProjectivePoint::from(a).is_on_curve()));
        assert!(!bool::from(ProjectivePoint::from(a).is_identity()));
        assert!(bool::from(ProjectivePoint::from(b).is_on_curve()));
        assert!(bool::from(ProjectivePoint::from(b).is_identity()));
    }

    #[test]
    fn test_doubling() {
        {
            let tmp = ProjectivePoint::identity().double();
            assert!(bool::from(tmp.is_identity()));
            assert!(bool::from(tmp.is_on_curve()));
        }
        {
            let tmp = ProjectivePoint::generator().double();
            assert!(!bool::from(tmp.is_identity()));
            assert!(bool::from(tmp.is_on_curve()));

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
    fn test_projective_addition() {
        {
            let a = ProjectivePoint::identity();
            let b = ProjectivePoint::identity();
            let c = a + b;
            assert!(bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
        }
        {
            let a = ProjectivePoint::identity();
            let mut b = ProjectivePoint::generator();
            {
                let z = Fp6 {
                    c0: Fp(0xa3c62f9770336022),
                    c1: Fp(0x69a1531152c92fc1),
                    c2: Fp(0x1636d9b90656a08a),
                    c3: Fp(0x635289066593aaf6),
                    c4: Fp(0x1e2178e3c54f6682),
                    c5: Fp(0xe15458e6d847f393),
                };

                b = ProjectivePoint {
                    x: b.x * z,
                    y: b.y * z,
                    z,
                };
            }
            let c = a + b;
            assert!(!bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
            assert!(c == ProjectivePoint::generator());
        }
        {
            let a = ProjectivePoint::generator().double_multi(2); // 4P
            let b = ProjectivePoint::generator().double(); // 2P
            let c = a + b;

            let mut d = ProjectivePoint::generator();
            for _ in 0..5 {
                d += ProjectivePoint::generator();
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
            let b = ProjectivePoint::identity();
            let c = a + b;
            assert!(bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
        }
        {
            let a = AffinePoint::identity();
            let mut b = ProjectivePoint::generator();
            {
                let z = Fp6 {
                    c0: Fp(0xa3c62f9770336022),
                    c1: Fp(0x69a1531152c92fc1),
                    c2: Fp(0x1636d9b90656a08a),
                    c3: Fp(0x635289066593aaf6),
                    c4: Fp(0x1e2178e3c54f6682),
                    c5: Fp(0xe15458e6d847f393),
                };

                b = ProjectivePoint {
                    x: b.x * z,
                    y: b.y * z,
                    z,
                };
            }
            let c = a + b;
            assert!(!bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
            assert!(c == ProjectivePoint::generator());
        }
        {
            let a = AffinePoint::identity();
            let mut b = ProjectivePoint::generator();
            {
                let z = Fp6 {
                    c0: Fp(0xa3c62f9770336022),
                    c1: Fp(0x69a1531152c92fc1),
                    c2: Fp(0x1636d9b90656a08a),
                    c3: Fp(0x635289066593aaf6),
                    c4: Fp(0x1e2178e3c54f6682),
                    c5: Fp(0xe15458e6d847f393),
                };

                b = ProjectivePoint {
                    x: b.x * z,
                    y: b.y * z,
                    z,
                };
            }
            let c = b + a;
            assert!(!bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
            assert!(c == ProjectivePoint::generator());
        }
        {
            let a = ProjectivePoint::generator().double_multi(2); // 4P
            let b = ProjectivePoint::generator().double(); // 2P
            let c = a + b;

            let mut d = ProjectivePoint::generator();
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
                let a = ProjectivePoint::random(&mut rng);
                let b = AffinePoint::random(&mut rng);

                // If this fails, we are very unlucky!
                assert_eq!(a.add_mixed(&b), a.add_mixed_unchecked(&b));
            }

            let a = ProjectivePoint::random(&mut rng);
            assert_ne!(a.add(&a), a.add_mixed_unchecked(&a.to_affine()));
        }
    }

    #[test]
    #[allow(clippy::eq_op)]
    fn test_projective_negation_and_subtraction() {
        let a = ProjectivePoint::generator().double();
        assert_eq!(a + (-a), ProjectivePoint::identity());
        assert_eq!(a + (-a), a - a);
    }

    #[test]
    fn test_affine_negation_and_subtraction() {
        let a = AffinePoint::generator();
        assert_eq!(ProjectivePoint::from(a) + (-a), ProjectivePoint::identity());
        assert_eq!(
            ProjectivePoint::from(a) + (-a),
            ProjectivePoint::from(a) - a
        );
    }

    #[test]
    fn test_projective_scalar_multiplication() {
        let mut rng = OsRng;
        let g = ProjectivePoint::generator();

        for _ in 0..100 {
            let a = Scalar::random(&mut rng);
            let b = Scalar::random(&mut rng);

            let c = a * b;

            assert_eq!((g * a) * b, g * c);
            assert_eq!(g * c, g.multiply_vartime(&c.to_bytes()));

            assert_eq!(g * c, &BASEPOINT_TABLE * c);
            assert_eq!(g * c, BASEPOINT_TABLE.multiply(&c.to_bytes()));
            assert_eq!(g * c, BASEPOINT_TABLE.multiply_vartime(&c.to_bytes()));

            let g = g * a;
            let basepoint_table = BasePointTable::create(&g);
            assert_eq!(g * b, &basepoint_table * b);
            assert_eq!(g * b, basepoint_table.multiply(&b.to_bytes()));
            assert_eq!(g * b, basepoint_table.multiply_vartime(&b.to_bytes()));
        }
    }

    #[test]
    fn test_projective_double_scalar_multiplication() {
        let mut rng = OsRng;
        let g = ProjectivePoint::generator() * Scalar::random(&mut rng);
        let h = ProjectivePoint::generator() * Scalar::random(&mut rng);

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
        }
    }

    #[test]
    fn test_projective_double_scalar_multiplication_with_basepoint() {
        let mut rng = OsRng;
        let g = ProjectivePoint::generator() * Scalar::random(&mut rng);
        let h = ProjectivePoint::generator();

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
    fn test_affine_scalar_multiplication() {
        let mut rng = OsRng;
        let g = AffinePoint::generator();

        for _ in 0..100 {
            let a = Scalar::random(&mut rng);
            let b = Scalar::random(&mut rng);

            let c = a * b;

            assert_eq!((g * a) * b, g * c);
            assert_eq!(g.multiply_vartime(&c.to_bytes()), g.multiply(&c.to_bytes()));
        }
    }

    #[test]
    fn test_affine_double_scalar_multiplication() {
        let mut rng = OsRng;
        let g = AffinePoint::from(ProjectivePoint::generator() * Scalar::random(&mut rng));
        let h = AffinePoint::from(ProjectivePoint::generator() * Scalar::random(&mut rng));

        for _ in 0..100 {
            let a = Scalar::random(&mut rng);
            let b = Scalar::random(&mut rng);

            assert_eq!(
                g.multiply_double(&h, &a.to_bytes(), &b.to_bytes()),
                AffinePoint::from(g * a + h * b)
            );
            assert_eq!(
                g.multiply_double_vartime(&h, &a.to_bytes(), &b.to_bytes()),
                AffinePoint::from(g * a + h * b)
            );
        }
    }

    #[test]
    fn test_clear_cofactor() {
        // the generator (and the identity) are always on the curve
        let generator = ProjectivePoint::generator();
        assert!(bool::from(generator.clear_cofactor().is_on_curve()));
        let id = ProjectivePoint::identity();
        assert!(bool::from(id.clear_cofactor().is_on_curve()));

        let point = ProjectivePoint {
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
        assert!(bool::from(a.is_on_curve()));
        assert!(!bool::from(a.is_torsion_free()));
        assert!(bool::from(AffinePoint::identity().is_torsion_free()));
        assert!(bool::from(AffinePoint::generator().is_torsion_free()));

        let a = ProjectivePoint::from(&a);
        assert!(bool::from(a.is_on_curve()));
        assert!(!bool::from(a.is_torsion_free()));
        assert!(bool::from(ProjectivePoint::identity().is_torsion_free()));
        assert!(bool::from(ProjectivePoint::generator().is_torsion_free()));
    }

    #[test]
    fn test_batch_normalize() {
        let a = ProjectivePoint::generator().double();
        let b = a.double();
        let c = b.double();

        for a_identity in (0..1).map(|n| n == 1) {
            for b_identity in (0..1).map(|n| n == 1) {
                for c_identity in (0..1).map(|n| n == 1) {
                    let mut v = [a, b, c];
                    if a_identity {
                        v[0] = ProjectivePoint::identity()
                    }
                    if b_identity {
                        v[1] = ProjectivePoint::identity()
                    }
                    if c_identity {
                        v[2] = ProjectivePoint::identity()
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

                    ProjectivePoint::batch_normalize(&v[..], &mut t[..]);

                    assert_eq!(&t[..], &expected[..]);
                }
            }
        }
    }

    #[test]
    fn test_zeroize() {
        use zeroize::Zeroize;

        let mut a = AffinePoint::generator();
        a.zeroize();
        assert_eq!(a, AffinePoint::identity());

        let mut a = ProjectivePoint::generator();
        a.zeroize();
        assert_eq!(a, ProjectivePoint::identity());
    }

    // POINT COMPRESSION
    // ================================================================================================

    #[test]
    fn test_point_compressed() {
        let mut rng = OsRng;
        // Random points
        for _ in 0..100 {
            let point = AffinePoint::random(&mut rng);
            let bytes = point.to_compressed();
            let point_decompressed = AffinePoint::from_compressed(&bytes).unwrap();
            assert_eq!(point, point_decompressed);

            let point = ProjectivePoint::random(&mut rng);
            let bytes = point.to_compressed();
            let point_decompressed = ProjectivePoint::from_compressed(&bytes).unwrap();
            assert_eq!(point, point_decompressed);
        }

        // Identity point
        {
            let bytes = AffinePoint::identity().to_compressed();
            let point_decompressed = AffinePoint::from_compressed(&bytes).unwrap();
            assert!(bool::from(point_decompressed.is_identity()));

            let bytes = ProjectivePoint::identity().to_compressed();
            let point_decompressed = ProjectivePoint::from_compressed(&bytes).unwrap();
            assert!(bool::from(point_decompressed.is_identity()));

            assert_eq!(
                ProjectivePoint::identity().to_compressed().to_bytes(),
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
            let bytes = point.to_compressed();
            let point_decompressed = AffinePoint::from_compressed(&bytes);
            assert!(bool::from(point_decompressed.is_none()));

            let point = ProjectivePoint::from(&point);
            let bytes = point.to_compressed();
            let point_decompressed = ProjectivePoint::from_compressed(&bytes);
            assert!(bool::from(point_decompressed.is_none()));
        }
        {
            let bytes = CompressedPoint([255; 49]);
            let point_decompressed = AffinePoint::from_compressed_unchecked(&bytes);
            assert!(bool::from(point_decompressed.is_none()));
        }
    }

    #[test]
    fn test_point_uncompressed() {
        let mut rng = OsRng;

        // Random points
        for _ in 0..100 {
            let point = AffinePoint::random(&mut rng);
            let bytes = point.to_uncompressed();
            let point_decompressed = AffinePoint::from_uncompressed(&bytes).unwrap();
            assert_eq!(point, point_decompressed);

            let point = ProjectivePoint::random(&mut rng);
            let bytes = point.to_uncompressed();
            let point_decompressed = ProjectivePoint::from_uncompressed(&bytes).unwrap();
            assert_eq!(point, point_decompressed);
        }

        // Identity point
        {
            let bytes = AffinePoint::identity().to_uncompressed();
            let point_decompressed = AffinePoint::from_uncompressed(&bytes).unwrap();
            assert!(bool::from(point_decompressed.is_identity()));

            let bytes = ProjectivePoint::identity().to_uncompressed();
            let point_decompressed = ProjectivePoint::from_uncompressed(&bytes).unwrap();
            assert!(bool::from(point_decompressed.is_identity()));

            assert_eq!(
                ProjectivePoint::identity().to_uncompressed().to_bytes(),
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
            let bytes = point.to_uncompressed();
            let point_decompressed = AffinePoint::from_uncompressed(&bytes);
            assert!(bool::from(point_decompressed.is_none()));

            let point = ProjectivePoint::from(&point);
            let bytes = point.to_uncompressed();
            let point_decompressed = ProjectivePoint::from_uncompressed(&bytes);
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
            let point_decompressed = AffinePoint::from_uncompressed_unchecked(&bytes);
            assert!(bool::from(point_decompressed.is_none()));

            let point_decompressed = ProjectivePoint::from_uncompressed_unchecked(&bytes);
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
            let point_decompressed = AffinePoint::from_uncompressed_unchecked(&bytes);
            assert!(bool::from(point_decompressed.is_none()));

            let point_decompressed = ProjectivePoint::from_uncompressed_unchecked(&bytes);
            assert!(bool::from(point_decompressed.is_none()));
        }
    }

    // SERDE SERIALIZATIOIN
    // ================================================================================================

    #[test]
    #[cfg(feature = "serialize")]
    fn test_serde_affine() {
        let mut rng = OsRng;
        let point = AffinePoint::random(&mut rng);
        let encoded = bincode::serialize(&point).unwrap();
        let parsed: AffinePoint = bincode::deserialize(&encoded).unwrap();
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
        assert!(bincode::deserialize::<AffinePoint>(&encoded[0..47]).is_err());

        let encoded = [255; 49];
        assert!(bincode::deserialize::<AffinePoint>(&encoded).is_err());
    }

    #[test]
    #[cfg(feature = "serialize")]
    fn test_serde_projective() {
        let mut rng = OsRng;
        let point = ProjectivePoint::random(&mut rng);
        let encoded = bincode::serialize(&point).unwrap();
        let parsed: ProjectivePoint = bincode::deserialize(&encoded).unwrap();
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
        assert!(bincode::deserialize::<ProjectivePoint>(&encoded[0..47]).is_err());

        let encoded = [255; 49];
        assert!(bincode::deserialize::<ProjectivePoint>(&encoded).is_err());
    }
}
