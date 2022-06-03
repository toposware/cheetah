// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This module provides an implementation of the cheetah curve in
//! affine coordinates.

use core::{
    fmt,
    hash::{Hash, Hasher},
    ops::{Add, Mul, Neg, Sub},
};

use super::B;

use crate::{CompressedPoint, UncompressedPoint};
use crate::{Fp, Fp6, Scalar};
use crate::{JacobianPoint, ModifiedJacobianPoint, ProjectivePoint};

use rand_core::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "serialize")]
use serde::{self, Deserialize, Deserializer, Serialize, Serializer};

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

impl<'a> From<&'a JacobianPoint> for AffinePoint {
    fn from(p: &'a JacobianPoint) -> AffinePoint {
        let z3 = p.z.square() * p.z;
        let zinv3 = z3.invert().unwrap_or(Fp6::zero());
        let x = p.x * zinv3 * p.z;
        let y = p.y * zinv3;

        let tmp = AffinePoint {
            x,
            y,
            infinity: Choice::from(0u8),
        };

        AffinePoint::conditional_select(&tmp, &AffinePoint::identity(), zinv3.ct_eq(&Fp6::zero()))
    }
}

impl From<JacobianPoint> for AffinePoint {
    fn from(p: JacobianPoint) -> AffinePoint {
        AffinePoint::from(&p)
    }
}

impl<'a> From<&'a ModifiedJacobianPoint> for AffinePoint {
    fn from(p: &'a ModifiedJacobianPoint) -> AffinePoint {
        let zinv4 = p.z4.invert().unwrap_or(Fp6::zero());
        let zinv3 = zinv4 * p.z;
        let x = p.x * zinv3 * p.z;
        let y = p.y * zinv3;

        let tmp = AffinePoint {
            x,
            y,
            infinity: Choice::from(0u8),
        };

        AffinePoint::conditional_select(&tmp, &AffinePoint::identity(), zinv4.ct_eq(&Fp6::zero()))
    }
}

impl From<ModifiedJacobianPoint> for AffinePoint {
    fn from(p: ModifiedJacobianPoint) -> AffinePoint {
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

impl<'a, 'b> Sub<&'b ProjectivePoint> for &'a AffinePoint {
    type Output = ProjectivePoint;

    #[inline]
    fn sub(self, rhs: &'b ProjectivePoint) -> ProjectivePoint {
        self + (-rhs)
    }
}

impl<'a, 'b> Add<&'b JacobianPoint> for &'a AffinePoint {
    type Output = JacobianPoint;

    #[inline]
    fn add(self, rhs: &'b JacobianPoint) -> JacobianPoint {
        rhs.add_mixed(self)
    }
}

impl<'a, 'b> Sub<&'b JacobianPoint> for &'a AffinePoint {
    type Output = JacobianPoint;

    #[inline]
    fn sub(self, rhs: &'b JacobianPoint) -> JacobianPoint {
        self + (-rhs)
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

impl<'a, 'b> Mul<&'b Scalar> for &'a AffinePoint {
    type Output = ProjectivePoint;

    fn mul(self, other: &'b Scalar) -> Self::Output {
        ProjectivePoint::from(self).multiply(&other.to_bytes())
    }
}

impl_binops_multiplicative_mixed!(AffinePoint, Scalar, ProjectivePoint);

// SERDE SERIALIZATION
// ================================================================================================

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

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::OsRng;

    #[test]
    fn test_is_on_curve() {
        assert!(bool::from(AffinePoint::identity().is_on_curve()));
        assert!(bool::from(AffinePoint::generator().is_on_curve()));

        let mut test = AffinePoint::generator();
        test.y = Fp6::zero();

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
    }

    #[test]
    fn test_zeroize() {
        use zeroize::Zeroize;

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
        }

        // Identity point
        {
            let bytes = AffinePoint::identity().to_compressed();
            let point_decompressed = AffinePoint::from_compressed(&bytes).unwrap();
            assert!(bool::from(point_decompressed.is_identity()));
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
        }

        // Identity point
        {
            let bytes = AffinePoint::identity().to_uncompressed();
            let point_decompressed = AffinePoint::from_uncompressed(&bytes).unwrap();
            assert!(bool::from(point_decompressed.is_identity()));
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
}
