// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This module provides an implementation of the STARK-friendly
//! cheetah curve over the sextic extension of the prime field Fp
//! of characteristic p = 2^62 + 2^56 + 2^55 + 1.

use core::{
    borrow::Borrow,
    fmt,
    iter::Sum,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::{Fp, Fp2, Fp6, Scalar};

use crate::LookupTable;

use group::ff::Field;
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
/// B = (1200866201009650596 * u + 1935817186716799185) * v^2
///         + (3999205700308519553 * u + 3518137720867787056) * v
///         + 2508413708960025374 * u + 1526905369741321712
pub const B: Fp6 = Fp6 {
    c0: Fp2 {
        c0: Fp(0x15ae48ea5d89311b),
        c1: Fp(0x17f1a248aadda1c0),
    },
    c1: Fp2 {
        c0: Fp(0x292d3259d8078df0),
        c1: Fp(0x32955e35847805dd),
    },
    c2: Fp2 {
        c0: Fp(0x40805fda8f3d4cf1),
        c1: Fp(0x30dd795db5ec864c),
    },
};

const B3: Fp6 = (&B).mul(&Fp6::new([3, 0, 0, 0, 0, 0]));

// cofactor = 639750783962362599710832166084722337
//          = 0x7b36261eef444af4e22bbe55a2bea1
const COFACTOR_BYTES: [u8; 15] = [
    161, 190, 162, 85, 190, 43, 226, 244, 74, 68, 239, 30, 38, 54, 123,
];

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

        if zinv == Fp6::zero() {
            AffinePoint::identity()
        } else {
            tmp
        }
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
    pub fn get_x(&self) -> Fp6 {
        self.x
    }

    /// Returns the y coordinate of this AffinePoint
    pub fn get_y(&self) -> Fp6 {
        self.y
    }

    /// Computes a random `AffinePoint` element
    pub fn random(mut rng: impl RngCore) -> Self {
        ProjectivePoint::random(&mut rng).into()
    }

    /// Returns a fixed generator of the curve in affine coordinates
    pub fn generator() -> AffinePoint {
        AffinePoint {
            x: Fp6 {
                c0: Fp2 {
                    c0: Fp(0xf6798582c92ece1),
                    c1: Fp(0x2b7c30a4c7d886c0),
                },
                c1: Fp2 {
                    c0: Fp(0x1269cdae98dc2fd0),
                    c1: Fp(0x11b78ef6c71c6132),
                },
                c2: Fp2 {
                    c0: Fp(0x3ac2244dfc47537),
                    c1: Fp(0x36dfeea4b9051daf),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp(0x334807e450d55e2f),
                    c1: Fp(0x200a54d42b84bd17),
                },
                c1: Fp2 {
                    c0: Fp(0x271af7bb20ab32e1),
                    c1: Fp(0x3df7b90927efc7ec),
                },
                c2: Fp2 {
                    c0: Fp(0xab8bbf4a53af6a0),
                    c1: Fp(0xe13dca26b2ac6ab),
                },
            },
            infinity: Choice::from(0u8),
        }
    }

    /// Outputs a compress byte representation of this `AffinePoint` element
    pub fn to_compressed(&self) -> [u8; 48] {
        // Strictly speaking, self.x is zero already when self.infinity is true, but
        // to guard against implementation mistakes we do not assume this.
        let mut res = Fp6::conditional_select(&self.x, &Fp6::zero(), self.infinity).to_bytes();

        // This point is in compressed form, so we set the most significant bit of the first Fp element.
        res[7] |= 1u8 << 7;

        // Is this point at infinity? If so, set the most significant bit of the second Fp element.
        res[15] |= u8::conditional_select(&0u8, &(1u8 << 7), self.infinity);

        // Is the y-coordinate the lexicographically largest of the two associated with the
        // x-coordinate? If so, the most significant bit of the third Fp element so long as this is not
        // the point at infinity.
        res[23] |= u8::conditional_select(
            &0u8,
            &(1u8 << 7),
            (!self.infinity) & self.y.lexicographically_largest(),
        );

        res
    }

    /// Outputs an uncompressed byte representation of this `AffinePoint` element
    /// It is twice larger than when calling `AffinePoint::to_compressed()`
    pub fn to_uncompressed(&self) -> [u8; 96] {
        let mut res = [0; 96];

        res[0..48].copy_from_slice(
            &Fp6::conditional_select(&self.x, &Fp6::zero(), self.infinity).to_bytes()[..],
        );
        res[48..96].copy_from_slice(
            &Fp6::conditional_select(&self.y, &Fp6::zero(), self.infinity).to_bytes()[..],
        );

        // Is this point at infinity? If so, set the most significant bit of the second Fp element.
        res[15] |= u8::conditional_select(&0u8, &(1u8 << 7), self.infinity);

        res
    }

    /// Attempts to deserialize an uncompressed element.
    pub fn from_uncompressed(bytes: &[u8; 96]) -> CtOption<Self> {
        Self::from_uncompressed_unchecked(bytes)
            .and_then(|p| CtOption::new(p, p.is_on_curve() & p.is_torsion_free()))
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is on the curve and not checking if it is in the correct subgroup.
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_uncompressed()` instead.
    pub fn from_uncompressed_unchecked(bytes: &[u8; 96]) -> CtOption<Self> {
        // Obtain the three flags
        let compression_flag_set = Choice::from((bytes[7] >> 7) & 1);
        let infinity_flag_set = Choice::from((bytes[15] >> 7) & 1);
        let sort_flag_set = Choice::from((bytes[23] >> 7) & 1);

        // Attempt to obtain the x-coordinate
        let x = {
            let mut tmp = [0; 48];
            tmp.copy_from_slice(&bytes[0..48]);

            // Mask away the flag bits
            tmp[7] &= 0b0111_1111;
            tmp[15] &= 0b0111_1111;
            tmp[23] &= 0b0111_1111;

            Fp6::from_bytes(&tmp)
        };

        // Attempt to obtain the y-coordinate
        let y = {
            let mut tmp = [0; 48];
            tmp.copy_from_slice(&bytes[48..96]);

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
                    ((!infinity_flag_set) | (infinity_flag_set & x.is_zero() & y.is_zero())) &
                    // The compression flag should not have been set, as this is an uncompressed element
                    (!compression_flag_set) &
                    // The sort flag should not have been set, as this is an uncompressed element
                    (!sort_flag_set),
                )
            })
        })
    }

    /// Attempts to deserialize a compressed element.
    pub fn from_compressed(bytes: &[u8; 48]) -> CtOption<Self> {
        // We already know the point is on the curve because this is established
        // by the y-coordinate recovery procedure in from_compressed_unchecked().

        Self::from_compressed_unchecked(bytes).and_then(|p| CtOption::new(p, p.is_torsion_free()))
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is in the correct subgroup.
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_compressed()` instead.
    pub fn from_compressed_unchecked(bytes: &[u8; 48]) -> CtOption<Self> {
        // Obtain the three flags
        let compression_flag_set = Choice::from((bytes[7] >> 7) & 1);
        let infinity_flag_set = Choice::from((bytes[15] >> 7) & 1);
        let sort_flag_set = Choice::from((bytes[23] >> 7) & 1);

        // Attempt to obtain the x-coordinate
        let x = {
            let mut tmp = [0; 48];
            tmp.copy_from_slice(&bytes[0..48]);

            // Mask away the flag bits
            tmp[7] &= 0b0111_1111;
            tmp[15] &= 0b0111_1111;
            tmp[23] &= 0b0111_1111;

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
            compression_flag_set & // Compression flag should be set
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
                        (!infinity_flag_set) & // Infinity flag should not be set
                        compression_flag_set, // Compression flag should be set
                    )
                })
            })
        })
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
    pub fn is_identity(&self) -> Choice {
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

    /// Multiplies by the curve cofactor
    pub fn clear_cofactor(&self) -> AffinePoint {
        let point: ProjectivePoint = self.into();

        point.clear_cofactor().into()
    }

    /// Returns true if this point is free of an $h$-torsion component.
    /// This should always return true unless an "unchecked" API was used.
    pub fn is_torsion_free(&self) -> Choice {
        const FQ_MODULUS_BYTES: [u8; 32] = [
            165, 106, 149, 48, 225, 174, 19, 30, 39, 43, 36, 53, 240, 100, 242, 242, 201, 236, 143,
            161, 185, 235, 182, 182, 124, 247, 149, 39, 117, 127, 51, 38,
        ];

        // Clear the r-torsion from the point and check if it is the identity
        self.multiply(&FQ_MODULUS_BYTES).is_identity()
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
        ProjectivePoint {
            x: self.x,
            y: -self.y,
            z: self.z,
        }
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
    a.mul(&B3)
}

impl ProjectivePoint {
    /// Returns the identity of the group: the point at infinity.
    pub fn identity() -> ProjectivePoint {
        ProjectivePoint {
            x: Fp6::zero(),
            y: Fp6::one(),
            z: Fp6::zero(),
        }
    }

    /// Returns the x coordinate of this ProjectivePoint
    pub fn get_x(&self) -> Fp6 {
        self.x
    }

    /// Returns the y coordinate of this ProjectivePoint
    pub fn get_y(&self) -> Fp6 {
        self.y
    }

    /// Returns the z coordinate of this ProjectivePoint
    pub fn get_z(&self) -> Fp6 {
        self.z
    }

    /// Returns a fixed generator of the curve in projective coordinates
    pub fn generator() -> ProjectivePoint {
        ProjectivePoint::from(AffinePoint::generator())
    }

    /// Outputs a compress byte representation of this `ProjectivePoint` element
    pub fn to_compressed(&self) -> [u8; 48] {
        AffinePoint::from(self).to_compressed()
    }

    /// Outputs an uncompressed byte representation of this `ProjectivePoint` element
    /// It is twice larger than when calling `ProjectivePoint::to_uncompress()`
    pub fn to_uncompressed(&self) -> [u8; 96] {
        AffinePoint::from(self).to_uncompressed()
    }

    /// Attempts to deserialize an uncompressed element.
    pub fn from_uncompressed(bytes: &[u8; 96]) -> CtOption<Self> {
        AffinePoint::from_uncompressed(bytes)
            .and_then(|p| CtOption::new(ProjectivePoint::from(p), 1.into()))
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is on the curve and not checking if it is in the correct subgroup.
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_uncompressed()` instead.
    pub fn from_uncompressed_unchecked(bytes: &[u8; 96]) -> CtOption<Self> {
        AffinePoint::from_uncompressed_unchecked(bytes)
            .and_then(|p| CtOption::new(ProjectivePoint::from(p), 1.into()))
    }

    /// Attempts to deserialize a compressed element.
    pub fn from_compressed(bytes: &[u8; 48]) -> CtOption<Self> {
        AffinePoint::from_compressed(bytes).map(ProjectivePoint::from)
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
        let a = a + x2;
        let a = a + z2;

        let b = self.y * self.z;
        let t3 = self.y * b;
        let c = t3 * self.x;

        let d = c.double();
        let c4 = d.double();
        let d = c4.double();

        let t = a.square();
        let d = t - d;

        let x3 = b * d;
        let x3 = x3.double();

        let y3 = c4 - d;
        let y3 = a * y3;
        let t3 = t3.square();

        let t3 = t3.double();
        let t3 = t3.double();
        let t3 = t3.double();

        let y3 = y3 - t3;
        let z3 = b.square();
        let z3 = z3 * b;

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

    /// Adds this point to another point.
    pub fn add(&self, rhs: &ProjectivePoint) -> ProjectivePoint {
        // Algorithm 1, https://eprint.iacr.org/2015/1060.pdf

        let t0 = self.x * rhs.x;
        let t1 = self.y * rhs.y;
        let t2 = self.z * rhs.z;

        let t3 = self.x + self.y;
        let t4 = rhs.x + rhs.y;
        let t3 = t3 * t4;

        let t4 = t0 + t1;
        let t3 = t3 - t4;
        let t4 = self.x + self.z;

        let t5 = rhs.x + rhs.z;
        let t4 = t4 * t5;
        let t5 = t0 + t2;

        let t4 = t4 - t5;
        let t5 = self.y + self.z;
        let x3 = rhs.y + rhs.z;

        let t5 = t5 * x3;
        let x3 = t1 + t2;
        let t5 = t5 - x3;

        let x3 = mul_by_3b(&t2);
        let z3 = x3 + t4;

        let x3 = t1 - z3;
        let z3 = t1 + z3;
        let y3 = x3 * z3;

        let t1 = t0.double();
        let t1 = t1 + t0;

        let t4 = mul_by_3b(&t4);
        let t1 = t1 + t2;
        let t2 = t0 - t2;

        let t4 = t4 + t2;
        let t0 = t1 * t4;

        let y3 = y3 + t0;
        let t0 = t5 * t4;
        let x3 = t3 * x3;

        let x3 = x3 - t0;
        let t0 = t3 * t1;
        let z3 = t5 * z3;

        let z3 = z3 + t0;

        ProjectivePoint {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    /// Adds this point to another point in the affine model.
    pub fn add_mixed(&self, rhs: &AffinePoint) -> ProjectivePoint {
        // Algorithm 2, https://eprint.iacr.org/2015/1060.pdf

        let t0 = self.x * rhs.x;
        let t1 = self.y * rhs.y;
        let t3 = rhs.x + rhs.y;

        let t4 = self.x + self.y;
        let t3 = t3 * t4;
        let t4 = t0 + t1;

        let t3 = t3 - t4;
        let t4 = rhs.x * self.z;
        let t4 = t4 + self.x;

        let t5 = rhs.y * self.z;
        let t5 = t5 + self.y;

        let x3 = mul_by_3b(&self.z);
        let z3 = x3 + t4;
        let x3 = t1 - z3;

        let z3 = t1 + z3;
        let y3 = x3 * z3;
        let t1 = t0.double();

        let t1 = t1 + t0;
        let t4 = mul_by_3b(&t4);

        let t1 = t1 + self.z;
        let t2 = t0 - self.z;

        let t4 = t4 + t2;
        let t0 = t1 * t4;
        let y3 = y3 + t0;

        let t0 = t5 * t4;
        let x3 = t3 * x3;
        let x3 = x3 - t0;

        let t0 = t3 * t1;
        let z3 = t5 * z3;
        let z3 = z3 + t0;

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
        let mut acc = ProjectivePoint::identity();
        let table = LookupTable::<16>::from(self);

        for digit in by
            .iter()
            .rev()
            .flat_map(|byte| (0..2).rev().map(move |i| (byte >> (i * 4)) & 15u8))
        {
            acc = acc.double_multi(4);
            acc = ProjectivePoint::conditional_select(
                &acc,
                &(acc + table.get_point(digit as i8)),
                ((digit != 0) as u8).into(),
            );
        }

        acc
    }

    /// Performs a projective scalar multiplication from `by`
    /// given as byte representation of a `Scalar` element.
    ///
    /// **This operation is variable time with respect
    /// to the scalar.** If the scalar is fixed,
    /// this operation is effectively constant time.
    pub fn multiply_vartime(&self, by: &[u8; 32]) -> ProjectivePoint {
        let mut acc = ProjectivePoint::identity();
        let table = LookupTable::<16>::from(self);

        for digit in by
            .iter()
            .rev()
            .flat_map(|byte| (0..2).rev().map(move |i| (byte >> (i * 4)) & 15u8))
        {
            acc = acc.double_multi(4);
            if digit != 0 {
                acc += table.get_point_vartime(digit as i8);
            }
        }

        acc
    }

    /// Performs the projective sum [`by_lhs` * `self` + `by_rhs` * `rhs`] with
    /// `by_lhs` and `by_rhs` given as byte representations of `Scalar` elements.
    #[inline]
    pub fn multiply_double(
        &self,
        rhs: &ProjectivePoint,
        by_lhs: &[u8; 32],
        by_rhs: &[u8; 32],
    ) -> ProjectivePoint {
        let mut acc = ProjectivePoint::identity();
        let table_lhs = LookupTable::<16>::from(self);
        let table_rhs = LookupTable::<16>::from(rhs);

        for (digit_lhs, digit_rhs) in by_lhs
            .iter()
            .rev()
            .flat_map(|byte| (0..2).rev().map(move |i| (byte >> (i * 4)) & 15u8))
            .zip(
                by_rhs
                    .iter()
                    .rev()
                    .flat_map(|byte| (0..2).rev().map(move |i| (byte >> (i * 4)) & 15u8)),
            )
        {
            acc = acc.double_multi(4);
            acc = ProjectivePoint::conditional_select(
                &acc,
                &(acc + table_lhs.get_point(digit_lhs as i8)),
                ((digit_lhs != 0) as u8).into(),
            );
            acc = ProjectivePoint::conditional_select(
                &acc,
                &(acc + table_rhs.get_point(digit_rhs as i8)),
                ((digit_rhs != 0) as u8).into(),
            );
        }

        acc
    }

    /// Performs the projective sum [`by_lhs` * `self` + `by_rhs` * `rhs`] with
    /// `by_lhs` and `by_rhs` given as byte representations of `Scalar` elements.
    ///
    /// **This operation is variable time with respect
    /// to the scalars.** If the scalars are fixed,
    /// this operation is effectively constant time.
    #[inline]
    pub fn multiply_double_vartime(
        &self,
        rhs: &ProjectivePoint,
        by_lhs: &[u8; 32],
        by_rhs: &[u8; 32],
    ) -> ProjectivePoint {
        let mut acc = ProjectivePoint::identity();
        let table_lhs = LookupTable::<16>::from(self);
        let table_rhs = LookupTable::<16>::from(rhs);

        for (digit_lhs, digit_rhs) in by_lhs
            .iter()
            .rev()
            .flat_map(|byte| (0..2).rev().map(move |i| (byte >> (i * 4)) & 15u8))
            .zip(
                by_rhs
                    .iter()
                    .rev()
                    .flat_map(|byte| (0..2).rev().map(move |i| (byte >> (i * 4)) & 15u8)),
            )
        {
            acc = acc.double_multi(4);
            if digit_lhs != 0 {
                acc += table_lhs.get_point_vartime(digit_lhs as i8);
            }
            if digit_rhs != 0 {
                acc += table_rhs.get_point_vartime(digit_rhs as i8);
            }
        }

        acc
    }

    /// Multiplies by the curve cofactor
    pub fn clear_cofactor(&self) -> ProjectivePoint {
        let mut acc = ProjectivePoint::identity();

        // This is a simple double-and-add implementation of point
        // multiplication, moving from most significant to least
        // significant bit of the scalar.
        for bit in COFACTOR_BYTES
            .iter()
            .rev()
            .flat_map(|byte| (0..8).rev().map(move |i| ((byte >> i) & 1u8) != 0))
        {
            acc = acc.double();
            acc = if bit { acc + self } else { acc };
        }

        acc
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

// GROUP TRAITS IMPLEMENTATION
// ================================================================================================

impl Group for ProjectivePoint {
    type Scalar = Scalar;

    fn random(mut rng: impl RngCore) -> Self {
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
                let p = p.unwrap().clear_cofactor();

                if bool::from(!p.is_identity()) {
                    return p.into();
                }
            }
        }
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
impl Serialize for AffinePoint {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        use serde::ser::SerializeTuple;
        let mut tup = serializer.serialize_tuple(48)?;
        for byte in self.to_compressed().iter() {
            tup.serialize_element(byte)?;
        }
        tup.end()
    }
}

#[cfg(feature = "serialize")]
impl Serialize for ProjectivePoint {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        use serde::ser::SerializeTuple;
        let mut tup = serializer.serialize_tuple(48)?;
        for byte in self.to_compressed().iter() {
            tup.serialize_element(byte)?;
        }
        tup.end()
    }
}

#[cfg(feature = "serialize")]
impl<'de> Deserialize<'de> for AffinePoint {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct AffinePointVisitor;

        impl<'de> Visitor<'de> for AffinePointVisitor {
            type Value = AffinePoint;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("a valid point in Ristretto format")
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<AffinePoint, A::Error>
            where
                A: serde::de::SeqAccess<'de>,
            {
                let mut bytes = [0u8; 48];
                for (i, byte) in bytes.iter_mut().enumerate() {
                    *byte = seq
                        .next_element()?
                        .ok_or_else(|| serde::de::Error::invalid_length(i, &"expected 48 bytes"))?;
                }

                let p = AffinePoint::from_compressed(&bytes);
                if bool::from(p.is_some()) {
                    Ok(p.unwrap())
                } else {
                    Err(serde::de::Error::custom("decompression failed"))
                }
            }
        }

        deserializer.deserialize_tuple(48, AffinePointVisitor)
    }
}

#[cfg(feature = "serialize")]
impl<'de> Deserialize<'de> for ProjectivePoint {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct ProjectivePointVisitor;

        impl<'de> Visitor<'de> for ProjectivePointVisitor {
            type Value = ProjectivePoint;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("48 bytes of data")
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<ProjectivePoint, A::Error>
            where
                A: serde::de::SeqAccess<'de>,
            {
                let mut bytes = [0u8; 48];
                for (i, byte) in bytes.iter_mut().enumerate() {
                    *byte = seq
                        .next_element()?
                        .ok_or_else(|| serde::de::Error::invalid_length(i, &"expected 48 bytes"))?;
                }

                let p = ProjectivePoint::from_compressed(&bytes);
                if bool::from(p.is_some()) {
                    Ok(p.unwrap())
                } else {
                    Err(serde::de::Error::custom("decompression failed"))
                }
            }
        }

        deserializer.deserialize_tuple(48, ProjectivePointVisitor)
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
            c0: Fp2 {
                c0: Fp::new(1780366954230337028),
                c1: Fp::new(1657393229643603266),
            },
            c1: Fp2 {
                c0: Fp::new(266397157796383323),
                c1: Fp::new(1882752790847015855),
            },
            c2: Fp2 {
                c0: Fp::new(3490694662347740613),
                c1: Fp::new(1934749607384548809),
            },
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
            c0: Fp2 {
                c0: Fp::new(1780366954230337028),
                c1: Fp::new(1657393229643603266),
            },
            c1: Fp2 {
                c0: Fp::new(266397157796383323),
                c1: Fp::new(1882752790847015855),
            },
            c2: Fp2 {
                c0: Fp::new(3490694662347740613),
                c1: Fp::new(1934749607384548809),
            },
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
            c0: Fp2 {
                c0: Fp::new(1780366954230337028),
                c1: Fp::new(1657393229643603266),
            },
            c1: Fp2 {
                c0: Fp::new(266397157796383323),
                c1: Fp::new(1882752790847015855),
            },
            c2: Fp2 {
                c0: Fp::new(3490694662347740613),
                c1: Fp::new(1934749607384548809),
            },
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
                        c0: Fp2 {
                            c0: Fp(0x1ba2d52806f212a),
                            c1: Fp(0x5e9353a4e8225c8),
                        },
                        c1: Fp2 {
                            c0: Fp(0x13e92423fef3bc2d),
                            c1: Fp(0x241081e7ae1db310),
                        },
                        c2: Fp2 {
                            c0: Fp(0x29f0073c3351026b),
                            c1: Fp(0x11233fe9eb7285c0),
                        },
                    },
                    y: Fp6 {
                        c0: Fp2 {
                            c0: Fp(0x3a19dfba18e15ed5),
                            c1: Fp(0x3691eb6949fca20b),
                        },
                        c1: Fp2 {
                            c0: Fp(0x3ea42cb9ad7430ab),
                            c1: Fp(0x1b840f91119a2eb3),
                        },
                        c2: Fp2 {
                            c0: Fp(0x1b94f8ccdafc47ba),
                            c1: Fp(0x19e92e12c3a9cfa),
                        },
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
                    c0: Fp2 {
                        c0: Fp::new(1780366954230337028),
                        c1: Fp::new(1657393229643603266),
                    },
                    c1: Fp2 {
                        c0: Fp::new(266397157796383323),
                        c1: Fp::new(1882752790847015855),
                    },
                    c2: Fp2 {
                        c0: Fp::new(3490694662347740613),
                        c1: Fp::new(1934749607384548809),
                    },
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
            let a = ProjectivePoint::generator().double().double(); // 4P
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
                    c0: Fp2 {
                        c0: Fp::new(1780366954230337028),
                        c1: Fp::new(1657393229643603266),
                    },
                    c1: Fp2 {
                        c0: Fp::new(266397157796383323),
                        c1: Fp::new(1882752790847015855),
                    },
                    c2: Fp2 {
                        c0: Fp::new(3490694662347740613),
                        c1: Fp::new(1934749607384548809),
                    },
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
                    c0: Fp2 {
                        c0: Fp::new(1780366954230337028),
                        c1: Fp::new(1657393229643603266),
                    },
                    c1: Fp2 {
                        c0: Fp::new(266397157796383323),
                        c1: Fp::new(1882752790847015855),
                    },
                    c2: Fp2 {
                        c0: Fp::new(3490694662347740613),
                        c1: Fp::new(1934749607384548809),
                    },
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
            let a = ProjectivePoint::generator().double().double(); // 4P
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
                c0: Fp2 {
                    c0: Fp::new(4088869953255777905),
                    c1: Fp::new(858439445101292708),
                },
                c1: Fp2 {
                    c0: Fp::new(3310975763830793951),
                    c1: Fp::new(2790117791195240420),
                },
                c2: Fp2 {
                    c0: Fp::new(3006738578199540275),
                    c1: Fp::new(2480895684612942729),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp::new(1351667951327929939),
                    c1: Fp::new(2705247942932011006),
                },
                c1: Fp2 {
                    c0: Fp::new(1639335439078005502),
                    c1: Fp::new(2264581412139512278),
                },
                c2: Fp2 {
                    c0: Fp::new(825746373732604926),
                    c1: Fp::new(3738666180738357964),
                },
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
                c0: Fp2 {
                    c0: Fp::new(4088869953255777905),
                    c1: Fp::new(858439445101292708),
                },
                c1: Fp2 {
                    c0: Fp::new(3310975763830793951),
                    c1: Fp::new(2790117791195240420),
                },
                c2: Fp2 {
                    c0: Fp::new(3006738578199540275),
                    c1: Fp::new(2480895684612942729),
                },
            },
            y: Fp6 {
                c0: Fp2 {
                    c0: Fp::new(1351667951327929939),
                    c1: Fp::new(2705247942932011006),
                },
                c1: Fp2 {
                    c0: Fp::new(1639335439078005502),
                    c1: Fp::new(2264581412139512278),
                },
                c2: Fp2 {
                    c0: Fp::new(825746373732604926),
                    c1: Fp::new(3738666180738357964),
                },
            },
            infinity: Choice::from(0u8),
        };
        assert!(bool::from(a.is_on_curve()));
        assert!(!bool::from(a.is_torsion_free()));
        assert!(bool::from(AffinePoint::identity().is_torsion_free()));
        assert!(bool::from(AffinePoint::generator().is_torsion_free()));
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
            let bytes = [
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            ];
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
            let bytes = [
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 255,
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            ];
            let point_decompressed = AffinePoint::from_uncompressed_unchecked(&bytes);
            assert!(bool::from(point_decompressed.is_none()));

            let point_decompressed = ProjectivePoint::from_uncompressed_unchecked(&bytes);
            assert!(bool::from(point_decompressed.is_none()));
        }
        {
            let bytes = [
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255,
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
            ];
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

        // Check that the encoding is 48 bytes exactly
        assert_eq!(encoded.len(), 48);

        // Check that the encoding itself matches the usual one
        assert_eq!(point, bincode::deserialize(&point.to_compressed()).unwrap());

        // Check that invalid encodings fail
        let point = AffinePoint::random(&mut rng);
        let mut encoded = bincode::serialize(&point).unwrap();
        encoded[47] = 127;
        assert!(bincode::deserialize::<AffinePoint>(&encoded).is_err());

        let encoded = bincode::serialize(&point).unwrap();
        assert!(bincode::deserialize::<AffinePoint>(&encoded[0..47]).is_err());
    }

    #[test]
    #[cfg(feature = "serialize")]
    fn test_serde_projective() {
        let mut rng = OsRng;
        let point = ProjectivePoint::random(&mut rng);
        let encoded = bincode::serialize(&point).unwrap();
        let parsed: ProjectivePoint = bincode::deserialize(&encoded).unwrap();
        assert_eq!(parsed, point);

        // Check that the encoding is 48 bytes exactly
        assert_eq!(encoded.len(), 48);

        // Check that the encoding itself matches the usual one
        assert_eq!(point, bincode::deserialize(&point.to_compressed()).unwrap());

        // Check that invalid encodings fail
        let point = ProjectivePoint::random(&mut rng);
        let mut encoded = bincode::serialize(&point).unwrap();
        encoded[47] = 127;
        assert!(bincode::deserialize::<ProjectivePoint>(&encoded).is_err());

        let encoded = bincode::serialize(&point).unwrap();
        assert!(bincode::deserialize::<ProjectivePoint>(&encoded[0..47]).is_err());
    }
}
