// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This module provides an implementation of the STARK-friendly
//! cheetah curve encoding.

use super::B;

use crate::AffinePoint;
use crate::Fp6;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "serialize")]
use core::fmt;
#[cfg(feature = "serialize")]
use serde::de::Visitor;
#[cfg(feature = "serialize")]
use serde::{self, Deserialize, Deserializer, Serialize, Serializer};

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
    pub(crate) fn from_affine(point: &AffinePoint) -> Self {
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
    pub(crate) fn to_affine(self) -> CtOption<AffinePoint> {
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
    pub(crate) fn from_affine(point: &AffinePoint) -> Self {
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
    /// The resulting point is not ensured to be on the curve.
    pub(crate) fn to_affine(self) -> CtOption<AffinePoint> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::OsRng;

    use crate::ProjectivePoint;

    // POINT COMPRESSION
    // ================================================================================================

    #[test]
    fn test_point_compressed() {
        let mut rng = OsRng;
        // Random points
        for _ in 0..100 {
            let point = AffinePoint::random(&mut rng);
            let compressed_point = point.to_compressed();
            let point_decompressed = compressed_point.to_affine().unwrap();
            assert_eq!(point, point_decompressed);
        }

        // Identity point
        {
            let point = AffinePoint::identity();
            let compressed_point = point.to_compressed();
            let point_decompressed = compressed_point.to_affine().unwrap();
            assert_eq!(point, point_decompressed);

            assert_eq!(
                compressed_point.to_bytes(),
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

            let compressed_point = point.to_compressed();
            let point_decompressed = compressed_point.to_affine();
            assert!(bool::from(point_decompressed.is_none()));
        }
    }

    #[test]
    fn test_point_uncompressed() {
        let mut rng = OsRng;

        // Random points
        for _ in 0..100 {
            let point = AffinePoint::random(&mut rng);
            let uncompressed_point = point.to_uncompressed();
            let point_decompressed = uncompressed_point.to_affine().unwrap();
            assert_eq!(point, point_decompressed);
        }

        // Identity point
        {
            let point = AffinePoint::identity();
            let uncompressed_point = point.to_uncompressed();
            let point_decompressed = uncompressed_point.to_affine().unwrap();
            assert_eq!(point, point_decompressed);

            assert_eq!(
                uncompressed_point.to_bytes(),
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
            let uncompressed_point = UncompressedPoint([
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 255,
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            ]);

            let point_decompressed = uncompressed_point.to_affine();
            assert!(bool::from(point_decompressed.is_none()));
        }
        {
            let uncompressed_point = UncompressedPoint([
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255,
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
                255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
            ]);

            let point_decompressed = uncompressed_point.to_affine();
            assert!(bool::from(point_decompressed.is_none()));
        }
    }

    // SERDE SERIALIZATIOIN
    // ================================================================================================

    #[test]
    #[cfg(feature = "serialize")]
    fn test_serde_projective() {
        let mut rng = OsRng;
        let point = ProjectivePoint::random(&mut rng);
        let compressed_point = point.to_compressed();
        let encoded = bincode::serialize(&compressed_point).unwrap();
        let parsed: CompressedPoint = bincode::deserialize(&encoded).unwrap();
        assert_eq!(parsed, compressed_point);

        // Check that the encoding is 49 bytes exactly
        assert_eq!(encoded.len(), 49);

        // Check that the encoding itself matches the usual one
        assert_eq!(
            compressed_point,
            bincode::deserialize(&compressed_point.0).unwrap()
        );

        let uncompressed_point = point.to_uncompressed();
        let encoded = bincode::serialize(&uncompressed_point).unwrap();
        let parsed: UncompressedPoint = bincode::deserialize(&encoded).unwrap();
        assert_eq!(parsed, uncompressed_point);

        // Check that the encoding is 97 bytes exactly
        assert_eq!(encoded.len(), 97);

        // Check that the encoding itself matches the usual one
        assert_eq!(
            uncompressed_point,
            bincode::deserialize(&uncompressed_point.0).unwrap()
        );

        // Check that invalid encodings fail
        let encoded = bincode::serialize(&compressed_point).unwrap();
        assert!(bincode::deserialize::<CompressedPoint>(&encoded[0..47]).is_err());
        let encoded = bincode::serialize(&uncompressed_point).unwrap();
        assert!(bincode::deserialize::<UncompressedPoint>(&encoded[0..47]).is_err());
    }
}
