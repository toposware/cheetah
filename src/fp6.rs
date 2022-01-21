// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This module implements arithmetic over the extension field Fp6,
//! defined with irreducible polynomial v^3 + v + 1.

use core::fmt::{self, Formatter};
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use group::ff::Field;
use rand_core::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "serialize")]
use serde::de::Visitor;
#[cfg(feature = "serialize")]
use serde::{self, Deserialize, Deserializer, Serialize, Serializer};

use crate::fp::Fp;
use crate::fp2::Fp2;
use crate::utils::square_assign_multi;

use crate::fp::TWO_ADICITY;

// 2^56 root of unity = 3537142827676469596*u + 1182629581807810213
//                    = 1992083664070583229*u + 2727688745413696580 in Montgomery form
const TWO_ADIC_ROOT_OF_UNITY_P6: Fp6 = Fp6 {
    c0: Fp2 {
        c0: Fp(0x25dab2764836dc44),
        c1: Fp(0x1ba54d89b7c923bd),
    },
    c1: Fp2::zero(),
    c2: Fp2::zero(),
};

const TWO_ADICITY_P6: u32 = TWO_ADICITY + 1;

#[derive(Copy, Clone)]
/// An element of the extension GF(p^6)
pub struct Fp6 {
    /// First coefficient
    pub c0: Fp2,
    /// Second coefficient
    pub c1: Fp2,
    /// Third coefficient
    pub c2: Fp2,
}

impl fmt::Debug for Fp6 {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{:?} + ({:?})*v + ({:?})*v^2", self.c0, self.c1, self.c2)
    }
}

impl Default for Fp6 {
    fn default() -> Self {
        Fp6::zero()
    }
}

impl zeroize::DefaultIsZeroes for Fp6 {}

impl From<Fp2> for Fp6 {
    fn from(f: Fp2) -> Self {
        Fp6 {
            c0: f,
            c1: Fp2::zero(),
            c2: Fp2::zero(),
        }
    }
}

impl From<[Fp2; 3]> for Fp6 {
    fn from(f: [Fp2; 3]) -> Self {
        Fp6 {
            c0: f[0],
            c1: f[1],
            c2: f[2],
        }
    }
}

impl From<Fp> for Fp6 {
    fn from(f: Fp) -> Self {
        Fp6::from(Fp2::from(f))
    }
}

impl From<[Fp; 6]> for Fp6 {
    fn from(f: [Fp; 6]) -> Self {
        Fp6 {
            c0: Fp2 { c0: f[0], c1: f[1] },
            c1: Fp2 { c0: f[2], c1: f[3] },
            c2: Fp2 { c0: f[4], c1: f[5] },
        }
    }
}

impl From<Fp6> for [Fp; 6] {
    fn from(f: Fp6) -> [Fp; 6] {
        [f.c0.c0, f.c0.c1, f.c1.c0, f.c1.c1, f.c2.c0, f.c2.c1]
    }
}

impl From<u64> for Fp6 {
    /// Converts a 64-bit value into a field element. If the value is greater than or equal to
    /// the field modulus, modular reduction is silently performed.
    fn from(value: u64) -> Self {
        Fp6::from(Fp::new(value))
    }
}

impl From<u32> for Fp6 {
    /// Converts a 32-bit value into a field element.
    fn from(value: u32) -> Self {
        Fp6::from(Fp::new(value as u64))
    }
}

impl From<u16> for Fp6 {
    /// Converts a 16-bit value into a field element.
    fn from(value: u16) -> Self {
        Fp6::from(Fp::new(value as u64))
    }
}

impl From<u8> for Fp6 {
    /// Converts an 8-bit value into a field element.
    fn from(value: u8) -> Self {
        Fp6::from(Fp::new(value as u64))
    }
}

impl ConstantTimeEq for Fp6 {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.c0.ct_eq(&other.c0) & self.c1.ct_eq(&other.c1) & self.c2.ct_eq(&other.c2)
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
        Fp6 {
            c0: Fp2::conditional_select(&a.c0, &b.c0, choice),
            c1: Fp2::conditional_select(&a.c1, &b.c1, choice),
            c2: Fp2::conditional_select(&a.c2, &b.c2, choice),
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
    /// The value is converted to Montgomery form by computing
    /// (a.R^0 * R^2) / R = a.R
    pub const fn new(value: [u64; 6]) -> Self {
        Fp6 {
            c0: Fp2 {
                c0: Fp::new(value[0]),
                c1: Fp::new(value[1]),
            },
            c1: Fp2 {
                c0: Fp::new(value[2]),
                c1: Fp::new(value[3]),
            },
            c2: Fp2 {
                c0: Fp::new(value[4]),
                c1: Fp::new(value[5]),
            },
        }
    }

    #[inline]
    /// The additive identity
    pub const fn zero() -> Self {
        Fp6 {
            c0: Fp2::zero(),
            c1: Fp2::zero(),
            c2: Fp2::zero(),
        }
    }

    #[inline]
    /// The multiplicative identity
    pub const fn one() -> Self {
        Fp6 {
            c0: Fp2::one(),
            c1: Fp2::zero(),
            c2: Fp2::zero(),
        }
    }

    /// Checks whether this element is zero or not
    pub fn is_zero(&self) -> Choice {
        self.c0.is_zero() & self.c1.is_zero() & self.c2.is_zero()
    }

    /// Returns whether or not this element is strictly lexicographically
    /// larger than its negation.
    #[inline]
    pub fn lexicographically_largest(&self) -> Choice {
        // If this element's c2 coefficient is lexicographically largest
        // then it is lexicographically largest. In the event
        // the c2 coefficient is zero and the c1 coefficient is
        // lexicographically largest, then this element is lexicographically
        // largest. Otherwise, in the event both the c2 and c1 coefficients
        // are zero and the c0 coefficient is lexicographically largest,
        // then this element is lexicographically largest.

        self.c2.lexicographically_largest()
            | (self.c2.is_zero() & self.c1.lexicographically_largest())
            | (self.c1.is_zero() & self.c0.lexicographically_largest())
    }

    #[inline]
    /// Computes the multiplication of two Fp6 elements
    pub const fn mul(&self, other: &Fp6) -> Fp6 {
        let aa = (&self.c0).mul(&other.c0);
        let bb = (&self.c1).mul(&other.c1);
        let cc = (&self.c2).mul(&other.c2);

        let ab_ab = (&self.c0).add(&self.c1);
        let tmp = (&other.c0).add(&other.c1);
        let ab_ab = (&ab_ab).mul(&tmp);

        let ac_ac = (&self.c0).add(&self.c2);
        let tmp = (&other.c0).add(&other.c2);
        let ac_ac = (&ac_ac).mul(&tmp);

        let bc_bc = (&self.c1).add(&self.c2);
        let tmp = (&other.c1).add(&other.c2);
        let bc_bc = (&bc_bc).mul(&tmp);

        let tmp = (&aa).add(&bb);
        let tmp = (&tmp).add(&cc);

        let c0 = (&tmp).sub(&bc_bc);

        let c1 = (&ab_ab).sub(&bc_bc);
        let c1 = (&c1).sub(&aa);

        let c2 = (&ac_ac).sub(&tmp);
        let c2 = (&c2).sub(&cc);
        let t2 = (&bb).double();
        let c2 = (&c2).add(&t2);

        Fp6 { c0, c1, c2 }
    }

    /// Computes the square of a field element
    #[inline]
    pub const fn square(&self) -> Self {
        let aa = (&self.c0).square();
        let bb = (&self.c1).square();
        let cc = (&self.c2).square();

        let ab_ab = (&self.c0).add(&self.c1);
        let ab_ab = (&ab_ab).square();

        let ac_ac = (&self.c0).add(&self.c2);
        let ac_ac = (&ac_ac).square();

        let bc_bc = (&self.c1).add(&self.c2);
        let bc_bc = (&bc_bc).square();

        let tmp = (&aa).add(&bb);
        let tmp = (&tmp).add(&cc);

        let c0 = (&tmp).sub(&bc_bc);

        let c1 = (&ab_ab).sub(&bc_bc);
        let c1 = (&c1).sub(&aa);

        let c2 = (&ac_ac).sub(&tmp);
        let c2 = (&c2).sub(&cc);
        let t2 = (&bb).double();
        let c2 = (&c2).add(&t2);

        Fp6 { c0, c1, c2 }
    }

    /// Computes the square root of this element, if it exists.
    pub fn sqrt(&self) -> CtOption<Self> {
        // Tonelli-Shank's algorithm for q mod 16 = 1
        // See https://eprint.iacr.org/2020/1497.pdf, page 3 for a
        // constant time specification of the algorithm.

        // We construct t and z from Fp6 seen as a direct extension
        // of Fp, and then use the isomorphism in `util.sage` to go
        // back to the tower extension.

        // Compute the progenitor y of self
        // y = self^((t - 1) // 2)
        //   = self^0x931696537b20000d7946712c80000083a6d75f80000002ae101c0000000007db0e00000000000c4
        let y = self.exp_vartime(&[
            0xb0e00000000000c4,
            0xe101c0000000007d,
            0x3a6d75f80000002a,
            0xd7946712c8000008,
            0x0931696537b20000,
            0x0000000000000000,
        ]);

        let mut s = self * y;
        let mut t = s * y;

        let mut z = TWO_ADIC_ROOT_OF_UNITY_P6;

        for k in (2..=TWO_ADICITY_P6).rev() {
            let mut b = t;

            square_assign_multi(&mut b, (k - 2) as usize);

            let new_s = s * z;
            s = Fp6::conditional_select(&new_s, &s, b.ct_eq(&Fp6::one()));
            z = z.square();
            let new_t = t * z;
            t = Fp6::conditional_select(&new_t, &t, b.ct_eq(&Fp6::one()));
        }

        CtOption::new(s, (s * s).ct_eq(self))
    }

    /// Computes the double of a field element
    #[inline]
    pub const fn double(&self) -> Self {
        Fp6 {
            c0: (&self.c0).double(),
            c1: (&self.c1).double(),
            c2: (&self.c2).double(),
        }
    }

    /// Computes the summation of two field elements
    #[inline]
    pub const fn add(&self, rhs: &Self) -> Self {
        Fp6 {
            c0: (&self.c0).add(&rhs.c0),
            c1: (&self.c1).add(&rhs.c1),
            c2: (&self.c2).add(&rhs.c2),
        }
    }

    /// Computes the difference of two field elements
    #[inline]
    pub const fn sub(&self, rhs: &Self) -> Self {
        Fp6 {
            c0: (&self.c0).sub(&rhs.c0),
            c1: (&self.c1).sub(&rhs.c1),
            c2: (&self.c2).sub(&rhs.c2),
        }
    }

    /// Computes the negation of a field element
    #[inline]
    pub const fn neg(&self) -> Self {
        Fp6 {
            c0: (&self.c0).neg(),
            c1: (&self.c1).neg(),
            c2: (&self.c2).neg(),
        }
    }

    /// Computes the multiplicative inverse of this field
    /// element, returning None in the case that this element
    /// is zero.
    pub fn invert(&self) -> CtOption<Self> {
        let c0_sq = self.c0.square();
        let c1_sq = self.c1.square();
        let c2_sq = self.c2.square();

        let inv = self.c0 * (c0_sq + c1_sq) - self.c1 * c1_sq
            + (self.c0 - self.c1 + self.c2) * c2_sq
            - (c0_sq.double() - (self.c0.double() + self.c0) * self.c1) * self.c2;

        let c0 = c0_sq + c1_sq - (self.c0.double() - self.c1) * self.c2 + c2_sq;
        let c1 = -(self.c0 * self.c1 + c2_sq);
        let c2 = c1_sq - self.c0 * self.c2 + c2_sq;

        inv.invert().map(|t| Fp6 {
            c0: c0 * t,
            c1: c1 * t,
            c2: c2 * t,
        })
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

    /// Outputs the internal representation as 6 64-bit limbs after Montgomery reduction
    pub const fn output_reduced_limbs(&self) -> [u64; 6] {
        [
            Fp::montgomery_reduce(self.c0.c0.0, 0).0,
            Fp::montgomery_reduce(self.c0.c1.0, 0).0,
            Fp::montgomery_reduce(self.c1.c0.0, 0).0,
            Fp::montgomery_reduce(self.c1.c1.0, 0).0,
            Fp::montgomery_reduce(self.c2.c0.0, 0).0,
            Fp::montgomery_reduce(self.c2.c1.0, 0).0,
        ]
    }

    /// Outputs the internal representation as 6 64-bit limbs without Montgomery reduction
    /// This is intended for uses like re-interpreting the type containing the internal value.
    pub const fn output_unreduced_limbs(&self) -> [u64; 6] {
        [
            self.c0.c0.0,
            self.c0.c1.0,
            self.c1.c0.0,
            self.c1.c1.0,
            self.c2.c0.0,
            self.c2.c1.0,
        ]
    }

    /// Converts an `Fp6` element into a byte representation in
    /// little-endian byte order.
    pub fn to_bytes(&self) -> [u8; 48] {
        let mut bytes = [0u8; 48];

        bytes[0..16].copy_from_slice(&self.c0.to_bytes());
        bytes[16..32].copy_from_slice(&self.c1.to_bytes());
        bytes[32..48].copy_from_slice(&self.c2.to_bytes());

        bytes
    }

    /// Attempts to convert a little-endian byte representation of
    /// a scalar into a `Fp6` element, failing if the input is not canonical.
    pub fn from_bytes(bytes: &[u8; 48]) -> CtOption<Self> {
        let mut array = [0u8; 16];

        array.copy_from_slice(&bytes[0..16]);
        let c0 = Fp2::from_bytes(&array);

        array.copy_from_slice(&bytes[16..32]);
        let c1 = Fp2::from_bytes(&array);

        array.copy_from_slice(&bytes[32..48]);
        let c2 = Fp2::from_bytes(&array);

        let is_some = c0.is_some() & c1.is_some() & c2.is_some();

        CtOption::new(
            Fp6 {
                c0: c0.unwrap_or(Fp2::zero()),
                c1: c1.unwrap_or(Fp2::zero()),
                c2: c2.unwrap_or(Fp2::zero()),
            },
            is_some,
        )
    }

    /// Constructs an element of `Fp6` without checking that it is
    /// canonical.
    pub const fn from_raw_unchecked(value: [u64; 6]) -> Self {
        Fp6 {
            c0: Fp2 {
                c0: Fp::from_raw_unchecked(value[0]),
                c1: Fp::from_raw_unchecked(value[1]),
            },
            c1: Fp2 {
                c0: Fp::from_raw_unchecked(value[2]),
                c1: Fp::from_raw_unchecked(value[3]),
            },
            c2: Fp2 {
                c0: Fp::from_raw_unchecked(value[4]),
                c1: Fp::from_raw_unchecked(value[5]),
            },
        }
    }
}

// FIELD TRAITS IMPLEMENTATION
// ================================================================================================

impl Field for Fp6 {
    fn random(mut rng: impl RngCore) -> Self {
        Fp6 {
            c0: Fp2::random(&mut rng),
            c1: Fp2::random(&mut rng),
            c2: Fp2::random(&mut rng),
        }
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
mod test {
    use super::*;
    use rand_core::OsRng;

    // DISPLAY
    // ================================================================================================

    #[test]
    fn test_debug() {
        assert_eq!(
            format!("{:?}", Fp6::zero()),
            "0 + 0*u + (0 + 0*u)*v + (0 + 0*u)*v^2"
        );
        assert_eq!(
            format!("{:?}", Fp6::one()),
            "1 + 0*u + (0 + 0*u)*v + (0 + 0*u)*v^2"
        );
        assert_eq!(
            format!("{:?}", Fp6::new([1, 2, 3, 4, 5, 6])),
            "1 + 2*u + (3 + 4*u)*v + (5 + 6*u)*v^2"
        );

        let a = Fp6::one().neg();
        assert_eq!(
            format!("{:?}", a),
            "4719772409484279808 + 0*u + (0 + 0*u)*v + (0 + 0*u)*v^2"
        );
    }

    #[test]
    fn test_output_reduced_limbs() {
        assert_eq!(
            format!("{:?}", Fp6::zero().output_reduced_limbs()),
            "[0, 0, 0, 0, 0, 0]"
        );
        assert_eq!(
            format!("{:?}", Fp6::one().output_reduced_limbs()),
            "[1, 0, 0, 0, 0, 0]"
        );
        let a = Fp6::one().neg();
        assert_eq!(
            format!("{:?}", a.output_reduced_limbs()),
            "[4719772409484279808, 0, 0, 0, 0, 0]"
        );
    }

    // BASIC ALGEBRA
    // ================================================================================================

    #[test]
    fn test_conditional_selection() {
        let a = Fp6 {
            c0: Fp2 {
                c0: Fp::one(),
                c1: Fp::new(2),
            },
            c1: Fp2 {
                c0: Fp::new(3),
                c1: Fp::new(4),
            },
            c2: Fp2 {
                c0: Fp::new(5),
                c1: Fp::new(6),
            },
        };
        let b = Fp6 {
            c0: Fp2 {
                c0: Fp::new(7),
                c1: Fp::new(8),
            },
            c1: Fp2 {
                c0: Fp::new(9),
                c1: Fp::new(10),
            },
            c2: Fp2 {
                c0: Fp::new(11),
                c1: Fp::new(12),
            },
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
                c0: Fp2 {
                    c0: Fp::one(),
                    c1: Fp::new(2),
                },
                c1: Fp2 {
                    c0: Fp::new(3),
                    c1: Fp::new(4),
                },
                c2: Fp2 {
                    c0: Fp::new(5),
                    c1: Fp::new(6),
                },
            },
            &Fp6 {
                c0: Fp2 {
                    c0: Fp::one(),
                    c1: Fp::new(2),
                },
                c1: Fp2 {
                    c0: Fp::new(3),
                    c1: Fp::new(4),
                },
                c2: Fp2 {
                    c0: Fp::new(5),
                    c1: Fp::new(6),
                },
            }
        ));

        assert!(!is_equal(
            &Fp6 {
                c0: Fp2 {
                    c0: Fp::new(2),
                    c1: Fp::new(2),
                },
                c1: Fp2 {
                    c0: Fp::new(3),
                    c1: Fp::new(4),
                },
                c2: Fp2 {
                    c0: Fp::new(5),
                    c1: Fp::new(6),
                },
            },
            &Fp6 {
                c0: Fp2 {
                    c0: Fp::one(),
                    c1: Fp::new(2),
                },
                c1: Fp2 {
                    c0: Fp::new(3),
                    c1: Fp::new(4),
                },
                c2: Fp2 {
                    c0: Fp::new(5),
                    c1: Fp::new(6),
                },
            }
        ));

        assert!(!is_equal(
            &Fp6 {
                c0: Fp2 {
                    c0: Fp::one(),
                    c1: Fp::new(2),
                },
                c1: Fp2 {
                    c0: Fp::new(4),
                    c1: Fp::new(4),
                },
                c2: Fp2 {
                    c0: Fp::new(5),
                    c1: Fp::new(6),
                },
            },
            &Fp6 {
                c0: Fp2 {
                    c0: Fp::one(),
                    c1: Fp::new(2),
                },
                c1: Fp2 {
                    c0: Fp::new(3),
                    c1: Fp::new(4),
                },
                c2: Fp2 {
                    c0: Fp::new(5),
                    c1: Fp::new(6),
                },
            }
        ));

        assert!(!is_equal(
            &Fp6 {
                c0: Fp2 {
                    c0: Fp::one(),
                    c1: Fp::new(2),
                },
                c1: Fp2 {
                    c0: Fp::new(3),
                    c1: Fp::new(4),
                },
                c2: Fp2 {
                    c0: Fp::new(6),
                    c1: Fp::new(6),
                },
            },
            &Fp6 {
                c0: Fp2 {
                    c0: Fp::one(),
                    c1: Fp::new(2),
                },
                c1: Fp2 {
                    c0: Fp::new(3),
                    c1: Fp::new(4),
                },
                c2: Fp2 {
                    c0: Fp::new(5),
                    c1: Fp::new(6),
                },
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
            c0: Fp2 {
                c0: Fp::new(2888944236454241314),
                c1: Fp::new(1900961609923933690),
            },
            c1: Fp2 {
                c0: Fp::new(4164134277002348062),
                c1: Fp::new(307346200243481763),
            },
            c2: Fp2 {
                c0: Fp::new(2163980486023891318),
                c1: Fp::new(3902966736580914535),
            },
        };
        let b = Fp6 {
            c0: Fp2 {
                c0: Fp::new(3125384363659457274),
                c1: Fp::new(3501855459985374183),
            },
            c1: Fp2 {
                c0: Fp::new(4665660136548904132),
                c1: Fp::new(926330806566119355),
            },
            c2: Fp2 {
                c0: Fp::new(3751072364085233784),
                c1: Fp::new(193238968045832370),
            },
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

        // (1817037221240944654 * u + 1972121079904647667) * v^2
        //      + (4179895386176115723 * u + 372058458154352541) * v
        //      + 799979765245758940 * u + 3901586259568969690
        // is not a quadratic residue in Fp6
        assert!(bool::from(
            Fp6 {
                c0: Fp2 {
                    c0: Fp::new(3901586259568969690),
                    c1: Fp::new(799979765245758940),
                },
                c1: Fp2 {
                    c0: Fp::new(372058458154352541),
                    c1: Fp::new(4179895386176115723),
                },
                c2: Fp2 {
                    c0: Fp::new(1972121079904647667),
                    c1: Fp::new(1817037221240944654),
                },
            }
            .sqrt()
            .is_none()
        ));
    }

    #[test]
    fn test_multiplication() {
        let a = Fp6 {
            c0: Fp2 {
                c0: Fp::one(),
                c1: Fp::new(2),
            },
            c1: Fp2 {
                c0: Fp::new(3),
                c1: Fp::new(4),
            },
            c2: Fp2 {
                c0: Fp::new(5),
                c1: Fp::new(6),
            },
        };
        let b = Fp6::one();
        let c = Fp6 {
            c0: Fp2 {
                c0: Fp::one(),
                c1: Fp::new(2),
            },
            c1: Fp2 {
                c0: Fp::new(3),
                c1: Fp::new(4),
            },
            c2: Fp2 {
                c0: Fp::new(5),
                c1: Fp::new(6),
            },
        };

        assert_eq!(a * b, c);

        let a = Fp6 {
            c0: Fp2 {
                c0: Fp::new(2888944236454241314),
                c1: Fp::new(1900961609923933690),
            },
            c1: Fp2 {
                c0: Fp::new(4164134277002348062),
                c1: Fp::new(307346200243481763),
            },
            c2: Fp2 {
                c0: Fp::new(2163980486023891318),
                c1: Fp::new(3902966736580914535),
            },
        };
        let b = Fp6 {
            c0: Fp2 {
                c0: Fp::new(2535997426787301359),
                c1: Fp::new(3971257432982727883),
            },
            c1: Fp2 {
                c0: Fp::new(4562449940607053508),
                c1: Fp::new(2369569849810406504),
            },
            c2: Fp2 {
                c0: Fp::new(4379339381077658533),
                c1: Fp::new(2874815267585070886),
            },
        };
        let c = Fp6 {
            c0: Fp2 {
                c0: Fp::new(3981622945204561084),
                c1: Fp::new(152859226758736115),
            },
            c1: Fp2 {
                c0: Fp::new(707484124120067374),
                c1: Fp::new(4240778718972633208),
            },
            c2: Fp2 {
                c0: Fp::new(136693932767651254),
                c1: Fp::new(174194258731295370),
            },
        };

        assert_eq!(a * b, c);
    }

    #[test]
    fn test_addition() {
        let a = Fp6 {
            c0: Fp2 {
                c0: Fp::one(),
                c1: Fp::new(2),
            },
            c1: Fp2 {
                c0: Fp::new(3),
                c1: Fp::new(4),
            },
            c2: Fp2 {
                c0: Fp::new(5),
                c1: Fp::new(6),
            },
        };
        let b = Fp6 {
            c0: Fp2 {
                c0: Fp::new(6),
                c1: Fp::new(5),
            },
            c1: Fp2 {
                c0: Fp::new(4),
                c1: Fp::new(3),
            },
            c2: Fp2 {
                c0: Fp::new(2),
                c1: Fp::one(),
            },
        };
        let c = Fp6 {
            c0: Fp2 {
                c0: Fp::new(7),
                c1: Fp::new(7),
            },
            c1: Fp2 {
                c0: Fp::new(7),
                c1: Fp::new(7),
            },
            c2: Fp2 {
                c0: Fp::new(7),
                c1: Fp::new(7),
            },
        };

        assert_eq!(a + b, c);
    }

    #[test]
    fn test_subtraction() {
        let a = Fp6 {
            c0: Fp2 {
                c0: Fp::new(6),
                c1: Fp::new(5),
            },
            c1: Fp2 {
                c0: Fp::new(4),
                c1: Fp::new(3),
            },
            c2: Fp2 {
                c0: Fp::new(2),
                c1: Fp::one(),
            },
        };
        let b = Fp6 {
            c0: Fp2 {
                c0: Fp::new(3),
                c1: Fp::new(3),
            },
            c1: Fp2 {
                c0: Fp::new(2),
                c1: Fp::new(2),
            },
            c2: Fp2 {
                c0: Fp::one(),
                c1: Fp::one(),
            },
        };
        let c = Fp6 {
            c0: Fp2 {
                c0: Fp::new(3),
                c1: Fp::new(2),
            },
            c1: Fp2 {
                c0: Fp::new(2),
                c1: Fp::new(1),
            },
            c2: Fp2 {
                c0: Fp::new(1),
                c1: Fp::new(0),
            },
        };

        assert_eq!(a - b, c);
    }

    #[test]
    fn test_negation() {
        let mut rng = OsRng;

        let a = Fp6 {
            c0: Fp2 {
                c0: Fp::one(),
                c1: Fp::new(2),
            },
            c1: Fp2 {
                c0: Fp::new(3),
                c1: Fp::new(4),
            },
            c2: Fp2 {
                c0: Fp::new(5),
                c1: Fp::new(6),
            },
        };
        let b = Fp6 {
            c0: -Fp2 {
                c0: Fp::one(),
                c1: Fp::new(2),
            },
            c1: -Fp2 {
                c0: Fp::new(3),
                c1: Fp::new(4),
            },
            c2: -Fp2 {
                c0: Fp::new(5),
                c1: Fp::new(6),
            },
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
            c0: Fp2 {
                c0: Fp::new(2888944236454241314),
                c1: Fp::new(1900961609923933690),
            },
            c1: Fp2 {
                c0: Fp::new(4164134277002348062),
                c1: Fp::new(307346200243481763),
            },
            c2: Fp2 {
                c0: Fp::new(2163980486023891318),
                c1: Fp::new(3902966736580914535),
            },
        };

        let b = Fp6 {
            c0: Fp2 {
                c0: Fp::new(411478125972401294),
                c1: Fp::new(2671067376479186442),
            },
            c1: Fp2 {
                c0: Fp::new(4339807523281293225),
                c1: Fp::new(4431608856960589867),
            },
            c2: Fp2 {
                c0: Fp::new(1673850528796443186),
                c1: Fp::new(533572733765271776),
            },
        };

        assert_eq!(a.invert().unwrap(), b);

        assert_eq!(Fp6::one().invert().unwrap(), Fp6::one());

        assert!(bool::from(Fp6::zero().invert().is_none()));
    }

    #[test]
    fn test_invert_is_pow() {
        let mut rng = OsRng;

        let p6_minus_2 = [
            0x88ffffffffffffff,
            0xfb61c00000000001,
            0x55c2038000000000,
            0x1074daebf0000000,
            0x01af28ce25900000,
            0x001262d2ca6f6400,
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
        let two_pow_56 = 1 << TWO_ADICITY_P6 as u64;
        assert_eq!(
            Fp6::one(),
            TWO_ADIC_ROOT_OF_UNITY_P6.exp(&[two_pow_56, 0, 0, 0, 0, 0])
        );
        assert_ne!(
            Fp6::one(),
            TWO_ADIC_ROOT_OF_UNITY_P6.exp(&[two_pow_56 - 1, 0, 0, 0, 0, 0])
        );
    }

    #[test]
    fn test_lexicographic_largest() {
        assert!(!bool::from(Fp6::zero().lexicographically_largest()));
        assert!(!bool::from(Fp6::one().lexicographically_largest()));
        assert!(bool::from(
            Fp6 {
                c0: Fp2 {
                    c0: Fp::new(4167072808029173087),
                    c1: Fp::new(1618085398724889560),
                },
                c1: Fp2 {
                    c0: Fp::new(4252814005348531461),
                    c1: Fp::new(1978937931976521722),
                },
                c2: Fp2 {
                    c0: Fp::new(733612977188794891),
                    c1: Fp::new(2521078467018751009),
                },
            }
            .lexicographically_largest()
        ));
        assert!(!bool::from(
            Fp6 {
                c0: Fp2 {
                    c0: Fp::new(444552187502873250),
                    c1: Fp::new(2993539596807156777),
                },
                c1: Fp2 {
                    c0: Fp::new(358810990183514876),
                    c1: Fp::new(2632687063555524615),
                },
                c2: Fp2 {
                    c0: Fp::new(3878012018343251446),
                    c1: Fp::new(2090546528513295328),
                },
            }
            .lexicographically_largest()
        ));
        assert!(bool::from(
            Fp6 {
                c0: Fp2 {
                    c0: Fp::new(444552187502873250),
                    c1: Fp::new(2993539596807156777),
                },
                c1: Fp2::zero(),
                c2: Fp2::zero(),
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
    fn test_from_raw_unchecked() {
        let mut element = Fp6::from_raw_unchecked([4287426845256712189, 0, 0, 0, 0, 0]);

        let element_normalized = Fp6::new([4287426845256712189, 0, 0, 0, 0, 0]);

        assert_eq!(element, Fp6::one());
        element *= &crate::fp::R2.into();

        assert!(element != Fp6::one());
        assert_eq!(element, element_normalized);
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

    #[test]
    fn test_from_fp2() {
        let mut rng = OsRng;
        let v0 = rng.next_u64();
        let v1 = rng.next_u64();
        let e = Fp2::new([v0, v1]);

        let e_fp6 = Fp6::new([v0, v1, 0, 0, 0, 0]);
        let array: [Fp2; 3] = [e, Fp2::zero(), Fp2::zero()];

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
                c0: Fp2 {
                    c0: Fp::new(3901586259568969690),
                    c1: Fp::new(799979765245758940),
                },
                c1: Fp2 {
                    c0: Fp::new(372058458154352541),
                    c1: Fp::new(4179895386176115723),
                },
                c2: Fp2 {
                    c0: Fp::new(1972121079904647667),
                    c1: Fp::new(1817037221240944654),
                },
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
                0, 0, 0, 0, 0, 0, 128, 65, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
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
                0, 0, 0, 0, 0, 0, 128, 65, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
            .unwrap(),
            -Fp6::one()
        );

        // Anything larger than M in one of the members is invalid
        assert!(bool::from(
            Fp6::from_bytes(&[
                2, 0, 0, 0, 0, 0, 128, 65, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
            .is_none()
        ));

        assert!(bool::from(
            Fp6::from_bytes(&[
                0, 0, 0, 0, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
            .is_none()
        ));

        assert!(bool::from(
            Fp6::from_bytes(&[
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
            .is_none()
        ));

        assert!(bool::from(
            Fp6::from_bytes(&[
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
            .is_none()
        ));

        assert!(bool::from(
            Fp6::from_bytes(&[
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
            .is_none()
        ));

        assert!(bool::from(
            Fp6::from_bytes(&[
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0
            ])
            .is_none()
        ));

        assert!(bool::from(
            Fp6::from_bytes(&[
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255
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
        let element = Fp6::random(&mut rng);
        let mut encoded = bincode::serialize(&element).unwrap();
        encoded[47] = 127;
        assert!(bincode::deserialize::<Fp6>(&encoded).is_err());

        let encoded = bincode::serialize(&element).unwrap();
        assert!(bincode::deserialize::<Fp6>(&encoded[0..47]).is_err());
    }
}
