// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This module implements arithmetic over the extension field Fp2,
//! defined with irreducible polynomial u^2 - 2x - 2.

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
use crate::fp::TWO_ADICITY;
use crate::utils::square_assign_multi;

// 2^56 root of unity = 3151433596914781827*u + 1568338812569497982
//                    = 1770206186931184045*u + 2949566222553095764 in Montgomery form
const TWO_ADIC_ROOT_OF_UNITY_P2: Fp2 = Fp2 {
    c0: Fp(0x28eef6da18c49e54),
    c1: Fp(0x18910925e73b61ad),
};

const TWO_ADICITY_P2: u32 = TWO_ADICITY + 1;

#[derive(Copy, Clone)]
/// An element of the extension GF(p^2)
pub struct Fp2 {
    /// First coefficient
    pub c0: Fp,
    /// Second coefficient
    pub c1: Fp,
}

impl fmt::Debug for Fp2 {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{:?} + {:?}*u", self.c0, self.c1)
    }
}

impl Default for Fp2 {
    fn default() -> Self {
        Fp2::zero()
    }
}

impl zeroize::DefaultIsZeroes for Fp2 {}

impl From<Fp> for Fp2 {
    fn from(f: Fp) -> Self {
        Fp2 {
            c0: f,
            c1: Fp::zero(),
        }
    }
}

impl From<[Fp; 2]> for Fp2 {
    fn from(f: [Fp; 2]) -> Self {
        Fp2 { c0: f[0], c1: f[1] }
    }
}

impl From<Fp2> for [Fp; 2] {
    fn from(f: Fp2) -> [Fp; 2] {
        [f.c0, f.c1]
    }
}

impl ConstantTimeEq for Fp2 {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.c0.ct_eq(&other.c0) & self.c1.ct_eq(&other.c1)
    }
}

impl Eq for Fp2 {}
impl PartialEq for Fp2 {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl ConditionallySelectable for Fp2 {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Fp2 {
            c0: Fp::conditional_select(&a.c0, &b.c0, choice),
            c1: Fp::conditional_select(&a.c1, &b.c1, choice),
        }
    }
}

impl<'a> Neg for &'a Fp2 {
    type Output = Fp2;

    #[inline]
    fn neg(self) -> Fp2 {
        self.neg()
    }
}

impl Neg for Fp2 {
    type Output = Fp2;

    #[inline]
    fn neg(self) -> Fp2 {
        -&self
    }
}

impl<'a, 'b> Sub<&'b Fp2> for &'a Fp2 {
    type Output = Fp2;

    #[inline]
    fn sub(self, rhs: &'b Fp2) -> Fp2 {
        self.sub(rhs)
    }
}

impl<'a, 'b> Add<&'b Fp2> for &'a Fp2 {
    type Output = Fp2;

    #[inline]
    fn add(self, rhs: &'b Fp2) -> Fp2 {
        self.add(rhs)
    }
}

impl<'a, 'b> Mul<&'b Fp2> for &'a Fp2 {
    type Output = Fp2;

    #[inline]
    fn mul(self, rhs: &'b Fp2) -> Fp2 {
        self.mul(rhs)
    }
}

impl_binops_additive!(Fp2, Fp2);
impl_binops_multiplicative!(Fp2, Fp2);

impl Fp2 {
    /// Creates a new field element from a `u64` value.
    /// The value is converted to Montgomery form by computing
    /// (a.R^0 * R^2) / R = a.R
    pub const fn new(value: [u64; 2]) -> Self {
        Fp2 {
            c0: Fp::new(value[0]),
            c1: Fp::new(value[1]),
        }
    }
    #[inline]
    /// The additive identity
    pub const fn zero() -> Self {
        Fp2 {
            c0: Fp::zero(),
            c1: Fp::zero(),
        }
    }

    #[inline]
    /// The multiplicative identity
    pub const fn one() -> Self {
        Fp2 {
            c0: Fp::one(),
            c1: Fp::zero(),
        }
    }

    /// Checks whether this element is zero or not
    pub fn is_zero(&self) -> Choice {
        self.c0.is_zero() & self.c1.is_zero()
    }

    /// Returns whether or not this element is strictly lexicographically
    /// larger than its negation.
    #[inline]
    pub fn lexicographically_largest(&self) -> Choice {
        // If this element's c1 coefficient is lexicographically largest
        // then it is lexicographically largest. Otherwise, in the event
        // the c1 coefficient is zero and the c0 coefficient is
        // lexicographically largest, then this element is lexicographically
        // largest.

        self.c1.lexicographically_largest()
            | (self.c1.is_zero() & self.c0.lexicographically_largest())
    }

    #[inline]
    /// Computes the multiplication of two Fp2 elements
    pub const fn mul(&self, other: &Fp2) -> Fp2 {
        let aa = (&self.c0).mul(&other.c0);
        let bb = (&self.c1).mul(&other.c1);

        let tmp = (&self.c0).sub(&self.c1);
        let tmp1 = (&other.c1).sub(&other.c0);
        let tmp = (&tmp).mul(&tmp1);

        let c0 = (&bb).double();
        let c0 = (&c0).add(&aa);

        let c1 = (&bb).add(&c0);
        let c1 = (&c1).add(&tmp);

        Fp2 { c0, c1 }
    }

    /// Computes the square of a field element
    #[inline]
    pub const fn square(&self) -> Self {
        let aa = (&self.c0).square();
        let bb = (&self.c1).square();

        let tmp = (&self.c0).sub(&self.c1);
        let tmp = (&tmp).square();

        let c0 = (&bb).double();
        let c0 = (&c0).add(&aa);

        let c1 = (&bb).add(&c0);
        let c1 = (&c1).sub(&tmp);

        Fp2 { c0, c1 }
    }

    /// Computes the square root of this element, if it exists.
    pub fn sqrt(&self) -> CtOption<Self> {
        // Tonelli-Shank's algorithm for q mod 16 = 1
        // See https://eprint.iacr.org/2020/1497.pdf, page 3 for a
        // constant time specification of the algorithm.

        // Compute the progenitor y of self
        // y = self^((t - 1) // 2)
        //   = self^0x86120000000000041
        let y = self.exp_vartime(&[0x6120000000000041, 0x0000000000000008]);

        let mut s = self * y;
        let mut t = s * y;

        let mut z = TWO_ADIC_ROOT_OF_UNITY_P2;

        for k in (2..=TWO_ADICITY_P2).rev() {
            let mut b = t;

            square_assign_multi(&mut b, (k - 2) as usize);

            let new_s = s * z;
            s = Fp2::conditional_select(&new_s, &s, b.ct_eq(&Fp2::one()));
            z = z.square();
            let new_t = t * z;
            t = Fp2::conditional_select(&new_t, &t, b.ct_eq(&Fp2::one()));
        }

        CtOption::new(s, (s * s).ct_eq(self))
    }

    /// Computes the double of a field element
    #[inline]
    pub const fn double(&self) -> Self {
        Fp2 {
            c0: (&self.c0).double(),
            c1: (&self.c1).double(),
        }
    }

    /// Computes the summation of two field elements
    #[inline]
    pub const fn add(&self, rhs: &Self) -> Self {
        Fp2 {
            c0: (&self.c0).add(&rhs.c0),
            c1: (&self.c1).add(&rhs.c1),
        }
    }

    /// Computes the difference of two field elements
    #[inline]
    pub const fn sub(&self, rhs: &Self) -> Self {
        Fp2 {
            c0: (&self.c0).sub(&rhs.c0),
            c1: (&self.c1).sub(&rhs.c1),
        }
    }

    /// Computes the negation of a field element
    #[inline]
    pub const fn neg(&self) -> Self {
        Fp2 {
            c0: (&self.c0).neg(),
            c1: (&self.c1).neg(),
        }
    }

    /// Computes the multiplicative inverse of this field
    /// element, returning None in the case that this element
    /// is zero.
    pub fn invert(&self) -> CtOption<Self> {
        let inv = self.c0.square() + self.c0.double() * self.c1 - self.c1.square().double();

        inv.invert().map(|t| Fp2 {
            c0: (self.c0 + self.c1.double()) * t,
            c1: (-self.c1) * t,
        })
    }

    /// Exponentiates `self` by `power`, where `power` is a
    /// little-endian order integer exponent.
    pub fn exp(self, by: &[u64; 2]) -> Self {
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
    pub fn exp_vartime(&self, by: &[u64; 2]) -> Self {
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

    /// Outputs the internal representation as 2 64-bit limbs after Montgomery reduction
    pub const fn output_reduced_limbs(&self) -> [u64; 2] {
        [
            Fp::montgomery_reduce(self.c0.0, 0).0,
            Fp::montgomery_reduce(self.c1.0, 0).0,
        ]
    }

    /// Outputs the internal representation as 2 64-bit limbs without Montgomery reduction
    /// This is intended for uses like re-interpreting the type containing the internal value.
    pub const fn output_unreduced_limbs(&self) -> [u64; 2] {
        [self.c0.0, self.c1.0]
    }

    /// Converts an `Fp2` element into a byte representation in
    /// little-endian byte order.
    pub fn to_bytes(&self) -> [u8; 16] {
        let mut bytes = [0u8; 16];

        bytes[0..8].copy_from_slice(&self.c0.to_bytes());
        bytes[8..16].copy_from_slice(&self.c1.to_bytes());

        bytes
    }

    /// Attempts to convert a little-endian byte representation of
    /// a scalar into a `Fp2` element, failing if the input is not canonical.
    pub fn from_bytes(bytes: &[u8; 16]) -> CtOption<Self> {
        let mut array = [0u8; 8];

        array.copy_from_slice(&bytes[0..8]);
        let c0 = Fp::from_bytes(&array);

        array.copy_from_slice(&bytes[8..16]);
        let c1 = Fp::from_bytes(&array);

        let is_some = c0.is_some() & c1.is_some();

        CtOption::new(
            Fp2 {
                c0: c0.unwrap_or(Fp::zero()),
                c1: c1.unwrap_or(Fp::zero()),
            },
            is_some,
        )
    }

    /// Constructs an element of `Fp2` without checking that it is
    /// canonical.
    pub const fn from_raw_unchecked(value: [u64; 2]) -> Self {
        Fp2 {
            c0: Fp::from_raw_unchecked(value[0]),
            c1: Fp::from_raw_unchecked(value[1]),
        }
    }
}

// FIELD TRAITS IMPLEMENTATION
// ================================================================================================

impl Field for Fp2 {
    fn random(mut rng: impl RngCore) -> Self {
        Fp2 {
            c0: Fp::random(&mut rng),
            c1: Fp::random(&mut rng),
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
impl Serialize for Fp2 {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        use serde::ser::SerializeTuple;
        let mut tup = serializer.serialize_tuple(16)?;
        for byte in self.to_bytes().iter() {
            tup.serialize_element(byte)?;
        }
        tup.end()
    }
}

#[cfg(feature = "serialize")]
impl<'de> Deserialize<'de> for Fp2 {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct Fp2Visitor;

        impl<'de> Visitor<'de> for Fp2Visitor {
            type Value = Fp2;

            fn expecting(&self, formatter: &mut Formatter) -> fmt::Result {
                formatter.write_str("a valid field element")
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<Fp2, A::Error>
            where
                A: serde::de::SeqAccess<'de>,
            {
                let mut bytes = [0u8; 16];
                for (i, byte) in bytes.iter_mut().enumerate() {
                    *byte = seq
                        .next_element()?
                        .ok_or_else(|| serde::de::Error::invalid_length(i, &"expected 16 bytes"))?;
                }
                let elem = Fp2::from_bytes(&bytes);
                if bool::from(elem.is_none()) {
                    Err(serde::de::Error::custom("decompression failed"))
                } else {
                    Ok(elem.unwrap())
                }
            }
        }

        deserializer.deserialize_tuple(16, Fp2Visitor)
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
        assert_eq!(format!("{:?}", Fp2::zero()), "0 + 0*u");
        assert_eq!(format!("{:?}", Fp2::one()), "1 + 0*u");
        assert_eq!(
            format!(
                "{:?}",
                Fp2 {
                    c0: Fp::zero(),
                    c1: Fp::one()
                }
            ),
            "0 + 1*u"
        );

        let a = Fp2::one().neg();
        assert_eq!(format!("{:?}", a), "4719772409484279808 + 0*u");
    }

    #[test]
    fn test_output_reduced_limbs() {
        assert_eq!(
            format!("{:?}", Fp2::zero().output_reduced_limbs()),
            "[0, 0]"
        );
        assert_eq!(format!("{:?}", Fp2::one().output_reduced_limbs()), "[1, 0]");
        let a = Fp2::one().neg();
        assert_eq!(
            format!("{:?}", a.output_reduced_limbs()),
            "[4719772409484279808, 0]"
        );
    }

    // BASIC ALGEBRA
    // ================================================================================================

    #[test]
    fn test_conditional_selection() {
        let a = Fp2 {
            c0: Fp::from_raw_unchecked(1),
            c1: Fp::from_raw_unchecked(2),
        };
        let b = Fp2 {
            c0: Fp::from_raw_unchecked(3),
            c1: Fp::from_raw_unchecked(4),
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
        fn is_equal(a: &Fp2, b: &Fp2) -> bool {
            let eq = a == b;
            let ct_eq = a.ct_eq(b);

            assert_eq!(eq, bool::from(ct_eq));

            eq
        }

        assert!(is_equal(
            &Fp2 {
                c0: Fp::from_raw_unchecked(1),
                c1: Fp::from_raw_unchecked(2),
            },
            &Fp2 {
                c0: Fp::from_raw_unchecked(1),
                c1: Fp::from_raw_unchecked(2),
            }
        ));

        assert!(!is_equal(
            &Fp2 {
                c0: Fp::from_raw_unchecked(2),
                c1: Fp::from_raw_unchecked(2),
            },
            &Fp2 {
                c0: Fp::from_raw_unchecked(1),
                c1: Fp::from_raw_unchecked(2),
            }
        ));

        assert!(!is_equal(
            &Fp2 {
                c0: Fp::from_raw_unchecked(1),
                c1: Fp::from_raw_unchecked(3),
            },
            &Fp2 {
                c0: Fp::from_raw_unchecked(1),
                c1: Fp::from_raw_unchecked(2),
            }
        ));

        assert!(bool::from(Fp2::default().is_zero()));
        assert!(!bool::from(Fp2::zero().ct_eq(&Fp2::one())));

        assert_eq!(Fp2::zero(), Fp2::new([0, 0]));
        assert_eq!(Fp2::one(), Fp2::new([1, 0]));
    }

    #[test]
    fn test_squaring() {
        let a = Fp2 {
            c0: Fp::new(1325350471228883001),
            c1: Fp::new(3555613217048211822),
        };
        let b = Fp2 {
            c0: Fp::new(4144385007586913737),
            c1: Fp::new(1006138256509718594),
        };

        assert_eq!(a.square(), b);
    }

    #[test]
    fn test_sqrt() {
        for _ in 0..100 {
            let a = Fp2::random(&mut OsRng).square();
            let b = a.sqrt().unwrap();
            assert_eq!(a, b.square());
        }

        assert_eq!(Fp2::zero().sqrt().unwrap(), Fp2::zero());
        assert_eq!(Fp2::one().sqrt().unwrap(), Fp2::one());

        // u + 2 is not a quadratic residue in Fp2
        assert!(bool::from(Fp2::new([2, 1]).sqrt().is_none()));
    }

    #[test]
    fn test_multiplication() {
        let a = Fp2 {
            c0: Fp::one(),
            c1: Fp::new(2),
        };
        let b = Fp2::one();
        let c = Fp2 {
            c0: Fp::one(),
            c1: Fp::new(2),
        };

        assert_eq!(a * b, c);

        let a = Fp2 {
            c0: Fp::new(702450889089710977),
            c1: Fp::new(131935822767345842),
        };
        let b = Fp2 {
            c0: Fp::new(71395483687675361),
            c1: Fp::new(1751573960731485260),
        };
        let c = Fp2 {
            c0: Fp::new(3870520324615741711),
            c1: Fp::new(2105308272415798137),
        };

        assert_eq!(a * b, c);
    }

    #[test]
    fn test_addition() {
        let a = Fp2 {
            c0: Fp::one(),
            c1: Fp::new(2),
        };
        let b = Fp2 {
            c0: Fp::new(2),
            c1: Fp::one(),
        };
        let c = Fp2 {
            c0: Fp::new(3),
            c1: Fp::new(3),
        };

        assert_eq!(a + b, c);
    }

    #[test]
    fn test_subtraction() {
        let a = Fp2 {
            c0: Fp::new(4),
            c1: Fp::new(3),
        };
        let b = Fp2 {
            c0: Fp::new(2),
            c1: Fp::one(),
        };
        let c = Fp2 {
            c0: Fp::new(2),
            c1: Fp::new(2),
        };

        assert_eq!(a - b, c);
    }

    #[test]
    fn test_negation() {
        let mut rng = OsRng;

        let a = Fp2 {
            c0: Fp::one(),
            c1: Fp::new(2),
        };
        let b = Fp2 {
            c0: -Fp::one(),
            c1: -Fp::new(2),
        };

        assert_eq!(-a, b);

        for _ in 0..100 {
            let a = Fp2::random(&mut rng);
            let b = -a;

            assert_eq!(a + b, Fp2::zero());
        }
    }

    #[test]
    fn test_inversion() {
        let a = Fp2 {
            c0: Fp::new(1325350471228883001),
            c1: Fp::new(3555613217048211822),
        };

        let b = Fp2 {
            c0: Fp::new(2842799053740762838),
            c1: Fp::new(2614512045708308582),
        };

        assert_eq!(a.invert().unwrap(), b);

        assert!(bool::from(Fp2::zero().invert().is_none()));
    }

    #[test]
    fn test_invert_is_pow() {
        let mut rng = OsRng;

        let p2_minus_2 = [0x82ffffffffffffff, 0x10c2400000000000];

        let mut r1 = Fp2::random(&mut rng);
        let mut r2 = r1;
        let mut r3 = r2;

        for _ in 0..100 {
            r1 = r1.invert().unwrap();
            r2 = r2.exp(&p2_minus_2);
            r3 = r3.exp_vartime(&p2_minus_2);

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
        let two_pow_56 = 1 << TWO_ADICITY_P2 as u64;
        assert_eq!(Fp2::one(), TWO_ADIC_ROOT_OF_UNITY_P2.exp(&[two_pow_56, 0]));
        assert_ne!(
            Fp2::one(),
            TWO_ADIC_ROOT_OF_UNITY_P2.exp(&[two_pow_56 - 1, 0])
        );
    }

    #[test]
    fn test_lexicographic_largest() {
        assert!(!bool::from(Fp2::zero().lexicographically_largest()));
        assert!(!bool::from(Fp2::one().lexicographically_largest()));
        assert!(bool::from(
            Fp2 {
                c0: Fp::new(2759934747484735216),
                c1: Fp::new(3406568298137380366),
            }
            .lexicographically_largest()
        ));
        assert!(!bool::from(
            Fp2 {
                c0: Fp::new(1851690248047311121),
                c1: Fp::new(1205056697394665971),
            }
            .lexicographically_largest()
        ));
        assert!(bool::from(
            Fp2 {
                c0: Fp::new(3406568298137380366),
                c1: Fp::zero(),
            }
            .lexicographically_largest()
        ));
    }

    #[test]
    fn test_zeroize() {
        use zeroize::Zeroize;

        let mut a = Fp2::one();
        a.zeroize();
        assert!(bool::from(a.is_zero()));
    }

    #[test]
    fn test_from_raw_unchecked() {
        let mut element = Fp2::from_raw_unchecked([4287426845256712189, 0]);

        let element_normalized = Fp2::new([4287426845256712189, 0]);

        assert_eq!(element, Fp2::one());
        element *= &crate::fp::R2.into();

        assert!(element != Fp2::one());
        assert_eq!(element, element_normalized);
    }

    #[test]
    fn test_from_fp() {
        let mut rng = OsRng;
        let v = rng.next_u64();
        let e = Fp::new(v);

        let e_fp2 = Fp2::new([v, 0]);
        let array: [Fp; 2] = [e, Fp::zero()];

        assert_eq!(e_fp2, e.into());
        assert_eq!(Fp2::from(array), e_fp2);
    }

    // FIELD TRAIT
    // ================================================================================================

    #[test]
    fn test_field_trait_methods() {
        assert_eq!(<Fp2 as Field>::zero(), Fp2::new([0, 0]));
        assert_eq!(<Fp2 as Field>::one(), Fp2::new([1, 0]));

        assert_eq!(
            bool::from(<Fp2 as Field>::zero().is_zero()),
            bool::from(Fp2::new([0, 0]).is_zero())
        );
        assert_eq!(
            bool::from(<Fp2 as Field>::one().is_zero()),
            bool::from(Fp2::new([1, 0]).is_zero())
        );

        let mut rng = OsRng;
        let e = Fp2::random(&mut rng).square();

        assert_eq!(<Fp2 as Field>::square(&e), e.square());
        assert_eq!(<Fp2 as Field>::double(&e), e.double());

        assert_eq!(<Fp2 as Field>::invert(&e).unwrap(), e.invert().unwrap());
        assert!(bool::from(<Fp2 as Field>::invert(&Fp2::zero()).is_none()));

        assert_eq!(<Fp2 as Field>::sqrt(&e).unwrap(), e.sqrt().unwrap());
        assert!(bool::from(
            <Fp2 as Field>::sqrt(&Fp2::new([2, 1])).is_none()
        ));
    }

    // SERIALIZATION / DESERIALIZATION
    // ================================================================================================

    #[test]
    fn test_to_bytes() {
        assert_eq!(
            Fp2::zero().to_bytes(),
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        );

        assert_eq!(
            Fp2::one().to_bytes(),
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        );

        assert_eq!(
            (-&Fp2::one()).to_bytes(),
            [0, 0, 0, 0, 0, 0, 128, 65, 0, 0, 0, 0, 0, 0, 0, 0]
        );
    }

    #[test]
    fn test_from_bytes() {
        let mut rng = OsRng;
        for _ in 0..100 {
            let a = Fp2::random(&mut rng);
            let bytes = a.to_bytes();
            assert_eq!(a, Fp2::from_bytes(&bytes).unwrap());
        }

        assert_eq!(
            Fp2::from_bytes(&[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]).unwrap(),
            Fp2::zero()
        );

        assert_eq!(
            Fp2::from_bytes(&[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]).unwrap(),
            Fp2::one()
        );

        // -1 should work
        assert_eq!(
            Fp2::from_bytes(&[0, 0, 0, 0, 0, 0, 128, 65, 0, 0, 0, 0, 0, 0, 0, 0]).unwrap(),
            -Fp2::one()
        );

        // Anything larger than M in one of the members is invalid
        assert!(bool::from(
            Fp2::from_bytes(&[1, 0, 0, 0, 0, 0, 128, 65, 0, 0, 0, 0, 0, 0, 0, 0]).is_none()
        ));

        assert!(bool::from(
            Fp2::from_bytes(&[0, 0, 0, 0, 255, 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0]).is_none()
        ));

        assert!(bool::from(
            Fp2::from_bytes(&[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255]).is_none()
        ));
    }

    #[test]
    #[cfg(feature = "serialize")]
    fn test_serde() {
        let mut rng = OsRng;
        let element = Fp2::random(&mut rng);
        let encoded = bincode::serialize(&element).unwrap();
        let parsed: Fp2 = bincode::deserialize(&encoded).unwrap();
        assert_eq!(parsed, element);

        // Check that the encoding is 16 bytes exactly
        assert_eq!(encoded.len(), 16);

        // Check that the encoding itself matches the usual one
        assert_eq!(element, bincode::deserialize(&element.to_bytes()).unwrap());

        // Check that invalid encodings fail
        let element = Fp2::random(&mut rng);
        let mut encoded = bincode::serialize(&element).unwrap();
        encoded[15] = 127;
        assert!(bincode::deserialize::<Fp2>(&encoded).is_err());

        let encoded = bincode::serialize(&element).unwrap();
        assert!(bincode::deserialize::<Fp2>(&encoded[0..15]).is_err());
    }
}
