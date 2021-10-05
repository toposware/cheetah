//! This module implements arithmetic over the extension field Fp2,
//! defined with irreducible polynomial u^2 - x - 1.

use core::fmt;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use rand_core::{CryptoRng, RngCore};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

use crate::fp::Fp;

const MODULUS_MINUS_ONE_DIV_TWO: [u64; 2] = [0x3fffc88000000000, 0x07fff22006042000];

const MODULUS_PLUS_ONE_DIV_TWO: [u64; 2] = [0x3fffc88000000001, 0x07fff22006042000];

#[derive(Copy, Clone)]
/// An element of the extension GF(p^2)
pub struct Fp2 {
    /// First coefficient
    pub c0: Fp,
    /// Second coefficient
    pub c1: Fp,
}

impl fmt::Debug for Fp2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
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
    pub const fn new(value: u64) -> Self {
        Fp2 {
            c0: Fp::new(value),
            c1: Fp::zero(),
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

    /// Generates a random element
    pub fn random(mut rng: impl CryptoRng + RngCore) -> Self {
        Fp2 {
            c0: Fp::random(&mut rng),
            c1: Fp::random(&mut rng),
        }
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
    // TODO: check how to reduce number of multiplications
    pub const fn mul(&self, other: &Fp2) -> Fp2 {
        let aa = (&self.c0).mul(&other.c0);
        let ab = (&self.c0).mul(&other.c1);

        let ba = (&self.c1).mul(&other.c0);
        let bb = (&self.c1).mul(&other.c1);

        let c0 = (&bb).add(&aa);

        let c1 = (&bb).add(&ab);
        let c1 = (&c1).add(&ba);

        Fp2 { c0, c1 }
    }

    /// Square this element
    #[inline]
    // TODO: check how to reduce number of multiplications
    pub const fn square(&self) -> Self {
        let aa = &self.c0.square();
        let ab = (&self.c0).mul(&self.c1);

        let bb = &self.c1.square();

        let c0 = aa.add(bb);

        let c1 = bb.add(&ab);
        let c1 = (&c1).add(&ab);

        Fp2 { c0, c1 }
    }

    /// Double this element
    pub const fn double(&self) -> Self {
        self.add(self)
    }

    /// Computes the square root of this element, if it exists.
    /// This operation is not constant time, as it requires going
    /// through a non-constant number of field elements until we reach
    /// a non-square.
    // TODO: Should rng be passed as argument?
    pub fn sqrt_vartime(&self) -> Option<Self> {
        // Cipolla's algorithm
        // https://en.wikipedia.org/wiki/Cipolla%27s_algorithm

        fn mul_in_fp4(a: &[Fp2; 2], b: &[Fp2; 2], beta: &Fp2) -> [Fp2; 2] {
            let v0 = (&a[0]).mul(&b[0]);
            let v1 = (&a[1]).mul(&b[1]);
            let c0 = (&v0).add(&(beta.mul(&v1)));
            let c1 = (&a[0]).add(&a[1]);
            let c1 = (&c1).mul(&b[0].add(&b[1]));
            let c1 = (&c1).sub(&v0);
            let c1 = (&c1).sub(&v1);

            [c0, c1]
        }

        fn exp_vartime_in_fp4(a: &[Fp2; 2], by: &[u64; 2], beta: &Fp2) -> [Fp2; 2] {
            let mut res = [Fp2::one(), Fp2::zero()];
            for e in by.iter().rev() {
                for i in (0..64).rev() {
                    res = mul_in_fp4(&res, &res, beta);

                    if ((*e >> i) & 1) == 1 {
                        res = mul_in_fp4(&res, a, beta);
                    }
                }
            }
            res
        }

        // check = self^((p^2 - 1) // 2)
        //       = self^0x7fff220060420003fffc88000000000
        let check = self.exp_vartime(&MODULUS_MINUS_ONE_DIV_TWO);

        if check != Fp2::one() {
            return None;
        }

        let mut a = Fp2::from(Fp::one().double());
        while (a.square() - self).exp_vartime(&MODULUS_MINUS_ONE_DIV_TWO) != -Fp2::one() {
            a += Fp2::one();
        }

        let beta = a.square() - self;

        let res = exp_vartime_in_fp4(&[a, Fp2::one()], &MODULUS_PLUS_ONE_DIV_TWO, &beta);
        if res[1] != Fp2::zero() {
            return None;
        }

        if bool::from(self.ct_eq(&res[0].square())) {
            Some(res[0])
        } else {
            None
        }
    }

    /// Add two elements together
    #[inline]
    pub const fn add(&self, rhs: &Self) -> Self {
        Fp2 {
            c0: (&self.c0).add(&rhs.c0),
            c1: (&self.c1).add(&rhs.c1),
        }
    }

    /// Substract two elements together
    #[inline]
    pub const fn sub(&self, rhs: &Self) -> Self {
        Fp2 {
            c0: (&self.c0).sub(&rhs.c0),
            c1: (&self.c1).sub(&rhs.c1),
        }
    }

    /// Negate `&self`
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
        let inv = self.c0.square() + self.c0 * self.c1 - self.c1.square();

        inv.invert().map(|t| Fp2 {
            c0: (self.c0 + self.c1) * t,
            c1: -self.c1 * t,
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
    pub const fn to_repr(&self) -> [u64; 2] {
        [
            Fp::montgomery_reduce(self.c0.0, 0).0,
            Fp::montgomery_reduce(self.c1.0, 0).0,
        ]
    }

    #[inline(always)]
    /// Normalizes the internal representation of a `Fp2` element
    pub fn normalize(&mut self) {
        self.c0.normalize();
        self.c1.normalize();
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::thread_rng;

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
    }

    #[test]
    fn test_squaring() {
        let a = Fp2 {
            c0: Fp::new(1373422459643663482),
            c1: Fp::new(1559760838287424487),
        };
        let b = Fp2 {
            c0: Fp::new(4508439400065818468),
            c1: Fp::new(3727974599246273355),
        };

        assert_eq!(a.square(), b);
    }

    #[test]
    fn test_sqrt() {
        for _ in 0..10 {
            let a = Fp2::random(&mut thread_rng()).square();
            let b = a.sqrt_vartime().unwrap();
            assert_eq!(a, b.square());
        }
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
            c0: Fp::new(1373422459643663482),
            c1: Fp::new(1559760838287424487),
        };
        let b = Fp2 {
            c0: Fp::new(4508439400065818468),
            c1: Fp::new(3727974599246273355),
        };
        let c = Fp2 {
            c0: Fp::new(3944078691460293030),
            c1: Fp::new(1622844793164244792),
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
        let a = Fp2 {
            c0: Fp::one(),
            c1: Fp::new(2),
        };
        let b = Fp2 {
            c0: -Fp::one(),
            c1: -Fp::new(2),
        };

        assert_eq!(-a, b);
    }

    #[test]
    fn test_inversion() {
        let a = Fp2 {
            c0: Fp::new(2077468337887652729),
            c1: Fp::new(3560716367442326177),
        };

        let b = Fp2 {
            c0: Fp::new(293161233865469014),
            c1: Fp::new(2017440449651245094),
        };

        assert_eq!(a.invert().unwrap(), b);

        assert!(bool::from(Fp2::zero().invert().is_none()));
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
}
