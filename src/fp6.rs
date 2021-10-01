//! This module implements arithmetic over the extension field Fp6,
//! defined with irreducible polynomial v^3 - 4v - 1.

use core::convert::TryInto;
use core::fmt;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use rand_core::{CryptoRng, RngCore};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

use crate::fp::Fp;
use crate::fp2::Fp2;

const MODULUS_MINUS_ONE_DIV_TWO: [u64; 6] = [
    0xbfff598000000000,
    0x953f2fe05a3de000,
    0xc22d66eed23dc5ea,
    0xb9538f1d34e3dc07,
    0xa978d997e2ea6efc,
    0x0007ffd6605a3d77,
];

const MODULUS_PLUS_ONE_DIV_TWO: [u64; 6] = [
    0xbfff598000000001,
    0x953f2fe05a3de000,
    0xc22d66eed23dc5ea,
    0xb9538f1d34e3dc07,
    0xa978d997e2ea6efc,
    0x0007ffd6605a3d77,
];

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
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
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

impl From<Fp> for Fp6 {
    fn from(f: Fp) -> Self {
        Fp6::from(Fp2::from(f))
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

    /// Generates a random element
    pub fn random(mut rng: impl CryptoRng + RngCore) -> Self {
        Fp6 {
            c0: Fp2::random(&mut rng),
            c1: Fp2::random(&mut rng),
            c2: Fp2::random(&mut rng),
        }
    }

    #[inline(always)]
    /// Multiply by nonresidue u.
    pub const fn mul_by_nonresidue(&self) -> Self {
        // Given a + bv + cv^2, this produces
        //     av + bv^2 + cv^3
        // but because v^3 = 4v + 1, we have
        //     c + (a + 4c)v + bv^2

        Fp6 {
            c0: self.c2,
            c1: (&(&self.c0).add(&Fp2 {
                c0: Fp::new(4),
                c1: Fp::zero(),
            }))
                .mul(&self.c2),
            c2: self.c1,
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

        self.c2.lexicographically_largest()
            | (self.c2.is_zero() & self.c1.lexicographically_largest())
            | (self.c1.is_zero() & self.c0.lexicographically_largest())
    }

    #[inline]
    /// Computes the multiplication of two Fp6 elements
    pub const fn mul(&self, other: &Fp6) -> Fp6 {
        let aa = (&self.c0).mul(&other.c0);
        let ab = (&self.c0).mul(&other.c1);
        let ac = (&self.c0).mul(&other.c2);

        let ba = (&self.c1).mul(&other.c0);
        let bb = (&self.c1).mul(&other.c1);
        let bc = (&self.c1).mul(&other.c2);

        let ca = (&self.c2).mul(&other.c0);
        let cb = (&self.c2).mul(&other.c1);
        let cc = (&self.c2).mul(&other.c2);

        let c0 = (&aa).add(&bc);
        let c0 = (&c0).add(&cb);

        let c1 = (&bc).add(&cb);
        let c1 = (&c1).mul(&Fp2 {
            c0: Fp::new(4),
            c1: Fp::zero(),
        });
        let c1 = (&c1).add(&ab);
        let c1 = (&c1).add(&ba);
        let c1 = (&c1).add(&cc);

        let c2 = (&cc).mul(&Fp2 {
            c0: Fp::new(4),
            c1: Fp::zero(),
        });
        let c2 = (&c2).add(&ac);
        let c2 = (&c2).add(&ca);
        let c2 = (&c2).add(&bb);

        Fp6 { c0, c1, c2 }
    }

    /// Square this element
    #[inline]
    pub const fn square(&self) -> Self {
        let aa = (&self.c0).mul(&self.c0);
        let ab = (&self.c0).mul(&self.c1);
        let ac = (&self.c0).mul(&self.c2);

        let bb = (&self.c1).mul(&self.c1);
        let bc = (&self.c1).mul(&self.c2);

        let cc = (&self.c2).mul(&self.c2);

        let c0 = (&aa).add(&bc);
        let c0 = (&c0).add(&bc);

        let c1 = (&bc).add(&bc);
        let c1 = (&c1).mul(&Fp2 {
            c0: Fp::new(4),
            c1: Fp::zero(),
        });
        let c1 = (&c1).add(&ab);
        let c1 = (&c1).add(&ab);
        let c1 = (&c1).add(&cc);

        let c2 = (&cc).mul(&Fp2 {
            c0: Fp::new(4),
            c1: Fp::zero(),
        });
        let c2 = (&c2).add(&ac);
        let c2 = (&c2).add(&ac);
        let c2 = (&c2).add(&bb);

        Fp6 { c0, c1, c2 }
    }

    /// Computes the square root of this element, if it exists.
    /// This operation is not constant time, as it requires going
    /// through a non-constant number of field elements until we reach
    /// a non-square.
    // TODO: Should rng be passed as argument?
    pub fn sqrt_vartime(&self) -> Option<Self> {
        // Cipolla's algorithm
        // https://en.wikipedia.org/wiki/Cipolla%27s_algorithm

        fn mul_in_fp12(a: &[Fp6; 2], b: &[Fp6; 2], beta: &Fp6) -> [Fp6; 2] {
            let v0 = (&a[0]).mul(&b[0]);
            let v1 = (&a[1]).mul(&b[1]);
            let c0 = (&v0).add(&(beta.mul(&v1)));
            let c1 = (&a[0]).add(&a[1]);
            let c1 = (&c1).mul(&b[0].add(&b[1]));
            let c1 = (&c1).sub(&v0);
            let c1 = (&c1).sub(&v1);

            [c0, c1]
        }

        fn exp_vartime_in_fp12(a: &[Fp6; 2], by: &[u64; 6], beta: &Fp6) -> [Fp6; 2] {
            let mut res = [Fp6::one(), Fp6::zero()];
            for e in by.iter().rev() {
                for i in (0..64).rev() {
                    res = mul_in_fp12(&res, &res, beta);

                    if ((*e >> i) & 1) == 1 {
                        res = mul_in_fp12(&res, a, beta);
                    }
                }
            }
            res
        }

        // check = self^((p^6 - 1) // 2)
        //       = self^0x7ffd6605a3d77a978d997e2ea6efcb9538f1d34e3dc07c22d66eed23dc5ea953f2fe05a3de000bfff598000000000
        let check = self.exp_vartime(&MODULUS_MINUS_ONE_DIV_TWO);

        if check != Fp6::one() {
            return None;
        }

        let mut a = Fp6::from(Fp2::one().double());
        while (a.square() - self).exp_vartime(&MODULUS_MINUS_ONE_DIV_TWO) != -Fp6::one() {
            a += Fp6::one();
        }

        let beta = a.square() - self;

        let res = exp_vartime_in_fp12(&[a, Fp6::one()], &MODULUS_PLUS_ONE_DIV_TWO, &beta);
        if res[1] != Fp6::zero() {
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
        Fp6 {
            c0: (&self.c0).add(&rhs.c0),
            c1: (&self.c1).add(&rhs.c1),
            c2: (&self.c2).add(&rhs.c2),
        }
    }

    /// Substract two elements together
    #[inline]
    pub const fn sub(&self, rhs: &Self) -> Self {
        Fp6 {
            c0: (&self.c0).sub(&rhs.c0),
            c1: (&self.c1).sub(&rhs.c1),
            c2: (&self.c2).sub(&rhs.c2),
        }
    }

    /// Negate `&self`
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
        let a2 = (&self.c0).square();
        let b2 = (&self.c1).square();
        let c2 = (&self.c2).square();

        let two = Fp2::one() + Fp2::one();
        let three = two + Fp2::one();
        let four = three + Fp2::one();
        let eight = four.double();
        let sixteen = four * four;

        let inv = a2 * self.c0 - four * self.c0 * b2 + self.c1 * b2 + sixteen * self.c0 * c2
            - four * self.c1 * c2
            + self.c2 * c2
            + (eight * a2 - three * self.c0 * self.c1) * self.c2;

        let c0 = a2 - four * b2 + (eight * self.c0 - self.c1) * self.c2 + sixteen * c2;
        let c1 = -self.c0 * self.c1 + c2;
        let c2 = b2 - self.c0 * self.c2 - four * c2;

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

    /// Outputs the internal representation as 4 64-bit limbs after Montgomery reduction
    pub const fn to_repr(&self) -> [u64; 6] {
        [
            Fp::montgomery_reduce(self.c0.c0.0, 0).0,
            Fp::montgomery_reduce(self.c0.c1.0, 0).0,
            Fp::montgomery_reduce(self.c1.c0.0, 0).0,
            Fp::montgomery_reduce(self.c1.c1.0, 0).0,
            Fp::montgomery_reduce(self.c2.c0.0, 0).0,
            Fp::montgomery_reduce(self.c2.c1.0, 0).0,
        ]
    }

    /// Converts a `Fp` element into a byte representation in
    /// little-endian byte order.
    pub fn to_bytes(&self) -> [u8; 48] {
        // Turn into canonical form by computing
        // (a.R) / R = a
        let tmp = self.to_repr();

        let mut res = [0u8; 48];
        res[0..8].copy_from_slice(&tmp[0].to_le_bytes());
        res[8..16].copy_from_slice(&tmp[1].to_le_bytes());
        res[16..24].copy_from_slice(&tmp[2].to_le_bytes());
        res[24..32].copy_from_slice(&tmp[3].to_le_bytes());
        res[32..40].copy_from_slice(&tmp[4].to_le_bytes());
        res[40..48].copy_from_slice(&tmp[5].to_le_bytes());

        res
    }

    /// Converts an array of bytes into an `Fp6` element
    pub fn from_bytes(bytes: &[u8; 48]) -> Fp6 {
        let mut res = Fp6::zero();

        let mut array: [u8; 16] = bytes[0..16].try_into().unwrap();
        res.c0 = Fp2::from_bytes(&array);
        array = bytes[16..32].try_into().unwrap();
        res.c1 = Fp2::from_bytes(&array);
        array = bytes[32..48].try_into().unwrap();
        res.c2 = Fp2::from_bytes(&array);

        res
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::thread_rng;

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
    }

    #[test]
    fn test_squaring() {
        let a = Fp6 {
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
        };
        let b = Fp6 {
            c0: Fp2 {
                c0: Fp::new(411871884247580752),
                c1: Fp::new(975422536930101490),
            },
            c1: Fp2 {
                c0: Fp::new(3660111879911998944),
                c1: Fp::new(3658765734269830397),
            },
            c2: Fp2 {
                c0: Fp::new(1336378171450403537),
                c1: Fp::new(3568041835212862283),
            },
        };

        assert_eq!(a.square(), b);
    }

    #[test]
    fn test_sqrt() {
        for _ in 0..10 {
            let a = Fp6::random(&mut thread_rng()).square();
            let b = a.sqrt_vartime().unwrap();
            assert_eq!(a, b.square());
        }
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
                c0: Fp::new(6),
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
                c0: Fp::new(6),
                c1: Fp::new(6),
            },
        };

        assert_eq!(a * b, c);

        let a = Fp6 {
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
        };
        let b = Fp6 {
            c0: Fp2 {
                c0: Fp::new(411871884247580752),
                c1: Fp::new(975422536930101490),
            },
            c1: Fp2 {
                c0: Fp::new(3660111879911998944),
                c1: Fp::new(3658765734269830397),
            },
            c2: Fp2 {
                c0: Fp::new(1336378171450403537),
                c1: Fp::new(3568041835212862283),
            },
        };
        let c = Fp6 {
            c0: Fp2 {
                c0: Fp::new(4100979495961732988),
                c1: Fp::new(815786812220801534),
            },
            c1: Fp2 {
                c0: Fp::new(4558525380807129266),
                c1: Fp::new(3579606074809773744),
            },
            c2: Fp2 {
                c0: Fp::new(3450239709028392840),
                c1: Fp::new(4151094398153102512),
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
    }

    #[test]
    fn test_inversion() {
        let a = Fp6 {
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
        };

        let b = Fp6 {
            c0: Fp2 {
                c0: Fp::new(250182847153364610),
                c1: Fp::new(1902439816108175911),
            },
            c1: Fp2 {
                c0: Fp::new(218641481061462196),
                c1: Fp::new(3629663639262592947),
            },
            c2: Fp2 {
                c0: Fp::new(376185455706256631),
                c1: Fp::new(854203883548072415),
            },
        };

        assert_eq!(a.invert().unwrap(), b);

        assert_eq!(Fp6::one().invert().unwrap(), Fp6::one());

        assert!(bool::from(Fp6::zero().invert().is_none()));
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
    fn test_bytes() {
        let mut rng = thread_rng();
        for _ in 0..100 {
            let a = Fp6::random(&mut rng);
            let bytes = a.to_bytes();
            assert_eq!(a, Fp6::from_bytes(&bytes));
        }
    }

    #[test]
    fn test_zeroize() {
        use zeroize::Zeroize;

        let mut a = Fp6::one();
        a.zeroize();
        assert!(bool::from(a.is_zero()));
    }
}
