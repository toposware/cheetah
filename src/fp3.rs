// Copyright (c) 2021-2023 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This module implements arithmetic over the extension field Fp3.
//! The implementation is minimal and not exposed through the public
//! API. It is only for internal use in Fp6 arithmetic operations.

use core::fmt::{self, Formatter};

use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

use crate::fp::reduce_u96;
use crate::fp::Fp;
use crate::fp6::Fp6;

use crate::fp::TWO_ADICITY;

const BETA: u128 = crate::fp::GENERATOR.0 as u128;

// 2^32 root of unity = 2800184025912956819
const TWO_ADIC_ROOT_OF_UNITY_P3: Fp3 = Fp3 {
    a0: Fp(2800184025912956819),
    a1: Fp::zero(),
    a2: Fp::zero(),
};

#[derive(Copy, Clone)]
/// An element of the extension GF(p^3)
pub(crate) struct Fp3 {
    /// First coefficient, lowest degree
    pub(crate) a0: Fp,
    /// Second coefficient
    pub(crate) a1: Fp,
    /// Third coefficient, highest degree
    pub(crate) a2: Fp,
}

impl fmt::Debug for Fp3 {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        let coeffs = [self.a0, self.a1, self.a2];

        let format_term = |coef: Fp, degree: usize| -> String {
            if coef == Fp::one() && degree > 0 {
                let exp = if degree > 1 {
                    format!("^{}", degree)
                } else {
                    "".to_string()
                };
                format!("u{}", exp)
            } else {
                let exp =  if degree > 0 {
                    format!("u{}", if degree > 1 { format!("^{}", degree) } else { "".to_string() })
                } else {
                    "".to_string()
                };
                format!("{}{}", coef, exp)
            }
        };

        let elem_rep = coeffs
            .iter()
            .enumerate()
            .filter_map(|(i, &coef)| {
                if coef == Fp::zero() {
                    None
                } else {
                    Some(format_term(coef, i))
                }
            })
            .collect::<Vec<String>>()
            .join(" + ");

        if *self == Fp3::zero() {
            return write!(f, "0"); // Handle the case where all coefficients are zero
        } 
        
        write!(f, "{}", elem_rep)?;

        Ok(())
    }
}

impl Default for Fp3 {
    fn default() -> Self {
        Self::zero()
    }
}

// When looking at the sextic extension as a towered one,
// i.e. Fp6 = Fp3[Y]/(Y^2 − γ) with γ = δ = 7, we use the
// lowest coefficient of the quadratic extension for
// conversion from Fp6 to Fp3 elements.
impl From<&Fp6> for Fp3 {
    fn from(f: &Fp6) -> Self {
        Self {
            a0: f.c0,
            a1: f.c2,
            a2: f.c4,
        }
    }
}

impl ConstantTimeEq for Fp3 {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.a0.ct_eq(&other.a0) & self.a1.ct_eq(&other.a1) & self.a2.ct_eq(&other.a2)
    }
}

impl Eq for Fp3 {}
impl PartialEq for Fp3 {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl ConditionallySelectable for Fp3 {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self {
            a0: Fp::conditional_select(&a.a0, &b.a0, choice),
            a1: Fp::conditional_select(&a.a1, &b.a1, choice),
            a2: Fp::conditional_select(&a.a2, &b.a2, choice),
        }
    }
}

impl Fp3 {
    #[inline]
    /// The additive identity
    pub(crate) const fn zero() -> Self {
        Self {
            a0: Fp::zero(),
            a1: Fp::zero(),
            a2: Fp::zero(),
        }
    }

    #[inline]
    /// The multiplicative identity
    pub(crate) const fn one() -> Self {
        Self {
            a0: Fp::one(),
            a1: Fp::zero(),
            a2: Fp::zero(),
        }
    }

    #[inline]
    /// Computes the multiplication of two Fp3 elements
    pub(crate) const fn mul(&self, other: &Fp3) -> Fp3 {
        let t00 = (&self.a0).mul(&other.a0).0 as u128;
        let t01 = (&self.a1).mul(&other.a1).0 as u128;
        let t02 = (&self.a2).mul(&other.a2).0 as u128;

        let s012 = (&self.a1).add(&self.a2);
        let tmp = (&other.a1).add(&other.a2);
        let s012 = (&s012).mul(&tmp).0 as u128;

        let s001 = (&self.a0).add(&self.a1);
        let tmp = (&other.a0).add(&other.a1);
        let s001 = (&s001).mul(&tmp).0 as u128;

        let s002 = (&self.a0).add(&self.a2);
        let tmp = (&other.a0).add(&other.a2);
        let s002 = (&s002).mul(&tmp).0 as u128;

        let d00 = t01 + t02;
        let d00 = s012 + 0x7fffffff800000008 - d00;
        let d00 = d00 * BETA;
        let d00 = d00 + t00;

        let d01 = t02 * BETA;
        let tmp = t00 + t01;
        let d01 = d01 + 0x7fffffff800000008 - tmp;
        let d01 = d01 + s001;

        let d02 = t00 + t02;
        let d02 = t01 + 0x7fffffff800000008 - d02;
        let d02 = d02 + s002;

        // Compute the final coordinates, reduced by the modulus
        let a0 = Fp(reduce_u96(d00));
        let a1 = Fp(reduce_u96(d01));
        let a2 = Fp(reduce_u96(d02));

        Fp3 { a0, a1, a2 }
    }

    /// Computes the square of a field element
    #[inline]
    pub(crate) const fn square(&self) -> Self {
        let t00 = (&self.a0).square().0 as u128;
        let t01 = (&self.a1).square().0 as u128;
        let t02 = (&self.a2).square().0 as u128;

        let s012 = (&self.a1).add(&self.a2);
        let s012 = (&s012).square().0 as u128;

        let s001 = (&self.a0).add(&self.a1);
        let s001 = (&s001).square().0 as u128;

        let s002 = (&self.a0).add(&self.a2);
        let s002 = (&s002).square().0 as u128;

        let d00 = t01 + t02;
        let d00 = s012 + 0x7fffffff800000008 - d00;
        let d00 = d00 * BETA;
        let d00 = d00 + t00;

        let d01 = t02 * BETA;
        let tmp = t00 + t01;
        let d01 = d01 + 0x7fffffff800000008 - tmp;
        let d01 = d01 + s001;

        let d02 = t00 + t02;
        let d02 = t01 + 0x7fffffff800000008 - d02;
        let d02 = d02 + s002;

        // Compute the final coordinates, reduced by the modulus
        let a0 = Fp(reduce_u96(d00));
        let a1 = Fp(reduce_u96(d01));
        let a2 = Fp(reduce_u96(d02));

        Fp3 { a0, a1, a2 }
    }

    /// Computes the square root of this element, if it exists.
    pub(crate) fn sqrt(&self) -> CtOption<Self> {
        // Tonelli-Shank's algorithm for q mod 16 = 1
        // See https://eprint.iacr.org/2020/1497.pdf, page 3 for a
        // constant time specification of the algorithm.

        // Compute the progenitor y of self
        // y = self^((t - 1) // 2)
        //   = self^0x7ffffffe80000002fffffffc80000002fffffffe
        let y = self.exp_vartime(&[0x80000002fffffffe, 0x80000002fffffffc, 0x000000007ffffffe]);

        let mut s = self.mul(&y);
        let mut t = s.mul(&y);

        let mut z = TWO_ADIC_ROOT_OF_UNITY_P3;

        for k in (2..=TWO_ADICITY).rev() {
            let mut b = t;

            for _ in 0..k - 2 {
                b = b.square();
            }

            let new_s = s.mul(&z);
            s = Fp3::conditional_select(&new_s, &s, b.ct_eq(&Fp3::one()));
            z = z.square();
            let new_t = t.mul(&z);
            t = Fp3::conditional_select(&new_t, &t, b.ct_eq(&Fp3::one()));
        }

        CtOption::new(s, (s.square()).ct_eq(self))
    }

    /// Exponentiates `self` by `power`, where `power` is a
    /// little-endian order integer exponent.
    ///
    /// **This operation is variable time with respect
    /// to the exponent.** If the exponent is fixed,
    /// this operation is effectively constant time.
    pub(crate) fn exp_vartime(&self, power: &[u64]) -> Self {
        let mut res = Self::one();
        for e in power.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();

                if ((*e >> i) & 1) == 1 {
                    res = res.mul(self);
                }
            }
        }
        res
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::{OsRng, RngCore};

    impl Fp3 {
        /// Generates a random canonical element
        fn random(mut rng: impl RngCore) -> Self {
            Self {
                a0: Fp::random(&mut rng),
                a1: Fp::random(&mut rng),
                a2: Fp::random(&mut rng),
            }
        }
    }

    // DISPLAY
    // ================================================================================================
    #[test]
    fn test_debug() {
        assert_eq!(format!("{:?}", Fp3::zero()), "0");
        assert_eq!(format!("{:?}", Fp3::one()), "1");

        let a = Fp3 {
            a0: Fp::new(0),
            a1: Fp::new(0),
            a2: Fp::new(7)
        };
        assert_eq!(format!("{:?}", a), "7u^2");

        let b = Fp3 {
            a0: Fp::new(1),
            a1: Fp::new(0),
            a2: Fp::new(11)
        };
        assert_eq!(format!("{:?}", b), "1 + 11u^2");

        let c = Fp3 {
            a0: Fp::new(1),
            a1: Fp::new(2),
            a2: Fp::new(1)
        };
        assert_eq!(format!("{:?}", c), "1 + 2u + u^2");
    }
    // BASIC ALGEBRA
    // ================================================================================================

    #[test]
    fn test_conditional_selection() {
        let a = Fp3 {
            a0: Fp::one(),
            a1: Fp::new(3),
            a2: Fp::new(5),
        };
        let b = Fp3 {
            a0: Fp::new(7),
            a1: Fp::new(9),
            a2: Fp::new(11),
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
        fn is_equal(a: &Fp3, b: &Fp3) -> bool {
            let eq = a == b;
            let ct_eq = a.ct_eq(b);

            assert_eq!(eq, bool::from(ct_eq));

            eq
        }

        assert!(is_equal(
            &Fp3 {
                a0: Fp::one(),
                a1: Fp::new(3),
                a2: Fp::new(5),
            },
            &Fp3 {
                a0: Fp::one(),
                a1: Fp::new(3),
                a2: Fp::new(5),
            }
        ));

        assert!(!is_equal(
            &Fp3 {
                a0: Fp::new(2),
                a1: Fp::new(3),
                a2: Fp::new(5),
            },
            &Fp3 {
                a0: Fp::one(),
                a1: Fp::new(3),
                a2: Fp::new(5),
            }
        ));

        assert!(!is_equal(
            &Fp3 {
                a0: Fp::one(),
                a1: Fp::new(2),
                a2: Fp::new(5),
            },
            &Fp3 {
                a0: Fp::one(),
                a1: Fp::new(3),
                a2: Fp::new(5),
            }
        ));

        assert!(!is_equal(
            &Fp3 {
                a0: Fp::one(),
                a1: Fp::new(3),
                a2: Fp::new(4),
            },
            &Fp3 {
                a0: Fp::one(),
                a1: Fp::new(3),
                a2: Fp::new(5),
            }
        ));

        assert!(!bool::from(Fp3::zero().ct_eq(&Fp3::one())));
    }

    #[test]
    fn test_squaring() {
        let a = Fp3 {
            a0: Fp::new(17095662164057750615),
            a1: Fp::new(10249841381153164628),
            a2: Fp::new(14023641212482863105),
        };
        let b = Fp3 {
            a0: Fp::new(11557468600072973281),
            a1: Fp::new(7488412806431150554),
            a2: Fp::new(14594457941668144169),
        };

        assert_eq!(a.square(), b);
    }

    #[test]
    fn test_sqrt() {
        for _ in 0..100 {
            let a = Fp3::random(&mut OsRng).square();
            let b = a.sqrt().unwrap();
            assert_eq!(a, b.square());
        }

        assert_eq!(Fp3::zero().sqrt().unwrap(), Fp3::zero());
        assert_eq!(Fp3::one().sqrt().unwrap(), Fp3::one());

        // u + 5
        // is not a quadratic residue in Fp3
        assert!(bool::from(
            Fp3 {
                a0: Fp::new(5),
                a1: Fp::one(),
                a2: Fp::zero(),
            }
            .sqrt()
            .is_none()
        ));
    }

    #[test]
    fn test_multiplication() {
        let a = Fp3 {
            a0: Fp::one(),
            a1: Fp::new(3),
            a2: Fp::new(5),
        };
        let b = Fp3::one();
        let c = Fp3 {
            a0: Fp::one(),
            a1: Fp::new(3),
            a2: Fp::new(5),
        };

        assert_eq!(a.mul(&b), c);

        let a = Fp3 {
            a0: Fp::new(4063754064363348009),
            a1: Fp::new(5263739465726006019),
            a2: Fp::new(6460939331453405671),
        };
        let b = Fp3 {
            a0: Fp::new(16715421830871665632),
            a1: Fp::new(11195742925840789795),
            a2: Fp::new(9795579848194115388),
        };
        let c = Fp3 {
            a0: Fp::new(3924021414384413108),
            a1: Fp::new(3448800052342118486),
            a2: Fp::new(14917403927156629509),
        };

        assert_eq!(a.mul(&b), c);
    }

    // ROOTS OF UNITY
    // ================================================================================================

    #[test]
    fn test_get_root_of_unity() {
        let two_pow_32 = 1 << TWO_ADICITY as u64;
        assert_eq!(
            Fp3::one(),
            TWO_ADIC_ROOT_OF_UNITY_P3.exp_vartime(&[two_pow_32, 0, 0,])
        );
        assert_ne!(
            Fp3::one(),
            TWO_ADIC_ROOT_OF_UNITY_P3.exp_vartime(&[two_pow_32 - 1, 0, 0,])
        );
    }
}
