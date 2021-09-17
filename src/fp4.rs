//! This module implements arithmetic over the extension field Fp4,
//! defined with irreducible polynomial u^4 + 3.

use core::fmt;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use rand_core::{CryptoRng, RngCore};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

use crate::fp::Fp;

/// U4 = -3
const U4: Fp = Fp::new(4611624995532046334);

#[derive(Copy, Clone)]
/// An element of the extension GF(p^4)
pub struct Fp4 {
    /// First coefficient
    pub c0: Fp,
    /// Second coefficient
    pub c1: Fp,
    /// Third coefficient
    pub c2: Fp,
    /// Fourth coefficient
    pub c3: Fp,
}

impl fmt::Debug for Fp4 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{:?} + {:?}*u + {:?}*u^2 + {:?}*u^3",
            self.c0, self.c1, self.c2, self.c3
        )
    }
}

impl Default for Fp4 {
    fn default() -> Self {
        Fp4::zero()
    }
}

impl zeroize::DefaultIsZeroes for Fp4 {}

impl From<Fp> for Fp4 {
    fn from(f: Fp) -> Self {
        Fp4 {
            c0: f,
            c1: Fp::zero(),
            c2: Fp::zero(),
            c3: Fp::zero(),
        }
    }
}

impl ConstantTimeEq for Fp4 {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.c0.ct_eq(&other.c0)
            & self.c1.ct_eq(&other.c1)
            & self.c2.ct_eq(&other.c2)
            & self.c3.ct_eq(&other.c3)
    }
}

impl Eq for Fp4 {}
impl PartialEq for Fp4 {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl ConditionallySelectable for Fp4 {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Fp4 {
            c0: Fp::conditional_select(&a.c0, &b.c0, choice),
            c1: Fp::conditional_select(&a.c1, &b.c1, choice),
            c2: Fp::conditional_select(&a.c2, &b.c2, choice),
            c3: Fp::conditional_select(&a.c3, &b.c3, choice),
        }
    }
}

impl<'a> Neg for &'a Fp4 {
    type Output = Fp4;

    #[inline]
    fn neg(self) -> Fp4 {
        self.neg()
    }
}

impl Neg for Fp4 {
    type Output = Fp4;

    #[inline]
    fn neg(self) -> Fp4 {
        -&self
    }
}

impl<'a, 'b> Sub<&'b Fp4> for &'a Fp4 {
    type Output = Fp4;

    #[inline]
    fn sub(self, rhs: &'b Fp4) -> Fp4 {
        self.sub(rhs)
    }
}

impl<'a, 'b> Add<&'b Fp4> for &'a Fp4 {
    type Output = Fp4;

    #[inline]
    fn add(self, rhs: &'b Fp4) -> Fp4 {
        self.add(rhs)
    }
}

impl<'a, 'b> Mul<&'b Fp4> for &'a Fp4 {
    type Output = Fp4;

    #[inline]
    fn mul(self, rhs: &'b Fp4) -> Fp4 {
        self.mul(rhs)
    }
}

impl_binops_additive!(Fp4, Fp4);
impl_binops_multiplicative!(Fp4, Fp4);

impl Fp4 {
    #[inline]
    /// The additive identity
    pub const fn zero() -> Self {
        Fp4 {
            c0: Fp::zero(),
            c1: Fp::zero(),
            c2: Fp::zero(),
            c3: Fp::zero(),
        }
    }

    #[inline]
    /// The multiplicative identity
    pub const fn one() -> Self {
        Fp4 {
            c0: Fp::one(),
            c1: Fp::zero(),
            c2: Fp::zero(),
            c3: Fp::zero(),
        }
    }

    /// Checks whether this element is zero or not
    pub fn is_zero(&self) -> Choice {
        self.c0.is_zero() & self.c1.is_zero() & self.c2.is_zero() & self.c3.is_zero()
    }

    /// Generates a random element
    pub fn random(mut rng: impl CryptoRng + RngCore) -> Self {
        Fp4 {
            c0: Fp::random(&mut rng),
            c1: Fp::random(&mut rng),
            c2: Fp::random(&mut rng),
            c3: Fp::random(&mut rng),
        }
    }

    #[inline(always)]
    /// Multiply by nonresidue u.
    pub const fn mul_by_nonresidue(&self) -> Self {
        // Given a + bu + cu^2 + du^3, this produces
        //     au + bu^2 + cu^3 + du^4
        // but because u^4 = -3, we have
        //     -3d + au + bu^2 + cu^3

        Fp4 {
            c0: (&U4).sub(&self.c3),
            c1: self.c0,
            c2: self.c1,
            c3: self.c2,
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

        self.c3.lexicographically_largest()
            | (self.c3.is_zero() & self.c2.lexicographically_largest())
            | (self.c2.is_zero() & self.c1.lexicographically_largest())
            | (self.c1.is_zero() & self.c0.lexicographically_largest())
    }

    #[inline]
    /// Computes the multiplication of two Fp4 elements
    pub const fn mul(&self, other: &Fp4) -> Fp4 {
        let aa = (&self.c0).mul(&other.c0);
        let ab = (&self.c0).mul(&other.c1);
        let ac = (&self.c0).mul(&other.c2);
        let ad = (&self.c0).mul(&other.c3);

        let ba = (&self.c1).mul(&other.c0);
        let bb = (&self.c1).mul(&other.c1);
        let bc = (&self.c1).mul(&other.c2);
        let bd = (&self.c1).mul(&other.c3);

        let ca = (&self.c2).mul(&other.c0);
        let cb = (&self.c2).mul(&other.c1);
        let cc = (&self.c2).mul(&other.c2);
        let cd = (&self.c2).mul(&other.c3);

        let da = (&self.c3).mul(&other.c0);
        let db = (&self.c3).mul(&other.c1);
        let dc = (&self.c3).mul(&other.c2);
        let dd = (&self.c3).mul(&other.c3);

        let c0 = (&bd).add(&cc);
        let c0 = (&c0).add(&db);
        let c0 = (&c0).mul(&U4);
        let c0 = (&c0).add(&aa);

        let c1 = (&cd).add(&dc);
        let c1 = (&c1).mul(&U4);
        let c1 = (&c1).add(&ab);
        let c1 = (&c1).add(&ba);

        let c2 = (&dd).mul(&U4);
        let c2 = (&c2).add(&ac);
        let c2 = (&c2).add(&ca);
        let c2 = (&c2).add(&bb);

        let c3 = (&ad).add(&da);
        let c3 = (&c3).add(&bc);
        let c3 = (&c3).add(&cb);

        Fp4 { c0, c1, c2, c3 }
    }

    /// Square this element
    #[inline]
    pub const fn square(&self) -> Self {
        let aa = (&self.c0).mul(&self.c0);
        let ab = (&self.c0).mul(&self.c1);
        let ac = (&self.c0).mul(&self.c2);
        let ad = (&self.c0).mul(&self.c3);

        let bb = (&self.c1).mul(&self.c1);
        let bc = (&self.c1).mul(&self.c2);
        let bd = (&self.c1).mul(&self.c3);

        let cc = (&self.c2).mul(&self.c2);
        let cd = (&self.c2).mul(&self.c3);

        let dd = (&self.c3).mul(&self.c3);

        let c0 = (&bd).add(&cc);
        let c0 = (&c0).add(&bd);
        let c0 = (&c0).mul(&U4);
        let c0 = (&c0).add(&aa);

        let c1 = (&cd).add(&cd);
        let c1 = (&c1).mul(&U4);
        let c1 = (&c1).add(&ab);
        let c1 = (&c1).add(&ab);

        let c2 = (&dd).mul(&U4);
        let c2 = (&c2).add(&ac);
        let c2 = (&c2).add(&ac);
        let c2 = (&c2).add(&bb);

        let c3 = (&ad).add(&ad);
        let c3 = (&c3).add(&bc);
        let c3 = (&c3).add(&bc);

        Fp4 { c0, c1, c2, c3 }
    }

    /// Add two elements together
    #[inline]
    pub const fn add(&self, rhs: &Self) -> Self {
        Fp4 {
            c0: (&self.c0).add(&rhs.c0),
            c1: (&self.c1).add(&rhs.c1),
            c2: (&self.c2).add(&rhs.c2),
            c3: (&self.c3).add(&rhs.c3),
        }
    }

    /// Substract two elements together
    #[inline]
    pub const fn sub(&self, rhs: &Self) -> Self {
        Fp4 {
            c0: (&self.c0).sub(&rhs.c0),
            c1: (&self.c1).sub(&rhs.c1),
            c2: (&self.c2).sub(&rhs.c2),
            c3: (&self.c3).sub(&rhs.c3),
        }
    }

    /// Negate `&self`
    #[inline]
    pub const fn neg(&self) -> Self {
        Fp4 {
            c0: (&self.c0).neg(),
            c1: (&self.c1).neg(),
            c2: (&self.c2).neg(),
            c3: (&self.c3).neg(),
        }
    }

    /// Computes the multiplicative inverse of this field
    /// element, returning None in the case that this element
    /// is zero.
    pub fn invert(&self) -> CtOption<Self> {
        let a2 = (&self.c0).square();
        let b2 = (&self.c1).square();
        let c2 = (&self.c2).square();
        let d2 = (&self.c3).square();

        let three = Fp::new(3);
        let two = Fp::new(2);

        let inv = a2.square() + three * b2.square()
            - three.double().double() * self.c0 * b2 * self.c2
            + three.double() * a2 * c2
            + three.square() * c2.square()
            + three.square() * three * d2.square()
            + three.square().double() * (b2 + two * self.c0 * self.c2) * d2
            + three.double().double() * (a2 * self.c1 - three * self.c1 * c2) * self.c3;

        let c0 = a2 * self.c0 - three * b2 * self.c2
            + three * self.c0 * c2
            + three.double() * self.c0 * self.c1 * self.c3
            + three.square() * self.c2 * d2;
        let c1 = -(a2 * self.c1 - three * self.c1 * c2
            + three.square() * d2 * self.c3
            + three * (b2 + two * self.c0 * self.c2) * self.c3);
        let c3 = -(b2 * self.c1 - two * self.c0 * self.c1 * self.c2
            + three * self.c1 * d2
            + (a2 - three * self.c2.square()) * self.c3);
        let c2 = (&self.c0) * b2 - a2 * self.c2 - three * c2 * self.c2
            + three.double() * self.c1 * self.c2 * self.c3
            - three * self.c0 * d2;

        inv.invert().map(|t| Fp4 {
            c0: c0 * t,
            c1: c1 * t,
            c2: c2 * t,
            c3: c3 * t,
        })
    }

    /// Exponentiates `self` by `power`, where `power` is a
    /// little-endian order integer exponent.
    pub fn exp(self, by: &[u64; 4]) -> Self {
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
    pub fn exp_vartime(&self, by: &[u64; 4]) -> Self {
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
    pub const fn to_repr(&self) -> [u64; 4] {
        [
            Fp::montgomery_reduce(self.c0.0, 0).0,
            Fp::montgomery_reduce(self.c1.0, 0).0,
            Fp::montgomery_reduce(self.c2.0, 0).0,
            Fp::montgomery_reduce(self.c3.0, 0).0,
        ]
    }
}

#[test]
fn test_conditional_selection() {
    let a = Fp4 {
        c0: Fp::from_raw_unchecked(1),
        c1: Fp::from_raw_unchecked(2),
        c2: Fp::from_raw_unchecked(3),
        c3: Fp::from_raw_unchecked(4),
    };
    let b = Fp4 {
        c0: Fp::from_raw_unchecked(5),
        c1: Fp::from_raw_unchecked(6),
        c2: Fp::from_raw_unchecked(7),
        c3: Fp::from_raw_unchecked(8),
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
    fn is_equal(a: &Fp4, b: &Fp4) -> bool {
        let eq = a == b;
        let ct_eq = a.ct_eq(&b);

        assert_eq!(eq, bool::from(ct_eq));

        eq
    }

    assert!(is_equal(
        &Fp4 {
            c0: Fp::from_raw_unchecked(1),
            c1: Fp::from_raw_unchecked(2),
            c2: Fp::from_raw_unchecked(3),
            c3: Fp::from_raw_unchecked(4),
        },
        &Fp4 {
            c0: Fp::from_raw_unchecked(1),
            c1: Fp::from_raw_unchecked(2),
            c2: Fp::from_raw_unchecked(3),
            c3: Fp::from_raw_unchecked(4),
        }
    ));

    assert!(!is_equal(
        &Fp4 {
            c0: Fp::from_raw_unchecked(2),
            c1: Fp::from_raw_unchecked(2),
            c2: Fp::from_raw_unchecked(3),
            c3: Fp::from_raw_unchecked(4),
        },
        &Fp4 {
            c0: Fp::from_raw_unchecked(1),
            c1: Fp::from_raw_unchecked(2),
            c2: Fp::from_raw_unchecked(3),
            c3: Fp::from_raw_unchecked(4),
        }
    ));

    assert!(!is_equal(
        &Fp4 {
            c0: Fp::from_raw_unchecked(1),
            c1: Fp::from_raw_unchecked(3),
            c2: Fp::from_raw_unchecked(3),
            c3: Fp::from_raw_unchecked(4),
        },
        &Fp4 {
            c0: Fp::from_raw_unchecked(1),
            c1: Fp::from_raw_unchecked(2),
            c2: Fp::from_raw_unchecked(3),
            c3: Fp::from_raw_unchecked(4),
        }
    ));

    assert!(!is_equal(
        &Fp4 {
            c0: Fp::from_raw_unchecked(1),
            c1: Fp::from_raw_unchecked(2),
            c2: Fp::from_raw_unchecked(4),
            c3: Fp::from_raw_unchecked(4),
        },
        &Fp4 {
            c0: Fp::from_raw_unchecked(1),
            c1: Fp::from_raw_unchecked(2),
            c2: Fp::from_raw_unchecked(3),
            c3: Fp::from_raw_unchecked(4),
        }
    ));

    assert!(!is_equal(
        &Fp4 {
            c0: Fp::from_raw_unchecked(1),
            c1: Fp::from_raw_unchecked(2),
            c2: Fp::from_raw_unchecked(3),
            c3: Fp::from_raw_unchecked(5),
        },
        &Fp4 {
            c0: Fp::from_raw_unchecked(1),
            c1: Fp::from_raw_unchecked(2),
            c2: Fp::from_raw_unchecked(3),
            c3: Fp::from_raw_unchecked(4),
        }
    ));
}

#[test]
fn test_squaring() {
    let a = Fp4 {
        c0: Fp::new(3221830727238336732),
        c1: Fp::new(979900501602246277),
        c2: Fp::new(3348828664716707940),
        c3: Fp::new(553377494921525747),
    };
    let b = Fp4 {
        c0: Fp::new(4296344331630153542),
        c1: Fp::new(3962027865415014205),
        c2: Fp::new(1295707009280005994),
        c3: Fp::new(175739025260361652),
    };

    assert_eq!(a.square(), b);
}

#[test]
fn test_multiplication() {
    let a = Fp4 {
        c0: Fp::new(1),
        c1: Fp::new(2),
        c2: Fp::new(3),
        c3: Fp::new(4),
    };
    let b = Fp4::one();
    let c = Fp4 {
        c0: Fp::new(1),
        c1: Fp::new(2),
        c2: Fp::new(3),
        c3: Fp::new(4),
    };

    assert_eq!(a * b, c);

    let a = Fp4 {
        c0: Fp::new(3221830727238336732),
        c1: Fp::new(979900501602246277),
        c2: Fp::new(3348828664716707940),
        c3: Fp::new(553377494921525747),
    };
    let b = Fp4 {
        c0: Fp::new(1512850102547057140),
        c1: Fp::new(3247687537609432777),
        c2: Fp::new(2002622220347860596),
        c3: Fp::new(1845869432939790873),
    };
    let c = Fp4 {
        c0: Fp::new(3911041241527762547),
        c1: Fp::new(1196892283469028245),
        c2: Fp::new(4235984110621735126),
        c3: Fp::new(495889613393377298),
    };

    assert_eq!(a * b, c);
}

#[test]
fn test_addition() {
    let a = Fp4 {
        c0: Fp::new(1),
        c1: Fp::new(2),
        c2: Fp::new(3),
        c3: Fp::new(4),
    };
    let b = Fp4 {
        c0: Fp::new(4),
        c1: Fp::new(3),
        c2: Fp::new(2),
        c3: Fp::new(1),
    };
    let c = Fp4 {
        c0: Fp::new(5),
        c1: Fp::new(5),
        c2: Fp::new(5),
        c3: Fp::new(5),
    };

    assert_eq!(a + b, c);
}

#[test]
fn test_subtraction() {
    let a = Fp4 {
        c0: Fp::new(9),
        c1: Fp::new(8),
        c2: Fp::new(7),
        c3: Fp::new(6),
    };
    let b = Fp4 {
        c0: Fp::new(5),
        c1: Fp::new(4),
        c2: Fp::new(3),
        c3: Fp::new(2),
    };
    let c = Fp4 {
        c0: Fp::new(4),
        c1: Fp::new(4),
        c2: Fp::new(4),
        c3: Fp::new(4),
    };

    assert_eq!(a - b, c);
}

#[test]
fn test_negation() {
    let a = Fp4 {
        c0: Fp::new(1),
        c1: Fp::new(2),
        c2: Fp::new(3),
        c3: Fp::new(4),
    };
    let b = Fp4 {
        c0: -Fp::new(1),
        c1: -Fp::new(2),
        c2: -Fp::new(3),
        c3: -Fp::new(4),
    };

    assert_eq!(-a, b);
}

#[test]
fn test_inversion() {
    let a = Fp4 {
        c0: Fp::new(3221830727238336732),
        c1: Fp::new(979900501602246277),
        c2: Fp::new(3348828664716707940),
        c3: Fp::new(553377494921525747),
    };

    let b = Fp4 {
        c0: Fp::new(31598296486133907),
        c1: Fp::new(1556075319026159673),
        c2: Fp::new(683494480746306747),
        c3: Fp::new(697610636998943008),
    };

    assert_eq!(a.invert().unwrap(), b);

    assert!(bool::from(Fp4::zero().invert().is_none()));
}

#[test]
fn test_lexicographic_largest() {
    assert!(!bool::from(Fp4::zero().lexicographically_largest()));
    assert!(!bool::from(Fp4::one().lexicographically_largest()));
    assert!(bool::from(
        Fp4 {
            c0: Fp::new(1389794268293709605),
            c1: Fp::new(3631724493929800060),
            c2: Fp::new(1262796330815338397),
            c3: Fp::new(4058247500610520590),
        }
        .lexicographically_largest()
    ));
    assert!(!bool::from(
        Fp4 {
            c0: Fp::new(3221830727238336732),
            c1: Fp::new(979900501602246277),
            c2: Fp::new(3348828664716707940),
            c3: Fp::new(553377494921525747),
        }
        .lexicographically_largest()
    ));
    assert!(bool::from(
        Fp4 {
            c0: Fp::new(0),
            c1: Fp::new(0),
            c2: Fp::new(4058247500610520590),
            c3: Fp::new(0),
        }
        .lexicographically_largest()
    ));
}

#[cfg(feature = "zeroize")]
#[test]
fn test_zeroize() {
    use zeroize::Zeroize;

    let mut a = Fp4::one();
    a.zeroize();
    assert!(bool::from(a.is_zero()));
}
