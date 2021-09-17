use core::{
    fmt::{self, Debug, Display, Formatter},
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::utils::{add64_with_carry, mul64_with_carry, sub64_with_carry};

use rand_core::{CryptoRng, RngCore};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

// CONSTANTS
// ================================================================================================

// Field modulus = 2^62 - 111 * 2^39 + 1
const M: Fp = Fp(4611624995532046337);

/// 2^64 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R: Fp = Fp(244091581366268);

/// 2^128 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R2: Fp = Fp(630444561284293700);

/// Two-adicity of the field: (p-1) % 2^39 = 0
const TWO_ADICITY: u32 = 39;

// 2^39 root of unity = 4421547261963328785
//                    = 117700978803869913 in Montgomery form
const TWO_ADIC_ROOT_OF_UNITY: Fp = Fp(117700978803869913);

/// -M^{-1} mod 2^64; this is used during element multiplication.
const U: u64 = 4611624995532046335;

// FIELD ELEMENT
// ================================================================================================

/// Represents a base field element.
///
/// Internal values are stored in Montgomery representation and can be in the range [0; 2M).
/// The backing type is `u64`.
#[derive(Copy, Clone, Eq, Default)]
pub struct Fp(pub(crate) u64);

impl Debug for Fp {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        let tmp = self.to_repr();
        write!(f, "{:?}", tmp)
    }
}

impl Display for Fp {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl ConstantTimeEq for Fp {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0.ct_eq(&other.0)
    }
}

impl PartialEq for Fp {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl ConditionallySelectable for Fp {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Fp(u64::conditional_select(&a.0, &b.0, choice))
    }
}

impl zeroize::DefaultIsZeroes for Fp {}

impl Fp {
    /// Creates a new field element from a `u64` value.
    /// The value is converted to Montgomery form by computing
    /// (a.R^0 * R^2) / R = a.R
    pub const fn new(value: u64) -> Self {
        (&Fp(value)).mul(&R2)
    }

    /// Returns zero, the additive identity.
    #[inline]
    pub const fn zero() -> Self {
        Fp(0)
    }

    /// Returns one, the multiplicative identity.
    #[inline]
    pub const fn one() -> Self {
        R
    }

    /// Checks whether `self` is zero or not
    pub fn is_zero(&self) -> Choice {
        self.ct_eq(&Fp::zero())
    }

    #[inline(always)]
    #[allow(clippy::too_many_arguments)]
    pub(crate) const fn montgomery_reduce(r0: u64, r1: u64) -> Self {
        let k = r0.wrapping_mul(U);
        let (_, carry) = mul64_with_carry(r0, k, M.0, 0);
        let (r1, _) = add64_with_carry(r1, 0, carry);

        // Result may be within M of the correct value, hence substracting the modulus
        (&Fp(r1)).sub(&M)
    }

    /// Computes the summation of two field elements
    #[inline]
    pub const fn add(&self, rhs: &Self) -> Self {
        let (d0, _) = add64_with_carry(self.0, rhs.0, 0);

        // Attempt to subtract the modulus, to ensure the value
        // is smaller than the modulus.
        (&Fp(d0)).sub(&M)
    }

    /// Computes the difference of two field elements
    #[inline]
    pub const fn sub(&self, rhs: &Self) -> Self {
        let (d0, borrow) = sub64_with_carry(self.0, rhs.0, 0);

        // If underflow occurred on the final limb,
        // borrow = 0xfff...fff, otherwise borrow = 0x000...000.
        let (d0, _) = add64_with_carry(d0, M.0 & borrow, 0);

        Fp(d0)
    }

    /// Computes the negation of a field element
    #[inline]
    pub const fn neg(&self) -> Self {
        // Subtract `self` from `MODULUS` to negate. Ignore the final
        // borrow because it cannot underflow; self is guaranteed to
        // be in the field.
        let (d0, _) = sub64_with_carry(M.0, self.0, 0);

        // `tmp` could be `MODULUS` if `self` was zero. Create a mask that is
        // zero if `self` was zero, and `u64::max_value()` if self was nonzero.
        let mask = ((self.0 == 0) as u64).wrapping_sub(1);

        Fp(d0 & mask)
    }

    /// Computes the multiplication of two field elements
    #[inline]
    pub const fn mul(&self, rhs: &Self) -> Self {
        // Schoolbook multiplication

        let (r0, r1) = mul64_with_carry(0, self.0, rhs.0, 0);

        Fp::montgomery_reduce(r0, r1)
    }

    /// Computes the square of a field element
    #[inline]
    pub const fn square(&self) -> Self {
        self.mul(&self)
    }

    /// Computes the double of a field element
    // Can be faster via bitshift
    #[inline]
    pub const fn double(&self) -> Self {
        self.add(self)
    }

    /// Outputs the internal representation as 4 64-bit limbs after Montgomery reduction
    pub const fn to_repr(&self) -> u64 {
        Fp::montgomery_reduce(self.0, 0).0
    }

    /// Converts a `Fp` element into a byte representation in
    /// little-endian byte order.
    pub fn to_bytes(&self) -> [u8; 8] {
        // Turn into canonical form by computing
        // (a.R) / R = a
        let tmp = Fp::montgomery_reduce(self.0, 0);

        tmp.0.to_le_bytes()
    }

    /// Returns whether or not this element is strictly lexicographically
    /// larger than its negation.
    pub fn lexicographically_largest(&self) -> Choice {
        // This can be determined by checking to see if the element is
        // larger than (p - 1) // 2. If we subtract by ((p - 1) // 2) + 1
        // and there is no underflow, then the element must be larger than
        // (p - 1) // 2.

        // First, because self is in Montgomery form we need to reduce it
        let tmp = Fp::montgomery_reduce(self.0, 0);

        let (_, borrow) = sub64_with_carry(tmp.0, 2305812497766023169, 0);

        // If the element was smaller, the subtraction will underflow
        // producing a borrow value of 0xffff...ffff, otherwise it will
        // be zero. We create a Choice representing true if there was
        // overflow (and so this element is not lexicographically larger
        // than its negation) and then negate it.

        !Choice::from((borrow as u8) & 1)
    }

    /// Exponentiates `self` by `power`, where `power` is a
    /// little-endian order integer exponent.
    pub fn exp(self, power: u64) -> Self {
        let mut res = Self::one();
        for i in (0..64).rev() {
            res = res.square();
            let mut tmp = res;
            tmp *= self;
            res.conditional_assign(&tmp, (((power >> i) & 1) as u8).into());
        }
        res
    }

    /// Exponentiates `self` by `power`, where `power` is a
    /// little-endian order integer exponent.
    ///
    /// **This operation is variable time with respect
    /// to the exponent.** If the exponent is fixed,
    /// this operation is effectively constant time.
    pub fn exp_vartime(&self, power: u64) -> Self {
        let mut res = Self::one();
        for i in (0..64).rev() {
            res = res.square();

            if ((power >> i) & 1) == 1 {
                res.mul_assign(self);
            }
        }
        res
    }

    /// Computes the multiplicative inverse of this element,
    /// failing if the element is zero.
    pub fn invert(self) -> CtOption<Self> {
        #[inline(always)]
        fn square_assign_multi(n: &mut Fp, num_times: usize) {
            for _ in 0..num_times {
                *n = n.square();
            }
        }
        // found using https://github.com/kwantam/addchain for M - 2
        let mut t2 = self.square(); //       1:   2
        let mut t1 = t2.square(); //         2:   4
        t1 *= &t2; //                        3:   6
        t2 = t1.square(); //                 4:   12
        t2 = t2.square(); //                 5:   24
        t2 *= &t1; //                        6:   30
        square_assign_multi(&mut t2, 2); //  8:   120
        t2 *= &t1; //                        9:   126
        square_assign_multi(&mut t2, 2); //  11:  504
        let t3 = t2 * self; //               12:  505
        t1 *= &t3; //                        13:  511
        let mut t0 = t1 * self; //           14:  512
        t2 = t0 * t1; //                     15:  1023
        t0 = t2.square(); //                 16:  2046
        square_assign_multi(&mut t0, 8); //  24:  523776
        t0 *= &t3; //                        25:  524281
        square_assign_multi(&mut t0, 14); // 39:  8589819904
        t0 *= &t2; //                        40:  8589820927
        square_assign_multi(&mut t0, 10); // 50:  8795976629248
        t0 *= &t2; //                        51:  8795976630271
        square_assign_multi(&mut t0, 10); // 61:  9007080069397504
        t0 *= &t2; //                        62:  9007080069398527
        square_assign_multi(&mut t0, 9); //  71:  4611624995532045824
        t0 *= &t1; //                        72:  4611624995532046335

        CtOption::new(t0, !self.ct_eq(&Self::zero()))
    }

    /// Computes the conjugate of a `Fp` element
    pub fn conjugate(&self) -> Self {
        Fp(self.0)
    }

    /// Computes a random `FieldElement` element
    pub fn random(mut rng: impl RngCore + CryptoRng) -> Self {
        Fp::new(rng.next_u64())
    }

    #[inline(always)]
    /// Normalizes the internal representation of a `Fp` element
    pub fn normalize(&self) -> Self {
        self * R2
    }

    /// Constructs an element of `Fp` without checking that it is
    /// canonical.
    pub const fn from_raw_unchecked(v: u64) -> Self {
        Fp(v)
    }

    /// Outputs a `Fp` element of multiplicative order equals to 2^n
    pub fn get_root_of_unity(n: u32) -> Self {
        assert!(n != 0, "cannot get root of unity for n = 0");
        assert!(n <= TWO_ADICITY, "order cannot exceed 2^{}", TWO_ADICITY);
        let power = 1u64 << (TWO_ADICITY - n);

        TWO_ADIC_ROOT_OF_UNITY.exp(power)
    }
}

// OVERLOADED OPERATORS
// ================================================================================================

impl<'a> Neg for &'a Fp {
    type Output = Fp;

    #[inline]
    fn neg(self) -> Fp {
        self.neg()
    }
}

impl Neg for Fp {
    type Output = Fp;

    #[inline]
    fn neg(self) -> Fp {
        -&self
    }
}

impl<'a, 'b> Sub<&'b Fp> for &'a Fp {
    type Output = Fp;

    #[inline]
    fn sub(self, rhs: &'b Fp) -> Fp {
        self.sub(rhs)
    }
}

impl<'a, 'b> Add<&'b Fp> for &'a Fp {
    type Output = Fp;

    #[inline]
    fn add(self, rhs: &'b Fp) -> Fp {
        self.add(rhs)
    }
}

impl<'a, 'b> Mul<&'b Fp> for &'a Fp {
    type Output = Fp;

    #[inline]
    fn mul(self, rhs: &'b Fp) -> Fp {
        self.mul(rhs)
    }
}

impl_binops_additive!(Fp, Fp);
impl_binops_multiplicative!(Fp, Fp);

impl Div for Fp {
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Fp {
        self.mul(&rhs.invert().unwrap_or(Self::zero()))
    }
}

impl DivAssign for Fp {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs
    }
}

// TYPE CONVERSIONS
// ================================================================================================

impl From<u64> for Fp {
    /// Converts a 64-bit value into a filed element. If the value is greater than or equal to
    /// the field modulus, modular reduction is silently preformed.
    fn from(value: u64) -> Self {
        Fp::new(value)
    }
}

impl From<u32> for Fp {
    /// Converts a 32-bit value into a filed element.
    fn from(value: u32) -> Self {
        Fp::new(value as u64)
    }
}

impl From<u16> for Fp {
    /// Converts a 16-bit value into a filed element.
    fn from(value: u16) -> Self {
        Fp::new(value as u64)
    }
}

impl From<u8> for Fp {
    /// Converts an 8-bit value into a filed element.
    fn from(value: u8) -> Self {
        Fp::new(value as u64)
    }
}

impl From<[u8; 8]> for Fp {
    /// Converts the value encoded in an array of 8 bytes into a field element. The bytes are
    /// assumed to encode the element in the canonical representation in little-endian byte order.
    /// If the value is greater than or equal to the field modulus, modular reduction is silently
    /// preformed.
    fn from(bytes: [u8; 8]) -> Self {
        let value = u64::from_le_bytes(bytes);
        Fp::new(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::thread_rng;

    const LARGEST: Fp = Fp(4611624995532046336);
    const TWO_POW_39: u64 = 549755813888;
    const TWO_POW_38: u64 = 274877906944;

    // DISPLAY
    // ================================================================================================

    #[test]
    fn test_debug() {
        assert_eq!(format!("{:?}", Fp::zero()), "0");
        assert_eq!(format!("{:?}", Fp::one()), "1");
        assert_eq!(format!("{:?}", R2), "244091581366268");
    }

    #[test]
    fn test_to_repr() {
        assert_eq!(format!("{:?}", Fp::zero().to_repr()), "0");
        assert_eq!(format!("{:?}", Fp::one().to_repr()), "1");
        assert_eq!(format!("{:?}", R2.to_repr()), "244091581366268");

        let r = Fp::from_raw_unchecked(244091581366268);
        assert!(format!("{:?}", r.0) != format!("{:?}", R.to_repr()));

        assert_eq!(format!("{:?}", r.to_repr()), format!("{:?}", R.to_repr()));
    }

    // BASIC ALGEBRA
    // ================================================================================================

    #[test]
    fn test_equality() {
        assert_eq!(Fp::default(), Fp::zero());
        assert_eq!(Fp::zero(), Fp::zero());
        assert_eq!(Fp::one(), Fp::one());

        assert!(bool::from(Fp::default().is_zero()));
        assert!(!bool::from(Fp::zero().ct_eq(&Fp::one())));

        assert!(Fp::zero() != Fp::one());
        assert!(Fp::one() != R2);
    }

    #[test]
    fn test_addition() {
        let mut tmp = LARGEST;
        tmp += &LARGEST;

        assert_eq!(tmp, Fp(4611624995532046335));

        assert_eq!(tmp, LARGEST.double());

        let mut tmp = LARGEST;
        tmp += &Fp(1);

        assert_eq!(tmp, Fp::zero());
    }

    #[test]
    fn test_subtraction() {
        let mut tmp = LARGEST;
        tmp -= &LARGEST;

        assert_eq!(tmp, Fp::zero());

        let mut tmp = Fp::zero();
        tmp -= &LARGEST;

        let mut tmp2 = M;
        tmp2 -= &LARGEST;

        assert_eq!(tmp, tmp2);
    }

    #[test]
    fn test_negation() {
        let tmp = -&LARGEST;

        assert_eq!(tmp, Fp(1));

        let tmp = -&Fp::zero();
        assert_eq!(tmp, Fp::zero());
        let tmp = -&Fp(1);
        assert_eq!(tmp, LARGEST);
    }

    #[test]
    fn test_multiplication() {
        let mut cur = LARGEST;

        for _ in 0..100 {
            let mut tmp = cur;
            tmp *= &cur;

            let mut tmp2 = Fp::zero();
            for b in cur
                .to_bytes()
                .iter()
                .rev()
                .flat_map(|byte| (0..8).rev().map(move |i| ((byte >> i) & 1u8) == 1u8))
            {
                let tmp3 = tmp2;
                tmp2.add_assign(&tmp3);

                if b {
                    tmp2.add_assign(&cur);
                }
            }

            assert_eq!(tmp, tmp2);

            cur.add_assign(&LARGEST);
        }
    }

    #[test]

    fn test_inversion() {
        assert!(bool::from(Fp::zero().invert().is_none()));
        assert_eq!(Fp::one().invert().unwrap(), Fp::one());
        assert_eq!((-&Fp::one()).invert().unwrap(), -&Fp::one());

        let mut tmp = Fp::random(&mut thread_rng());

        for _ in 0..100 {
            let mut tmp2 = tmp.invert().unwrap();
            tmp2.mul_assign(&tmp);

            assert_eq!(tmp2, Fp::one());

            tmp.add_assign(&Fp::random(&mut thread_rng()));
        }
    }

    #[test]
    fn test_squaring() {
        let mut cur = LARGEST;

        for _ in 0..100 {
            let mut tmp = cur;
            let pow2 = tmp.exp(2);
            tmp = tmp.square();

            let mut tmp2 = Fp::zero();
            for b in cur
                .to_bytes()
                .iter()
                .rev()
                .flat_map(|byte| (0..8).rev().map(move |i| ((byte >> i) & 1u8) == 1u8))
            {
                let tmp3 = tmp2;
                tmp2.add_assign(&tmp3);

                if b {
                    tmp2.add_assign(&cur);
                }
            }

            assert_eq!(tmp, tmp2);
            assert_eq!(tmp, pow2);

            cur.add_assign(&LARGEST);
        }
    }

    #[test]
    fn test_invert_is_pow() {
        let q_minus_2 = 4611624995532046335;

        let mut r1 = R;
        let mut r2 = R;
        let mut r3 = R;

        for _ in 0..100 {
            r1 = r1.invert().unwrap();
            r2 = r2.exp(q_minus_2);
            r3 = r3.exp_vartime(q_minus_2);

            assert_eq!(r1, r2);
            assert_eq!(r2, r3);
            // Add R so we check something different next time around
            r1.add_assign(&R);
            r2 = r1;
            r3 = r1;
        }
    }

    #[test]
    fn test_conjugate() {
        let a = Fp::random(&mut thread_rng());
        let b = a.conjugate();
        assert_eq!(a, b);
    }

    // ROOTS OF UNITY
    // ================================================================================================

    #[test]
    fn test_get_root_of_unity() {
        let root_39 = Fp::get_root_of_unity(39);
        assert_eq!(TWO_ADIC_ROOT_OF_UNITY, root_39);
        assert_eq!(Fp::one(), root_39.exp(TWO_POW_39));

        let root_38 = Fp::get_root_of_unity(38);
        let expected = root_39.exp(2);
        assert_eq!(expected, root_38);
        assert_eq!(Fp::one(), root_38.exp(TWO_POW_38));
    }

    // SERIALIZATION / DESERIALIZATION
    // ================================================================================================

    #[test]
    fn test_to_bytes() {
        assert_eq!(Fp::zero().to_bytes(), [0, 0, 0, 0, 0, 0, 0, 0,]);

        assert_eq!(Fp::one().to_bytes(), [1, 0, 0, 0, 0, 0, 0, 0,]);

        assert_eq!(R2.to_bytes(), [252, 255, 255, 255, 255, 221, 0, 0]);

        assert_eq!((-&Fp::one()).to_bytes(), [0, 0, 0, 0, 128, 200, 255, 63]);
    }

    #[test]
    fn test_lexicographically_largest() {
        let a = Fp::from_raw_unchecked(244091581366268);

        let b = Fp::from_raw_unchecked(4611380903950680069);

        assert_eq!(a.square(), b.square());
        assert!(!bool::from(a.lexicographically_largest()));
        assert!(bool::from(b.lexicographically_largest()));
    }

    #[test]
    fn test_zeroize() {
        use zeroize::Zeroize;

        let mut a = Fp::one();
        a.zeroize();
        assert_eq!(a, Fp::zero());
    }

    // INITIALIZATION
    // ================================================================================================

    #[test]
    fn test_from_int() {
        let n = 42u8;
        let element = Fp::from(n);

        assert_eq!(element, Fp::from(n as u16));
        assert_eq!(element, Fp::from(n as u32));
        assert_eq!(element, Fp::from(n as u64));
    }

    #[test]
    fn test_from_raw_unchecked() {
        let mut element = Fp::from_raw_unchecked(244091581366268);

        let element_normalized = Fp::new(244091581366268);

        assert_eq!(element, Fp::one());
        element = element.normalize();

        assert!(element != Fp::one());
        assert_eq!(element, element_normalized);
    }
}
