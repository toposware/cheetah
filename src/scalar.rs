use core::{
    convert::{TryFrom, TryInto},
    fmt::{self, Debug, Display, Formatter},
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::utils::{add64_with_carry, mul64_with_carry, sub64_with_carry};

use bitvec::{order::Lsb0, slice::BitSlice};
use rand_core::{CryptoRng, RngCore};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "serialize")]
use serde::de::Visitor;
#[cfg(feature = "serialize")]
use serde::{self, Deserialize, Deserializer, Serialize, Serializer};

// CONSTANTS
// ================================================================================================

// Field modulus = 452288908625799774674098781249295805571160305918807263004273270801344763473
const M: Scalar = Scalar([
    0x02f096739c941651,
    0x2ef8715df33bd3e0,
    0x46b09ca43418b591,
    0x00fffc8804831564,
]);

/// 2^256 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R: Scalar = Scalar([
    0x0f698c636be9af00,
    0x078ea20cc42c1ffd,
    0x4f635bcbe74a6ed1,
    0x000377fb7cea9bb9,
]);

/// 2^512 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R2: Scalar = Scalar([
    0xb598f7a893b7ab5a,
    0xced8a71a39a3d46e,
    0x997959903d1e6b70,
    0x005fa27e825555e8,
]);

/// 2^768 mod M; this is used during element inversion.
pub(crate) const R3: Scalar = Scalar([
    0xf116e1df1d4f7994,
    0x7894e4d7b5c9d436,
    0xaed9f68f025ae2a7,
    0x00154013fc66b6c6,
]);

/// -M^{-1} mod 2^64; this is used during element multiplication.
const U: u64 = 0xc0000a221cbb0d4f;

// SCALAR FIELD ELEMENT
// ================================================================================================

/// Represents a scalar field element.
///
/// Internal values are stored in their canonical form in the range [0, M).
/// The backing type is `[u64; 4]`.
#[derive(Copy, Clone, Eq)]
pub struct Scalar(pub(crate) [u64; 4]);

impl Default for Scalar {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

impl Debug for Scalar {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        let tmp = self.to_bytes();
        write!(f, "0x")?;
        for &b in tmp.iter().rev() {
            write!(f, "{:02x}", b)?;
        }
        Ok(())
    }
}

impl Display for Scalar {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl ConstantTimeEq for Scalar {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0[0].ct_eq(&other.0[0])
            & self.0[1].ct_eq(&other.0[1])
            & self.0[2].ct_eq(&other.0[2])
            & self.0[3].ct_eq(&other.0[3])
    }
}

impl PartialEq for Scalar {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl ConditionallySelectable for Scalar {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Scalar([
            u64::conditional_select(&a.0[0], &b.0[0], choice),
            u64::conditional_select(&a.0[1], &b.0[1], choice),
            u64::conditional_select(&a.0[2], &b.0[2], choice),
            u64::conditional_select(&a.0[3], &b.0[3], choice),
        ])
    }
}

impl zeroize::DefaultIsZeroes for Scalar {}

impl Scalar {
    /// Creates a new field element from a [u64; 4] value.
    /// The value is converted to Montgomery form by computing
    /// (a.R^0 * R^2) / R = a.R
    pub const fn new(value: [u64; 4]) -> Self {
        (&Scalar(value)).mul(&R2)
    }

    /// Returns zero, the additive identity.
    #[inline]
    pub const fn zero() -> Self {
        Scalar([0, 0, 0, 0])
    }

    /// Returns one, the multiplicative identity.
    #[inline]
    pub const fn one() -> Self {
        R
    }

    /// Checks whether `self` is zero or not
    pub fn is_zero(&self) -> Choice {
        self.ct_eq(&Scalar::zero())
    }

    #[inline(always)]
    #[allow(clippy::too_many_arguments)]
    const fn montgomery_reduce(
        r0: u64,
        r1: u64,
        r2: u64,
        r3: u64,
        r4: u64,
        r5: u64,
        r6: u64,
        r7: u64,
    ) -> Self {
        let k = r0.wrapping_mul(U);
        let (_, carry) = mul64_with_carry(r0, k, M.0[0], 0);
        let (r1, carry) = mul64_with_carry(r1, k, M.0[1], carry);
        let (r2, carry) = mul64_with_carry(r2, k, M.0[2], carry);
        let (r3, carry) = mul64_with_carry(r3, k, M.0[3], carry);
        let (r4, carry2) = add64_with_carry(r4, 0, carry);

        let k = r1.wrapping_mul(U);
        let (_, carry) = mul64_with_carry(r1, k, M.0[0], 0);
        let (r2, carry) = mul64_with_carry(r2, k, M.0[1], carry);
        let (r3, carry) = mul64_with_carry(r3, k, M.0[2], carry);
        let (r4, carry) = mul64_with_carry(r4, k, M.0[3], carry);
        let (r5, carry2) = add64_with_carry(r5, carry2, carry);

        let k = r2.wrapping_mul(U);
        let (_, carry) = mul64_with_carry(r2, k, M.0[0], 0);
        let (r3, carry) = mul64_with_carry(r3, k, M.0[1], carry);
        let (r4, carry) = mul64_with_carry(r4, k, M.0[2], carry);
        let (r5, carry) = mul64_with_carry(r5, k, M.0[3], carry);
        let (r6, carry2) = add64_with_carry(r6, carry2, carry);

        let k = r3.wrapping_mul(U);
        let (_, carry) = mul64_with_carry(r3, k, M.0[0], 0);
        let (r4, carry) = mul64_with_carry(r4, k, M.0[1], carry);
        let (r5, carry) = mul64_with_carry(r5, k, M.0[2], carry);
        let (r6, carry) = mul64_with_carry(r6, k, M.0[3], carry);
        let (r7, _) = add64_with_carry(r7, carry2, carry);

        // Result may be within M of the correct value, hence substracting the modulus
        (&Scalar([r4, r5, r6, r7])).sub(&M)
    }

    /// Computes the summation of two scalar elements
    #[inline]
    pub const fn add(&self, rhs: &Self) -> Self {
        let (d0, carry) = add64_with_carry(self.0[0], rhs.0[0], 0);
        let (d1, carry) = add64_with_carry(self.0[1], rhs.0[1], carry);
        let (d2, carry) = add64_with_carry(self.0[2], rhs.0[2], carry);
        let (d3, _) = add64_with_carry(self.0[3], rhs.0[3], carry);

        // Attempt to subtract the modulus, to ensure the value
        // is smaller than the modulus.
        (&Scalar([d0, d1, d2, d3])).sub(&M)
    }

    /// Computes the difference of two scalar elements
    #[inline]
    pub const fn sub(&self, rhs: &Self) -> Self {
        let (d0, borrow) = sub64_with_carry(self.0[0], rhs.0[0], 0);
        let (d1, borrow) = sub64_with_carry(self.0[1], rhs.0[1], borrow);
        let (d2, borrow) = sub64_with_carry(self.0[2], rhs.0[2], borrow);
        let (d3, borrow) = sub64_with_carry(self.0[3], rhs.0[3], borrow);

        // If underflow occurred on the final limb,
        // borrow = 0xfff...fff, otherwise borrow = 0x000...000.
        let (d0, carry) = add64_with_carry(d0, M.0[0] & borrow, 0);
        let (d1, carry) = add64_with_carry(d1, M.0[1] & borrow, carry);
        let (d2, carry) = add64_with_carry(d2, M.0[2] & borrow, carry);
        let (d3, _) = add64_with_carry(d3, M.0[3] & borrow, carry);

        Scalar([d0, d1, d2, d3])
    }

    /// Computes the negation of a scalar element
    #[inline]
    pub const fn neg(&self) -> Self {
        // Subtract `self` from `MODULUS` to negate. Ignore the final
        // borrow because it cannot underflow; self is guaranteed to
        // be in the field.
        let (d0, borrow) = sub64_with_carry(M.0[0], self.0[0], 0);
        let (d1, borrow) = sub64_with_carry(M.0[1], self.0[1], borrow);
        let (d2, borrow) = sub64_with_carry(M.0[2], self.0[2], borrow);
        let (d3, _) = sub64_with_carry(M.0[3], self.0[3], borrow);

        // `tmp` could be `MODULUS` if `self` was zero. Create a mask that is
        // zero if `self` was zero, and `u64::max_value()` if self was nonzero.
        let mask = (((self.0[0] | self.0[1] | self.0[2] | self.0[3]) == 0) as u64).wrapping_sub(1);

        Scalar([d0 & mask, d1 & mask, d2 & mask, d3 & mask])
    }

    /// Computes the multiplication of two scalar elements
    #[inline]
    pub const fn mul(&self, rhs: &Self) -> Self {
        // Schoolbook multiplication

        let (r0, carry) = mul64_with_carry(0, self.0[0], rhs.0[0], 0);
        let (r1, carry) = mul64_with_carry(0, self.0[0], rhs.0[1], carry);
        let (r2, carry) = mul64_with_carry(0, self.0[0], rhs.0[2], carry);
        let (r3, r4) = mul64_with_carry(0, self.0[0], rhs.0[3], carry);

        let (r1, carry) = mul64_with_carry(r1, self.0[1], rhs.0[0], 0);
        let (r2, carry) = mul64_with_carry(r2, self.0[1], rhs.0[1], carry);
        let (r3, carry) = mul64_with_carry(r3, self.0[1], rhs.0[2], carry);
        let (r4, r5) = mul64_with_carry(r4, self.0[1], rhs.0[3], carry);

        let (r2, carry) = mul64_with_carry(r2, self.0[2], rhs.0[0], 0);
        let (r3, carry) = mul64_with_carry(r3, self.0[2], rhs.0[1], carry);
        let (r4, carry) = mul64_with_carry(r4, self.0[2], rhs.0[2], carry);
        let (r5, r6) = mul64_with_carry(r5, self.0[2], rhs.0[3], carry);

        let (r3, carry) = mul64_with_carry(r3, self.0[3], rhs.0[0], 0);
        let (r4, carry) = mul64_with_carry(r4, self.0[3], rhs.0[1], carry);
        let (r5, carry) = mul64_with_carry(r5, self.0[3], rhs.0[2], carry);
        let (r6, r7) = mul64_with_carry(r6, self.0[3], rhs.0[3], carry);

        Scalar::montgomery_reduce(r0, r1, r2, r3, r4, r5, r6, r7)
    }

    /// Computes the square of a scalar element
    #[inline]
    pub const fn square(&self) -> Self {
        let (r1, carry) = mul64_with_carry(0, self.0[0], self.0[1], 0);
        let (r2, carry) = mul64_with_carry(0, self.0[0], self.0[2], carry);
        let (r3, r4) = mul64_with_carry(0, self.0[0], self.0[3], carry);

        let (r3, carry) = mul64_with_carry(r3, self.0[1], self.0[2], 0);
        let (r4, r5) = mul64_with_carry(r4, self.0[1], self.0[3], carry);

        let (r5, r6) = mul64_with_carry(r5, self.0[2], self.0[3], 0);

        let r7 = r6 >> 63;
        let r6 = (r6 << 1) | (r5 >> 63);
        let r5 = (r5 << 1) | (r4 >> 63);
        let r4 = (r4 << 1) | (r3 >> 63);
        let r3 = (r3 << 1) | (r2 >> 63);
        let r2 = (r2 << 1) | (r1 >> 63);
        let r1 = r1 << 1;

        let (r0, carry) = mul64_with_carry(0, self.0[0], self.0[0], 0);
        let (r1, carry) = add64_with_carry(0, r1, carry);
        let (r2, carry) = mul64_with_carry(r2, self.0[1], self.0[1], carry);
        let (r3, carry) = add64_with_carry(0, r3, carry);
        let (r4, carry) = mul64_with_carry(r4, self.0[2], self.0[2], carry);
        let (r5, carry) = add64_with_carry(0, r5, carry);
        let (r6, carry) = mul64_with_carry(r6, self.0[3], self.0[3], carry);
        let (r7, _) = add64_with_carry(0, r7, carry);

        Scalar::montgomery_reduce(r0, r1, r2, r3, r4, r5, r6, r7)
    }

    /// Computes the double of a scalar element
    // Can be faster via bitshift
    #[inline]
    pub const fn double(&self) -> Self {
        self.add(self)
    }

    /// Attempts to convert a little-endian byte representation of
    /// a scalar into a `Scalar`, failing if the input is not canonical.
    pub fn from_bytes(bytes: &[u8; 32]) -> CtOption<Scalar> {
        let mut tmp = Scalar([0, 0, 0, 0]);

        tmp.0[0] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[0..8]).unwrap());
        tmp.0[1] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[8..16]).unwrap());
        tmp.0[2] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[16..24]).unwrap());
        tmp.0[3] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[24..32]).unwrap());

        // Try to subtract the modulus
        let (_, borrow) = sub64_with_carry(tmp.0[0], M.0[0], 0);
        let (_, borrow) = sub64_with_carry(tmp.0[1], M.0[1], borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[2], M.0[2], borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[3], M.0[3], borrow);

        // If the element is smaller than MODULUS then the
        // subtraction will underflow, producing a borrow value
        // of 0xffff...ffff. Otherwise, it'll be zero.
        let is_some = (borrow as u8) & 1;

        // Convert to Montgomery form by computing
        // (a.R^0 * R^2) / R = a.R
        tmp *= &R2;

        CtOption::new(tmp, Choice::from(is_some))
    }

    /// Convert a little-endian bit sequence into a Scalar element
    ///
    /// **This operation is variable time with respect
    /// to the bit slice.** If the slice is fixed,
    /// this operation is effectively constant time.
    pub fn from_bits_vartime(bit_slice: &BitSlice<Lsb0, u8>) -> Scalar {
        assert_eq!(bit_slice.len(), 256);

        let mut result = Scalar::zero();
        for i in 0..256 {
            result = result.double();
            if bit_slice[255 - i] {
                result += Scalar::one();
            }
        }

        result
    }

    /// Outputs the internal representation as 4 64-bit limbs after Montgomery reduction
    pub const fn to_repr(&self) -> [u64; 4] {
        Scalar::montgomery_reduce(self.0[0], self.0[1], self.0[2], self.0[3], 0, 0, 0, 0).0
    }

    /// Converts a `Scalar` element into a byte representation in
    /// little-endian byte order.
    pub fn to_bytes(&self) -> [u8; 32] {
        // Turn into canonical form by computing
        // (a.R) / R = a
        let tmp = Scalar::montgomery_reduce(self.0[0], self.0[1], self.0[2], self.0[3], 0, 0, 0, 0);

        let mut res = [0; 32];
        res[0..8].copy_from_slice(&tmp.0[0].to_le_bytes());
        res[8..16].copy_from_slice(&tmp.0[1].to_le_bytes());
        res[16..24].copy_from_slice(&tmp.0[2].to_le_bytes());
        res[24..32].copy_from_slice(&tmp.0[3].to_le_bytes());

        res
    }

    /// Converts a 512-bit little endian integer into
    /// a `Scalar` by reducing by the modulus.
    pub fn from_bytes_wide(bytes: &[u8; 64]) -> Self {
        Scalar::from_u512([
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[0..8]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[8..16]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[16..24]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[24..32]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[32..40]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[40..48]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[48..56]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[56..64]).unwrap()),
        ])
    }

    fn from_u512(limbs: [u64; 8]) -> Self {
        // We reduce an arbitrary 512-bit number by decomposing it into two 256-bit digits
        // with the higher bits multiplied by 2^256. Thus, we perform two reductions
        //
        // 1. the lower bits are multiplied by R^2, as normal
        // 2. the upper bits are multiplied by R^2 * 2^256 = R^3
        //
        // and computing their sum in the field. It remains to see that arbitrary 256-bit
        // numbers can be placed into Montgomery form safely using the reduction. The
        // reduction works so long as the product is less than R=2^256 multiplied by
        // the modulus. This holds because for any `c` smaller than the modulus, we have
        // that (2^256 - 1)*c is an acceptable product for the reduction. Therefore, the
        // reduction always works so long as `c` is in the field; in this case it is either the
        // constant `R2` or `R3`.
        let d0 = Scalar([limbs[0], limbs[1], limbs[2], limbs[3]]);
        let d1 = Scalar([limbs[4], limbs[5], limbs[6], limbs[7]]);

        // Convert to Montgomery form
        d0 * R2 + d1 * R3
    }

    /// Returns whether or not this element is strictly lexicographically
    /// larger than its negation.
    pub fn lexicographically_largest(&self) -> Choice {
        // This can be determined by checking to see if the element is
        // larger than (p - 1) // 2. If we subtract by ((p - 1) // 2) + 1
        // and there is no underflow, then the element must be larger than
        // (p - 1) // 2.

        // First, because self is in Montgomery form we need to reduce it
        let tmp = Scalar::montgomery_reduce(self.0[0], self.0[1], self.0[2], self.0[3], 0, 0, 0, 0);

        let (_, borrow) = sub64_with_carry(tmp.0[0], 0x0178_4b39_ce4a_0b29, 0);
        let (_, borrow) = sub64_with_carry(tmp.0[1], 0x977c_38ae_f99d_e9f0, borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[2], 0x2358_4e52_1a0c_5ac8, borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[3], 0x007f_fe44_0241_8ab2, borrow);

        // If the element was smaller, the subtraction will underflow
        // producing a borrow value of 0xffff...ffff, otherwise it will
        // be zero. We create a Choice representing true if there was
        // overflow (and so this element is not lexicographically larger
        // than its negation) and then negate it.

        !Choice::from((borrow as u8) & 1)
    }

    /// Exponentiates `self` by `power`, where `power` is a
    /// little-endian order integer exponent.
    pub fn exp(self, power: &[u64; 4]) -> Self {
        let mut res = Self::one();
        for e in power.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();
                let mut tmp = res;
                tmp *= self;
                res.conditional_assign(&tmp, (((*e >> i) & 1) as u8).into());
            }
        }
        res
    }

    /// Exponentiates `self` by `power`, where `power` is a
    /// little-endian order integer exponent.
    ///
    /// **This operation is variable time with respect
    /// to the exponent.** If the exponent is fixed,
    /// this operation is effectively constant time.
    pub fn exp_vartime(&self, power: &[u64; 4]) -> Self {
        let mut res = Self::one();
        for e in power.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();

                if ((*e >> i) & 1) == 1 {
                    res.mul_assign(self);
                }
            }
        }
        res
    }

    /// Computes the multiplicative inverse of this element,
    /// failing if the element is zero.
    pub fn invert(self) -> CtOption<Self> {
        #[inline(always)]
        fn square_assign_multi(n: &mut Scalar, num_times: usize) {
            for _ in 0..num_times {
                *n = n.square();
            }
        }
        // found using https://github.com/kwantam/addchain for M - 2

        let mut t1 = self.square(); //       1:   2
        let mut t0 = t1 * self; //           2:   3
        let mut t11 = t1.square(); //        3:   4
        let mut t10 = t0.square(); //        4:   6
        let mut t5 = t11.square(); //        5:   8
        let t20 = t10 * t0; //               6:   9
        let mut t2 = t10.square(); //        7:   12
        let t17 = t20 * t11; //              8:   13
        t0 = t5 * t10; //                    9:   14
        t1 = t20 * t10; //                   10:  15
        let t14 = t20 * t5; //               11:  17
        let t13 = t0 * t20; //               12:  23
        let t8 = t14 * t0; //                13:  31
        let t3 = t13 * t0; //                14:  37
        let t4 = t8 * t5; //                 15:  39
        let t18 = t3 * t11; //               16:  41
        let t19 = t8 * t2; //                17:  43
        let t7 = t4 * t5; //                 18:  47
        let t16 = t3 * t2; //                19:  49
        let t9 = t7 * t0; //                 20:  61
        let t6 = t9 * t0; //                 21:  75
        t2 = t6 * t0; //                     22:  89
        t11 = t2 * t11; //                   23:  93
        let t15 = t11 * t0; //               24:  107
        let t12 = t15 * t10; //              25:  113
        t5 = t15 * t5; //                    26:  115
        t10 = t15 * t0; //                   27:  121
        let t21 = t12 * t0; //               28:  127
        t0 = t21.square(); //                29:  254
        square_assign_multi(&mut t0, 6); //  35:  16256
        t0 *= &t21; //                       36:  16383
        square_assign_multi(&mut t0, 7); //  43:  2097024
        t0 *= &t14; //                       44:  2097041
        square_assign_multi(&mut t0, 12); // 56:  8589479936
        t0 *= &t20; //                       57:  8589479945
        square_assign_multi(&mut t0, 11); // 68:  17591254927360
        t0 *= &t16; //                       69:  17591254927409
        square_assign_multi(&mut t0, 7); //  76:  2251680630708352
        t0 *= &t19; //                       77:  2251680630708395
        square_assign_multi(&mut t0, 7); //  84:  288215120730674560
        t0 *= &t14; //                       85:  288215120730674577
        square_assign_multi(&mut t0, 10); // 95:  295132283628210766848
        t0 *= &t15; //                       96:  295132283628210766955
        square_assign_multi(&mut t0, 10); // 106: 302215458435287825361920
        t0 *= &t4; //                        107: 302215458435287825361959
        square_assign_multi(&mut t0, 8); //  115: 77367157359433683292661504
        t0 *= &t18; //                       116: 77367157359433683292661545
        square_assign_multi(&mut t0, 8); //  124: 19805992284015022922921355520
        t0 *= &t17; //                       125: 19805992284015022922921355533
        square_assign_multi(&mut t0, 11); // 136: 40562672197662766946142936131584
        t0 *= &t16; //                       137: 40562672197662766946142936131633
        square_assign_multi(&mut t0, 8); //  145: 10384044082601668338212591649698048
        t0 *= &t15; //                       146: 10384044082601668338212591649698155
        square_assign_multi(&mut t0, 7); //  153: 1329157642573013547291211731161363840
        t0 *= &t14; //                       154: 1329157642573013547291211731161363857
        square_assign_multi(&mut t0, 7); //  161: 170132178249345734053275101588654573696
        t0 *= &t13; //                       162: 170132178249345734053275101588654573719
        square_assign_multi(&mut t0, 6); //  168: 10888459407958126979409606501673892718016
        t0 *= &t8; //                        169: 10888459407958126979409606501673892718047
        square_assign_multi(&mut t0, 11); // 180: 22299564867498244053830874115428132286560256
        t0 *= &t12; //                       181: 22299564867498244053830874115428132286560369
        square_assign_multi(&mut t0, 8); //  189: 5708688606079550477780703773549601865359454464
        t0 *= &t11; //                       190: 5708688606079550477780703773549601865359454557
        square_assign_multi(&mut t0, 7); //  197: 730712141578182461155930083014349038766010183296
        t0 *= &t10; //                       198: 730712141578182461155930083014349038766010183417
        square_assign_multi(&mut t0, 6); //  204: 46765577061003677513979525312918338481024651738688
        t0 *= &t4; //                        205: 46765577061003677513979525312918338481024651738727
        square_assign_multi(&mut t0, 7); //  212: 5985993863808470721789379240053547325571155422557056
        t0 *= &t9; //                        213: 5985993863808470721789379240053547325571155422557117
        square_assign_multi(&mut t0, 7); //  220: 766207214567484252389040542726854057673107894087310976
        t0 *= &t8; //                        221: 766207214567484252389040542726854057673107894087311007
        square_assign_multi(&mut t0, 17); // 238: 100428312027789295929136322016294215047329597893812028309504
        t0 *= &t7; //                        239: 100428312027789295929136322016294215047329597893812028309551
        square_assign_multi(&mut t0, 11); // 250: 205677183032912478062871187489370552416931016486527033977960448
        t0 *= &t6; //                        251: 205677183032912478062871187489370552416931016486527033977960523
        square_assign_multi(&mut t0, 9); //  260: 105306717712851188768190047994557722837468680441101841396715787776
        t0 *= &t5; //                        261: 105306717712851188768190047994557722837468680441101841396715787891
        square_assign_multi(&mut t0, 6); //  267: 6739629933622476081164163071651694261597995548230517849389810425024
        t0 *= &t4; //                        268: 6739629933622476081164163071651694261597995548230517849389810425063
        square_assign_multi(&mut t0, 8); //  276: 1725345263007353876778025746342833730969086860347012569443791468816128
        t0 *= &t3; //                        277: 1725345263007353876778025746342833730969086860347012569443791468816165
        square_assign_multi(&mut t0, 12); // 289: 7067014197278121479282793457020246962049379779981363484441769856271011840
        t0 *= &t2; //                        290: 7067014197278121479282793457020246962049379779981363484441769856271011929
        square_assign_multi(&mut t0, 6); //  296: 452288908625799774674098781249295805571160305918807263004273270801344763456
        t0 *= &t1; //                        297: 452288908625799774674098781249295805571160305918807263004273270801344763471 = M-2

        CtOption::new(t0, !self.ct_eq(&Self::zero()))
    }

    /// Computes a random `Scalar` element
    pub fn random(mut rng: impl RngCore + CryptoRng) -> Self {
        let mut buf = [0; 64];
        rng.fill_bytes(&mut buf);
        Self::from_bytes_wide(&buf)
    }

    #[allow(unused)]
    /// Constructs a `Scalar` element without checking that it is
    /// canonical.
    pub const fn from_raw_unchecked(v: [u64; 4]) -> Self {
        Scalar(v)
    }
}

// OVERLOADED OPERATORS
// ================================================================================================

impl<'a> Neg for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn neg(self) -> Scalar {
        self.neg()
    }
}

impl Neg for Scalar {
    type Output = Scalar;

    #[inline]
    fn neg(self) -> Scalar {
        -&self
    }
}

impl<'a, 'b> Sub<&'b Scalar> for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn sub(self, rhs: &'b Scalar) -> Scalar {
        self.sub(rhs)
    }
}

impl<'a, 'b> Add<&'b Scalar> for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn add(self, rhs: &'b Scalar) -> Scalar {
        self.add(rhs)
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn mul(self, rhs: &'b Scalar) -> Scalar {
        self.mul(rhs)
    }
}

impl_binops_additive!(Scalar, Scalar);
impl_binops_multiplicative!(Scalar, Scalar);

impl Div for Scalar {
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Scalar {
        self.mul(&rhs.invert().unwrap_or(Self::zero()))
    }
}

impl DivAssign for Scalar {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs
    }
}

// TYPE CONVERSIONS
// ================================================================================================

impl From<u128> for Scalar {
    /// Converts a 128-bit value into a field element. If the value is greater than or equal to
    /// the field modulus, modular reduction is silently preformed.
    fn from(value: u128) -> Self {
        let value_high: u64 = (value >> 64).try_into().unwrap();
        let value_low: u64 = (value & (u64::MAX as u128)).try_into().unwrap();
        Scalar::new([value_low, value_high, 0, 0])
    }
}

impl From<u64> for Scalar {
    /// Converts a 64-bit value into a field element.
    fn from(value: u64) -> Self {
        Scalar([value, 0, 0, 0]) * R2
    }
}

impl From<u32> for Scalar {
    /// Converts a 32-bit value into a field element.
    fn from(value: u32) -> Self {
        Scalar([value as u64, 0, 0, 0]) * R2
    }
}

impl From<u16> for Scalar {
    /// Converts a 16-bit value into a field element.
    fn from(value: u16) -> Self {
        Scalar([value as u64, 0, 0, 0]) * R2
    }
}

impl From<u8> for Scalar {
    /// Converts an 8-bit value into a field element.
    fn from(value: u8) -> Self {
        Scalar([value as u64, 0, 0, 0]) * R2
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bitvec::view::AsBits;
    use rand::thread_rng;

    const LARGEST: Scalar = Scalar([
        0x02f096739c941650,
        0x2ef8715df33bd3e0,
        0x46b09ca43418b591,
        0x00fffc8804831564,
    ]);

    // DISPLAY
    // ================================================================================================

    #[test]
    fn test_debug() {
        assert_eq!(
            format!("{:?}", Scalar::zero()),
            "0x0000000000000000000000000000000000000000000000000000000000000000"
        );
        assert_eq!(
            format!("{:?}", Scalar::one()),
            "0x0000000000000000000000000000000000000000000000000000000000000001"
        );
        assert_eq!(
            format!("{:?}", R2),
            "0x000377fb7cea9bb94f635bcbe74a6ed1078ea20cc42c1ffd0f698c636be9af00"
        );
    }

    #[test]
    fn test_to_repr() {
        assert_eq!(format!("{:?}", Scalar::zero().to_repr()), "[0, 0, 0, 0]");
        assert_eq!(format!("{:?}", Scalar::one().to_repr()), "[1, 0, 0, 0]");
        assert_eq!(
            format!("{:?}", R2.to_repr()),
            "[1110573141763665664, 544550780672942077, 5720516883007565521, 976346946378681]"
        );
    }

    // BASIC ALGEBRA
    // ================================================================================================

    #[test]
    fn test_equality() {
        assert_eq!(Scalar::default(), Scalar::zero());
        assert_eq!(Scalar::zero(), Scalar::zero());
        assert_eq!(Scalar::one(), Scalar::one());

        assert!(bool::from(Scalar::default().is_zero()));
        assert!(!bool::from(Scalar::zero().ct_eq(&Scalar::one())));

        assert!(Scalar::zero() != Scalar::one());
        assert!(Scalar::one() != R2);
    }

    #[test]
    fn test_addition() {
        let mut tmp = LARGEST;
        tmp += &LARGEST;

        assert_eq!(
            tmp,
            Scalar([
                0x02f096739c94164f,
                0x2ef8715df33bd3e0,
                0x46b09ca43418b591,
                0x00fffc8804831564,
            ])
        );

        assert_eq!(tmp, LARGEST.double());

        let mut tmp = LARGEST;
        tmp += &Scalar([1, 0, 0, 0]);

        assert_eq!(tmp, Scalar::zero());
    }

    #[test]
    fn test_subtraction() {
        let mut tmp = LARGEST;
        tmp -= &LARGEST;

        assert_eq!(tmp, Scalar::zero());

        let mut tmp = Scalar::zero();
        tmp -= &LARGEST;

        let mut tmp2 = M;
        tmp2 -= &LARGEST;

        assert_eq!(tmp, tmp2);
    }

    #[test]
    fn test_negation() {
        let tmp = -&LARGEST;

        assert_eq!(tmp, Scalar([1, 0, 0, 0]));

        let tmp = -&Scalar::zero();
        assert_eq!(tmp, Scalar::zero());
        let tmp = -&Scalar([1, 0, 0, 0]);
        assert_eq!(tmp, LARGEST);
    }

    #[test]
    fn test_multiplication() {
        let mut cur = LARGEST;

        for _ in 0..100 {
            let mut tmp = cur;
            tmp *= &cur;

            let mut tmp2 = Scalar::zero();
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
        assert!(bool::from(Scalar::zero().invert().is_none()));
        assert_eq!(Scalar::one().invert().unwrap(), Scalar::one());
        assert_eq!((-&Scalar::one()).invert().unwrap(), -&Scalar::one());

        let mut tmp = Scalar::random(&mut thread_rng());

        for _ in 0..100 {
            let mut tmp2 = tmp.invert().unwrap();
            tmp2.mul_assign(&tmp);

            assert_eq!(tmp2, Scalar::one());

            tmp.add_assign(&Scalar::random(&mut thread_rng()));
        }
    }

    #[test]
    fn test_squaring() {
        let mut cur = LARGEST;

        for _ in 0..100 {
            let mut tmp = cur;
            let pow2 = tmp.exp(&[2, 0, 0, 0]);
            tmp = tmp.square();

            let mut tmp2 = Scalar::zero();
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
        let q_minus_2 = [
            0x02f096739c94164f,
            0x2ef8715df33bd3e0,
            0x46b09ca43418b591,
            0x00fffc8804831564,
        ];

        let mut r1 = R;
        let mut r2 = R;
        let mut r3 = R;

        for _ in 0..100 {
            r1 = r1.invert().unwrap();
            r2 = r2.exp(&q_minus_2);
            r3 = r3.exp_vartime(&q_minus_2);

            assert_eq!(r1, r2);
            assert_eq!(r2, r3);
            // Add R so we check something different next time around
            r1.add_assign(&R);
            r2 = r1;
            r3 = r1;
        }
    }

    // SERIALIZATION / DESERIALIZATION
    // ================================================================================================

    #[test]
    fn test_to_bytes() {
        assert_eq!(
            Scalar::zero().to_bytes(),
            [
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ]
        );

        assert_eq!(
            Scalar::one().to_bytes(),
            [
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ]
        );

        assert_eq!(
            R2.to_bytes(),
            [
                0, 175, 233, 107, 99, 140, 105, 15, 253, 31, 44, 196, 12, 162, 142, 7, 209, 110,
                74, 231, 203, 91, 99, 79, 185, 155, 234, 124, 251, 119, 3, 0
            ]
        );

        assert_eq!(
            (-&Scalar::one()).to_bytes(),
            [
                80, 22, 148, 156, 115, 150, 240, 2, 224, 211, 59, 243, 93, 113, 248, 46, 145, 181,
                24, 52, 164, 156, 176, 70, 100, 21, 131, 4, 136, 252, 255, 0
            ]
        );
    }

    #[test]
    fn test_from_bytes() {
        assert_eq!(
            Scalar::from_bytes(&[
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ])
            .unwrap(),
            Scalar::zero()
        );

        assert_eq!(
            Scalar::from_bytes(&[
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ])
            .unwrap(),
            Scalar::one()
        );

        assert_eq!(
            Scalar::from_bytes(&[
                0, 175, 233, 107, 99, 140, 105, 15, 253, 31, 44, 196, 12, 162, 142, 7, 209, 110,
                74, 231, 203, 91, 99, 79, 185, 155, 234, 124, 251, 119, 3, 0
            ])
            .unwrap(),
            R2
        );

        // -1 should work
        assert_eq!(
            Scalar::from_bytes(&[
                80, 22, 148, 156, 115, 150, 240, 2, 224, 211, 59, 243, 93, 113, 248, 46, 145, 181,
                24, 52, 164, 156, 176, 70, 100, 21, 131, 4, 136, 252, 255, 0
            ])
            .unwrap(),
            -Scalar::one(),
        );

        // M is invalid
        assert!(bool::from(
            Scalar::from_bytes(&[
                81, 22, 148, 156, 115, 150, 240, 2, 224, 211, 59, 243, 93, 113, 248, 46, 145, 181,
                24, 52, 164, 156, 176, 70, 100, 21, 131, 4, 136, 252, 255, 0
            ])
            .is_none()
        ));

        // Anything larger than the M is invalid
        assert!(bool::from(
            Scalar::from_bytes(&[
                82, 22, 148, 156, 115, 150, 240, 2, 224, 211, 59, 243, 93, 113, 248, 46, 145, 181,
                24, 52, 164, 156, 176, 70, 100, 21, 131, 4, 136, 252, 255, 0
            ])
            .is_none()
        ));
        assert!(bool::from(
            Scalar::from_bytes(&[
                81, 22, 148, 156, 115, 150, 240, 2, 224, 211, 59, 243, 93, 113, 248, 46, 145, 181,
                24, 52, 164, 156, 176, 70, 100, 21, 131, 4, 136, 252, 255, 1
            ])
            .is_none()
        ));
    }

    #[test]
    fn test_from_bits_vartime() {
        let bytes = Scalar::zero().to_bytes();
        assert_eq!(
            Scalar::from_bits_vartime(bytes.as_bits::<Lsb0>()),
            Scalar::zero()
        );

        let bytes = Scalar::one().to_bytes();
        assert_eq!(
            Scalar::from_bits_vartime(bytes.as_bits::<Lsb0>()),
            Scalar::one()
        );

        let bytes = R2.to_bytes();
        assert_eq!(Scalar::from_bits_vartime(bytes.as_bits::<Lsb0>()), R2);

        // -1 should work
        let bytes = (-Scalar::one()).to_bytes();
        assert_eq!(
            Scalar::from_bits_vartime(bytes.as_bits::<Lsb0>()),
            -Scalar::one()
        );

        // Modulus results in Scalar::zero()
        let bytes = [
            81, 22, 148, 156, 115, 150, 240, 2, 224, 211, 59, 243, 93, 113, 248, 46, 145, 181, 24,
            52, 164, 156, 176, 70, 100, 21, 131, 4, 136, 252, 255, 0,
        ];
        assert_eq!(
            Scalar::from_bits_vartime(bytes.as_bits::<Lsb0>()),
            Scalar::zero()
        );
    }

    #[test]
    fn test_from_u512_zero() {
        assert_eq!(
            Scalar::zero(),
            Scalar::from_u512([M.0[0], M.0[1], M.0[2], M.0[3], 0, 0, 0, 0])
        );
    }

    #[test]
    fn test_from_u512_r() {
        assert_eq!(R, Scalar::from_u512([1, 0, 0, 0, 0, 0, 0, 0]));
    }

    #[test]
    fn test_from_u512_r2() {
        assert_eq!(R2, Scalar::from_u512([0, 0, 0, 0, 1, 0, 0, 0]));
    }

    #[test]
    fn test_from_u512_max() {
        let max_u64 = 0xffff_ffff_ffff_ffff;
        assert_eq!(
            R3 - R,
            Scalar::from_u512([
                max_u64, max_u64, max_u64, max_u64, max_u64, max_u64, max_u64, max_u64
            ])
        );
    }

    #[test]
    fn test_from_bytes_wide_r2() {
        assert_eq!(
            R2,
            Scalar::from_bytes_wide(&[
                0, 175, 233, 107, 99, 140, 105, 15, 253, 31, 44, 196, 12, 162, 142, 7, 209, 110,
                74, 231, 203, 91, 99, 79, 185, 155, 234, 124, 251, 119, 3, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            ])
        );
    }

    #[test]
    fn test_from_bytes_wide_negative_one() {
        assert_eq!(
            -&Scalar::one(),
            Scalar::from_bytes_wide(&[
                80, 22, 148, 156, 115, 150, 240, 2, 224, 211, 59, 243, 93, 113, 248, 46, 145, 181,
                24, 52, 164, 156, 176, 70, 100, 21, 131, 4, 136, 252, 255, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            ])
        );
    }

    #[test]
    fn test_from_bytes_wide_maximum() {
        assert_eq!(
            Scalar([
                0xe1ad557bb165ca94,
                0x710642caf19db439,
                0x5f769ac31b1073d6,
                0x0011c8187f7c1b0d,
            ]),
            Scalar::from_bytes_wide(&[0xff; 64])
        );
    }

    #[test]
    fn test_lexicographically_largest() {
        // a = 103745092170913632337395666750586378382978886801660409715261768932739876504
        let a = Scalar([
            0x0bdea22ae33b8697,
            0x05aff6ecd268cc65,
            0x9400646c57a18ebb,
            0x00d1058d04b33452,
        ]);

        // b = 348543816454886142336703114498709427188181419117146853289011501868604886969
        let b = Scalar([
            0xf711f448b9588fba,
            0x29487a7120d3077a,
            0xb2b03837dc7726d6,
            0x002ef6faffcfe111,
        ]);

        assert_eq!(a.square(), b.square());
        assert!(!bool::from(a.lexicographically_largest()));
        assert!(bool::from(b.lexicographically_largest()));
    }

    #[test]
    fn test_zeroize() {
        use zeroize::Zeroize;

        let mut a = Scalar::one();
        a.zeroize();
        assert_eq!(a, Scalar::zero());
    }

    // INITIALIZATION
    // ================================================================================================

    #[test]
    fn test_from_int() {
        let n = 42u8;
        let element = Scalar::from(n);

        assert_eq!(element, Scalar::from(n as u16));
        assert_eq!(element, Scalar::from(n as u32));
        assert_eq!(element, Scalar::from(n as u64));
        assert_eq!(element, Scalar::from(n as u128));
    }

    // SERDE SERIALIZATIOIN
    // ================================================================================================

    #[test]
    #[cfg(feature = "serialize")]
    fn test_serde_scalar() {
        let mut rng = thread_rng();
        let scalar = Scalar::random(&mut rng);
        let encoded = bincode::serialize(&scalar).unwrap();
        let parsed: Scalar = bincode::deserialize(&encoded).unwrap();
        assert_eq!(parsed, scalar);

        // Check that the encoding is 32 bytes exactly
        assert_eq!(encoded.len(), 32);

        // Check that the encoding itself matches the usual one
        assert_eq!(scalar, bincode::deserialize(&scalar.to_bytes()).unwrap());

        // Check that invalid encodings fail
        let scalar = Scalar::random(&mut rng);
        let mut encoded = bincode::serialize(&scalar).unwrap();
        encoded[31] = 127;
        assert!(bincode::deserialize::<Scalar>(&encoded).is_err());

        let encoded = bincode::serialize(&scalar).unwrap();
        assert!(bincode::deserialize::<Scalar>(&encoded[0..31]).is_err());
    }
}
