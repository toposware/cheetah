use core::{
    convert::{TryFrom, TryInto},
    fmt::{self, Debug, Display, Formatter},
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::utils::{add64_with_carry, mul64_with_carry, square_assign_multi, sub64_with_carry};

use bitvec::{order::Lsb0, slice::BitSlice};
use group::ff::{Field, PrimeField};
use rand_core::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "serialize")]
use serde::de::Visitor;
#[cfg(feature = "serialize")]
use serde::{self, Deserialize, Deserializer, Serialize, Serializer};

// CONSTANTS
// ================================================================================================

// Field modulus = 19059060964990286085476939486711779813624322081295683206882804078282896682513
const M: Scalar = Scalar([
    0xfe94da7e25f1ae11,
    0x4880a66664b93701,
    0xeb604b4894046cd6,
    0x2a230bd593a1a279,
]);

// 2^256 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R: Scalar = Scalar([
    0x0882e10b1c55eb9a,
    0x4cfc1999a3a8b5f4,
    0x7bbe3c4c87e572fa,
    0x032db8fe8a363124,
]);

// 2^512 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R2: Scalar = Scalar([
    0xa2fad5b0ff7248f5,
    0x50d020a7612c2ebd,
    0x8528734b9835d288,
    0x165ef986d0bcb87e,
]);

// 2^768 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R3: Scalar = Scalar([
    0xe9d7977686156655,
    0x473cd715cf34db3f,
    0xd8d14d270e51010f,
    0x03e4cc1ae7863573,
]);

// Multiplicative generator g of order q-1
// g = 5
//   = 0xfe49cf8b30ef5b66ab72d7ea77b3ee380ec8000324b8dc42a8e65378dad9a02 in Montgomery form
const GENERATOR: Scalar = Scalar([
    0x2a8e65378dad9a02,
    0x80ec8000324b8dc4,
    0x6ab72d7ea77b3ee3,
    0x0fe49cf8b30ef5b6,
]);

// Two-adicity of the field: (q-1) % 2^4 = 0
const TWO_ADICITY: u32 = 4;

// 2^4 root of unity = 0xb89c7221d6ffc2bc0f548e825b15c6337727a16c3dca7b87df3bedcfea25fee
//                   = 0xaaa534ded054530efa0e08b631ae066ca4a0afe638ae93c21d02c60c8519ff0 in Montgomery form
const TWO_ADIC_ROOT_OF_UNITY: Scalar = Scalar([
    0x21d02c60c8519ff0,
    0xca4a0afe638ae93c,
    0xefa0e08b631ae066,
    0x0aaa534ded054530,
]);

// -M^{-1} mod 2^64; this is used during element multiplication.
const U: u64 = 11876065156671798543;

// SCALAR FIELD ELEMENT
// ================================================================================================

/// Represents a scalar field element.
///
/// Internal values are stored in Montgomery representation.
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

    /// Computes the square root of this element, if it exists.
    pub fn sqrt(&self) -> CtOption<Self> {
        // Tonelli-Shank's algorithm for q mod 16 = 1

        // w = self^((t - 1) // 2)
        //   = self^0x151185eac9d0d13cf5b025a44a02366b24405333325c9b80ff4a6d3f12f8d70
        let w = self.exp_vartime(&[
            0xfa5369f897c6b844,
            0x2202999992e4dc07,
            0xad812d225011b359,
            0x208c2f564e8689e7,
        ]);

        let mut v = TWO_ADICITY;
        let mut x = self * w;
        let mut b = x * w;

        // Initialize z as the 2^S root of unity.
        let mut z = TWO_ADIC_ROOT_OF_UNITY;

        for max_v in (1..=TWO_ADICITY).rev() {
            let mut k = 1;
            let mut tmp = b.square();
            let mut j_less_than_v: Choice = 1.into();

            for j in 2..max_v {
                let tmp_is_one = tmp.ct_eq(&Scalar::one());
                let squared = Scalar::conditional_select(&tmp, &z, tmp_is_one).square();
                tmp = Scalar::conditional_select(&squared, &tmp, tmp_is_one);

                let new_z = Scalar::conditional_select(&z, &squared, tmp_is_one);
                j_less_than_v &= !j.ct_eq(&v);
                k = u32::conditional_select(&j, &k, tmp_is_one);
                z = Scalar::conditional_select(&z, &new_z, j_less_than_v);
            }

            let result = x * z;
            x = Scalar::conditional_select(&result, &x, b.ct_eq(&Scalar::one()));
            z = z.square();
            b *= z;
            v = k;
        }

        CtOption::new(x, (x * x).ct_eq(self))
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
    pub const fn output_limbs(&self) -> [u64; 4] {
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

        let (_, borrow) = sub64_with_carry(tmp.0[0], 0xff4a_6d3f_12f8_d708, 0);
        let (_, borrow) = sub64_with_carry(tmp.0[1], 0x2440_5333_325c_9b80, borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[2], 0xf5b0_25a4_4a02_366b, borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[3], 0x1511_85ea_c9d0_d13c, borrow);

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
    pub fn invert(&self) -> CtOption<Self> {
        // found using https://github.com/kwantam/addchain for M - 2
        let mut t0 = self.square(); //       1:   2
        let t15 = t0 * self; //              2:   3
        let t8 = t15 * t0; //                3:   5
        let t2 = t8 * t0; //                 4:   7
        let t5 = t2 * t0; //                 5:   9
        let t13 = t5 * t0; //                6:   11
        let t3 = t13 * t0; //                7:   13
        let t1 = t3 * t0; //                 8:   15
        let t12 = t1 * t0; //                9:   17
        let t6 = t12 * t0; //                10:  19
        let t14 = t6 * t0; //                11:  21
        let t10 = t14 * t0; //               12:  23
        let t11 = t10 * t0; //               13:  25
        let t7 = t11 * t0; //                14:  27
        let t9 = t7 * t0; //                 15:  29
        let t4 = t9 * t0; //                 16:  31
        t0 = t14.square(); //                17:  42
        square_assign_multi(&mut t0, 3); //  20:  336
        t0 *= self; //                       21:  337
        square_assign_multi(&mut t0, 5); //  26:  10784
        t0 *= t15; //                        27:  10787
        square_assign_multi(&mut t0, 9); //  36:  5522944
        t0 *= t10; //                        37:  5522967
        square_assign_multi(&mut t0, 5); //  42:  176734944
        t0 *= t14; //                        43:  176734965
        square_assign_multi(&mut t0, 6); //  49:  11311037760
        t0 *= t11; //                        50:  11311037785
        square_assign_multi(&mut t0, 7); //  57:  1447812836480
        t0 *= t9; //                         58:  1447812836509
        square_assign_multi(&mut t0, 8); //  66:  370640086146304
        t0 *= t3; //                         67:  370640086146317
        square_assign_multi(&mut t0, 8); //  75:  94883862053457152
        t0 *= t6; //                         76:  94883862053457171
        square_assign_multi(&mut t0, 5); //  81:  3036283585710629472
        t0 *= t11; //                        82:  3036283585710629497
        square_assign_multi(&mut t0, 5); //  87:  97161074742740143904
        t0 *= t9; //                         88:  97161074742740143933
        square_assign_multi(&mut t0, 6); //  94:  6218308783535369211712
        t0 *= t7; //                         95:  6218308783535369211739
        square_assign_multi(&mut t0, 10); // 105: 6367548194340218072820736
        t0 *= t5; //                         106: 6367548194340218072820745
        square_assign_multi(&mut t0, 5); //  111: 203761542218886978330263840
        t0 *= t3; //                         112: 203761542218886978330263853
        square_assign_multi(&mut t0, 7); //  119: 26081477404017533226273773184
        t0 *= t12; //                        120: 26081477404017533226273773201
        square_assign_multi(&mut t0, 5); //  125: 834607276928561063240760742432
        t0 *= t8; //                         126: 834607276928561063240760742437
        square_assign_multi(&mut t0, 12); // 138: 3418551406299386115034156001021952
        t0 *= t12; //                        139: 3418551406299386115034156001021969
        square_assign_multi(&mut t0, 4); //  143: 54696822500790177840546496016351504
        t0 *= t13; //                        144: 54696822500790177840546496016351515
        square_assign_multi(&mut t0, 6); //  150: 3500596640050571381794975745046496960
        t0 *= t3; //                         151: 3500596640050571381794975745046496973
        square_assign_multi(&mut t0, 6); //  157: 224038184963236568434878447682975806272
        t0 *= t11; //                        158: 224038184963236568434878447682975806297
        square_assign_multi(&mut t0, 7); //  165: 28676887675294280759664441303420903206016
        t0 *= t12; //                        166: 28676887675294280759664441303420903206033
        square_assign_multi(&mut t0, 10); // 176: 29365132979501343497896387894703004882977792
        t0 *= t8; //                         177: 29365132979501343497896387894703004882977797
        square_assign_multi(&mut t0, 7); //  184: 3758737021376171967730737650521984625021158016
        t0 *= t11; //                        185: 3758737021376171967730737650521984625021158041
        square_assign_multi(&mut t0, 5); //  190: 120279584684037502967383604816703508000677057312
        t0 *= t6; //                         191: 120279584684037502967383604816703508000677057331
        square_assign_multi(&mut t0, 7); //  198: 15395786839556800379825101416538049024086663338368
        t0 *= t11; //                        199: 15395786839556800379825101416538049024086663338393
        square_assign_multi(&mut t0, 7); //  206: 1970660715463270448617612981316870275083092907314304
        t0 *= t10; //                        207: 1970660715463270448617612981316870275083092907314327
        square_assign_multi(&mut t0, 7); //  214: 252244571579298617423054461608559395210635892136233856
        t0 *= t6; //                         215: 252244571579298617423054461608559395210635892136233875
        square_assign_multi(&mut t0, 4); //  219: 4035913145268777878768871385736950323370174274179742000
        t0 *= t2; //                         220: 4035913145268777878768871385736950323370174274179742007
        square_assign_multi(&mut t0, 12); // 232: 16531100243020914191437297195978548524524233827040223260672
        t0 *= t4; //                         233: 16531100243020914191437297195978548524524233827040223260703
        square_assign_multi(&mut t0, 5); //  238: 528995207776669254125993510271313552784775482465287144342496
        t0 *= t9; //                         239: 528995207776669254125993510271313552784775482465287144342525
        square_assign_multi(&mut t0, 5); //  244: 16927846648853416132031792328682033689112815438889188618960800
        t0 *= t8; //                         245: 16927846648853416132031792328682033689112815438889188618960805
        square_assign_multi(&mut t0, 7); //  252: 2166764371053237264900069418071300312206440376177816143226983040
        t0 *= t7; //                         253: 2166764371053237264900069418071300312206440376177816143226983067
        square_assign_multi(&mut t0, 6); //  259: 138672919747407184953604442756563219981212184075380233166526916288
        t0 *= t6; //                         260: 138672919747407184953604442756563219981212184075380233166526916307
        square_assign_multi(&mut t0, 4); //  264: 2218766715958514959257671084105011519699394945206083730664430660912
        t0 *= t1; //                         265: 2218766715958514959257671084105011519699394945206083730664430660927
        square_assign_multi(&mut t0, 7); //  272: 284002139642689914784981898765441474521522552986378717525047124598656
        t0 *= t5; //                         273: 284002139642689914784981898765441474521522552986378717525047124598665
        square_assign_multi(&mut t0, 6); //  279: 18176136937132154546238841520988254369377443391128237921603015974314560
        t0 *= t4; //                         280: 18176136937132154546238841520988254369377443391128237921603015974314591
        square_assign_multi(&mut t0, 7); //  287: 2326545527952915781918571714686496559280312754064414453965186044712267648
        t0 *= t3; //                         288: 2326545527952915781918571714686496559280312754064414453965186044712267661
        square_assign_multi(&mut t0, 4); //  292: 37224728447246652510697147434983944948485004065030631263442976715396282576
        t0 *= t2; //                         293: 37224728447246652510697147434983944948485004065030631263442976715396282583
        square_assign_multi(&mut t0, 9); //  302: 19059060964990286085476939486711779813624322081295683206882804078282896682496
        t0 *= t1; //                         303: 19059060964990286085476939486711779813624322081295683206882804078282896682511 = M - 2

        CtOption::new(t0, !self.ct_eq(&Self::zero()))
    }

    /// Computes the conjugate of a `Scalar` element
    pub fn conjugate(&self) -> Self {
        Scalar(self.0)
    }

    #[inline(always)]
    /// Normalizes the internal representation of an `Fp` element
    pub fn normalize(&mut self) {
        *self *= &R2
    }

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

// FIELD TRAITS IMPLEMENTATION
// ================================================================================================

impl Field for Scalar {
    fn random(mut rng: impl RngCore) -> Self {
        let mut buf = [0; 64];
        rng.fill_bytes(&mut buf);
        Self::from_bytes_wide(&buf)
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

impl PrimeField for Scalar {
    type Repr = [u8; 32];

    fn from_repr(r: Self::Repr) -> CtOption<Self> {
        Self::from_bytes(&r)
    }

    fn to_repr(&self) -> Self::Repr {
        self.to_bytes()
    }

    fn is_odd(&self) -> Choice {
        (self.to_bytes()[0] & 1).ct_eq(&1)
    }

    const NUM_BITS: u32 = 254;
    const CAPACITY: u32 = Self::NUM_BITS - 1;

    fn multiplicative_generator() -> Self {
        GENERATOR
    }

    const S: u32 = TWO_ADICITY;

    fn root_of_unity() -> Self {
        TWO_ADIC_ROOT_OF_UNITY
    }
}

// SERDE SERIALIZATION
// ================================================================================================

#[cfg(feature = "serialize")]
impl Serialize for Scalar {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        use serde::ser::SerializeTuple;
        let mut tup = serializer.serialize_tuple(32)?;
        for byte in self.to_bytes().iter() {
            tup.serialize_element(byte)?;
        }
        tup.end()
    }
}

#[cfg(feature = "serialize")]
impl<'de> Deserialize<'de> for Scalar {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct ScalarVisitor;

        impl<'de> Visitor<'de> for ScalarVisitor {
            type Value = Scalar;

            fn expecting(&self, formatter: &mut Formatter) -> fmt::Result {
                formatter.write_str("a valid field element")
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<Scalar, A::Error>
            where
                A: serde::de::SeqAccess<'de>,
            {
                let mut bytes = [0u8; 32];
                for (i, byte) in bytes.iter_mut().enumerate() {
                    *byte = seq
                        .next_element()?
                        .ok_or_else(|| serde::de::Error::invalid_length(i, &"expected 32 bytes"))?;
                }
                let elem = Scalar::from_bytes(&bytes);
                if bool::from(elem.is_none()) {
                    Err(serde::de::Error::custom("decompression failed"))
                } else {
                    Ok(elem.unwrap())
                }
            }
        }

        deserializer.deserialize_tuple(32, ScalarVisitor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bitvec::view::AsBits;
    use rand::thread_rng;

    const LARGEST: Scalar = Scalar([
        0xfe94da7e25f1ae10,
        0x4880a66664b93701,
        0xeb604b4894046cd6,
        0x2a230bd593a1a279,
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
            "0x032db8fe8a3631247bbe3c4c87e572fa4cfc1999a3a8b5f40882e10b1c55eb9a"
        );
    }

    #[test]
    fn test_output_limbs() {
        assert_eq!(
            format!("{:?}", Scalar::zero().output_limbs()),
            "[0, 0, 0, 0]"
        );
        assert_eq!(
            format!("{:?}", Scalar::one().output_limbs()),
            "[1, 0, 0, 0]"
        );
        assert_eq!(
            format!("{:?}", R2.output_limbs()),
            "[613299937112091546, 5547336988680041972, 8916630611635303162, 229042559445774628]"
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
                0xfe94da7e25f1ae0f,
                0x4880a66664b93701,
                0xeb604b4894046cd6,
                0x2a230bd593a1a279,
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
        let mut rng = thread_rng();

        let tmp = -&LARGEST;

        assert_eq!(tmp, Scalar([1, 0, 0, 0]));

        let tmp = -&Scalar::zero();
        assert_eq!(tmp, Scalar::zero());
        let tmp = -&Scalar([1, 0, 0, 0]);
        assert_eq!(tmp, LARGEST);

        for _ in 0..100 {
            let a = Scalar::random(&mut rng);
            let b = -a;

            assert_eq!(a + b, Scalar::zero());
        }
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
            0xfe94da7e25f1ae0f,
            0x4880a66664b93701,
            0xeb604b4894046cd6,
            0x2a230bd593a1a279,
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

    #[test]
    fn test_conjugate() {
        let a = Scalar::random(&mut thread_rng());
        let b = a.conjugate();
        assert_eq!(a, b);
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
                154, 235, 85, 28, 11, 225, 130, 8, 244, 181, 168, 163, 153, 25, 252, 76, 250, 114,
                229, 135, 76, 60, 190, 123, 36, 49, 54, 138, 254, 184, 45, 3
            ]
        );

        assert_eq!(
            (-&Scalar::one()).to_bytes(),
            [
                16, 174, 241, 37, 126, 218, 148, 254, 1, 55, 185, 100, 102, 166, 128, 72, 214, 108,
                4, 148, 72, 75, 96, 235, 121, 162, 161, 147, 213, 11, 35, 42
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
                154, 235, 85, 28, 11, 225, 130, 8, 244, 181, 168, 163, 153, 25, 252, 76, 250, 114,
                229, 135, 76, 60, 190, 123, 36, 49, 54, 138, 254, 184, 45, 3
            ])
            .unwrap(),
            R2
        );

        // -1 should work
        assert_eq!(
            Scalar::from_bytes(&[
                16, 174, 241, 37, 126, 218, 148, 254, 1, 55, 185, 100, 102, 166, 128, 72, 214, 108,
                4, 148, 72, 75, 96, 235, 121, 162, 161, 147, 213, 11, 35, 42
            ])
            .unwrap(),
            -Scalar::one(),
        );

        // M is invalid
        assert!(bool::from(
            Scalar::from_bytes(&[
                17, 174, 241, 37, 126, 218, 148, 254, 1, 55, 185, 100, 102, 166, 128, 72, 214, 108,
                4, 148, 72, 75, 96, 235, 121, 162, 161, 147, 213, 11, 35, 42
            ])
            .is_none()
        ));

        // Anything larger than the M is invalid
        assert!(bool::from(
            Scalar::from_bytes(&[
                18, 174, 241, 37, 126, 218, 148, 254, 1, 55, 185, 100, 102, 166, 128, 72, 214, 108,
                4, 148, 72, 75, 96, 235, 121, 162, 161, 147, 213, 11, 35, 42
            ])
            .is_none()
        ));
        assert!(bool::from(
            Scalar::from_bytes(&[
                15, 174, 241, 37, 126, 218, 148, 254, 2, 55, 185, 100, 102, 166, 128, 72, 214, 108,
                4, 148, 72, 75, 96, 235, 121, 162, 161, 147, 213, 11, 35, 42
            ])
            .is_none()
        ));
        assert!(bool::from(
            Scalar::from_bytes(&[
                15, 174, 241, 37, 126, 218, 148, 254, 1, 55, 185, 100, 102, 166, 128, 72, 214, 108,
                4, 148, 72, 75, 96, 235, 121, 162, 161, 147, 213, 11, 35, 255
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
            17, 174, 241, 37, 126, 218, 148, 254, 1, 55, 185, 100, 102, 166, 128, 72, 214, 108, 4,
            148, 72, 75, 96, 235, 121, 162, 161, 147, 213, 11, 35, 42,
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
                154, 235, 85, 28, 11, 225, 130, 8, 244, 181, 168, 163, 153, 25, 252, 76, 250, 114,
                229, 135, 76, 60, 190, 123, 36, 49, 54, 138, 254, 184, 45, 3, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
        );
    }

    #[test]
    fn test_from_bytes_wide_negative_one() {
        assert_eq!(
            -&Scalar::one(),
            Scalar::from_bytes_wide(&[
                16, 174, 241, 37, 126, 218, 148, 254, 1, 55, 185, 100, 102, 166, 128, 72, 214, 108,
                4, 148, 72, 75, 96, 235, 121, 162, 161, 147, 213, 11, 35, 42, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
        );
    }

    #[test]
    fn test_from_bytes_wide_maximum() {
        assert_eq!(
            Scalar([
                0xe154b66b69bf7abb,
                0xfa40bd7c2b8c254b,
                0x5d1310da866b8e14,
                0x00b7131c5d50044f,
            ]),
            Scalar::from_bytes_wide(&[0xff; 64])
        );
    }

    #[test]
    fn test_lexicographically_largest() {
        // a = 3194249248449702161625907824974691224153811289682184296505134093839210143169
        let a = Scalar([
            0x1a27928f89f45dc1,
            0xb0a11736719d96b6,
            0x0a4011f33f25c7af,
            0x070fe189977f18a3,
        ]);

        // b = 15864811716540583923851031661737088589470510791613498910377669984443686539344
        let b = Scalar([
            0xe46d47ee9bfd5050,
            0x97df8f2ff31ba04b,
            0xe120395554dea526,
            0x23132a4bfc2289d6,
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

    #[test]
    fn test_from_raw_unchecked() {
        let mut element = Scalar::from_raw_unchecked([
            613299937112091546,
            5547336988680041972,
            8916630611635303162,
            229042559445774628,
        ]);

        let element_normalized = Scalar::new([
            613299937112091546,
            5547336988680041972,
            8916630611635303162,
            229042559445774628,
        ]);

        assert_eq!(element, Scalar::one());
        element.normalize();

        assert!(element != Scalar::one());
        assert_eq!(element, element_normalized);
    }

    // SERDE SERIALIZATIOIN
    // ================================================================================================

    #[test]
    #[cfg(feature = "serialize")]
    fn test_serde_field() {
        let mut rng = thread_rng();
        let element = Scalar::random(&mut rng);
        let encoded = bincode::serialize(&element).unwrap();
        let parsed: Scalar = bincode::deserialize(&encoded).unwrap();
        assert_eq!(parsed, element);

        // Check that the encoding is 32 bytes exactly
        assert_eq!(encoded.len(), 32);

        // Check that the encoding itself matches the usual one
        assert_eq!(element, bincode::deserialize(&element.to_bytes()).unwrap());

        // Check that invalid encodings fail
        let element = Scalar::random(&mut rng);
        let mut encoded = bincode::serialize(&element).unwrap();
        encoded[31] = 127;
        assert!(bincode::deserialize::<Scalar>(&encoded).is_err());

        let encoded = bincode::serialize(&element).unwrap();
        assert!(bincode::deserialize::<Scalar>(&encoded[0..31]).is_err());
    }
}
