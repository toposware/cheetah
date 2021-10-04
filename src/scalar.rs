use core::{
    convert::{TryFrom, TryInto},
    fmt::{self, Debug, Display, Formatter},
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::utils::{add64_with_carry, mul64_with_carry, sub64_with_carry};

use bitvec::{order::Lsb0, slice::BitSlice};
use rand_core::{CryptoRng, RngCore};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

// CONSTANTS
// ================================================================================================

// Field modulus = 9618866709266846897190681717726035074731103115011964004359729618342275068175251479314708351071139295228317576333
const M: Scalar = Scalar([
    0x6210a856e8e2ac8d,
    0x2774ffd1eeca0208,
    0x833c7cca5f5735d0,
    0x72a71e3a69c7b80f,
    0x52f1b32fc5d4ddf9,
    0x000fffacc0b47aef,
]);

/// 2^384 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R: Scalar = Scalar([
    0xf57a9171d5373000,
    0xb002e1135fdf79de,
    0x38335a0a8ca2fd88,
    0x8e1c5963847f07cc,
    0xe4cd03a2b22068d5,
    0x000533f4b8510ad0,
]);

/// 2^768 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R2: Scalar = Scalar([
    0xb5abf62b949211ef,
    0xf6fe56a92fd697c6,
    0x97b8009b9c1a408b,
    0x964d02bf42d11f0d,
    0xcf9bd5ee5563ea90,
    0x0003cb51dd9af85b,
]);

/// 2^1152 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R3: Scalar = Scalar([
    0xbffff69417dd9d74,
    0x400dbe34b209590a,
    0x6339965a283a9add,
    0x0082f3ee34e8f243,
    0x3bc0c639696414ad,
    0x000a96be6cadda40,
]);

/// -M^{-1} mod 2^64; this is used during element multiplication.
const U: u64 = 17813468779381197243;

// SCALAR FIELD ELEMENT
// ================================================================================================

/// Represents a scalar field element.
///
/// Internal values are stored in their canonical form in the range [0, M).
/// The backing type is `[u64; 6]`.
#[derive(Copy, Clone, Eq)]
pub struct Scalar(pub(crate) [u64; 6]);

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
            & self.0[4].ct_eq(&other.0[4])
            & self.0[5].ct_eq(&other.0[5])
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
            u64::conditional_select(&a.0[4], &b.0[4], choice),
            u64::conditional_select(&a.0[5], &b.0[5], choice),
        ])
    }
}

impl zeroize::DefaultIsZeroes for Scalar {}

impl Scalar {
    /// Creates a new field element from a [u64; 6] value.
    /// The value is converted to Montgomery form by computing
    /// (a.R^0 * R^2) / R = a.R
    pub const fn new(value: [u64; 6]) -> Self {
        (&Scalar(value)).mul(&R2)
    }

    /// Returns zero, the additive identity.
    #[inline]
    pub const fn zero() -> Self {
        Scalar([0, 0, 0, 0, 0, 0])
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
        r8: u64,
        r9: u64,
        r10: u64,
        r11: u64,
    ) -> Self {
        let k = r0.wrapping_mul(U);
        let (_, carry) = mul64_with_carry(r0, k, M.0[0], 0);
        let (r1, carry) = mul64_with_carry(r1, k, M.0[1], carry);
        let (r2, carry) = mul64_with_carry(r2, k, M.0[2], carry);
        let (r3, carry) = mul64_with_carry(r3, k, M.0[3], carry);
        let (r4, carry) = mul64_with_carry(r4, k, M.0[4], carry);
        let (r5, carry) = mul64_with_carry(r5, k, M.0[5], carry);
        let (r6, carry2) = add64_with_carry(r6, 0, carry);

        let k = r1.wrapping_mul(U);
        let (_, carry) = mul64_with_carry(r1, k, M.0[0], 0);
        let (r2, carry) = mul64_with_carry(r2, k, M.0[1], carry);
        let (r3, carry) = mul64_with_carry(r3, k, M.0[2], carry);
        let (r4, carry) = mul64_with_carry(r4, k, M.0[3], carry);
        let (r5, carry) = mul64_with_carry(r5, k, M.0[4], carry);
        let (r6, carry) = mul64_with_carry(r6, k, M.0[5], carry);
        let (r7, carry2) = add64_with_carry(r7, carry2, carry);

        let k = r2.wrapping_mul(U);
        let (_, carry) = mul64_with_carry(r2, k, M.0[0], 0);
        let (r3, carry) = mul64_with_carry(r3, k, M.0[1], carry);
        let (r4, carry) = mul64_with_carry(r4, k, M.0[2], carry);
        let (r5, carry) = mul64_with_carry(r5, k, M.0[3], carry);
        let (r6, carry) = mul64_with_carry(r6, k, M.0[4], carry);
        let (r7, carry) = mul64_with_carry(r7, k, M.0[5], carry);
        let (r8, carry2) = add64_with_carry(r8, carry2, carry);

        let k = r3.wrapping_mul(U);
        let (_, carry) = mul64_with_carry(r3, k, M.0[0], 0);
        let (r4, carry) = mul64_with_carry(r4, k, M.0[1], carry);
        let (r5, carry) = mul64_with_carry(r5, k, M.0[2], carry);
        let (r6, carry) = mul64_with_carry(r6, k, M.0[3], carry);
        let (r7, carry) = mul64_with_carry(r7, k, M.0[4], carry);
        let (r8, carry) = mul64_with_carry(r8, k, M.0[5], carry);
        let (r9, carry2) = add64_with_carry(r9, carry2, carry);

        let k = r4.wrapping_mul(U);
        let (_, carry) = mul64_with_carry(r4, k, M.0[0], 0);
        let (r5, carry) = mul64_with_carry(r5, k, M.0[1], carry);
        let (r6, carry) = mul64_with_carry(r6, k, M.0[2], carry);
        let (r7, carry) = mul64_with_carry(r7, k, M.0[3], carry);
        let (r8, carry) = mul64_with_carry(r8, k, M.0[4], carry);
        let (r9, carry) = mul64_with_carry(r9, k, M.0[5], carry);
        let (r10, carry2) = add64_with_carry(r10, carry2, carry);

        let k = r5.wrapping_mul(U);
        let (_, carry) = mul64_with_carry(r5, k, M.0[0], 0);
        let (r6, carry) = mul64_with_carry(r6, k, M.0[1], carry);
        let (r7, carry) = mul64_with_carry(r7, k, M.0[2], carry);
        let (r8, carry) = mul64_with_carry(r8, k, M.0[3], carry);
        let (r9, carry) = mul64_with_carry(r9, k, M.0[4], carry);
        let (r10, carry) = mul64_with_carry(r10, k, M.0[5], carry);
        let (r11, _) = add64_with_carry(r11, carry2, carry);

        // Result may be within M of the correct value, hence substracting the modulus
        (&Scalar([r6, r7, r8, r9, r10, r11])).sub(&M)
    }

    /// Computes the summation of two scalar elements
    #[inline]
    pub const fn add(&self, rhs: &Self) -> Self {
        let (d0, carry) = add64_with_carry(self.0[0], rhs.0[0], 0);
        let (d1, carry) = add64_with_carry(self.0[1], rhs.0[1], carry);
        let (d2, carry) = add64_with_carry(self.0[2], rhs.0[2], carry);
        let (d3, carry) = add64_with_carry(self.0[3], rhs.0[3], carry);
        let (d4, carry) = add64_with_carry(self.0[4], rhs.0[4], carry);
        let (d5, _) = add64_with_carry(self.0[5], rhs.0[5], carry);

        // Attempt to subtract the modulus, to ensure the value
        // is smaller than the modulus.
        (&Scalar([d0, d1, d2, d3, d4, d5])).sub(&M)
    }

    /// Computes the difference of two scalar elements
    #[inline]
    pub const fn sub(&self, rhs: &Self) -> Self {
        let (d0, borrow) = sub64_with_carry(self.0[0], rhs.0[0], 0);
        let (d1, borrow) = sub64_with_carry(self.0[1], rhs.0[1], borrow);
        let (d2, borrow) = sub64_with_carry(self.0[2], rhs.0[2], borrow);
        let (d3, borrow) = sub64_with_carry(self.0[3], rhs.0[3], borrow);
        let (d4, borrow) = sub64_with_carry(self.0[4], rhs.0[4], borrow);
        let (d5, borrow) = sub64_with_carry(self.0[5], rhs.0[5], borrow);

        // If underflow occurred on the final limb,
        // borrow = 0xfff...fff, otherwise borrow = 0x000...000.
        let (d0, carry) = add64_with_carry(d0, M.0[0] & borrow, 0);
        let (d1, carry) = add64_with_carry(d1, M.0[1] & borrow, carry);
        let (d2, carry) = add64_with_carry(d2, M.0[2] & borrow, carry);
        let (d3, carry) = add64_with_carry(d3, M.0[3] & borrow, carry);
        let (d4, carry) = add64_with_carry(d4, M.0[4] & borrow, carry);
        let (d5, _) = add64_with_carry(d5, M.0[5] & borrow, carry);

        Scalar([d0, d1, d2, d3, d4, d5])
    }

    /// Computes the negation of a scalar element
    #[inline]
    pub const fn neg(&self) -> Self {
        // Subtract `self` from `M` to negate. Ignore the final
        // borrow because it cannot underflow; self is guaranteed to
        // be in the field.
        let (d0, borrow) = sub64_with_carry(M.0[0], self.0[0], 0);
        let (d1, borrow) = sub64_with_carry(M.0[1], self.0[1], borrow);
        let (d2, borrow) = sub64_with_carry(M.0[2], self.0[2], borrow);
        let (d3, borrow) = sub64_with_carry(M.0[3], self.0[3], borrow);
        let (d4, borrow) = sub64_with_carry(M.0[4], self.0[4], borrow);
        let (d5, _) = sub64_with_carry(M.0[5], self.0[5], borrow);

        // `tmp` could be `M` if `self` was zero. Create a mask that is
        // zero if `self` was zero, and `u64::max_value()` if self was nonzero.
        let mask = (((self.0[0] | self.0[1] | self.0[2] | self.0[3] | self.0[4] | self.0[5]) == 0)
            as u64)
            .wrapping_sub(1);

        Scalar([
            d0 & mask,
            d1 & mask,
            d2 & mask,
            d3 & mask,
            d4 & mask,
            d5 & mask,
        ])
    }

    /// Computes the multiplication of two scalar elements
    #[inline]
    pub const fn mul(&self, rhs: &Self) -> Self {
        // Schoolbook multiplication

        let (t0, carry) = mul64_with_carry(0, self.0[0], rhs.0[0], 0);
        let (t1, carry) = mul64_with_carry(0, self.0[0], rhs.0[1], carry);
        let (t2, carry) = mul64_with_carry(0, self.0[0], rhs.0[2], carry);
        let (t3, carry) = mul64_with_carry(0, self.0[0], rhs.0[3], carry);
        let (t4, carry) = mul64_with_carry(0, self.0[0], rhs.0[4], carry);
        let (t5, t6) = mul64_with_carry(0, self.0[0], rhs.0[5], carry);

        let (t1, carry) = mul64_with_carry(t1, self.0[1], rhs.0[0], 0);
        let (t2, carry) = mul64_with_carry(t2, self.0[1], rhs.0[1], carry);
        let (t3, carry) = mul64_with_carry(t3, self.0[1], rhs.0[2], carry);
        let (t4, carry) = mul64_with_carry(t4, self.0[1], rhs.0[3], carry);
        let (t5, carry) = mul64_with_carry(t5, self.0[1], rhs.0[4], carry);
        let (t6, t7) = mul64_with_carry(t6, self.0[1], rhs.0[5], carry);

        let (t2, carry) = mul64_with_carry(t2, self.0[2], rhs.0[0], 0);
        let (t3, carry) = mul64_with_carry(t3, self.0[2], rhs.0[1], carry);
        let (t4, carry) = mul64_with_carry(t4, self.0[2], rhs.0[2], carry);
        let (t5, carry) = mul64_with_carry(t5, self.0[2], rhs.0[3], carry);
        let (t6, carry) = mul64_with_carry(t6, self.0[2], rhs.0[4], carry);
        let (t7, t8) = mul64_with_carry(t7, self.0[2], rhs.0[5], carry);

        let (t3, carry) = mul64_with_carry(t3, self.0[3], rhs.0[0], 0);
        let (t4, carry) = mul64_with_carry(t4, self.0[3], rhs.0[1], carry);
        let (t5, carry) = mul64_with_carry(t5, self.0[3], rhs.0[2], carry);
        let (t6, carry) = mul64_with_carry(t6, self.0[3], rhs.0[3], carry);
        let (t7, carry) = mul64_with_carry(t7, self.0[3], rhs.0[4], carry);
        let (t8, t9) = mul64_with_carry(t8, self.0[3], rhs.0[5], carry);

        let (t4, carry) = mul64_with_carry(t4, self.0[4], rhs.0[0], 0);
        let (t5, carry) = mul64_with_carry(t5, self.0[4], rhs.0[1], carry);
        let (t6, carry) = mul64_with_carry(t6, self.0[4], rhs.0[2], carry);
        let (t7, carry) = mul64_with_carry(t7, self.0[4], rhs.0[3], carry);
        let (t8, carry) = mul64_with_carry(t8, self.0[4], rhs.0[4], carry);
        let (t9, t10) = mul64_with_carry(t9, self.0[4], rhs.0[5], carry);

        let (t5, carry) = mul64_with_carry(t5, self.0[5], rhs.0[0], 0);
        let (t6, carry) = mul64_with_carry(t6, self.0[5], rhs.0[1], carry);
        let (t7, carry) = mul64_with_carry(t7, self.0[5], rhs.0[2], carry);
        let (t8, carry) = mul64_with_carry(t8, self.0[5], rhs.0[3], carry);
        let (t9, carry) = mul64_with_carry(t9, self.0[5], rhs.0[4], carry);
        let (t10, t11) = mul64_with_carry(t10, self.0[5], rhs.0[5], carry);

        Self::montgomery_reduce(t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11)
    }

    /// Computes the square of a scalar element
    #[inline]
    pub const fn square(&self) -> Self {
        let (t1, carry) = mul64_with_carry(0, self.0[0], self.0[1], 0);
        let (t2, carry) = mul64_with_carry(0, self.0[0], self.0[2], carry);
        let (t3, carry) = mul64_with_carry(0, self.0[0], self.0[3], carry);
        let (t4, carry) = mul64_with_carry(0, self.0[0], self.0[4], carry);
        let (t5, t6) = mul64_with_carry(0, self.0[0], self.0[5], carry);

        let (t3, carry) = mul64_with_carry(t3, self.0[1], self.0[2], 0);
        let (t4, carry) = mul64_with_carry(t4, self.0[1], self.0[3], carry);
        let (t5, carry) = mul64_with_carry(t5, self.0[1], self.0[4], carry);
        let (t6, t7) = mul64_with_carry(t6, self.0[1], self.0[5], carry);

        let (t5, carry) = mul64_with_carry(t5, self.0[2], self.0[3], 0);
        let (t6, carry) = mul64_with_carry(t6, self.0[2], self.0[4], carry);
        let (t7, t8) = mul64_with_carry(t7, self.0[2], self.0[5], carry);

        let (t7, carry) = mul64_with_carry(t7, self.0[3], self.0[4], 0);
        let (t8, t9) = mul64_with_carry(t8, self.0[3], self.0[5], carry);

        let (t9, t10) = mul64_with_carry(t9, self.0[4], self.0[5], 0);

        let t11 = t10 >> 63;
        let t10 = (t10 << 1) | (t9 >> 63);
        let t9 = (t9 << 1) | (t8 >> 63);
        let t8 = (t8 << 1) | (t7 >> 63);
        let t7 = (t7 << 1) | (t6 >> 63);
        let t6 = (t6 << 1) | (t5 >> 63);
        let t5 = (t5 << 1) | (t4 >> 63);
        let t4 = (t4 << 1) | (t3 >> 63);
        let t3 = (t3 << 1) | (t2 >> 63);
        let t2 = (t2 << 1) | (t1 >> 63);
        let t1 = t1 << 1;

        let (t0, carry) = mul64_with_carry(0, self.0[0], self.0[0], 0);
        let (t1, carry) = add64_with_carry(t1, 0, carry);
        let (t2, carry) = mul64_with_carry(t2, self.0[1], self.0[1], carry);
        let (t3, carry) = add64_with_carry(t3, 0, carry);
        let (t4, carry) = mul64_with_carry(t4, self.0[2], self.0[2], carry);
        let (t5, carry) = add64_with_carry(t5, 0, carry);
        let (t6, carry) = mul64_with_carry(t6, self.0[3], self.0[3], carry);
        let (t7, carry) = add64_with_carry(t7, 0, carry);
        let (t8, carry) = mul64_with_carry(t8, self.0[4], self.0[4], carry);
        let (t9, carry) = add64_with_carry(t9, 0, carry);
        let (t10, carry) = mul64_with_carry(t10, self.0[5], self.0[5], carry);
        let (t11, _) = add64_with_carry(t11, 0, carry);

        Self::montgomery_reduce(t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11)
    }

    /// Computes the double of a scalar element
    // Can be faster via bitshift
    #[inline]
    pub const fn double(&self) -> Self {
        self.add(self)
    }

    /// Attempts to convert a little-endian byte representation of
    /// a scalar into a `Scalar`, failing if the input is not canonical.
    pub fn from_bytes(bytes: &[u8; 48]) -> CtOption<Scalar> {
        let mut tmp = Scalar([0, 0, 0, 0, 0, 0]);

        tmp.0[0] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[0..8]).unwrap());
        tmp.0[1] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[8..16]).unwrap());
        tmp.0[2] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[16..24]).unwrap());
        tmp.0[3] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[24..32]).unwrap());
        tmp.0[4] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[32..40]).unwrap());
        tmp.0[5] = u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[40..48]).unwrap());

        // Try to subtract the modulus
        let (_, borrow) = sub64_with_carry(tmp.0[0], M.0[0], 0);
        let (_, borrow) = sub64_with_carry(tmp.0[1], M.0[1], borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[2], M.0[2], borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[3], M.0[3], borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[4], M.0[4], borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[5], M.0[5], borrow);

        // If the element is smaller than M then the
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
        assert_eq!(bit_slice.len(), 384);

        let mut result = Scalar::zero();
        for i in 0..384 {
            result = result.double();
            if bit_slice[383 - i] {
                result += Scalar::one();
            }
        }

        result
    }

    /// Outputs the internal representation as 6 64-bit limbs after Montgomery reduction
    pub const fn to_repr(&self) -> [u64; 6] {
        Scalar::montgomery_reduce(
            self.0[0], self.0[1], self.0[2], self.0[3], self.0[4], self.0[5], 0, 0, 0, 0, 0, 0,
        )
        .0
    }

    /// Converts a `Scalar` element into a byte representation in
    /// little-endian byte order.
    pub fn to_bytes(&self) -> [u8; 48] {
        // Turn into canonical form by computing
        // (a.R) / R = a
        let tmp = Scalar::montgomery_reduce(
            self.0[0], self.0[1], self.0[2], self.0[3], self.0[4], self.0[5], 0, 0, 0, 0, 0, 0,
        );

        let mut res = [0; 48];
        res[0..8].copy_from_slice(&tmp.0[0].to_le_bytes());
        res[8..16].copy_from_slice(&tmp.0[1].to_le_bytes());
        res[16..24].copy_from_slice(&tmp.0[2].to_le_bytes());
        res[24..32].copy_from_slice(&tmp.0[3].to_le_bytes());
        res[32..40].copy_from_slice(&tmp.0[4].to_le_bytes());
        res[40..48].copy_from_slice(&tmp.0[5].to_le_bytes());

        res
    }

    /// Converts a 768-bit little endian integer into
    /// a `Scalar` by reducing by the modulus.
    pub fn from_bytes_wide(bytes: &[u8; 96]) -> Self {
        Scalar::from_u768([
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[0..8]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[8..16]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[16..24]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[24..32]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[32..40]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[40..48]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[48..56]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[56..64]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[64..72]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[72..80]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[80..88]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[88..96]).unwrap()),
        ])
    }

    fn from_u768(limbs: [u64; 12]) -> Self {
        // We reduce an arbitrary 768-bit number by decomposing it into two 384-bit digits
        // with the higher bits multiplied by 2^384. Thus, we perform two reductions
        //
        // 1. the lower bits are multiplied by R^2, as normal
        // 2. the upper bits are multiplied by R^2 * 2^384 = R^3
        //
        // and computing their sum in the field. It remains to see that arbitrary 384-bit
        // numbers can be placed into Montgomery form safely using the reduction. The
        // reduction works so long as the product is less than R=2^384 multiplied by
        // the modulus. This holds because for any `c` smaller than the modulus, we have
        // that (2^384 - 1)*c is an acceptable product for the reduction. Therefore, the
        // reduction always works so long as `c` is in the field; in this case it is either the
        // constant `R2` or `R3`.
        let d0 = Scalar([limbs[0], limbs[1], limbs[2], limbs[3], limbs[4], limbs[5]]);
        let d1 = Scalar([limbs[6], limbs[7], limbs[8], limbs[9], limbs[10], limbs[11]]);

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
        let tmp = Scalar::montgomery_reduce(
            self.0[0], self.0[1], self.0[2], self.0[3], self.0[4], self.0[5], 0, 0, 0, 0, 0, 0,
        );

        let (_, borrow) = sub64_with_carry(tmp.0[0], 0x3108_542b_7471_5647, 0);
        let (_, borrow) = sub64_with_carry(tmp.0[1], 0x13ba_7fe8_f765_0104, borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[2], 0xc19e_3e65_2fab_9ae8, borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[3], 0xb953_8f1d_34e3_dc07, borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[4], 0xa978_d997_e2ea_6efc, borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[5], 0x0007_ffd6_605a_3d77, borrow);

        // If the element was smaller, the subtraction will underflow
        // producing a borrow value of 0xffff...ffff, otherwise it will
        // be zero. We create a Choice representing true if there was
        // overflow (and so this element is not lexicographically larger
        // than its negation) and then negate it.

        !Choice::from((borrow as u8) & 1)
    }

    /// Exponentiates `self` by `power`, where `power` is a
    /// little-endian order integer exponent.
    pub fn exp(self, power: &[u64; 6]) -> Self {
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
    pub fn exp_vartime(&self, power: &[u64; 6]) -> Self {
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
        let t1 = self.square(); //          1    2
        let t7 = t1 * self; //              2:   3
        let mut t0 = t1.square(); //        3:   4
        let t4 = t0 * t7; //                4:   7
        let t8 = t4 * t1; //                5:   9
        let t1 = t4 * t0; //                6:   11
        let t13 = t8 * t0; //               7:   13
        let t10 = t1 * t0; //               8:   15
        let t14 = t13 * t0; //              9:   17
        let t12 = t10 * t0; //              10:  19
        let t3 = t14 * t0; //               11:  21
        let t6 = t12 * t0; //               12:  23
        let t2 = t3 * t0; //                13:  25
        let t15 = t6 * t0; //               14:  27
        let t9 = t2 * t0; //                15:  29
        let t11 = t15 * t0; //              16:  31
        t0 = t11.square(); //               17:  62
        square_assign_multi(&mut t0, 4); // 21:  992
        t0 *= &t11; //                      22:  1023
        square_assign_multi(&mut t0, 5); // 27:  32736
        t0 *= &t9; //                       28:  32765
        square_assign_multi(&mut t0, 3); // 31:  262120
        t0 *= &t7; //                       32:  262123
        square_assign_multi(&mut t0, 4); // 36:  4193968
        t0 *= &t7; //                       37:  4193971
        square_assign_multi(&mut t0, 10); //47:  4294626304
        t0 *= &t1; //                       48:  4294626315
        square_assign_multi(&mut t0, 6); // 54:  274856084160
        t0 *= &t14; //                      55:  274856084177
        square_assign_multi(&mut t0, 5); // 60:  8795394693664
        t0 *= &t9; //                       61:  8795394693693
        square_assign_multi(&mut t0, 6); // 67:  562905260396352
        t0 *= &t9; //                       68:  562905260396381
        square_assign_multi(&mut t0, 5); // 73:  18012968332684192
        t0 *= &t9; //                       74:  18012968332684221
        square_assign_multi(&mut t0, 5); // 79:  576414986645895072
        t0 *= &t8; //                       80:  576414986645895081
        square_assign_multi(&mut t0, 5); // 85:  18445279572668642592
        t0 *= &t10; //                      86:  18445279572668642607
        square_assign_multi(&mut t0, 8); // 94:  4721991570603172507392
        t0 *= &t15; //                      95:  4721991570603172507419
        square_assign_multi(&mut t0, 7); // 102: 604414921037206080949632
        t0 *= &t2; //                       103: 604414921037206080949657
        square_assign_multi(&mut t0, 6); // 109: 38682554946381189180778048
        t0 *= &t11; //                      110: 38682554946381189180778079
        square_assign_multi(&mut t0, 5); // 115: 1237841758284198053784898528
        t0 *= &t14; //                      116: 1237841758284198053784898545
        square_assign_multi(&mut t0, 6); // 122: 79221872530188675442233506880
        t0 *= &t9; //                       123: 79221872530188675442233506909
        square_assign_multi(&mut t0, 6); // 129: 5070199841932075228302944442176
        t0 *= &t12; //                      130: 5070199841932075228302944442195
        square_assign_multi(&mut t0, 6); // 136: 324492789883652814611388444300480
        t0 *= &t9; //                       137: 324492789883652814611388444300509
        square_assign_multi(&mut t0, 5); // 142: 10383769276276890067564430217616288
        t0 *= &t11; //                      143: 10383769276276890067564430217616319
        square_assign_multi(&mut t0, 7); // 150: 1329122467363441928648247067854888832
        t0 *= &t6; //                       151: 1329122467363441928648247067854888855
        square_assign_multi(&mut t0, 7); // 158: 170127675822520566866975624685425773440
        t0 *= &t3; //                       159: 170127675822520566866975624685425773461
        square_assign_multi(&mut t0, 5); // 164: 5444085626320658139743219989933624750752
        t0 *= &t4; //                       165: 5444085626320658139743219989933624750759
        square_assign_multi(&mut t0, 7); // 172: 696842960169044241887132158711503968097152
        t0 *= &t10; //                      173: 696842960169044241887132158711503968097167
        square_assign_multi(&mut t0, 8); // 181: 178391797803275325923105832630145015832874752
        t0 *= &t9; //                       182: 178391797803275325923105832630145015832874781
        square_assign_multi(&mut t0, 6); // 188: 11417075059409620859078773288329281013303985984
        t0 *= &t13; //                      189: 11417075059409620859078773288329281013303985997
        square_assign_multi(&mut t0, 5); // 194: 365346401901107867490520745226536992425727551904
        t0 *= &t4; //                       195: 365346401901107867490520745226536992425727551911
        square_assign_multi(&mut t0, 7); // 202: 46764339443341807038786655388996735030493126644608
        t0 *= &t10; //                      203: 46764339443341807038786655388996735030493126644623
        square_assign_multi(&mut t0, 4); // 207: 748229431093468912620586486223947760487890026313968
        t0 *= &t4; //                       208: 748229431093468912620586486223947760487890026313975
        square_assign_multi(&mut t0, 12); //220: 3064747749758848666093922247573290026958397547782041600
        t0 *= &t11; //                      221: 3064747749758848666093922247573290026958397547782041631
        square_assign_multi(&mut t0, 7); // 228: 392287711969132629260022047689381123450674886116101328768
        t0 *= &t7; //                       229: 392287711969132629260022047689381123450674886116101328771
        square_assign_multi(&mut t0, 6); // 235: 25106413566024488272641411052120391900843192711430485041344
        t0 *= &t10; //                      236: 25106413566024488272641411052120391900843192711430485041359
        square_assign_multi(&mut t0, 8); // 244: 6427241872902268997796201229342820326615857334126204170587904
        t0 *= &t11; //                      245: 6427241872902268997796201229342820326615857334126204170587935
        square_assign_multi(&mut t0, 7); // 252: 822686959731490431717913757355881001806829738768154133835255680
        t0 *= &t2; //                       253: 822686959731490431717913757355881001806829738768154133835255705
        square_assign_multi(&mut t0, 5); // 258: 26325982711407693814973240235388192057818551640580932282728182560
        t0 *= &t8; //                       259: 26325982711407693814973240235388192057818551640580932282728182569
        square_assign_multi(&mut t0, 6); // 265: 1684862893530092404158287375064844291700387304997179666094603684416
        t0 *= &t11; //                      266: 1684862893530092404158287375064844291700387304997179666094603684447
        square_assign_multi(&mut t0, 6); // 272: 107831225185925913866130392004150034668824787519819498630054635804608
        t0 *= &t3; //                       273: 107831225185925913866130392004150034668824787519819498630054635804629
        square_assign_multi(&mut t0, 5); // 278: 3450599205949629243716172544132801109402393200634223956161748345748128
        t0 *= &t2; //                       279: 3450599205949629243716172544132801109402393200634223956161748345748153
        square_assign_multi(&mut t0, 5); // 284: 110419174590388135798917521412249635500876582420295166597175947063940896
        t0 *= &t3; //                       285: 110419174590388135798917521412249635500876582420295166597175947063940917
        square_assign_multi(&mut t0, 4); // 289: 1766706793446210172782680342595994168014025318724722665554815153023054672
        t0 *= &t13; //                      290: 1766706793446210172782680342595994168014025318724722665554815153023054685
        square_assign_multi(&mut t0, 11); //301: 3618215512977838433858929341636596056092723852748232019056261433391215994880
        t0 *= &t12; //                      302: 3618215512977838433858929341636596056092723852748232019056261433391215994899
        square_assign_multi(&mut t0, 5); // 307: 115782896415290829883485738932371073794967163287943424609800365868518911836768
        t0 *= &t6; //                       308: 115782896415290829883485738932371073794967163287943424609800365868518911836791
        square_assign_multi(&mut t0, 6); // 314: 7410105370578613112543087291671748722877898450428379175027223415585210357554624
        t0 *= &t12; //                      315: 7410105370578613112543087291671748722877898450428379175027223415585210357554643
        square_assign_multi(&mut t0, 5); // 320: 237123371858515619601378793333495959132092750413708133600871149298726731441748576
        t0 *= &t11; //                      321: 237123371858515619601378793333495959132092750413708133600871149298726731441748607
        square_assign_multi(&mut t0, 5); // 326: 7587947899472499827244121386671870692226968013238660275227876777559255406135955424
        t0 *= &t9; //                       327: 7587947899472499827244121386671870692226968013238660275227876777559255406135955453
        square_assign_multi(&mut t0, 7); // 334: 971257331132479977887247537493999448605051905694548515229168227527584691985402297984
        t0 *= &t10; //                      335: 971257331132479977887247537493999448605051905694548515229168227527584691985402297999
        square_assign_multi(&mut t0, 6); // 341: 62160469192478718584783842399615964710723321964451104974666766561765420287065747071936
        t0 *= &t9; //                       342: 62160469192478718584783842399615964710723321964451104974666766561765420287065747071965
        square_assign_multi(&mut t0, 4); // 346: 994567507079659497356541478393855435371573151431217679594668264988246724593051953151440
        t0 *= &t8; //                       347: 994567507079659497356541478393855435371573151431217679594668264988246724593051953151449
        square_assign_multi(&mut t0, 2); // 349: 3978270028318637989426165913575421741486292605724870718378673059952986898372207812605796
        t0 *= &self; //                     350: 3978270028318637989426165913575421741486292605724870718378673059952986898372207812605797
        square_assign_multi(&mut t0, 8); // 358: 1018437127249571325293098473875307965820490907065566903904940303347964645983285200027084032
        t0 *= &self; //                     359: 1018437127249571325293098473875307965820490907065566903904940303347964645983285200027084033
        square_assign_multi(&mut t0, 6); // 365: 65179976143972564818758302328019709812511418052196281849916179414269737342930252801733378112
        t0 *= &self; //                     366: 65179976143972564818758302328019709812511418052196281849916179414269737342930252801733378113
        square_assign_multi(&mut t0, 6); // 372: 4171518473214244148400531348993261428000730755340562038394635482513263189947536179310936199232
        t0 *= &t7; //                       373: 4171518473214244148400531348993261428000730755340562038394635482513263189947536179310936199235
        square_assign_multi(&mut t0, 4); // 377: 66744295571427906374408501583892182848011692085448992614314167720212211039160578868974979187760
        t0 *= &self; //                     378: 66744295571427906374408501583892182848011692085448992614314167720212211039160578868974979187761
        square_assign_multi(&mut t0, 5); // 383: 2135817458285693003981072050684549851136374146734367763658053367046790753253138523807199334008352
        t0 *= &self; //                     384: 2135817458285693003981072050684549851136374146734367763658053367046790753253138523807199334008353
        square_assign_multi(&mut t0, 9); // 393: 1093538538642274818038308889950489523781823563127996294992923323927956865665606924189286059012276736
        t0 *= &t3; //                       394: 1093538538642274818038308889950489523781823563127996294992923323927956865665606924189286059012276757
        square_assign_multi(&mut t0, 9); // 403: 559891731784844706835614151654650636176293664321534103036376741851113915220790745184914462214285699584
        t0 *= &t3; //                       404: 559891731784844706835614151654650636176293664321534103036376741851113915220790745184914462214285699605
        square_assign_multi(&mut t0, 5); // 409: 17916535417115030618739652852948820357641397258289091297164055739235645287065303845917262790857142387360
        t0 *= &t6; //                       410: 17916535417115030618739652852948820357641397258289091297164055739235645287065303845917262790857142387383
        square_assign_multi(&mut t0, 2); // 412: 71666141668460122474958611411795281430565589033156365188656222956942581148261215383669051163428569549532
        t0 *= &self; //                     413: 71666141668460122474958611411795281430565589033156365188656222956942581148261215383669051163428569549533
        square_assign_multi(&mut t0, 6); // 419: 4586633066781447838397351130354898011556197698122007372073998269244325193488717784554819274459428451170112
        t0 *= &t4; //                       420: 4586633066781447838397351130354898011556197698122007372073998269244325193488717784554819274459428451170119
        square_assign_multi(&mut t0, 8); // 428: 1174178065096050646629721889370853890958386610719233887250943556926547249533111752846033734261613683499550464
        t0 *= &t3; //                       429: 1174178065096050646629721889370853890958386610719233887250943556926547249533111752846033734261613683499550485
        square_assign_multi(&mut t0, 6); // 435: 75147396166147241384302200919734649021336743086030968784060387643299023970119152182146158992743275743971231040
        t0 *= &t2; //                       436: 75147396166147241384302200919734649021336743086030968784060387643299023970119152182146158992743275743971231065
        square_assign_multi(&mut t0, 7); // 443: 9618866709266846897190681717726035074731103115011964004359729618342275068175251479314708351071139295228317576320
        t0 *= &t1; //                       444: 9618866709266846897190681717726035074731103115011964004359729618342275068175251479314708351071139295228317576331 = M - 2

        CtOption::new(t0, !self.ct_eq(&Self::zero()))
    }

    /// Computes a random `Scalar` element
    pub fn random(mut rng: impl RngCore + CryptoRng) -> Self {
        let mut buf = [0; 96];
        rng.fill_bytes(&mut buf);
        Self::from_bytes_wide(&buf)
    }

    #[allow(unused)]
    /// Constructs a `Scalar` element without checking that it is
    /// canonical.
    pub const fn from_raw_unchecked(v: [u64; 6]) -> Self {
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
        Scalar::new([value_low, value_high, 0, 0, 0, 0])
    }
}

impl From<u64> for Scalar {
    /// Converts a 64-bit value into a field element.
    fn from(value: u64) -> Self {
        Scalar([value, 0, 0, 0, 0, 0]) * R2
    }
}

impl From<u32> for Scalar {
    /// Converts a 32-bit value into a field element.
    fn from(value: u32) -> Self {
        Scalar([value as u64, 0, 0, 0, 0, 0]) * R2
    }
}

impl From<u16> for Scalar {
    /// Converts a 16-bit value into a field element.
    fn from(value: u16) -> Self {
        Scalar([value as u64, 0, 0, 0, 0, 0]) * R2
    }
}

impl From<u8> for Scalar {
    /// Converts an 8-bit value into a field element.
    fn from(value: u8) -> Self {
        Scalar([value as u64, 0, 0, 0, 0, 0]) * R2
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bitvec::view::AsBits;
    use rand::thread_rng;

    const LARGEST: Scalar = Scalar([
        0x6210a856e8e2ac8c,
        0x2774ffd1eeca0208,
        0x833c7cca5f5735d0,
        0x72a71e3a69c7b80f,
        0x52f1b32fc5d4ddf9,
        0x000fffacc0b47aef,
    ]);

    // DISPLAY
    // ================================================================================================

    #[test]
    fn test_debug() {
        assert_eq!(
            format!("{:?}", Scalar::zero()),
            "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
        );
        assert_eq!(
            format!("{:?}", Scalar::one()),
            "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001"
        );
        assert_eq!(
            format!("{:?}", R2),
            "0x000533f4b8510ad0e4cd03a2b22068d58e1c5963847f07cc38335a0a8ca2fd88b002e1135fdf79def57a9171d5373000"
        );
    }

    #[test]
    fn test_to_repr() {
        assert_eq!(
            format!("{:?}", Scalar::zero().to_repr()),
            "[0, 0, 0, 0, 0, 0]"
        );
        assert_eq!(
            format!("{:?}", Scalar::one().to_repr()),
            "[1, 0, 0, 0, 0, 0]"
        );
        assert_eq!(
            format!("{:?}", R2.to_repr()),
            "[17688610404545540096, 12682946973957847518, 4049679491291872648, 10240157936693217228, 16486837808181307605, 1464501040909008]"
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
                0x6210a856e8e2ac8b,
                0x2774ffd1eeca0208,
                0x833c7cca5f5735d0,
                0x72a71e3a69c7b80f,
                0x52f1b32fc5d4ddf9,
                0x000fffacc0b47aef,
            ])
        );

        assert_eq!(tmp, LARGEST.double());

        let mut tmp = LARGEST;
        tmp += &Scalar([1, 0, 0, 0, 0, 0]);

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

        assert_eq!(tmp, Scalar([1, 0, 0, 0, 0, 0]));

        let tmp = -&Scalar::zero();
        assert_eq!(tmp, Scalar::zero());
        let tmp = -&Scalar([1, 0, 0, 0, 0, 0]);
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
            let pow2 = tmp.exp(&[2, 0, 0, 0, 0, 0]);
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
            0x6210a856e8e2ac8b,
            0x2774ffd1eeca0208,
            0x833c7cca5f5735d0,
            0x72a71e3a69c7b80f,
            0x52f1b32fc5d4ddf9,
            0x000fffacc0b47aef,
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
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ]
        );

        assert_eq!(
            Scalar::one().to_bytes(),
            [
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ]
        );

        assert_eq!(
            R2.to_bytes(),
            [
                0, 48, 55, 213, 113, 145, 122, 245, 222, 121, 223, 95, 19, 225, 2, 176, 136, 253,
                162, 140, 10, 90, 51, 56, 204, 7, 127, 132, 99, 89, 28, 142, 213, 104, 32, 178,
                162, 3, 205, 228, 208, 10, 81, 184, 244, 51, 5, 0
            ]
        );

        assert_eq!(
            (-&Scalar::one()).to_bytes(),
            [
                140, 172, 226, 232, 86, 168, 16, 98, 8, 2, 202, 238, 209, 255, 116, 39, 208, 53,
                87, 95, 202, 124, 60, 131, 15, 184, 199, 105, 58, 30, 167, 114, 249, 221, 212, 197,
                47, 179, 241, 82, 239, 122, 180, 192, 172, 255, 15, 0
            ]
        );
    }

    #[test]
    fn test_from_bytes() {
        assert_eq!(
            Scalar::from_bytes(&[
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
            .unwrap(),
            Scalar::zero()
        );

        assert_eq!(
            Scalar::from_bytes(&[
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
            .unwrap(),
            Scalar::one()
        );

        assert_eq!(
            Scalar::from_bytes(&[
                0, 48, 55, 213, 113, 145, 122, 245, 222, 121, 223, 95, 19, 225, 2, 176, 136, 253,
                162, 140, 10, 90, 51, 56, 204, 7, 127, 132, 99, 89, 28, 142, 213, 104, 32, 178,
                162, 3, 205, 228, 208, 10, 81, 184, 244, 51, 5, 0
            ])
            .unwrap(),
            R2
        );

        // -1 should work
        assert_eq!(
            Scalar::from_bytes(&[
                140, 172, 226, 232, 86, 168, 16, 98, 8, 2, 202, 238, 209, 255, 116, 39, 208, 53,
                87, 95, 202, 124, 60, 131, 15, 184, 199, 105, 58, 30, 167, 114, 249, 221, 212, 197,
                47, 179, 241, 82, 239, 122, 180, 192, 172, 255, 15, 0
            ])
            .unwrap(),
            -Scalar::one(),
        );

        // M is invalid
        assert!(bool::from(
            Scalar::from_bytes(&[
                141, 172, 226, 232, 86, 168, 16, 98, 8, 2, 202, 238, 209, 255, 116, 39, 208, 53,
                87, 95, 202, 124, 60, 131, 15, 184, 199, 105, 58, 30, 167, 114, 249, 221, 212, 197,
                47, 179, 241, 82, 239, 122, 180, 192, 172, 255, 15, 0
            ])
            .is_none()
        ));

        // Anything larger than the M is invalid
        assert!(bool::from(
            Scalar::from_bytes(&[
                142, 172, 226, 232, 86, 168, 16, 98, 8, 2, 202, 238, 209, 255, 116, 39, 208, 53,
                87, 95, 202, 124, 60, 131, 15, 184, 199, 105, 58, 30, 167, 114, 249, 221, 212, 197,
                47, 179, 241, 82, 239, 122, 180, 192, 172, 255, 15, 0
            ])
            .is_none()
        ));
        assert!(bool::from(
            Scalar::from_bytes(&[
                141, 172, 226, 232, 86, 168, 16, 98, 8, 2, 202, 238, 209, 255, 116, 39, 208, 53,
                87, 95, 202, 124, 60, 131, 15, 184, 199, 105, 58, 30, 167, 114, 249, 221, 212, 197,
                47, 179, 241, 82, 239, 122, 180, 192, 172, 255, 255, 255
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
            141, 172, 226, 232, 86, 168, 16, 98, 8, 2, 202, 238, 209, 255, 116, 39, 208, 53, 87,
            95, 202, 124, 60, 131, 15, 184, 199, 105, 58, 30, 167, 114, 249, 221, 212, 197, 47,
            179, 241, 82, 239, 122, 180, 192, 172, 255, 15, 0,
        ];
        assert_eq!(
            Scalar::from_bits_vartime(bytes.as_bits::<Lsb0>()),
            Scalar::zero()
        );
    }

    #[test]
    fn test_from_u768_zero() {
        assert_eq!(
            Scalar::zero(),
            Scalar::from_u768([M.0[0], M.0[1], M.0[2], M.0[3], M.0[4], M.0[5], 0, 0, 0, 0, 0, 0])
        );
    }

    #[test]
    fn test_from_u768_r() {
        assert_eq!(R, Scalar::from_u768([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]));
    }

    #[test]
    fn test_from_u768_r2() {
        assert_eq!(R2, Scalar::from_u768([0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]));
    }

    #[test]
    fn test_from_u768_max() {
        let max_u64 = 0xffff_ffff_ffff_ffff;
        assert_eq!(
            R3 - R,
            Scalar::from_u768([
                max_u64, max_u64, max_u64, max_u64, max_u64, max_u64, max_u64, max_u64, max_u64,
                max_u64, max_u64, max_u64
            ])
        );
    }

    #[test]
    fn test_from_bytes_wide_r2() {
        assert_eq!(
            R2,
            Scalar::from_bytes_wide(&[
                0, 48, 55, 213, 113, 145, 122, 245, 222, 121, 223, 95, 19, 225, 2, 176, 136, 253,
                162, 140, 10, 90, 51, 56, 204, 7, 127, 132, 99, 89, 28, 142, 213, 104, 32, 178,
                162, 3, 205, 228, 208, 10, 81, 184, 244, 51, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0,
            ])
        );
    }

    #[test]
    fn test_from_bytes_wide_negative_one() {
        assert_eq!(
            -&Scalar::one(),
            Scalar::from_bytes_wide(&[
                140, 172, 226, 232, 86, 168, 16, 98, 8, 2, 202, 238, 209, 255, 116, 39, 208, 53,
                87, 95, 202, 124, 60, 131, 15, 184, 199, 105, 58, 30, 167, 114, 249, 221, 212, 197,
                47, 179, 241, 82, 239, 122, 180, 192, 172, 255, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
        );
    }

    #[test]
    fn test_from_bytes_wide_maximum() {
        assert_eq!(
            Scalar([
                0xca85652242a66d74,
                0x900add215229df2b,
                0x2b063c4f9b979d54,
                0x72669a8ab069ea77,
                0x56f3c296b743abd7,
                0x000562c9b45ccf6f,
            ]),
            Scalar::from_bytes_wide(&[0xff; 96])
        );
    }

    #[test]
    fn test_lexicographically_largest() {
        // a = 3270721541412937332678431542494659405168602257955891291263700876781669849818040027313410211252413088299582187994
        let a = Scalar::new([
            0x96230c7e33c769da,
            0x9d952e94df9317e9,
            0x4f918f79f9ae9b6a,
            0x17613e631e56f5b1,
            0xd80e0b7e3e1ebc99,
            0x000570a8fc1ff080,
        ]);

        // b = 6348145167853909564512250175231375669562500857056072713096028741560605218357211452001298139818726206928735388339
        let b = Scalar::new([
            0xcbed9bd8b51b42b3,
            0x89dfd13d0f36ea1e,
            0x33aaed5065a89a65,
            0x5b45dfd74b70c25e,
            0x7ae3a7b187b62160,
            0x000a8f03c4948a6e,
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
}
