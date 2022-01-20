// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This module provides an implementation of the scalar field Fq
//! of the STARK-friendly cheetah curve, with characteristic
//! q = 0x26337f752795f77cb6b6ebb9a18fecc9f2f264f035242b271e13aee130956aa5.

use core::{
    convert::{TryFrom, TryInto},
    fmt::{self, Debug, Display, Formatter},
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::utils::{
    add64_with_carry, mul64_with_carry, shl64_by_u32_with_carry, square_assign_multi,
    sub64_with_carry,
};

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

// Field modulus = 10230194559610033405039867617070259612247645045591847851798073552054039295467
const M: Scalar = Scalar([
    0x91776f9d741b35eb,
    0xa0a87b464e777106,
    0x9677863fe04d28d8,
    0x169e15bdd5149539,
]);

// 2^256 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R: Scalar = Scalar([
    0xbfde343c02d4aee7,
    0x18c2b3faa0de24b7,
    0x88dd3b415caf3eb1,
    0x073510d7d81d9686,
]);

// 2^512 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R2: Scalar = Scalar([
    0xeaa8a49e7a002ce6,
    0x7876e0d513a7b125,
    0xb0d96d540c4c36bb,
    0x05ee76ff16800566,
]);

// 2^768 mod M; this is used for conversion of elements into Montgomery representation.
pub(crate) const R3: Scalar = Scalar([
    0xbe225474b993a745,
    0x8ca6ee9b17fe89b7,
    0x0a4fd990d07c8929,
    0x0e9fba567142babb,
]);

// Multiplicative generator g of order q-1
// g = 3
//   = 0x159f32878858c3939a97b1c4160dbc134a481befe29a6e273f9a9cb4087e0cb5 in Montgomery form
const GENERATOR: Scalar = Scalar([
    0x3f9a9cb4087e0cb5,
    0x4a481befe29a6e27,
    0x9a97b1c4160dbc13,
    0x159f32878858c393,
]);

// Two-adicity of the field: (q-1) % 2^1 = 0
const TWO_ADICITY: u32 = 1;

// 2^1 root of unity = 0x169e15bdd51495399677863fe04d28d8a0a87b464e77710691776f9d741b35ea
//                   = 0xf6904e5fcf6feb30d9a4afe839dea2787e5c74bad994c4ed1993b6171468704 in Montgomery form
const TWO_ADIC_ROOT_OF_UNITY: Scalar = Scalar([
    0xd1993b6171468704,
    0x87e5c74bad994c4e,
    0x0d9a4afe839dea27,
    0x0f6904e5fcf6feb3,
]);

// -M^{-1} mod 2^64; this is used during element multiplication.
const U: u64 = 18156396707220272445;

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

        // The result may be within M of the correct value,
        // hence substracting the modulus
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
        // Subtract `self` from `M` to negate. Ignore the final
        // borrow because it cannot underflow; self is guaranteed to
        // be in the field.
        let (d0, borrow) = sub64_with_carry(M.0[0], self.0[0], 0);
        let (d1, borrow) = sub64_with_carry(M.0[1], self.0[1], borrow);
        let (d2, borrow) = sub64_with_carry(M.0[2], self.0[2], borrow);
        let (d3, _) = sub64_with_carry(M.0[3], self.0[3], borrow);

        // `tmp` could be `M` if `self` was zero. Create a mask that is
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
        // Tonelli-Shanks algorithm for M mod 4 = 3
        // Using a custom chain for faster exponentiation
        // found using https://github.com/kwantam/addchain for (M + 1) / 4

        let t4 = self.square(); //              1: 2
        let t1 = t4 * self; //                  2: 3
        let t11 = t1 * t4; //                   3: 5
        let t8 = t11 * t4; //                   4: 7
        let t13 = t8 * t4; //                   5: 9
        let t12 = t13 * t4; //                  6: 11
        let t3 = t12 * t4; //                   7: 13
        let t2 = t3 * t4; //                    8: 15
        let t7 = t2 * t4; //                    9: 17
        let t14 = t7 * t4; //                  10: 19
        let t10 = t14 * t4; //                 11: 21
        let mut t0 = t12.square(); //          12: 22
        let t15 = t10 * t4; //                 13: 23
        let t9 = t0 * t1; //                   14: 25
        let t4 = t0 * t11; //                  15: 27
        let t5 = t0 * t8; //                   16: 29
        let t6 = t0 * t13; //                  17: 31
        t0 *= t15; //                          18: 45
        square_assign_multi(&mut t0, 6); //    24: 2880
        t0 *= t2; //                           25: 2895
        square_assign_multi(&mut t0, 9); //    34: 1482240
        t0 *= t10; //                          35: 1482261
        square_assign_multi(&mut t0, 5); //    40: 47432352
        t0 *= t15; //                          41: 47432375
        square_assign_multi(&mut t0, 5); //    46: 1517836000
        t0 *= t15; //                          47: 1517836023
        square_assign_multi(&mut t0, 6); //    53: 97141505472
        t0 *= t10; //                          54: 97141505493
        square_assign_multi(&mut t0, 6); //    60: 6217056351552
        t0 *= t11; //                          61: 6217056351557
        square_assign_multi(&mut t0, 6); //    67: 397891606499648
        t0 *= t13; //                          68: 397891606499657
        square_assign_multi(&mut t0, 4); //    72: 6366265703994512
        t0 *= t11; //                          73: 6366265703994517
        square_assign_multi(&mut t0, 5); //    78: 203720502527824544
        t0 *= t8; //                           79: 203720502527824551
        square_assign_multi(&mut t0, 7); //    86: 26076224323561542528
        t0 *= t9; //                           87: 26076224323561542553
        square_assign_multi(&mut t0, 6); //    93: 1668878356707938723392
        t0 *= t9; //                           94: 1668878356707938723417
        square_assign_multi(&mut t0, 2); //    96: 6675513426831754893668
        t0 *= t1; //                           97: 6675513426831754893671
        square_assign_multi(&mut t0, 5); //   102: 213616429658616156597472
        t0 *= t2; //                          103: 213616429658616156597487
        square_assign_multi(&mut t0, 6); //   109: 13671451498151434022239168
        t0 *= t1; //                          110: 13671451498151434022239171
        square_assign_multi(&mut t0, 8); //   118: 3499891583526767109693227776
        t0 *= t6; //                          119: 3499891583526767109693227807
        square_assign_multi(&mut t0, 4); //   123: 55998265336428273755091644912
        t0 *= t2; //                          124: 55998265336428273755091644927
        square_assign_multi(&mut t0, 11); //  135: 114684447409005104650427688810496
        t0 *= t14; //                         136: 114684447409005104650427688810515
        square_assign_multi(&mut t0, 5); //   141: 3669902317088163348813686041936480
        t0 *= t13; //                         142: 3669902317088163348813686041936489
        square_assign_multi(&mut t0, 6); //   148: 234873748293642454324075906683935296
        t0 *= t7; //                          149: 234873748293642454324075906683935313
        square_assign_multi(&mut t0, 4); //   153: 3757979972698279269185214506942965008
        t0 *= t12; //                         154: 3757979972698279269185214506942965019
        square_assign_multi(&mut t0, 6); //   160: 240510718252689873227853728444349761216
        t0 *= t11; //                         161: 240510718252689873227853728444349761221
        square_assign_multi(&mut t0, 10); //  171: 246282975490754430185322217927014155490304
        t0 *= t10; //                         172: 246282975490754430185322217927014155490325
        square_assign_multi(&mut t0, 8); //   180: 63048441725633134127442487789315623805523200
        t0 *= t2; //                          181: 63048441725633134127442487789315623805523215
        square_assign_multi(&mut t0, 5); //   186: 2017550135220260292078159609258099961776742880
        t0 *= t3; //                          187: 2017550135220260292078159609258099961776742893
        square_assign_multi(&mut t0, 8); //   195: 516492834616386634772008859970073590214846180608
        t0 *= t9; //                          196: 516492834616386634772008859970073590214846180633
        square_assign_multi(&mut t0, 5); //   201: 16527770707724372312704283519042354886875077780256
        t0 *= t8; //                          202: 16527770707724372312704283519042354886875077780263
        square_assign_multi(&mut t0, 7); //   209: 2115554650588719656026148290437421425520009955873664
        t0 *= t5; //                          210: 2115554650588719656026148290437421425520009955873693
        square_assign_multi(&mut t0, 5); //   215: 67697748818839028992836745293997485616640318587958176
        t0 *= t4; //                          216: 67697748818839028992836745293997485616640318587958203
        square_assign_multi(&mut t0, 5); //   221: 2166327962202848927770775849407919539732490194814662496
        t0 *= t7; //                          222: 2166327962202848927770775849407919539732490194814662513
        square_assign_multi(&mut t0, 9); //   231: 1109159916647858651018637234896854804343034979745107206656
        t0 *= t3; //                          232: 1109159916647858651018637234896854804343034979745107206669
        square_assign_multi(&mut t0, 7); //   239: 141972469330925907330385566066797414955908477407373722453632
        t0 *= t7; //                          240: 141972469330925907330385566066797414955908477407373722453649
        square_assign_multi(&mut t0, 6); //   246: 9086238037179258069144676228275034557178142554071918237033536
        t0 *= t5; //                          247: 9086238037179258069144676228275034557178142554071918237033565
        square_assign_multi(&mut t0, 5); //   252: 290759617189736258212629639304801105829700561730301383585074080
        t0 *= t4; //                          253: 290759617189736258212629639304801105829700561730301383585074107
        square_assign_multi(&mut t0, 6); //   259: 18608615500143120525608296915507270773100835950739288549444742848
        t0 *= t6; //                          260: 18608615500143120525608296915507270773100835950739288549444742879
        square_assign_multi(&mut t0, 7); //   267: 2381902784018319427277862005184930658956907001694628934328927088512
        t0 *= t5; //                          268: 2381902784018319427277862005184930658956907001694628934328927088541
        square_assign_multi(&mut t0, 6); //   274: 152441778177172443345783168331835562173242048108456251797051333666624
        t0 *= t5; //                          275: 152441778177172443345783168331835562173242048108456251797051333666653
        square_assign_multi(&mut t0, 10); //  285: 156100380853424581986081964371799615665399857263059201840180565674652672
        t0 *= t4; //                          286: 156100380853424581986081964371799615665399857263059201840180565674652699
        square_assign_multi(&mut t0, 6); //   292: 9990424374619173247109245719795175402585590864835788917771556203177772736
        t0 *= t3; //                          293: 9990424374619173247109245719795175402585590864835788917771556203177772749
        square_assign_multi(&mut t0, 5); //   298: 319693579987813543907495863033445612882738907674745245368689798501688727968
        t0 *= t2; //                          299: 319693579987813543907495863033445612882738907674745245368689798501688727983
        square_assign_multi(&mut t0, 3); //   302: 2557548639902508351259966904267564903061911261397961962949518388013509823864
        t0 *= t1; //                          303: 2557548639902508351259966904267564903061911261397961962949518388013509823867 = (M + 1) / 4

        CtOption::new(t0, (t0 * t0).ct_eq(self))
    }

    /// Computes the double of a scalar element
    #[inline]
    pub const fn double(&self) -> Self {
        let (d0, carry) = shl64_by_u32_with_carry(self.0[0], 1, 0);
        let (d1, carry) = shl64_by_u32_with_carry(self.0[1], 1, carry);
        let (d2, carry) = shl64_by_u32_with_carry(self.0[2], 1, carry);
        let (d3, _carry) = shl64_by_u32_with_carry(self.0[3], 1, carry);

        // Attempt to subtract the modulus, to ensure the value
        // is smaller than the modulus.
        (&Scalar([d0, d1, d2, d3])).sub(&M)
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

        // If the element is smaller than M then the
        // subtraction will underflow, producing a borrow value
        // of 0xffff...ffff. Otherwise, it'll be zero.
        let is_some = (borrow as u8) & 1;

        // Convert to Montgomery form by computing
        // (a.R^0 * R^2) / R = a.R
        tmp *= &R2;

        CtOption::new(tmp, Choice::from(is_some))
    }

    /// Converts a little-endian bit sequence into a Scalar element
    pub fn from_bits(bit_slice: &BitSlice<Lsb0, u8>) -> Scalar {
        assert_eq!(bit_slice.len(), 256);

        let mut result = Scalar::zero();
        for i in (0..256).rev() {
            result = result.double();
            let tmp = Scalar::conditional_select(
                &Scalar::zero(),
                &Scalar::one(),
                Choice::from(bit_slice[i] as u8),
            );
            result += tmp;
        }

        result
    }

    /// Outputs the internal representation as 4 64-bit limbs after Montgomery reduction
    pub const fn output_reduced_limbs(&self) -> [u64; 4] {
        Scalar::montgomery_reduce(self.0[0], self.0[1], self.0[2], self.0[3], 0, 0, 0, 0).0
    }

    /// Outputs the internal representation as 4 64-bit limbs without Montgomery reduction
    /// This is intended for uses like re-interpreting the type containing the internal value.
    pub const fn output_unreduced_limbs(&self) -> [u64; 4] {
        self.0
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
        // and computing their sum in the field.
        // The reduction works so long as the product is less than R=2^256 multiplied by
        // the modulus. This holds because for any `c` smaller than the modulus, we have
        // that (2^256 - 1)*c is an acceptable product for the reduction. Therefore, the
        // reduction always works so long as `c` is in the field; in this case it is either the
        // constant `R2` or `R3`.
        let d0 = Scalar([limbs[0], limbs[1], limbs[2], limbs[3]]);
        let d1 = Scalar([limbs[4], limbs[5], limbs[6], limbs[7]]);

        // Convert to Montgomery form
        d0 * R2 + d1 * R3
    }

    /// Converts a `Scalar` element given as byte representation into a radix-16
    /// representation, where each resulting coefficient is in [-8;8).
    ///
    /// The resulting decomposition `[a_0, ..., a_63]` is such that
    /// `sum(a_j * 2^(j * 4)) == a`.
    pub(crate) fn bytes_to_radix_16(bytes: &[u8; 32]) -> [i8; 64] {
        let mut result = [0i8; 64];

        // Convert from bytes to radix-16
        for i in 0..32 {
            result[2 * i] = (bytes[i] & 0xf) as i8;
            result[2 * i + 1] = ((bytes[i] >> 4) & 0xf) as i8;
        }

        // Shift every coefficients from [0; 16) to [-8; 8)
        for i in 0..63 {
            let carry = (result[i] + 8) >> 4;
            result[i] -= carry << 4;
            result[i + 1] += carry;
        }

        result
    }

    /// Returns whether or not this element is strictly lexicographically
    /// larger than its negation.
    pub fn lexicographically_largest(&self) -> Choice {
        // This can be determined by checking to see if the element is
        // larger than (M - 1) // 2. If we subtract by ((M - 1) // 2) + 1
        // and there is no underflow, then the element must be larger than
        // (M - 1) // 2.

        // First, because self is in Montgomery form we need to reduce it
        let tmp = Scalar::montgomery_reduce(self.0[0], self.0[1], self.0[2], self.0[3], 0, 0, 0, 0);
        let (_, borrow) = sub64_with_carry(tmp.0[0], 0x48bbb7ceba0d9af6, 0);
        let (_, borrow) = sub64_with_carry(tmp.0[1], 0x50543da3273bb883, borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[2], 0xcb3bc31ff026946c, borrow);
        let (_, borrow) = sub64_with_carry(tmp.0[3], 0x0b4f0adeea8a4a9c, borrow);

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
        let t4 = self.square(); //    1: 2
        let t14 = t4 * self; //    2: 3
        let t11 = t14 * t4; //    3: 5
        let t8 = t11 * t4; //    4: 7
        let t1 = t8 * t4; //    5: 9
        let t12 = t1 * t4; //    6: 11
        let t3 = t12 * t4; //    7: 13
        let t2 = t3 * t4; //    8: 15
        let t7 = t2 * t4; //    9: 17
        let t13 = t7 * t4; //   10: 19
        let t10 = t13 * t4; //   11: 21
        let mut t0 = t12.square(); //   12: 22
        let t15 = t10 * t4; //   13: 23
        let t9 = t0 * t14; //   14: 25
        let t4 = t0 * t11; //   15: 27
        let t5 = t0 * t8; //   16: 29
        let t6 = t0 * t1; //   17: 31
        t0 *= t15; //   18: 45
        square_assign_multi(&mut t0, 6); //   24: 2880
        t0 *= t2; //   25: 2895
        square_assign_multi(&mut t0, 9); //   34: 1482240
        t0 *= t10; //   35: 1482261
        square_assign_multi(&mut t0, 5); //   40: 47432352
        t0 *= t15; //   41: 47432375
        square_assign_multi(&mut t0, 5); //   46: 1517836000
        t0 *= t15; //   47: 1517836023
        square_assign_multi(&mut t0, 6); //   53: 97141505472
        t0 *= t10; //   54: 97141505493
        square_assign_multi(&mut t0, 6); //   60: 6217056351552
        t0 *= t11; //   61: 6217056351557
        square_assign_multi(&mut t0, 6); //   67: 397891606499648
        t0 *= t1; //   68: 397891606499657
        square_assign_multi(&mut t0, 4); //   72: 6366265703994512
        t0 *= t11; //   73: 6366265703994517
        square_assign_multi(&mut t0, 5); //   78: 203720502527824544
        t0 *= t8; //   79: 203720502527824551
        square_assign_multi(&mut t0, 7); //   86: 26076224323561542528
        t0 *= t9; //   87: 26076224323561542553
        square_assign_multi(&mut t0, 6); //   93: 1668878356707938723392
        t0 *= t9; //   94: 1668878356707938723417
        square_assign_multi(&mut t0, 2); //   96: 6675513426831754893668
        t0 *= t14; //   97: 6675513426831754893671
        square_assign_multi(&mut t0, 5); //   102: 213616429658616156597472
        t0 *= t2; //  103: 213616429658616156597487
        square_assign_multi(&mut t0, 6); //   109: 13671451498151434022239168
        t0 *= t14; //  110: 13671451498151434022239171
        square_assign_multi(&mut t0, 8); //   118: 3499891583526767109693227776
        t0 *= t6; //  119: 3499891583526767109693227807
        square_assign_multi(&mut t0, 4); //   123: 55998265336428273755091644912
        t0 *= t2; //  124: 55998265336428273755091644927
        square_assign_multi(&mut t0, 11); //  135: 114684447409005104650427688810496
        t0 *= t13; //                         136: 114684447409005104650427688810515
        square_assign_multi(&mut t0, 5); //   141: 3669902317088163348813686041936480
        t0 *= t1; //                          142: 3669902317088163348813686041936489
        square_assign_multi(&mut t0, 6); //   148: 234873748293642454324075906683935296
        t0 *= t7; //                          149: 234873748293642454324075906683935313
        square_assign_multi(&mut t0, 4); //   153: 3757979972698279269185214506942965008
        t0 *= t12; //                         154: 3757979972698279269185214506942965019
        square_assign_multi(&mut t0, 6); //   160: 240510718252689873227853728444349761216
        t0 *= t11; //                         161: 240510718252689873227853728444349761221
        square_assign_multi(&mut t0, 10); //  171: 246282975490754430185322217927014155490304
        t0 *= t10; //                         172: 246282975490754430185322217927014155490325
        square_assign_multi(&mut t0, 8); //   180: 63048441725633134127442487789315623805523200
        t0 *= t2; //                          181: 63048441725633134127442487789315623805523215
        square_assign_multi(&mut t0, 5); //   186: 2017550135220260292078159609258099961776742880
        t0 *= t3; //                          187: 2017550135220260292078159609258099961776742893
        square_assign_multi(&mut t0, 8); //   195: 516492834616386634772008859970073590214846180608
        t0 *= t9; //                          196: 516492834616386634772008859970073590214846180633
        square_assign_multi(&mut t0, 5); //   201: 16527770707724372312704283519042354886875077780256
        t0 *= t8; //                          202: 16527770707724372312704283519042354886875077780263
        square_assign_multi(&mut t0, 7); //   209: 2115554650588719656026148290437421425520009955873664
        t0 *= t5; //                          210: 2115554650588719656026148290437421425520009955873693
        square_assign_multi(&mut t0, 5); //   215: 67697748818839028992836745293997485616640318587958176
        t0 *= t4; //                          216: 67697748818839028992836745293997485616640318587958203
        square_assign_multi(&mut t0, 5); //   221: 2166327962202848927770775849407919539732490194814662496
        t0 *= t7; //                          222: 2166327962202848927770775849407919539732490194814662513
        square_assign_multi(&mut t0, 9); //   231: 1109159916647858651018637234896854804343034979745107206656
        t0 *= t3; //                          232: 1109159916647858651018637234896854804343034979745107206669
        square_assign_multi(&mut t0, 7); //   239: 141972469330925907330385566066797414955908477407373722453632
        t0 *= t7; //                          240: 141972469330925907330385566066797414955908477407373722453649
        square_assign_multi(&mut t0, 6); //   246: 9086238037179258069144676228275034557178142554071918237033536
        t0 *= t5; //                          247: 9086238037179258069144676228275034557178142554071918237033565
        square_assign_multi(&mut t0, 5); //   252: 290759617189736258212629639304801105829700561730301383585074080
        t0 *= t4; //                          253: 290759617189736258212629639304801105829700561730301383585074107
        square_assign_multi(&mut t0, 6); //   259: 18608615500143120525608296915507270773100835950739288549444742848
        t0 *= t6; //                          260: 18608615500143120525608296915507270773100835950739288549444742879
        square_assign_multi(&mut t0, 7); //   267: 2381902784018319427277862005184930658956907001694628934328927088512
        t0 *= t5; //                          268: 2381902784018319427277862005184930658956907001694628934328927088541
        square_assign_multi(&mut t0, 6); //   274: 152441778177172443345783168331835562173242048108456251797051333666624
        t0 *= t5; //                          275: 152441778177172443345783168331835562173242048108456251797051333666653
        square_assign_multi(&mut t0, 10); //  285: 156100380853424581986081964371799615665399857263059201840180565674652672
        t0 *= t4; //                          286: 156100380853424581986081964371799615665399857263059201840180565674652699
        square_assign_multi(&mut t0, 6); //   292: 9990424374619173247109245719795175402585590864835788917771556203177772736
        t0 *= t3; //                          293: 9990424374619173247109245719795175402585590864835788917771556203177772749
        square_assign_multi(&mut t0, 5); //   298: 319693579987813543907495863033445612882738907674745245368689798501688727968
        t0 *= t2; //                          299: 319693579987813543907495863033445612882738907674745245368689798501688727983
        square_assign_multi(&mut t0, 5); //   304: 10230194559610033405039867617070259612247645045591847851798073552054039295456
        t0 *= t1; //                          305: 10230194559610033405039867617070259612247645045591847851798073552054039295465

        CtOption::new(t0, !self.ct_eq(&Self::zero()))
    }

    /// Constructs a `Scalar` element without checking that it is
    /// canonical.
    pub const fn from_raw_unchecked(v: [u64; 4]) -> Self {
        Scalar(v)
    }

    /// Outputs a `Scalar` element of multiplicative order equals to 2^n
    pub fn get_root_of_unity(n: u32) -> Self {
        assert!(n != 0, "cannot get root of unity for n = 0");
        assert!(n <= TWO_ADICITY, "order cannot exceed 2^{}", TWO_ADICITY);
        let power = 1u64 << (TWO_ADICITY - n);

        TWO_ADIC_ROOT_OF_UNITY.exp(&[power, 0, 0, 0])
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
    /// the field modulus, modular reduction is silently performed.
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
    use rand_core::OsRng;

    const LARGEST: Scalar = Scalar([
        0x91776f9d741b35ea,
        0xa0a87b464e777106,
        0x9677863fe04d28d8,
        0x169e15bdd5149539,
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
            "0x073510d7d81d968688dd3b415caf3eb118c2b3faa0de24b7bfde343c02d4aee7"
        );
    }

    #[test]
    fn test_output_reduced_limbs() {
        assert_eq!(
            format!("{:?}", Scalar::zero().output_reduced_limbs()),
            "[0, 0, 0, 0]"
        );
        assert_eq!(
            format!("{:?}", Scalar::one().output_reduced_limbs()),
            "[1, 0, 0, 0]"
        );
        assert_eq!(
            format!("{:?}", R2.output_reduced_limbs()),
            "[13825545338424176359, 1784186291414246583, 9862103910925156017, 519339851260991110]"
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
                0x91776f9d741b35e9,
                0xa0a87b464e777106,
                0x9677863fe04d28d8,
                0x169e15bdd5149539,
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
        let mut rng = OsRng;

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

        let mut tmp = Scalar::random(&mut OsRng);

        for _ in 0..100 {
            let mut tmp2 = tmp.invert().unwrap();
            tmp2.mul_assign(&tmp);

            assert_eq!(tmp2, Scalar::one());

            tmp.add_assign(&Scalar::random(&mut OsRng));
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
    fn test_sqrt() {
        for _ in 0..100 {
            let a = Scalar::random(&mut OsRng).square();
            let b = a.sqrt().unwrap();
            assert_eq!(a, b.square());
        }

        assert_eq!(Scalar::zero().sqrt().unwrap(), Scalar::zero());
        assert_eq!(Scalar::one().sqrt().unwrap(), Scalar::one());

        // 3 is not a quadratic residue in Scalar
        assert!(bool::from(Scalar::new([3, 0, 0, 0]).sqrt().is_none()));
    }

    #[test]
    fn test_invert_is_pow() {
        let q_minus_2 = [
            0x91776f9d741b35e9,
            0xa0a87b464e777106,
            0x9677863fe04d28d8,
            0x169e15bdd5149539,
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

    // ROOTS OF UNITY
    // ================================================================================================

    #[test]
    fn test_get_root_of_unity() {
        let root_1 = Scalar::get_root_of_unity(1);
        assert_eq!(TWO_ADIC_ROOT_OF_UNITY, root_1);
        assert_eq!(Scalar::one(), root_1.exp(&[2, 0, 0, 0]));
    }

    #[test]
    fn test_lexicographically_largest() {
        // a = 4745189257869841308931709460885349048933591656118912969124572917019716039591
        let a = Scalar([
            0x1f4bddc6a36defc4,
            0x115fd9f0c807adaf,
            0xb3a334c8abab41d6,
            0x00de74ceb4eec9c4,
        ]);

        // b = 5485005301740192096108158156184910563314053389472934882673500635034323255876
        let b = Scalar([
            0x722b91d6d0ad4627,
            0x8f48a155866fc357,
            0xe2d4517734a1e702,
            0x15bfa0ef2025cb74,
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

        assert_eq!(element, Scalar::new([n as u64, 0, 0, 0]));
    }

    #[test]
    fn test_from_raw_unchecked() {
        let mut element = Scalar::from_raw_unchecked([
            0xbfde343c02d4aee7,
            0x18c2b3faa0de24b7,
            0x88dd3b415caf3eb1,
            0x073510d7d81d9686,
        ]);

        let element_normalized = Scalar::new([
            0xbfde343c02d4aee7,
            0x18c2b3faa0de24b7,
            0x88dd3b415caf3eb1,
            0x073510d7d81d9686,
        ]);

        assert_eq!(element, Scalar::one());
        element *= &R2;

        assert!(element != Scalar::one());
        assert_eq!(element, element_normalized);
    }

    // FIELD TRAIT
    // ================================================================================================

    #[test]
    fn test_field_trait_methods() {
        assert_eq!(<Scalar as Field>::zero(), Scalar::new([0, 0, 0, 0]));
        assert_eq!(<Scalar as Field>::one(), Scalar::new([1, 0, 0, 0]));

        assert_eq!(
            bool::from(<Scalar as Field>::zero().is_zero()),
            bool::from(Scalar::new([0, 0, 0, 0]).is_zero())
        );
        assert_eq!(
            bool::from(<Scalar as Field>::one().is_zero()),
            bool::from(Scalar::new([1, 0, 0, 0]).is_zero())
        );

        let mut rng = OsRng;
        let e = Scalar::random(&mut rng).square();

        assert_eq!(<Scalar as Field>::square(&e), e.square());
        assert_eq!(<Scalar as Field>::double(&e), e.double());

        assert_eq!(<Scalar as Field>::invert(&e).unwrap(), e.invert().unwrap());
        assert!(bool::from(
            <Scalar as Field>::invert(&Scalar::zero()).is_none()
        ));

        assert_eq!(<Scalar as Field>::sqrt(&e).unwrap(), e.sqrt().unwrap());
        assert!(bool::from(
            <Scalar as Field>::sqrt(&Scalar::new([3, 0, 0, 0])).is_none()
        ));
    }

    #[test]
    fn test_primefield_trait_methods() {
        assert_eq!(
            Scalar::from_repr([
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ])
            .unwrap(),
            Scalar::zero()
        );
        assert_eq!(
            Scalar::from_repr([
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ])
            .unwrap(),
            Scalar::one()
        );

        let mut rng = OsRng;
        let e = Scalar::random(&mut rng).square();

        assert_eq!(Scalar::from_repr(Scalar::to_repr(&e)).unwrap(), e);

        assert_eq!(
            Scalar::get_root_of_unity(<Scalar as PrimeField>::S),
            <Scalar as PrimeField>::root_of_unity()
        )
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
                231, 174, 212, 2, 60, 52, 222, 191, 183, 36, 222, 160, 250, 179, 194, 24, 177, 62,
                175, 92, 65, 59, 221, 136, 134, 150, 29, 216, 215, 16, 53, 7
            ]
        );

        assert_eq!(
            (-&Scalar::one()).to_bytes(),
            [
                234, 53, 27, 116, 157, 111, 119, 145, 6, 113, 119, 78, 70, 123, 168, 160, 216, 40,
                77, 224, 63, 134, 119, 150, 57, 149, 20, 213, 189, 21, 158, 22
            ]
        );
    }

    #[test]
    fn test_from_bytes() {
        let mut rng = OsRng;
        for _ in 0..100 {
            let a = Scalar::random(&mut rng);
            let bytes = a.to_bytes();
            assert_eq!(a, Scalar::from_bytes(&bytes).unwrap());
        }

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
                231, 174, 212, 2, 60, 52, 222, 191, 183, 36, 222, 160, 250, 179, 194, 24, 177, 62,
                175, 92, 65, 59, 221, 136, 134, 150, 29, 216, 215, 16, 53, 7
            ])
            .unwrap(),
            R2
        );

        // -1 should work
        assert_eq!(
            Scalar::from_bytes(&[
                234, 53, 27, 116, 157, 111, 119, 145, 6, 113, 119, 78, 70, 123, 168, 160, 216, 40,
                77, 224, 63, 134, 119, 150, 57, 149, 20, 213, 189, 21, 158, 22
            ])
            .unwrap(),
            -Scalar::one(),
        );

        // M is invalid
        assert!(bool::from(
            Scalar::from_bytes(&[
                235, 53, 27, 116, 157, 111, 119, 145, 6, 113, 119, 78, 70, 123, 168, 160, 216, 40,
                77, 224, 63, 134, 119, 150, 57, 149, 20, 213, 189, 21, 158, 22
            ])
            .is_none()
        ));

        // Anything larger than M is invalid
        assert!(bool::from(
            Scalar::from_bytes(&[
                236, 53, 27, 116, 157, 111, 119, 145, 6, 113, 119, 78, 70, 123, 168, 160, 216, 40,
                77, 224, 63, 134, 119, 150, 57, 149, 20, 213, 189, 21, 158, 22
            ])
            .is_none()
        ));
        assert!(bool::from(
            Scalar::from_bytes(&[
                234, 54, 27, 116, 157, 111, 119, 145, 6, 113, 119, 78, 70, 123, 168, 160, 216, 40,
                77, 224, 63, 134, 119, 150, 57, 149, 20, 213, 189, 21, 158, 22
            ])
            .is_none()
        ));
        assert!(bool::from(
            Scalar::from_bytes(&[
                0, 0, 0, 116, 157, 111, 119, 145, 6, 113, 119, 78, 70, 123, 168, 160, 216, 40, 77,
                224, 63, 134, 119, 150, 57, 149, 20, 213, 189, 255, 255, 255
            ])
            .is_none()
        ));
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
                231, 174, 212, 2, 60, 52, 222, 191, 183, 36, 222, 160, 250, 179, 194, 24, 177, 62,
                175, 92, 65, 59, 221, 136, 134, 150, 29, 216, 215, 16, 53, 7, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
        );
    }

    #[test]
    fn test_from_bytes_wide_negative_one() {
        assert_eq!(
            -&Scalar::one(),
            Scalar::from_bytes_wide(&[
                234, 53, 27, 116, 157, 111, 119, 145, 6, 113, 119, 78, 70, 123, 168, 160, 216, 40,
                77, 224, 63, 134, 119, 150, 57, 149, 20, 213, 189, 21, 158, 22, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            ])
        );
    }

    #[test]
    fn test_from_bytes_wide_maximum() {
        assert_eq!(
            Scalar([
                0xfe442038b6bef85e,
                0x73e43aa0772064ff,
                0x81729e4f73cd4a78,
                0x076aa97e99252434,
            ]),
            Scalar::from_bytes_wide(&[0xff; 64])
        );
    }

    #[test]
    fn test_from_bits() {
        let bytes = Scalar::zero().to_bytes();
        assert_eq!(Scalar::from_bits(bytes.as_bits::<Lsb0>()), Scalar::zero());

        let bytes = Scalar::one().to_bytes();
        assert_eq!(Scalar::from_bits(bytes.as_bits::<Lsb0>()), Scalar::one());

        let bytes = R2.to_bytes();
        assert_eq!(Scalar::from_bits(bytes.as_bits::<Lsb0>()), R2);

        // -1 should work
        let bytes = (-Scalar::one()).to_bytes();
        assert_eq!(Scalar::from_bits(bytes.as_bits::<Lsb0>()), -Scalar::one());

        // Modulus results in Scalar::zero()
        let bytes = [
            235, 53, 27, 116, 157, 111, 119, 145, 6, 113, 119, 78, 70, 123, 168, 160, 216, 40, 77,
            224, 63, 134, 119, 150, 57, 149, 20, 213, 189, 21, 158, 22,
        ];
        assert_eq!(Scalar::from_bits(bytes.as_bits::<Lsb0>()), Scalar::zero());
    }

    #[test]
    #[cfg(feature = "serialize")]
    fn test_serde() {
        let mut rng = OsRng;
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
