// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::{AffinePoint, BasePointTable, NafLookupTable, ProjectivePoint};
use crate::{Fp, Fp6, Scalar};

use core::ops::Mul;
use subtle::Choice;

/// A fixed shift point of the curve in projective coordinates
/// The point has been generated from the Simplified Shallue-van
/// de Woestijne-Ulas method for hashing to elliptic curves in
/// Short Weierstrass form, applied on the integer value of the
/// binary encoding of the string "Cheetah - Shift point".
///
/// See https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve
/// for more.
pub const SHIFT_POINT: ProjectivePoint = ProjectivePoint {
    x: Fp6 {
        c0: Fp(0xab95fb0d188ff2b1),
        c1: Fp(0xf38263a8c954493e),
        c2: Fp(0x1bc968d6739028b7),
        c3: Fp(0x4f67684fe6b27d1f),
        c4: Fp(0x6504f8bbba9077b4),
        c5: Fp(0x860838116d76e545),
    },
    y: Fp6 {
        c0: Fp(0x56988918d5b02cf3),
        c1: Fp(0xc9fe3ef9d7652c54),
        c2: Fp(0x4923e90c4d6f4538),
        c3: Fp(0xc747ac3f82369f27),
        c4: Fp(0x6e74c106972c71f5),
        c5: Fp(0xcf190938f1c64bb5),
    },
    z: Fp6::one(),
};

lazy_static! {
    /// A hardcoded `BasePointTable` for the generator of the Cheetah
    /// curve, to allow for efficient single scalar multiplication.
    pub static ref BASEPOINT_TABLE: BasePointTable =
        BasePointTable::create(&ProjectivePoint::generator());

    /// A hardcoded `NafLookupTable` for odd multiples of the generator
    /// of the Cheetah curve, to allow for efficient double scalar
    /// multiplications.
    pub static ref ODD_MULTIPLES_BASEPOINT: NafLookupTable::<64> =
        NafLookupTable::<64>::from(&ProjectivePoint::generator());

    /// The opposite of 2^4.SHIFT_POINT. This is used when multiplying a
    /// `BasePointTable` with a `Scalar` for faster scalar multiplication.
    pub static ref MINUS_SHIFT_POINT_POW_4: AffinePoint = AffinePoint {
        x: Fp6 {
            c0: Fp(0x59f135aa8369eea0),
            c1: Fp(0xdfcadde02ca0362b),
            c2: Fp(0x99e498bacad59dc0),
            c3: Fp(0xc0ae607c0e4289fe),
            c4: Fp(0xa3c9733bbefe411c),
            c5: Fp(0x2c0e26555061fbbf),
            },
        y: Fp6 {
            c0: Fp(0xfb19a0cfab8e245c),
            c1: Fp(0x23f1f9b307ac72e1),
            c2: Fp(0xb9d43296c72c79f2),
            c3: Fp(0xa430f22e5903080c),
            c4: Fp(0x2ca860527297557e),
            c5: Fp(0x9e05765dcd66172c),
        },
        infinity: Choice::from(0u8),
    };

    /// The opposite of 2^256.SHIFT_POINT. This is used when multiplying
    /// an `AffinePoint` or a `ProjectivePoint` with a `Scalar` for faster
    /// scalar multiplication.
    pub static ref MINUS_SHIFT_POINT_POW_256: AffinePoint = AffinePoint {
        x: Fp6 {
            c0: Fp(0xbf055af286583ccc),
            c1: Fp(0xee33e2330c725cdd),
            c2: Fp(0xa1bdadd9ea138f24),
            c3: Fp(0xf97d628d13c4e145),
            c4: Fp(0xb6389323811c4747),
            c5: Fp(0xed92b921e0c7e575),
        },
        y: Fp6 {
            c0: Fp(0xcad959f8181d46a4),
            c1: Fp(0x2d434d5931a7c485),
            c2: Fp(0x39483eba519e8d5f),
            c3: Fp(0xeadcd00fc1d23ee),
            c4: Fp(0x6f5350dfabc6b254),
            c5: Fp(0x286e25dbcc086611),
        },
        infinity: Choice::from(0u8),
    };
}

impl<'a, 'b> Mul<&'b Scalar> for &'a BASEPOINT_TABLE {
    type Output = ProjectivePoint;

    fn mul(self, scalar: &'b Scalar) -> Self::Output {
        self.multiply(&scalar.to_bytes())
    }
}

impl_binops_multiplicative_mixed!(BASEPOINT_TABLE, Scalar, ProjectivePoint);
