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

    /// An array of points corresponding to -2^i.SHIFT_POINT, for i ranging
    /// from 0 to 256 included. This is used when multiplying an `AffinePoint`,
    /// a `ProjectivePoint` or a `BasePointTable` with a `Scalar` for faster
    /// scalar multiplication.
    pub static ref MINUS_SHIFT_POINT_ARRAY: [AffinePoint; 257] =
        generate_opposite_powers(&SHIFT_POINT);
}

fn generate_opposite_powers(point: &ProjectivePoint) -> [AffinePoint; 257] {
    let mut array = [point.neg(); 257];
    for i in 1..257 {
        array[i] = array[i - 1].double();
    }

    let mut array_affine = [AffinePoint::identity(); 257];
    ProjectivePoint::batch_normalize(&array, &mut array_affine);

    array_affine
}

impl<'a, 'b> Mul<&'b Scalar> for &'a BASEPOINT_TABLE {
    type Output = ProjectivePoint;

    fn mul(self, scalar: &'b Scalar) -> Self::Output {
        self.multiply(&scalar.to_bytes())
    }
}

impl_binops_multiplicative_mixed!(BASEPOINT_TABLE, Scalar, ProjectivePoint);
