// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::{BasePointTable, NafLookupTable, ProjectivePoint, Scalar};
use core::ops::Mul;

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
}

impl<'a, 'b> Mul<&'b Scalar> for &'a BASEPOINT_TABLE {
    type Output = ProjectivePoint;

    fn mul(self, scalar: &'b Scalar) -> Self::Output {
        self.multiply(&scalar.to_bytes())
    }
}

impl_binops_multiplicative_mixed!(BASEPOINT_TABLE, Scalar, ProjectivePoint);
