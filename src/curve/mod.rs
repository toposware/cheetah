// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This module provides an implementation of the STARK-friendly
//! cheetah curve over the sextic extension of the prime field Fp
//! of characteristic p = 2^64 - 2^32 + 1.

mod affine;
mod encoding;
mod jacobian;
mod projective;

pub use affine::AffinePoint;
pub use encoding::{CompressedPoint, UncompressedPoint};
pub use jacobian::{JacobianPoint, ModifiedJacobianPoint};
pub use projective::ProjectivePoint;

// HELPER METHODS
// ================================================================================================

use crate::fp::{reduce_u96, GENERATOR};
use crate::{Fp, Fp6};

//  A = 1
/// B = u + 395
pub const B: Fp6 = Fp6 {
    c0: Fp(395),
    c1: Fp::one(),
    c2: Fp::zero(),
    c3: Fp::zero(),
    c4: Fp::zero(),
    c5: Fp::zero(),
};

pub(crate) const B3: Fp6 = (&B).mul(&Fp6::new([3, 0, 0, 0, 0, 0]));

#[inline(always)]
pub(crate) const fn mul_by_3b(a: &Fp6) -> Fp6 {
    let aa = (&a.c0).mul(&B3.c0).0 as u128;
    let ab = a.c0.0 as u128;
    let ab = ab + (ab << 1);

    let ba = (&a.c1).mul(&B3.c0).0 as u128;
    let bb = a.c1.0 as u128;
    let bb = bb + (bb << 1);

    let ca = (&a.c2).mul(&B3.c0).0 as u128;
    let cb = a.c2.0 as u128;
    let cb = cb + (cb << 1);

    let da = (&a.c3).mul(&B3.c0).0 as u128;
    let db = a.c3.0 as u128;
    let db = db + (db << 1);

    let ea = (&a.c4).mul(&B3.c0).0 as u128;
    let eb = a.c4.0 as u128;
    let eb = eb + (eb << 1);

    let fa = (&a.c5).mul(&B3.c0).0 as u128;
    let fb = a.c5.0 as u128;
    let fb = fb + (fb << 1);

    let c0 = fb * GENERATOR.0 as u128;
    let c0 = c0 + aa;
    let c0 = Fp(reduce_u96(c0));

    let c1 = ab + ba;
    let c1 = Fp(reduce_u96(c1));

    let c2 = bb + ca;
    let c2 = Fp(reduce_u96(c2));

    let c3 = cb + da;
    let c3 = Fp(reduce_u96(c3));

    let c4 = db + ea;
    let c4 = Fp(reduce_u96(c4));

    let c5 = eb + fa;
    let c5 = Fp(reduce_u96(c5));

    Fp6 {
        c0,
        c1,
        c2,
        c3,
        c4,
        c5,
    }
}
