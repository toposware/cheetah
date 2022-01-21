// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This crate provides an implementation of the STARK-friendly elliptic
//! curve E(Fp6): y^2 = x^3 + x + B with
//! p = 2^62 + 2^56 + 2^55 + 1
//!   = 4719772409484279809
//! B = (1200866201009650596 * u + 1935817186716799185) * v^2
//!     + (3999205700308519553 * u + 3518137720867787056) * v
//!     + 2508413708960025374 * u + 1526905369741321712
//! where u^2 - 2u - 2 = 0
//! and v^3 + v + 1 = 0.
//!
//! Portions of the codebase have been adapted from https://github.com/zkcrypto/bls12_381

#![no_std]
#![cfg_attr(docsrs, feature(doc_cfg))]
// Catch documentation errors caused by code changes.
#![deny(broken_intra_doc_links)]
#![deny(missing_debug_implementations)]
#![deny(missing_docs)]
#![deny(unsafe_code)]
#![allow(clippy::many_single_char_names)]

#[cfg(test)]
#[macro_use]
extern crate std;

// Re-export group (which re-exports ff)
pub use group;

#[macro_use]
mod utils;

mod constants;
mod lookup;

mod curve;
mod fp;
mod fp2;
mod fp6;
mod scalar;

pub use scalar::Scalar;

pub use fp::Fp;
pub use fp2::Fp2;
pub use fp6::Fp6;

pub use constants::BASEPOINT_TABLE;
pub use lookup::{BasePointTable, LookupTable};

pub use curve::{AffinePoint, ProjectivePoint, B};
