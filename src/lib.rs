// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This crate provides an implementation of the STARK-friendly elliptic
//! curve E(Fp6): y^2 = x^3 + x + B with
//! p = 2^64 - 2^32 + 1
//!   = 18446744069414584321
//! B = 8751648879042009540*u^5 + 9152663345541292493*u^4 + 11461655964319546493*u^3
//!         + 18139339604299112321*u^2 + 11484739320139982777*u + 13239264519837388909
//! where u^6 - 7 = 0
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
