//! This crate provides an implementation of the STARK-friendly elliptic
//! curve E(Fp): y^2 = x^3 + x + B with
//! p = 2^251 + 17 * 2^192 + 1
//!   = 3618502788666131213697322783095070105623107215331596699973092056135872020481
//! B = 3141592653589793238462643383279502884197169399375105820974944592307816406665
//!

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

#[macro_use]
mod utils;

mod fp;
mod fp4;
mod scalar;

pub use scalar::Scalar;

pub use fp::Fp;
pub use fp4::Fp4;

// pub use curve::{AffinePoint, ProjectivePoint};
