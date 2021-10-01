//! This crate provides an implementation of the STARK-friendly elliptic
//! curve E(Fp4): y^2 = x^3 + x + B with
//! p = 2^62 - 111 * 2^39 + 1
//!   = 4611624995532046337
//! B = 2708278037369052277*u^3 + 489710895200713542*u^2 + 3456610074177457817*u + 1669244588749562658
//! where u^4 + 3 = 0.
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

// mod curve;
mod fp;
mod fp2;
mod fp4;
mod fp6;
mod scalar;

pub use scalar::Scalar;

pub use fp::Fp;
pub use fp2::Fp2;
pub use fp4::Fp4;
pub use fp6::Fp6;

// pub use curve::{AffinePoint, ProjectivePoint};
