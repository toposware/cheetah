//! This crate provides an implementation of the STARK-friendly elliptic
//! curve E(Fp6): y^2 = x^3 + x + B with
//! p = 2^62 - 111 * 2^39 + 1
//!   = 4611624995532046337
//! B = (4198454498232167043*u + 1236584124511164798)*v^2 + (109675342022917428*u + 2744078705171299893)*v + 4212198449936917436*u + 2883908786857436727
//! where u^2 - u -1 = 0
//! and v^3 - v - 2 = 0.
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

// Re-export group (which re-exports ff)
pub use group;

#[macro_use]
mod utils;

mod curve;
mod fp;
mod fp2;
mod fp6;
mod scalar;

pub use scalar::Scalar;

pub use fp::Fp;
pub use fp2::Fp2;
pub use fp6::Fp6;

pub use curve::{AffinePoint, ProjectivePoint, B};
