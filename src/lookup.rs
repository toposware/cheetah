// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This module implements lookup tables for points in Projective coordinates.
//!
//! Adapted from https://github.com/RustCrypto/elliptic-curves

use crate::{AffinePoint, ProjectivePoint, Scalar};

use core::ops::Mul;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};
use zeroize::Zeroize;

/// A lookup table for storing the precomputed values `[p, 2p, 3p, ..., Np]`
/// in projective coordinates.
#[derive(Clone, Copy, Debug)]
pub struct LookupTable<const N: usize>(pub(crate) [ProjectivePoint; N]);

impl<const N: usize> Default for LookupTable<N> {
    fn default() -> Self {
        Self([ProjectivePoint::default(); N])
    }
}

impl<const N: usize> From<ProjectivePoint> for LookupTable<N> {
    fn from(p: ProjectivePoint) -> Self {
        let mut points = [p; N];
        for i in 1..N {
            points[i] = p + points[i - 1];
        }

        Self(points)
    }
}

impl<const N: usize> From<&ProjectivePoint> for LookupTable<N> {
    fn from(p: &ProjectivePoint) -> Self {
        let mut points = [*p; N];
        for i in 1..N {
            points[i] = p + points[i - 1];
        }

        Self(points)
    }
}

impl<const N: usize> From<AffinePoint> for LookupTable<N> {
    fn from(p: AffinePoint) -> Self {
        let mut points = [ProjectivePoint::from(&p); N];
        for i in 1..N {
            points[i] = p + points[i - 1];
        }

        Self(points)
    }
}

impl<const N: usize> From<&AffinePoint> for LookupTable<N> {
    fn from(p: &AffinePoint) -> Self {
        let mut points = [ProjectivePoint::from(p); N];
        for i in 1..N {
            points[i] = p + points[i - 1];
        }

        Self(points)
    }
}

impl<const N: usize> Zeroize for LookupTable<N> {
    fn zeroize(&mut self) {
        for point in self.0.iter_mut() {
            point.zeroize();
        }
    }
}

impl<const N: usize> LookupTable<N> {
    /// Given an `i8` x value, returns x.P.
    // To do so, we first compute |x|.P, and then conditionally
    // negate the result based on the sign of x.
    pub(crate) fn get_point(&self, x: i8) -> ProjectivePoint {
        debug_assert!(x >= -(N as i8));
        debug_assert!(x <= N as i8);

        // Compute xabs = |x|
        let xmask = x >> 7;
        let xabs = (x + xmask) ^ xmask;

        // Get an array element in constant time
        let mut t = ProjectivePoint::identity();
        for j in 1..N + 1 {
            let c = (xabs as u8).ct_eq(&(j as u8));
            t.conditional_assign(&self.0[j - 1], c);
        }
        // Now t == |x| * p.

        let neg_mask = Choice::from((xmask & 1) as u8);
        t.conditional_assign(&-t, neg_mask);
        // Now t == x * p.

        t
    }

    /// Given an `i8` x value, returns x.P.
    // To do so, we first compute |x|.P, and then conditionally
    // negate the result based on the sign of x.
    pub(crate) fn get_point_vartime(&self, x: i8) -> ProjectivePoint {
        debug_assert!(x >= -(N as i8));
        debug_assert!(x <= N as i8);

        // Compute xabs = |x|
        let xmask = x >> 7;
        let xabs = (x + xmask) ^ xmask;

        // Get an array element
        let mut t = ProjectivePoint::identity();
        for j in 1..N + 1 {
            if xabs as usize == j {
                t = self.0[j - 1];
            }
        }
        // Now t == |x| * p.

        if xmask & 1 == 1 {
            t = -t;
        }
        // Now t == x * p.

        t
    }
}

/// A list of `LookupTable` of multiples of a point `p` to perform
/// efficient scalar multiplication with the Pippenger's algorith,
/// https://cr.yp.to/papers/pippenger.pdfm.
///
/// Creating a base point table is costly, and hence should be
/// done in the purpose of being used more than once.
#[derive(Clone, Debug)]
pub struct BasePointTable(pub [LookupTable<8>; 32]);

impl BasePointTable {
    /// Returns a precomputed table of multiples of a given point.
    pub fn create(basepoint: &ProjectivePoint) -> Self {
        let mut table = BasePointTable([LookupTable::default(); 32]);
        let mut point = *basepoint;
        for i in 0..32 {
            table.0[i] = LookupTable::from(&point);
            point = point.double_multi(8);
        }

        table
    }

    /// Get the basepoint of this table.
    pub fn get_basepoint(&self) -> ProjectivePoint {
        self.0[0].get_point(1)
    }

    /// Get the basepoint of this table.
    pub fn get_basepoint_vartime(&self) -> ProjectivePoint {
        self.0[0].get_point_vartime(1)
    }

    /// Performs a projective scalar multiplication from `by`
    /// given as byte representation of a `Scalar` element, by
    /// using internally the Pippenger's algorithm.
    #[inline]
    pub fn multiply(&self, scalar: &[u8; 32]) -> ProjectivePoint {
        let a = Scalar::bytes_to_radix_16(scalar);

        let tables = &self.0;
        let mut acc = ProjectivePoint::identity();

        for i in (0..64).filter(|x| x % 2 == 1) {
            acc += &tables[i / 2].get_point(a[i]);
        }

        acc = acc.double_multi(4);

        for i in (0..64).filter(|x| x % 2 == 0) {
            acc += &tables[i / 2].get_point(a[i]);
        }

        acc
    }

    /// Performs a projective scalar multiplication from `by`
    /// given as byte representation of a `Scalar` element, by
    /// using internally the Pippenger's algorithm.
    ///
    /// **This operation is variable time with respect
    /// to the scalar.** If the scalar is fixed,
    /// this operation is effectively constant time.
    #[inline]
    pub fn multiply_vartime(&self, scalar: &[u8; 32]) -> ProjectivePoint {
        let a = Scalar::bytes_to_radix_16(scalar);

        let tables = &self.0;
        let mut acc = ProjectivePoint::identity();

        for i in (0..64).filter(|x| x % 2 == 1) {
            acc += &tables[i / 2].get_point_vartime(a[i]);
        }

        acc = acc.double_multi(4);

        for i in (0..64).filter(|x| x % 2 == 0) {
            acc += &tables[i / 2].get_point_vartime(a[i]);
        }

        acc
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a BasePointTable {
    type Output = ProjectivePoint;

    fn mul(self, scalar: &'b Scalar) -> Self::Output {
        self.multiply(&scalar.to_bytes())
    }
}

impl<'a> Mul<Scalar> for &'a BasePointTable {
    type Output = ProjectivePoint;

    #[inline]
    fn mul(self, scalar: Scalar) -> ProjectivePoint {
        self.multiply(&scalar.to_bytes())
    }
}
