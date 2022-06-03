// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! This module implements NAF lookup tables for points in Affine coordinates.

use crate::{AffinePoint, JacobianPoint, ProjectivePoint};

use zeroize::Zeroize;

/// A NAF lookup table for storing the precomputed values
/// `[P, 3P, 5P, ..., (2N-1)P]` in affine coordinates.
#[derive(Clone, Copy, Debug)]
pub struct NafLookupTable<const N: usize>(pub(crate) [AffinePoint; N]);

impl<const N: usize> From<ProjectivePoint> for NafLookupTable<N> {
    fn from(p: ProjectivePoint) -> Self {
        let mut points = [p; N];
        let double = p.double();
        for i in 1..N {
            points[i] = double + points[i - 1];
        }

        let mut points_affine = [AffinePoint::identity(); N];
        ProjectivePoint::batch_normalize(&points, &mut points_affine);

        Self(points_affine)
    }
}

impl<const N: usize> From<&ProjectivePoint> for NafLookupTable<N> {
    fn from(p: &ProjectivePoint) -> Self {
        let mut points = [*p; N];
        let double = p.double();
        for i in 1..N {
            points[i] = double + points[i - 1];
        }

        let mut points_affine = [AffinePoint::identity(); N];
        ProjectivePoint::batch_normalize(&points, &mut points_affine);

        Self(points_affine)
    }
}

impl<const N: usize> From<JacobianPoint> for NafLookupTable<N> {
    fn from(p: JacobianPoint) -> Self {
        let mut points = [p; N];
        let double = p.double();
        for i in 1..N {
            points[i] = double + points[i - 1];
        }

        let mut points_affine = [AffinePoint::identity(); N];
        JacobianPoint::batch_normalize(&points, &mut points_affine);

        Self(points_affine)
    }
}

impl<const N: usize> From<&JacobianPoint> for NafLookupTable<N> {
    fn from(p: &JacobianPoint) -> Self {
        let mut points = [*p; N];
        let double = p.double();
        for i in 1..N {
            points[i] = double + points[i - 1];
        }

        let mut points_affine = [AffinePoint::identity(); N];
        JacobianPoint::batch_normalize(&points, &mut points_affine);

        Self(points_affine)
    }
}

impl<const N: usize> From<AffinePoint> for NafLookupTable<N> {
    fn from(p: AffinePoint) -> Self {
        Self::from(ProjectivePoint::from(&p))
    }
}

impl<const N: usize> From<&AffinePoint> for NafLookupTable<N> {
    fn from(p: &AffinePoint) -> Self {
        Self::from(ProjectivePoint::from(p))
    }
}

impl<const N: usize> Zeroize for NafLookupTable<N> {
    fn zeroize(&mut self) {
        for point in self.0.iter_mut() {
            point.zeroize();
        }
    }
}

impl<const N: usize> NafLookupTable<N> {
    /// Given an `usize` x value, returns x.P.
    ///
    /// **NOTE**: The provided value x **MUST** be odd and smaller
    /// than 2N for a `NafLookupTable` of N elements.
    pub(crate) fn get_point(&self, x: usize) -> AffinePoint {
        debug_assert!(x & 1 == 1);
        debug_assert!(x < 2 * N);

        self.0[x / 2]
    }
}
