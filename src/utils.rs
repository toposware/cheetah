// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use group::ff::Field;

/// Compute a + b + carry, returning the result and the new carry over.
#[inline(always)]
pub const fn add64_with_carry(a: u64, b: u64, carry: u64) -> (u64, u64) {
    let ret = (a as u128) + (b as u128) + (carry as u128);
    (ret as u64, (ret >> 64) as u64)
}

/// Compute a - (b + borrow), returning the result and the new borrow.
#[inline(always)]
pub const fn sub64_with_carry(a: u64, b: u64, borrow: u64) -> (u64, u64) {
    let ret = (a as u128).wrapping_sub((b as u128) + ((borrow >> 63) as u128));
    (ret as u64, (ret >> 64) as u64)
}

/// Compute a + (b * c) + carry, returning the result and the new carry over.
#[inline(always)]
pub const fn mul64_with_carry(a: u64, b: u64, c: u64, carry: u64) -> (u64, u64) {
    let ret = (a as u128) + ((b as u128) * (c as u128)) + (carry as u128);
    (ret as u64, (ret >> 64) as u64)
}

/// Compute (a << b) + carry, returning the result and the new carry over.
#[inline(always)]
pub const fn shl64_by_u32_with_carry(a: u64, b: u32, carry: u64) -> (u64, u64) {
    let ret = ((a as u128) << (b as u128)) + (carry as u128);
    (ret as u64, (ret >> 64) as u64)
}

#[inline(always)]
pub(crate) fn square_assign_multi<F: Field>(n: &mut F, num_times: usize) {
    for _ in 0..num_times {
        *n = n.square();
    }
}

macro_rules! impl_add_binop_specify_output {
    ($lhs:ident, $rhs:ident, $output:ident) => {
        impl<'b> Add<&'b $rhs> for $lhs {
            type Output = $output;

            #[inline]
            fn add(self, rhs: &'b $rhs) -> $output {
                &self + rhs
            }
        }

        impl<'a> Add<$rhs> for &'a $lhs {
            type Output = $output;

            #[inline]
            fn add(self, rhs: $rhs) -> $output {
                self + &rhs
            }
        }

        impl Add<$rhs> for $lhs {
            type Output = $output;

            #[inline]
            fn add(self, rhs: $rhs) -> $output {
                &self + &rhs
            }
        }
    };
}

macro_rules! impl_sub_binop_specify_output {
    ($lhs:ident, $rhs:ident, $output:ident) => {
        impl<'b> Sub<&'b $rhs> for $lhs {
            type Output = $output;

            #[inline]
            fn sub(self, rhs: &'b $rhs) -> $output {
                &self - rhs
            }
        }

        impl<'a> Sub<$rhs> for &'a $lhs {
            type Output = $output;

            #[inline]
            fn sub(self, rhs: $rhs) -> $output {
                self - &rhs
            }
        }

        impl Sub<$rhs> for $lhs {
            type Output = $output;

            #[inline]
            fn sub(self, rhs: $rhs) -> $output {
                &self - &rhs
            }
        }
    };
}

macro_rules! impl_binops_additive_specify_output {
    ($lhs:ident, $rhs:ident, $output:ident) => {
        impl_add_binop_specify_output!($lhs, $rhs, $output);
        impl_sub_binop_specify_output!($lhs, $rhs, $output);
    };
}

macro_rules! impl_binops_multiplicative_mixed {
    ($lhs:ident, $rhs:ident, $output:ident) => {
        impl<'b> Mul<&'b $rhs> for $lhs {
            type Output = $output;

            #[inline]
            fn mul(self, rhs: &'b $rhs) -> $output {
                &self * rhs
            }
        }

        impl<'a> Mul<$rhs> for &'a $lhs {
            type Output = $output;

            #[inline]
            fn mul(self, rhs: $rhs) -> $output {
                self * &rhs
            }
        }

        impl Mul<$rhs> for $lhs {
            type Output = $output;

            #[inline]
            fn mul(self, rhs: $rhs) -> $output {
                &self * &rhs
            }
        }
    };
}

macro_rules! impl_binops_additive {
    ($lhs:ident, $rhs:ident) => {
        impl_binops_additive_specify_output!($lhs, $rhs, $lhs);

        impl SubAssign<$rhs> for $lhs {
            #[inline]
            fn sub_assign(&mut self, rhs: $rhs) {
                *self = &*self - &rhs;
            }
        }

        impl AddAssign<$rhs> for $lhs {
            #[inline]
            fn add_assign(&mut self, rhs: $rhs) {
                *self = &*self + &rhs;
            }
        }

        impl<'b> SubAssign<&'b $rhs> for $lhs {
            #[inline]
            fn sub_assign(&mut self, rhs: &'b $rhs) {
                *self = &*self - rhs;
            }
        }

        impl<'b> AddAssign<&'b $rhs> for $lhs {
            #[inline]
            fn add_assign(&mut self, rhs: &'b $rhs) {
                *self = &*self + rhs;
            }
        }
    };
}

macro_rules! impl_binops_multiplicative {
    ($lhs:ident, $rhs:ident) => {
        impl_binops_multiplicative_mixed!($lhs, $rhs, $lhs);

        impl MulAssign<$rhs> for $lhs {
            #[inline]
            fn mul_assign(&mut self, rhs: $rhs) {
                *self = &*self * &rhs;
            }
        }

        impl<'b> MulAssign<&'b $rhs> for $lhs {
            #[inline]
            fn mul_assign(&mut self, rhs: &'b $rhs) {
                *self = &*self * rhs;
            }
        }
    };
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_add_with_carry() {
        use crate::utils::add64_with_carry;
        {
            let a = 0u64;
            let b = 1u64;
            let (c, carry) = add64_with_carry(a, b, 0);

            assert_eq!(b, c);
            assert_eq!(carry, 0);
        }
        {
            let a = u64::MAX;
            let b = 1u64;
            let (c, carry) = add64_with_carry(a, b, 1);

            assert_eq!(c, b);
            assert_eq!(carry, 1);
        }
        {
            let a = u64::MAX;
            let b = 1u64;
            let (c, carry) = add64_with_carry(a, b, u64::MAX);

            assert_eq!(c, u64::MAX);
            assert_eq!(carry, 1);
        }
    }

    #[test]
    fn test_sub_with_carry() {
        use crate::utils::sub64_with_carry;
        {
            let a = 1u64;
            let b = 0u64;
            let (c, carry) = sub64_with_carry(a, b, 0);

            assert_eq!(a, c);
            assert_eq!(carry, 0);
        }
        {
            let a = u64::MAX;
            let b = u64::MAX;
            let (c, carry) = sub64_with_carry(a, b, 1);

            assert_eq!(c, 0);
            assert_eq!(carry, 0);
        }
        {
            let a = 1u64;
            let b = u64::MAX;
            let (c, carry) = sub64_with_carry(a, b, 0);

            assert_eq!(c, 2);
            assert_eq!(carry, u64::MAX);
        }
    }

    #[test]
    fn test_mul_with_carry() {
        use crate::utils::mul64_with_carry;
        {
            let a = 0u64;
            let b = 1u64;
            let c = 1u64;
            let (d, carry) = mul64_with_carry(a, b, c, 0);

            assert_eq!(b, d);
            assert_eq!(carry, 0);
        }
        {
            let a = u64::MAX;
            let b = 1u64;
            let c = 1u64;
            let (d, carry) = mul64_with_carry(a, b, c, 0);

            assert_eq!(d, 0);
            assert_eq!(carry, 1);
        }
        {
            let a = u64::MAX;
            let b = 42u64;
            let c = 1u64;
            let (d, carry) = mul64_with_carry(a, b, c, 1);

            assert_eq!(d, b);
            assert_eq!(carry, 1);
        }
    }

    #[test]
    fn test_shl_with_carry() {
        use crate::utils::shl64_by_u32_with_carry;
        {
            let a = 1u64;
            let b = 1u32;
            let c = 1u64;
            let (d, carry) = shl64_by_u32_with_carry(a, b, c);

            assert_eq!(d, 3);
            assert_eq!(carry, 0);
        }
        {
            let a = u64::MAX;
            let b = 1u32;
            let c = 1u64;
            let (d, carry) = shl64_by_u32_with_carry(a, b, c);

            assert_eq!(d, u64::MAX);
            assert_eq!(carry, 1);
        }
        {
            let a = u64::MAX;
            let b = 64u32;
            let c = 1u64;
            let (d, carry) = shl64_by_u32_with_carry(a, b, c);

            assert_eq!(d, 1);
            assert_eq!(carry, u64::MAX);
        }
    }

    #[test]
    fn test_square_assign_multi() {
        use crate::fp::Fp;
        use crate::utils::square_assign_multi;
        use group::ff::Field;
        use rand_core::OsRng;

        let mut rng = OsRng;
        {
            let mut e = Fp::random(&mut rng);
            let e_copy = e;
            square_assign_multi(&mut e, 1);
            assert_eq!(e, e_copy.square());
        }
        {
            let mut e = Fp::random(&mut rng);
            let e_copy = e;
            square_assign_multi(&mut e, 2);
            assert_eq!(e, e_copy.square().square());
        }
    }
}
