use core::{
    borrow::Borrow,
    fmt,
    iter::Sum,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::{fp::Fp, fp4::Fp4, scalar::Scalar};

use rand_core::{CryptoRng, RngCore};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

impl_binops_additive!(ProjectivePoint, AffinePoint);
impl_binops_additive_specify_output!(AffinePoint, ProjectivePoint, ProjectivePoint);

// A = 1
// B = 2708278037369052277*u^3 + 489710895200713542*u^2 + 3456610074177457817*u + 1669244588749562658
pub const B: Fp4 = Fp4 {
    c0: Fp::new(1669244588749562658),
    c1: Fp::new(3456610074177457817),
    c2: Fp::new(489710895200713542),
    c3: Fp::new(2708278037369052277),
};

/// An affine point
#[derive(Copy, Clone, Debug)]
pub struct AffinePoint {
    pub(crate) x: Fp4,
    pub(crate) y: Fp4,
    infinity: Choice,
}

impl Default for AffinePoint {
    fn default() -> AffinePoint {
        AffinePoint::identity()
    }
}

impl zeroize::DefaultIsZeroes for AffinePoint {}

impl fmt::Display for AffinePoint {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl<'a> From<&'a ProjectivePoint> for AffinePoint {
    fn from(p: &'a ProjectivePoint) -> AffinePoint {
        let zinv = p.z.invert().unwrap_or(Fp4::zero());
        let x = p.x * zinv;
        let y = p.y * zinv;

        let tmp = AffinePoint {
            x,
            y,
            infinity: Choice::from(0u8),
        };

        if zinv == Fp4::zero() {
            AffinePoint::identity()
        } else {
            tmp
        }
    }
}

impl From<ProjectivePoint> for AffinePoint {
    fn from(p: ProjectivePoint) -> AffinePoint {
        AffinePoint::from(&p)
    }
}

impl ConstantTimeEq for AffinePoint {
    fn ct_eq(&self, other: &Self) -> Choice {
        // The only cases in which two points are equal are
        // 1. infinity is set on both
        // 2. infinity is not set on both, and their coordinates are equal

        (self.infinity & other.infinity)
            | ((!self.infinity)
                & (!other.infinity)
                & self.x.ct_eq(&other.x)
                & self.y.ct_eq(&other.y))
    }
}

impl ConditionallySelectable for AffinePoint {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        AffinePoint {
            x: Fp4::conditional_select(&a.x, &b.x, choice),
            y: Fp4::conditional_select(&a.y, &b.y, choice),
            infinity: Choice::conditional_select(&a.infinity, &b.infinity, choice),
        }
    }
}

impl Eq for AffinePoint {}
impl PartialEq for AffinePoint {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl<'a> Neg for &'a AffinePoint {
    type Output = AffinePoint;

    #[inline]
    fn neg(self) -> AffinePoint {
        AffinePoint {
            x: self.x,
            y: Fp4::conditional_select(&-self.y, &Fp4::one(), self.infinity),
            infinity: self.infinity,
        }
    }
}

impl Neg for AffinePoint {
    type Output = AffinePoint;

    #[inline]
    fn neg(self) -> AffinePoint {
        -&self
    }
}

impl<'a, 'b> Add<&'b ProjectivePoint> for &'a AffinePoint {
    type Output = ProjectivePoint;

    #[inline]
    fn add(self, rhs: &'b ProjectivePoint) -> ProjectivePoint {
        rhs.add_mixed(self)
    }
}

impl<'a, 'b> Add<&'b AffinePoint> for &'a ProjectivePoint {
    type Output = ProjectivePoint;

    #[inline]
    fn add(self, rhs: &'b AffinePoint) -> ProjectivePoint {
        self.add_mixed(rhs)
    }
}

impl<'a, 'b> Sub<&'b ProjectivePoint> for &'a AffinePoint {
    type Output = ProjectivePoint;

    #[inline]
    fn sub(self, rhs: &'b ProjectivePoint) -> ProjectivePoint {
        self + (-rhs)
    }
}

impl<'a, 'b> Sub<&'b AffinePoint> for &'a ProjectivePoint {
    type Output = ProjectivePoint;

    #[inline]
    fn sub(self, rhs: &'b AffinePoint) -> ProjectivePoint {
        self + (-rhs)
    }
}

impl<T> Sum<T> for ProjectivePoint
where
    T: Borrow<ProjectivePoint>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Self::identity(), |acc, item| acc + item.borrow())
    }
}

impl AffinePoint {
    /// Returns the identity point in affine coordinates
    pub fn identity() -> AffinePoint {
        AffinePoint {
            x: Fp4::zero(),
            y: Fp4::one(),
            infinity: Choice::from(1u8),
        }
    }

    /// Returns the x coordinate of this AffinePoint
    pub fn get_x(&self) -> Fp4 {
        self.x
    }

    /// Returns the y coordinate of this AffinePoint
    pub fn get_y(&self) -> Fp4 {
        self.y
    }

    /// Computes a random `AffinePoint` element
    pub fn random(mut rng: impl RngCore + CryptoRng) -> Self {
        loop {
            let x = Fp4::random(&mut rng);
            let flip_sign = rng.next_u32() % 2 != 0;

            // Obtain the corresponding y-coordinate given x as y = sqrt(x^3 + x + B)
            if let Some(p) = ((x.square() * x) + x + B)
                .sqrt_vartime()
                .map(|y| AffinePoint {
                    x,
                    y: if flip_sign { -y } else { y },
                    infinity: 0.into(),
                })
            {
                if bool::from(!p.is_identity()) {
                    return p;
                }
            }
        }
    }

    /// Returns a fixed generator of the curve in affine coordinates
    pub fn generator() -> AffinePoint {
        AffinePoint {
            x: Fp4 {
                c0: Fp::from_raw_unchecked(225412627959290117),
                c1: Fp::from_raw_unchecked(498792661968344066),
                c2: Fp::from_raw_unchecked(3164519460772836424),
                c3: Fp::from_raw_unchecked(3347389561768906505),
            },
            y: Fp4 {
                c0: Fp::from_raw_unchecked(3607967135093124863),
                c1: Fp::from_raw_unchecked(4039667423912778261),
                c2: Fp::from_raw_unchecked(3737455456557202021),
                c3: Fp::from_raw_unchecked(1753374771055551422),
            },
            infinity: Choice::from(0u8),
        }
    }

    /// Outputs a compress byte representation of this `AffinePoint` element
    pub fn to_compressed(&self) -> [u8; 32] {
        // Strictly speaking, self.x is zero already when self.infinity is true, but
        // to guard against implementation mistakes we do not assume this.
        let mut res = Fp4::conditional_select(&self.x, &Fp4::zero(), self.infinity).to_bytes();

        // This point is in compressed form, so we set the most significant bit.
        res[31] |= 1u8 << 7;

        // Is the y-coordinate the lexicographically largest of the two associated with the
        // x-coordinate? If so, set the third-most significant bit so long as this is not
        // the point at infinity.
        res[31] |= u8::conditional_select(
            &0u8,
            &(1u8 << 6),
            (!self.infinity) & self.y.lexicographically_largest(),
        );

        res
    }

    /// Outputs an uncompressed byte representation of this `AffinePoint` element
    /// It is twice larger than when calling `AffinePoint::to_compressed()`
    pub fn to_uncompressed(&self) -> [u8; 64] {
        let mut res = [0; 64];

        res[0..32].copy_from_slice(
            &Fp4::conditional_select(&self.x, &Fp4::zero(), self.infinity).to_bytes()[..],
        );
        res[32..64].copy_from_slice(
            &Fp4::conditional_select(&self.y, &Fp4::zero(), self.infinity).to_bytes()[..],
        );

        res
    }

    /// Attempts to deserialize an uncompressed element.
    pub fn from_uncompressed(bytes: &[u8; 64]) -> CtOption<Self> {
        Self::from_uncompressed_unchecked(bytes).and_then(|p| CtOption::new(p, p.is_on_curve()))
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is on the curve and not checking if it is in the correct subgroup.
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_uncompressed()` instead.
    pub fn from_uncompressed_unchecked(bytes: &[u8; 64]) -> CtOption<Self> {
        // Obtain the three flags from the start of the byte sequence
        let compression_flag_set = Choice::from((bytes[63] >> 7) & 1);
        let sort_flag_set = Choice::from((bytes[63] >> 6) & 1);

        // Attempt to obtain the x-coordinate
        let x = {
            let mut tmp = [0; 32];
            tmp.copy_from_slice(&bytes[0..32]);

            Fp4::from_bytes(&tmp)
        };

        // Attempt to obtain the y-coordinate
        let y = {
            let mut tmp = [0; 32];
            tmp.copy_from_slice(&bytes[32..64]);

            // Mask away the flag bits
            tmp[31] &= 0b0011_1111;

            Fp4::from_bytes(&tmp)
        };

        let infinity_flag = x.is_zero() & y.is_zero();
        // Create a point representing this value
        let p = AffinePoint::conditional_select(
            &AffinePoint {
                x,
                y,
                infinity: infinity_flag,
            },
            &AffinePoint::identity(),
            infinity_flag,
        );

        CtOption::new(
            p,
            // The compression flag should not have been set, as this is an uncompressed element
            (!compression_flag_set) &
            // The sort flag should not have been set, as this is an uncompressed element
            (!sort_flag_set),
        )
    }

    /// Attempts to deserialize a compressed element.
    pub fn from_compressed(bytes: &[u8; 32]) -> Option<Self> {
        // We already know the point is on the curve because this is established
        // by the y-coordinate recovery procedure.

        // Obtain the three flags from the start of the byte sequence
        let compression_flag_set = Choice::from((bytes[31] >> 7) & 1);
        let sort_flag_set = Choice::from((bytes[31] >> 6) & 1);

        // Attempt to obtain the x-coordinate
        let x = {
            let mut tmp = [0; 32];
            tmp.copy_from_slice(&bytes[0..32]);

            // Mask away the flag bits
            tmp[31] &= 0b0011_1111;

            Fp4::from_bytes(&tmp)
        };

        // Return the identity assuming
        // the x-coordinate is zero and the sort bit is not set.
        //
        // Otherwise, return a recovered point (assuming the correct
        // y-coordinate can be found) so long as the infinity flag
        // was not set.
        if bool::from(
            compression_flag_set & // Compression flag should be set
                (!sort_flag_set) & // Sort flag should not be set
                x.is_zero(),
        ) {
            Some(AffinePoint::identity())
        } else {
            None
        }
        .or_else(|| {
            // Recover a y-coordinate given x by y = sqrt(x^3 + x + B)
            ((x.square() * x) + x + B).sqrt_vartime().and_then(|y| {
                // Switch to the correct y-coordinate if necessary.
                let y =
                    Fp4::conditional_select(&y, &-y, y.lexicographically_largest() ^ sort_flag_set);

                if bool::from(
                    compression_flag_set, // Compression flag should be set
                ) {
                    Some(AffinePoint {
                        x,
                        y,
                        infinity: Choice::from(0u8),
                    })
                } else {
                    None
                }
            })
        })
    }

    #[allow(unused)]
    /// Constructs an `AffinePoint` element without checking that it is a valid point.
    /// Assumes the coordinates do not represent the infinity point.
    pub fn from_raw_coordinates(elems: [Fp4; 2]) -> Self {
        AffinePoint {
            x: elems[0],
            y: elems[1],
            infinity: Choice::from(0u8),
        }
    }

    /// Returns true if this element is the identity (the point at infinity).
    #[inline]
    pub fn is_identity(&self) -> Choice {
        self.infinity
    }

    /// Returns true if this point is on the curve. This should always return
    /// true unless an "unchecked" API was used.
    pub fn is_on_curve(&self) -> Choice {
        // y^2 - x^3 - x ?= B

        (self.y.square() - self.x * (self.x.square() + Fp4::one())).ct_eq(&B) | self.infinity
    }

    /// Performs an affine scalar multiplication from `by`
    /// given as byte representation of a `Scalar` element.
    #[inline]
    pub fn multiply(&self, by: &[u8; 32]) -> AffinePoint {
        let mut acc = ProjectivePoint::identity();

        // This is a simple double-and-add implementation of point
        // multiplication, moving from most significant to least
        // significant bit of the scalar.
        //
        // We skip the first 3 leading bits because it's always unset for Fq
        // elements.
        for bit in by
            .iter()
            .rev()
            .flat_map(|byte| (0..8).rev().map(move |i| ((byte >> i) & 1u8) != 0))
            .skip(7)
        {
            acc = acc.double();
            acc = if bit { acc + self } else { acc };
        }

        acc.into()
    }
}

/// A projective point
#[derive(Copy, Clone, Debug)]
pub struct ProjectivePoint {
    pub(crate) x: Fp4,
    pub(crate) y: Fp4,
    pub(crate) z: Fp4,
}

impl Default for ProjectivePoint {
    fn default() -> ProjectivePoint {
        ProjectivePoint::identity()
    }
}

impl zeroize::DefaultIsZeroes for ProjectivePoint {}

impl fmt::Display for ProjectivePoint {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl<'a> From<&'a AffinePoint> for ProjectivePoint {
    fn from(p: &'a AffinePoint) -> ProjectivePoint {
        ProjectivePoint {
            x: p.x,
            y: p.y,
            z: Fp4::conditional_select(&Fp4::one(), &Fp4::zero(), p.infinity),
        }
    }
}

impl From<AffinePoint> for ProjectivePoint {
    fn from(p: AffinePoint) -> ProjectivePoint {
        ProjectivePoint::from(&p)
    }
}

impl ConstantTimeEq for ProjectivePoint {
    fn ct_eq(&self, other: &Self) -> Choice {
        // Is (xz, yz, z) equal to (x'z', y'z', z') when converted to affine?

        let x1 = self.x * other.z;
        let x2 = other.x * self.z;

        let y1 = self.y * other.z;
        let y2 = other.y * self.z;

        let self_is_zero = self.z.is_zero();
        let other_is_zero = other.z.is_zero();

        (self_is_zero & other_is_zero) // Both point at infinity
            | ((!self_is_zero) & (!other_is_zero) & x1.ct_eq(&x2) & y1.ct_eq(&y2))
        // Neither point at infinity, coordinates are the same
    }
}

impl ConditionallySelectable for ProjectivePoint {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        ProjectivePoint {
            x: Fp4::conditional_select(&a.x, &b.x, choice),
            y: Fp4::conditional_select(&a.y, &b.y, choice),
            z: Fp4::conditional_select(&a.z, &b.z, choice),
        }
    }
}

impl Eq for ProjectivePoint {}
impl PartialEq for ProjectivePoint {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl<'a> Neg for &'a ProjectivePoint {
    type Output = ProjectivePoint;

    #[inline]
    fn neg(self) -> ProjectivePoint {
        ProjectivePoint {
            x: self.x,
            y: -self.y,
            z: self.z,
        }
    }
}

impl Neg for ProjectivePoint {
    type Output = ProjectivePoint;

    #[inline]
    fn neg(self) -> ProjectivePoint {
        -&self
    }
}

impl<'a, 'b> Add<&'b ProjectivePoint> for &'a ProjectivePoint {
    type Output = ProjectivePoint;

    #[inline]
    fn add(self, rhs: &'b ProjectivePoint) -> ProjectivePoint {
        self.add(rhs)
    }
}

impl<'a, 'b> Sub<&'b ProjectivePoint> for &'a ProjectivePoint {
    type Output = ProjectivePoint;

    #[inline]
    fn sub(self, rhs: &'b ProjectivePoint) -> ProjectivePoint {
        self + (-rhs)
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a ProjectivePoint {
    type Output = ProjectivePoint;

    fn mul(self, other: &'b Scalar) -> Self::Output {
        self.multiply(&other.to_bytes())
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a AffinePoint {
    type Output = ProjectivePoint;

    fn mul(self, other: &'b Scalar) -> Self::Output {
        ProjectivePoint::from(self).multiply(&other.to_bytes())
    }
}

impl_binops_additive!(ProjectivePoint, ProjectivePoint);
impl_binops_multiplicative!(ProjectivePoint, Scalar);
impl_binops_multiplicative_mixed!(AffinePoint, Scalar, ProjectivePoint);

#[inline(always)]
fn mul_by_3b(a: Fp4) -> Fp4 {
    let b3 = B + B + B;
    a * b3
}

impl ProjectivePoint {
    /// Returns the identity of the group: the point at infinity.
    pub fn identity() -> ProjectivePoint {
        ProjectivePoint {
            x: Fp4::zero(),
            y: Fp4::one(),
            z: Fp4::zero(),
        }
    }

    /// Returns the x coordinate of this ProjectivePoint
    pub fn get_x(&self) -> Fp4 {
        self.x
    }

    /// Returns the y coordinate of this ProjectivePoint
    pub fn get_y(&self) -> Fp4 {
        self.y
    }

    /// Returns the z coordinate of this ProjectivePoint
    pub fn get_z(&self) -> Fp4 {
        self.z
    }

    /// Computes a random `ProjectivePoint` element
    pub fn random(mut rng: impl RngCore + CryptoRng) -> Self {
        AffinePoint::random(&mut rng).into()
    }

    /// Returns a fixed generator of the curve in projective coordinates
    pub fn generator() -> ProjectivePoint {
        ProjectivePoint::from(AffinePoint::generator())
    }

    /// Outputs a compress byte representation of this `ProjectivePoint` element
    pub fn to_compressed(&self) -> [u8; 32] {
        AffinePoint::from(self).to_compressed()
    }

    /// Outputs an uncompressed byte representation of this `ProjectivePoint` element
    /// It is twice larger than when calling `ProjectivePoint::to_uncompress()`
    pub fn to_uncompressed(&self) -> [u8; 64] {
        AffinePoint::from(self).to_uncompressed()
    }

    /// Attempts to deserialize an uncompressed element.
    pub fn from_uncompressed(bytes: &[u8; 64]) -> CtOption<Self> {
        AffinePoint::from_uncompressed(bytes)
            .and_then(|p| CtOption::new(ProjectivePoint::from(p), 1.into()))
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is on the curve and not checking if it is in the correct subgroup.
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_uncompressed()` instead.
    pub fn from_uncompressed_unchecked(bytes: &[u8; 64]) -> CtOption<Self> {
        AffinePoint::from_uncompressed_unchecked(bytes)
            .and_then(|p| CtOption::new(ProjectivePoint::from(p), 1.into()))
    }

    /// Attempts to deserialize a compressed element.
    pub fn from_compressed(bytes: &[u8; 32]) -> Option<Self> {
        AffinePoint::from_compressed(bytes).map(ProjectivePoint::from)
    }

    #[allow(unused)]
    /// Constructs a `ProjectivePoint` element without checking that it is a valid point.
    pub const fn from_raw_coordinates(elems: [Fp4; 3]) -> Self {
        ProjectivePoint {
            x: elems[0],
            y: elems[1],
            z: elems[2],
        }
    }

    /// Computes the doubling of this point.
    pub fn double(&self) -> ProjectivePoint {
        // Algorithm 3, https://eprint.iacr.org/2015/1060.pdf

        let t0 = self.x.square();
        let t1 = self.y.square();
        let t2 = self.z.square();

        let t3 = self.x * self.y;
        let t3 = t3 + t3;
        let z3 = self.x * self.z;

        let z3 = z3 + z3;
        let y3 = mul_by_3b(t2);

        let y3 = z3 + y3;
        let x3 = t1 - y3;
        let y3 = t1 + y3;

        let y3 = x3 * y3;
        let x3 = t3 * x3;
        let z3 = mul_by_3b(z3);

        let t3 = t0 - t2;

        let t3 = t3 + z3;
        let z3 = t0 + t0;
        let t0 = z3 + t0;

        let t0 = t0 + t2;
        let t0 = t0 * t3;
        let y3 = y3 + t0;

        let t2 = self.y * self.z;
        let t2 = t2 + t2;
        let t0 = t2 * t3;

        let x3 = x3 - t0;
        let z3 = t2 * t1;
        let z3 = z3 + z3;

        let z3 = z3 + z3;

        let tmp = ProjectivePoint {
            x: x3,
            y: y3,
            z: z3,
        };

        ProjectivePoint::conditional_select(&tmp, &ProjectivePoint::identity(), self.is_identity())
    }

    /// Adds this point to another point.
    pub fn add(&self, rhs: &ProjectivePoint) -> ProjectivePoint {
        // Algorithm 1, https://eprint.iacr.org/2015/1060.pdf

        let t0 = self.x * rhs.x;
        let t1 = self.y * rhs.y;
        let t2 = self.z * rhs.z;

        let t3 = self.x + self.y;
        let t4 = rhs.x + rhs.y;
        let t3 = t3 * t4;

        let t4 = t0 + t1;
        let t3 = t3 - t4;
        let t4 = self.x + self.z;

        let t5 = rhs.x + rhs.z;
        let t4 = t4 * t5;
        let t5 = t0 + t2;

        let t4 = t4 - t5;
        let t5 = self.y + self.z;
        let x3 = rhs.y + rhs.z;

        let t5 = t5 * x3;
        let x3 = t1 + t2;
        let t5 = t5 - x3;

        let x3 = mul_by_3b(t2);
        let z3 = x3 + t4;

        let x3 = t1 - z3;
        let z3 = t1 + z3;
        let y3 = x3 * z3;

        let t1 = t0 + t0;
        let t1 = t1 + t0;

        let t4 = mul_by_3b(t4);
        let t1 = t1 + t2;
        let t2 = t0 - t2;

        let t4 = t4 + t2;
        let t0 = t1 * t4;

        let y3 = y3 + t0;
        let t0 = t5 * t4;
        let x3 = t3 * x3;

        let x3 = x3 - t0;
        let t0 = t3 * t1;
        let z3 = t5 * z3;

        let z3 = z3 + t0;

        ProjectivePoint {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    /// Adds this point to another point in the affine model.
    pub fn add_mixed(&self, rhs: &AffinePoint) -> ProjectivePoint {
        // Algorithm 2, https://eprint.iacr.org/2015/1060.pdf

        let t0 = self.x * rhs.x;
        let t1 = self.y * rhs.y;
        let t3 = rhs.x + rhs.y;

        let t4 = self.x + self.y;
        let t3 = t3 * t4;
        let t4 = t0 + t1;

        let t3 = t3 - t4;
        let t4 = rhs.x * self.z;
        let t4 = t4 + self.x;

        let t5 = rhs.y * self.z;
        let t5 = t5 + self.y;

        let x3 = mul_by_3b(self.z);
        let z3 = x3 + t4;
        let x3 = t1 - z3;

        let z3 = t1 + z3;
        let y3 = x3 * z3;
        let t1 = t0 + t0;

        let t1 = t1 + t0;
        let t4 = mul_by_3b(t4);

        let t1 = t1 + self.z;
        let t2 = t0 - self.z;

        let t4 = t4 + t2;
        let t0 = t1 * t4;
        let y3 = y3 + t0;

        let t0 = t5 * t4;
        let x3 = t3 * x3;
        let x3 = x3 - t0;

        let t0 = t3 * t1;
        let z3 = t5 * z3;
        let z3 = z3 + t0;

        let tmp = ProjectivePoint {
            x: x3,
            y: y3,
            z: z3,
        };

        ProjectivePoint::conditional_select(&tmp, self, rhs.is_identity())
    }

    /// Performs a projective scalar multiplication from `by`
    /// given as byte representation of a `Scalar` element
    pub fn multiply(&self, by: &[u8; 32]) -> ProjectivePoint {
        let mut acc = ProjectivePoint::identity();

        // This is a simple double-and-add implementation of point
        // multiplication, moving from most significant to least
        // significant bit of the scalar.
        //
        // We skip the first 3 leading bits because it's always unset for Fq
        // elements.
        for bit in by
            .iter()
            .rev()
            .flat_map(|byte| (0..8).rev().map(move |i| ((byte >> i) & 1u8) != 0))
            .skip(7)
        {
            acc = acc.double();
            acc = if bit { acc + self } else { acc };
        }

        acc
    }

    /// Converts a batch of `ProjectivePoint` elements into `AffinePoint` elements. This
    /// function will panic if `p.len() != q.len()`.
    pub fn batch_normalize(p: &[Self], q: &mut [AffinePoint]) {
        assert_eq!(p.len(), q.len());

        let mut acc = Fp4::one();
        for (p, q) in p.iter().zip(q.iter_mut()) {
            // We use the `x` field of `AffinePoint` to store the product
            // of previous z-coordinates seen.
            q.x = acc;

            // We will end up skipping all identities in p
            acc = Fp4::conditional_select(&(acc * p.z), &acc, p.is_identity());
        }

        // This is the inverse, as all z-coordinates are nonzero and the ones
        // that are not are skipped.
        acc = acc.invert().unwrap();

        for (p, q) in p.iter().rev().zip(q.iter_mut().rev()) {
            let skip = p.is_identity();

            // Compute tmp = 1/z`
            let tmp = q.x * acc;

            // Cancel out z-coordinate in denominator of `acc`
            acc = Fp4::conditional_select(&(acc * p.z), &acc, skip);

            // Set the coordinates to the correct value
            q.x = p.x * tmp;
            q.y = p.y * tmp;
            q.infinity = Choice::from(0);

            *q = AffinePoint::conditional_select(q, &AffinePoint::identity(), skip);
        }
    }

    /// Returns true if this element is the identity (the point at infinity).
    #[inline]
    pub fn is_identity(&self) -> Choice {
        self.z.is_zero()
    }

    /// Returns true if this point is on the curve. This should always return
    /// true unless an "unchecked" API was used.
    pub fn is_on_curve(&self) -> Choice {
        // Y^2 Z = X^3 + X Z^2 + b Z^3

        (self.y.square() * self.z).ct_eq(
            &(self.x.square() * self.x + self.x * self.z.square() + self.z.square() * self.z * B),
        ) | (self.z.is_zero())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::thread_rng;

    #[test]
    fn test_is_on_curve() {
        assert!(bool::from(AffinePoint::identity().is_on_curve()));
        assert!(bool::from(AffinePoint::generator().is_on_curve()));
        assert!(bool::from(ProjectivePoint::identity().is_on_curve()));
        assert!(bool::from(ProjectivePoint::generator().is_on_curve()));

        let z = Fp4 {
            c0: Fp::from_raw_unchecked(602829550959103537),
            c1: Fp::from_raw_unchecked(2650934129500922588),
            c2: Fp::from_raw_unchecked(2639669314250633790),
            c3: Fp::from_raw_unchecked(4605125144701140866),
        };

        let gen = AffinePoint::generator();
        let mut test = ProjectivePoint {
            x: gen.x * z,
            y: gen.y * z,
            z,
        };

        assert!(bool::from(test.is_on_curve()));

        test.x = z;
        assert!(!bool::from(test.is_on_curve()));
    }

    #[test]
    #[allow(clippy::eq_op)]
    fn test_affine_point_equality() {
        let a = AffinePoint::generator();
        let b = AffinePoint::identity();
        let c = AffinePoint::default();

        assert!(a == a);
        assert!(b == b);
        assert!(b == c);
        assert!(a != b);
        assert!(b != a);

        assert!(bool::from(b.is_identity()));
        assert!(!bool::from(a.ct_eq(&b)));
    }

    #[test]
    #[allow(clippy::eq_op)]
    fn test_projective_point_equality() {
        let a = ProjectivePoint::generator();
        let b = ProjectivePoint::identity();
        let c = ProjectivePoint::default();

        assert!(a == a);
        assert!(b == b);
        assert!(b == c);
        assert!(a != b);
        assert!(b != a);

        assert!(bool::from(b.is_identity()));
        assert!(!bool::from(a.ct_eq(&b)));

        let z = Fp4 {
            c0: Fp::from_raw_unchecked(602829550959103537),
            c1: Fp::from_raw_unchecked(2650934129500922588),
            c2: Fp::from_raw_unchecked(2639669314250633790),
            c3: Fp::from_raw_unchecked(4605125144701140866),
        };

        let mut c = ProjectivePoint {
            x: a.x * z,
            y: a.y * z,
            z,
        };
        assert!(bool::from(c.is_on_curve()));

        assert!(a == c);
        assert!(b != c);
        assert!(c == a);
        assert!(c != b);

        c.y = -c.y;
        assert!(bool::from(c.is_on_curve()));

        assert!(a != c);
        assert!(b != c);
        assert!(c != a);
        assert!(c != b);

        c.y = -c.y;
        c.x = z;
        assert!(!bool::from(c.is_on_curve()));
        assert!(a != b);
        assert!(a != c);
        assert!(b != c);
    }

    #[test]
    fn test_projective_to_affine() {
        let a = ProjectivePoint::generator();
        let b = ProjectivePoint::identity();

        assert!(bool::from(AffinePoint::from(a).is_on_curve()));
        assert!(!bool::from(AffinePoint::from(a).is_identity()));
        assert!(bool::from(AffinePoint::from(b).is_on_curve()));
        assert!(bool::from(AffinePoint::from(b).is_identity()));

        let z = Fp4 {
            c0: Fp::from_raw_unchecked(602829550959103537),
            c1: Fp::from_raw_unchecked(2650934129500922588),
            c2: Fp::from_raw_unchecked(2639669314250633790),
            c3: Fp::from_raw_unchecked(4605125144701140866),
        };

        let c = ProjectivePoint {
            x: a.x * z,
            y: a.y * z,
            z,
        };

        assert_eq!(AffinePoint::from(c), AffinePoint::generator());
    }

    #[test]
    fn test_affine_to_projective() {
        let a = AffinePoint::generator();
        let b = AffinePoint::identity();

        assert!(bool::from(ProjectivePoint::from(a).is_on_curve()));
        assert!(!bool::from(ProjectivePoint::from(a).is_identity()));
        assert!(bool::from(ProjectivePoint::from(b).is_on_curve()));
        assert!(bool::from(ProjectivePoint::from(b).is_identity()));
    }

    #[test]
    fn test_doubling() {
        {
            let tmp = ProjectivePoint::identity().double();
            assert!(bool::from(tmp.is_identity()));
            assert!(bool::from(tmp.is_on_curve()));
        }
        {
            let tmp = ProjectivePoint::generator().double();
            assert!(!bool::from(tmp.is_identity()));
            assert!(bool::from(tmp.is_on_curve()));

            assert_eq!(
                AffinePoint::from(tmp),
                AffinePoint {
                    x: Fp4 {
                        c0: Fp::from_raw_unchecked(318189894563717415),
                        c1: Fp::from_raw_unchecked(3693971075707880361),
                        c2: Fp::from_raw_unchecked(1364015372176290066),
                        c3: Fp::from_raw_unchecked(3144734462061498524),
                    },
                    y: Fp4 {
                        c0: Fp::from_raw_unchecked(2071941157988472139),
                        c1: Fp::from_raw_unchecked(930813477875358905),
                        c2: Fp::from_raw_unchecked(1458973203145121080),
                        c3: Fp::from_raw_unchecked(2111343693983758130),
                    },
                    infinity: Choice::from(0u8),
                }
            );
        }
    }

    #[test]
    fn test_projective_addition() {
        {
            let a = ProjectivePoint::identity();
            let b = ProjectivePoint::identity();
            let c = a + b;
            assert!(bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
        }
        {
            let a = ProjectivePoint::identity();
            let mut b = ProjectivePoint::generator();
            {
                let z = Fp4 {
                    c0: Fp::from_raw_unchecked(2071941157988472139),
                    c1: Fp::from_raw_unchecked(930813477875358905),
                    c2: Fp::from_raw_unchecked(1458973203145121080),
                    c3: Fp::from_raw_unchecked(2111343693983758130),
                };

                b = ProjectivePoint {
                    x: b.x * z,
                    y: b.y * z,
                    z,
                };
            }
            let c = a + b;
            assert!(!bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
            assert!(c == ProjectivePoint::generator());
        }
        {
            let a = ProjectivePoint::generator().double().double(); // 4P
            let b = ProjectivePoint::generator().double(); // 2P
            let c = a + b;

            let mut d = ProjectivePoint::generator();
            for _ in 0..5 {
                d += ProjectivePoint::generator();
            }
            assert!(!bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
            assert!(!bool::from(d.is_identity()));
            assert!(bool::from(d.is_on_curve()));
            assert_eq!(c, d);
        }
    }

    #[test]
    fn test_mixed_addition() {
        {
            let a = AffinePoint::identity();
            let b = ProjectivePoint::identity();
            let c = a + b;
            assert!(bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
        }
        {
            let a = AffinePoint::identity();
            let mut b = ProjectivePoint::generator();
            {
                let z = Fp4 {
                    c0: Fp::from_raw_unchecked(2071941157988472139),
                    c1: Fp::from_raw_unchecked(930813477875358905),
                    c2: Fp::from_raw_unchecked(1458973203145121080),
                    c3: Fp::from_raw_unchecked(2111343693983758130),
                };

                b = ProjectivePoint {
                    x: b.x * z,
                    y: b.y * z,
                    z,
                };
            }
            let c = a + b;
            assert!(!bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
            assert!(c == ProjectivePoint::generator());
        }
        {
            let a = AffinePoint::identity();
            let mut b = ProjectivePoint::generator();
            {
                let z = Fp4 {
                    c0: Fp::from_raw_unchecked(2071941157988472139),
                    c1: Fp::from_raw_unchecked(930813477875358905),
                    c2: Fp::from_raw_unchecked(1458973203145121080),
                    c3: Fp::from_raw_unchecked(2111343693983758130),
                };

                b = ProjectivePoint {
                    x: b.x * z,
                    y: b.y * z,
                    z,
                };
            }
            let c = b + a;
            assert!(!bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
            assert!(c == ProjectivePoint::generator());
        }
        {
            let a = ProjectivePoint::generator().double().double(); // 4P
            let b = ProjectivePoint::generator().double(); // 2P
            let c = a + b;

            let mut d = ProjectivePoint::generator();
            for _ in 0..5 {
                d += AffinePoint::generator();
            }
            assert!(!bool::from(c.is_identity()));
            assert!(bool::from(c.is_on_curve()));
            assert!(!bool::from(d.is_identity()));
            assert!(bool::from(d.is_on_curve()));
            assert_eq!(c, d);
        }
    }

    #[test]
    #[allow(clippy::eq_op)]
    fn test_projective_negation_and_subtraction() {
        let a = ProjectivePoint::generator().double();
        assert_eq!(a + (-a), ProjectivePoint::identity());
        assert_eq!(a + (-a), a - a);
    }

    #[test]
    fn test_affine_negation_and_subtraction() {
        let a = AffinePoint::generator();
        assert_eq!(ProjectivePoint::from(a) + (-a), ProjectivePoint::identity());
        assert_eq!(
            ProjectivePoint::from(a) + (-a),
            ProjectivePoint::from(a) - a
        );
    }

    #[test]
    fn test_projective_scalar_multiplication() {
        let g = ProjectivePoint::generator();
        let a = Scalar::new([
            0xef427d940c471145,
            0xf9d1c30637e9f84d,
            0x843a5b754596e86b,
            0x05b910f89b6b601c,
        ]);
        let b = Scalar::new([
            0xcdf47d5adc756906,
            0x381699324f082566,
            0x725be442943c3f0f,
            0x0701db10daaec421,
        ]);
        let c = a * b;

        assert_eq!((g * a) * b, g * c);

        for _ in 0..100 {
            let a = Scalar::random(&mut thread_rng());
            let b = Scalar::random(&mut thread_rng());
            let c = a * b;

            assert_eq!((g * a) * b, g * c);
        }
    }

    #[test]
    fn test_affine_scalar_multiplication() {
        let g = AffinePoint::generator();
        let a = Scalar::new([
            0xb951ca4b11baeb8c,
            0xbd8bccd724d2d460,
            0x3520dbe0f992ab40,
            0x02a7506357d39b4e,
        ]);
        let b = Scalar::new([
            0x80996fb6c25f0316,
            0xa518a33400a43fdd,
            0x8e456b2de42d5671,
            0x0401b958b504dd68,
        ]);
        let c = a * b;

        assert_eq!(AffinePoint::from(g * a) * b, g * c);

        for _ in 0..100 {
            let a = Scalar::random(&mut thread_rng());
            let b = Scalar::random(&mut thread_rng());
            let c = a * b;

            assert_eq!((g * a) * b, g * c);
        }
    }

    #[test]
    fn test_batch_normalize() {
        let a = ProjectivePoint::generator().double();
        let b = a.double();
        let c = b.double();

        for a_identity in (0..1).map(|n| n == 1) {
            for b_identity in (0..1).map(|n| n == 1) {
                for c_identity in (0..1).map(|n| n == 1) {
                    let mut v = [a, b, c];
                    if a_identity {
                        v[0] = ProjectivePoint::identity()
                    }
                    if b_identity {
                        v[1] = ProjectivePoint::identity()
                    }
                    if c_identity {
                        v[2] = ProjectivePoint::identity()
                    }

                    let mut t = [
                        AffinePoint::identity(),
                        AffinePoint::identity(),
                        AffinePoint::identity(),
                    ];
                    let expected = [
                        AffinePoint::from(v[0]),
                        AffinePoint::from(v[1]),
                        AffinePoint::from(v[2]),
                    ];

                    ProjectivePoint::batch_normalize(&v[..], &mut t[..]);

                    assert_eq!(&t[..], &expected[..]);
                }
            }
        }
    }

    #[test]
    fn test_zeroize() {
        use zeroize::Zeroize;

        let mut a = AffinePoint::generator();
        a.zeroize();
        assert_eq!(a, AffinePoint::identity());

        let mut a = ProjectivePoint::generator();
        a.zeroize();
        assert_eq!(a, ProjectivePoint::identity());
    }

    // POINT COMPRESSION
    // ================================================================================================

    #[test]
    fn test_point_compressed() {
        let mut rng = thread_rng();
        // Random points
        for _ in 0..100 {
            let point = AffinePoint::random(&mut rng);
            let bytes = point.to_compressed();
            let point_decompressed = AffinePoint::from_compressed(&bytes).unwrap();
            assert_eq!(point, point_decompressed);

            let point = ProjectivePoint::random(&mut rng);
            let bytes = point.to_compressed();
            let point_decompressed = ProjectivePoint::from_compressed(&bytes).unwrap();
            assert_eq!(point, point_decompressed);
        }

        // Identity point
        {
            let bytes = AffinePoint::identity().to_compressed();
            let point_decompressed = AffinePoint::from_compressed(&bytes).unwrap();
            assert!(bool::from(point_decompressed.is_identity()));

            let bytes = ProjectivePoint::identity().to_compressed();
            let point_decompressed = ProjectivePoint::from_compressed(&bytes).unwrap();
            assert!(bool::from(point_decompressed.is_identity()));
        }

        // Invalid points
        {
            let point = AffinePoint {
                x: Fp4::zero(),
                y: Fp4::zero(),
                infinity: 0.into(),
            };
            let bytes = point.to_compressed();
            let point_decompressed = AffinePoint::from_compressed(&bytes);
            assert_eq!(point_decompressed.unwrap(), AffinePoint::identity());

            let point = ProjectivePoint::from(&point);
            let bytes = point.to_compressed();
            let point_decompressed = ProjectivePoint::from_compressed(&bytes);
            assert_eq!(point_decompressed.unwrap(), ProjectivePoint::identity());
        }
    }

    #[test]
    fn test_point_uncompressed() {
        let mut rng = thread_rng();

        // Random points
        for _ in 0..100 {
            let point = AffinePoint::random(&mut rng);
            let bytes = point.to_uncompressed();
            let point_decompressed = AffinePoint::from_uncompressed(&bytes).unwrap();
            assert_eq!(point, point_decompressed);

            let point = ProjectivePoint::random(&mut rng);
            let bytes = point.to_uncompressed();
            let point_decompressed = ProjectivePoint::from_uncompressed(&bytes).unwrap();
            assert_eq!(point, point_decompressed);
        }

        // Identity point
        {
            let bytes = AffinePoint::identity().to_uncompressed();
            let point_decompressed = AffinePoint::from_uncompressed(&bytes).unwrap();
            assert!(bool::from(point_decompressed.is_identity()));

            let bytes = ProjectivePoint::identity().to_uncompressed();
            let point_decompressed = ProjectivePoint::from_uncompressed(&bytes).unwrap();
            assert!(bool::from(point_decompressed.is_identity()));
        }

        // Invalid points
        {
            let point = AffinePoint {
                x: Fp4::zero(),
                y: Fp4::zero(),
                infinity: 0.into(),
            };
            let bytes = point.to_uncompressed();
            let point_decompressed = AffinePoint::from_uncompressed(&bytes);
            assert_eq!(point_decompressed.unwrap(), AffinePoint::identity());

            let point = ProjectivePoint::from(&point);
            let bytes = point.to_uncompressed();
            let point_decompressed = ProjectivePoint::from_uncompressed(&bytes);
            assert_eq!(point_decompressed.unwrap(), ProjectivePoint::identity());
        }
    }
}
