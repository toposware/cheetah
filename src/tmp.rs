pub fn multiply(&self, by: &[u8; 32]) -> AffinePoint {
    let mut acc = ProjectivePoint::identity();

    let mut table = [*self; 17];
    for i in 2..17 {
        table[i] = (ProjectivePoint::from(table[i - 1]) + self).into();
    }

    // This is a simple double-and-add implementation of point
    // multiplication, moving from most significant to least
    // significant bit of the scalar.
    //
    // We skip the first two leading bits because they are always unset for Fq
    // elements.
    for bit in by
        .iter()
        .rev()
        .flat_map(|byte| (0..2).rev().map(move |i| (byte >> (i * 4)) & 15u8))
    {
        acc = acc.double();
        acc = acc.double();
        acc = acc.double();
        acc = acc.double();
        acc = ProjectivePoint::conditional_select(
            &acc,
            &(acc + table[bit as usize]),
            ((bit != 0) as u8).into(),
        );
    }

    acc.into()
}

/// Adds this point to another point different from +/- self.
pub fn add_different(&self, rhs: &ProjectivePoint) -> ProjectivePoint {
    // Algorithm 1, https://eprint.iacr.org/2015/1060.pdf
    let t2 = self.z * rhs.z;

    let s1 = self.y * rhs.z;
    let s2 = self.z * rhs.y;

    let u1 = self.x * rhs.z;
    let u2 = self.z * rhs.x;

    let a = s2 - s1;
    let b = u2 - u1;

    let b2 = b.square();
    let b3 = b2 * b;

    let c = a.square();
    let c = c * t2;
    let c = c - b3;

    let t3 = b2 * u1;

    let t = t3.double();
    let c = c - t;

    let x3 = b * c;

    let y3 = t3 - c;
    let y3 = a * y3;
    let t = b3 * s1;
    let y3 = y3 - t;

    let z3 = b3 * t2;

    ProjectivePoint {
        x: x3,
        y: y3,
        z: z3,
    }
}

/// Adds this point to another point different from +/- self.
pub fn add_mixed_different(&self, rhs: &AffinePoint) -> ProjectivePoint {
    // Algorithm 1, https://eprint.iacr.org/2015/1060.pdf

    let s2 = self.z * rhs.y;
    let u2 = self.z * rhs.x;

    let a = s2 - self.y;
    let b = u2 - self.x;

    let b2 = b.square();
    let b3 = b2 * b;

    let c = a.square();
    let c = c * self.z;
    let c = c - b3;

    let t3 = b2 * self.x;

    let t = t3.double();
    let c = c - t;

    let x3 = b * c;

    let y3 = t3 - c;
    let y3 = a * y3;
    let t = b3 * self.y;
    let y3 = y3 - t;

    let z3 = b3 * self.z;

    let tmp = ProjectivePoint {
        x: x3,
        y: y3,
        z: z3,
    };

    ProjectivePoint::conditional_select(&tmp, self, rhs.is_identity())
}

/// TODO
pub fn multiply_lookup_8(&self, by: &[u8; 32]) -> ProjectivePoint {
    let mut acc = ProjectivePoint::identity();
    let table = LookupTable::<8>::from(self);

    // This is a simple double-and-add implementation of point
    // multiplication, moving from most significant to least
    // significant bit of the scalar.
    //
    // We skip the first two leading bits because they are always unset for Fq
    // elements.
    for digit in by
        .iter()
        .rev()
        .flat_map(|byte| (0..4).rev().map(move |i| (byte >> (i * 2)) & 7u8))
    {
        acc = acc.double();
        acc = acc.double();
        acc = ProjectivePoint::conditional_select(
            &acc,
            &(acc + table.get_from_index(digit as usize)),
            ((digit != 0) as u8).into(),
        );
    }

    acc.into()
}

/// TODO
pub fn multiply_lookup_8_vartime(&self, by: &[u8; 32]) -> ProjectivePoint {
    let mut acc = ProjectivePoint::identity();
    let table = LookupTable::<8>::from(self);

    // This is a simple double-and-add implementation of point
    // multiplication, moving from most significant to least
    // significant bit of the scalar.
    //
    // We skip the first two leading bits because they are always unset for Fq
    // elements.
    for digit in by
        .iter()
        .rev()
        .flat_map(|byte| (0..4).rev().map(move |i| (byte >> (i * 2)) & 7u8))
    {
        acc = acc.double();
        acc = acc.double();
        if digit != 0 {
            acc += table.get_from_index(digit as usize);
        }
    }

    acc.into()
}

/// TODO
pub fn multiply_lookup_16(&self, by: &[u8; 32]) -> ProjectivePoint {
    let mut acc = ProjectivePoint::identity();
    let table = LookupTable::<16>::from(self);

    // This is a simple double-and-add implementation of point
    // multiplication, moving from most significant to least
    // significant bit of the scalar.
    //
    // We skip the first two leading bits because they are always unset for Fq
    // elements.
    for digit in by
        .iter()
        .rev()
        .flat_map(|byte| (0..2).rev().map(move |i| (byte >> (i * 4)) & 15u8))
    {
        acc = acc.double();
        acc = acc.double();
        acc = acc.double();
        acc = acc.double();
        acc = ProjectivePoint::conditional_select(
            &acc,
            &(acc + table.get_from_index(digit as usize)),
            ((digit != 0) as u8).into(),
        );
    }

    acc.into()
}

/// TODO
pub fn multiply_lookup_16_vartime(&self, by: &[u8; 32]) -> ProjectivePoint {
    let mut acc = ProjectivePoint::identity();
    let table = LookupTable::<16>::from(self);

    // This is a simple double-and-add implementation of point
    // multiplication, moving from most significant to least
    // significant bit of the scalar.
    //
    // We skip the first two leading bits because they are always unset for Fq
    // elements.
    for digit in by
        .iter()
        .rev()
        .flat_map(|byte| (0..2).rev().map(move |i| (byte >> (i * 4)) & 15u8))
    {
        acc = acc.double();
        acc = acc.double();
        acc = acc.double();
        acc = acc.double();
        if digit != 0 {
            acc += table.get_from_index(digit as usize);
        }
    }
    acc.into()
}
