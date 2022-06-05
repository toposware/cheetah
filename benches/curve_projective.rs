// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use rand_core::OsRng;

#[macro_use]
extern crate criterion;

use criterion::black_box;
use criterion::Criterion;

extern crate cheetah;

use cheetah::Scalar;
use cheetah::{AffinePoint, ProjectivePoint};

use cheetah::{BasePointTable, BASEPOINT_TABLE};

static BATCH_SIZES: [u32; 5] = [1, 10, 100, 1000, 10000];

fn criterion_benchmark(c: &mut Criterion) {
    let mut rng = OsRng;
    let p = ProjectivePoint::random(&mut rng);
    let p2 = ProjectivePoint::random(&mut rng);
    let p2_affine = AffinePoint::from(&p2);
    let pow = Scalar::new([
        0xb283b88ede33b702,
        0x4460d73106ed6dad,
        0x6c2f3c993870f926,
        0x71acba2667ccf016,
    ])
    .to_bytes();
    let pow2 = Scalar::new([
        0xafd872bc45b8fbd0,
        0xe714cd570760ccad,
        0x4ec25207c1649a1a,
        0x6fa45dca4626832f,
    ])
    .to_bytes();

    let uncompressed_encoding = p.to_uncompressed();
    let compressed_encoding = p.to_compressed();

    c.bench_function("Projective addition", |bench| {
        bench.iter(|| black_box(p) + black_box(p2))
    });

    c.bench_function("Projective mixed addition", |bench| {
        bench.iter(|| black_box(p) + black_box(p2_affine))
    });

    c.bench_function("Projective unchecked mixed addition", |bench| {
        bench.iter(|| black_box(p).add_mixed_unchecked(black_box(&p2_affine)))
    });

    c.bench_function("Projective subtraction", |bench| {
        bench.iter(|| black_box(p) - black_box(p2))
    });

    c.bench_function("Projective doubling", |bench| {
        bench.iter(|| black_box(p).double())
    });

    c.bench_function("Projective unchecked doubling", |bench| {
        bench.iter(|| black_box(p).double_unchecked())
    });

    c.bench_function(
        "Projective scalar multiplication (variable base)",
        |bench| bench.iter(|| ProjectivePoint::multiply(black_box(&p), black_box(&pow))),
    );

    c.bench_function(
        "Projective scalar multiplication (variable base) - variable time",
        |bench| bench.iter(|| ProjectivePoint::multiply_vartime(black_box(&p), black_box(&pow))),
    );

    c.bench_function(
        "Projective double scalar multiplication (variable base)",
        |bench| {
            bench.iter(|| {
                ProjectivePoint::multiply_double(
                    black_box(&p),
                    black_box(&p2),
                    black_box(&pow),
                    black_box(&pow2),
                )
            })
        },
    );

    c.bench_function(
        "Projective double scalar multiplication (variable base) - variable time",
        |bench| {
            bench.iter(|| {
                ProjectivePoint::multiply_double_vartime(
                    black_box(&p),
                    black_box(&p2),
                    black_box(&pow),
                    black_box(&pow2),
                )
            })
        },
    );

    c.bench_function(
        "Projective double scalar multiplication with basepoint - variable time",
        |bench| {
            bench.iter(|| {
                ProjectivePoint::multiply_double_with_basepoint_vartime(
                    black_box(&p),
                    black_box(&pow),
                    black_box(&pow2),
                )
            })
        },
    );

    c.bench_function("Projective scalar multiplication (basepoint)", |bench| {
        bench.iter(|| BasePointTable::multiply(black_box(&BASEPOINT_TABLE), black_box(&pow)))
    });

    c.bench_function(
        "Projective scalar multiplication (basepoint) - variable time",
        |bench| {
            bench.iter(|| {
                BasePointTable::multiply_vartime(black_box(&BASEPOINT_TABLE), black_box(&pow))
            })
        },
    );

    c.bench_function("Projective basepoint table creation", |bench| {
        bench.iter(|| BasePointTable::from(black_box(&p)))
    });

    c.bench_function("Projective uncompressed encoding", |bench| {
        bench.iter(|| ProjectivePoint::to_uncompressed(black_box(&p)))
    });

    c.bench_function("Projective uncompressed decoding", |bench| {
        bench.iter(|| ProjectivePoint::from_uncompressed(black_box(&uncompressed_encoding)))
    });

    c.bench_function("Projective compressed encoding", |bench| {
        bench.iter(|| ProjectivePoint::to_compressed(black_box(&p)))
    });

    c.bench_function("Projective compressed decoding", |bench| {
        bench.iter(|| ProjectivePoint::from_compressed(black_box(&compressed_encoding)))
    });

    c.bench_function("Projective curve check", |bench| {
        bench.iter(|| ProjectivePoint::is_on_curve(black_box(&p)))
    });

    c.bench_function("Projective subgroup check", |bench| {
        bench.iter(|| ProjectivePoint::is_torsion_free(black_box(&p)))
    });

    c.bench_function("Projective points equality", |bench| {
        bench.iter(|| black_box(p) == black_box(p))
    });

    c.bench_function("Projective conversion to affine", |bench| {
        bench.iter(|| AffinePoint::from(black_box(&p)))
    });

    let batch_str = "Projective batch normalize to affine - ".to_string();
    for &batch_size in BATCH_SIZES.iter() {
        let name = batch_str.clone() + &batch_size.to_string();
        c.bench_function(&name, |bench| {
            let projective_points = vec![ProjectivePoint::random(&mut rng); batch_size as usize];
            let mut affine_points = vec![AffinePoint::identity(); batch_size as usize];
            bench.iter(|| {
                ProjectivePoint::batch_normalize(
                    black_box(&projective_points),
                    black_box(&mut affine_points),
                )
            })
        });
    }
}

criterion_group!(
    name = benches;
    config = Criterion::default();
    targets = criterion_benchmark);
criterion_main!(benches);
