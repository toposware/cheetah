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

use cheetah::Fp2;
use group::ff::Field;

fn criterion_benchmark(c: &mut Criterion) {
    let mut rng = OsRng;

    c.bench_function("fp2 add", |bench| {
        let x = Fp2::random(&mut rng);
        let y = Fp2::random(&mut rng);
        bench.iter(|| black_box(x) + black_box(y))
    });

    c.bench_function("fp2 sub", |bench| {
        let x = Fp2::random(&mut rng);
        let y = Fp2::random(&mut rng);
        bench.iter(|| black_box(x) - black_box(y))
    });

    c.bench_function("fp2 double", |bench| {
        let x = Fp2::random(&mut rng);
        bench.iter(|| black_box(x).double())
    });

    c.bench_function("fp2 mul", |bench| {
        let x = Fp2::random(&mut rng);
        let y = Fp2::random(&mut rng);
        bench.iter(|| black_box(x) * black_box(y))
    });

    c.bench_function("fp2 square", |bench| {
        let x = Fp2::random(&mut rng);
        bench.iter(|| black_box(x).square())
    });

    c.bench_function("fp2 square_from_mul", |bench| {
        let x = Fp2::random(&mut rng);
        bench.iter(|| black_box(x) * black_box(x))
    });

    c.bench_function("fp2 exp", |bench| {
        let x = Fp2::random(&mut rng);
        let y = Fp2::random(&mut rng).output_reduced_limbs();
        bench.iter(|| Fp2::exp(black_box(x), black_box(&y)))
    });

    c.bench_function("fp2 invert", |bench| {
        let x = Fp2::random(&mut rng);
        bench.iter(|| Fp2::invert(black_box(&x)))
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default();
    targets = criterion_benchmark);
criterion_main!(benches);
