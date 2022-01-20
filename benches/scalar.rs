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
use group::ff::Field;

fn criterion_benchmark(c: &mut Criterion) {
    let mut rng = OsRng;

    c.bench_function("scalar add", |bench| {
        let x = Scalar::random(&mut rng);
        let y = Scalar::random(&mut rng);
        bench.iter(|| black_box(x) + black_box(y))
    });

    c.bench_function("scalar sub", |bench| {
        let x = Scalar::random(&mut rng);
        let y = Scalar::random(&mut rng);
        bench.iter(|| black_box(x) - black_box(y))
    });

    c.bench_function("scalar double", |bench| {
        let x = Scalar::random(&mut rng);
        bench.iter(|| black_box(x).double())
    });

    c.bench_function("scalar mul", |bench| {
        let x = Scalar::random(&mut rng);
        let y = Scalar::random(&mut rng);
        bench.iter(|| black_box(x) * black_box(y))
    });

    c.bench_function("scalar square", |bench| {
        let x = Scalar::random(&mut rng);
        bench.iter(|| black_box(x).square())
    });

    c.bench_function("scalar square_from_mul", |bench| {
        let x = Scalar::random(&mut rng);
        bench.iter(|| black_box(x) * black_box(x))
    });

    c.bench_function("scalar exp", |bench| {
        let x = Scalar::random(&mut rng);
        let y = Scalar::random(&mut rng).output_reduced_limbs();
        bench.iter(|| Scalar::exp(black_box(x), black_box(&y)))
    });

    c.bench_function("scalar invert", |bench| {
        let x = Scalar::random(&mut rng);
        bench.iter(|| Scalar::invert(black_box(&x)))
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default();
    targets = criterion_benchmark);
criterion_main!(benches);
