// Copyright (c) 2021-2023 Toposware, Inc.
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

fn criterion_benchmark(c: &mut Criterion) {
    let mut rng = OsRng;
    let x = Scalar::random(&mut rng);
    let x_squared = x.square();
    let x_bytes = x.to_bytes();
    let y = Scalar::random(&mut rng);
    let pow = y.output_reduced_limbs();

    c.bench_function("Scalar add", |bench| {
        bench.iter(|| black_box(x) + black_box(y))
    });

    c.bench_function("Scalar sub", |bench| {
        bench.iter(|| black_box(x) - black_box(y))
    });

    c.bench_function("Scalar double", |bench| {
        bench.iter(|| black_box(x).double())
    });

    c.bench_function("Scalar mul", |bench| {
        bench.iter(|| black_box(x) * black_box(y))
    });

    c.bench_function("Scalar square", |bench| {
        bench.iter(|| black_box(x).square())
    });

    c.bench_function("Scalar square_from_mul", |bench| {
        bench.iter(|| black_box(x) * black_box(x))
    });

    c.bench_function("Scalar sqrt", |bench| {
        bench.iter(|| Scalar::sqrt(black_box(&x_squared)))
    });

    c.bench_function("Scalar exp", |bench| {
        bench.iter(|| Scalar::exp(black_box(&x), black_box(&pow)))
    });

    c.bench_function("Scalar invert", |bench| {
        bench.iter(|| Scalar::invert(black_box(&x)))
    });

    c.bench_function("Scalar encoding", |bench| {
        bench.iter(|| Scalar::to_bytes(black_box(&x)))
    });

    c.bench_function("Scalar decoding", |bench| {
        bench.iter(|| Scalar::from_bytes(black_box(&x_bytes)))
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default();
    targets = criterion_benchmark);
criterion_main!(benches);
