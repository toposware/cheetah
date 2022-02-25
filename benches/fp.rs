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

use cheetah::Fp;

fn criterion_benchmark(c: &mut Criterion) {
    let mut rng = OsRng;
    let x = Fp::random(&mut rng);
    let x_squared = x.square();
    let x_bytes = x.to_bytes();
    let y = Fp::random(&mut rng);
    let pow = y.output_internal();

    c.bench_function("Fp add", |bench| bench.iter(|| black_box(x) + black_box(y)));

    c.bench_function("Fp sub", |bench| bench.iter(|| black_box(x) - black_box(y)));

    c.bench_function("Fp double", |bench| bench.iter(|| black_box(x).double()));

    c.bench_function("Fp mul", |bench| bench.iter(|| black_box(x) * black_box(y)));

    c.bench_function("Fp square", |bench| bench.iter(|| black_box(x).square()));

    c.bench_function("Fp square_from_mul", |bench| {
        bench.iter(|| black_box(x) * black_box(x))
    });

    c.bench_function("Fp sqrt", |bench| {
        bench.iter(|| Fp::sqrt(black_box(&x_squared)))
    });

    c.bench_function("Fp exp", |bench| {
        bench.iter(|| Fp::exp(black_box(x), black_box(pow)))
    });

    c.bench_function("Fp invert", |bench| {
        bench.iter(|| Fp::invert(black_box(&x)))
    });

    c.bench_function("Fp encoding", |bench| {
        bench.iter(|| Fp::to_bytes(black_box(&x)))
    });

    c.bench_function("Fp decoding", |bench| {
        bench.iter(|| Fp::from_bytes(black_box(&x_bytes)))
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default();
    targets = criterion_benchmark);
criterion_main!(benches);
