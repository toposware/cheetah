// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use rand_core::OsRng;
use rand_core::RngCore;

#[macro_use]
extern crate criterion;

use criterion::black_box;
use criterion::Criterion;

extern crate cheetah;

use cheetah::Fp6;

fn criterion_benchmark(c: &mut Criterion) {
    let mut rng = OsRng;
    let x = Fp6::random(&mut rng);
    let x_squared = x.square();
    let x_bytes = x.to_bytes();
    let y = Fp6::random(&mut rng);
    let pow = y.output_internal();
    let y32 = rng.next_u32();

    c.bench_function("Fp6 add", |bench| {
        bench.iter(|| black_box(x) + black_box(y))
    });

    c.bench_function("Fp6 sub", |bench| {
        bench.iter(|| black_box(x) - black_box(y))
    });

    c.bench_function("Fp6 double", |bench| bench.iter(|| black_box(x).double()));

    c.bench_function("Fp6 mul u32", |bench| {
        bench.iter(|| black_box(x).mul_by_u32(black_box(y32)))
    });

    c.bench_function("Fp6 mul", |bench| {
        bench.iter(|| black_box(x) * black_box(y))
    });

    c.bench_function("Fp6 square", |bench| bench.iter(|| black_box(x).square()));

    c.bench_function("Fp6 square_from_mul", |bench| {
        bench.iter(|| black_box(x) * black_box(x))
    });

    c.bench_function("Fp6 sqrt", |bench| {
        bench.iter(|| Fp6::sqrt(black_box(&x_squared)))
    });

    c.bench_function("Fp6 exp", |bench| {
        bench.iter(|| Fp6::exp(black_box(&x), black_box(&pow)))
    });

    c.bench_function("Fp6 invert", |bench| {
        bench.iter(|| Fp6::invert(black_box(&x)))
    });

    c.bench_function("Fp6 encoding", |bench| {
        bench.iter(|| Fp6::to_bytes(black_box(&x)))
    });

    c.bench_function("Fp6 decoding", |bench| {
        bench.iter(|| Fp6::from_bytes(black_box(&x_bytes)))
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default();
    targets = criterion_benchmark);
criterion_main!(benches);
