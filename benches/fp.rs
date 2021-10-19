use rand::thread_rng;

#[macro_use]
extern crate criterion;

use criterion::black_box;
use criterion::Criterion;

extern crate cheetah;

use cheetah::Fp;
use group::ff::Field;

fn criterion_benchmark(c: &mut Criterion) {
    let mut rng = thread_rng();

    c.bench_function("fp add", |bench| {
        let x = Fp::random(&mut rng);
        let y = Fp::random(&mut rng);
        bench.iter(|| black_box(x) + black_box(y))
    });

    c.bench_function("fp sub", |bench| {
        let x = Fp::random(&mut rng);
        let y = Fp::random(&mut rng);
        bench.iter(|| black_box(x) - black_box(y))
    });

    c.bench_function("fp mul", |bench| {
        let x = Fp::random(&mut rng);
        let y = Fp::random(&mut rng);
        bench.iter(|| black_box(x) * black_box(y))
    });

    c.bench_function("fp square", |bench| {
        let x = Fp::random(&mut rng);
        bench.iter(|| black_box(x).square())
    });

    c.bench_function("fp square_from_mul", |bench| {
        let x = Fp::random(&mut rng);
        bench.iter(|| black_box(x) * black_box(x))
    });

    c.bench_function("fp exp", |bench| {
        let x = Fp::random(&mut rng);
        let y = Fp::random(&mut rng).to_repr();
        bench.iter(|| Fp::exp(black_box(x), black_box(y)))
    });

    c.bench_function("fp invert", |bench| {
        let x = Fp::random(&mut rng);
        bench.iter(|| Fp::invert(black_box(&x)))
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default();
    targets = criterion_benchmark);
criterion_main!(benches);
