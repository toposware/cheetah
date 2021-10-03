use rand::thread_rng;

#[macro_use]
extern crate criterion;

use criterion::black_box;
use criterion::Criterion;

extern crate cheetah;

use cheetah::Fp4;

fn criterion_benchmark(c: &mut Criterion) {
    let mut rng = thread_rng();

    c.bench_function("fp4 add", |bench| {
        let x = Fp4::random(&mut rng);
        let y = Fp4::random(&mut rng);
        bench.iter(|| black_box(x) + black_box(y))
    });

    c.bench_function("fp4 sub", |bench| {
        let x = Fp4::random(&mut rng);
        let y = Fp4::random(&mut rng);
        bench.iter(|| black_box(x) - black_box(y))
    });

    c.bench_function("fp4 mul", |bench| {
        let x = Fp4::random(&mut rng);
        let y = Fp4::random(&mut rng);
        bench.iter(|| black_box(x) * black_box(y))
    });

    c.bench_function("fp4 square", |bench| {
        let x = Fp4::random(&mut rng);
        bench.iter(|| black_box(x).square())
    });

    c.bench_function("fp4 square_from_mul", |bench| {
        let x = Fp4::random(&mut rng);
        bench.iter(|| black_box(x) * black_box(x))
    });

    c.bench_function("fp4 exp", |bench| {
        let x = Fp4::random(&mut rng);
        let y = Fp4::random(&mut rng).to_repr();
        bench.iter(|| Fp4::exp(black_box(x), black_box(&y)))
    });

    c.bench_function("fp4 invert", |bench| {
        let x = Fp4::random(&mut rng);
        bench.iter(|| Fp4::invert(black_box(&x)))
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default();
    targets = criterion_benchmark);
criterion_main!(benches);
