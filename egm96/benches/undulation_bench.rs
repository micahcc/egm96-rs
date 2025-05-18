use std::hint::black_box;

use egm96::egm96_compute_altitude_offset;

use criterion::{criterion_group, criterion_main, Criterion};

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("cts undulation", |b| {
        b.iter(|| {
            egm96_compute_altitude_offset(black_box(33.4818), black_box(-117.556));
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
