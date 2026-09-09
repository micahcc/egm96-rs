use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};

fn criterion_benchmark(c: &mut Criterion) {
    let test_points = [
        (33.4818, -117.556),  // US West Coast
        (0.0, 0.0),           // Equator / Prime Meridian
        (89.0, 0.0),          // Near Pole
        (-33.8688, 151.2093), // Southern Hemisphere
    ];

    c.bench_function("cts autoselect", |b| {
        b.iter(|| {
            egm96::egm96_altitude_offset(black_box(33.4818), black_box(-117.556));
        })
    });

    c.bench_function("cts harmonics", |b| {
        b.iter(|| {
            egm96::egm96_compute_altitude_offset(black_box(33.4818), black_box(-117.556));
        })
    });

    c.bench_function("cts harmonics multi-point", |b| {
        b.iter(|| {
            for &(lat, lon) in &test_points {
                egm96::egm96_compute_altitude_offset(black_box(lat), black_box(lon));
            }
        })
    });

    #[cfg(feature = "raster_15_min")]
    c.bench_function("cts 15 min", |b| {
        b.iter(|| {
            egm96::egm96_raster_15_min_altitude_offset(black_box(33.4818), black_box(-117.556));
        })
    });

    #[cfg(feature = "raster_5_min")]
    c.bench_function("cts 5 min", |b| {
        b.iter(|| {
            egm96::egm96_raster_5_min_altitude_offset(black_box(33.4818), black_box(-117.556));
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
