use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion, Throughput};

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

    let mut harmonics_group = c.benchmark_group("harmonics multi-point");
    harmonics_group.throughput(Throughput::Elements(test_points.len() as u64));
    harmonics_group.bench_function("cts harmonics multi-point", |b| {
        b.iter(|| {
            for &(lat, lon) in &test_points {
                egm96::egm96_compute_altitude_offset(black_box(lat), black_box(lon));
            }
        })
    });
    harmonics_group.finish();

    #[cfg(feature = "raster_15_min")]
    {
        let mut r15_group = c.benchmark_group("raster 15 min multi-point");
        r15_group.throughput(Throughput::Elements(test_points.len() as u64));
        r15_group.bench_function("cts 15 min multi-point", |b| {
            b.iter(|| {
                for &(lat, lon) in &test_points {
                    egm96::egm96_raster_15_min_altitude_offset(black_box(lat), black_box(lon));
                }
            })
        });
        r15_group.finish();

        c.bench_function("cts 15 min", |b| {
            b.iter(|| {
                egm96::egm96_raster_15_min_altitude_offset(black_box(33.4818), black_box(-117.556));
            })
        });
    }

    #[cfg(feature = "raster_5_min")]
    {
        let mut r5_group = c.benchmark_group("raster 5 min multi-point");
        r5_group.throughput(Throughput::Elements(test_points.len() as u64));
        r5_group.bench_function("cts 5 min multi-point", |b| {
            b.iter(|| {
                for &(lat, lon) in &test_points {
                    egm96::egm96_raster_5_min_altitude_offset(black_box(lat), black_box(lon));
                }
            })
        });
        r5_group.finish();

        c.bench_function("cts 5 min", |b| {
            b.iter(|| {
                egm96::egm96_raster_5_min_altitude_offset(black_box(33.4818), black_box(-117.556));
            })
        });
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
