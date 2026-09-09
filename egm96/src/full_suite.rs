use super::*;
use std::f64::consts::FRAC_PI_2;

#[test]
fn test_radgra_at_equator() {
    let lat = 0.0;
    let lon = 0.0;
    let mut rlat = 0.0;
    let mut gr = 0.0;
    let mut re = 0.0;
    radgra(lat, lon, &mut rlat, &mut gr, &mut re);
    assert!((rlat).abs() < 1e-9); // Geocentric latitude should be close to 0 at equator
    assert!((re - 6378137.0).abs() < 1e-9); // Geocentric radius should be close to semi-major axis
    assert!((gr - 9.7803253359).abs() < 1e-9); // Normal gravity at equator
}

#[test]
fn test_radgra_at_pole() {
    let lat = FRAC_PI_2;
    let lon = 0.0;
    let mut rlat = 0.0;
    let mut gr = 0.0;
    let mut re = 0.0;
    radgra(lat, lon, &mut rlat, &mut gr, &mut re);
    assert!((rlat - FRAC_PI_2).abs() < 1e-9); // Geocentric latitude should be close to pi/2 at pole
    assert!((re - 6356752.314245).abs() < 1e-6); // Geocentric radius at pole (semi-minor axis)
    assert!((gr - 9.8321863685).abs() < 1e-5); // Normal gravity at pole
}

#[test]
fn test_undulation_at_locations() {
    struct Check {
        lat: f64,
        lon: f64,
        geoid: f64,
    }

    let checks = [
        //Houston       :
        Check {
            lat: 29.7604,
            lon: -95.3698,
            geoid: -28.41,
        },
        //San Antonio   :
        Check {
            lat: 29.4241,
            lon: -98.4936,
            geoid: -26.52,
        },
        //San Diego     :
        Check {
            lat: 32.7157,
            lon: -117.1611,
            geoid: -35.22,
        },
        //Dallas        :
        Check {
            lat: 32.7767,
            lon: -96.797,
            geoid: -27.34,
        },
        //San Jose      :
        Check {
            lat: 37.3382,
            lon: -121.8863,
            geoid: -32.37,
        },
        //Los Angeles   :
        Check {
            lat: 34.0522,
            lon: -118.2437,
            geoid: -35.17,
        },
        //New York      :
        Check {
            lat: 40.7128,
            lon: -74.006,
            geoid: -32.73,
        },
        //San Francisco :
        Check {
            lat: 37.7749,
            lon: -122.4194,
            geoid: -32.17,
        },
        //Chicago       :
        Check {
            lat: 41.8781,
            lon: -87.6298,
            geoid: -33.93,
        },
        //London        :
        Check {
            lat: 51.5074,
            lon: 0.1278,
            geoid: 45.78,
        },
        //Paris         :
        Check {
            lat: 48.8566,
            lon: 2.3522,
            geoid: 44.61,
        },
        //Tokyo          :
        Check {
            lat: 35.6895,
            lon: 139.6917,
            geoid: 36.71,
        },
        //Philadelphia  :
        Check {
            lat: 40.05,
            lon: -75.45,
            geoid: -34.32,
        },
        //Phoenix       :
        Check {
            lat: 33.4484,
            lon: -112.074,
            geoid: -30.25,
        },
        //null island
        Check {
            lat: 0.0,
            lon: 0.0,
            geoid: 17.22,
        },
    ];

    for check in checks {
        let computed = egm96_compute_altitude_offset(check.lat, check.lon);
        let expected = check.geoid;
        let err = (computed - expected).abs();
        if err.is_nan() || err.is_infinite() || err > 0.5 {
            panic!(
                "Lat: {}, Lon: {}, Expected: {expected}, Computed: {computed}",
                check.lat, check.lon
            );
        }
    }
}

#[test]
fn test_wrap_degrees() {
    assert_eq!(wrap_degrees(0.0), 0.0);
    assert_eq!(wrap_degrees(179.8), 179.8);
    assert_eq!(wrap_degrees(-179.0), -179.0);
    assert_eq!(wrap_degrees(-181.0), 179.0);
    assert_eq!(wrap_degrees(190.0), -170.0);
    assert_eq!(wrap_degrees(-190.0), 170.0);
    assert_eq!(wrap_degrees(-190.0 - 360.0), 170.0);
    assert_eq!(wrap_degrees(360.0), 0.0);
    assert_eq!(wrap_degrees(540.0), -180.0);
    assert_eq!(wrap_degrees(-540.0), -180.0);
    assert_eq!(wrap_degrees(1000.0), -80.0);
    assert_eq!(wrap_degrees(-1000.0), 80.0);
}

#[cfg(feature = "raster_5_min")]
#[test]
fn test_5min_at_locations() {
    let _ = env_logger::builder().is_test(true).try_init();

    struct Check {
        lat: f64,
        lon: f64,
        geoid: f64,
    }

    let checks = [
        //Houston       :
        Check {
            lat: 29.7604,
            lon: -95.3698,
            geoid: -28.41,
        },
        //San Antonio   :
        Check {
            lat: 29.4241,
            lon: -98.4936,
            geoid: -26.52,
        },
        //San Diego     :
        Check {
            lat: 32.7157,
            lon: -117.1611,
            geoid: -35.22,
        },
        //Dallas        :
        Check {
            lat: 32.7767,
            lon: -96.797,
            geoid: -27.34,
        },
        //San Jose      :
        Check {
            lat: 37.3382,
            lon: -121.8863,
            geoid: -32.37,
        },
        //Los Angeles   :
        Check {
            lat: 34.0522,
            lon: -118.2437,
            geoid: -35.17,
        },
        //New York      :
        Check {
            lat: 40.7128,
            lon: -74.006,
            geoid: -32.73,
        },
        //San Francisco :
        Check {
            lat: 37.7749,
            lon: -122.4194,
            geoid: -32.17,
        },
        //Chicago       :
        Check {
            lat: 41.8781,
            lon: -87.6298,
            geoid: -33.93,
        },
        //London        :
        Check {
            lat: 51.5074,
            lon: 0.1278,
            geoid: 45.78,
        },
        //Paris         :
        Check {
            lat: 48.8566,
            lon: 2.3522,
            geoid: 44.61,
        },
        //Tokyo          :
        Check {
            lat: 35.6895,
            lon: 139.6917,
            geoid: 36.71,
        },
        //Philadelphia  :
        Check {
            lat: 40.05,
            lon: -75.45,
            geoid: -34.32,
        },
        //Phoenix       :
        Check {
            lat: 33.4484,
            lon: -112.074,
            geoid: -30.25,
        },
        //null island
        Check {
            lat: 0.0,
            lon: 0.0,
            geoid: 17.22,
        },
    ];

    for check in checks {
        let computed = egm96_raster_5_min_altitude_offset(check.lat, check.lon);
        let expected = check.geoid;
        let err = (computed - expected).abs();
        if err.is_nan() || err.is_infinite() || err > 0.5 {
            panic!(
                "Lat: {}, Lon: {}, Expected: {expected}, Computed: {computed}",
                check.lat, check.lon
            );
        }
    }
}

#[test]
#[cfg(feature = "raster_15_min")]
fn test_15min_at_locations() {
    let _ = env_logger::builder().is_test(true).try_init();

    struct Check {
        lat: f64,
        lon: f64,
        geoid: f64,
    }

    let checks = [
        //Houston       :
        Check {
            lat: 29.7604,
            lon: -95.3698,
            geoid: -28.41,
        },
        //San Antonio   :
        Check {
            lat: 29.4241,
            lon: -98.4936,
            geoid: -26.52,
        },
        //San Diego     :
        Check {
            lat: 32.7157,
            lon: -117.1611,
            geoid: -35.22,
        },
        //Dallas        :
        Check {
            lat: 32.7767,
            lon: -96.797,
            geoid: -27.34,
        },
        //San Jose      :
        Check {
            lat: 37.3382,
            lon: -121.8863,
            geoid: -32.37,
        },
        //Los Angeles   :
        Check {
            lat: 34.0522,
            lon: -118.2437,
            geoid: -35.17,
        },
        //New York      :
        Check {
            lat: 40.7128,
            lon: -74.006,
            geoid: -32.73,
        },
        //San Francisco :
        Check {
            lat: 37.7749,
            lon: -122.4194,
            geoid: -32.17,
        },
        //Chicago       :
        Check {
            lat: 41.8781,
            lon: -87.6298,
            geoid: -33.93,
        },
        //London        :
        Check {
            lat: 51.5074,
            lon: 0.1278,
            geoid: 45.78,
        },
        //Paris         :
        Check {
            lat: 48.8566,
            lon: 2.3522,
            geoid: 44.61,
        },
        //Tokyo
        Check {
            lat: 35.355,
            lon: 139.895,
            geoid: 47303.0 * 0.003 - 108.0,
        },
        //Philadelphia  :
        Check {
            lat: 40.05,
            lon: -75.45,
            geoid: -34.32,
        },
        //Phoenix       :
        Check {
            lat: 33.4484,
            lon: -112.074,
            geoid: -30.25,
        },
        //null island
        Check {
            lat: 0.0,
            lon: 0.0,
            geoid: 17.22,
        },
    ];

    for (i, check) in checks.iter().enumerate() {
        let computed = egm96_raster_15_min_altitude_offset(check.lat, check.lon);
        let expected = check.geoid;
        let err = (computed - expected).abs();
        if err.is_nan() || err.is_infinite() || err > 0.5 {
            panic!(
                "{i}, Lat: {}, Lon: {}, Expected: {expected}, Computed: {computed}",
                check.lat, check.lon
            );
        }
    }
}

// In these tests, we will use a simple 3x3 image.
// We define WIDTH = 3 and HEIGHT = 3.
// The pixel array is provided row-major order.
// For example, we fill the pixels with various values:
// Row 0: 100, 110, 120
// Row 1: 130, 140, 150
// Row 2: 160, 170, 180
//
// Also note that the interpolation function applies:
// value => value * SCALE + OFFSET, where SCALE = 0.003 and OFFSET = -108.
//
// Therefore:
// For pixel 100: 100 * 0.003 - 108 = -107.7
// For pixel 110: 110 * 0.003 - 108 = -107.67
// For pixel 120: 120 * 0.003 - 108 = -107.64
// For pixel 130: 130 * 0.003 - 108 = -107.61
// For pixel 140: 140 * 0.003 - 108 = -107.58
// For pixel 150: 150 * 0.003 - 108 = -107.55
// For pixel 160: 160 * 0.003 - 108 = -107.52
// For pixel 170: 170 * 0.003 - 108 = -107.49
// For pixel 180: 180 * 0.003 - 108 = -107.46

const WIDTH: usize = 3;
const HEIGHT: usize = 3;

// Define our test image pixels.
const PIXELS: [u16; WIDTH * HEIGHT] = [100, 110, 120, 130, 140, 150, 160, 170, 180];

// We'll use x_start = 0.0 and y_start = 0.0 with steps of 1.0.
const X_START: f64 = 0.0;
const Y_START: f64 = 0.0;
const X_STEP: f64 = 1.0;
const Y_STEP: f64 = 1.0;

// A helper function to compare two floating-point numbers within an epsilon.
fn approx_eq(a: f64, b: f64, epsilon: f64) -> bool {
    (a - b).abs() < epsilon
}

#[test]
fn test_exact_top_left() {
    // Testing a point that exactly maps to the top-left pixel.
    // lat = 0, lon = 0; therefore, x = 0, y = 0.
    let lat = 0.0;
    let lon = 0.0;
    let result = interpolate::<WIDTH, HEIGHT>(lat, lon, X_START, Y_START, X_STEP, Y_STEP, &PIXELS);
    // Expected value is the value of pixel at (0,0): 100 * 0.003 - 108 = -107.7
    let expected = 100.0 * 0.003 - 108.0;
    assert!(
        approx_eq(result, expected, 1e-6),
        "Expected {}, got {}",
        expected,
        result
    );
}

#[test]
fn test_exact_bottom_right() {
    // Testing a point that exactly maps to the bottom-right pixel.
    // For a 3x3 image, bottom-right pixel is at (2,2).
    let lat = 2.0;
    let lon = 2.0;
    let result = interpolate::<WIDTH, HEIGHT>(lat, lon, X_START, Y_START, X_STEP, Y_STEP, &PIXELS);
    // Expected value is the value of pixel at (2,2): 180 * 0.003 - 108 = -107.46
    let expected = 180.0 * 0.003 - 108.0;
    assert!(
        approx_eq(result, expected, 1e-6),
        "Expected {}, got {}",
        expected,
        result
    );
}

#[test]
fn test_center_interpolation() {
    // Testing a point that lies in the exact center of the top-left 2x2 block.
    // For lat = 0.5 and lon = 0.5, x, y = 0.5, so dx = 0.5 and dy = 0.5.
    //
    // The four neighboring pixel values:
    // top_left: (0,0) => 100 * 0.003 - 108 = -107.7
    // top_right: (1,0) => 110 * 0.003 - 108 = -107.67
    // bottom_left: (0,1) => 130 * 0.003 - 108 = -107.61
    // bottom_right: (1,1) => 140 * 0.003 - 108 = -107.58
    //
    // Interpolating:
    // top = -107.7 + 0.5 * ( -107.67 - (-107.7) ) = -107.7 + 0.5 * 0.03 = -107.685
    // bottom = -107.61 + 0.5 * ( -107.58 - (-107.61) ) = -107.61 + 0.5 * 0.03 = -107.595
    // result = top + 0.5 * (bottom - top) = -107.685 + 0.5 * (0.09) = -107.685 + 0.045 = -107.64
    let lat = 0.5;
    let lon = 0.5;
    let result = interpolate::<WIDTH, HEIGHT>(lat, lon, X_START, Y_START, X_STEP, Y_STEP, &PIXELS);
    let expected = -107.64; // Based on the interpolation above.
    assert!(
        approx_eq(result, expected, 1e-6),
        "Expected {}, got {}",
        expected,
        result
    );
}

#[test]
fn test_right_edge_clamping() {
    // Testing a point on the right boundary that may have an x coordinate equal to WIDTH - 1.
    // For lon = 2.0 and lat = 0.0, x = 2.0 and y = 0.0.
    // The surrounding indices become:
    // x0 = 2, x1 = 3 (clamped to 2); y0 = 0, y1 = 1.
    // All x values used in interpolation will refer to column index 2.
    let lat = 0.0;
    let lon = 2.0;
    let result = interpolate::<WIDTH, HEIGHT>(lat, lon, X_START, Y_START, X_STEP, Y_STEP, &PIXELS);
    // Expected value is the pixel at (2,0): 120 * 0.003 - 108 = -107.64.
    let expected = 120.0 * 0.003 - 108.0;
    assert!(
        approx_eq(result, expected, 1e-6),
        "Expected {}, got {}",
        expected,
        result
    );
}

#[test]
fn test_negative_coordinates_clamping() {
    // Testing a point with negative lat and lon.
    // For lat = -0.5 and lon = -0.5, the computed x and y are -0.5.
    // x0 = -1 and y0 = -1. After clamping, they become 0.
    // x1 and y1 are then also clamped to 0.
    //
    // Thus, all four neighbors are pixel (0,0).
    let lat = -0.5;
    let lon = -0.5;
    let result = interpolate::<WIDTH, HEIGHT>(lat, lon, X_START, Y_START, X_STEP, Y_STEP, &PIXELS);
    // Expected value is pixel at (0,0): 100 * 0.003 - 108 = -107.7.
    let expected = 100.0 * 0.003 - 108.0;
    assert!(
        approx_eq(result, expected, 1e-6),
        "Expected {}, got {}",
        expected,
        result
    );
}

#[test]
fn test_bottom_edge_interpolation() {
    // Testing a point on the lower edge where y is exactly at HEIGHT-1.
    // For lat = 2.0 and lon = 1.5, x = 1.5, y = 2.0.
    // x0 = 1, x1 = 2; y0 = 2, y1 = 3 (clamped to 2).
    // In this case the top and bottom rows are the same (row 2).
    let lat = 2.0;
    let lon = 1.5;
    let result = interpolate::<WIDTH, HEIGHT>(lat, lon, X_START, Y_START, X_STEP, Y_STEP, &PIXELS);
    // The interpolation is only in x.
    // For row 2, the pixels at x=1 and x=2 are:
    // pixel at (1,2): 170 * 0.003 - 108 = -107.49
    // pixel at (2,2): 180 * 0.003 - 108 = -107.46
    // x factor = 0.5, so expected = -107.49 + 0.5 * (-107.46 + 107.49) = -107.49 + 0.5 * 0.03 = -107.475
    let expected = -107.475;
    assert!(
        approx_eq(result, expected, 1e-6),
        "Expected {}, got {}",
        expected,
        result
    );
}

#[test]
fn test_extreme_lat_lon() {
    // Poles
    assert!(egm96_compute_altitude_offset(90.0, 0.0).is_finite());
    assert!(egm96_compute_altitude_offset(-90.0, 0.0).is_finite());

    // Antimeridian
    assert!(egm96_compute_altitude_offset(10.0, 180.0).is_finite());
    assert!(egm96_compute_altitude_offset(10.0, -180.0).is_finite());

    // Multiple wraps
    let a = egm96_compute_altitude_offset(10.0, 30.0);
    let b = egm96_compute_altitude_offset(10.0, 390.0);
    assert!((a - b).abs() < 1e-9);

    // Slightly out of bounds (lat clamp)
    let c = egm96_compute_altitude_offset(95.0, 0.0);
    let d = egm96_compute_altitude_offset(90.0, 0.0);
    assert!((c - d).abs() < 1e-9);
}

#[test]
fn test_nan_inputs() {
    assert!(egm96_compute_altitude_offset(f64::NAN, 0.0).is_nan());
    assert!(egm96_compute_altitude_offset(0.0, f64::NAN).is_nan());
}

#[test]
fn test_infinite_inputs() {
    assert!(egm96_compute_altitude_offset(f64::INFINITY, 0.0).is_finite());
    assert!(egm96_compute_altitude_offset(0.0, f64::INFINITY).is_nan());
}

#[test]
fn test_longitude_symmetry() {
    let a = egm96_compute_altitude_offset(45.0, 30.0);
    let b = egm96_compute_altitude_offset(45.0, 390.0);
    assert!((a - b).abs() < 1e-9);
}

#[test]
fn test_antipodal_points() {
    let a = egm96_compute_altitude_offset(10.0, 20.0);
    let b = egm96_compute_altitude_offset(-10.0, 200.0);
    assert!(a.is_finite());
    assert!(b.is_finite());
}

#[cfg(feature = "raster_5_min")]
#[test]
fn test_raster_5min_edges() {
    // First pixel
    let a = egm96_raster_5_min_altitude_offset(-90.0, 0.0);
    assert!(a.is_finite());

    // Last pixel
    let b = egm96_raster_5_min_altitude_offset(90.0, 359.999);
    assert!(b.is_finite());

    // Meridian crossing
    let c = egm96_raster_5_min_altitude_offset(10.0, -0.0001);
    let d = egm96_raster_5_min_altitude_offset(10.0, 359.9999);
    assert!((c - d).abs() < 1e-3);
}

#[cfg(feature = "raster_15_min")]
#[test]
fn test_raster_15min_edges() {
    let a = egm96_raster_15_min_altitude_offset(-90.0, 0.0);
    assert!(a.is_finite());

    let b = egm96_raster_15_min_altitude_offset(90.0, 359.999);
    assert!(b.is_finite());

    let c = egm96_raster_15_min_altitude_offset(10.0, -0.0001);
    let d = egm96_raster_15_min_altitude_offset(10.0, 359.9999);
    assert!((c - d).abs() < 1e-3);
}

#[test]
fn test_multithread_scratch() {
    use std::thread;

    let handles: Vec<_> = (0..16)
        .map(|_| thread::spawn(|| egm96_compute_altitude_offset(45.0, 45.0)))
        .collect();

    for h in handles {
        assert!(h.join().unwrap().is_finite());
    }
}

#[test]
fn test_random_fuzz() {
    for _ in 0..500 {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .subsec_nanos() as f64;

        let lat = (nanos / 1e9) * 180.0 - 90.0;
        let lon = (nanos / 1e9) * 360.0 - 180.0;

        let v = egm96_compute_altitude_offset(lat, lon);
        assert!(v.is_finite());
    }
}

#[test]
fn test_performance_sanity() {
    use std::time::Instant;

    let start = Instant::now();
    for _ in 0..100 {
        let _ = egm96_compute_altitude_offset(40.0, -74.0);
    }
    let elapsed = start.elapsed();

    assert!(elapsed.as_millis() < 2000);
}

#[test]
fn test_lat_clamping() {
    let a = egm96_compute_altitude_offset(90.0, 0.0);
    let b = egm96_compute_altitude_offset(120.0, 0.0);
    assert!((a - b).abs() < 1e-9);
}

#[test]
fn test_lon_wrapping() {
    let a = egm96_compute_altitude_offset(10.0, 20.0);
    let b = egm96_compute_altitude_offset(10.0, 380.0);
    assert!((a - b).abs() < 1e-9);
}

#[test]
fn test_feature_fallback() {
    let _a = egm96_altitude_offset(10.0, 20.0);
    let _b = egm96_compute_altitude_offset(10.0, 20.0);

    #[cfg(not(any(feature = "raster_5_min", feature = "raster_15_min")))]
    assert!((a - b).abs() < 1e-9);
}
