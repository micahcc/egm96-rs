#![allow(clippy::needless_return)]
#![allow(clippy::missing_safety_doc)]

/// Takes a lat, lon and returns the undulation, this is uses the most accurate
/// method available, depending on build options.
/// - If 5 arc-minute raster is available at build time it will use it
/// - If 15 arc-minute raster is available at build time it will use it
/// - Otherwise will use spherical harmonic model that is slower but uses very little memory
#[unsafe(no_mangle)]
pub unsafe extern "C" fn egm96_altitude_offset(lat: f64, lon: f64) -> f64 {
    return egm96::egm96_altitude_offset(lat, lon);
}

/// Takes a lat, lon and returns the undulation
#[unsafe(no_mangle)]
pub unsafe extern "C" fn egm96_compute_altitude_offset(lat: f64, lon: f64) -> f64 {
    return egm96::egm96_compute_altitude_offset(lat, lon);
}

/// Uses 15 arc-minute raster, lookup with interpolation
#[cfg(feature = "raster_15_min")]
#[allow(unused)]
#[unsafe(no_mangle)]
pub unsafe extern "C" fn egm96_raster_15_min_altitude_offset(lat: f64, lon: f64) -> f64 {
    return egm96::egm96_raster_15_min_altitude_offset(lat, lon);
}

/// Uses 5 arc-minute raster, lookup with interpolation
#[allow(unused)]
#[cfg(feature = "raster_5_min")]
#[unsafe(no_mangle)]
pub unsafe extern "C" fn egm96_raster_5_min_altitude_offset(lat: f64, lon: f64) -> f64 {
    return egm96::egm96_raster_5_min_altitude_offset(lat, lon);
}
