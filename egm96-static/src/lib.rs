#![allow(clippy::needless_return)]
#![allow(clippy::missing_safety_doc)]

/// Takes a lat, lon and returns the undulation
#[unsafe(no_mangle)]
pub unsafe extern "C" fn egm96_compute_altitude_offset(lat: f64, lon: f64) -> f64 {
    return egm96::egm96_compute_altitude_offset(lat, lon);
}
