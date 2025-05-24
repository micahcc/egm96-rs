mod algo;
mod egm96_data;

pub use algo::egm96_compute_altitude_offset;
pub use algo::egm96_compute_altitude_offset_harmonics;

#[cfg(feature = "raster_15_min")]
pub use algo::egm96_compute_altitude_offset_raster_15_min;

#[cfg(feature = "raster_5_min")]
pub use algo::egm96_compute_altitude_offset_raster_5_min;
