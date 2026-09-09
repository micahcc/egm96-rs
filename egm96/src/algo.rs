#![allow(clippy::all)]
#![allow(clippy::needless_return)]

/*
 * Copyright (c) 2006 D.Ineiev <ineiev@yahoo.co.uk>
 * Copyright (c) 2020 Emeric Grange <emeric.grange@gmail.com>
 * Copyright (c) 2025 Micah Chambers <micahc.vt@gmail.com>
 *
 * This software is provided 'as-is', without any express or implied warranty.
 * In no event will the authors be held liable for any damages arising from
 * the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 * claim that you wrote the original software. If you use this software
 * in a product, an acknowledgment in the product documentation would be
 * appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 * misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 **/

/*
 * This program is designed for the calculation of a geoid undulation at a point
 * whose latitude and longitude is specified.
 *
 * This program is designed to be used with the constants of EGM96 and those of
 * the WGS84(g873) system. The undulation will refer to the WGS84 ellipsoid.
 *
 * It's designed to use the potential coefficient model EGM96 and a set of
 * spherical harmonic coefficients of a correction term.
 * The correction term is composed of several different components, the primary
 * one being the conversion of a height anomaly to a geoid undulation.
 * The principles of this procedure were initially described in the paper:
 * - use of potential coefficient models for geoid undulation determination using
 * a spherical harmonic representation of the height anomaly/geoid undulation
 * difference by R.H. Rapp, Journal of Geodesy, 1996.
 *
 * This program is a modification of the program described in the following report:
 * - a fortran program for the computation of gravimetric quantities from high
 * degree spherical harmonic expansions, Richard H. Rapp, report 334, Department
 * of Geodetic Science and Surveying, the Ohio State University, Columbus, 1982
 **/

use std::cell::RefCell;
use std::f64::consts::PI;
use std::sync::OnceLock;

use crate::egm96_data::EGM96_DATA;

/****************************************************************************/

// Maximum degree and orders of harmonic coefficients
const NMAX: usize = 360;
const NMAX1: usize = 361;
const N361: usize = 361;
// Size of correction and harmonic coefficients arrays (361*181)
const COEFFS: usize = 65341;

#[cfg(feature = "raster_15_min")]
const EGM96_15_BYTES: &[u8] = include_bytes!(concat!(env!("OUT_DIR"), "/egm96-15.png"));

#[cfg(feature = "raster_5_min")]
const EGM96_5_BYTES: &[u8] = include_bytes!(concat!(env!("OUT_DIR"), "/egm96-5.png"));

/***************************************************************************/

struct Egm96Scratch {
    p: Box<[f64; COEFFS + 1]>,
    sinml: Box<[f64; N361 + 1]>,
    cosml: Box<[f64; N361 + 1]>,
    rleg: Box<[f64; N361 + 1]>,
    rlnn: Box<[f64; N361 + 1]>,
}

impl Egm96Scratch {
    fn new() -> Self {
        Self {
            p: Box::new([0.0; COEFFS + 1]),
            sinml: Box::new([0.0; N361 + 1]),
            cosml: Box::new([0.0; N361 + 1]),
            rleg: Box::new([0.0; N361 + 1]),
            rlnn: Box::new([0.0; N361 + 1]),
        }
    }
}

thread_local! {
    static SCRATCH: RefCell<Egm96Scratch> = RefCell::new(Egm96Scratch::new());
}

/// Compute sine and cosine values for the given longitude
fn dscml(rlon: f64, sinml: &mut [f64; N361 + 1], cosml: &mut [f64; N361 + 1]) {
    let a = rlon.sin();
    let b = rlon.cos();

    sinml[1] = a;
    cosml[1] = b;
    sinml[2] = 2.0 * b * a;
    cosml[2] = 2.0 * b * b - 1.0;

    for m in 3..=NMAX {
        sinml[m] = 2.0 * b * sinml[m - 1] - sinml[m - 2];
        cosml[m] = 2.0 * b * cosml[m - 1] - cosml[m - 2];
    }
}

/// Compute height undulation based on coefficients
fn hundu(
    p: &[f64; COEFFS + 1],
    sinml: &[f64; N361 + 1],
    cosml: &[f64; N361 + 1],
    gr: f64,
    re: f64,
) -> f64 {
    // WGS 84 gravitational constant in m^3/s^2 (mass of Earth's atmosphere included)
    const GM: f64 = 0.3986004418e15;
    // WGS 84 datum surface equatorial radius
    const AE: f64 = 6378137.0;

    let ar = AE / re;
    let mut arn = ar;
    let mut ac = 0.0;
    let mut a = 0.0;

    let mut k = 3;
    for n in 2..=NMAX {
        arn *= ar;
        k += 1;
        let mut sum = p[k] * EGM96_DATA[k][2] as f64;
        let mut sumc = p[k] * EGM96_DATA[k][0] as f64;

        for m in 1..=n {
            k += 1;
            let tempc = EGM96_DATA[k][0] as f64 * cosml[m] + EGM96_DATA[k][1] as f64 * sinml[m];
            let temp = EGM96_DATA[k][2] as f64 * cosml[m] + EGM96_DATA[k][3] as f64 * sinml[m];
            sumc += p[k] * tempc;
            sum += p[k] * temp;
        }
        ac += sumc;
        a += sum * arn;
    }
    ac += EGM96_DATA[1][0] as f64
        + (p[2] * EGM96_DATA[2][0] as f64)
        + (p[3] * (EGM96_DATA[3][0] as f64 * cosml[1] + EGM96_DATA[3][1] as f64 * sinml[1]));

    // Add haco = ac/100 to convert height anomaly on the ellipsoid to the undulation
    // Add -0.53m to make undulation refer to the WGS84 ellipsoid

    ((a * GM) / (gr * re)) + (ac / 100.0) - 0.53
}

/// Computes geocentric distance, geocentric latitude, and approximate normal gravity
fn radgra(lat: f64, lon: f64, rlat: &mut f64, gr: &mut f64, re: &mut f64) {
    const A: f64 = 6378137.0;
    const E2: f64 = 0.00669437999013;
    const GEQT: f64 = 9.7803253359;
    const K: f64 = 0.00193185265246;
    let t1 = lat.sin().powi(2);
    let n = A / (1.0 - (E2 * t1)).sqrt();
    let t2 = n * lat.cos();
    let x = t2 * lon.cos();
    let y = t2 * lon.sin();
    let z = (n * (1.0 - E2)) * lat.sin();

    *re = (x * x + y * y + z * z).sqrt(); // compute the geocentric radius
    *rlat = (z / (x * x + y * y).sqrt()).atan(); // compute the geocentric latitude
    *gr = GEQT * (1.0 + (K * t1)) / (1.0 - (E2 * t1)).sqrt(); // compute normal gravity (m/sec²)
}

/// Compute the geoid undulation from the EGM96 model
fn undulation(lat: f64, lon: f64) -> f64 {
    static DRTS_DIRT: OnceLock<([f64; 1301], [f64; 1301])> = OnceLock::new();
    let (drts, dirt) = DRTS_DIRT.get_or_init(|| {
        let nmax2p = (2 * NMAX) + 1;
        let mut drts = [0.0; 1301];
        let mut dirt = [0.0; 1301];
        for n in 1..=nmax2p {
            drts[n] = (n as f64).sqrt();
            dirt[n] = 1.0 / drts[n];
        }
        (drts, dirt)
    });

    SCRATCH.with(|scratch| {
        let mut s = scratch.borrow_mut();
        let Egm96Scratch {
            p,
            sinml,
            cosml,
            rleg,
            rlnn,
        } = &mut *s;

        let mut rlat = 0.0;
        let mut gr = 0.0;
        let mut re = 0.0;

        radgra(lat, lon, &mut rlat, &mut gr, &mut re);
        rlat = (PI / 2.0) - rlat;
        let cothet = rlat.cos();
        let sithet = rlat.sin();

        rlnn[1] = 1.0;
        rlnn[2] = sithet * drts[3];
        for j in 1..=NMAX1 {
            let m = j - 1;
            let m1 = m + 1;
            for n1 in 3..=m1 {
                let n = n1 - 1;
                let n2 = 2 * n;
                rlnn[n1] = drts[n2 + 1] * dirt[n2] * sithet * rlnn[n];
            }
        }

        for j in 1..=NMAX1 {
            let m = j - 1;
            let m1 = m + 1;
            let m2 = m + 2;
            let m3 = m + 3;

            if m == 0 {
                rleg[1] = 1.0;
                rleg[2] = cothet * drts[3];
            } else if m == 1 {
                rleg[2] = rlnn[2];
                rleg[3] = drts[5] * cothet * rleg[2];
            }
            rleg[m1] = rlnn[m1];

            if m2 <= NMAX1 {
                rleg[m2] = drts[m1 * 2 + 1] * cothet * rleg[m1];
                for n1 in m3..=NMAX1 {
                    let n = n1 - 1;
                    if (!m == 0 && n < 2) || (m == 1 && n < 3) {
                        continue;
                    }
                    let n2 = 2 * n;
                    rleg[n1] = drts[n2 + 1]
                        * dirt[n + m]
                        * dirt[n - m]
                        * (drts[n2 - 1] * cothet * rleg[n1 - 1]
                            - drts[n + m - 1] * drts[n - m - 1] * dirt[n2 - 3] * rleg[n1 - 2]);
                }
            }

            for i in j..=NMAX1 {
                p[((i - 1) * i) / 2 + m + 1] = rleg[i];
            }
        }
        dscml(lon, sinml, cosml);

        hundu(p, sinml, cosml, gr, re)
    })
}

fn wrap_degrees(mut degrees: f64) -> f64 {
    degrees += 180.0;
    degrees = degrees.rem_euclid(360.0);
    degrees - 180.0
}

pub fn egm96_compute_altitude_offset(lat: f64, lon: f64) -> f64 {
    let lon = wrap_degrees(lon);
    let lat = lat.clamp(-90.0, 90.0);
    undulation(lat.to_radians(), lon.to_radians())
}

#[allow(unused)]
fn interpolate<const WIDTH: usize, const HEIGHT: usize>(
    lat: f64,
    lon: f64,
    x_start: f64,
    y_start: f64,
    x_step: f64,
    y_step: f64,
    pixels: &[u16],
) -> f64 {
    const SCALE: f64 = 0.003;
    const OFFSET: f64 = -108.0;

    // X_geo = GT(0) + X_pixel * GT(1)
    // Y_geo = GT(3) + Y_line * GT(5)

    // X_pixel = (X_geo - GT(0)) / GT(1)
    // Y_line  = (Y_geo - GT(3)) / GT(5)

    let x = (lon - x_start) / x_step;
    let y = (lat - y_start) / y_step;

    // Determine the integer coordinates surrounding the point.
    let x0 = x.floor() as isize;
    let y0 = y.floor() as isize;
    let x1 = x0 + 1;
    let y1 = y0 + 1;

    // Clamp indices so they remain within the image bounds.
    // Convert the floating-point location differences to factors.
    let x0_clamped = x0.clamp(0, (WIDTH - 1) as isize) as usize;
    let y0_clamped = y0.clamp(0, (HEIGHT - 1) as isize) as usize;
    let x1_clamped = x1.clamp(0, (WIDTH - 1) as isize) as usize;
    let y1_clamped = y1.clamp(0, (HEIGHT - 1) as isize) as usize;

    // Compute the fractional part (distance between the point and the floor indices).
    let dx = (x - x0 as f64).clamp(0.0, 1.0);
    let dy = (y - y0 as f64).clamp(0.0, 1.0);

    // Retrieve the values at the four neighboring pixels.
    let top_left = pixels[y0_clamped * WIDTH + x0_clamped] as f64 * SCALE + OFFSET;
    let top_right = pixels[y0_clamped * WIDTH + x1_clamped] as f64 * SCALE + OFFSET;
    let bottom_left = pixels[y1_clamped * WIDTH + x0_clamped] as f64 * SCALE + OFFSET;
    let bottom_right = pixels[y1_clamped * WIDTH + x1_clamped] as f64 * SCALE + OFFSET;

    // Interpolate in the x direction on the top and bottom rows.
    let top = top_left + dx * (top_right - top_left);
    let bottom = bottom_left + dx * (bottom_right - bottom_left);

    // interpolate in the y direction between the top and bottom interpolated values.
    top + dy * (bottom - top)
}

fn load_image<const WIDTH: usize, const HEIGHT: usize>(bytes: &[u8]) -> Vec<u16> {
    let decoder = png::Decoder::new(bytes);
    let mut reader = decoder.read_info().expect("Failed to check info");
    // Allocate the output buffer.
    let mut buf = vec![0; reader.output_buffer_size()];
    // Read the next frame. An Atiff might contain multiple frames.
    let info = reader.next_frame(&mut buf).expect("Failed to get frame");

    // Grab the bytes of the image.
    buf.truncate(info.buffer_size());
    assert!(buf.len() == WIDTH * HEIGHT * 2);

    let mut out = vec![0; WIDTH * HEIGHT];
    for row in 0..HEIGHT {
        for col in 0..WIDTH {
            let index = 2 * (col + row * WIDTH);
            out[row * WIDTH + col] =
                u16::from_be_bytes(buf[index..(index + 2)].try_into().expect("pair"));
        }
    }

    out
}

#[cfg(feature = "raster_5_min")]
pub fn egm96_raster_5_min_altitude_offset(lat: f64, lon: f64) -> f64 {
    const WIDTH: usize = 4320;
    const HEIGHT: usize = 2161;
    static IMAGE: OnceLock<Vec<u16>> = OnceLock::new();
    let image = IMAGE.get_or_init(|| load_image::<WIDTH, HEIGHT>(EGM96_5_BYTES));

    let mut lon = wrap_degrees(lon);
    if lon < 0.0 {
        lon += 360.0;
    }

    let lat = lat.clamp(-90.0, 90.0);
    // from https://gdal.org/en/stable/tutorials/geotransforms_tut.html
    // X_geo = GT(0) + X_pixel * GT(1) + Y_line * GT(2)
    // Y_geo = GT(3) + X_pixel * GT(4) + Y_line * GT(5)
    // GT(0) = -0.04166666666666666
    // GT(1) = 0.08333333333333333
    // GT(2) = 0
    // GT(3) = 90.04166666666666666
    // GT(4) = 0
    // GT(5) = -0.08333333333333333
    interpolate::<WIDTH, HEIGHT>(
        lat,
        lon,
        -0.04166666666666666,
        90.04166666666666666,
        0.08333333333333333,
        -0.08333333333333333,
        image,
    )
}

#[cfg(feature = "raster_15_min")]
pub fn egm96_raster_15_min_altitude_offset(lat: f64, lon: f64) -> f64 {
    const WIDTH: usize = 1440;
    const HEIGHT: usize = 721;
    static IMAGE: OnceLock<Vec<u16>> = OnceLock::new();
    let image = IMAGE.get_or_init(|| load_image::<WIDTH, HEIGHT>(EGM96_15_BYTES));

    let mut lon = wrap_degrees(lon);
    if lon < 0.0 {
        lon += 360.0;
    }

    let lat = lat.clamp(-90.0, 90.0);

    // from https://gdal.org/en/stable/tutorials/geotransforms_tut.html
    // X_geo = GT(0) + X_pixel * GT(1) + Y_line * GT(2)
    // Y_geo = GT(3) + X_pixel * GT(4) + Y_line * GT(5)
    // GT(0) = -0.125
    // GT(1) = 0.25
    // GT(2) = 0
    // GT(3) = 90.12500000000000000
    // GT(4) = 0
    // GT(5) = -0.25000000000000000
    interpolate::<WIDTH, HEIGHT>(lat, lon, -0.125, 90.125, 0.25, -0.25, image)
}

/// Public function to compute altitude offset using EGM96 model
pub fn egm96_altitude_offset(lat: f64, lon: f64) -> f64 {
    #[cfg(feature = "raster_5_min")]
    {
        // highest res
        egm96_raster_5_min_altitude_offset(lat, lon)
    }

    #[cfg(all(feature = "raster_15_min", not(feature = "raster_5_min")))]
    {
        // medium res
        egm96_raster_15_min_altitude_offset(lat, lon)
    }

    // slow, but little memory
    #[cfg(all(not(feature = "raster_15_min"), not(feature = "raster_5_min")))]
    {
        egm96_compute_altitude_offset(lat, lon)
    }
}

#[cfg(test)]
#[path = "full_suite.rs"]
mod full_suite;
