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

//use std::f64::consts::PI;
//use std::f64::sqrt;
//use std::f64::sin;
//use std::f64::cos;
//
//use crate::egm96_data::EGM96_DATA;
//
//* ************************************************************************** */
//
//const _COEFFS: usize = 65341; //!< Size of correction and harmonic coefficients arrays (361*181)
//const _NMAX: usize = 360;   //!< Maximum degree and orders of harmonic coefficients.
//const _361: usize = 361;
//
//* ************************************************************************** */
//
//fn hundu(
//    p: &[f64],
//    sinml: &[f64],
//    cosml: &[f64],
//    gr: f64,
//    re: f64,
//) -> f64 {
//    // WGS 84 gravitational constant in m³/s² (mass of Earth’s atmosphere included)
//    const GM: f64 = 0.3986004418e15;
//    // WGS 84 datum surface equatorial radius
//    const AE: f64 = 6378137.0;
//
//    let ar = AE / re;
//    let mut arn = ar;
//    let mut ac = 0.0;
//    let mut a = 0.0;
//
//    let mut k = 3;
//    for n in 2..=_NMAX {
//        arn *= ar;
//        k += 1;
//        let mut sum = p[k] * EGM96_DATA[k][2];
//        let mut sumc = p[k] * EGM96_DATA[k][0];
//
//        for m in 1..=n {
//            k += 1;
//            let tempc = EGM96_DATA[k][0] * cosml[m] + EGM96_DATA[k][1] * sinml[m];
//            let temp = EGM96_DATA[k][2] * cosml[m] + EGM96_DATA[k][3] * sinml[m];
//            sumc += p[k] * tempc;
//            sum += p[k] * temp;
//        }
//        ac += sumc;
//        a += sum * arn;
//    }
//    ac += EGM96_DATA[1][0] + (p[2] * EGM96_DATA[2][0]) + (p[3] * (EGM96_DATA[3][0] * cosml[1] + EGM96_DATA[3][1] * sinml[1]));
//
//    // Add haco = ac/100 to convert height anomaly on the ellipsoid to the undulation
//    // Add -0.53m to make undulation refer to the WGS84 ellipsoid
//
//    (a * GM) / (gr * re) + (ac / 100.0) - 0.53
//}
//
//fn dscml(rlon: f64, sinml: &mut [f64], cosml: &mut [f64]) {
//    let a = sin(rlon);
//    let b = cos(rlon);
//
//    sinml[1] = a;
//    cosml[1] = b;
//    sinml[2] = 2.0 * b * a;
//    cosml[2] = 2.0 * b * b - 1.0;
//
//    for m in 3..=_NMAX {
//        sinml[m] = 2.0 * b * sinml[m - 1] - sinml[m - 2];
//        cosml[m] = 2.0 * b * cosml[m - 1] - cosml[m - 2];
//    }
//}
//
//*!
// * \param m: order.
// * \param theta: Colatitude (radians).
// * \param rleg: Normalized legendre function.
// *
// * This subroutine computes all normalized legendre function in 'rleg'.
// * The dimensions of array 'rleg' must be at least equal to nmax+1.
// * All calculations are in double precision.
// *
// * Original programmer: Oscar L. Colombo, Dept. of Geodetic Science the Ohio State University, August 1980.
// * ineiev: I removed the derivatives, for they are never computed here.
// */
//fn legfdn(m: usize, theta: f64, rleg: &mut [f64]) {
//    static mut DRTS: [f64; 1302] = [0.0; 1302];
//    static mut DIRT: [f64; 1302] = [0.0; 1302];
//    static mut RLNN: [f64; _361 + 1] = [0.0; _361 + 1];
//    static mut IR: i32 = 0; // TODO 'ir' must be set to zero before the first call to this sub.
//
//    let nmax1 = _NMAX + 1;
//    let nmax2p = (2 * _NMAX) + 1;
//    let m1 = m + 1;
//    let m2 = m + 2;
//    let m3 = m + 3;
//
//    unsafe {
//        if IR == 0 {
//            IR = 1;
//            for n in 1..=nmax2p {
//                DRTS[n] = sqrt(n as f64);
//                DIRT[n] = 1.0 / DRTS[n];
//            }
//        }
//
//        let cothet = cos(theta);
//        let sithet = sin(theta);
//
//        // compute the legendre functions
//        RLNN[1] = 1.0;
//        RLNN[2] = sithet * DRTS[3];
//        for n1 in 3..=m1 {
//            let n = n1 - 1;
//            let n2 = 2 * n;
//            RLNN[n1] = DRTS[n2 + 1] * DIRT[n2] * sithet * RLNN[n];
//        }
//
//        match m {
//            1 => {
//                rleg[2] = RLNN[2];
//                rleg[3] = DRTS[5] * cothet * rleg[2];
//            }
//            0 => {
//                rleg[1] = 1.0;
//                rleg[2] = cothet * DRTS[3];
//            }
//            _ => {}
//        }
//        rleg[m1] = RLNN[m1];
//
//        if m2 <= nmax1 {
//            rleg[m2] = DRTS[m1 * 2 + 1] * cothet * rleg[m1];
//            if m3 <= nmax1 {
//                for n1 in m3..=nmax1 {
//                    let n = n1 - 1;
//                    if (m == 0 && n < 2) || (m == 1 && n < 3) {
//                        continue;
//                    }
//                    let n2 = 2 * n;
//                    rleg[n1] = DRTS[n2 + 1] * DIRT[n + m] * DIRT[n - m]
//                        * (DRTS[n2 - 1] * cothet * rleg[n1 - 1]
//                            - DRTS[n + m - 1] * DRTS[n - m - 1] * DIRT[n2 - 3] * rleg[n1 - 2]);
//                }
//            }
//        }
//    }
//}
//
//*!
// * \param lat: Latitude in radians.
// * \param lon: Longitude in radians.
// * \param re: Geocentric radius.
// * \param rlat: Geocentric latitude.
// * \param gr: Normal gravity (m/sec²).
// *
// * This subroutine computes geocentric distance to the point, the geocentric
// * latitude, and an approximate value of normal gravity at the point based the
// * constants of the WGS84(g873) system are used.
// */
//fn radgra(lat: f64, lon: f64, rlat: &mut f64, gr: &mut f64, re: &mut f64) {
//    const A: f64 = 6378137.0;
//    const E2: f64 = 0.00669437999013;
//    const GEQT: f64 = 9.7803253359;
//    const K: f64 = 0.00193185265246;
//    let t1 = sin(lat).powi(2);
//    let n = A / sqrt(1.0 - (E2 * t1));
//    let t2 = n * cos(lat);
//    let x = t2 * cos(lon);
//    let y = t2 * sin(lon);
//    let z = (n * (1.0 - E2)) * sin(lat);
//
//    *re = sqrt((x * x) + (y * y) + (z * z)); // compute the geocentric radius
//    *rlat = (z / sqrt((x * x) + (y * y))).atan(); // compute the geocentric latitude
//    *gr = GEQT * (1.0 + (K * t1)) / sqrt(1.0 - (E2 * t1)); // compute normal gravity (m/sec²)
//}
//
//*!
// * \brief Compute the geoid undulation from the EGM96 potential coefficient model, for a given latitude and longitude.
// * \param lat: Latitude in radians.
// * \param lon: Longitude in radians.
// * \return The geoid undulation / altitude offset (in meters).
// */
//fn undulation(lat: f64, lon: f64) -> f64 {
//    let mut p = [0.0; _COEFFS + 1];
//    let mut sinml = [0.0; _361 + 1];
//    let mut cosml = [0.0; _361 + 1];
//    let mut rleg = [0.0; _361 + 1];
//
//    let mut rlat = 0.0;
//    let mut gr = 0.0;
//    let mut re = 0.0;
//    let nmax1 = _NMAX + 1;
//
//    // compute the geocentric latitude, geocentric radius, normal gravity
//    radgra(lat, lon, &mut rlat, &mut gr, &mut re);
//    let rlat_colat = (PI / 2.0) - rlat;
//
//    for j in 1..=nmax1 {
//        let m = j - 1;
//        legfdn(m, rlat_colat, &mut rleg);
//        for i in j..=nmax1 {
//            p[(((i - 1) * i) / 2) + m + 1] = rleg[i];
//        }
//    }
//    dscml(lon, &mut sinml, &mut cosml);
//
//    hundu(&p, &sinml, &cosml, gr, re)
//}
//
//* ************************************************************************** */
//
//pub fn egm96_compute_altitude_offset(lat_deg: f64, lon_deg: f64) -> f64 {
//    const RAD: f64 = PI / 180.0;
//    undulation(lat_deg * RAD, lon_deg * RAD)
//}
//
//* ************************************************************************** */
//
//#[cfg(test)]
//mod tests {
//    use super::*;
//    use std::f64::consts::FRAC_PI_2;
//
//    // Example usage and a simple test (you might need more comprehensive tests)
//    #[test]
//    fn test_undulation_at_zero() {
//        // Approximate value at 0, 0 should be around 47-51m based on EGM96
//        let lat = 0.0;
//        let lon = 0.0;
//        let undulation_value = undulation(lat, lon);
//        println!("Undulation at (0, 0): {}", undulation_value);
//        assert!((undulation_value - 49.0).abs() < 5.0); // Allow for some tolerance
//    }
//
//    #[test]
//    fn test_egm96_compute_altitude_offset() {
//        let lat_deg = 0.0;
//        let lon_deg = 0.0;
//        let offset = egm96_compute_altitude_offset(lat_deg, lon_deg);
//        println!("Altitude offset at (0, 0): {}", offset);
//        assert!((offset - 49.0).abs() < 5.0);
//    }
//
//    #[test]
//    fn test_radgra_at_equator() {
//        let lat = 0.0;
//        let lon = 0.0;
//        let mut rlat = 0.0;
//        let mut gr = 0.0;
//        let mut re = 0.0;
//        radgra(lat, lon, &mut rlat, &mut gr, &mut re);
//        assert!((rlat).abs() < 1e-9); // Geocentric latitude should be close to 0 at equator
//        assert!((re - 6378137.0).abs() < 1e-9); // Geocentric radius should be close to semi-major axis
//        assert!((gr - 9.7803253359).abs() < 1e-9); // Normal gravity at equator
//    }
//
//    #[test]
//    fn test_radgra_at_pole() {
//        let lat = FRAC_PI_2;
//        let lon = 0.0;
//        let mut rlat = 0.0;
//        let mut gr = 0.0;
//        let mut re = 0.0;
//        radgra(lat, lon, &mut rlat, &mut gr, &mut re);
//        assert!((rlat - FRAC_PI_2).abs() < 1e-9); // Geocentric latitude should be close to pi/2 at pole
//        assert!((re - 6356752.314245).abs() < 1e-6); // Geocentric radius at pole (semi-minor axis)
//        assert!((gr - 9.8321863685).abs() < 1e-9); // Normal gravity at pole
//    }
//}
