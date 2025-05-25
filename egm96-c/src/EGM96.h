#pragma once

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/**
 * Takes a lat, lon and returns the undulation, this is uses the most accurate
 * method available, depending on build options.
 * - If 5 arc-minute raster is available at build time it will use it
 * - If 15 arc-minute raster is available at build time it will use it
 * - Otherwise will use spherical harmonic model that is slower but uses very little memory
 */
double egm96_altitude_offset(double lat, double lon);

/**
 * Takes a lat, lon and returns the undulation
 */
double egm96_compute_altitude_offset(double lat, double lon);

/**
 * Uses 15 arc-minute raster, lookup with interpolation
 */
double egm96_raster_15_min_altitude_offset(double lat, double lon);

/**
 * Uses 5 arc-minute raster, lookup with interpolation
 */
double egm96_raster_5_min_altitude_offset(double lat, double lon);

#ifdef __cplusplus
}  // extern "C"
#endif  // __cplusplus
