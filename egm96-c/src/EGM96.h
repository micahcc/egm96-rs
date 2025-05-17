#pragma once

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/**
 * Takes a lat, lon and returns the undulation
 */
double egm96_compute_altitude_offset(double lat, double lon);

#ifdef __cplusplus
}  // extern "C"
#endif  // __cplusplus
