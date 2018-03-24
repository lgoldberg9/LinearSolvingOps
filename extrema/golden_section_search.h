#ifndef __GOLDEN_SECTION_SEARCH_H__
#define __GOLDEN_SECTION_SEARCH_H__

#include <stdint.h>
#include <stdbool.h>
#include <math.h>

void set_tolerance(double new_tol);

void set_max_iters(double new_max_iter);

double find_univariate_min(double lower_bound, double upper_bound, double (*extrema_func)(double));

double find_univariate_max(double lower_bound, double upper_bound, double (*extrema_func)(double));

#endif
