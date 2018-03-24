#include "golden_section_search.h"
#include <stdio.h>

#define GOLDEN_RATIO (1 + sqrt(5)) / 2
#define REVERSE_GR GOLDEN_RATIO - 1

static double TOL = 1E-6;
static uint32_t MAX_ITER = 1000;

void set_tolerance(double new_tol) {
  TOL = new_tol;
}

void set_max_iters(double new_max_iter) {
  MAX_ITER = new_max_iter;
}

static double midpoint(double a, double b) {
  return a + (b - a) / 2;
}

static bool decide_min(double a, double b) {
  return a < b;
}

static bool decide_max(double a, double b) {
  return a > b;
}

/**
 * 
 *
 *
 */
static double golden_section_search(double domain_lower_bound, double domain_upper_bound, double (*extrema_func)(double), bool (*comparator)(double, double)) {

  // Set an iteration counter to eventually terminate from the loop
  uint32_t iter = 0;

  double max_dist = domain_upper_bound - domain_lower_bound;
  double lower_guestimate = domain_lower_bound + (1 - (REVERSE_GR)) * max_dist;
  double upper_guestimate = domain_lower_bound + (REVERSE_GR) * max_dist;
  
  while (fabs(upper_guestimate - lower_guestimate) >= TOL && iter < MAX_ITER) {
    if (comparator(extrema_func(lower_guestimate), extrema_func(upper_guestimate))) {
      // Supposing the extrema is on the 'left' side of the interval
      domain_upper_bound = upper_guestimate;
      upper_guestimate = lower_guestimate;
      max_dist = domain_upper_bound - domain_lower_bound;
      lower_guestimate = domain_lower_bound + (1 - (REVERSE_GR)) * max_dist;
    } else {
      // Supposing the extrema is on the 'right' side of the interval
      domain_lower_bound = lower_guestimate;
      lower_guestimate = upper_guestimate;
      max_dist = domain_upper_bound - domain_lower_bound;
      upper_guestimate = domain_lower_bound + (REVERSE_GR) * max_dist;
    }
    iter++;
  }
  return midpoint(domain_upper_bound, domain_lower_bound);
}

double find_univariate_min(double lower_bound, double upper_bound, double (*extrema_func)(double)) {
  return golden_section_search(lower_bound, upper_bound, extrema_func, decide_min);
}

double find_univariate_max(double lower_bound, double upper_bound, double (*extrema_func)(double)) {
  return golden_section_search(lower_bound, upper_bound, extrema_func, decide_max);
}
