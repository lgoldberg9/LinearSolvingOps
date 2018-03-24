#include "golden_section_search.h"

#define GOLDEN_RATIO (1 + sqrt(5)) / 2;

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
  return a <= b;
}

static bool decide_max(double a, double b) {
  return a >= b;
}

static double golden_section_search(double lower_bound, double upper_bound, double (*extrema_func)(double), bool (*comparator)(double, double)) {

  double inward_step = (upper_bound - lower_bound) / GOLDEN_RATIO;
  double candidate_lower_bound = upper_bound - inward_step;
  double candidate_upper_bound = lower_bound + inward_step;
  
  while (fabs(candidate_upper_bound - candidate_lower_bound) >= TOL) {
    if (comparator(extrema_func(candidate_lower_bound), extrema_func(candidate_upper_bound))) {
      upper_bound = candidate_upper_bound;
    } else {
      lower_bound = candidate_lower_bound;
    }
    inward_step = (upper_bound - lower_bound) / GOLDEN_RATIO;
    candidate_lower_bound = upper_bound - inward_step;
    candidate_upper_bound = lower_bound + inward_step;
  }
  return midpoint(upper_bound, lower_bound);
}

double find_univariate_min(double lower_bound, double upper_bound, double (*extrema_func)(double)) {
  return golden_section_search(lower_bound, upper_bound, extrema_func, decide_min);
}

double find_univariate_max(double lower_bound, double upper_bound, double (*extrema_func)(double)) {
  return golden_section_search(lower_bound, upper_bound, extrema_func, decide_max);
}
