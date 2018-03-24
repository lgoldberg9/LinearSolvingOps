#include <stdio.h>

#include "golden_section_search.h"

double func1(double x) {
  return (x - 2) * (x - 2);
}

double func2(double x) {
  return 4 - x * x;
}

int main(void) {
  
  double lower_bound = 0;
  double upper_bound = 4;
  
  printf("min: %lf\n", find_univariate_min(lower_bound, upper_bound, func1));
  printf("max: %lf\n", find_univariate_max(-1, 3, func2));

  return 0;
  
}
