#include <stdio.h>

#include "golden_section_search.h"

double func1(double x) {
  return (x - 2) * (x - 2);
}

int main(void) {
  
  double lower_bound = 0;
  double upper_bound = 4;
  
  printf("min: %lf\n", find_univariate_min(lower_bound, upper_bound, func1));

  return 0;
  
}
