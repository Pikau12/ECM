#include <math.h>
#include <stdio.h>

double my_exp(double x, long int accuracy) {
  double result = 1.0; 
  double a = 1;
  for (long int n = 1; n <= accuracy; ++n) {
    a = a*x/n;
    result += a;
  }

  return result;
}

int main() {
  long int number, accuracy;
  accuracy = 4000000000;
  number = 1;
  printf("my_exp = %f \n", my_exp(number, accuracy));
  return 0;
}
