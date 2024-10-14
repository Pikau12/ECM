#include <cmath>
#include <iostream>

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
  std::cout << my_exp(number, accuracy) << std::endl;
  return 0;
}
