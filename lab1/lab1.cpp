#include <cmath>
#include <iostream>

using namespace std;

double my_exp(double x, int accuracy) {
  double result = 1.0; 
  double a = 1;
  for (int n = 1; n <= accuracy; ++n) {
    a = a*x/n;
    result += a;
  }

  return result;
}

int main() {
  int number, accuracy;
/*
  cout << "Please input number: ";
  cin >> number;
  cout << "Please input Required accuracy ";
  cin >> accuracy;
*/
  accuracy = 1500000000;// 1 500 000 000
  number = 1;
  cout << "my_exp = " << my_exp(number, accuracy) << endl;
  return 0;
}