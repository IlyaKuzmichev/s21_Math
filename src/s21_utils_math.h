#ifndef SRC_S21_UTILS_MATH_H_
#define SRC_S21_UTILS_MATH_H_

#define S21_EPSILON 1e-15
#define S21_SQRT_E 1.648721270700128102982932876141575206L

long double s21_trunc(double x);
long double asin_first_quater(double x);
long double s21_fast_int_pow(long double base, unsigned long long exp);
long double s21_fast_pow(long double base, long double exp);

#endif  // SRC_S21_UTILS_MATH_H_