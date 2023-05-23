#include "s21_math.h"

#include "s21_utils_math.h"

union ieee754double {
  double number;
  struct {
    unsigned long long mantissa : 52;
    unsigned exponent : 11;
    unsigned sign : 1;
  };
};

static const unsigned int taylor_limitation = 500000;

static const int EXPONENT_ALIGNMENT = 1023;
static const int MANTISSA_SIZE = 52;
static const unsigned long long POSITIVE_MASK = 0xFFFFFFFFFFFFFFFF;

int s21_abs(int x) {
  if (x < 0) {
    x *= -1;
  }
  return x;
}

long double s21_fabs(double x) {
  long double result = (long double)x;
  if (result < 0) {
    result *= -1.;
  }
  return result;
}

long double s21_ceil(double x) {
  long double truncated = s21_trunc(x);
  if (x > 0 && truncated >= 0 && truncated != (long double)x) {
    truncated++;
  }
  return truncated;
}

long double s21_floor(double x) {
  long double truncated = s21_trunc(x);
  if (x < 0 && truncated <= 0 && truncated != (long double)x) {
    truncated--;
  }
  return truncated;
}

long double s21_exp(double x) {
  long double result = 0;
  if (s21_is_nan(x)) {
    result = S21_NAN;
  } else if (s21_is_inf(x)) {
    result = x < 0 ? 0 : S21_INF;
  }
  long double taylor_series_element = 1.;

  long double denom = 0.;
  size_t i = 0;
  long double fx = s21_fabs(x);
  while (i < taylor_limitation &&
         !(s21_is_inf(x) || s21_is_inf(result) || s21_is_nan(result)) &&
         taylor_series_element > S21_EPSILON) {
    result += taylor_series_element;
    i++;
    denom++;
    taylor_series_element *= fx / denom;
  }
  if (!s21_is_inf(x)) {
    result += taylor_series_element;
  }
  if (x < 0 && result != 0) {
    result = 1. / result;
  }
  return result;
}

long double s21_log(double x) {
  long double result = 0.;

  if (x < 0 || s21_is_nan(x)) {
    result = S21_NAN;
  } else if (s21_is_inf(x)) {
    result = x;
  } else if (x < 1e-30) {
    result = -S21_INF;
  } else {
    while (x >= S21_E) {
      x /= S21_E;
      result += 1.;
    }
    if (x > 2.0L) {
      x /= S21_SQRT_E;
      result += 0.5;
    }
    while (x < 0.5L) {
      x *= S21_E;
      result -= 1.;
    }

    x--;
    long double numerator = x;
    long double denominator = 1;
    long double taylor_series_elem = x;
    for (unsigned int i = 0;
         i < taylor_limitation && s21_fabs(taylor_series_elem) > S21_EPSILON;
         i++) {
      result += taylor_series_elem;
      denominator++;
      numerator *= -x;
      taylor_series_elem = numerator / denominator;
    }
    result += taylor_series_elem;
  }
  return result;
}

long double s21_fast_int_pow(long double base, unsigned long long exp) {
  long double result = 1.;
  while (exp > 0) {
    if (exp & 1) {
      result *= base;
    }
    base *= base;
    exp >>= 1;
  }
  return result;
}

long double s21_fast_pow(long double base, long double exp) {
  long double result = 1.;
  while (exp > (long double)ULLONG_MAX) {
    result *= s21_fast_int_pow(base, ULLONG_MAX);
    exp -= (long double)ULLONG_MAX;
  }
  result *= s21_fast_int_pow(base, exp);
  return result;
}

long double s21_pow(double base, double exp) {
  long double x = (long double)base;
  long double y = (long double)exp;
  long double result = 0.;
  if (0. == y) {
    result = 1.;
  } else if (1. == y) {
    result = x;
  } else if (0. == x) {
    result = y < 0 ? S21_INF : 0;
  } else if (1. == x) {
    result = 1.;
  } else if (-1 == y) {
    result = 1. / x;
  } else if (s21_is_nan(x) || s21_is_nan(y)) {
    result = S21_NAN;
  } else if (y == -S21_INF) {
    result = 0.;
  } else if (s21_is_inf(x) || y == S21_INF) {
    result = S21_INF;
  } else {
    long double absExp = s21_fabs(y);
    long double expIntPart = s21_trunc(absExp);
    long double expFractPart = absExp - expIntPart;
    result = s21_exp(expFractPart * s21_log(x)) * s21_fast_pow(x, expIntPart);
    if (y < 0) {
      result = 1. / result;
    }
  }
  return result;
}

long double s21_trunc(double x) {
  union ieee754double x_repr = {0};
  x_repr.number = x;
  int exponent = x_repr.exponent - EXPONENT_ALIGNMENT;
  if (exponent < 0) {
    x_repr.number = 0;
  } else {
    int fractional_part = MANTISSA_SIZE - exponent;
    if (fractional_part > 0) {
      x_repr.mantissa &= POSITIVE_MASK << fractional_part;
    }
  }
  return (long double)x_repr.number;
}

long double s21_fmod(double x, double y) {
  long double result = S21_NAN;

  if (y && !(s21_is_nan(x)) && !(s21_is_nan(y)) && !s21_is_inf(x) &&
      !s21_is_inf(y)) {
    long double round_result = s21_trunc(x / y);
    result = x - round_result * y;
  } else if (!(x != x) && !s21_is_inf(x) && s21_is_inf(y)) {
    result = x;
  }
  return result;
}

long double s21_sin(double x) {
  long double result = S21_NAN;

  if (!s21_is_nan(x) && !s21_is_inf(x)) {
    long double taylor_member = 0;
    int sign = 1;

    if (x < 0) {
      sign = -1;
    }
    x = s21_fmod(s21_fabs(x), 2 * S21_PI);
    if (x > S21_PI) {
      x -= S21_PI;
      sign *= -1;
    }
    if (x > S21_PI / 2) {
      x = S21_PI - x;
    }
    result = taylor_member = x;
    for (unsigned long int i = 3; s21_fabs(taylor_member) > S21_EPSILON;
         i += 2) {
      taylor_member *= (-1) * x * x / (long double)(i * (i - 1));
      result += taylor_member;
    }
    result *= sign;
  }
  return result;
}

long double s21_cos(double x) {
  long double result = S21_NAN;

  if (!s21_is_nan(x) && !s21_is_inf(x)) {
    long double taylor_member = 0;
    int sign = 1;
    x = s21_fmod(s21_fabs(x), 2 * S21_PI);
    if (x > S21_PI) {
      x -= S21_PI;
      sign *= -1;
    }
    if (x > S21_PI / 2) {
      x = S21_PI - x;
      sign *= -1;
    }
    result = 1 - x * x / 2.;
    taylor_member = -x * x / 2.;
    for (unsigned long int i = 4; s21_fabs(taylor_member) > S21_EPSILON;
         i += 2) {
      taylor_member *= (-1) * x * x / (long double)(i * (i - 1));
      result += taylor_member;
    }
    result *= sign;
  }
  return result;
}

long double s21_tan(double x) {
  long double result = S21_NAN;

  if (!s21_is_nan(x) && !s21_is_inf(x)) {
    if (s21_fabs(s21_cos(x)) < S21_EPSILON) {
      result = s21_sin(x) * S21_INF;
    } else {
      result = s21_sin(x) / s21_cos(x);
    }
  }
  return result;
}

long double asin_first_quater(double x) {
  long double result = x;
  long double taylor_member = x;

  for (unsigned long int i = 2; s21_fabs(taylor_member) > S21_EPSILON; i += 2) {
    taylor_member *= i * (i - 1) * (i - 1) * x * x / (i * i * (i + 1));
    result += taylor_member;
  }
  return result;
}

long double s21_asin(double x) {
  long double result = S21_NAN;

  if (!s21_is_nan(x) && x >= -1 && x <= 1) {
    if (1 - s21_fabs(x) < S21_EPSILON) {
      if (x < 0) {
        result = -S21_PI / 2;
      } else {
        result = S21_PI / 2;
      }
    } else {
      int sign = (x > 0) ? 1 : -1;
      if (s21_fabs(x) < 0.5) {
        result = asin_first_quater(s21_fabs(x));
      } else {
        result = S21_PI / 2 - asin_first_quater(s21_sqrt(1 - x * x));
      }
      result *= sign;
    }
  }
  return result;
}

long double s21_acos(double x) {
  long double result = S21_NAN;

  if (!s21_is_nan(x) && x >= -1 && x <= 1) {
    result = S21_PI / 2 - s21_asin(x);
  }
  return result;
}

long double s21_atan(double x) {
  long double result = S21_NAN;

  if (!s21_is_nan(x)) {
    if (s21_is_inf(x)) {
      if (x > 0) {
        result = S21_PI / 2;
      } else {
        result = -1 * S21_PI / 2;
      }
    } else if (s21_fabs(x) == 1) {
      result = x * S21_PI / 4;
    } else {
      int reverse_flag = 0;
      int sign = 1;
      if (s21_fabs(x) > 1) {
        reverse_flag = 1;
        x = 1 / x;
        if (x < 0) {
          sign = -1;
        }
      }
      long double taylor_member = x;
      result = x;
      for (unsigned long int i = 3; s21_fabs(taylor_member) > S21_EPSILON;
           i += 2) {
        taylor_member *= (-1) * x * x * (i - 2) / i;
        result += taylor_member;
      }
      if (reverse_flag) {
        result = sign * S21_PI / 2 - result;
      }
    }
  }
  return result;
}

long double s21_sqrt(double x) {
  long double result = S21_NAN;

  if (x == 0.) {
    result = 0;
  }
  if (x > 0 && s21_is_inf(x)) {
    result = S21_INF;
  }
  if (x > 0 && !s21_is_inf(x) && !s21_is_nan(x)) {
    long double prev_result = 0;
    result = x / 2;
    do {
      prev_result = result;
      result = (result + x / result) / 2;
    } while (s21_fabs(result - prev_result) > S21_EPSILON);
  }
  return result;
}
