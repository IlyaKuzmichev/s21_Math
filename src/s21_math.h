/* ==========================================================================*/
/* Copyright 2022 ©
 *
 * This file contains Original Code created by
 * Ilya Kuzmichev aka wilmerno,
 * Kirill Safin aka tabathae.
 *
 * The Original Code and all software developed in the process of
 * participation on learning by experimental programming educational method.
 * The whole methodology was developed and distributed by
 * Autonomous non-profit organization «School 21» (ANO «School 21»).
 *
 * Redistribution and use of this file, its parts, or entire project
 * are permitted by confirmation of its original creators.
 */
/* ==========================================================================*/

#ifndef SRC_S21_MATH_H_
#define SRC_S21_MATH_H_

#include <float.h>
#include <limits.h>
#include <stddef.h>

#define S21_PI 3.14159265358979323846
#define S21_NAN (0.0f / 0.0f)
#define S21_INF (1.0f / 0.0f)
#define S21_E 2.718281828459045235360287471352662500L
#define s21_is_inf(x) (x == S21_INF || x == -S21_INF)
#define s21_is_nan(x) (x != x)

int s21_abs(int x);
long double s21_ceil(double x);
long double s21_floor(double x);
long double s21_exp(double x);
long double s21_log(double x);
long double s21_fabs(double x);
long double s21_pow(double base, double exp);
long double s21_sin(double x);
long double s21_cos(double x);
long double s21_tan(double x);
long double s21_acos(double x);
long double s21_asin(double x);
long double s21_atan(double x);
long double s21_fmod(double x, double y);
long double s21_sqrt(double x);

#endif  // SRC_S21_MATH_H_
