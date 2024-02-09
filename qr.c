#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

int64_t mul_scale(int64_t a, int64_t b);
// int32_t div_scale(int32_t a, int32_t b);

void print_matrix_float (int m, int n, double A[m][n]) {
  for (int8_t i = 0; i < m; i++) {
    printf("[ ");
    for (int8_t j = 0; j < n; j++)
      printf("%9.3f ", A[i][j]);
    printf("]\n");
  }
}

void print_matrix_int32 (int m, int n, int64_t A[m][n]) {
  for (int8_t i = 0; i < m; i++) {
    printf("[ ");
    for (int8_t j = 0; j < n; j++)
      printf("%6lld, ", A[i][j]);
    printf("];\n");
  }
}

double product_float (double *a, double *b, int n) {
  double ret = 0;
  for (int8_t k = 0; k < n; k++)
    ret += a[k] * b[k];
  return ret;
}

double norm_float (double *a, int n) {
  double ret = 0;
  for (int8_t k = 0; k < n; k++)
    ret += a[k] * a[k];
  return sqrt(ret);
}

void transpose_float (int n, double a[][n], int col, double *ret) {
  for (int i = 0; i < n; i++)
    ret[i] = a[i][col];
}

void QR_float(int m, int n, double A[m][n], double Q[m][m], double R[n][n]) {
  double Q_t[m];
  double w_float[m];
  double w_norm = 0;
  for (int8_t j = 0; j < n; j++) {
    for (int8_t k = 0; k < m; k++)
      w_float[k] = A[k][j];
    for (int8_t i = 0; i < j; i++) {
      transpose_float(n, Q, i, Q_t);
      R[i][j] = product_float(Q_t, w_float, n);
      for (int8_t k = 0; k < m; k++)
        w_float[k] -= R[i][j] * Q[k][i];
    }
    w_norm = norm_float(w_float, n);
    for (int8_t k = 0; k < m; k++)
      Q[k][j] = w_float[k] / w_norm;
    R[j][j] = w_norm;
  }
}

void matrix_mul_float (int n, double a[][n], double b[][n], double c[][n]) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      c[i][j] = 0;
      for (int k = 0; k < n; k++)
        c[i][j] += a[i][k] * b[k][j];
    }
  }
}

void matrix_mul_scale (int n, int64_t a[][n], int64_t b[][n], int64_t c[][n]) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      c[i][j] = 0;
      for (int k = 0; k < n; k++)
        c[i][j] += mul_scale(a[i][k], b[k][j]);
    }
  }
}

int64_t mul_scale (int64_t a, int64_t b) {
  return (a * b) >> 8;
}

int32_t div_scale(int32_t a, int32_t b) {
  return (a << 8) / b;
}

int32_t my_sqrt(int64_t data) {
  int64_t l = 1, r = data >> 1;
  while (l <= r) {
    int64_t mid = (l + r) / 2;
    int64_t tmp = mid * mid;
    if (tmp == data)
      return mid;
    else if (tmp < data)
      l = mid + 1;
    else
      r = mid - 1;
  }
  return l;
}

int32_t norm_scale(int m, int64_t *w) {
  int32_t tmp = 0;
  for (int8_t k = 0; k < m; k++)
    tmp += mul_scale(w[k], w[k]);
  double tmp_sqrt = sqrt((double)tmp / 256);
  int32_t tmp_sqrt_int = my_sqrt(tmp);
  tmp_sqrt *= 256;
  tmp_sqrt_int *= 16;
  // printf("%lf, %d\n", tmp_sqrt, tmp_sqrt_int);
  // return (int32_t)tmp_sqrt;
  return tmp_sqrt_int;
}

int64_t product (int n, int64_t *a, int64_t *b) {
  int64_t ret = 0;
  for (int8_t k = 0; k < n; k++)
    ret += mul_scale(a[k], b[k]);
  return ret;
}

void transpose (int n, int64_t a[][n], int8_t col, int64_t *ret) {
  for (int8_t i = 0; i < n; i++)
    ret[i] = a[i][col];
}

void QR_fixed (int m, int n, int64_t A_scale[m][n], int64_t Q_scale[m][m], int64_t R_scale[n][n]) {
  int64_t w[m];
  int64_t Q_t[m];
  int64_t w_norm;
  for (int8_t j = 0; j < n; j++) {
    for (int8_t k = 0; k < m; k++)
      w[k] = A_scale[k][j];
    for (int8_t i = 0; i < j; i++) {
      transpose(n, Q_scale, i, Q_t);
      R_scale[i][j] = product(n, Q_t, w); // Q is << 8, product will >> 8, so R will be correct
      for (int8_t k = 0; k < m; k++)
        w[k] -= mul_scale(R_scale[i][j], Q_scale[k][i]);
    }
    w_norm = norm_scale(m, w);
    for (int8_t k = 0; k < m; k++)
      Q_scale[k][j] = div_scale(w[k], w_norm); // Q will << 8
    R_scale[j][j] = w_norm;
  }
}

int main() {
  int m = 4;
  int n = 4;
  double A_float[4][4] = {{356,  -8,  -4,   6}, 
                          { -8, 364,  -2,  -4}, 
                          {  4,  -6, 356,  -8}, 
                          {  2,   4,  -8, 364}};
  double Q_float[4][4] = {0};
  double R_float[4][4] = {0};

  printf("A_float = \n");
  print_matrix_float(m, n, A_float);

  QR_float(m, n, A_float, Q_float, R_float);

  printf("Q_float = \n");
  print_matrix_float(m, n, Q_float);

  printf("R_float = \n");
  print_matrix_float(m, n, R_float);

  double result_float[4][4];
  matrix_mul_float(m, Q_float, R_float, result_float);

  printf("QR_float = \n");
  print_matrix_float(m, n, result_float);

  printf("----------------------------------\n");

  m = 4;
  n = 4;

  int32_t ch0[2] = {364, -4};
  int32_t ch1[2] = {0, 4};
  int32_t ch2[2] = {0, -2};
  int32_t ch3[2] = {358, 0};

  int32_t A[4][4] = {{ch0[0], ch1[0], -ch0[1], -ch1[1]}, 
                     {ch2[0], ch3[0], -ch2[1], -ch3[1]}, 
                     {ch0[1], ch1[1], ch0[0], ch1[0]}, 
                     {ch2[1], ch3[1], ch2[0], ch3[0]}};

  int64_t A_scale[m][n];
  int64_t Q_scale[m][n];
  int64_t R_scale[m][n];
  memset(Q_scale, 0, sizeof(Q_scale));
  memset(R_scale, 0, sizeof(Q_scale));

  for (int8_t i = 0; i < m; i++)
    for (int8_t j = 0; j < n; j++)
      A_scale[i][j] = A[i][j];
  
  printf("A_scale = \n");
  print_matrix_int32(m, n, A_scale);

  QR_fixed(m, n, A_scale, Q_scale, R_scale);
 
  printf("Q_scale = \n");
  print_matrix_int32(m, n, Q_scale);

  printf("R_scale = \n");
  print_matrix_int32(m, n, R_scale);

  int64_t result[4][4];
  matrix_mul_scale(m, Q_scale, R_scale, result);

  printf("QR_scale = \n");
  print_matrix_int32(m, n, result);

  return 0;
}
