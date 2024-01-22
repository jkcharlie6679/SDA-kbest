#include <stdio.h>
#include <stdint.h>
#include <string.h>

typedef struct c16_t{
  int16_t r;
  int16_t i;
} c16_t;

typedef struct c32_t{
  int32_t r;
  int32_t i;
} c32_t;



int8_t xRange[3][8] = {
  {-1, 1, 0, 0, 0, 0, 0, 0},
  {-3, -1, 1, 3, 0, 0, 0, 0},
  {-7, -5, -3, -1, 1, 3, 5, 7}
};

int32_t mul_scale(int32_t a, int32_t b);

void print_matrix_int16 (int m, int n, int16_t A[m][n]) {
  for (int8_t i = 0; i < m; i++) {
    printf("[ ");
    for (int8_t j = 0; j < n; j++)
      printf("%6d, ", A[i][j]);
    printf("];\n");
  }
}


void print_matrix_int32 (int m, int n, int32_t A[m][n]) {
  for (int8_t i = 0; i < m; i++) {
    printf("[ ");
    for (int8_t j = 0; j < n; j++)
      printf("%6d, ", A[i][j]);
    printf("];\n");
  }
}

void matrix_mul_scale (int n, int16_t a[][n], int16_t b[][n], int16_t c[][n]) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      c[i][j] = 0;
      for (int k = 0; k < n; k++)
        c[i][j] += mul_scale(a[i][k], b[k][j]);
    }
  }
}

int32_t mul_scale (int32_t a, int32_t b) {
  return (a * b) >> 8;
}

int32_t div_scale(int32_t a, int32_t b) {
  return (a << 8) / b;
}

int64_t my_sqrt(int64_t data) {
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

int32_t norm_scale(int m, int32_t *w) {
  int32_t tmp = 0;
  for (int8_t k = 0; k < m; k++)
    tmp += mul_scale(w[k], w[k]);
  int32_t tmp_sqrt_int = my_sqrt(tmp);
  tmp_sqrt_int *= 16;
  return tmp_sqrt_int;
}

int32_t product (int n, int32_t *a, int32_t *b) {
  int32_t ret = 0;
  for (int8_t k = 0; k < n; k++)
    ret += mul_scale(a[k], b[k]);
  return ret;
}

void transpose (int n, int16_t a[][n], int8_t col, int32_t *ret) {
  for (int8_t i = 0; i < n; i++)
    ret[i] = a[i][col];
}

void QR_fixed (int m, int n, int16_t A_scale[m][n], int16_t Q_scale[m][m], int32_t R_scale[n][n]) {
  int32_t w[m];
  int32_t Q_t[m];
  int32_t w_norm;
  for (int8_t j = 0; j < n; j++) {
    for (int8_t k = 0; k < m; k++)
      w[k] = A_scale[k][j];
    for (int8_t i = 0; i < j; i++) {
      transpose(n, Q_scale, i, Q_t);
      R_scale[i][j] = product(n, Q_t, w); // Q is << 8, product will >> 8, so R will be correct
      for (int8_t k = 0; k < m; k++)
        w[k] -= mul_scale(R_scale[i][j], (int32_t)Q_scale[k][i]);
    }
    w_norm = norm_scale(m, w);
    for (int8_t k = 0; k < m; k++)
      Q_scale[k][j] = div_scale(w[k], w_norm); // Q will << 8
    R_scale[j][j] = w_norm;
  }
}

int16_t first_norm (c16_t a, c16_t b) {
  int16_t res;
  int32_t a32 = 0, b32 = 0;
  a32 = a.r * a.r + a.i * a.i;
  b32 = b.r * b.r + b.i * b.i;
  res = my_sqrt(a32 + b32);
  return res;
}

c16_t cmplx_mul (c16_t a, c16_t b) {
  c16_t result = {0, 0};
  result.r = (a.r * b.r - a.i * b.i) >> 8;
  result.i = (a.r * b.i + a.i * b.r) >> 8;
  return result;
}

c16_t cmplx_add (c16_t a, c16_t b) {
  c16_t result = {0};
  result.r = a.r + b.r;
  result.i = a.i + b.i;
  return  result;
}
c16_t cmplx_sub (c16_t a, c16_t b) {
  c16_t result = {0};
  result.r = a.r - b.r;
  result.i = a.i - b.i;
  return  result;
}
c16_t cmplx_scale (c16_t a, int b) {
  if (b >= 0) {
    a.r = (int16_t)(((int32_t)a.r * b) >> 8);
    a.i = (int16_t)(((int32_t)a.i * b) >> 8);
  } 
  return a;
}

int binarySearch(int64_t *a, int8_t *res, int l, int r, int target) {
  if (r <= l)
    return (a[target] > a[res[l]] ? l + 1 : l);
  int mid = (l + r) / 2;
  
  if (a[target] == a[res[mid]])
    return mid + 1;
  if (a[target] > a[res[mid]])
    return binarySearch(a, res, mid + 1, r, target);
  return binarySearch(a, res, l, mid - 1, target);
}

void sort_ped (int k_size, int64_t *a, int8_t *res) {
  int j = 0;
  res[0] = 0;
  for (int k = 1; k < k_size; k++) {
    j = k - 1;
    int selected = k;
    int loc = binarySearch(a, res, 0, j, selected);
    while (j >= loc) {
      res[j + 1] = res[j];
      j--;
    }
    res[j + 1] = selected;
  }
}

int main() {

  int m = 4;
  int n = 4;
  int k_best = 4;
  int8_t order = 6;

  c16_t ch0 = {360, 0};
  c16_t ch1 = {4, 0};
  c16_t ch2 = {2, -4};
  c16_t ch3 = {358, -4};

  c16_t data_in[2] = {{572/12, -3123/12}, {-3145/12, 4429/12}};

  // ZF
  c16_t y_zf[2] = {0};
  c16_t tmp_ch = {0};
  c16_t zf0 = cmplx_mul(ch3, data_in[0]);
  tmp_ch.r = -ch1.r;
  tmp_ch.i = -ch1.i;
  c16_t zf1 = cmplx_mul(tmp_ch, data_in[1]);
  tmp_ch.r = -ch2.r;
  tmp_ch.i = -ch2.i;
  c16_t zf2 = cmplx_mul(tmp_ch, data_in[0]);
  c16_t zf3 = cmplx_mul(ch0, data_in[1]);
  
  y_zf[0] = cmplx_add(zf0, zf1);
  y_zf[1] = cmplx_add(zf2, zf3);

  printf("y_zf[0] (%d, %d)\n", y_zf[0].r, y_zf[0].i);
  printf("y_zf[1] (%d, %d)\n", y_zf[1].r, y_zf[1].i);


  int16_t raduis = 0; 
  c16_t data_yzf_0 = cmplx_sub(data_in[0], y_zf[0]);
  c16_t data_yzf_1 = cmplx_sub(data_in[1], y_zf[1]);
  switch (order) {
    case 2:
      data_yzf_0 = cmplx_scale(data_yzf_0, 362);
      data_yzf_1 = cmplx_scale(data_yzf_1, 362);
      raduis = first_norm(data_yzf_0, data_yzf_1);
      data_in[0] = cmplx_scale(data_in[0], 362);
      data_in[1] = cmplx_scale(data_in[1], 362);
      break;
    case 4:
      data_yzf_0 = cmplx_scale(data_yzf_0, 810);
      data_yzf_1 = cmplx_scale(data_yzf_1, 810);
      raduis = first_norm(data_yzf_0, data_yzf_1);
      data_in[0] = cmplx_scale(data_in[0], 810);
      data_in[1] = cmplx_scale(data_in[1], 810);
      break;
    case 6:
      data_yzf_0 = cmplx_scale(data_yzf_0, 1659);
      data_yzf_1 = cmplx_scale(data_yzf_1, 1659);
      raduis = first_norm(data_yzf_0, data_yzf_1);
      data_in[0] = cmplx_scale(data_in[0], 1659);
      data_in[1] = cmplx_scale(data_in[1], 1659);
      break;
    default:
      printf("Modulation order %d is not support.\n", order);
  }
  printf("data_in[0]: (%d, %d)\n", data_in[0].r, data_in[0].i);
  printf("data_in[1]: (%d, %d)\n", data_in[1].r, data_in[1].i);
  printf("raduis: %d\n", raduis);
  
  int16_t H[4][4] = {{ch0.r, ch1.r, -ch0.i, -ch1.i}, 
                     {ch2.r, ch3.r, -ch2.i, -ch3.i}, 
                     {ch0.i, ch1.i,  ch0.r,  ch1.r}, 
                     {ch2.i, ch3.i,  ch2.r,  ch3.r}};
  int16_t Q[4][4] = {0};
  int32_t R[4][4] = {0};

  QR_fixed(m, n, H, Q, R);
  printf("Q = \n");
  print_matrix_int16(m, n, Q);
  printf("R = \n");
  print_matrix_int32(m, n, R);
  int xRangeLen = 2 * xRange[order/2 -1 ][0] / -2 + 1;
  int k_size = 2 * (xRange[order/2 -1 ][0] / -2 + 1);
  printf("%d\n", k_size);

  int16_t level = 0;
  int16_t xHat[4][32] = {0};
  int32_t dp[32] = {0};
  int64_t ped[32] = {0};
  int16_t py_extend[4][32];

  for (int i = 0; i < 8; i++) {
    xHat[level][i] = xRange[order / 2 - 1][i];
  }
 
  int16_t data_in16[4] = {data_in[0].r, data_in[1].r, data_in[0].i, data_in[1].i};
  int16_t yR[4] = {0}; 

  for (int i = 0; i < 4; i++) {
    int32_t tmp = 0;
    for (int j = 0; j < 4; j++) {
      tmp += (int32_t)Q[j][i] * (int32_t)data_in16[j] / 256;
    }
    yR[i] = (int16_t)(tmp);
    printf("yR[%d]: %d\n", i, yR[i]);
  }
  
  for (int i = 0; i < k_size; i++)
    dp[i] = raduis;
  for (int i = 0; i < m; i++) 
    for (int j = 0; j < k_size; j++)
      py_extend[i][j] = yR[i];
  int dimp = 3;
  int dim = 2;
  for (int i = 0; i < 4; i++) {
    int k_best_cur = k_size;
    if (k_best_cur > k_best)
      k_best_cur = k_best;
    int64_t tmp_w[32] = {INT64_MAX};
    for (int k = 0; k < k_size; k++) {
      tmp_w[k] = (int64_t)py_extend[dimp][k] - (int64_t)R[dimp][dimp] * xHat[level][k]; 
      tmp_w[k] = tmp_w[k] * tmp_w[k];
    }
    int64_t tmp_dp[32] = {INT64_MAX};
    for (int k = 0; k < k_size; k++) {
      tmp_dp[k] = (int64_t)dp[k] * dp[k] - tmp_w[k]; 
      tmp_dp[k] = my_sqrt(tmp_dp[k]);
      ped[k] += tmp_w[k];
      printf("%lld, ", tmp_w[k]);
    }
    printf("\n");
    int8_t res[32];
    sort_ped(k_size, ped, res);
    printf("fin: %d\n", res[0]);
    if (i == 3) {
      // for (int l = 0; l < 4; l++) {
      //   for (int k = 0; k < 32; k++) {
      //     printf("%2d, ", xHat[l][k]);
      //   }
      //   printf("\n");
      // }
      printf("Ans: ");
      for (int l = 3; l >= 0; l--)
        printf("%2d, ", xHat[l][res[0]]);
      break;
    }
    k_size = k_best_cur * xRangeLen;
    // for (int k = 0; k < k_size; k++) {
    //   printf("%d, ", res[k]);
    // }
    int16_t xHat_tmp[4][32] = {0};
    // int32_t dp_tmp[32] = {0};
    int64_t ped_tmp[32] = {0};
    int16_t py_extend_tmp[4][32];
    memcpy(xHat_tmp, xHat, sizeof(xHat));
    // memcpy(dp_tmp, dp, sizeof(dp));
    memcpy(ped_tmp, ped, sizeof(ped));
    memcpy(py_extend_tmp, py_extend, sizeof(py_extend));

    for (int k = 0; k < k_best_cur; k++) {
      for (int j = 0; j <= dim; j++) {
        py_extend_tmp[j][res[k]] = py_extend_tmp[j][res[k]] - R[j][dimp] * xHat_tmp[level][res[k]];
      }
      for (int j = 0; j < xRangeLen; j++) {
        for (int l = 0; l < 4; l++) {
          py_extend[l][xRangeLen * k + j] = py_extend_tmp[l][res[k]];
        }
        xHat[level][xRangeLen * k + j] = xHat_tmp[level][res[k]];
        xHat[level + 1][xRangeLen * k + j] = xRange[order / 2 - 1][j];
        ped[xRangeLen * k + j] = ped_tmp[res[k]];
        dp[xRangeLen * k + j] = tmp_dp[res[k]]; 
        // printf("%3d, ", xHat[level][xRangeLen * k + j]);
      }
      // printf("\n");
      // for (int j = 0; j < xRangeLen; j++)
      //   printf("%3d, ", xHat[level + 1][xRangeLen * k + j]);
      // printf("\n");
    }
    for (int l = 0; l < 32; l++) {
      printf("%6d, ", ped[l]);
      if (l % 16 == 15)
        printf("\n");
    }
    // for (int l = 0; l < 4; l++) {
    //   for (int s = 0; s < 32; s++) {
    //     printf("%6d, ", py_extend[l][s]);
    //     if (s % 16 == 15)
    //       printf("\n");
    //   }
    // }
    //     printf("\n");
    for (int l = 0; l < 4; l++) {
      for (int k = 0; k < 32; k++) {
        printf("%2d, ", xHat[l][k]);
      }
      printf("\n");
    }
    level++;
    dim--;
    dimp--;
  }

  return 0;
}
