#include <stdio.h>
#include <stdint.h>
#include <string.h>

#define SQRT2  362
#define SQRT10 810
#define SQRT42 1659

typedef struct c16_t{
  int16_t r;
  int16_t i;
} c16_t;

int8_t xRange[3][8] = {
  {-1,  1,  0,  0, 0, 0, 0, 0},
  {-3, -1,  1,  3, 0, 0, 0, 0},
  {-7, -5, -3, -1, 1, 3, 5, 7},
};

int16_t QPSK[4] = {3, 2, 1, 0};

int16_t QAM16[16] = {15, 14, 10, 11,
                     13, 12,  8,  9,
                      5,  4,  0,  1,
                      7,  6,  2,  3};

int16_t QAM64[64] = {63, 62, 58, 59, 43, 42, 46, 47,
                     61, 60, 56, 57, 41, 40, 44, 45,
                     53, 52, 48, 49, 33, 32, 36, 37,
                     55, 54, 50, 51, 35, 34, 38, 39,
                     23, 22, 18, 19,  3,  2,  6,  7,
                     21, 20, 16, 17,  1,  0,  4,  5,
                     29, 28, 24, 25,  9,  8, 12, 13,
                     31, 30, 26, 27, 13, 12, 14, 15};

void print_matrix_int16 (int m, int n, int16_t A[m][n]) {
  for (int8_t i = 0; i < m; i++) {
    printf("[ ");
    for (int8_t j = 0; j < n; j++)
      printf("%6d, ", A[i][j]);
    printf("];\n");
  }
}

int32_t mul_scale (int32_t a, int32_t b) {
  return (a * b) >> 8;
}

int32_t div_scale(int32_t a, int32_t b) {
  return (a << 8) / b;
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

void QR_fixed (int m, int n, int16_t A_scale[m][n], int16_t Q_scale[m][m], int16_t R_scale[n][n]) {
  int32_t w[m];
  int32_t Q_t[m];
  int32_t w_norm;
  for (int8_t j = 0; j < n; j++) {
    for (int8_t k = 0; k < m; k++)
      w[k] = A_scale[k][j];
    for (int8_t i = 0; i < j; i++) {
      transpose(n, Q_scale, i, Q_t);
      R_scale[i][j] = product(n, Q_t, w); // Hq is << 8, product will >> 8, so Hr will be correct
      for (int8_t k = 0; k < m; k++)
        w[k] -= mul_scale(R_scale[i][j], (int32_t)Q_scale[k][i]);
    }
    w_norm = norm_scale(m, w);
    for (int8_t k = 0; k < m; k++)
      Q_scale[k][j] = div_scale(w[k], w_norm); // Hq will << 8
    R_scale[j][j] = w_norm;
  }
}

int16_t raduis_norm (c16_t a, c16_t b) {
  int16_t sorted_idx;
  int32_t a32 = 0, b32 = 0;
  a32 = a.r * a.r + a.i * a.i;
  b32 = b.r * b.r + b.i * b.i;
  sorted_idx = my_sqrt(a32 + b32);
  return sorted_idx;
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

int binary_search(int64_t *a, int8_t *sorted_idx, int l, int r, int target) {
  if (r <= l)
    return (a[target] > a[sorted_idx[l]] ? l + 1 : l);
  int mid = (l + r) / 2;
  if (a[target] == a[sorted_idx[mid]])
    return mid + 1;
  if (a[target] > a[sorted_idx[mid]])
    return binary_search(a, sorted_idx, mid + 1, r, target);
  return binary_search(a, sorted_idx, l, mid - 1, target);
}

void sort_ped (int k_size, int64_t *a, int8_t *sorted_idx) {
  int j = 0;
  sorted_idx[0] = 0;
  for (int k = 1; k < k_size; k++) {
    j = k - 1;
    int selected = k;
    int loc = binary_search(a, sorted_idx, 0, j, selected);
    while (j >= loc) {
      sorted_idx[j + 1] = sorted_idx[j];
      j--;
    }
    sorted_idx[j + 1] = selected;
  }
}

int map_iq_to_qam_idx (int8_t xRangeLen, int8_t *range, int amplitude, int phase) {
  int amplitude_index = 0, phase_index = 0;
  
  for (amplitude_index = 0; amplitude_index < xRangeLen; amplitude_index++) {
    if (range[amplitude_index] == amplitude) {
      break;
    }
  }
  
  for (phase_index = 0; phase_index < xRangeLen; phase_index++) {
    if (range[phase_index] == phase) {
      break;
    }
  }
  return amplitude_index * xRangeLen + phase_index;
}

void my_SDA_kbest (int numLayer,
                   int numAnt,
                   c16_t *ch,
                   c16_t *data_in,
                   int order,
                   int k_best) 
{
  for (int i = 0; i < 2 && order == 6; i++) {
    data_in[i].r /= 12;
    data_in[i].i /= 12;
  }
  c16_t ch0 = ch[0];
  c16_t ch1 = ch[1];
  c16_t ch2 = ch[2];
  c16_t ch3 = ch[3];
  // ZF H^-1 * y
  // c16_t y_zf[2] = {0};
  // c16_t tmp_ch = {0};
  // c16_t zf0 = cmplx_mul(ch3, data_in[0]);
  // tmp_ch.r = -ch1.r;
  // tmp_ch.i = -ch1.i;
  // c16_t zf1 = cmplx_mul(tmp_ch, data_in[1]);
  // tmp_ch.r = -ch2.r;
  // tmp_ch.i = -ch2.i;
  // c16_t zf2 = cmplx_mul(tmp_ch, data_in[0]);
  // c16_t zf3 = cmplx_mul(ch0, data_in[1]);
  // y_zf[0] = cmplx_add(zf0, zf1);
  // y_zf[1] = cmplx_add(zf2, zf3);

  // c16_t data_yzf_0 = cmplx_sub(data_in[0], y_zf[0]);
  // c16_t data_yzf_1 = cmplx_sub(data_in[1], y_zf[1]);

  switch (order) {
    case 2:
      // data_yzf_0 = cmplx_scale(data_yzf_0, SQRT2);
      // data_yzf_1 = cmplx_scale(data_yzf_1, SQRT2);
      data_in[0] = cmplx_scale(data_in[0], SQRT2);
      data_in[1] = cmplx_scale(data_in[1], SQRT2);
      break;
    case 4:
      // data_yzf_0 = cmplx_scale(data_yzf_0, SQRT10);
      // data_yzf_1 = cmplx_scale(data_yzf_1, SQRT10);
      data_in[0] = cmplx_scale(data_in[0], SQRT10);
      data_in[1] = cmplx_scale(data_in[1], SQRT10);
      break;
    case 6:
      // data_yzf_0 = cmplx_scale(data_yzf_0, SQRT42);
      // data_yzf_1 = cmplx_scale(data_yzf_1, SQRT42);
      data_in[0] = cmplx_scale(data_in[0], SQRT42);
      data_in[1] = cmplx_scale(data_in[1], SQRT42);
      break;
    default:
      printf("Modulation order %d is not support.\n", order);
  }
  // int32_t raduis = raduis_norm(data_yzf_0, data_yzf_1);

  int16_t H[4][4] = {{ch0.r, ch1.r, -ch0.i, -ch1.i}, 
                     {ch2.r, ch3.r, -ch2.i, -ch3.i}, 
                     {ch0.i, ch1.i,  ch0.r,  ch1.r}, 
                     {ch2.i, ch3.i,  ch2.r,  ch3.r}};
  int16_t Hq[4][4] = {0};
  int16_t Hr[4][4] = {0};

  QR_fixed(numLayer, numAnt, H, Hq, Hr);

  int xRangeLen = 2 * xRange[order/2 -1 ][0] / -2 + 1;
  int k_size = xRangeLen;

  int16_t level = 0;
  int16_t xHat[4][64] = {0};
  // int32_t dp[64] = {0};      // radius
  int64_t ped[64] = {0};     // add for DP
  int16_t py[4][64] = {0};   // Store Hq * Y -> yR

  for (int i = 0; i < 8; i++)  // init xHat
    xHat[level][i] = xRange[order / 2 - 1][i];
 
  int16_t data_in16[4] = {data_in[0].r, data_in[1].r, data_in[0].i, data_in[1].i};
  int16_t yR[4] = {0}; 

  for (int i = 0; i < 4; i++) { // Hq^-1 * Y
    int32_t tmp = 0;
    for (int j = 0; j < 4; j++) {
      tmp += (int32_t)Hq[j][i] * (int32_t)data_in16[j] / 256;
    }
    yR[i] = (int16_t)tmp;
  }
  
  // for (int i = 0; i < k_size; i++) // init dp
  //   dp[i] = raduis;

  for (int i = 0; i < numLayer; i++) // init py
    for (int j = 0; j < k_size; j++)
      py[i][j] = yR[i];
  
  int dimp = numLayer - 1, dim = numLayer - 2;

  for (int layers = 0; layers < numLayer; layers++) {
    int k_best_cur = (k_size > k_best) ? k_best : k_size;

    int64_t tmp_w[64] = {INT64_MAX};
    for (int k = 0; k < k_size; k++) { // (py - Hr * xHat)^2
      tmp_w[k] = ((int64_t)py[dimp][k] - (int64_t)Hr[dimp][dimp] * xHat[level][k]); 
      if (tmp_w[k] < 0)
        tmp_w[k] *= -1;
      ped[k] += tmp_w[k];
    }

    // int64_t tmp_dp[64] = {INT64_MAX};
    // for (int k = 0; k < k_size; k++) { // ped = radius - tmp_w
    //   tmp_dp[k] = (int64_t)dp[k] - tmp_w[k]; 
    //   if (tmp_dp[k] < 0)
    //     tmp_dp[k] = 0;
    //   ped[k] += tmp_w[k];
    // }

    int8_t sorted_idx[64];
    sort_ped(k_size, ped, sorted_idx); // sort for finding kbest

    if (layers == numLayer - 1) { // handle the anser
      printf("Ans: (%d, %dj), (%d, %dj)\n", xHat[3][sorted_idx[0]], xHat[1][sorted_idx[0]], xHat[2][sorted_idx[0]], xHat[0][sorted_idx[0]]);
      
      // int l0_idx = map_iq_to_qam_idx(xRangeLen, xRange[order / 2 - 1], xHat[3][sorted_idx[0]], xHat[1][sorted_idx[0]]);
      // int l1_idx = map_iq_to_qam_idx(xRangeLen, xRange[order / 2 - 1], xHat[2][sorted_idx[0]], xHat[0][sorted_idx[0]]);
      // switch (order) {
      //   case 2:
      //     l0_idx = QPSK[l0_idx];
      //     l1_idx = QPSK[l1_idx];
      //     break;
      //   case 4:
      //     l0_idx = QAM16[l0_idx];
      //     l1_idx = QAM16[l1_idx];
      //     break;
      //   case 6:
      //     l0_idx = QAM64[l0_idx];
      //     l1_idx = QAM64[l1_idx];
      //     break;
      //   default:
      //     break;
      // }
      // break;
    }
    
    int16_t xHat_tmp[4][64] = {0};
    int64_t ped_tmp[64] = {0};
    int16_t py_extend_tmp[4][64];
    memcpy(xHat_tmp, xHat, sizeof(xHat));
    memcpy(ped_tmp, ped, sizeof(ped));
    memcpy(py_extend_tmp, py, sizeof(py));

    for (int k = 0; k < k_best_cur; k++) {
      for (int j = 0; j <= dim; j++) {
        py_extend_tmp[j][sorted_idx[k]] = py_extend_tmp[j][sorted_idx[k]] - Hr[j][dimp] * xHat_tmp[level][sorted_idx[k]];
      }
      for (int j = 0; j < xRangeLen; j++) {
        for (int l = 0; l < 4; l++) {
          py[l][xRangeLen * k + j] = py_extend_tmp[l][sorted_idx[k]];
        }
        // set add the new xHat
        xHat[level][xRangeLen * k + j] = xHat_tmp[level][sorted_idx[k]];
        xHat[level + 1][xRangeLen * k + j] = xRange[order / 2 - 1][j];
        // set the new ped
        ped[xRangeLen * k + j] = ped_tmp[sorted_idx[k]];
        // set the new dp (radius)
        // dp[xRangeLen * k + j] = tmp_dp[sorted_idx[k]]; 
      }
    }
    k_size = k_best_cur * xRangeLen; // set the new k_size
    level++;
    dim--;
    dimp--;
  }
}

int main() {

  int numLayer = 2;
  int numAnt = 2;
  int k_best = 8;
  int8_t order = 4;
  
  c16_t ch[] = {
    {334, -30},
    {40, -6},
    {-6, -4},
    {378, 0},
  };
  c16_t data_in[2] = {{591, -544}, {-223, 93}};

  my_SDA_kbest(2 * numLayer, 2 * numAnt, ch, data_in, order, k_best);
  
  return 0;
}
