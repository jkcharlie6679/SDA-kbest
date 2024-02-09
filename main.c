#include <stdio.h>
#include <stdint.h>
#include <string.h>

// #define DEBUG

#define SQRT2  362
#define SQRT10 810
#define SQRT42 1659

typedef struct c16_t{
  int16_t r;
  int16_t i;
} c16_t;

typedef struct c32_t{
  int32_t r;
  int32_t i;
} c32_t;

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

int64_t my_sqrt(uint64_t data) {
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
    if (j == n - 1)
      w_norm = -w_norm;
    for (int8_t k = 0; k < m; k++)
      Q_scale[k][j] = div_scale(w[k], w_norm); // Hq will << 8
    R_scale[j][j] = w_norm;
  }
}

int64_t radius_norm (c32_t a, c32_t b) {
  int64_t sorted_idx;
  uint64_t a64 = 0, b64 = 0;
  a64 = (int64_t)a.r * a.r + (int64_t)a.i * a.i;
  b64 = (int64_t)b.r * b.r + (int64_t)b.i * b.i;
  uint64_t tmp = a64 + b64;
  sorted_idx = my_sqrt(tmp);
  return sorted_idx;
}

c32_t cmplx_mul (c16_t a, c16_t b) {
  c32_t result = {0, 0};
  result.r = ((int32_t)a.r * (int32_t)b.r - (int32_t)a.i * (int32_t)b.i);
  result.i = ((int32_t)a.r * (int32_t)b.i + (int32_t)a.i * (int32_t)b.r);
  return result;
}

c32_t cmplx_add (c32_t a, c32_t b) {
  c32_t result = {0};
  result.r = a.r + b.r;
  result.i = a.i + b.i;
  return  result;
}
c32_t cmplx_sub (c16_t a, c32_t b) {
  c32_t result = {0};
  result.r = (int32_t)a.r - b.r;
  result.i = (int32_t)a.i - b.i;
  return  result;
}

c32_t cmplx32_scale (c32_t a, int b) {
  if (b >= 0) {
    a.r = (((int32_t)a.r * b) >> 8);
    a.i = (((int32_t)a.i * b) >> 8);
  } 
  return a;
}

c16_t cmplx16_scale (c16_t a, int b) {
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

void my_SDA_kbest (int numLayer,
                   int numAnt,
                   c16_t ch[numLayer][numAnt],
                   c16_t *data_in,
                   int order,
                   int k_best) 
{
  // ZF H^-1 * y
  c32_t y_zf[2] = {0};
  c16_t tmp_ch = {0};
  c32_t zf0 = cmplx_mul(ch[1][1], data_in[0]);
  tmp_ch.r = -ch[0][1].r;
  tmp_ch.i = -ch[0][1].i;
  c32_t zf1 = cmplx_mul(tmp_ch, data_in[1]);
  tmp_ch.r = -ch[1][0].r;
  tmp_ch.i = -ch[1][0].i;
  c32_t zf2 = cmplx_mul(tmp_ch, data_in[0]);
  c32_t zf3 = cmplx_mul(ch[0][0], data_in[1]);
  y_zf[0] = cmplx_add(zf0, zf1);
  y_zf[1] = cmplx_add(zf2, zf3);

  c32_t data_yzf_0 = cmplx_sub(data_in[0], y_zf[0]);
  c32_t data_yzf_1 = cmplx_sub(data_in[1], y_zf[1]);
  int64_t det = (int64_t)ch[0][0].r * ch[1][1].r - (int64_t)ch[0][0].i * ch[1][1].i - (int64_t)ch[0][1].r * ch[1][0].r + (int64_t)ch[0][1].i * ch[1][0].i;

  switch (order) {
    case 2:
      data_yzf_0 = cmplx32_scale(data_yzf_0, SQRT2);
      data_yzf_1 = cmplx32_scale(data_yzf_1, SQRT2);
      data_in[0] = cmplx16_scale(data_in[0], SQRT2);
      data_in[1] = cmplx16_scale(data_in[1], SQRT2);
      break;
    case 4:
      data_yzf_0 = cmplx32_scale(data_yzf_0, SQRT10);
      data_yzf_1 = cmplx32_scale(data_yzf_1, SQRT10);
      data_in[0] = cmplx16_scale(data_in[0], SQRT10);
      data_in[1] = cmplx16_scale(data_in[1], SQRT10);
      break;
    case 6:
      data_yzf_0 = cmplx32_scale(data_yzf_0, SQRT42);
      data_yzf_1 = cmplx32_scale(data_yzf_1, SQRT42);
      data_in[0] = cmplx16_scale(data_in[0], SQRT42);
      data_in[1] = cmplx16_scale(data_in[1], SQRT42);
      break;
    default:
      printf("Modulation order %d is not support.\n", order);
  }

  int64_t radius = radius_norm(data_yzf_0, data_yzf_1);
  radius /= my_sqrt(det);
  radius *= radius;

  int16_t layer_ext = numLayer << 1;
  int16_t ant_ext = numAnt << 1;

  int16_t H[layer_ext][ant_ext];
  for (int l = 0; l < layer_ext; l++) {
    for (int a = 0; a < ant_ext; a++) {
      if (l < numLayer && a < numAnt) {
        H[l][a] = ch[l % numLayer][a % numAnt].r;
      } else if (l < numLayer&& a >= numAnt) {
        H[l][a] = -ch[l % numLayer][a % numAnt].i;
      } else if (l >= numLayer && a < numAnt) {
        H[l][a] = ch[l % numLayer][a % numAnt].i;
      } else {
        H[l][a] = ch[l % numLayer][a % numAnt].r;
      }
    }
  }

  int16_t Hq[layer_ext][ant_ext];
  int16_t Hr[layer_ext][ant_ext];
  memset(Hq, 0, sizeof(Hq));
  memset(Hr, 0, sizeof(Hr));

  QR_fixed(layer_ext, ant_ext, H, Hq, Hr);

#ifdef DEBUG
  printf("Hq = \n");
  print_matrix_int16(layer_ext, ant_ext, Hq);
  printf("\n");
  printf("Hr = \n");
  print_matrix_int16(layer_ext, ant_ext, Hr);
  printf("\n");
#endif
  
  int xRangeLen = 2 * xRange[order/2 -1 ][0] / -2 + 1;
  int k_size = xRangeLen;

  int16_t level = 0;
  int16_t xHat[layer_ext][k_best * xRangeLen];
  int16_t data_in16[layer_ext];
  int64_t dp[k_best * xRangeLen];             // radius
  int64_t py[layer_ext][k_best * xRangeLen];   // Store Hq * Y -> yR
  int8_t sorted_idx[k_best * xRangeLen];
  
  memset(xHat, 0, sizeof(xHat));
  memset(dp, 0, sizeof(dp));
  memset(py, 0, sizeof(py));

  for (int i = 0; i < xRangeLen; i++)  // init xHat
    xHat[level][i] = xRange[order / 2 - 1][i];

  for (int l = 0; l < layer_ext; l++) {
    if (l < numLayer)
      data_in16[l] = data_in[l % numLayer].r;
    else
      data_in16[l] = data_in[l % numLayer].i;
  }

  int64_t yR[4];
#ifdef DEBUG
  printf("yR = [");
#endif
  for (int l = 0; l < layer_ext; l++) { // Hq^-1 * Y
    int64_t tmp = 0;
    for (int a = 0; a < ant_ext; a++) {
      tmp += (int64_t)Hq[a][l] * (int64_t)data_in16[a] / 256;
    }
    yR[l] = tmp;
#ifdef DEBUG
    printf("%lld, ",  yR[l]);
#endif
  }
#ifdef DEBUG
  printf("]\n");
#endif
  
  for (int i = 0; i < k_size; i++) // init dp
    dp[i] = radius;

  for (int i = 0; i < layer_ext; i++) // init py
    for (int j = 0; j < k_size; j++)
      py[i][j] = yR[i];
  
  int dimp = layer_ext - 1, dim = layer_ext - 2;

  for (int layers = 0; layers < layer_ext; layers++) {
    int k_best_cur = (k_size > k_best) ? k_best : k_size;
#ifdef DEBUG
    printf("------ l = %d -------\n", layers);
    printf("xHat = [\n");
    for (int i = 0; i <= level; i++) {
      for (int j = 0; j < k_size; j++)
        printf("%2d, ", xHat[i][j]);
      printf("\n");
    }
    printf("]\n");
#endif
    int64_t tmp_w[k_best * xRangeLen];
    memset(tmp_w, INT64_MAX, sizeof(tmp_w));
    for (int k = 0; k < k_size; k++) { // (py - Hr * xHat)^2
      tmp_w[k] = ((int64_t)py[dimp][k] - (int64_t)Hr[dimp][dimp] * xHat[level][k]); 
      tmp_w[k] *= tmp_w[k];
    }
#ifdef DEBUG
    printf("tmp_w = [\n");
    for (int j = 0; j < k_size; j++)
      printf("%lld, ", tmp_w[j]);
    printf("]\n");
    printf("dp_B = [\n");
    for (int j = 0; j < k_size; j++)
      printf("%lld, ", dp[j]);
    printf("]\n");
#endif
    int eff_count = k_size;
    for (int k = 0; k < k_size; k++) { // radius - tmp_w
      dp[k] -= tmp_w[k]; 
      if (dp[k] < 0) {
        eff_count--;
        dp[k] = 0;
      }
    }
    if (eff_count < k_best_cur)
      k_best_cur = eff_count;

#ifdef DEBUG
    printf("dp_A = [\n");
    for (int j = 0; j < k_size; j++)
      printf("%lld, ", dp[j]);
    printf("]\n");
#endif

    sort_ped(k_size, dp, sorted_idx); // sort for finding kbest

    if (layers == layer_ext - 1) { // handle the anser
      printf("Ans: (%d, %dj), (%d, %dj)\n", xHat[3][sorted_idx[k_size - 1]], xHat[1][sorted_idx[k_size - 1]], xHat[2][sorted_idx[k_size - 1]], xHat[0][sorted_idx[k_size - 1]]);
      break;
    }
    
    int16_t xHat_tmp[layer_ext][k_best * xRangeLen];
    int64_t dp_tmp[k_best * xRangeLen];
    int64_t py_tmp[layer_ext][k_best * xRangeLen];
    memcpy(xHat_tmp, xHat, sizeof(xHat));
    memcpy(dp_tmp, dp, sizeof(dp));
    memcpy(py_tmp, py, sizeof(py));
    for (int k = 0; k < k_best_cur; k++) {
      int sel_idx = sorted_idx[k_size - k - 1];
      for (int l = 0; l <= dim; l++)
        py_tmp[l][sel_idx] = py_tmp[l][sel_idx] - Hr[l][dimp] * xHat_tmp[level][sel_idx];

      for (int j = 0; j < xRangeLen; j++) {
        for (int l = 0; l < layer_ext; l++)
          py[l][xRangeLen * k + j] = py_tmp[l][sel_idx];

        // set add the new xHat
        for (int i = 0; i < level + 1; i++) {
          xHat[i][xRangeLen * k + j] = xHat_tmp[i][sel_idx];
        }
        xHat[level + 1][xRangeLen * k + j] = xRange[order / 2 - 1][j];

        // set the new dp (radius)
        dp[xRangeLen * k + j] = dp_tmp[sel_idx];
      }
    }
    k_size = k_best_cur * xRangeLen; // set the new k_size

#ifdef DEBUG
    printf("dp_C = [\n");
    for (int j = 0; j < k_size; j++)
      printf("%lld, ", dp[j]);
    printf("]\n");
#endif

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
  
  c16_t ch[2][2] = {
    // {{384, 46},{-58, 18}},
    // {{16, -48},{334, 38}},
  };
  c16_t data_in[2] = {{419, -652}, {-29, 178}};

  my_SDA_kbest(numLayer, numAnt, ch, data_in, order, k_best);
  
  return 0;
}
