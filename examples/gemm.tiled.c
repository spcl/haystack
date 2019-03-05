int main() {
  #define T 32
  #define D 1024

  float A[D][D];
  float B[D][D];
  float C[D][D];
  float beta, alpha;
  
#pragma scop
  for (int i = 0; i < D/T; i++) {
    for (int j = 0; j < D/T; j++) {
      for (int ii = i*T; ii < i*T + T; ii++)
        for (int jj = j*T; jj < j*T + T; jj++)
          C[ii][jj] *= beta;
      for (int k = 0; k < D/T; k++) 
        for (int ii = i*T; ii < i*T + T; ii++)
          for (int kk = k*T; kk < k*T + T; kk++)
            for (int jj = j*T; jj < j*T + T; jj++)
              C[ii][jj] += alpha * A[ii][kk] * B[kk][jj];
    }
  }
#pragma endscop
}
