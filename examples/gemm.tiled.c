int main() {
  float A[128][128];
  float B[128][128];
  float C[128][128];
  float beta, alpha;
  
  // TODO test simpler formulation?

  #define T 64
  #define D 1024
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
  // for (int i = 0; i < D; i+=T) {
  //   for (int j = 0; j < D; j+=T) {
  //     for (int ii = i; ii < i + T; ii++)
  //       for (int jj = j; jj < j + T; jj++)
  //         C[ii][jj] *= beta;
  //     for (int k = 0; k < D; k+=T) 
  //       for (int ii = i; ii < i + T; ii++)
  //         for (int kk = k; kk < k + T; kk++)
  //           for (int jj = j; jj < j + T; jj++)
  //             C[ii][jj] += alpha * A[ii][kk] * B[kk][jj];
  //   }
  // }
#pragma endscop
}
