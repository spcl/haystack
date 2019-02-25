int main() {
  float A[128][128];
  float B[128][128];
  float C[128][128];
  float alpha, beta, tmp;

  #define D 1024
#pragma scop
  for (int i = 0; i < D; i++) {
    for (int j = 0; j < D; j++)
      C[i][j] *= beta;
    for (int k = 0; k < D; k++) 
      for (int j = 0; j < D; j++)
        C[i][j] += alpha * A[i][k] * B[k][j];
  }
#pragma endscop
}

