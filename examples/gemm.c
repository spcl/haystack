int main() {
  #define D 1024

  float A[D][D];
  float B[D][D];
  float C[D][D];
  float alpha, beta, tmp;
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

