int main() {
  int N1, N2, N3;
  float A[N1][N3];
  float B[N3][N2];
  float C[N1][N2];
  float alpha, beta, tmp;
  
#pragma scop
  for (int i = 0; i < N1; i++) {
    for (int j = 0; j < N2; j++)
S0:   C[i][j] *= beta;
    for (int k = 0; k < N3; k++) 
      for (int j = 0; j < N2; j++) {
S1:     C[i][j] += alpha * A[i][k] * B[k][j]; 
    }
  }
#pragma endscop
}