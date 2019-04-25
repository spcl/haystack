int main() {
  int N1, N2;
  float A[N2][N2];
  float B[N2][N2];

#pragma scop
  for (int t = 0; t < N1; t++) {
    for (int i = 1; i < N2 - 1; i++)
      for (int j = 1; j < N2 - 1; j++)
S0:     B[i][j] = 0.2 * (A[i][j] + A[i][j - 1] + A[i][1 + j] + A[1 + i][j] + A[i - 1][j]);
    for (int i = 1; i < N2 - 1; i++)
      for (int j = 1; j < N2 - 1; j++)
S1:     A[i][j] = 0.2 * (B[i][j] + B[i][j - 1] + B[i][1 + j] + B[1 + i][j] + B[i - 1][j]);
  }
#pragma endscop
}