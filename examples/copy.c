int main() {
  int i, j, k;
  int N1, N2, N3;
  float A[N1][N2][N3];
  float B[N1][N2][N3];
  float T[N3][N2][N1];

#pragma scop
  for (i = 0; i < N1; i++)
    for (j = 0; j < N2; j++)
      for (k = 0; k < N3; k++)
S0:     T[k][j][i] = A[i][j][k];
  for (i = 0; i < N1; i++)
    for (j = 0; j < N2; j++)
      for (k = 0; k < N3; k++)
S1:     B[i][j][k] = T[k][j][i];
#pragma endscop
}
