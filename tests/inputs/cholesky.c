int main() {
  int N;
  float A[N][N];

#pragma scop
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < i; j++) {
      for (int k = 0; k < j; k++) {
S0:     A[i][j] -= A[i][k] * A[j][k];
      }
S1:   A[i][j] /= A[j][j]; 
    }
    for (int k = 0; k < i; k++) {
S2:   A[i][i] -= A[i][k] * A[i][k];
    }
S3: A[i][i] = A[i][i];
  }
#pragma endscop
}
