int main() {
  int i, j, k;
  int N;
  float A[N][N];
  float sqrt;

#pragma scop
  for (i = 0; i < N; i++) {
    // j<i
    for (j = 0; j < i; j++) {
      for (k = 0; k < j; k++) {
S0:     A[i][j] -= A[i][k] * A[j][k];
      }
S1:   A[i][j] /= A[j][j]; 
    }
    // i==j case
    for (k = 0; k < i; k++) {
S2:   A[i][i] -= A[i][k] * A[i][k];
    }
S3: A[i][i] = sqrt * A[i][i];
  }
#pragma endscop
}
