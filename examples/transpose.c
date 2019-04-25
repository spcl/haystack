int main() {
  #define D 1024

  float A[D][D];
  float B[D][D];
#pragma scop
  for (int i = 0; i < D; i++) 
    for (int j = 0; j < D; j++) {
      B[i][j] = A[j][i];
  }
#pragma endscop
}

