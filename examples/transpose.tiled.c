int main() {
  #define T 256
  #define D 1024

  float A[D][D];
  float B[D][D];
#pragma scop
  for (int i = 0; i < D/T; i++) 
    for (int j = 0; j < D/T; j++) 
      for (int ii = i*T; ii < i*T + T; ii++)
        for (int jj = j*T; jj < j*T + T; jj++) {
          B[ii][jj] = A[jj][ii];
  }
#pragma endscop
}

