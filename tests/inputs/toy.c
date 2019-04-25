int main() {
  int N;
  float A[N];

#pragma scop
  for (int i = 0; i < N; i++) {
S0: A[i] = 0;
S1: A[N-1-i] = 1;
    if (i < N/2)
S2:   A[2*i] = 2;
  }
#pragma endscop
}
