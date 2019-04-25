int main() {
  #define I 64
  #define J 16
  #define K 64
  
  float tmp1, tmp2;
  
  float A[K][J][I];
  float B[K][J][I];
  float C[K][J][I];
  float D[K][J][I];
  float E[K][J][I];
  float F[K][J][I];
  float G[K][J][I];
  float H[K][J][I];
  float L[K][J][I];
  float M[K][J][I];
  float N[K][J][I];
  float O[K][J][I];
  float P[K][J][I];
  float Q[K][J][I];
  float R[K][J][I];
  float S[K][J][I];
  float T[K][J][I];
  float U[K][J][I];
  float V[K][J][I];

#pragma scop
  #define TK1 1
  #define TJ1 4

  for(int k=0; k<K/TK1; k++) 
    for(int j=0; j<J/TJ1; j++) { 
      for(int kk=max(1, k*TK1); kk<min(k*TK1+TK1, K-1)+1; kk++) 
        for(int jj=max(3, j*TJ1); jj<min(j*TJ1+TJ1, J-3)+1; jj++) 
          for(int i=3; i<I-3+1; i++) 
            G[kk][jj][i] = A[kk][jj][i] * F[kk][jj][i] + (1.0 - A[kk][jj][i]) * F[kk-1][jj][i];

      for(int kk=max(1, k*TK1); kk<min(k*TK1+TK1, K-1); kk++) 
        for(int jj=max(3, j*TJ1); jj<min(j*TJ1+TJ1, J-3)+1; jj++) 
          for(int i=3; i<I-3+1; i++) 
            H[kk][jj][i] = G[kk+1][jj][i] - G[kk][jj][i];

      for(int kk=max(1, k*TK1); kk<min(k*TK1+TK1, K-1); kk++) 
        for(int jj=max(3, j*TJ1); jj<min(j*TJ1+TJ1, J-3); jj++) 
          for(int i=3; i<I-3; i++)
            L[kk][jj][i] = (F[kk][jj][i+1] - F[kk][jj][i]) + (H[kk][jj][i+1] + H[kk][jj][i]) * 0.5 * 
                ((E[kk+1][jj][i] + E[kk][jj][i]) - (E[kk+1][jj][i+1] + E[kk][jj][i+1])) /
                ((E[kk+1][jj][i] - E[kk][jj][i]) + (E[kk+1][jj][i+1] - E[kk][jj][i+1]));

      for(int kk=max(1, k*TK1); kk<min(k*TK1+TK1, K-1); kk++) 
        for(int jj=max(3, j*TJ1); jj<min(j*TJ1+TJ1, J-3); jj++) 
          for(int i=3; i<I-3; i++) 
            M[kk][jj][i] = (F[kk][jj+1][i] - F[kk][jj][i]) + (H[kk][jj+1][i] + H[kk][jj][i]) * 0.5 * 
                ((E[kk+1][jj][i] + E[kk][jj][i]) - (E[kk+1][jj+1][i] + E[kk][jj+1][i])) /
                ((E[kk+1][jj][i] - E[kk][jj][i]) + (E[kk+1][jj+1][i] - E[kk][jj+1][i]));

      for(int kk=max(1, k*TK1); kk<min(k*TK1+TK1, K-1); kk++) 
        for(int jj=max(3, j*TJ1); jj<min(j*TJ1+TJ1, J-3); jj++)
          for(int i=3; i<I-3; i++) 
            N[kk][jj][i] =  P[kk][jj][i] + 0.01 * (O[kk][jj][i] - L[kk][jj][i] *
                (2.0 / (B[kk][jj][i+1] + B[kk][jj][i])));

      for(int kk=max(1, k*TK1); kk<min(k*TK1+TK1, K-1); kk++) 
        for(int jj=max(3, j*TJ1); jj<min(j*TJ1+TJ1, J-3); jj++) 
          for(int i=3; i<I-3; i++) 
            Q[kk][jj][i] =  S[kk][jj][i] + 0.01 * (R[kk][jj][i] - M[kk][jj][i] *
                (2.0 / (B[kk][jj+1][i] + B[kk][jj][i])));
    }

  #define TK2 4
  #define TJ2 4

  for(int k=0; k<K/TK2; k++) 
    for(int j=0; j<J/TJ2; j++) { 
      for(int kk=max(1, k*TK2); kk<min(k*TK2+TK2, K-1)+1; kk++) 
        for(int jj=max(3, j*TJ2); jj<min(j*TJ2+TJ2, J-3); jj++) 
          for(int i=3-1; i<I-3; i++) {
            tmp1 = 0.5 * (A[kk][jj][i] + A[kk][jj][i+1]);
            T[kk][jj][i] = tmp1 * N[kk][jj][i] + (1.0 - tmp1) * N[kk-1][jj][i];
          }

      for(int kk=max(1, k*TK2); kk<min(k*TK2+TK2, K-1)+1; kk++) 
        for(int jj=max(3, j*TJ2)-1; jj<min(j*TJ2+TJ2, J-3); jj++) 
          for(int i=3; i<I-3; i++) {
            tmp2 = 0.5 * (A[kk][jj][i] + A[kk][jj+1][i]);
            U[kk][jj][i] = tmp2 * Q[kk][jj][i] + (1.0 - tmp2) * Q[kk-1][jj][i];
          }

      for(int kk=max(1, k*TK2); kk<min(k*TK2+TK2, K-1); kk++) 
        for(int jj=max(3, j*TJ2); jj<min(j*TJ2+TJ2, J-3); jj++) 
          for(int i=3; i<I-3; i++) 
            V[kk][jj][i] =
                0.1 * ((T[kk+1][jj][i] - T[kk][jj][i]) * C[kk][jj][i] + N[kk][jj][i]) + 
                0.1 * ((T[kk+1][jj][i-1] - T[kk][jj][i-1]) * C[kk][jj][i-1] - N[kk][jj][i-1]) + 
                0.2 * ((U[kk+1][jj][i] - U[kk][jj][i]) * D[kk][jj][i] + Q[kk][jj][i]) + 
                0.3 * ((U[kk+1][jj-1][i] - U[kk][jj-1][i]) * D[kk][jj-1][i] - Q[kk][jj-1][i]);
  }
#pragma endscop
}