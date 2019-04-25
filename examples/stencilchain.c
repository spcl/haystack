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
  for(int k=1; k<K-1; k++) 
    for(int j=3; j<J-3; j++) 
      for(int i=3; i<I-3; i++) 
        G[k][j][i] = A[k][j][i] * F[k][j][i] + (1.0 - A[k][j][i]) * F[k-1][j][i];

  for(int k=1; k<K-1; k++) 
    for(int j=3; j<J-3; j++) 
      for(int i=3; i<I-3; i++) 
        H[k][j][i] = G[k+1][j][i] - G[k][j][i];

  for(int k=1; k<K-1; k++) 
    for(int j=3; j<J-3; j++) 
      for(int i=3; i<I-3; i++)
        L[k][j][i] = (F[k][j][i+1] - F[k][j][i]) + (H[k][j][i+1] + H[k][j][i]) * 0.5 * 
            ((E[k+1][j][i] + E[k][j][i]) - (E[k+1][j][i+1] + E[k][j][i+1])) /
            ((E[k+1][j][i] - E[k][j][i]) + (E[k+1][j][i+1] - E[k][j][i+1]));

  for(int k=1; k<K-1; k++) 
    for(int j=3; j<J-3; j++) 
      for(int i=3; i<I-3; i++) 
        M[k][j][i] = (F[k][j+1][i] - F[k][j][i]) + (H[k][j+1][i] + H[k][j][i]) * 0.5 * 
            ((E[k+1][j][i] + E[k][j][i]) - (E[k+1][j+1][i] + E[k][j+1][i])) /
            ((E[k+1][j][i] - E[k][j][i]) + (E[k+1][j+1][i] - E[k][j+1][i]));

  for(int k=1; k<K-1; k++) 
    for(int j=3; j<J-3; j++) 
      for(int i=3; i<I-3; i++) 
        N[k][j][i] =  P[k][j][i] + 0.01 * (O[k][j][i] - L[k][j][i] *
            (2.0 / (B[k][j][i+1] + B[k][j][i])));

  for(int k=1; k<K-1; k++) 
    for(int j=3; j<J-3; j++) 
      for(int i=3; i<I-3; i++) 
        Q[k][j][i] =  S[k][j][i] + 0.01 * (R[k][j][i] - M[k][j][i] *
            (2.0 / (B[k][j+1][i] + B[k][j][i])));

  for(int k=1; k<K-1; k++) 
    for(int j=3; j<J-3; j++) 
      for(int i=3; i<I-3; i++) {
        tmp1 = 0.5 * (A[k][j][i] + A[k][j][i+1]);
        T[k][j][i] = tmp1 * N[k][j][i] + (1.0 - tmp1) * N[k-1][j][i];
      }

  for(int k=1; k<K-1; k++) 
    for(int j=3; j<J-3; j++) 
      for(int i=3; i<I-3; i++) {
        tmp2 = 0.5 * (A[k][j][i] + A[k][j+1][i]);
        U[k][j][i] = tmp2 * Q[k][j][i] + (1.0 - tmp2) * Q[k-1][j][i];
      }

  for(int k=1; k<K-1; k++) 
    for(int j=3; j<J-3; j++) 
      for(int i=3; i<I-3; i++) 
        V[k][j][i] =
            0.1 * ((T[k+1][j][i] - T[k][j][i]) * C[k][j][i] + N[k][j][i]) + 
            0.1 * ((T[k+1][j][i-1] - T[k][j][i-1]) * C[k][j][i-1] - N[k][j][i-1]) + 
            0.2 * ((U[k+1][j][i] - U[k][j][i]) * D[k][j][i] + Q[k][j][i]) + 
            0.3 * ((U[k+1][j-1][i] - U[k][j-1][i]) * D[k][j-1][i] - Q[k][j-1][i]);
#pragma endscop
}