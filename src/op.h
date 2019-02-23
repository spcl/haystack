#include <isl/isl-noexceptions.h>

#ifndef ISL_OPERATORS
#define ISL_OPERATORS

inline isl::pw_aff operator+(isl::pw_aff A, isl::pw_aff B) { return A.add(B); }

inline isl::pw_aff operator+(isl::val V, isl::pw_aff A) {
  isl::pw_aff AV(A.domain(), V);
  return A.add(AV);
}

inline isl::pw_aff operator+(isl::pw_aff A, isl::val V) { return V + A; }

inline isl::pw_aff operator+(int i, isl::pw_aff A) {
  isl::ctx ctx = A.get_ctx();
  return A + isl::val(ctx, i);
}

inline isl::pw_aff operator+(isl::pw_aff A, int i) { return i + A; }

inline isl::val operator+(isl::val A, isl::val B) {
  return A.mul(B);
}

inline isl::pw_aff operator*(isl::pw_aff A, isl::pw_aff B) { return A.mul(B); }

inline isl::pw_aff operator*(isl::val V, isl::pw_aff A) {
  isl::pw_aff AV(A.domain(), V);
  return A.mul(AV);
}

inline isl::pw_aff operator*(isl::pw_aff A, isl::val V) { return V * A; }

inline isl::pw_aff operator*(int i, isl::pw_aff A) {
  isl::ctx ctx = A.get_ctx();
  return A * isl::val(ctx, i);
}

inline isl::pw_aff operator*(isl::pw_aff A, int i) { return i * A; }

inline isl::val operator*(isl::val A, isl::val B) {
  return A.mul(B);
}

inline isl::pw_aff operator-(isl::pw_aff A, isl::pw_aff B) { return A.sub(B); }

inline isl::pw_aff operator-(isl::val V, isl::pw_aff A) {
  isl::pw_aff AV(A.domain(), V);
  return AV - A;
}

inline isl::pw_aff operator-(isl::pw_aff A, isl::val V) {
  isl::pw_aff AV(A.domain(), V);
  return A - AV;
}

inline isl::pw_aff operator-(int i, isl::pw_aff A) {
  isl::ctx ctx = A.get_ctx();
  return isl::val(ctx, i) - A;
}

inline isl::pw_aff operator-(isl::pw_aff A, int i) {
  isl::ctx ctx = A.get_ctx();
  return A - isl::val(ctx, i);
}

inline isl::pw_aff operator/(isl::pw_aff A, isl::pw_aff B) {
  return A.tdiv_q(B);
}

inline isl::pw_aff operator/(isl::val V, isl::pw_aff A) {
  isl::pw_aff AV(A.domain(), V);
  return AV / A;
}

inline isl::pw_aff operator/(isl::pw_aff A, isl::val V) {
  isl::pw_aff AV(A.domain(), V);
  return A / AV;
}

inline isl::pw_aff operator/(int i, isl::pw_aff A) {
  isl::ctx ctx = A.get_ctx();
  return isl::val(ctx, i) / A;
}

inline isl::pw_aff operator/(isl::pw_aff A, int i) {
  isl::ctx ctx = A.get_ctx();
  return A / isl::val(ctx, i);
}

inline isl::pw_aff operator%(isl::pw_aff A, isl::pw_aff B) {
  return A.tdiv_r(B);
}

inline isl::pw_aff operator%(isl::val V, isl::pw_aff A) {
  isl::pw_aff AV(A.domain(), V);
  return AV % A;
}

inline isl::pw_aff operator%(isl::pw_aff A, isl::val V) {
  isl::pw_aff AV(A.domain(), V);
  return A % AV;
}

inline isl::pw_aff operator%(int i, isl::pw_aff A) {
  isl::ctx ctx = A.get_ctx();
  return isl::val(ctx, i) % A;
}

inline isl::pw_aff operator%(isl::pw_aff A, int i) {
  isl::ctx ctx = A.get_ctx();
  return A % isl::val(ctx, i);
}

inline isl::set operator==(isl::pw_aff A, int i) {
  return A.eq_set(0 * A + i);
}

inline isl::set operator==(isl::pw_aff A, isl::pw_aff B) {
  return A.eq_set(B);
}

inline isl::set operator>=(isl::pw_aff A, isl::pw_aff B) {
  return A.ge_set(B);
}

inline isl::set operator&&(isl::set A, isl::set B) {
  return A.intersect(B);
}

inline isl::set operator||(isl::set A, isl::set B) {
  return A.unite(B);
}

#endif
