#define main main_c3c
#include "c3c.cu"
#undef main
static Prog hand15(){ Prog P = seed();
  P.op[2] = NOP; P.op[5] = FMA_; P.d[5] = 13; P.a[5] = 0; P.b[5] = 12; P.c[5] = 4;
  P.op[6] = FMA_; P.d[6] = 13; P.a[6] = 0; P.b[6] = 13; P.c[6] = 10; return P; }
int main(){
  Backend B; STRIDE = 1;
  Prog S = seed(), H = hand15(), T = hand15();
  // « toutes astuces connues » : réduction d'argument en 3 opérations (fma)
  memset(T.opA, 0, MAXA);
  T.opA[0] = FMA_; T.dA[0] = 7; T.aA[0] = 0; T.bA[0] = 2; T.cA[0] = 1;   // tk = fma(x, 128, MAGIC)
  T.opA[1] = SUB;  T.dA[1] = 8; T.aA[1] = 7; T.bA[1] = 1;                // kf = tk - MAGIC
  T.opA[2] = FMA_; T.dA[2] = 10; T.aA[2] = 8; T.bA[2] = 4; T.cA[2] = 0;  // h = fma(kf, -2^-7, x)
  T.otk = 7; T.oh = 10;
  const char* nm[3] = {"graine (5 + 16)", "main, coeur à 15 (5 + 15)", "main, toutes astuces (3 + 15)"};
  Prog* Ps[3] = {&S, &H, &T};
  for (int q = 0; q < 3; q++){ unsigned w, s; std::vector<uint32_t> f; B.full(*Ps[q], w, s, f);
    printf("%-32s : %d opérations vivantes, %u erreurs, %u renvois (2^32 entrées)\n", nm[q], live_count(*Ps[q]), w, s); }
  FILE* f = fopen("hand_progs.cuh", "w");
  fprintf(f, "#pragma once\n#include \"fexp2.h\"\n%s%s", to_c(H, "hand15_fast").c_str(), to_c(T, "handall_fast").c_str()); fclose(f);
}
