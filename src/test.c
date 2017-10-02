#include "slater_condon.h"
#include <stdio.h>


#define BYTE_TO_BINARY_PATTERN "%c%c%c%c%c%c%c%c"
#define BYTE_TO_BINARY(byte)  \
  (byte & 0x0080 ? '1' : '0'), \
  (byte & 0x0040 ? '1' : '0'), \
  (byte & 0x0020 ? '1' : '0'), \
  (byte & 0x0010 ? '1' : '0'), \
  (byte & 0x0008 ? '1' : '0'), \
  (byte & 0x0004 ? '1' : '0'), \
  (byte & 0x0002 ? '1' : '0'), \
  (byte & 0x0001 ? '1' : '0') 


int main(int agrc, char** argv)
{
   bucket_t        N_int = (bucket_t) 1;
   determinant_t   d1[alpha_beta][1];
   determinant_t   d2[alpha_beta][1];
   exc_number_t    exc;
   orbital_t       holes[8*sizeof(determinant_t)];
   orbital_t       particles[8*sizeof(determinant_t)];
   orbital_t       iorb;



   d1[alpha][0] = 0x03ff;    d1[beta][0] = 0x0000;
   d2[alpha][0] = 0x0cff;    d2[beta][0] = 0x0000;

   exc = exc_degree(N_int, d1, d2);
   printf("Excitation degree: %d\n", (int) exc);

   exc = get_holes(N_int, d1[alpha], d2[alpha], holes);
   printf("Number of holes : %d\n", (int) exc);
   for (iorb=(orbital_t)0 ; iorb < exc ; iorb++)
      printf("%d ", (int) holes[iorb]);
   printf("\n");

   exc = get_particles(N_int, d1[alpha], d2[alpha], particles);
   printf("Number of particles: %d\n", (int) exc);
   for (iorb=(orbital_t)0 ; iorb < exc ; iorb++)
      printf("%d ", (int) particles[iorb]);
   printf("\n");

   exc = get_holes_particles(N_int, d1[alpha], d2[alpha], holes, particles);
   printf("Number of holes/particles : %d\n", (int) exc);
   for (iorb=(orbital_t)0 ; iorb < exc ; iorb++)
      printf("%d ", (int) holes[iorb]);
   printf("\n");
   for (iorb=(orbital_t)0 ; iorb < exc ; iorb++)
      printf("%d ", (int) particles[iorb]);
   printf("\n");

   return 0;
}
