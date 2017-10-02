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
   determinant_t   d1[1];
   determinant_t   d2[1];
   exc_number_t    exc;
   orbital_t       holes[2];
   orbital_t       particles[2];
   orbital_t       iorb;

  printf("NORB_PER_INT_SHIFT: %d\n", (int) NORB_PER_INT_SHIFT);
  printf("NORB_PER_INT      : %d\n", (int) NORB_PER_INT);



   d1[0] = 0x03ff;  
   d2[0] = 0x0cff;  

   exc = exc_degree(N_int, d1, d2);
   printf("Excitation degree: %d\n", (int) exc);

   exc = get_holes(N_int, d1, d2, holes);
   printf("Number of holes : %d\n", (int) exc);
   for (iorb=(orbital_t)0 ; iorb < exc ; iorb++)
      printf("%d ", (int) holes[iorb]);
   printf("\n");

   exc = get_particles(N_int, d1, d2, particles);
   printf("Number of particles: %d\n", (int) exc);
   for (iorb=(orbital_t)0 ; iorb < exc ; iorb++)
      printf("%d ", (int) particles[iorb]);
   printf("\n");

   exc = get_holes_particles(N_int, d1, d2, holes, particles);
   printf("Number of holes/particles : %d\n", (int) exc);
   for (iorb=(orbital_t)0 ; iorb < exc ; iorb++)
      printf("%d ", (int) holes[iorb]);
   printf("\n");
   for (iorb=(orbital_t)0 ; iorb < exc ; iorb++)
      printf("%d ", (int) particles[iorb]);
   printf("\n");

   exc = get_holes_particles(N_int, d1, d2, holes, particles);
   printf("Number of holes/particles : %d\n", (int) exc);
   for (iorb=(orbital_t)0 ; iorb < exc ; iorb++)
      printf("%d ", (int) holes[iorb]);
   printf("\n");
   for (iorb=(orbital_t)0 ; iorb < exc ; iorb++)
      printf("%d ", (int) particles[iorb]);
   printf("\n");

   return 0;
}
