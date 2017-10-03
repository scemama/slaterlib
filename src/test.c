#include "slater_condon.h"
#include <stdio.h>
#include <stdlib.h>
#include <err.h>




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

int test_trailz()
{
  unsigned int k;
  determinant_t x;
  int test_ok;

  test_ok = 0;
  for (k=0 ; k<NORB_PER_INT ; k++)
  {
    x = ((determinant_t) 1) << k;
    if ( trailz_simple(x) != k ) test_ok = -1;
    if ( trailz_simple(x) != trailz(x) ) {
      test_ok = -1;
      errx(1,"test_trailz: %llu  %u  %u\n", x, trailz_simple(x), popcnt(x));
    }
  }
  x = ((determinant_t) 0);
  if ( trailz_simple(x) != NONE ) test_ok = -1;
  return test_ok;
}

determinant_t random_det()
{
  return
    (((determinant_t) rand() <<  0) & 0x000000000000FFFFull) | 
    (((determinant_t) rand() << 16) & 0x00000000FFFF0000ull) | 
    (((determinant_t) rand() << 32) & 0x0000FFFF00000000ull) |
    (((determinant_t) rand() << 48) & 0xFFFF000000000000ull);
}

int test_popcnt()
{
  unsigned int k;
  determinant_t x;
  int test_ok;

  test_ok = 0;
  for (k=0 ; k<1000000 ; k++)
  {
    x = random_det();
    if ( popcnt_simple(x) != popcnt(x) ) {
      test_ok = -1;
      errx(1,"test_popcnt: %llu  %u  %u\n", x, popcnt_simple(x), popcnt(x));
    }
  }
  return test_ok;
}





int main(int agrc, char** argv)
{
   bucket_t        N_int = (bucket_t) 1;
   determinant_t   d1[1];
   determinant_t   d2[1];
   exc_number_t    exc;
   orbital_t       holes[2];
   orbital_t       particles[2];
   orbital_t       iorb;

  test_trailz();
  test_popcnt();
    



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





