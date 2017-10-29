#include "slater_condon.h"

/* Print a determinant for debugging purposes */
void debug_det(bucket_t N_int, determinant_t d[N_int])
{
  orbital_t list[N_int * NORB_PER_INT];
  orbital_t k;
  bucket_t  l;
  unsigned int nmax;

  
  to_orbital_list(N_int, d, list);

  nmax = 0;
  for (l=0 ; l<N_int ; l++)
  {
    nmax += popcnt(d[l]);
    for (k=0 ; k<NORB_PER_INT ; k++)
    {
      if ( (d[l] & ((determinant_t)1) << k) == 0 )
        printf("0");
      else
        printf("1");
    }
    printf(" ");
  }
  printf(":\n");
  for (k=0 ; k<nmax; k++)
    printf("%d ",list[k]);
  printf("\n");
}


