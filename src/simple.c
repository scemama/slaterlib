#include "slater_condon.h"


/* Returns the excitation degree between two determinants.
 */

exc_number_t exc_degree_simple(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int])
{
  bucket_t       i;
  exc_number_t   result;

  result = (exc_number_t) 0;
  for (i=(bucket_t) 0 ; i<N_int ; i++) {
    result += (exc_number_t) ( popcnt(d1[i]^d2[i]) );
  }
  return result/2;
}


unsigned int trailz_simple(determinant_t x)
{
  unsigned int k;
  if ( x == (determinant_t) 0 ) {
    k= NONE;
  } else {
    k=0;
    while ( (k < NORB_PER_INT) &&
            ( ( (~x) & (((unsigned long long)1) << k)) != ((unsigned long long)0) ) )
      k++;
  }
  return k;
}


unsigned int popcnt_simple(determinant_t x)
{
  unsigned int k, i;
  k=0;
  for (i=0 ; i<NORB_PER_INT ; i++)
  {
    if ( (x & (((determinant_t) 1) << i)) != 0 )
      k++;
  }
  return k;
}


