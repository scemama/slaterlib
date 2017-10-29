#include "slater_condon.h"
#include "debug.h"

/* Naive implementation of trailz
 */
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


/* Naive implementation of popcnt
 */
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




/* Returns the excitation degree between two determinants.
 */
exc_number_t exc_degree_simple(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int])
{
  orbital_t      i,j,norb,jmin;
  exc_number_t   result;

  orbital_t  list1[N_int*NORB_PER_INT];
  orbital_t  list2[N_int*NORB_PER_INT];

  norb = to_orbital_list(N_int, d1, list1);
  if (norb != to_orbital_list(N_int, d2, list2)) {
    printf("N_orb: %d %d\n", norb, to_orbital_list(N_int, d2, list2));
    debug_det(N_int,d1);
    debug_det(N_int,d2);
    errx(1,"exc_degree_simple: Wrong number of orbitals in determinants");
  }

  result = (exc_number_t) norb;
  jmin = (orbital_t) 0;
  for (i=(orbital_t)0 ; i<norb ; i++) {
    for (j=jmin ; j<norb ; j++) {
      if (list2[j] > list1[i])
        break;
      else if (list1[i] == list2[j]) {
        jmin = j;
        result--;
        break;
      }
    }
  }
  return result;
}


exc_number_t get_holes_simple(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int],
                        orbital_t holes[2])
{
  orbital_t      i,j,norb;
  orbital_t      k;

  orbital_t  list1[N_int*NORB_PER_INT];
  orbital_t  list2[N_int*NORB_PER_INT];
  int found;

  holes[0] = 0;
  holes[1] = 0;
  norb = to_orbital_list(N_int, d1, list1);
  if (norb != to_orbital_list(N_int, d2, list2)) {
    printf("N_orb: %d %d\n", norb, to_orbital_list(N_int, d2, list2));
    debug_det(N_int,d1);
    debug_det(N_int,d2);
    errx(1,"get_holes_simple: Wrong number of orbitals in determinants");
  }

  k = (orbital_t) 0;
  for (i=(orbital_t)0 ; i<norb ; i++) {
    found=0;
    for (j=(orbital_t)0 ; j<norb ; j++) {
      if (list1[i] == list2[j]) {
        found=1;
        break;
      }
    }
    if (found == 0) {
        holes[k] = (orbital_t) list1[i];
        k++;
    }
  }
  return k;
}



exc_number_t get_particles_simple(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int],
                        orbital_t particles[2])
{
  orbital_t      i,j,norb;
  orbital_t      k;

  orbital_t  list1[N_int*NORB_PER_INT];
  orbital_t  list2[N_int*NORB_PER_INT];
  int found;

  particles[0] = 0;
  particles[1] = 0;
  norb = to_orbital_list(N_int, d1, list1);
  if (norb != to_orbital_list(N_int, d2, list2)) {
    printf("N_orb: %d %d\n", norb, to_orbital_list(N_int, d2, list2));
    debug_det(N_int,d1);
    debug_det(N_int,d2);
    errx(1,"get_particles_simple: Wrong number of orbitals in determinants");
  }

  k = (orbital_t) 0;
  for (i=(orbital_t)0 ; i<norb ; i++) {
    found=0;
    for (j=(orbital_t)0 ; j<norb ; j++) {
      if (list2[i] == list1[j]) {
        found=1;
        break;
      }
    }
    if (found == 0) {
        particles[k] = (orbital_t) list2[i];
        k++;
    }
  }
  return k;
}

