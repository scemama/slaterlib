#include "slater_condon.h"
#include "slater_condon_simple.h"
#include "debug.h"
#include <string.h>
#include <assert.h>

#define MAXTESTS 10000



unsigned int random_int()
{
  unsigned int result;
  result = (rand() & (NORB_PER_INT-1));
  return result;
}


determinant_t random_det(orbital_t N_orb)
{
  orbital_t r;
  determinant_t result;

  result = (determinant_t) 0;

  while (popcnt(result) < N_orb)
  {
    r = random_int();
    result |= ( ((determinant_t)1) << r);
  }
  return result;
}


unsigned int random_single_exc(bucket_t N_int,
                    determinant_t d1[N_int],
                    determinant_t d2[N_int])
{
  orbital_t list_occ[N_int*NORB_PER_INT];
  orbital_t list_virt[N_int*NORB_PER_INT];
  bucket_t k;
  orbital_t N_orb, particle, pos, pos_particle, hole;
  unsigned int nperm;

  for (k=(bucket_t)0 ; k<N_int ; k++)
    d2[k] = ~d1[k];
  N_orb = to_orbital_list(N_int, d2, list_virt);
  assert (N_orb > 0);

  pos_particle = rand() % N_orb;
  particle = list_virt[ pos_particle ];

  N_orb = to_orbital_list(N_int, d1, list_occ);
  assert (N_orb > 0);

  pos = pos_particle;
  while (pos == pos_particle)
    pos  = rand() % N_orb;

  assert(pos != pos_particle);

  hole = list_occ[pos];
  list_occ[pos] = particle;
  of_orbital_list(N_int, N_orb, d2, list_occ);
  
  nperm = 0;
  if (particle > hole) {
    while (pos < N_orb-1) {
      pos++;
      if (list_occ[pos] < particle) nperm++;
    }
  } else {
    while (pos > 0) {
      pos--;
      if (list_occ[pos] > particle) nperm++;
    }
  }
  return nperm;
}


unsigned int random_double_exc(bucket_t N_int,
                    determinant_t d1[N_int],
                    determinant_t d2[N_int])
{
  determinant_t d15[N_int];
  unsigned int failed,i;
  unsigned int nperm, nperm2;

  failed = 1;
  nperm = random_single_exc(N_int,d1,d15);
  while (failed) {
    nperm2 = random_single_exc(N_int,d15,d2);
    for (i=1 ; i<N_int ; i++) {
      if (d2[i] != d1[i]) {
        failed = 0;
        break;
      }
    }
  }
  return nperm + nperm2;
}


/* Tests */

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
      err(1,"test_trailz: %llu  %u  %u\n", x, trailz_simple(x), popcnt(x));
    }
  }
  x = ((determinant_t) 0);
  if ( trailz_simple(x) != NONE ) test_ok = -1;
  return test_ok;
}


int test_popcnt()
{
  unsigned int itest;
  determinant_t x;
  int test_ok;

  test_ok = 0;
  for (itest=0 ; itest<MAXTESTS ; itest++)
  {
    x = random_det( (orbital_t) (random_int()));
    if ( popcnt_simple(x) != popcnt(x) ) {
      test_ok = -1;
      itest = MAXTESTS;
      err(1,"test_popcnt: %llu  %u  %u\n", x, popcnt_simple(x), popcnt(x));
    }
  }
  return test_ok;
}


int test_exc_degree(bucket_t N_int)
{
  unsigned int    itest;
  orbital_t       r;
  int             i,j,k,l,m;
  determinant_t   d1[N_int];
  determinant_t   d2[N_int];
  int             test_ok;

  test_ok = 0;
  for (itest=0 ; itest<MAXTESTS/10 ; itest++)
  {
    /* Excitation on one int */
    for (k=0 ; k<N_int ; k++)
    {
      r = random_int();
      r = r < 4 ? 4 : r-2;  // Be sure that we have at least two holes and particles
      memset((determinant_t*) d1, 0, N_int*sizeof(determinant_t));
      d1[k] = random_det(r);
      for (j=0 ; j<N_int ; j++)
      {
        memset((determinant_t*) d2, 0, N_int*sizeof(determinant_t));
        d2[j] = random_det(r);
        if (exc_degree(N_int, d1, d2) != exc_degree_simple(N_int, d1, d2))
        {
          debug_det(N_int, d1);
          debug_det(N_int, d2);
          test_ok = -1;
          itest = MAXTESTS;
          err(1,"Failure in exc_degree on same int : %d != %d", (int) exc_degree(N_int, d1, d2), (int) exc_degree_simple(N_int, d1, d2));
        }
      }
    }

    /* Excitation on two ints */
    for (k=0 ; k<N_int-1 ; k++)
    {
      r = random_int();
      r = r < 4 ? 4 : r-2;
      memset((determinant_t*) d1, 0, N_int*sizeof(determinant_t));
      d1[k] = random_det(r);
      for (j=k+1 ; j<N_int ; j++)
      {
        memset((determinant_t*) d1+k+1, 0, (N_int-k-1)*sizeof(determinant_t));
        d1[j] = random_det(r-(r/2));
        for (m=0 ; m<N_int-1 ; m++)
        {
          memset((determinant_t*) d2, 0, N_int*sizeof(determinant_t));
          d2[m] = random_det(r);
          for (i=m+1 ; i<N_int ; i++)
          {
            memset((determinant_t*) d2+m+1, 0, (N_int-m-1)*sizeof(determinant_t));
            d2[i] = random_det(r-(r/2));
            if (exc_degree(N_int, d1, d2) != exc_degree_simple(N_int, d1, d2))
            {
              debug_det(N_int, d1);
              debug_det(N_int, d2);
              test_ok = -1;
              itest = MAXTESTS;
              err(1,"Failure in exc_degree on 2 ints : %d != %d", (int) exc_degree(N_int, d1, d2), (int) exc_degree_simple(N_int, d1, d2));
            }
          }
        }
      }
    }

    /* Excitation on three ints */
    for (k=0 ; k<N_int ; k++)
    {
      r = random_int();
      r = r < 4 ? 4 : r-2;
      for (l=0 ; l<N_int ; l++) d1[l] = random_det(r);
      d1[k] = (determinant_t)0;
      for (j=0 ; j<N_int ; j++)
      {
        for (l=0 ; l<N_int ; l++) d2[l] = random_det(r);
        d2[j] = (determinant_t)0;
        if (exc_degree(N_int, d1, d2) != exc_degree_simple(N_int, d1, d2))
        {
          debug_det(N_int, d1);
          debug_det(N_int, d2);
          test_ok = -1;
          itest = MAXTESTS;
          err(1,"Failure in exc_degree on 3 ints : %d != %d", (int) exc_degree(N_int, d1, d2), (int) exc_degree_simple(N_int, d1, d2));
        }
      }
    }

    /* Excitation on four ints */
    for (k=0 ; k<N_int ; k++)
    {
      r = random_int();
      r = r < 4 ? 4 : r-2;
      for (l=0 ; l<N_int ; l++) {
        d1[l] = random_det(r);
        d2[l] = random_det(r);
      }
      if (exc_degree(N_int, d1, d2) != exc_degree_simple(N_int, d1, d2))
      {
        debug_det(N_int, d1);
        debug_det(N_int, d2);
        test_ok = -1;
        itest = MAXTESTS;
        err(1,"Failure in exc_degree on 3 ints : %d != %d", (int) exc_degree(N_int, d1, d2), (int) exc_degree_simple(N_int, d1, d2));
      }
    }
  }
  return test_ok;
}


int test_to_of_list(bucket_t N_int)
{
  unsigned int    itest;
  orbital_t       r;
  determinant_t   d1[N_int];
  determinant_t   d2[N_int];
  orbital_t       N_orb;
  orbital_t       list[N_int*NORB_PER_INT];
  bucket_t l;
  int test_ok;

  test_ok = 0;
  for (itest=0 ; itest<MAXTESTS ; itest++)
  {
    r = random_int();
    r = r < 4 ? 4 : r-2;
    for (l=(bucket_t) 0 ; l<N_int ; l++)
      d1[l] = random_det(r);
    N_orb = to_orbital_list(N_int, d1, list);
    of_orbital_list(N_int, N_orb, d2, list);
    if (exc_degree(N_int, d1, d2) != 0)
    {
        debug_det(N_int, d1);
        debug_det(N_int, d2);
        itest = MAXTESTS;
        test_ok = -1;
        err(1,"Failure in to/of_list : %d != 0", (int) exc_degree(N_int, d1, d2) );
    }
  }
  return test_ok;
}


int test_get_holes_particles(bucket_t N_int)
{
  int test_ok;
  unsigned int    itest;
  determinant_t   d1[N_int];
  determinant_t   d2[N_int];
  orbital_t       r;
  bucket_t        l;
  exc_number_t    exc, exc2;
  orbital_t       holes[2], holes2[2];
  orbital_t       particles[2], particles2[2];

  test_ok = 0;
  for (itest=0 ; itest<MAXTESTS ; itest++)
  {
    r = random_int();
    r = r < 4 ? 4 : r-2;
    for (l=(bucket_t)0 ; l<N_int ; l++)
      d1[l] = random_det(r);


    /* Single excitations */
    random_single_exc(N_int,d1,d2);

    exc = get_holes(N_int, d1, d2, holes);
    if (exc != 1)
    {
        debug_det(N_int, d1);
        debug_det(N_int, d2);
        test_ok = -1;
        itest = MAXTESTS;
        err(1,"Failure in random_single_exc: %d", (int) exc_degree(N_int, d1, d2) );
    }

    exc2 = get_holes_simple(N_int, d1, d2, holes2);
    if ( (exc2 != exc) || (holes[0] != holes2[0]) )
    {
        debug_det(N_int, d1);
        debug_det(N_int, d2);
        test_ok = -1;
        itest = MAXTESTS;
        err(1,"Failure in get_holes: %d %d", (int) holes[0], (int) holes2[0]);
    }

    exc2 = get_particles(N_int, d1, d2, particles);
    if (exc2 != exc)
    {
        debug_det(N_int, d1);
        debug_det(N_int, d2);
        test_ok = -1;
        itest = MAXTESTS;
        err(1,"Failure in get_particles: %d", (int) exc_degree(N_int, d1, d2) );
    }

    exc2 = get_particles_simple(N_int, d1, d2, particles2);
    if ( (exc2 != exc) || (particles[0] != particles2[0]) )
    {
        debug_det(N_int, d1);
        debug_det(N_int, d2);
        test_ok = -1;
        itest = MAXTESTS;
        err(1,"Failure in get_particles: %d, %d=%d", (int) exc, (int) particles[0], (int) particles2[0]);
    }

    exc2 = get_holes_particles(N_int,d1,d2,holes2,particles2);
    if ( (exc2 != exc) || (particles[0] != particles2[0]) || (holes[0] != holes2[0]) )
    {
        debug_det(N_int, d1);
        debug_det(N_int, d2);
        test_ok = -1;
        itest = MAXTESTS;
        err(1,"Failure in get_holes_particles: %d, %d=%d %d=%d", (int) exc,
          (int) holes[0], (int) holes2[0], (int) particles[0], (int) particles2[0]);
    }


    /* Double excitations */
    while (exc == 1) {
      random_double_exc(N_int,d1,d2);
      exc = get_holes(N_int, d1, d2, holes);
    }

    if ( (exc != 2) )
    {
        debug_det(N_int, d1);
        debug_det(N_int, d2);
        test_ok = -1;
        itest = MAXTESTS;
        err(1,"Failure in random_double_exc: %d", (int) exc_degree(N_int, d1, d2) );
    }

    exc2 = get_holes_simple(N_int, d1, d2, holes2);
    if ( (exc2 != exc) || (holes[0] != holes2[0]) || (holes[1] != holes2[1]) ) 
    {
        debug_det(N_int, d1);
        debug_det(N_int, d2);
        test_ok = -1;
        itest = MAXTESTS;
        err(1,"Failure in get_holes: %d, %d=%d, %d=%d", (int) exc2, 
          (int) holes[0], (int) holes2[0], (int) holes[1], (int) holes2[1]);
    }

    exc2 = get_particles(N_int, d1, d2, particles);
    if ( (exc2 != exc) )
    {
        debug_det(N_int, d1);
        debug_det(N_int, d2);
        test_ok = -1;
        itest = MAXTESTS;
        err(1,"Failure in get_particles: %d", (int) exc_degree(N_int, d1, d2));
    }

    exc2 = get_particles_simple(N_int, d1, d2, particles2);
    if ( (exc2 != exc) || (particles[0] != particles2[0]) || (particles[1] != particles2[1]) ) 
    {
        debug_det(N_int, d1);
        debug_det(N_int, d2);
        test_ok = -1;
        itest = MAXTESTS;
        err(1,"Failure in get_particles: %d,  %d=%d %d=%d", (int) exc2, 
          (int) particles[0], (int) particles2[0], (int) particles[1], (int) particles2[1]);
    }

    exc2 = get_holes_particles(N_int,d1,d2,holes2,particles2);
    if ( (exc2 != exc) || (particles[0] != particles2[0]) || (holes[0] != holes2[0]) 
                       || (particles[1] != particles2[1]) || (holes[1] != holes2[1]) )
    {
        debug_det(N_int, d1);
        debug_det(N_int, d2);
        test_ok = -1;
        itest = MAXTESTS;
        err(1,"Failure in get_holes_particles: %d, %d=%d %d=%d %d=%d %d=%d", (int) exc,
          (int) holes[0], (int) holes2[0], (int) particles[0], (int) particles2[0],
          (int) holes[1], (int) holes2[1], (int) particles[1], (int) particles2[1]);
    }
  }
  return test_ok;
}


int test_nperm(bucket_t N_int)
{
  int test_ok;
  unsigned int    itest;
  determinant_t   d1[N_int];
  determinant_t   d2[N_int];
  orbital_t       r;
  bucket_t        l;
  orbital_t       holes[2], particles[2];
  unsigned int    nperm, nperm2;

  test_ok = 0;
  for (itest=0 ; itest<MAXTESTS ; itest++)
  {
    r = random_int();
    r = r < 4 ? 4 : r-2;
    for (l=(bucket_t)0 ; l<N_int ; l++)
      d1[l] = random_det(r);


    /* Single excitations */
    nperm  = random_single_exc(N_int,d1,d2);
    get_holes_particles(N_int, d1, d2, holes, particles);
    nperm2 = get_nperm_single(N_int, d1, d2, holes, particles);
    if (nperm != nperm2)
    {
        debug_det(N_int, d1);
        debug_det(N_int, d2);
        test_ok = -1;
        itest = MAXTESTS;
        printf("%d -> %d\n",holes[0],particles[0]);
        err(1,"Failure in n_perm_single: %d %d", (int) nperm, (int) nperm2);
    }
  }
  return test_ok;
}


/* Main */

#include <time.h>
int main(int agrc, char** argv)
{
  bucket_t        N_int = 4;
  int             test_ok;

srand(time(NULL));
  test_ok = 0;
 
  printf("NORB_PER_INT_SHIFT: %d\n", (int) NORB_PER_INT_SHIFT);
  printf("NORB_PER_INT      : %d\n", (int) NORB_PER_INT);

#define TEST(X,Y) \
  if (test_##X (Y) != 0) { \
    printf("\t\t[ Failed ]\n"); \
    test_ok = -1; } \
  else \
    printf("\t\t[   OK   ]\n")

  printf("Trailz                "); TEST(trailz,);
  printf("Popcnt                "); TEST(popcnt,);
  printf("exc_degree            "); TEST(exc_degree,N_int);
  printf("to/of_list            "); TEST(to_of_list,N_int);
  printf("get_holes/particles   "); TEST(get_holes_particles,N_int);
  printf("n_perm                "); TEST(nperm,N_int);

  return test_ok;
}





