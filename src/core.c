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





/* Returns the excitation degree between two determinants.
 */
/*@
   requires
      ValidDeterminants(d1,d2,N_int);
   ensures
      0 <= \result < (N_int*NORB_PER_INT)/2; 
   assigns \nothing;
*/ 
exc_number_t exc_degree(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int])
{
  bucket_t       i;
  exc_number_t   result;

  result = (exc_number_t) 0;
  for (i=(bucket_t) 0 ; i<N_int ; i++) {
    result += (exc_number_t) ( popcnt(d1[i]^d2[i]) );
  }
/*@ assert result % 2 == 0; */

  return (result >> 1);
}




/* Returns the number of holes in the d1 -> d2 excitation.
 * `holes` is the list of orbital indices.
 */
/*@ 
    requires
      ValidDeterminants(d1,d2,N_int) && \valid(holes+(0..1));

    behavior same_det:
      ensures \result == 0;
      assigns \nothing;

    behavior single_exc:
      ensures \result == 1;
      assigns holes[0];

    behavior double_exc:
      ensures \result == 2;
      assigns *(holes+(0..1));

    complete behaviors;
    disjoint behaviors;
*/ 
exc_number_t get_holes(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int],
                        orbital_t holes[2])
{
  exc_number_t    n_holes;
  bucket_t        i;
  orbital_t       shift;
  determinant_t   tmp;
  unsigned int    k;
  unsigned int    pos;

  k = 0;
  shift = ORBITAL_SHIFT;
/*@ assert ORBITAL_SHIFT == 0 || ORBITAL_SHIFT == 1; */

/*@ loop invariant (ORBITAL_SHIFT <= shift < NORB_PER_INT*N_int+ORBITAL_SHIFT) ; */
  for (i=(bucket_t) 0 ; i<N_int ; i++)
  {
      tmp = (d1[i]^d2[i]) & d1[i];

/*@   loop invariant (0 < pos <= NORB_PER_INT); */
      while (tmp != (determinant_t) 0) 
      {
          pos = trailz(tmp);
          holes[k] = ( (orbital_t) pos) + shift;
          tmp ^= ( ((determinant_t) 1) << pos);
          k++;
/*@       assert (0 < k <= 2); */
      }
      shift += NORB_PER_INT;
  }
  n_holes = (exc_number_t) k;
  return n_holes;
}



/* Returns the number of particles in the d1 -> d2 excitation.
 * `particles` is the list of orbital indices.
 */
/*@ 
    requires
      ValidDeterminants(d1,d2,N_int) && \valid(particles+(0..1));

    behavior same_det:
      ensures \result == 0;
      assigns \nothing;

    behavior single_exc:
      ensures \result == 1;
      assigns particles[0];

    behavior double_exc:
      ensures \result == 2;
      assigns *(particles+(0..1));

    complete behaviors;
    disjoint behaviors;
*/ 
exc_number_t get_particles(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int],
                        orbital_t particles[2])
{
  exc_number_t    n_particles;
  bucket_t        i;
  orbital_t       shift;
  determinant_t   tmp;
  unsigned int    k;
  unsigned int    pos;

  k = 0;
  shift = ORBITAL_SHIFT;
/*@ assert ORBITAL_SHIFT == 0 || ORBITAL_SHIFT == 1; */

/*@ loop invariant (ORBITAL_SHIFT <= shift < NORB_PER_INT*N_int+ORBITAL_SHIFT) ; */
  for (i=(bucket_t) 0 ; i<N_int ; i++)
  {
      tmp = (d1[i]^d2[i]) & d2[i];

/*@   loop invariant (0 < pos <= NORB_PER_INT); */
      while (tmp != (determinant_t) 0) 
      {
          pos = trailz(tmp);
          particles[k] = ( (orbital_t) pos) + shift;
          tmp ^= ( ((determinant_t) 1) << pos);
          k++;
/*@       assert (0 < k <= 2); */
      }
      shift += NORB_PER_INT;
  }
  n_particles = (exc_number_t) k;
  return n_particles;
}





/* Returns the number of holes or particles in the d1 -> d2 excitation.
 * `holes` and `particles` are lists of orbital indices.
 */
/*@ 
    requires
      ValidDeterminants(d1,d2,N_int) &&
      \valid(particles+(0..1)) &&
      \valid(holes+(0..1));

    behavior same_det:
      ensures \result == 0;
      assigns \nothing;

    behavior single_exc:
      ensures \result == 1;
      assigns particles[0], holes[0];

    behavior double_exc:
      ensures (\result == 2) && (holes[0] < holes[1]) && (particles[0] < particles[1]);
      assigns *(particles+(0..1)), *(holes+(0..1));

    complete behaviors;
    disjoint behaviors;
*/ 
exc_number_t get_holes_particles(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int],
                        orbital_t holes[2],    
                        orbital_t particles[2])
{
  exc_number_t    n_particles;
  bucket_t        i;
  unsigned int    k_holes, k_particles;
  unsigned int    shift;
  unsigned int    pos;
  determinant_t   tmp;
  determinant_t   tmp_xor;

  k_holes = 0;
  k_particles = 0;
  shift = ORBITAL_SHIFT;
/*@ assert ORBITAL_SHIFT == 0 || ORBITAL_SHIFT == 1; */

/*@ loop invariant (ORBITAL_SHIFT <= shift < NORB_PER_INT*N_int+ORBITAL_SHIFT) ; */
  for (i=(bucket_t) 0 ; i<N_int ; i++)
  {
      tmp_xor = (d1[i]^d2[i]);

      tmp = tmp_xor & d1[i];
/*@   loop invariant (0 < pos <= NORB_PER_INT); */
      while (tmp != (determinant_t) 0) 
      {
          pos = trailz(tmp);
          holes[k_holes] = (orbital_t) (pos + shift);
          tmp ^= ( ((determinant_t) 1) << pos);
          k_holes++;
/*@       assert (0 < k_holes <= 2); */
      }

      tmp = tmp_xor & d2[i];
/*@   loop invariant (0 < pos <= NORB_PER_INT); */
      while (tmp != (determinant_t) 0) 
      {
          pos = trailz(tmp);
          particles[k_particles] = (orbital_t) (pos + shift);
          tmp ^= ( ((determinant_t) 1) << pos);
          k_particles++;
/*@       assert (0 < k_particles <= 2); */
      }

      shift += NORB_PER_INT;
  }
  /*@ assert k_particles == k_holes; */
  n_particles = (exc_number_t) k_particles;
  return n_particles;
}




/* Returns the number of permutations obtained from the single
 * excitation hole[0] -> particle[0]
 */
/*@ 
    requires
      ValidDeterminants(d1,d2,N_int) &&
      \valid_read(particle)          &&
      \valid_read(hole);
    assigns \nothing;
    ensures \result >= 0;
*/ 
unsigned int get_nperm_single(bucket_t N_int,
                 determinant_t d1[N_int],
                 determinant_t d2[N_int],
                 orbital_t* hole,    
                 orbital_t* particle)
{
  orbital_t  high, low;
  bucket_t  j,k,l;
  unsigned int m,n;
  unsigned int nperm;
  determinant_t mask[N_int];

  high = ( (*particle > *hole) ? *particle : *hole    ) - ORBITAL_SHIFT;
  low  = ( (*particle > *hole) ? *hole     : *particle) - ORBITAL_SHIFT + 1;
/*@ assert ORBITAL_SHIFT <= low < high <= MAXORB+ORBITAL_SHIFT; */

  k = (bucket_t) (high >> NORB_PER_INT_SHIFT );
/*@ assert k == high / NORB_PER_INT;  */

  j = (bucket_t) (low  >> NORB_PER_INT_SHIFT );
/*@ assert j == low / NORB_PER_INT;  */

  m = (unsigned int) (high & (NORB_PER_INT - 1) );
/*@ assert m == high % NORB_PER_INT;  */

  n = (unsigned int) (low  & (NORB_PER_INT - 1) );
/*@ assert n == low % NORB_PER_INT;  */

/*@ assert (j<=k) ; */

  for (l=j ; l<k ; l++)
      mask[l] = ~((determinant_t)0);
  mask[k] = ( ((determinant_t)1) << m) - (determinant_t)1 ;
  mask[j] &= ( ~(((determinant_t)1) << n) + ((determinant_t)1) );

  nperm = (unsigned int) 0;
  for (l=j ; l<=k ; l++)
     nperm += popcnt( d1[l] & mask[l] ); 

  return nperm;
}



/* Returns the number of permutations obtained from the double
 * excitation hole[0] -> particle[0] and hole[1] -> particle[1].
 */
/*@ 
    requires
      ValidDeterminants(d1,d2,N_int) &&
      \valid_read(particle+(0..1))   &&
      \valid_read(hole+(0..1));
    assigns \nothing;
    ensures \result >= 0;
*/ 
unsigned int get_nperm_double(bucket_t N_int,
                 determinant_t d1[N_int],
                 determinant_t d2[N_int],
                 orbital_t hole[2],    
                 orbital_t particle[2])
{
  unsigned int nperm;

  nperm  = get_nperm_single(N_int, d1, d2, &hole[0], &particle[0]);
  nperm += get_nperm_single(N_int, d1, d2, &hole[1], &particle[1]);

  if ( (hole[1] < particle[0]) || (hole[0] > particle[1]) )
    nperm++;

  return nperm;
}



/* Returns the number of occupied orbitals, and an ordered list of occupied orbitals.
 */
/*  
    requires
      ValidDeterminant(d1,N_int) &&
      \valid(hole+(0..(N_int-1)*NORB_PER_INT));
    assigns (hole+(0..(N_int-1)*NORB_PER_INT));
*/ 
orbital_t to_orbital_list(bucket_t N_int,
                     determinant_t d1[N_int],
                     orbital_t* list)
{
  bucket_t        i;
  orbital_t       shift;
  determinant_t   tmp;
  unsigned int    k;
  unsigned int    pos;

  k = 0;
  shift = ORBITAL_SHIFT;
/*@ assert ORBITAL_SHIFT == 0 || ORBITAL_SHIFT == 1; */

/*@ loop invariant (ORBITAL_SHIFT <= shift < NORB_PER_INT*N_int+ORBITAL_SHIFT) ; */
  for (i=(bucket_t) 0 ; i<N_int ; i++)
  {
      tmp = d1[i];

/*@   loop invariant (0 < pos <= NORB_PER_INT); */
      while (tmp != (determinant_t) 0) 
      {
          pos = trailz(tmp);
          list[k] = ( (orbital_t) pos) + shift;
          tmp ^= ( ((determinant_t) 1) << pos);
          k++;
/*@       assert (0 < k <= N_int*NORB_PER_INT); */
      }
      shift += NORB_PER_INT;
  }
  return k;
}







void of_orbital_list(bucket_t N_int,
                     orbital_t N_orb, 
                     determinant_t d1[N_int],
                     orbital_t list[N_orb])
{
  bucket_t        i;
  orbital_t       pos, iorb;
  unsigned int    k;

  for (i=(bucket_t) 0 ; i<N_int ; i++) {
      d1[i] = (determinant_t) 0;
  }

  for (pos=((orbital_t)0) ; pos < N_orb ; pos++)
  {
    iorb = list[pos] - ORBITAL_SHIFT;
    i = (bucket_t) (iorb >> NORB_PER_INT_SHIFT );
    k = (unsigned int) (iorb & (NORB_PER_INT - 1) );
    d1[i] |= ((determinant_t)1) << k;
  }

}



