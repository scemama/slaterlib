#include "slater_condon.h"


/* Returns the excitation degree between two determinants.
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
  return (result >> 1);
}




/* Returns the number of holes in the d1 -> d2 excitation.
 * `holes` is the list of orbital indices.
 */
exc_number_t get_holes(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int],
                        orbital_t holes[2])
{
  exc_number_t    n_holes;
  bucket_t        i;
  unsigned int    k;
  unsigned int    shift;
  unsigned int    pos;
  determinant_t   tmp;

  k = 0;
  shift = ORBITAL_SHIFT;

  for (i=(bucket_t) 0 ; i<N_int ; i++)
  {
      tmp = (d1[i]^d2[i]) & d1[i];

      while (tmp != (determinant_t) 0) 
      {
          pos = trailz(tmp);
          holes[k] = (orbital_t) (pos + shift);
          tmp ^= ( ((determinant_t) 1) << pos);
          k++;

      }
      shift += (unsigned int) (8*sizeof(determinant_t));

  }
  n_holes = (exc_number_t) k;
  return n_holes;
}



/* Returns the number of particles in the d1 -> d2 excitation.
 * `particles` is the list of orbital indices.
 */
exc_number_t get_particles(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int],
                        orbital_t particles[2])
{
  exc_number_t    n_particles;
  bucket_t        i;
  unsigned int    k;
  unsigned int    shift;
  unsigned int    pos;
  determinant_t   tmp;

  k = 0;
  shift = ORBITAL_SHIFT;

  for (i=(bucket_t) 0 ; i<N_int ; i++)
  {
      tmp = (d1[i]^d2[i]) & d2[i];

      while (tmp != (determinant_t) 0) 
      {
          pos = trailz(tmp);
          particles[k] = (orbital_t) (pos + shift);
          tmp ^= ( ((determinant_t) 1) << pos);
          k++;

      }
      shift += (unsigned int) (8*sizeof(determinant_t));

  }
  n_particles = (exc_number_t) k;
  return n_particles;
}





/* Returns the number of holes or particles in the d1 -> d2 excitation.
 * `holes` and `particles` are lists of orbital indices.
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

  for (i=(bucket_t) 0 ; i<N_int ; i++)
  {
      tmp_xor = (d1[i]^d2[i]);

      tmp = tmp_xor & d1[i];
      while (tmp != (determinant_t) 0) 
      {
          pos = trailz(tmp);
          holes[k_holes] = (orbital_t) (pos + shift);
          tmp ^= ( ((determinant_t) 1) << pos);
          k_holes++;
      }

      tmp = tmp_xor & d2[i];
      while (tmp != (determinant_t) 0) 
      {
          pos = trailz(tmp);
          particles[k_particles] = (orbital_t) (pos + shift);
          tmp ^= ( ((determinant_t) 1) << pos);
          k_particles++;
      }

      shift += (unsigned int) (8*sizeof(determinant_t));
  }
  n_particles = (exc_number_t) k_particles;
  /*@ assert k_particles == k_holes; */
  return n_particles;
}




/* Returns the number of permutations obtained from the single
 * excitation hole[0] -> particle[0]
 */
unsigned int get_nperm_single(bucket_t N_int,
                 determinant_t d1[N_int],
                 determinant_t d2[N_int],
                 orbital_t hole[2],    
                 orbital_t particle[2])
{
  orbital_t  high, low;
  bucket_t  j,k,l;
  unsigned short m,n;
  unsigned int nperm;
  determinant_t mask[N_int];

  high = ( (particle[0] > hole[0]) ? particle[0] : hole[0]    ) - ORBITAL_SHIFT;
  low  = ( (particle[0] > hole[0]) ? hole[0]     : particle[0]) - ORBITAL_SHIFT;

  k = (bucket_t) (high >> NORB_PER_INT_SHIFT );
  j = (bucket_t) (low  >> NORB_PER_INT_SHIFT );
  m = (unsigned short) (high & (NORB_PER_INT - 1) );
  n = (unsigned short) (low  & (NORB_PER_INT - 1) );

  for (l=j ; l<k ; l++)
      mask[l] = ~((determinant_t)0);

  mask[k] = ( ((determinant_t)1) << m) - ((determinant_t)1) ;
  mask[j] &= ~(((determinant_t)1) << (n+1)) + ((determinant_t)1) ;

  nperm = (unsigned int) 0;
  for (l=j ; l<=k ; l++)
     nperm += popcnt( d1 & mask[l] ); 

  return nperm;
}



/* Returns the number of permutations obtained from the double
 * excitation hole[0] -> particle[0] and hole[1] -> particle[1].
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



