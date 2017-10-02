#include "slater_condon.h"


/* Returns the excitation degree between two determinants.
 */

exc_number_t exc_degree(bucket_t N_int,
                        determinant_t d1[alpha_beta][N_int],
                        determinant_t d2[alpha_beta][N_int])
{
  bucket_t       i;
  exc_number_t   result;

  result = (exc_number_t) 0;
  for (i=(bucket_t) 0 ; i<N_int ; i++) {
    result += (exc_number_t) ( popcnt(d1[alpha][i]^d2[alpha][i])
                             + popcnt(d1[beta ][i]^d2[beta ][i]) );
  }
  return (result >> 1);
}




/* Returns the number of holes in the d1 -> d2 excitation.
 * `holes` is the list of orbital indices.
 */
exc_number_t get_holes(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int],
                        orbital_t*    holes)
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
                        orbital_t*    particles)
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
                        orbital_t*    holes,    
                        orbital_t*    particles)
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
