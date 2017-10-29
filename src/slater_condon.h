#ifndef SLATER_CONDON_H
#define SLATER_CONDON_H

#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <assert.h>

#define ORBITAL_SHIFT 1
#define INT_SIZE      64


#define NONE   (exc_number_t) 123456
#define MAXORB (orbital_t) 10000

typedef unsigned int         bucket_t;
typedef unsigned int         exc_number_t;
typedef unsigned int         orbital_t;

typedef enum { alpha=0, beta=1, alpha_beta=2 } spin_t;
typedef enum { phase_p=0, phase_m=1 } phase_t;
/* static double phase_double[2] = { 1., -1. };
 */

typedef struct {
  orbital_t    holes[2];
  orbital_t    particles[2];
  exc_number_t exc_degree;
  phase_t      phase;
} excitation_operator_t;
  

/* Popcount and trailz */
#if INT_SIZE == 64

  typedef unsigned long long   determinant_t;

  /*@ assigns \result;
      assigns \result \from x_0; */
  extern int __builtin_popcountll (unsigned long long x_0);
  #define popcnt(X) __builtin_popcountll((unsigned long long) X)

  /*@ assigns \result;
      assigns \result \from x_0; */
  extern int __builtin_ctzll (unsigned long long x_0);
  #define trailz(X) __builtin_ctzll((unsigned long long) X)

#elif INT_SIZE == 32

  /*@ assigns \result;
      assigns \result \from x_0; */
  extern int __builtin_popcountl (unsigned long x_0);
  #define popcnt(X) __builtin_popcountl((unsigned long) X)

  /*@ assigns \result;
      assigns \result \from x_0; */
  extern int __builtin_ctzl(unsigned long x_0);
  #define trailz(X) __builtin_ctzl((unsigned long) X)

#elif INT_SIZE == 16

  /*@ assigns \result;
      assigns \result \from x_0; */
  extern int __builtin_popcount (unsigned int x_0);
  #define popcnt(X) __builtin_popcount((unsigned int) X)

  /*@ assigns \result;
      assigns \result \from x_0; */
  extern int __builtin_ctz (unsigned int x_0);
  #define trailz(X) __builtin_ctz((unsigned int) X)

#else

 #error("Invalid INT_SIZE")

#endif

#define NORB_PER_INT ( (orbital_t) (8*sizeof(determinant_t) ) )
#define NORB_PER_INT_SHIFT ( trailz( NORB_PER_INT ) )



/* Optimized functions */

exc_number_t exc_degree(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int]);

exc_number_t get_holes(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int],
                        orbital_t holes[2]);

exc_number_t get_particles(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int],
                        orbital_t particles[2]);

exc_number_t get_holes_particles(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int],
                        orbital_t holes[2],    
                        orbital_t particles[2]);

unsigned int get_nperm_single(bucket_t N_int,
                 determinant_t d1[N_int],
                 determinant_t d2[N_int],
                 orbital_t hole[2],
                 orbital_t particle[2]);

unsigned int get_nperm_double(bucket_t N_int,
                 determinant_t d1[N_int],
                 determinant_t d2[N_int],
                 orbital_t hole[2],
                 orbital_t particle[2]);


orbital_t to_orbital_list(bucket_t N_int,
                     determinant_t d1[N_int],
                     orbital_t* list);

void of_orbital_list(bucket_t N_int,
                     orbital_t N_orb, 
                     determinant_t d1[N_int],
                     orbital_t list[N_orb]);


/*@
predicate
  ValidDeterminant (determinant_t * a, integer N_int) = 
      \valid_read(a+(0..N_int-1)) &&
      0 < N_int < MAXORB/NORB_PER_INT;
predicate
  ValidDeterminants (determinant_t * a, determinant_t * b, integer N_int) =
      \valid_read(a+(0..N_int-1)) &&
      \valid_read(b+(0..N_int-1)) &&
      0 < N_int < (MAXORB/NORB_PER_INT);
*/

#endif
