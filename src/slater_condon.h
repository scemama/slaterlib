#ifndef SLATER_CONDON_H
#define SLATER_CONDON_H

#define ORBITAL_SHIFT 1

/*@ assigns \result;
    assigns \result \from x_0; */
extern int __builtin_ctzll (unsigned long long x_0);

/*@ assigns \result;
    assigns \result \from x_0; */
extern int __builtin_popcountll (unsigned long long x_0);

#define popcnt(X) __builtin_popcountll((unsigned long long) X)
#define trailz(X) __builtin_ctzll((unsigned long long) X)

#define NONE  (exc_number_t) 123456

typedef unsigned int         bucket_t;
typedef unsigned int         exc_number_t;
typedef unsigned int         orbital_t;

typedef enum { alpha=0, beta=1, alpha_beta=2 } spin_t;

typedef unsigned long long   determinant_t;

#define NORB_PER_INT (orbital_t) (8*sizeof(determinant_t) )
#define NORB_PER_INT_SHIFT trailz( (unsigned long long) NORB_PER_INT )

#define STRINGIFY(x) XSTRINGIFY(x)
#define XSTRINGIFY(x) #x

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

#endif
