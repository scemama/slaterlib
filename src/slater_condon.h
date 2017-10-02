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


exc_number_t exc_degree(bucket_t N_int,
                        determinant_t d1[alpha_beta][N_int],
                        determinant_t d2[alpha_beta][N_int]);

exc_number_t get_holes(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int],
                        orbital_t*    holes);

exc_number_t get_particles(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int],
                        orbital_t*    particles);

exc_number_t get_holes_particles(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int],
                        orbital_t*    holes,    
                        orbital_t*    particles);

#endif
