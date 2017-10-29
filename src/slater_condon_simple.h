#ifndef SLATER_CONDON_SIMPLE_H
#define SLATER_CONDON_SIMPLE_H

#include "slater_condon.h"

unsigned int trailz_simple(determinant_t x);
unsigned int popcnt_simple(determinant_t x);

exc_number_t exc_degree_simple(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int]);

exc_number_t get_holes_simple(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int],
                        orbital_t holes[2]);

exc_number_t get_particles_simple(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int],
                        orbital_t particles[2]);

/*
exc_number_t get_holes_particles_simple(bucket_t N_int,
                        determinant_t d1[N_int],
                        determinant_t d2[N_int],
                        orbital_t holes[2],    
                        orbital_t particles[2]);

unsigned int get_nperm_single_simple(bucket_t N_int,
                 determinant_t d1[N_int],
                 determinant_t d2[N_int],
                 orbital_t hole[2],
                 orbital_t particle[2]);

unsigned int get_nperm_double_simple(bucket_t N_int,
                 determinant_t d1[N_int],
                 determinant_t d2[N_int],
                 orbital_t hole[2],
                 orbital_t particle[2]);
*/


#endif
