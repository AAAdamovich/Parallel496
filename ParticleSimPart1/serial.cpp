/*  Antony Adamovich, Dimitra Deliopoulos, Mahmoud Gudarzi
*    serial.cpp for XSEDE hw2 "Particle Simulator" - Serial version
*    West Chester University - CSC 496 - Dr. Lihn B. Ngo
*    Created: 04-MAR-2020 - Last Edited: 05-APR-2020 by Antony Adamovich
*    Please see https://github.com/AAAdamovich/Parallel496
*       for version tracking
*    Definitions
*       Cell: A subdivision of the simulation area, as a square
*       Frontier: An area close a cell's border where particles from 
*       two neighboring cells could potentially interact
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "common.h"

int main( int argc, char **argv )
{    
    // !! Testing Value !! Must be altered dynamically, NOT YET IMPLEMENTED
    int divisions = 3;
    // Corresponding to a MxM matrix where SIZE_DIV = M

    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;
    // The length of the square that defines the area in which the particles will be moving around in
    // Duplicate value from "common.cpp"
    double size;
    // The length of the squares that correspond to cells in the matrix subdivision of the simulation area
    double subsize;

    // Command line option finding
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    // Particle Initialization
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    // Obtain size information from common
    size = set_size( n );
    subsize = size / ((double)divisions);
    init_particles( n, particles );
    
    // Time when simulation begins
    double simulation_time = read_timer( );
	
    // How long the simulation is to run is based on NSTEPS
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        davg = 0.0;
        dmin = 1.0;
        
        // Application of near-particle forces (acceleration)
        for( int i = 0; i < n; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
            for (int j = 0; j < n; j++ )
				apply_force( particles[i], particles[j],&dmin,&davg,&navg);
        }
 
        // Move each particle
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );		

        if( find_option( argc, argv, "-no" ) == -1 )
        {
            if (navg) {
                absavg +=  davg/navg;
                nabsavg++;
            }
            if (dmin < absmin) 
                absmin = dmin;
            // Prints particle positions to a file based on SAVEFREQ
            if( fsave && (step % SAVEFREQ) == 0 )
                save( fsave, n, particles );
        }
    }
    // Simulation concluded, record running time and print
    simulation_time = read_timer( ) - simulation_time;
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    // If -no flag is not set, simulation correctness will be checked
    if( find_option( argc, argv, "-no" ) == -1 )
    {
        if (nabsavg) 
            absavg /= nabsavg;
        // Correctness check for simulation:
        //  - The minimum distance absmin between 2 particles during the run of the simulation
        //  - A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
        //  - A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
        //  - The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
        printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
        if (absmin < 0.4)
            printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
        if (absavg < 0.8)
            printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    // Print Summary Time
    if(fsum) 
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
    // Cleanup / Memory Management
    if(fsum)
        fclose( fsum );    
    free( particles );
    if(fsave)
        fclose( fsave );
    
    return 0;
}