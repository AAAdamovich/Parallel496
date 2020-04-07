/*  Antony Adamovich, Dimitra Deliopoulos, Mahmoud Gudarzi
*    serial.cpp for XSEDE hw2 "Particle Simulator" - Serial version
*    West Chester University - CSC 496 - Dr. Lihn B. Ngo
*    Created: 04-MAR-2020 - Last Edited: 06-APR-2020 by Antony Adamovich
*    Please see https://github.com/AAAdamovich/Parallel496
*       for version tracking
*    Definitions:
*       Field: The entire square area where particles will be interacting, 
*       the total simulation area
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

enum bins {C, N, NE, E, SE, S, SW, W, NW, NUM_BINS};

int get_square(particle_t &particle, int divisions, double subsize){
    int xcoordinate = (int)(particle.x / subsize);
    int ycoordinate = (int)(particle.y / subsize);
    int location = ycoordinate + (divisions * xcoordinate);
    // Opportunity for optimization for near-wall cases
    if(particle.x <= (subsize * (double)xcoordinate) + CUTOFF){
        // Particle is in bin NW, N or NE
        if(particle.y <= (subsize * (double)ycoordinate) + CUTOFF){
            // In bin NW
            return (location * NUM_BINS) + NW;
        }
        if(particle.y >= (subsize * (double)(ycoordinate + 1)) - CUTOFF){
            // In bin NE
            return (location * NUM_BINS) + NE;
        }
        // Otherwise, particle in bin N
        return (location * NUM_BINS) + N;
    }
    // Particle is in bin W, C, E, SW, S, or SE
    if(particle.x >= (subsize * (double)(xcoordinate + 1)) - CUTOFF){
        // Particle is in bin SW, S, or SE
        if(particle.y <= (subsize * (double)ycoordinate) + CUTOFF){
            // In bin SW
            return (location * NUM_BINS) + SW;
        }
        if(particle.y >= (subsize * (double)(ycoordinate + 1)) - CUTOFF){
            // In bin SE
            return (location * NUM_BINS) + SE;
        }
        // Otherwise, particle in bin S
        return (location * NUM_BINS) + S;
    }
    // Particle is in bin W, C, or E
    if(particle.y <= (subsize * (double)ycoordinate) + CUTOFF){
        // In bin W
        return (location * NUM_BINS) + W;
    }
    if(particle.y >= (subsize * (double)(ycoordinate + 1)) - CUTOFF){
        // In bin E
        return (location * NUM_BINS) + E;
    }
    // Options exhausted, particle in bin C
    return (location * NUM_BINS) + C;
}


int main( int argc, char **argv )
{    
    
    // !! Testing Value !! Must be altered dynamically, NOT YET IMPLEMENTED
    int divisions = 15;
    // Corresponding to a MxM matrix where divisions = M
    
    /* ==== BIN SETUP ====
    *  
    *  NW - 8  |  N - 1  |  NE - 2
    * ---------|---------|---------
    *   W - 7  |  C - 0  |  E - 3
    * ---------|---------|---------
    *  SW - 6  |  S - 5  |  SE - 4
    *
    * Hence, NINE (9) is the total number of bins as NUM_BINS = 9
    */
    


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
    // Set divisions dynamically
    divisions = (int) sqrt(n);
    
    subsize = size / ((double)divisions);
    init_particles( n, particles );
    

    // 2D master field array that contains all particles
    std::vector< std::vector<particle_t> > field (divisions * divisions * NUM_BINS, std::vector<particle_t>());

    // Time when simulation begins
    double simulation_time = read_timer( );
    
	// "Sort" particles into master field array
    for(int i = 0; i < n; i++){
        field[get_square(particles[i], divisions, subsize)].push_back(particles[i]);
    }
    
    // DEBUGGING PRINTS
    /*
      printf("subsize: %lf  subsize*2: %lf\n", subsize, subsize * 2.0);
      for(int i = 0; i < field.size(); i++){
          if(field[i].size() == 0){
             printf("i: %d  EMPTY\n", i);
          }
          for(int j = 0; j < field[i].size(); j++){
              printf("i: %d  j:  %d  X: %lf  Y: %lf\n", i, j, field[i][j].x, field[i][j].y);
          }
       }
    */
    // How long the simulation is to run is based on NSTEPS
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        davg = 0.0;
        dmin = 1.0;
        
        // === PARTICLE INTERACTIONS WITHIN SAME CELL ===
        
        // Go by cells, cell % NUM_BINS would be its "location" in the grid
        for(int cell = 0; cell < field.size(); cell += NUM_BINS){
            // Go by the bins within cells
            for(int bin = cell; bin < cell + NUM_BINS; bin++){
                // Now we start looking at particles within cells
                for(int i = 0; i < field[bin].size(); i++){
                    // Reset acceleration values in prep for calculation
                    (field[bin][i]).ax = 0;
                    (field[bin][i]).ay = 0;
                    
                    // Particle interactions within the same bin
                    for(int j = 0; j < i; j++){
                        apply_force( field[bin][i], field[bin][j],&dmin,&davg,&navg);
                    }
                    //printf("cell: %d  bin:  %d  i: %d  XAccel: %lf\n", cell, bin, i, (field[bin][i]).ax);
                    
                    // Particle interactions within the same cell, different bins
                    for(int localbin = bin + 1; localbin < cell + NUM_BINS; localbin++){
                        for(int j = 0; j < field[localbin].size(); j++){
                            apply_force( field[bin][i], field[localbin][j],&dmin,&davg,&navg);
                        }
                    }
                    //printf("cell: %d  bin:  %d  i: %d  XAccel: %lf\n", cell, bin, i, (field[bin][i]).ax);
                }
            }
        }
        
        // TODO === PARTICLE INTERACTIONS ALONG CELL BORDERS ===
        
        // Go by cells, cell % NUM_BINS would be its "location" in the grid
        for(int cell = 0; cell < field.size(); cell += NUM_BINS){
            // Go by the bins within cells, skip cell 0 or 'C'
            for(int bin = cell + 1; bin < cell + NUM_BINS; bin++){
                switch(bin % NUM_BINS){
                    case N:
                        for(int i = 0; i < field[bin].size(); i++){
                            if((cell - (divisions * NUM_BINS)) + SE >= 0){
                                for(int j = 0; j < field[(cell - (divisions * NUM_BINS)) + SE].size(); j++){
                                    apply_force( field[bin][i], field[(cell - (divisions * NUM_BINS)) + SE][j],&dmin,&davg,&navg);
                                    //printf("1\n");
                                }
                            }
                            if((cell - (divisions * NUM_BINS)) + S >= 0){
                                for(int j = 0; j < field[(cell - (divisions * NUM_BINS)) + S].size(); j++){
                                    apply_force( field[bin][i], field[(cell - (divisions * NUM_BINS)) + S][j],&dmin,&davg,&navg);
                                    //printf("2\n");
                                }
                            }
                            if((cell - (divisions * NUM_BINS)) + SW >= 0){
                                for(int j = 0; j < field[(cell - (divisions * NUM_BINS)) + SW].size(); j++){
                                    apply_force( field[bin][i], field[(cell - (divisions * NUM_BINS)) + SW][j],&dmin,&davg,&navg);
                                    //printf("3\n");
                                }
                            }
                        }
                        break;
                        
                    case NE:
                        for(int i = 0; i < field[bin].size(); i++){
                            if((cell - (divisions * NUM_BINS)) + SE >= 0){
                                for(int j = 0; j < field[(cell - (divisions * NUM_BINS)) + SE].size(); j++){
                                    apply_force( field[bin][i], field[(cell - (divisions * NUM_BINS)) + SE][j],&dmin,&davg,&navg);
                                    //printf("4\n");
                                }
                            }
                            if((cell - (divisions * NUM_BINS)) + S >= 0){
                                for(int j = 0; j < field[(cell - (divisions * NUM_BINS)) + S].size(); j++) {
                                    apply_force( field[bin][i], field[(cell - (divisions * NUM_BINS)) + S][j],&dmin,&davg,&navg);
                                    //printf("5\n");
                            }}
                            if((cell - (divisions * NUM_BINS)) + NUM_BINS + SW >= 0 && ((cell / NUM_BINS) % divisions != 2)){
                                for(int j = 0; j < field[(cell - (divisions * NUM_BINS)) + NUM_BINS + SW].size(); j++) {
                                    apply_force( field[bin][i], field[(cell - (divisions * NUM_BINS)) + NUM_BINS + SW][j],&dmin,&davg,&navg);
                                    //printf("6\n");
                            }}
                            if((cell / NUM_BINS) % divisions != (divisions - 1)){
                                for(int j = 0; j < field[cell + NUM_BINS + NW].size(); j++) {
                                    apply_force( field[bin][i], field[cell + NUM_BINS + NW][j],&dmin,&davg,&navg);
                                    //printf("7\n");
                            }}
                            if((cell / NUM_BINS) % divisions != (divisions - 1)){
                                for(int j = 0; j < field[cell + NUM_BINS + W].size(); j++){
                                    apply_force( field[bin][i], field[cell + NUM_BINS + W][j],&dmin,&davg,&navg);
                                    //printf("8\n");
                                }
                            }
                        }
                        break;
                        
                    case E:
                        for(int i = 0; i < field[bin].size(); i++){
                            if((cell / NUM_BINS) % divisions != (divisions - 1)){
                                for(int j = 0; j < field[cell + NUM_BINS + NW].size(); j++) {
                                    apply_force( field[bin][i], field[cell + NUM_BINS + NW][j],&dmin,&davg,&navg);
                                    //printf("9\n");
                            }}
                            if((cell / NUM_BINS) % divisions != (divisions - 1)){
                                for(int j = 0; j < field[cell + NUM_BINS + W].size(); j++) {
                                    apply_force( field[bin][i], field[cell + NUM_BINS + W][j],&dmin,&davg,&navg);
                                    //printf("10\n");
                            }}
                            if((cell / NUM_BINS) % divisions != (divisions - 1)){
                                for(int j = 0; j < field[cell + NUM_BINS + SW].size(); j++) {
                                    apply_force( field[bin][i], field[cell + NUM_BINS + SW][j],&dmin,&davg,&navg);
                                    //printf("11\n");
                            }}
                        }
                        break;
                        
                        
                    case SE:
                        for(int i = 0; i < field[bin].size(); i++){
                            
                            if((cell / NUM_BINS) % divisions != (divisions - 1)){
                                for(int j = 0; j < field[cell + NUM_BINS + W].size(); j++) {
                                    apply_force( field[bin][i], field[cell + NUM_BINS + W][j],&dmin,&davg,&navg);
                                //printf("12\n");
                            }}
                            if((cell / NUM_BINS) % divisions != (divisions - 1)){
                                for(int j = 0; j < field[cell + NUM_BINS + SW].size(); j++) {
                                    apply_force( field[bin][i], field[cell + NUM_BINS + SW][j],&dmin,&davg,&navg);
                                //printf("13\n");
                            }}
                            
                            if((cell / NUM_BINS) % divisions != (divisions - 1) && (cell + (divisions * NUM_BINS)) + NUM_BINS < field.size()){
                                for(int j = 0; j < field[(cell + (divisions * NUM_BINS)) + NUM_BINS + NW].size(); j++) {
                                    apply_force( field[bin][i], field[(cell + (divisions * NUM_BINS)) + NUM_BINS + NW][j],&dmin,&davg,&navg);
                                //printf("14\n");
                        }}
                            
                            if(cell + (divisions * NUM_BINS) < field.size()){
                                for(int j = 0; j < field[(cell + (divisions * NUM_BINS)) + N].size(); j++) {
                                    apply_force( field[bin][i], field[(cell + (divisions * NUM_BINS)) + N][j],&dmin,&davg,&navg);
                                //printf("15\n");
                            }}
                            if(cell + (divisions * NUM_BINS) < field.size()){
                                for(int j = 0; j < field[(cell + (divisions * NUM_BINS)) + NE].size(); j++) {
                                    apply_force( field[bin][i], field[(cell + (divisions * NUM_BINS)) + NE][j],&dmin,&davg,&navg);
                                //printf("16\n");
                            }}
                        }
                    break;
                    
                    case S:
                    for(int i = 0; i < field[bin].size(); i++){
                        if(cell + (divisions * NUM_BINS) < field.size()){
                            for(int j = 0; j < field[(cell + (divisions * NUM_BINS)) + NW].size(); j++) 
                                apply_force( field[bin][i], field[(cell + (divisions * NUM_BINS)) + NW][j],&dmin,&davg,&navg);
                        }
                        if(cell + (divisions * NUM_BINS) < field.size()){
                            for(int j = 0; j < field[(cell + (divisions * NUM_BINS)) + N].size(); j++) 
                                apply_force( field[bin][i], field[(cell + (divisions * NUM_BINS)) + N][j],&dmin,&davg,&navg);
                        }
                        if(cell + (divisions * NUM_BINS) < field.size()){
                            for(int j = 0; j < field[(cell + (divisions * NUM_BINS)) + NE].size(); j++) 
                                apply_force( field[bin][i], field[(cell + (divisions * NUM_BINS)) + NE][j],&dmin,&davg,&navg);
                        }
                    }
                    break;
                    
                    case SW:
                    for(int i = 0; i < field[bin].size(); i++){
                        
                        if((cell / NUM_BINS) % divisions != 0){
                            for(int j = 0; j < field[cell - NUM_BINS + E].size(); j++) 
                                apply_force( field[bin][i], field[cell - NUM_BINS + E][j],&dmin,&davg,&navg);
                        }
                        if((cell / NUM_BINS) % divisions != 0){
                            for(int j = 0; j < field[cell - NUM_BINS + SE].size(); j++) 
                                apply_force( field[bin][i], field[cell - NUM_BINS + SE][j],&dmin,&davg,&navg);
                        }
                        if((cell / NUM_BINS) % divisions != 0 && cell + (divisions * NUM_BINS) - NUM_BINS < field.size()){
                            for(int j = 0; j < field[cell + (divisions * NUM_BINS) - NUM_BINS + NE].size(); j++) 
                                apply_force( field[bin][i], field[cell + (divisions * NUM_BINS) - NUM_BINS + NE][j],&dmin,&davg,&navg);
                        }
                        if(cell + (divisions * NUM_BINS) < field.size()){
                            for(int j = 0; j < field[(cell + (divisions * NUM_BINS)) + N].size(); j++) 
                                apply_force( field[bin][i], field[(cell + (divisions * NUM_BINS)) + N][j],&dmin,&davg,&navg);
                        }
                        if(cell + (divisions * NUM_BINS) < field.size()){
                            for(int j = 0; j < field[(cell + (divisions * NUM_BINS)) + NW].size(); j++) 
                                apply_force( field[bin][i], field[(cell + (divisions * NUM_BINS)) + NW][j],&dmin,&davg,&navg);
                        }
                        
                    }
                    break;
                    case W:
                    for(int i = 0; i < field[bin].size(); i++){
                        
                        if((cell / NUM_BINS) % divisions != 0){
                            for(int j = 0; j < field[cell - NUM_BINS + NE].size(); j++) 
                                apply_force( field[bin][i], field[cell - NUM_BINS + NE][j],&dmin,&davg,&navg);
                        }
                        
                        if((cell / NUM_BINS) % divisions != 0){
                            for(int j = 0; j < field[cell - NUM_BINS + E].size(); j++) 
                                apply_force( field[bin][i], field[cell - NUM_BINS + E][j],&dmin,&davg,&navg);
                        }
                        if((cell / NUM_BINS) % divisions != 0){
                            for(int j = 0; j < field[cell - NUM_BINS + SE].size(); j++) 
                                apply_force( field[bin][i], field[cell - NUM_BINS + SE][j],&dmin,&davg,&navg);
                        }
                        
                    }
                    break;
                    case NW:
                    for(int i = 0; i < field[bin].size(); i++){
                            if((cell - (divisions * NUM_BINS)) + SW >= 0){
                                for(int j = 0; j < field[(cell - (divisions * NUM_BINS)) + SW].size(); j++){
                                    apply_force( field[bin][i], field[(cell - (divisions * NUM_BINS)) + SW][j],&dmin,&davg,&navg);
                                    //printf("Force applied in %d\n", cell);
                                }
                            }
                            if((cell - (divisions * NUM_BINS)) + S >= 0){
                                for(int j = 0; j < field[(cell - (divisions * NUM_BINS)) + S].size(); j++) {
                                    apply_force( field[bin][i], field[(cell - (divisions * NUM_BINS)) + S][j],&dmin,&davg,&navg);
                                    //printf("Force applied in %d\n", cell);
                            }}
                            
                            if((cell - (divisions * NUM_BINS)) - NUM_BINS + SE >= 0 && (cell / NUM_BINS) % divisions != 0){
                                for(int j = 0; j < field[(cell - (divisions * NUM_BINS)) - NUM_BINS + SE].size(); j++) {
                                    apply_force( field[bin][i], field[(cell - (divisions * NUM_BINS)) - NUM_BINS + SE][j],&dmin,&davg,&navg);
                                    //printf("Force applied in %d\n", cell);
                            }}
                            if((cell / NUM_BINS) % divisions != 0){
                                for(int j = 0; j < field[cell - NUM_BINS + NE].size(); j++) {
                                    apply_force( field[bin][i], field[cell - NUM_BINS + NE][j],&dmin,&davg,&navg);
                                    //printf("Force applied in %d\n", cell);
                            }}
                            if((cell / NUM_BINS) % divisions != 0){
                                for(int j = 0; j < field[cell - NUM_BINS + E].size(); j++){
                                    apply_force( field[bin][i], field[cell - NUM_BINS + E][j],&dmin,&davg,&navg);
                                    //printf("Force applied in %d\n", cell);
                                }
                            }
                        }
                    break;
                }    
            }
        }
               
            
        
        
        
        
        // Starter Code
    /*    
          for( int i = 0; i < n; i++ )
          {
            particles[i].ax = particles[i].ay = 0;
              for (int j = 0; j < n; j++ ){
				  apply_force( particles[i], particles[j],&dmin,&davg,&navg);
                }
          }
    */   

        // Particle movement alt.
        int my_index;
        for(int i = 0; i < field.size(); i++){
            for(int j = 0; j < field[i].size(); j++){
                //printf("AccelX: %lf  AccelY: %lf\n", field[i][j].ax, field[i][j].ay);
                move(field[i][j]);
                my_index = get_square(field[i][j], divisions, subsize);
                if(my_index != i){
                    // Particle has moved outside of its bin into another one, must be updated
                    //printf("WHERE DID MY LOST LAMB GO?! i: %d myindex: %d\n", i, my_index);
                    // Add particle to its new home
                    field[my_index].push_back(field[i][j]);
                    // Remove duplicate particle from old location
                    (field[i]).erase(field[i].begin() + j);
                }
            }
        }
 
 /*
        // Move each particle
         for( int i = 0; i < n; i++ ){
             
             move( particles[i] );
         } 
            		
*/
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