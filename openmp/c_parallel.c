
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#include <omp.h>        // Use OpenMP

#include "my_include.h" // Some local functions I'll use


int main()
{


  /*

    ALWAYS SET THE NUMBER OF THREADS TO USE
    If not, your process will fight to use all of the available CPU. 
    Make sure you know how many the computer you are working on has.
    Generally, multiples of 2 are good numbers to use. Sometimes 
    it is forced to be a power of two, which is better than a multiple
    of two. Just don't use all, or even most if there is a remote 
    chance someone else will use the same computer.

    If unbounded, it will decrease the run time of every other process
    people are trying to run, and is the like #2 reason people
    will hate you during this project. Don't make me hate you.
    Multithread responsibly.

  */


  int N_threads = 4 ;

  omp_set_num_threads( N_threads ) ; // Function to set the number of threads to use






  // First let's look at no parallelization...

  for ( int i = 0 ; i < 10 ; ++i )
  {
    int thread_num  = omp_get_thread_num () ;
    int num_threads = omp_get_num_threads() ;

    printf("No parallel thread number is %i of %i\n", thread_num, num_threads );
  }
  printf("\n");



  // Then add some! I'll explain below

  #pragma omp parallel for
  for ( int i = 0 ; i < 10 ; ++i )
  {
    int thread_num  = omp_get_thread_num () ;
    int num_threads = omp_get_num_threads() ;
    if (omp_in_parallel())
    printf("   Parallel thread number is %i of %i\n", thread_num, num_threads );
  }
  printf("\n");

  /*
    0 is the master thread, and with no parallelization there is only 1 thread, so
    the first loop prints 0 and 1. When parallelized, you see threads 0 to OMP_NUM_THREADS-1. 
    Additionally, the program will wait until every parallel thread finishes before
    moving on in the code. Magical shit.


    #pragma omp parallel    is the basis of most of the parallelization you will be
    using. It sets the next region it encounters as parallel


    for indicates a for loop, which you will most likely be using as well.
    Some other useful options are:

    num_threads( var )           <- Set number of threads, or can set globally with omp_set_num_threads()

    default(shared)              <- Sets variables to default to being shared
    default(private)             <- Sets variables to default to being private

    private( var1, var2, ... )   <- Sets var1, var2, etc to private
    shared( var1, var2, ... )    <- Sets var1, var2, etc to shared



    So, what is shared and private? Yup, we are getting into riddles. Basically, 
    private means each thread gets it's own copy of a variable to play with, and
    shared allows everything to access and change the same one. Be very careful, as
    this is one of the easiest things to fuck up on. 
  */





  // Example of private/shared
  {
    int x = 0 ;
    int y = 0 ;


    // Lets try this out, guess what will happen:
    #pragma omp parallel for private( x ) shared ( y ) num_threads(N_threads) 
    for ( int i = 0 ; i < 10 ; ++i )
    {
      x = 4 + x ;
      y = 2 + y ;

      int z = i * x ;

      printf("i =%2i,  Thread:%3i,   Private ( x = 4 + x ):%3i,  Shared ( y = 2 + y ):%3i,  Local (z=i*x):%3i\n", 
             i , 
             omp_get_thread_num(),
             x , 
             y ,
             z );
    }
    printf("\n");



    /*
      Before you try to understand it, run it a few times. You'll notice numbers 
      are changing. Even if you think you know what is going on, guess what? 
      Parallelization is weird. You will notice for the same threads, the private 
      variable (x) is doing all sorts of crazy and random shit! This is because each 
      thread got it's own x, and in creating the variable for each thread it was 
      left uninitialized. 

      We can fix this by using firstprivate, which will treat x as private, but use 
      the value it had when the parallel region began
    */

    x = 0 ; // Reset these...
    y = 0 ;


    // With the fix, see if you can track whats being done in each thread to understand what's going on

    #pragma omp parallel for firstprivate( x ) shared ( y ) num_threads(N_threads) 
    for ( int i = 0 ; i < 10 ; ++i )
    {
      x = 4 + x ;
      y = 2 + y ;

      int z = i * x ;

      printf("i =%2i,  Thread:%3i,   Private ( x = 4 + x ):%3i,  Shared ( y = 2 + y ):%3i,  Local (z=i*x):%3i\n", 
             i , 
             omp_get_thread_num(),
             x , 
             y ,
             z );
    }
    printf("\n");



  }






  /*
    Let's do a fun example! Ever want to estimate pi? Yes you have.



    If you take the ratio of the area of a circle to square with equal diameter and length, you get

    pi R^2 / D^2 = pi (D/2)^2 / D^2 = pi/4 D^2/D^2 = pi/4




    If you thus throw a random distribution onto this area, the fraction that end up in pi is simply

    4 * c/T

    where c is the number that landed in the circle, and T is the total number thrown in. This is the
    kind of problem where parallelization is super useful, as we want to do a shitton of calculations
    that are independant of eachother

  */

  srand(time(NULL))         ; // Sets a random seed for psudo RNG based on execution time

  int  N_probes = 1000000   ; // Number of darts we are throwing

  int *N_circle = generateIntArray( N_threads ) ; // A nice work around for counting in parallelization
                                                  // Each thread can count it's amount individually
                                                  // Avoids the effects of shared variables





  #pragma omp parallel for shared(N_circle)       // We can share array, since each thread gets it's own index
  for ( int i = 0 ; i < N_probes ; ++i )
  {

    int thread_num = omp_get_thread_num() ;

    double x  = randVal( 0.0, 1.0 ) ; // Since these are defined inside the loop, no need for private
    double y  = randVal( 0.0, 1.0 ) ; // This is my own function, for info on random number generation, see google

    double r2 = x*x + y*y           ; // r^2, if it's in the circle, r^2 < 1


    // If our probe landed in the circle, count it
    if ( r2 < 1.0 ) 
      N_circle[ thread_num ] = N_circle[ thread_num ] + 1 ;
  }




  // Since we counted by thread, need to sum count in each thread for total

  double N_inCircle = 0 ;

  for ( int i = 0 ; i < N_threads ; ++i )
    N_inCircle += N_circle[i] ;



  // Don't need array anymore, free dat shit

  free( N_circle );




  // This answer will change every time due to how I set the seed, not due to parallelization

  printf("Pi is approximately: %10.6f\n", 4*N_inCircle/N_probes );





  return 0 ;
}
