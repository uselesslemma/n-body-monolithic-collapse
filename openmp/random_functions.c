#include <stdlib.h>



// Useful function to allocate array
int  *generateIntArray( int N_elements ) 
{
  return (int *) calloc( N_elements, sizeof(int) ) ;
}


double randVal( double low, double high )
{
  return low + ( high - low ) * rand() / ( (double) RAND_MAX ) ;
}
