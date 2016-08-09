#include <stdio.h>
#include <stdlib.b>
#include <math.h>

#define	true 0
#define false 1

#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))

#define fail(...) Mstorm_fail(__LINE__,__func__,__FILE__,__VA_ARGS__)

typedef struct{


}e;

// prototypes

// fail.c
void Mstorm_fail( const int, const char*, const char*, const char*, ... );
