#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "jutil.h"

#define	true 0
#define false 1

#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))

#define fail(...) Mstorm_fail(__LINE__,__func__,__FILE__,__VA_ARGS__)

/*
typedef struct{
    int Lm; // Number of interior grid points in the XI-direction
    int Mm; // Number of internal grid points in the ETA-direction
    int LBi; // lower bound I-dimension
    int UBi; // upper bound I-dimension
    int LBj; // lower bound J-dimension
    int UBj; // upper bound J-dimension
    int LBij; // lower bound MIN(I,J)-dimension
    int UBij; // upper bound MAX(I,J)-dimension

    int edge[4][4]; // boundary edges I- or J-indices

    int Istr; // starting tile I-direction
    int Iend;  // ending   tile I-direction
    int Jstr;  // starting tile J-direction
    int Jend;  // ending   tile J-direction
    int IstrR; // starting tile I-direction (RHO)
    int IendR; // ending   tile I-direction (RHO)
    int IstrU; // starting tile I-direction (U)
    int JstrR; // starting tile J-direction (RHO)
    int JendR; // ending   tile J-direction (RHO)
    int JstrV; // starting tile J-direction (V)
    int IstrB; // starting obc I-direction (RHO,V)
    int IendB; // ending   obc I-direction (RHO,V)
    int IstrM; // starting obc I-direction (PSI,U)
    int JstrB; // starting obc J-direction (RHO,U)
    int JendB; // ending   obc J-direction (RHO,U)
    int JstrM; // starting obc J-direction (PSI,V)

    int Istrm3;    // starting I-halo, Istr-3
    int Istrm2;    // starting I-halo, Istr-2
    int Istrm1;    // starting I-halo, Istr-1
    int IstrUm2;   // starting I-halo, IstrU-2
    int IstrUm1;   // starting I-halo, IstrU-1
    int Iendp1;    // ending   I-halo, Iend+1
    int Iendp2;    // ending   I-halo, Iend+2
    int Iendp2i;   // ending   I-halo, Iend+2
    int Iendp3;    // ending   I-halo, Iend+3
    int Jstrm3;    // starting J-halo, Jstr-3
    int Jstrm2;    // starting J-halo, Jstr-2
    int Jstrm1;    // starting J-halo, Jstr-1
    int JstrVm2;   // starting J-halo, JstrV-2
    int JstrVm1;   // starting J-halo, JstrV-1
    int Jendp1;    // ending   J-halo, Jend+1
    int Jendp2;    // ending   J-halo, Jend+2
    int Jendp2i;   // ending   J-halo, Jend+2
    int Jendp3;    // ending   J-halo, Jend+3

    // ?? what are these??
    //int Imin(:,:,:)  // starting ghost I-direction
    //int Imax(:,:,:)  // ending   ghost I-direction
    //int Jmin(:,:,:)  // starting ghost J-direction
    //int Jmax(:,:,:)  // ending   ghost J-direction

}T_BOUNDS;
*/

typedef struct{

    int time_ref;

    int nffiles;



    char    **frcname;
    char    *grid_name;


    float   ntimes;
    float   dt;
    float   rdrg2;
    float   dcrit;

    float   dstart;

    float   run_time;

}e;

// prototypes

void initialize( e* );
void malloc_arrays( e* );

// fail.c
void Mstorm_fail( const int, const char*, const char*, const char*, ... );
