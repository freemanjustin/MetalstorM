// gridr
//
// freeman.justin@gmail.com

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#ifndef __jutil_H__
#define __jutil_H__

typedef struct{
    float   x;
    float   y;
}pts;

#define pi  (M_PI)
#define earthradius (6356750.52)
#define deg2rad (M_PI/180.0)
#define rad2deg (180.0/M_PI)


/* function prototypes */

/* memory allocation functions */
/* deallocate memory with free() */

extern double**** malloc4d_double( int, int, int, int );

extern double*** malloc3d_double( int, int, int);

extern float*** malloc3d_float(int, int, int );

extern int ***malloc3d_int( int, int, int );

extern short int ***malloc3d_short( int, int, int );

extern signed char ***malloc3d_char(int, int, int );

extern double **malloc2d_double( int, int );

extern int **malloc2d_int( int, int );

extern float **malloc2d_float( int, int );

extern signed char **malloc2d_signedchar( int, int );

extern char **malloc2d_char( int, int );

extern short int **malloc2d_short( int, int );

extern char *lr_pack( char* );
extern char *strip_slash_n( char* );
extern void get_command_line_arg_as_string( char**, char* );
extern void get_command_line_arg_as_int( int*, char* );

/* math utils */
/* round_double is named to avoid potential namespace clashes with the c99 function round() */
/* c99 targets should use round() rather than this one */
extern double round_double( double );

extern double j_fmod(double, double);

extern void julday(int, int, int, int, int, double, double*);
extern void caldat(double, int*, int*, int*, int*, int*, double*);


int pnpoly(pts *polygon, int npol, pts p);
double spheriq_dist(double lon1, double lat1, double lon2, double lat2, int debug);
float point_to_line_distance(float x, float y, float x1, float y1, float x2, float y2);

#endif /*__jutil_H__ */
