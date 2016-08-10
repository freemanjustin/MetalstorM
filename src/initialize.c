#include "MetalstorM.h"

void initialize(e *E){

  // read input parameter file
  // for now I am just hardcoding the things we need
  // to bootstrap MetalstorM into existence
  E->Lm = 142;
  E->Mm = 115;
  E->dt = 6.0;
  E->ntimes = 12;
  E->rdrg2 = 1.0e-3;
  E->dcrit = 0.1;
  E->dstart =  15006.0;  // days
  E->time_ref = 19700101;
  E->grid_name = malloc(strlen("tiny.nc")+1);
  strcpy(E->grid_name, "tiny.nc");
  E->nffiles = 2;
  E->frcname = (char**)malloc(2);
  E->frcname[0] = (char*)malloc(strlen("pair_0001.nc")+1);
  strcpy(E->frcname[0], "pair_0001.nc");
  E->frcname[1] = (char*)malloc(strlen("strs_0001.nc")+1);
  strcpy(E->frcname[1], "strs_0001.nc");

  // set the rest of the stuff we need here

  // malloc memory so we can do some real work
  malloc_arrays(E);

  // Initialize nonlinear model state variables
  //init(E);	// jNOTE: fix the name of this function to be something more meaningful...

  // set the total simulation run time in seconds
  E->run_time=E->dt*E->ntimes;

}
