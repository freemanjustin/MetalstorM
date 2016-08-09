#include "MetalstorM.h"

void initialize(e *E){

  // read input parameter file
  // for now I am just hardcoding the things we need 
  // to bootstrap MetalstorM into existence
  E->dt = 6.0;
  E->ntimes = 12;
  // set the rest of the stuff we need here

  // malloc memory so we can do some real work
  malloc_arrays(E); 

  // Initialize nonlinear model state variables
  init(E);	// jNOTE: fix the name of this function to be something more meaningful...

  // set the total simulation run time in seconds
  E->run_time=E->dt*E->ntimes;

}
