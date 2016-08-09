#include "MetalstorM.h"

void run(e *E){

  int t;

  for(t=0;t<E->ntimes;t++){
    // determine the current time
    get_time_info(E); //Set time indices and time clock --> see sub time_string()

    // read in required data from the input netcdf files
    get_data(E);

    // if required process the input data and time interpolare between the data
    // snapshots
    set_data(E);

    // initialize the free surface
    ini_zeta(E);

    // initialize the other fields
    ini_fields(E);

    // compute diagnostics
    // diag(E); <-- not doing this yet

    // Set vertical boundary conditions
    set_vbc(E);

    // write out fields to netcdf file
    output(E); //

    // Solve the vertically integrated primitive equations for the
    // free-surface and momentum components.
    //
    // Set time indices for predictor step. The PREDICTOR_2D_STEP switch
    // it is assumed to be false before the first time-step.
    
    // TODO...insert some code from main2d here 193-207

    // Predictor step - Advance barotropic equations using 2D time-step
    step2d(E);

    // Set time indices for corrector step.
    // TODO ... insert some code from main2d 220-228

    // Corrector step - Apply 2D time-step corrector scheme.
    step2d(E);

  }
}
