#include "MetalstorM.h"

int main(int argc,char **argv)
{
  // Everything
  e *E;

  E = malloc(sizeof(e));

  initialize(E);
  run(E);
  finalize(E);

  return 0;
}
