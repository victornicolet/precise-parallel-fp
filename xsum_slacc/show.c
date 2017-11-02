/* SHOW A FLOATING POINT VALUE GIVEN AS AN ARGUMENT IN BINARY. */

#include "pbinary.h"
#include <stdlib.h>
#include <stdio.h>

int main (int argc, char **argv)
{
  int i;
  for (i = 1; i < argc; i++)
  { pbinary_double(atof(argv[i]));
    printf("\n");
  }
  return 0;
}
