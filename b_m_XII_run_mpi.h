#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <sys/stat.h>
#include <mpi.h>

#include "JRand_dyn.h"
#include "model_kmc.h"
#include "kinetic_monte_carlo.h"
#include "bound_misbound_XII.h"

#define NUMBER_PARTICLES 1000
