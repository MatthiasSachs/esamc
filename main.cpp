//
//  main.cpp
//  esamc
//
//  Created by Matthias Sachs on 14/06/2018.
//  Copyright © 2018 Matthias Sachs. All rights reserved.
//

#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include "expokit_translation2.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include "test_exponential.hpp"
#include "test_DPD_mBAOAB.hpp"


/*
 Linked cell alogrithm in c++ for efficient integration of markovian and non-markovian DPD systems using BAOAB in combination with Krylov sub-space methods
 */


//#include "domains.hpp"

//#include "grid.hpp"

int main(int argc, const char * argv[]) {
    //test_Verlet();
    test_mBAOAB();
//     test_mBAOAB_RR3();
    //test_sparseExp();
}
