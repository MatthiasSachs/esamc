//
//  test_DPD_mBAOAB.hpp
//  linkedCellMD
//
//  Created by Matthias Sachs on 10/06/2018.
//  Copyright © 2018 Matthias Sachs. All rights reserved.
//

#ifndef test_DPD_mBAOAB_hpp
#define test_DPD_mBAOAB_hpp

#include <stdio.h>
#include <iostream>
#include "grid.hpp"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include "expokit_translation2.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include "domains.hpp"
#include "potentials.hpp"
#include "outputshedulers.hpp"
#include "dynamics.hpp"
#include "langevin_DPD.hpp"
#include "langevin_DPD_literature.hpp"
void test_Verlet();
void test_mBAOAB();


#endif /* test_DPD_mBAOAB_hpp */
