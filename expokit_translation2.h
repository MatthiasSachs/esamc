//
//  expokit_translation2.h
//  linkedCellMD
// 
//  Created by Matthias Sachs on 28/04/2018.
//  Copyright © 2018 Matthias Sachs. All rights reserved.
//

#ifndef expokit_translation2_h
#define expokit_translation2_h

#ifdef __cplusplus
extern "C" {
#endif
#include <gsl/gsl_blas.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>
#include <math.h>
#include <stdio.h>
//#include "expohelp.h"
    
//void dgexpv2( double t, const gsl_spmatrix *A, const gsl_vector *v, double tol, size_t m, gsl_vector *w, int *err, double *hump);
void dgexpv2( double t, const gsl_spmatrix *A, const gsl_vector *v, double tol, size_t m, gsl_vector *w, double *hump);
void printMatrix(gsl_matrix *A);
void printspMatrix(const gsl_spmatrix *A);
void printVector(gsl_vector *v);
int test_non_zero(const gsl_vector *v, double tol);
#ifdef __cplusplus
}
#endif


#endif /* expokit_translation2_h */
