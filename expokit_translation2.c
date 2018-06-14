//
//  expokit_translation2.c
//  linkedCellMD
//
//  Created by Matthias Sachs on 28/04/2018.
//  Copyright © 2018 Matthias Sachs. All rights reserved.
//

#include "expokit_translation2.h"

#define my_max(x,y) ((x) >= (y)) ? (x) : (y)

int msign(double d){
    return d<-GSL_DBL_EPSILON?-1:d>GSL_DBL_EPSILON;
}
/*
double dot_spmatrixcol_vec(const gsl_spmatrix *A, const size_t col_index, const gsl_vector *p){
    for (size_t k; k<A->nz && A->i; k++) {
        <#statements#>
    }
}
 */

int test_non_zero(const gsl_vector *v, double tol){
    int counter = 0;
    for (int i =0; i < v->size; i++) {
        if(fabs(v->data[i]) > tol){
            counter++;
        }
    }
    return counter;
}
void printMatrix(gsl_matrix *A){
    /*char str[80];
    strcpy(str, "\n ");
    strcat(str, name);
    strcat(str, " = \n");
    printf("%s", str);*/
    //printf("\n");
    for (int i = 0; i < A->size1; i++){
        for (int j = 0; j < A->size2; j++)
            printf( "%f ", gsl_matrix_get(A, i, j));
        printf("\n");
    }
    
}

void printspMatrix(const gsl_spmatrix *A){
    /*char str[80];
     strcpy(str, "\n ");
     strcat(str, name);
     strcat(str, " = \n");
     printf("%s", str);*/
    //printf("\n");
    gsl_matrix *A_dens = gsl_matrix_alloc(A->size1, A->size2);
    gsl_spmatrix_sp2d(A_dens, A);
    printMatrix(A_dens);
    gsl_matrix_free(A_dens);
    
}
void printVector(gsl_vector *v){
    /*char str[80];
     strcpy(str, "\n ");
     strcat(str, name);
     strcat(str, " = \n");
     printf("%s", str);*/
    //printf("\n");
    for (int i = 0; i < v->size; i++){
        printf( "%f ", v->data[i]);
        
    }
    printf("\n");
}
double spmatrix_infnorm(const gsl_spmatrix *A){
    double colsum = 0.0;
    double colsum_max = 0.0;
    size_t i = A->i[0];
    for(size_t k=0; k<A->nz; k++){
        if (i == A->i[k] ) {
            colsum += fabs(A->data[k]);
        }
        else{
            if (colsum > colsum_max) {
                colsum_max = colsum;
            }
            i = A->i[k];
            colsum = fabs(A->data[k]);
        }
    }
    return colsum_max;
}

void dgexpv2( double t, const gsl_spmatrix *A, const gsl_vector *v, double tol, size_t m, gsl_vector *w, double *hump){
    if (A->size1 != A->size2) {
        printf("Matrix A must be square matrix");
    }
    if (test_non_zero(v, 0.001) > 0) {
        size_t n = A->size1;
        double anorm = spmatrix_infnorm(A);
        int mxrej = 10;  double btol  = 1.0e-7;
        double gamma = 0.9; double delta = 1.2;
        size_t mb    = m; double t_out   = fabs(t);
        size_t nstep = 0; double t_new   = 0;
        double t_now = 0; double s_error = 0;
        double rndoff= anorm*GSL_DBL_EPSILON;
        
        double err_loc, avnorm; size_t mx; //added when translating
        
        size_t k1 = 2; double xm = 1.0/m; double normv = gsl_blas_dnrm2(v); double beta = normv;
        double fact = pow(((m+1)/exp(1.0)),(m+1))*sqrt(2*M_PI*(m+1));
        t_new = (1/anorm)*pow((fact*tol)/(4*beta*anorm),xm);
        double s = pow(10,floor(log10(t_new))-1); t_new = ceil(t_new/s)*s;
        double sgn = msign(t);
        nstep = 0;
        gsl_blas_dcopy (v, w); // w= v
        *hump = normv;
        gsl_matrix *V = gsl_matrix_calloc(n, m+1);
        gsl_matrix *H = gsl_matrix_calloc (m+2, m+2);
        gsl_matrix *Hscaled = gsl_matrix_calloc (m+2, m+2);
        gsl_matrix *F = gsl_matrix_calloc (m+2, m+2);
        gsl_matrix_view Fview, Vview, Hview, Hviewscaled;
        gsl_vector_view Fcol;
        gsl_vector *p = gsl_vector_calloc(n);
        while (t_now < t_out){
            nstep = nstep + 1;
            double t_step = fmin( t_out-t_now,t_new );
            gsl_matrix_set_zero(V);
            gsl_matrix_set_zero(H);
            gsl_matrix_set_col(V, 0, w);
            gsl_matrix_scale(V, 1.0/beta);//V(:,1) = (1/beta)*w;
            //printf("\n V = \n");
            //printMatrix(V); OK
            gsl_vector_view colVj;
            gsl_vector_view colVi;
            double result;
            double s;
            for (size_t j=0; j<m; j++) {
                
                colVj = gsl_matrix_column (V, j);
                //printf("\n V = \n");
                //printMatrix(V);
                gsl_spblas_dgemv(CblasNoTrans, 1.0, A, &colVj.vector, 0.0, p); //p = A*V(:,j);
                //printf("\n p = \n");
                //printVector(p);
                for (int i=0; i<=j; i++) {
                    colVi = gsl_matrix_column (V, i);
                    gsl_blas_ddot(&colVi.vector, p, &result);
                    gsl_matrix_set(H, i, j, result); //H(i,j) = V(:,i)'*p;
                    gsl_blas_daxpy(-result, &colVi.vector, p); //p = p-H(i,j)*V(:,i);
                    //printf("\n p = \n");
                    //printVector(p);
                }
                //printf("\n H = \n");
                //printMatrix(H);
                s =  gsl_blas_dnrm2(p); //s = norm(p);
                if (s < btol) {
                    k1 = 0;
                    mb = j+1;
                    t_step = t_out-t_now;
                    break;
                }
                gsl_matrix_set(H, j+1, j, s);//H(j+1,j) = s;
                gsl_vector_scale(p, 1.0/s); // p not used again before reassignment. So scaling ok
                gsl_matrix_set_col(V, j+1, p); // V(:,j+1) = (1/s)*p;
            }
            //printf("\n H = \n");
            //printMatrix(H);
            if (k1 != 0) {
                gsl_matrix_set(H, m+1, m, 1.0);//H(m+2,m+1) = 1;
                //printf("\n H = \n");
                //printMatrix(H);
                colVj = gsl_matrix_column (V, m);
                gsl_spblas_dgemv(CblasNoTrans, 1.0, A, &colVj.vector, 0.0, p);
                avnorm = gsl_blas_dnrm2(p); //avnorm = norm(A*V(:,m+1));
            }
            //OK
            int ireject = 0;
            while (ireject <= mxrej) {
                mx = mb + k1;
                Hview = gsl_matrix_submatrix (H, 0, 0, mx, mx);
                Hviewscaled =gsl_matrix_submatrix (Hscaled, 0, 0, mx, mx);
                gsl_matrix_memcpy (&Hviewscaled.matrix, &Hview.matrix);
                gsl_matrix_scale(&Hviewscaled.matrix, sgn*t_step);
                Fview = gsl_matrix_submatrix (F, 0, 0, mx, mx);
                gsl_linalg_exponential_ss(&Hviewscaled.matrix, &Fview.matrix, .00001);//F = expm(sgn*t_step*H(1:mx,1:mx));
                //printf("\n Fview = \n");
                //printMatrix(&Fview.matrix);
                if (k1 == 0) {
                    err_loc = btol;
                    break;
                }
                else{
                    //printf("Fview=\n");
                    //printMatrix(&Fview.matrix);
                    double phi1 = fabs( beta* gsl_matrix_get(&Fview.matrix, m, 0));//phi1 = abs( beta*F(m+1,1) );
                    double phi2 = fabs( beta*avnorm * gsl_matrix_get(&Fview.matrix, m+1,0));//phi2 = abs( beta*F(m+2,1) * avnorm );
                    if (phi1 > 10*phi2) {
                        err_loc = phi2;
                        xm = 1/m;
                    }
                    else{
                        if (phi1 > phi2) {
                            err_loc = (phi1*phi2)/(phi1-phi2);
                            xm = 1/m;
                        }
                        else{
                            err_loc = phi1;
                            xm = 1/(m-1);
                        }
                    }
                }
                if (err_loc <= delta * t_step*tol) {
                    break;
                }
                else{
                    t_step = gamma * t_step * pow(t_step*tol/err_loc,xm);//t_step = gamma * t_step * (t_step*tol/err_loc)^xm;
                    s = pow(10,(floor(log10(t_step))-1));//s = 10^(floor(log10(t_step))-1);
                    t_step = ceil(t_step/s) * s; //t_step = ceil(t_step/s) * s;
                
                    if (ireject == mxrej) {
                        printf("The requested tolerance is too high.\n");
                        printf("Matrix A:\n");
                        printspMatrix(A);                    
                        printf("Vector v:\n");
                        printVector(v);
                    }
                    ireject = ireject + 1;
                }
            }
            if ( (int)k1-1 > 0)
                mx = mb + k1-1;
            else
                mx = mb;
            Fview = gsl_matrix_submatrix (F, 0, 0, mx, mx);
            //printf("Fview=\n");
            //printMatrix(&Fview.matrix);
            Fcol = gsl_matrix_column(&Fview.matrix, 0);
            Vview =gsl_matrix_submatrix(V, 0,0, V->size1, mx);
            gsl_blas_dgemv(CblasNoTrans, beta, &Vview.matrix, &Fcol.vector, 0.0, w); //w = V(:,1:mx)*(beta*F(1:mx,1));
            beta = gsl_blas_dnrm2(w);
            *hump = my_max(*hump,beta);
            
            t_now = t_now + t_step;
            t_new = gamma * t_step * pow(t_step*tol/err_loc, xm);
            s = pow(10, (floor(log10(t_new))-1));
            t_new = ceil(t_new/s) * s;
            err_loc = my_max(err_loc,rndoff);
            s_error = s_error + err_loc;

        }
        gsl_matrix_free(V);
        gsl_matrix_free(H);
        gsl_matrix_free(Hscaled);
        gsl_matrix_free(F);
        gsl_vector_free(p);
    }
    else{
        gsl_vector_set_zero(w);
    }
}
