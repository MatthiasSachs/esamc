//
//  test_exponential.cpp
//  linkedCellMD
//
//  Created by Matthias Sachs on 31/05/2018.
//  Copyright © 2018 Matthias Sachs. All rights reserved.
//

#include "test_exponential.hpp"

void test_sparseExp(){
    // insert code here...
    /*
     std::cout << "Hello, World!\n";
     int Np= 10;
     int sdim = 2;
     double cutoff = 1.0;
     int Ni = 10;
     const gsl_rng_type * T;
     gsl_rng * r;
     gsl_rng_env_setup();
     T = gsl_rng_default;
     r = gsl_rng_alloc (T);
     
     
     double L = 5.0;
     gsl_vector *L_vec = gsl_vector_alloc(sdim);
     gsl_vector_set_all (L_vec, L);
     gsl_vector_int *nc_vec = gsl_vector_int_alloc(sdim);
     for (int d=0; d<L_vec->size; d++)
     {
     nc_vec->data[d] = floor(L_vec->data[d]/cutoff);
     }
     
     ParticleSystem *ps = new ParticleSystem(Np, sdim, Ni, L_vec->data);
     for (int i=0; i<Np*sdim; i++){
     ps->position->data[i] = gsl_rng_uniform(r)*L;
     std::cout << "value: " << ps->position->data[i] <<".\n";
     }
     LcGrid *lcgrid = new LcGrid(ps,nc_vec->data);
     //std::vector<Particle> =
     //Particle *particle;
     //ParticleList *particleList;
     Particle **particle_array = (Particle**) malloc(Np*sizeof(Particle*));
     ParticleList **particleList_array =(ParticleList**) malloc(Np*sizeof(ParticleList*));
     for (int i=0; i<Np; i++)
     {
     //particle = (Particle) malloc(sizeof(Particle));
     particle_array[i] = new Particle(&(ps->position->data[i * ps->position->tda]),
     &(ps->momentum->data[i * ps->momentum->tda]),
     &(ps->mass->data[i]),&(ps->force->data[i * ps->force->tda]),
     &(ps->U->data[i]),
     i);
     //particleList = (ParticleList) malloc(sizeof(ParticleList));
     particleList_array[i] = new ParticleList(particle_array[i]);
     lcgrid->addParticleList(particleList_array[i]);
     //lcgrid->findIndex( &(ps->position->data[i * ps->position->tda]), ic_vec->data);
     //linear_index = lcgrid->index(ic_vec->data);
     //lcgrid->addParticle(&particle);
     if( lcgrid->stamp->max_index > Np){
     printf("Pointer error 1");
     }
     
     }
     if( lcgrid->stamp->max_index > Np){
     printf("Pointer error 2");
     }
     lcgrid->sortParticles();
     
     double k_stiffness = 5.0;
     double r_cutoff = 1.0;
     DPDPot *potential = new DPDPot(sdim, k_stiffness, r_cutoff);
     ps->addPotential(potential);
     
     double stepsize = .01;
     
     long int nsample = 10000;
     int modprnt = 1;
     BufferedOutputSheduler *outp = new BufferedOutputSheduler(nsample, modprnt, ps);
     
     VelocityVerlet *sampler = new VelocityVerlet(stepsize, outp, ps, lcgrid);
     
     time_t tstart, tend;
     tstart = time(0);
     sampler->sample();
     tend = time(0);
     
     
     for (int i=0; i<Np*sdim; i++){
     std::cout << "Final value: " << ps->position->data[i] <<".\n";
     }
     std::cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< ".\n";
     gsl_rng_free (r);
     
     */
    int n = 100;
    int m = 10;
    double tol = 0.01;
    double t = -1.0;

    int ideg = 6; int extra = 2;
    int liwsp = m +2;
    int *iwsp = (int*) malloc(liwsp * sizeof(int));
    int itrace = 1;
    int iflag = 0;
    gsl_vector *v = gsl_vector_alloc(n); gsl_vector_set_all(v, 1.0);
    //v->data[0]=1.0;
    gsl_vector *w = gsl_vector_alloc(n); gsl_vector_set_all(w, 1.0);
    gsl_vector *wdense = gsl_vector_alloc(n); gsl_vector_set_all(w, 1.0);
    gsl_matrix *densM = gsl_matrix_alloc(n, n); gsl_matrix *expdensM = gsl_matrix_alloc(n, n);
    gsl_matrix_set_identity(densM);
    gsl_matrix_scale(densM, -1);
    gsl_spmatrix *spM = gsl_spmatrix_alloc(n, n);
    gsl_spmatrix_d2sp(spM, densM);
    //gsl_spmatrix_set(spM, 0, 1, -1.0);
    //gsl_spmatrix_sp2d(densM, spM);

    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            //gsl_spmatrix_set(spM, i, j, gsl_rng_uniform (r));
            if(gsl_rng_uniform (r)>.1)
                gsl_spmatrix_set(spM, i, j, 1.0);
            }
    gsl_rng_uniform (r);
    printf("\n densM = \n");
    gsl_spmatrix_sp2d(densM, spM);
    printMatrix(densM);
    double hump;
    //dgexpv( n, m, t, v->data, w->data,tol, spM, anorm, wsp, lwsp, iwsp, liwsp, itrace, &iflag );
    dgexpv2( t, spM, v, tol, m, w, &hump);
    /*
     printf("\n A = \n");
     for (int i = 0; i < n; i++){
     for (int j = 0; j < n; j++)
     printf( "%f ", gsl_matrix_get(densM, i, j));
     printf("\n");
     }
     */
    printf("\n INPUT v = \n");
    for (int i = 0; i < n; i++)
    printf( "%f ", gsl_vector_get(v, i));
    printf("\n");

    printf("\n OUTPUT w = \n");
    for (int i = 0; i < n; i++)
    printf( "%f ", gsl_vector_get(w, i));
    printf("\n");

    gsl_matrix_scale(densM, t);
    gsl_linalg_exponential_ss(densM, expdensM, GSL_MODE_PREC(GSL_PREC_DOUBLE));
    gsl_blas_dgemv(CblasNoTrans, 1.0, expdensM, v, 0.0, wdense);
    printf("Dense computation");
    printf("\n exp(A) = \n");
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            printf( "%f ", gsl_matrix_get(expdensM, i, j));
            printf("\n");
            }

    printf("\n INPUT v = \n");
    for (int i = 0; i < n; i++)
    printf( "%f ", gsl_vector_get(v, i));
    printf("\n");

    printf("\n OUTPUT wdense = \n");
    for (int i = 0; i < n; i++)
    printf( "%f ", gsl_vector_get(wdense, i));
    printf("\n");

    printf("\n DIFF = \n");
    for (int i = 0; i < n; i++)
    printf( "%f ", gsl_vector_get(w, i) - gsl_vector_get(wdense, i));
    printf("\n");

}
