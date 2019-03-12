//
//  test_DPD_mBAOAB.cpp
//  linkedCellMD
//
//  Created by Matthias Sachs on 10/06/2018.
//  Copyright © 2018 Matthias Sachs. All rights reserved.
//

#include "test_DPD_mBAOAB.hpp"

void test_Verlet()
{
    
    std::cout << "Hello, World!\n";
    int Np= 5;
    int sdim = 2;
    double cutoff = 1.0;
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    

    double L = 5.0;
    double a = 5.0;
    gsl_vector *L_vec = gsl_vector_alloc(sdim);
    gsl_vector_set_all (L_vec, L);
    gsl_vector_int *nc_vec = gsl_vector_int_alloc(sdim);
    for (int d=0; d<L_vec->size; d++)
    {
        nc_vec->data[d] = floor(L_vec->data[d]/cutoff);
    }
    Torus *domain = new Torus(sdim, L_vec->data);
    //RRn *domain = new RRn(sdim, L_vec->data);
    ParticleSystem *ps = new ParticleSystem(Np, domain);
    /*
    for (int i=0; i<Np*sdim; i++){
        ps->position->data[i] = gsl_rng_uniform(r)*a;
        std::cout << "value: " << ps->position->data[i] <<".\n";
    }
     */
    LcGrid *lcgrid = new LcGrid(ps,nc_vec->data);
    Particle **particle_array = (Particle**) malloc(Np*sizeof(Particle*));
    ParticleList **particleList_array =(ParticleList**) malloc(Np*sizeof(ParticleList*));
    for (int i=0; i<Np; i++)
    {
        //particle = (Particle) malloc(sizeof(Particle));
        particle_array[i] = new Particle(&(ps->position->data[i * ps->position->tda]),
                                         &(ps->momentum->data[i * ps->momentum->tda]),
                                         &(ps->mass->data[i]),&(ps->force->data[i * ps->force->tda]),
                                         &(ps->U->data[i]),
                                         &(ps->laplace->data[i]),
                                         i);
        //particleList = (ParticleList) malloc(sizeof(ParticleList));
        particleList_array[i] = new ParticleList(particle_array[i]);
        lcgrid->addParticleList(particleList_array[i]);
        //lcgrid->findIndex( &(ps->position->data[i * ps->position->tda]), ic_vec->data);
        //linear_index = lcgrid->index(ic_vec->data);
        //lcgrid->addParticle(&particle);
        /*
        if( lcgrid->stamp->max_index > Np){
            printf("Pointer error 1");
        }
        */
    }
    /*
    if( lcgrid->stamp->max_index > Np){
        printf("Pointer error 2");
    }
     */
    lcgrid->sortParticles();
    
    //lcgrid->printState();
    
    double k_stiffness = 5.0;
    double r_cutoff = 1.0;
    //MorsePot *potential = new MorsePot(sdim, ps, 1.0, 1.0, 1.0);
    //HarmonicPairPot *potential = new HarmonicPairPot(sdim, ps, 1.0);
    DPDPot *potential = new DPDPot(sdim, ps, k_stiffness, r_cutoff);
    //LJPot *potential = new LJPot(sdim, ps,  .001, 1.0);
    gsl_vector *center = gsl_vector_calloc(sdim);
    double stiffness = 1.0;
    HarmonicPotential *epotential = new HarmonicPotential(sdim, ps, center->data, stiffness);
    ps->addPotential(potential);
    //ps->addPotential(epotential);
    double stepsize = .0001;
    
    long int nsample = 10000;
    int modprnt = 1;
    BufferedOutputShedulerU *outp = new BufferedOutputShedulerU(nsample, modprnt, ps);
    
    SingleVarOT outpt1 = SingleVarOT(ps->position, ps, "position");
    outp->addOutputTask(&outpt1);
    SingleVarOT outpt2 = SingleVarOT(ps->U, ps, "potential");
    outp->addOutputTask(&outpt2);
    SingleVarOT outpt3 = SingleVarOT(ps->momentum, ps, "momentum");
    outp->addOutputTask(&outpt3);
    
    VelocityVerlet *sampler = new VelocityVerlet(stepsize, outp, ps, lcgrid);
    
    time_t tstart, tend;
    tstart = time(0);
    sampler->sample();
    tend = time(0);
    
//     outp->h5write("/Users/msachs2/Documents/Code/outputs/fastMD_output/test/testfile.h5");
    outp->h5write("/home/xshang/Codes/2019_GLE_DPD/esamc/matlab/testfile.h5");
    /*for (int i=0; i<Np*sdim; i++){
        std::cout << "Final value: " << ps->position->data[i] <<".\n";
    }*/
    std::cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< ".\n";
    gsl_rng_free (r);
     
}

void test_mBAOAB()
{
    double Tk_B = 1.0;
    int Np = 500;
    int sdim = 3;
    double density = 3.0;
    double cutoff = 1.0;
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    
    double L = pow( Np/density, 1.0/sdim);
    double a = L;
    gsl_vector *L_vec = gsl_vector_alloc(sdim);
    gsl_vector_set_all (L_vec, L);
    gsl_vector_int *nc_vec = gsl_vector_int_alloc(sdim);
    for (int d=0; d<L_vec->size; d++)
    {
        nc_vec->data[d] = floor(L_vec->data[d]/cutoff);
    }
    
    
    Torus *domain = new Torus(sdim, L_vec->data);
//     RRn *domain = new RRn(sdim, L_vec->data);

    ParticleSystem *ps = new ParticleSystem(Np, domain);

/******************
 Initialization of Linked cell algorithm and allocation of memory
 *******************************/
     for (int i=0; i<Np*sdim; i++){
         ps->position->data[i] = gsl_rng_uniform(ps->r)*a;
     //std::cout << "value: " << ps->position->data[i] <<".\n";
     }
    ps->initMomentumDPD(1.0/Tk_B);
    
    LcGrid *lcgrid = new LcGrid(ps,nc_vec->data);
    Particle **particle_array = (Particle**) malloc(Np*sizeof(Particle*));
    ParticleList **particleList_array =(ParticleList**) malloc(Np*sizeof(ParticleList*));
    
    
    for (int i=0; i<Np; i++)
    {
        particle_array[i] = new Particle(&(ps->position->data[i * ps->position->tda]),
                                         &(ps->momentum->data[i * ps->momentum->tda]),
                                         &(ps->mass->data[i]),&(ps->force->data[i * ps->force->tda]),
                                         &(ps->U->data[i]),
                                         &(ps->laplace->data[i]),
                                         i);

        particleList_array[i] = new ParticleList(particle_array[i]);
        lcgrid->addParticleList(particleList_array[i]);
    }

    lcgrid->sortParticles();
/*************************************************************************/
    //lcgrid->printState();
    
    double k_stiffness = 25.0;
    double r_cutoff = 1.0;
    DPDPot *potential = new DPDPot(sdim, ps, k_stiffness, r_cutoff);
    ps->addPotential(potential);
    double stepsize = .05*pow(1.15, 7.0);
    double Time = 1000;
    
    long int nsample = floor(Time/stepsize);
    int modprnt = 1;
    BufferedOutputShedulerU *outp = new BufferedOutputShedulerU(nsample, modprnt, ps);
    
    SingleVarOT outpt1 = SingleVarOT(ps->position, ps, "position");
    outp->addOutputTask(&outpt1);
    SingleVarOT outpt2 = SingleVarOT(ps->U, ps, "potential");
    outp->addOutputTask(&outpt2);
    SingleVarOT outpt3 = SingleVarOT(ps->momentum, ps, "momentum");
    outp->addOutputTask(&outpt3);
    LaplaceOT outpt4 = LaplaceOT(ps);
    outp->addOutputTask(&outpt4);
    SingleVarOT outpt5 = SingleVarOT(ps->force, ps, "force");
    outp->addOutputTask(&outpt5);
    
    double gamma_friction = 4.5;
    DPD_InteractionTerm interaction_term =  DPD_InteractionTerm(ps, lcgrid, r_cutoff, gamma_friction);
    DPD_Tensor ft = DPD_Tensor(ps,
                                           lcgrid,
                                           &interaction_term,
                                           r_cutoff);
    
    int n_substeps = 1;
    int m_exp = 20;
    double tol_exp = 1E-2;
    Langevin_DPD_m_krylovABOBA *sampler  = new Langevin_DPD_m_krylovABOBA(stepsize, outp, ps, lcgrid, Tk_B, &ft, n_substeps, m_exp, tol_exp);
    
    /*
     Other samplers:
     
    VelocityVerlet *sampler = new VelocityVerlet(stepsize, outp, ps, lcgrid);
    Langevin_DPD_m_krylovBAOAB *sampler  = new Langevin_DPD_m_krylovBAOAB(stepsize, outp, ps, lcgrid, Tk_B, &ft, n_substeps, m_exp, tol_exp);
    Langevin_DPD_m_vv *sampler  = new Langevin_DPD_m_vv(stepsize, outp, ps, lcgrid, Tk_B, &ft, n_substeps);
    Langevin_scalar_BAOAB *sampler = new Langevin_scalar_BAOAB(stepsize, outp, ps, lcgrid, Tk_B, gamma_friction );
    */
    
    time_t tstart, tend;
    tstart = time(0);
    sampler->sample();
    tend = time(0);
    
//     outp->h5write("/Users/msachs2/Documents/Code/outputs/fastMD_output/test/testfile.h5");
    outp->h5write("/home/xshang/Codes/2019_GLE_DPD/esamc/matlab/testfile.h5");

    std::cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< ".\n";
    gsl_rng_free (r);
    
}

void test_mBAOAB_RR3()
{
    
    int Np = 500;//75;
    int sdim = 3;
    double density = 3.0;
    double cutoff = 1.0;
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    
    double L = pow( Np/density, 1.0/sdim); 
    double a = L;
    gsl_vector *L_vec = gsl_vector_alloc(sdim);
    gsl_vector_set_all (L_vec, L);
    gsl_vector_int *nc_vec = gsl_vector_int_alloc(sdim);
    for (int d=0; d<L_vec->size; d++)
    {
        nc_vec->data[d] = floor(L_vec->data[d]/cutoff);
    }
    
    
    //Torus *domain = new Torus(sdim, L_vec->data);
    RRn *domain = new RRn(sdim, L_vec->data);
    //RRn *domain = new RRn(sdim, L_vec->data);
    
    ParticleSystem *ps = new ParticleSystem(Np, domain);
    
    /******************
     Initialization of Linked cell algorithm and allocation of memory
     *******************************/
    for (int i=0; i<Np*sdim; i++){
        ps->position->data[i] = gsl_rng_uniform(r)*a;
        //std::cout << "value: " << ps->position->data[i] <<".\n";
    }
    
    LcGrid *lcgrid = new LcGrid(ps,nc_vec->data);
    Particle **particle_array = (Particle**) malloc(Np*sizeof(Particle*));
    ParticleList **particleList_array =(ParticleList**) malloc(Np*sizeof(ParticleList*));
    
    
    for (int i=0; i<Np; i++)
    {
        particle_array[i] = new Particle(&(ps->position->data[i * ps->position->tda]),
                                         &(ps->momentum->data[i * ps->momentum->tda]),
                                         &(ps->mass->data[i]),&(ps->force->data[i * ps->force->tda]),
                                         &(ps->U->data[i]),
                                         &(ps->laplace->data[i]),
                                         i);
        
        particleList_array[i] = new ParticleList(particle_array[i]);
        lcgrid->addParticleList(particleList_array[i]);
    }
    
    lcgrid->sortParticles();
    /*************************************************************************/
    //lcgrid->printState();
    
    double k_stiffness = 25.0;
    double r_cutoff = 1.0;
    //MorsePot *potential = new MorsePot(sdim, ps, 1.0, 1.0, 1.0);
    //HarmonicPairPot *potential = new HarmonicPairPot(sdim, ps, 1.0);
    DPDPot *potential = new DPDPot(sdim, ps, k_stiffness, r_cutoff);
    //LJPot *potential = new LJPot(sdim, ps,  .001, 1.0);
    gsl_vector *center = gsl_vector_calloc(sdim);
    double stiffness = .5;
    HarmonicPotential *epotential = new HarmonicPotential(sdim, ps, center->data, stiffness);
    
    ps->addPotential(potential);
    ps->addPotential(epotential);
    double stepsize = .05;
    
    long int nsample = 1000;
    int modprnt = 1;
    BufferedOutputShedulerU *outp = new BufferedOutputShedulerU(nsample, modprnt, ps);
    
    SingleVarOT outpt1 = SingleVarOT(ps->position, ps, "position");
    outp->addOutputTask(&outpt1);
    SingleVarOT outpt2 = SingleVarOT(ps->U, ps, "potential");
    outp->addOutputTask(&outpt2);
    SingleVarOT outpt3 = SingleVarOT(ps->momentum, ps, "momentum");
    outp->addOutputTask(&outpt3);
    LaplaceOT outpt4 = LaplaceOT(ps);
    outp->addOutputTask(&outpt4);
    SingleVarOT outpt5 = SingleVarOT(ps->force, ps, "force");
    outp->addOutputTask(&outpt5);
    
    double gamma_friction = 4.5;
    double Tk_B = 1.0;
    DPD_InteractionTerm interaction_term =  DPD_InteractionTerm(ps, lcgrid, r_cutoff, gamma_friction);
    DPD_Tensor ft = DPD_Tensor(ps,
                               lcgrid,
                               &interaction_term,
                               r_cutoff);
    
    //VelocityVerlet *sampler = new VelocityVerlet(stepsize, outp, ps, lcgrid);
    int n_substeps = 1;
    //Langevin_DPD_m_vv *sampler  = new Langevin_DPD_m_vv(stepsize, outp, ps, lcgrid, Tk_B, &ft, n_substeps);
    int m_exp = 20;
    double tol_exp = 1E-2;
    Langevin_DPD_m_krylovBAOAB *sampler  = new Langevin_DPD_m_krylovBAOAB(stepsize, outp, ps, lcgrid, Tk_B, &ft, n_substeps, m_exp, tol_exp);
    //VelocityVerlet *sampler = new VelocityVerlet(stepsize, outp, ps, lcgrid);
//     Langevin_scalar_BAOAB *sampler = new Langevin_scalar_BAOAB(stepsize, outp, ps, lcgrid, Tk_B, gamma_friction );
    time_t tstart, tend;
    tstart = time(0);
    sampler->sample();
    tend = time(0);
    
//     outp->h5write("/Users/msachs2/Documents/Code/outputs/fastMD_output/test/testfile.h5");
    outp->h5write("/home/xshang/Codes/2019_GLE_DPD/esamc/matlab/testfile.h5");
    
    std::cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< ".\n";
    gsl_rng_free (r);
    
}


