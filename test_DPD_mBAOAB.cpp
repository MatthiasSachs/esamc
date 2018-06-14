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
    ParticleSystem *ps = new ParticleSystem(Np, sdim, domain);
    /*
    for (int i=0; i<Np*sdim; i++){
        ps->position->data[i] = gsl_rng_uniform(r)*a;
        std::cout << "value: " << ps->position->data[i] <<".\n";
    }
     */
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
    
    long int nsample = 1000;
    int modprnt = 100;
    BufferedOutputShedulerU *outp = new BufferedOutputShedulerU(nsample, modprnt, ps);
    
    OutputTask outpt1 = OutputTask(ps->position, "position");
    outp->addOutputTask(&outpt1);
    OutputTask outpt2 = OutputTask(ps->U, "potential");
    outp->addOutputTask(&outpt2);
    OutputTask outpt3 = OutputTask(ps->momentum, "momentum");
    outp->addOutputTask(&outpt3);
    
    VelocityVerlet *sampler = new VelocityVerlet(stepsize, outp, ps, lcgrid);
    
    time_t tstart, tend;
    tstart = time(0);
    sampler->sample();
    tend = time(0);
    
    outp->h5write("/Users/msachs2/Documents/Code/outputs/fastMD_output/test/testfile.h5");
    /*for (int i=0; i<Np*sdim; i++){
        std::cout << "Final value: " << ps->position->data[i] <<".\n";
    }*/
    std::cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< ".\n";
    gsl_rng_free (r);
     
}

void test_mBAOAB()
{
    
    std::cout << "Hello, World!\n";
    int Np=500;
    int sdim = 3;
    double cutoff = 1.0;
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    
    double L = 10.0;
    double a = 10.0;
    gsl_vector *L_vec = gsl_vector_alloc(sdim);
    gsl_vector_set_all (L_vec, L);
    gsl_vector_int *nc_vec = gsl_vector_int_alloc(sdim);
    for (int d=0; d<L_vec->size; d++)
    {
        nc_vec->data[d] = floor(L_vec->data[d]/cutoff);
    }
    Torus *domain = new Torus(sdim, L_vec->data);
    //RRn *domain = new RRn(sdim, L_vec->data);
    ParticleSystem *ps = new ParticleSystem(Np, sdim, domain);
    
     for (int i=0; i<Np*sdim; i++){
     ps->position->data[i] = gsl_rng_uniform(r)*a;
     //std::cout << "value: " << ps->position->data[i] <<".\n";
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

        
    }

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
    double stepsize = .5;
    
    long int nsample = 100;
    int modprnt = 1;
    BufferedOutputShedulerU *outp = new BufferedOutputShedulerU(nsample, modprnt, ps);
    
    OutputTask outpt1 = OutputTask(ps->position, "position");
    outp->addOutputTask(&outpt1);
    OutputTask outpt2 = OutputTask(ps->U, "potential");
    outp->addOutputTask(&outpt2);
    OutputTask outpt3 = OutputTask(ps->momentum, "momentum");
    outp->addOutputTask(&outpt3);

    double gamma_friction =5.0;
    double Tk_B=1.0;
    DPD_InteractionTerm interaction_term =  DPD_InteractionTerm(ps, lcgrid, r_cutoff, gamma_friction);
    DPD_Tensor ft = DPD_Tensor(ps,
                                           lcgrid,
                                           &interaction_term,
                                           r_cutoff);
    //VelocityVerlet *sampler = new VelocityVerlet(stepsize, outp, ps, lcgrid);
    int n_substeps = 5;
    int m_exp = 20;
    double tol_exp = 1E-4;
    Langevin_DPD_mBAOAB *sampler  = new Langevin_DPD_mBAOAB(stepsize, outp, ps, lcgrid, Tk_B, &ft, n_substeps, m_exp, tol_exp);
    
    time_t tstart, tend;
    tstart = time(0);
    sampler->sample();
    tend = time(0);
    
    outp->h5write("/Users/msachs2/Documents/Code/outputs/fastMD_output/test/testfile.h5");

    std::cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< ".\n";
    gsl_rng_free (r);
    
}
