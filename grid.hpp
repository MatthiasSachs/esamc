//
//  grid.hpp
//  linkedCellMD
//
//  Created by Matthias Sachs on 26/12/2017.
//  Copyright © 2017 Matthias Sachs. All rights reserved.
//

#ifndef grid_hpp
#define grid_hpp
#include <stdio.h>
//#include <util.h>
#include <iostream>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spmatrix.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <vector>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "expokit_translation2.h"

#define DEFAULT 0
#define DPD 1

class PairPotential;
class Stamp;
class LcGrid;
class ParticleSystem;
class ExternalPotential;
class OutputSheduler;

class Domain{
public:
    int sdim;
    ParticleSystem *ps;
    double *L;
    Domain( int sdim_a);
    virtual void apply_boundary();
    virtual void comp_rel_position(double *rel_pos, double *pos1, double *pos2);
    virtual double comp_rel_position(double *pos1, double *pos2, int d);
};


class ParticleSystem {
    /*
     Particle system class: Stores all relevant information about the system including the system state simulation time
     */
public:
    int sdim;  // Number of dimensions of spacial domain
    int Np; // Number of simulated particles
    int dim; //Total number of dimensions of simulated system
    gsl_matrix *position;  // matrix storing position of particles
    gsl_matrix *momentum; // matrix storing momentum of particles
    gsl_vector *mass; // vector storing mass of particles
    gsl_matrix *force; // matrix storing force
    gsl_vector *U; //vector storing the contribution to the total potential energy of each particle
    gsl_vector *laplace; //vector storing the contribution of each particle to the total laplacian of the system
    /* Vector views of the above variables as */
    gsl_vector_view position_as_vec;
    gsl_vector_view momentum_as_vec;
    gsl_vector_view force_as_vec;
    
    
    gsl_spmatrix *rel_distance; // relative distence between particles
    gsl_spmatrix *rel_position;  // relative position between particles
    Domain *domain; //Domain on which particles are simulated (e.g. sdim-dimensional torus, R^sdim
    
    PairPotential *ppotential = NULL;           //array of pair potentials
    ExternalPotential *epotential = NULL;       //array of external potentials
    LcGrid *lcgrid;                             //Linked cell grid used for the computation of forces

    const gsl_rng_type * T; //Random number environment variables
    gsl_rng * r; //Random number generator for this particle system. Should be used for any generation of random numbers related to this system
    
    ParticleSystem(int Np_a, Domain *domain_a); // standard constructor
    
    void apply_boundary();                      // Resolves boundary conditions (calls the respective function of in domain attribute)
    
    void comp_rel_position(double *rel_pos, double *pos1, double *pos2); // Computes relative position between the position of two particles. Output is stored in the array *rel_pos
    
    double comp_rel_position(double *pos1, double *pos2, int d); // Returns relative position in space dimension d between the position of two particles.
    
    void addPotential(PairPotential *ppotential); // Adds *ppotential to the array of pair potentials
    
    void addPotential(ExternalPotential *epotential); //Adds *epotential to the array of external potential
    
    void initMomentum(double beta);
    void initMomentumDPD(double beta);
    void centerMomentum();
    void scaleMomentum(double beta);
};

/*
 Base class for particles used in the LinkedCell algorithm for the computation of short range forces.
 */
class Particle{
public:
    double *position;
    double *momentum;
    double *mass;
    double *force;
    double *U;
    double *laplace;
    int id;
    Particle(double *position_a, double *momentum_a, double *mass_a, double *force_a, double *U_a, double *laplace, int id_a);
};

/*
 Base class for particleLists used in the LinkedCell algorithm for the computation of short range forces.
 */
class ParticleList {
    public:
    Particle *p;
    ParticleList *next;
    ParticleList(Particle *p);
    ParticleList(double *position, double *momentum, double *mass, double *force, double *U, double *laplace, int id);
};
typedef ParticleList* Cell;
/*
 Linked cell grid implementing the linked cell algorithm in 1,2, or 3 spacial dimension.
 */
class LcGrid{
public:
    int sdim;
    int *nc;
    double *L;
    int nc_all;
    
    Cell *cells;
    ParticleSystem *ps;
    Stamp *stamp;
    
    LcGrid(ParticleSystem *ps,int *nc_a);
    /* computes forces of all particles */
    void compForce();
    /* computes potential enrgy of all particles  */
    void compPotential();
    
    /* computes Laplacian of the potential  */
    void compLaplacePotential();
    
    /* computes relative position + distance between particles. Ensures that only values corresponding to particle pairs with distance < cutoff are added to the matrices rel_position and rel_distance */
    void compRelPos(gsl_spmatrix *rel_position,gsl_spmatrix *rel_distance, double cutoff); //
    
    /* Helper functions used in the linked cell algorithm: */
    void findIndex( double* position, int *index );
    void gridIndex(int *ic);
    int sub2lin(int *ic);
    void lin2sub(int *ic, int lin_index);
    void sortParticles();
    void addParticleList(ParticleList *pl);
    static void insertList(ParticleList **root_list, ParticleList *i);
    static void deleteList(ParticleList **q);
    
    
    /* Usefull for debugging purposes: */
    void printState();
    void printStateForce();
    void printStateNeighors();
    void printStateNeighorsList();
    void printStateNeighorsListDistance();


};
/*
Base class used in linked cell algorithm to determine neighbouring cells
 */
class Stamp{
public:
    int sdim;
    int min_index;
    int max_index;
    Stamp(int sdim_a);
    
    void getNeighbour(int *kc, const int *ic, const int lin_index);
};

/*
 Base class for external potentials
 */
class ExternalPotential {
public:
    int sdim;
    ParticleSystem *ps;
    ExternalPotential(size_t sdim_a, ParticleSystem *ps_a);
    virtual void comp_force(Particle *pi);
    virtual void comp_pot(Particle *pi);
    virtual void comp_laplace(Particle *pi);
};



/*
 Base class for pair potentials
 */
class PairPotential {
public:
    int sdim;
    ParticleSystem *ps;
    PairPotential(size_t sdim_a, ParticleSystem *ps_a);
    virtual void comp_force(Particle *pi, Particle *pj);
    virtual void comp_pot(Particle *pi, Particle *pj);
    virtual void comp_laplace(Particle *pi, Particle *pj);
};

/*
 Base class for an output task. See outpusheduler.hpp for derived classes
 */
class OutputTask {
public:
    std::string variableName;
    size_t size;
    ParticleSystem *ps;
    OutputTask(ParticleSystem *ps, std::string variableName);
    virtual void comp_output(double *outputTraj);
};




/*
 Base class for an outputscheduler . See outpusheduler.chh for derived classes
 */
class OutputSheduler{
public:
    long int nsample;
    int modprnt;
    size_t traj_length;
    std::vector<OutputTask*> outputTasks;
    ParticleSystem *ps;
    virtual void feed(long int t);
    void addOutputTask(OutputTask outputTask);
};

    
#endif
