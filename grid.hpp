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
#include <util.h>
#include <iostream>
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
    
public:
    int sdim;
    int Np;
    int dim;
    gsl_matrix *position;
    gsl_matrix *momentum;
    gsl_vector *mass;
    gsl_matrix *force;
    gsl_vector *U;
    gsl_vector_view position_as_vec;
    gsl_vector_view momentum_as_vec;
    gsl_vector_view force_as_vec;
    gsl_spmatrix *rel_distance;
    gsl_spmatrix *rel_position;
    Domain *domain;
    
    ParticleSystem(int Np_a, int sdim_a, Domain *domain_a);
    void apply_boundary();
    void comp_rel_position(double *rel_pos, double *pos1, double *pos2);
    double comp_rel_position(double *pos1, double *pos2, int d);
    void addPotential(PairPotential *ppotential);
    void addPotential(ExternalPotential *epotential);
    PairPotential *ppotential = NULL;
    ExternalPotential *epotential = NULL;
    LcGrid *lcgrid;
};


class Particle{
public:
    double *position;
    double *momentum;
    double *mass;
    double *force;
    double *U;
    int id;
    Particle(double *position_a, double *momentum_a, double *mass_a, double *force_a, double *U_a, int id_a);
};


class ParticleList {
    public:
    Particle *p;
    ParticleList *next;
    ParticleList(Particle *p);
    ParticleList(double *position, double *momentum, double *mass, double *force, double *U, int id);
};
typedef ParticleList* Cell;

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
    
    /* computes relative position + distance between particles. Ensures that only values corresponding to particle pairs with distance < cutoff are added to the matrices rel_position and rel_distance */
    void compRelPos(gsl_spmatrix *rel_position,gsl_spmatrix *rel_distance, double cutoff); //
    
    /* Used in the linked cell algorithm: */
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

class Stamp{
public:
    int sdim;
    int min_index;
    int max_index;
    Stamp(int sdim_a);
    
    void getNeighbour(int *kc, const int *ic, const int lin_index);
};


class ExternalPotential {
public:
    int sdim;
    ParticleSystem *ps;
    ExternalPotential(size_t sdim_a, ParticleSystem *ps_a);
    virtual void comp_force(Particle *pi);
    virtual void comp_pot(Particle *pi);
};



class PairPotential {
public:
    int sdim;
    ParticleSystem *ps;
    PairPotential(size_t sdim_a, ParticleSystem *ps_a);
    virtual void comp_force(Particle *pi, Particle *pj);
    virtual void comp_pot(Particle *pi, Particle *pj);
};


class OutputTask{
public:
    double *variable;
    std::string variableName;
    double *target;
    size_t size;
    OutputTask(gsl_vector *data, std::string variableName);
    OutputTask(gsl_matrix *data, std::string variableName_a);
};


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
