//
//  langevin_DPD.hpp
//  linkedCellMD
//
//  Created by Matthias Sachs on 13/06/2018.
//  Copyright © 2018 Matthias Sachs. All rights reserved.
//

#ifndef langevin_DPD_hpp
#define langevin_DPD_hpp

#include <stdio.h>
#include "grid.hpp"
#include "dynamics.hpp"

class InteractionTerm;


class Friction_Tensor{
public:
    ParticleSystem *ps;
    LcGrid *lcgrid;
    Friction_Tensor(ParticleSystem *ps_a,
                    LcGrid *lcgrid_a);
    virtual void update_Gamma();
    virtual void update_Sigma();
};


class DPD_Tensor : public Friction_Tensor{
public:
    
    double cutoff;
    gsl_spmatrix * weights;
    gsl_spmatrix * weights_2sqrt;
    gsl_spmatrix * Gamma;
    gsl_spmatrix * Sigma;
    gsl_matrix *Diagonal;
    gsl_matrix *Noise_Matrix;
    gsl_vector_view Noise_Matrix_as_vec;
    DPD_Tensor(ParticleSystem *ps_a,
               LcGrid *lcgrid_a,
               InteractionTerm *interactionTerm_a, double cutoff);
    InteractionTerm *interactionTerm;
    // Before calling any of the functions below make sure that ps->rel_distance and ps->rel_position is updated
    void update_weights();
    // Before calling any of the functions below make sure that this->weights is updated
    void update_Gamma();
    void update_Noise();
    // Before calling any of the functions below make sure that this->Gamma is updated
    void update_Sigma();
    
    void update_FD(); // calls update_weights() and update_Gamma() in the right order (including the update of ps->rel_distance and ps->rel_position).
    const gsl_rng_type * T;
    gsl_rng * r;
};



class InteractionTerm{
public:
    InteractionTerm(ParticleSystem *ps_a, LcGrid *lcgrid_a);
    virtual void update_weights(gsl_spmatrix *weights);
    ParticleSystem *ps;
    LcGrid *lcgrid;
};


class DPD_InteractionTerm : public InteractionTerm{
public:
    DPD_InteractionTerm(ParticleSystem *ps_a, LcGrid *lcgrid_a, double cutoff, double gamma);
    void update_weights(gsl_spmatrix *weights);
    double cutoff;
    double gamma;
    double weight_function(double distance);
};



class Langevin_DPD : public Langevin{
public:
    Langevin_DPD( double stepsize_a, OutputSheduler *outp_a,
                 ParticleSystem *ps_a,
                 LcGrid *lcgrid_a, double Tk_B_a, DPD_Tensor *ft);
    DPD_Tensor *ft;
    void traverse();
};

class Langevin_DPD_m : public Langevin_DPD{
public:
    Langevin_DPD_m( double stepsize_a, OutputSheduler *outp_a,
                          ParticleSystem *ps_a,
                          LcGrid *lcgrid_a, double Tk_B_a, DPD_Tensor *ft, size_t n_subsetps);
    size_t n_substeps;
    double substepsize;
    gsl_vector *momentum_copy;
    
    void O_step(double stepsize_factor);
    void traverse();
};


class Langevin_DPD_m_krylov : public Langevin_DPD_m{
public:
    Langevin_DPD_m_krylov( double stepsize_a, OutputSheduler *outp_a,
                   ParticleSystem *ps_a,
                   LcGrid *lcgrid_a, double Tk_B_a, DPD_Tensor *ft, size_t n_subsetps, int m_exp, double tol_exp);
    
    int m_exp;
    double tol_exp;

    void O_step(double stepsize_factor);
    void traverse();
};


class Langevin_DPD_m_krylovBAOAB : public Langevin_DPD_m_krylov{
public:
    Langevin_DPD_m_krylovBAOAB( double stepsize_a, OutputSheduler *outp_a,
                        ParticleSystem *ps_a,
                        LcGrid *lcgrid_a, double Tk_B_a, DPD_Tensor *ft, size_t n_subsetps, int m_exp, double tol_exp);
    
    void traverse();
};


#endif /* langevin_DPD_hpp */
