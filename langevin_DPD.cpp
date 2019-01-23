//
//  langevin_DPD.cpp
//  linkedCellMD
//
//  Created by Matthias Sachs on 13/06/2018.
//  Copyright © 2018 Matthias Sachs. All rights reserved.
//

#include "langevin_DPD.hpp"

Friction_Tensor::Friction_Tensor(ParticleSystem *ps_a,
                                 LcGrid *lcgrid_a){
    this->ps = ps_a;
    this->lcgrid = lcgrid_a;
}
void Friction_Tensor::update_Gamma(){}
void Friction_Tensor::update_Sigma(){}


DPD_Tensor::DPD_Tensor(ParticleSystem *ps_a,
                       LcGrid *lcgrid_a,
                       InteractionTerm *interactionTerm_a, double cutoff)
:Friction_Tensor::Friction_Tensor(ps_a,lcgrid_a)
{
    this->cutoff = cutoff;
    this->ps = ps_a;
    this->lcgrid = lcgrid_a;
    this->Diagonal = gsl_matrix_calloc(this->ps->Np,this->ps->sdim*this->ps->sdim);
    this->interactionTerm = interactionTerm_a;
    this->weights = gsl_spmatrix_alloc(this->ps->Np, this->ps->Np);
    this->weights_2sqrt = gsl_spmatrix_alloc(this->ps->Np, this->ps->Np);
    this->Gamma =  gsl_spmatrix_alloc(this->ps->dim, this->ps->dim);
    this->Noise_Matrix = gsl_matrix_calloc(this->ps->Np, this->ps->sdim);
    this->Noise_Matrix_as_vec = gsl_vector_view_array(this->Noise_Matrix->data, this->Noise_Matrix->size1*this->Noise_Matrix->size2);
    
    gsl_rng_env_setup();
    this->T = gsl_rng_default;
    this->r = gsl_rng_alloc (this->T);
}
void DPD_Tensor::update_FD(){
    this->lcgrid->compRelPos(this->ps->rel_position, this->ps->rel_distance, this->cutoff);
    this->update_weights();
    this->update_Gamma();
}
void DPD_Tensor::update_weights(){
    this->interactionTerm->update_weights(this->weights);
    for (size_t k =0; k < this->weights->nz; k++) {
        gsl_spmatrix_set(this->weights_2sqrt, this->weights->i[k], this->weights->p[k], sqrt(2.0*this->weights->data[k]));
    }
}

void DPD_Tensor::update_Gamma(){
    gsl_spmatrix_set_zero(this->Gamma);
    gsl_matrix_set_zero(this->Diagonal);
    gsl_spmatrix *rd = this->ps->rel_distance;
    gsl_spmatrix *rp = this->ps->rel_position;
    double val;
    double weight;
    double distance;
    size_t id1,id2;
    int sdim = this->ps->sdim;
    for (size_t k =0; k < rd->nz; k++) {
        id1 = rd->i[k];
        id2 = rd->p[k];
        distance = rd->data[k];
        weight = gsl_spmatrix_get(this->weights, id1, id2);
        //add the off-diagonal terms to Gamma
        if (distance > 0.0){
            for (size_t d1 =0; d1 < this->ps->sdim; d1++){
                for (size_t d2 =0; d2 < this->ps->sdim; d2++){
                    val = -weight   * gsl_spmatrix_get(rp, sdim*id1+d1, sdim*id2+d1)
                    * gsl_spmatrix_get(rp, sdim*id1+d2, sdim*id2+d2)
                    /(distance*distance);
                    gsl_spmatrix_set(this->Gamma, sdim*id1+d1, sdim*id2+d2, -val);
                    //compute the diagonal terms
                    this->Diagonal->data[id1 * this->Diagonal->tda + sdim*d1+d2] += val;
                    
                }
            }
        }
    }
    //add the diagonal terms to Gamma
    for (id1 =0; id1 < this->ps->Np; id1++) {
        for (size_t d1 =0; d1 < this->ps->sdim; d1++){
            for (size_t d2 =0; d2 < this->ps->sdim; d2++){
                val = this->Diagonal->data[id1 * this->Diagonal->tda + sdim*d1+d2];
                if (val != 0.0){
                    gsl_spmatrix_set(this->Gamma, sdim*id1+d1, sdim*id1+d2, val);
                }
            }
        }
    }
}

void DPD_Tensor::update_Noise(){
    gsl_spmatrix *rd = this->ps->rel_distance;
    gsl_spmatrix *rp = this->ps->rel_position;
    double val;
    double weight;
    double distance;
    double rand;
    size_t id1,id2;
    int sdim = this->ps->sdim;
    gsl_matrix_set_zero(this->Noise_Matrix);
    for (size_t k =0; k < rd->nz; k++) {
        id1 = rd->i[k];
        id2 = rd->p[k];
        distance = rd->data[k];
        weight = gsl_spmatrix_get(this->weights_2sqrt, id1, id2);
        if (id1<id2 && distance>0) {
            rand = gsl_ran_gaussian (this->r, 1.0);
            for (int d=0; d < this->ps->sdim ; d++) {
                val =gsl_spmatrix_get(rp, sdim*id1+d, sdim*id2+d)/distance * weight * rand;
                this->Noise_Matrix->data[id1 * this->Noise_Matrix->tda + d] += val;
                this->Noise_Matrix->data[id2 * this->Noise_Matrix->tda + d] += -val;            }
            
        }
    }
}
void DPD_Tensor::update_Sigma(){};


InteractionTerm::InteractionTerm(ParticleSystem *ps_a, LcGrid *lcgrid_a){
    this->ps = ps_a;
    this->lcgrid = lcgrid_a;
}
void InteractionTerm::update_weights(gsl_spmatrix *weights){};


DPD_InteractionTerm::DPD_InteractionTerm(ParticleSystem *ps_a, LcGrid *lcgrid_a, double cutoff_a, double gamma_a) : InteractionTerm(ps_a, lcgrid_a){
    this->cutoff = cutoff_a;
    this->gamma = gamma_a;
}
double DPD_InteractionTerm::weight_function(double distance){
    double val =0.0;
    if (distance < cutoff) {
        val = gamma*(1.0 - distance/cutoff);
    }
    return val;
}
void DPD_InteractionTerm::update_weights(gsl_spmatrix *weights){
    gsl_spmatrix_set_zero(weights);
    size_t id1,id2;
    gsl_spmatrix *rd = this->ps->rel_distance;
    for (size_t k =0; k < rd->nz; k++) {
        id1 = rd->i[k];
        id2 = rd->p[k];
        gsl_spmatrix_set(weights, id1, id2, this->weight_function(rd->data[k]));
    }
}



Langevin_DPD::Langevin_DPD( double stepsize_a, OutputSheduler *outp_a, ParticleSystem *ps_a, LcGrid *lcgrid_a, double Tk_B_a, DPD_Tensor *ft_a) : Langevin::Langevin( stepsize_a, outp_a, ps_a,lcgrid_a,Tk_B_a){
    this->ft = ft_a;
}
void Langevin_DPD::traverse(){};

Langevin_DPD_m::Langevin_DPD_m( double stepsize_a, OutputSheduler *outp_a, ParticleSystem *ps_a, LcGrid *lcgrid_a, double Tk_B_a, DPD_Tensor *ft_a, size_t n_substeps_a) : Langevin_DPD::Langevin_DPD( stepsize_a, outp_a, ps_a,lcgrid_a,Tk_B_a, ft_a){
    this->n_substeps = n_substeps_a;
    this->substepsize = (stepsize_a/(double)n_substeps_a);
    this->momentum_copy = gsl_vector_calloc(this->ps->dim);
}
void Langevin_DPD_m::traverse(){};

Langevin_DPD_m_krylov::Langevin_DPD_m_krylov( double stepsize_a, OutputSheduler *outp_a, ParticleSystem *ps_a, LcGrid *lcgrid_a, double Tk_B_a, DPD_Tensor *ft_a, size_t n_substeps_a, int m_exp, double tol_exp) : Langevin_DPD_m::Langevin_DPD_m( stepsize_a, outp_a, ps_a,lcgrid_a,Tk_B_a, ft_a, n_substeps_a)
{
    this->m_exp = m_exp;
    this->tol_exp = tol_exp;
}

void Langevin_DPD_m_krylov::traverse(){};

void Langevin_DPD_m_krylov::O_step(double stepsize_factor){
    double hs = this->substepsize * stepsize_factor;
    double var_sqrt = sqrt(hs*this->Tk_B);
    double hump;
    
    this->ft->update_FD();
    
    for (int k =0; k < this->n_substeps; k++) {
        // D-part
        gsl_vector_memcpy(this->momentum_copy, &this->ps->momentum_as_vec.vector);
        dgexpv2( .5 * hs, this->ft->Gamma, this->momentum_copy, this->tol_exp, this->m_exp, &this->ps->momentum_as_vec.vector, &hump);
        // F part
        this->ft->update_Noise();
        gsl_blas_daxpy(var_sqrt, &this->ft->Noise_Matrix_as_vec.vector, &this->ps->momentum_as_vec.vector);
        // D-part
        gsl_vector_memcpy(this->momentum_copy, &this->ps->momentum_as_vec.vector);
        dgexpv2( .5 * hs, this->ft->Gamma, this->momentum_copy, this->tol_exp, this->m_exp, &this->ps->momentum_as_vec.vector, &hump);
    }
}



Langevin_DPD_m_krylovBAOAB::Langevin_DPD_m_krylovBAOAB( double stepsize_a, OutputSheduler *outp_a, ParticleSystem *ps_a, LcGrid *lcgrid_a, double Tk_B_a, DPD_Tensor *ft_a, size_t n_substeps_a, int m_exp, double tol_exp): Langevin_DPD_m_krylov::Langevin_DPD_m_krylov( stepsize_a, outp_a, ps_a, lcgrid_a, Tk_B_a, ft_a, n_substeps_a, m_exp, tol_exp){};
void Langevin_DPD_m_krylovBAOAB::traverse(){
    this->B_step(.5 );
    this->A_step(.5 );
    this->ps->apply_boundary();
    this->lcgrid->compForce();
    this->O_step(1.0);
    this->A_step(.5 );
    this->ps->apply_boundary();
    this->lcgrid->compForce();
    this->B_step(.5 );
}
