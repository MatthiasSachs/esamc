//
//  langevin_DPD_literature.cpp
//  esamc
//
//  Created by Matthias Sachs on 20/06/2018.
//  Copyright © 2018 Matthias Sachs. All rights reserved.
//

#include "langevin_DPD_literature.hpp"


Langevin_DPD_m_vv::Langevin_DPD_m_vv( double stepsize_a, OutputSheduler *outp_a, ParticleSystem *ps_a,LcGrid *lcgrid_a, double Tk_B_a, DPD_Tensor *ft_a, size_t n_subsetps_a) : Langevin_DPD_m::Langevin_DPD_m( stepsize_a, outp_a, ps_a,lcgrid_a,Tk_B_a, ft_a ,n_subsetps_a){};

void Langevin_DPD_m_vv::O_step(double stepsize_factor){
    this->O_step(stepsize_factor, true);
}

void Langevin_DPD_m_vv::O_step(double stepsize_factor, bool need_FD_update){
    double hs = this->substepsize * stepsize_factor;
    double var_sqrt = sqrt(hs*this->Tk_B);
    if (need_FD_update) {
            this->ft->update_FD();
    }
    for (int k =0; k < this->n_substeps; k++) {
        // D-part
        gsl_vector_memcpy(this->momentum_copy, &this->ps->momentum_as_vec.vector);
        gsl_spblas_dgemv(CblasNoTrans, hs, this->ft->Gamma, this->momentum_copy, 1.0, &this->ps->momentum_as_vec.vector);
        // F part
        this->ft->update_Noise();
        gsl_blas_daxpy(var_sqrt, &this->ft->Noise_Matrix_as_vec.vector, &this->ps->momentum_as_vec.vector);
    }
    
}
void Langevin_DPD_m_vv::traverse(){
    this->O_step(.5, false);
    this->B_step(.5);
    this->A_step(1.0);
    this->ps->apply_boundary();
    this->lcgrid->compForce();
    this->O_step(.5, true);
    this->B_step(.5);
}
