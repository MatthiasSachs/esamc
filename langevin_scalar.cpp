//
//  langevin_scalar.cpp
//  esamc
//
//  Created by Matthias Sachs on 05/02/2019.
//  Copyright © 2019 Matthias Sachs. All rights reserved.
//

#include "langevin_scalar.hpp"


Langevin_scalar::Langevin_scalar( double stepsize_a, OutputSheduler *outp_a,
                                 ParticleSystem *ps_a, LcGrid *lcgrid_a, double Tk_B_a, double gamma): Langevin::Langevin( stepsize_a, outp_a, ps_a,lcgrid_a,Tk_B_a){
    this->gamma = gamma;
    this->alpha1 = exp(- .5 * this->stepsize * this->gamma);
    this->zeta1 = sqrt(1.0 - this->alpha1* this->alpha1) * sqrt(this->Tk_B);
    this->alpha2 = exp(-this->stepsize * this->gamma);
    this->zeta2 = sqrt(1.0 - this->alpha2* this->alpha2) * sqrt(this->Tk_B);

};

void Langevin_scalar::O_step(double alpha, double zeta){
    double* momentum = this->ps->momentum->data;
    size_t size = this->ps->momentum->size1 * this->ps->momentum->size2;
    for (size_t i = 0; i < size ; i++) {
        momentum[i] *= alpha;
        momentum[i] += zeta * gsl_ran_gaussian (this->ps->r, 1.0);
    }
}
void Langevin_scalar::traverse(){};

Langevin_scalar_BAOAB::Langevin_scalar_BAOAB( double stepsize_a, OutputSheduler *outp_a,
                                             ParticleSystem *ps_a, LcGrid *lcgrid_a, double Tk_B_a, double gamma) : Langevin_scalar( stepsize_a, outp_a, ps_a,lcgrid_a, Tk_B_a, gamma){};
void Langevin_scalar_BAOAB::traverse(){
    this->B_step(.5);
    this->A_step(.5);
    this->O_step(this->alpha2,this->zeta2);
    this->A_step(.5);
    this->ps->apply_boundary();
    this->lcgrid->compForce();
    this->B_step(.5 );
};


