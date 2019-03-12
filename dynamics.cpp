//
//  dynamics.cpp
//  linkedCellMD
//
//  Created by Matthias Sachs on 13/06/2018.
//  Copyright © 2018 Matthias Sachs. All rights reserved.
//

#include "dynamics.hpp"

HamIntegrator::HamIntegrator(){};
HamIntegrator::HamIntegrator(double stepsize_a,
                             OutputSheduler *outp_a,
                             ParticleSystem *ps_a,
                             LcGrid *lcgrid_a){
    this->stepsize = stepsize_a;
    this->outp = outp_a;
    this->ps = ps_a;
    this->lcgrid = lcgrid_a;
}

void HamIntegrator::traverse(){};

void HamIntegrator::sample(){
    for (this->t = 0; this->t < this->outp->nsample; this->t++) {
        this->traverse();
        this->outp->feed(this->t);
        
        force_timer.increment_time();
        tensor_timer.increment_time();
        dstep_timer.increment_time();
        fstep_timer.increment_time();
    }
}

void HamIntegrator::A_step(double stepsize_factor){
    double* position = this->ps->position->data;
    double* momentum = this->ps->momentum->data;
    double h = this->stepsize * stepsize_factor;
    size_t size = this->ps->position->size1 * this->ps->position->size2;
    // mass not considered here (Only works for m[i] =1, i =1,...,N_p )
    for (size_t i = 0; i < size ; i++) {
        position[i] += h * momentum[i];
    }
}
void HamIntegrator::B_step(double stepsize_factor){
    double* momentum = this->ps->momentum->data;
    double* force = this->ps->force->data;
    double h = this->stepsize * stepsize_factor;
    size_t size = this->ps->position->size1 * this->ps->position->size2;
    for (size_t i = 0; i < size ; i++) {
        momentum[i] += h * force[i];
    }
}


VelocityVerlet::VelocityVerlet(double stepsize_a,
                               OutputSheduler *outp_a,
                               ParticleSystem *ps_a,
                               LcGrid *lcgrid_a){
    this->stepsize = stepsize_a;
    this->outp = outp_a;
    this->ps = ps_a;
    this->lcgrid = lcgrid_a;
}
void VelocityVerlet::traverse(){
    this->B_step(.5);
    this->A_step(1.0);
    this->ps->apply_boundary();
    this->lcgrid->compForce();
    //printf("%ld\n",this->t);
    //this->lcgrid->printState();
    //this->lcgrid->printStateForce();
    
    //this->lcgrid->printStateNeighors();
    //this->lcgrid->printStateNeighorsList();
    //this->lcgrid->printStateNeighorsListDistance();
    //printf("Force: \n");
    //printMatrix(this->ps->force);
    this->B_step(.5 );
};

Thermostat::Thermostat( double stepsize_a, OutputSheduler *outp_a,
                       ParticleSystem *ps_a,
                       LcGrid *lcgrid_a, double Tk_B_a): HamIntegrator::HamIntegrator( stepsize_a, outp_a, ps_a,lcgrid_a){
    this->Tk_B = Tk_B_a;
};

void Thermostat::traverse(){};

Langevin::Langevin( double stepsize_a, OutputSheduler *outp_a,
                   ParticleSystem *ps_a,
                   LcGrid *lcgrid_a, double Tk_B_a) : Thermostat::Thermostat( stepsize_a, outp_a, ps_a,lcgrid_a,Tk_B_a){
};
void Langevin::traverse(){};
void Langevin::O_step(){
    printf("Warning: O_step not implemented");
}
