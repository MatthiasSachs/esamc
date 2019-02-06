//
//  potentials.cpp
//  linkedCellMD
//
//  Created by Matthias Sachs on 13/06/2018.
//  Copyright © 2018 Matthias Sachs. All rights reserved.
//

#include "potentials.hpp"
#define sqr(x) ((x)*(x))

HarmonicPotential::HarmonicPotential(size_t sdim_a, ParticleSystem *ps_a, double *center_a, double stiffness_a)  : ExternalPotential::ExternalPotential(sdim_a, ps_a){
    this->center = center_a;
    this->stiffness = stiffness_a;
}
void HarmonicPotential::comp_force(Particle *pi){
    for (int d =0; d < ps->sdim; d++) {
        pi->force[d] += -this->stiffness * this->ps->comp_rel_position(this->center, pi->position, d);
    }
};
void HarmonicPotential::comp_pot(Particle *pi){
    double r2, dist;
    r2= 0;
    for (int d =0; d < ps->sdim; d++) {
        dist = this->ps->comp_rel_position(this->center, pi->position, d);
        r2 += sqr(dist);
    }
    *(pi->U) += .5 * stiffness * r2;
};

void HarmonicPotential::comp_laplace(Particle *pi){
    *(pi->laplace) += stiffness * ps->sdim;
};


HarmonicPairPot::HarmonicPairPot(int sdim_a, ParticleSystem *ps_a, double stiffness) : PairPotential(sdim_a, ps_a) {
    this->stiffness = stiffness;
}

void HarmonicPairPot::comp_force(Particle *pi, Particle *pj){
    for (int d = 0; d < this->sdim; d++){
        pi->force[d] += -this->stiffness * this->ps->comp_rel_position(pj->position, pi->position, d);
    }
}

void HarmonicPairPot::comp_pot(Particle *pi, Particle *pj){
    double r2,dist;
    r2 = 0;
    for (int d = 0; d < this->sdim; d++){
        dist = this->ps->comp_rel_position(pj->position, pi->position, d);
        r2 += sqr(dist);
    }
    *(pi->U) += .5 * .5 * stiffness * r2;
}

void HarmonicPairPot::comp_laplace(Particle *pi, Particle *pj){
    *(pi->laplace) += .5 * stiffness;
}

LJPot::LJPot(int sdim_a, ParticleSystem *ps_a, double epsilon_a, double sigma_a) : PairPotential(sdim_a, ps_a) {
    this->epsilon = epsilon_a;
    this->sigma = sigma_a;
}

void LJPot::comp_force(Particle *pi, Particle *pj){
    double r2,dist;
    r2 = 0;
    for (int d = 0; d < this->sdim; d++){
        dist = this->ps->comp_rel_position(pi->position, pj->position, d);
        r2 += sqr(dist);
    }
    double s = sqr(sigma) / r2;
    s = sqr(s) * s;
    double f = 24 * epsilon * s / r2 * (1 - 2 * s);
    for (int d=0; d<this->sdim; d++)
        pi->force[d] += f * this->ps->comp_rel_position(pi->position, pj->position, d);
}

void LJPot::comp_pot(Particle *pi, Particle *pj){
    double r2,dist;
    r2 = 0;
    for (int d = 0; d < this->sdim; d++){
        dist = this->ps->comp_rel_position(pi->position, pj->position, d);
        r2 += sqr(dist);
    }
    double s = sqr(sigma) / r2;
    s = sqr(s) * s;
    *(pi->U) += .5 * 4.0 * epsilon * ( sqr(s) - s);
}
void LJPot::comp_laplace(Particle *pi, Particle *pj){
    printf("Warning: function LJPot::comp_laplace(Particle *pi, Particle *pj) not implemented");
};

MorsePot::MorsePot(int sdim_a, ParticleSystem *ps_a, double D_a, double r_e_a, double a_a) : PairPotential(sdim_a, ps_a){
    this->D = D_a;
    this->r_e = r_e_a;
    this->a = a_a;
}

void MorsePot::comp_force(Particle *pi, Particle *pj){
    
    double r1, r2, dist;
    r2=0.0;
    for (int d = 0; d < this->sdim; d++){
        dist = this->ps->comp_rel_position(pi->position, pj->position, d);
        r2 += sqr(dist);
    }
    r1 = sqrt(r2);
    double cexpm = exp(-this->a*(r1 - this->r_e));
    double cf = -2.0 * this->D * this->a * (cexpm * (1.0 - cexpm) / r1);
    for (int d=0; d<this->sdim; d++)
        pi->force[d] += -cf * (pj->position[d] - pi->position[d]);
}

void MorsePot::comp_pot(Particle *pi, Particle *pj){
    double r1, r2;
    r2=0.0;
    for (int d = 0; d < this->sdim; d++)
        r2 += sqr(pj->position[d] - pi->position[d]);
    r1 = sqrt(r2);
    double ctemp = (1 - exp(-this->a * (r1 - this->r_e)));
    *(pi->U) += .5 * this->D * sqr(ctemp);
}
void MorsePot::comp_laplace(Particle *pi, Particle *pj){
    printf("Warning: function MorsePot::comp_laplace(Particle *pi, Particle *pj) not implemented");
};

DPDPot::DPDPot(int sdim_a, ParticleSystem *ps_a, double k_stiffness_a, double r_cutoff_a) : PairPotential(sdim_a, ps_a){
    
    this->k_stiffness = k_stiffness_a;
    this->r_cutoff = r_cutoff_a;
    this->r_cutoff2 = sqr(r_cutoff);
}

void DPDPot::comp_force(Particle *pi, Particle *pj){
    double r1, r2, dist,factor;
    r2=0.0;
    for (int d = 0; d < this->sdim; d++){
        dist = this->ps->comp_rel_position(pi->position, pj->position, d);
        r2 += sqr(dist);
    }
    if (r2 < r_cutoff2 && r2 > 0)
    {
        r1 = sqrt(r2);
        factor = -k_stiffness * (r_cutoff - r1)/(r_cutoff);
        for (int d = 0; d < this->sdim; d++)
            pi->force[d] += factor * this->ps->comp_rel_position(pi->position, pj->position, d)/r1;
    }
}
void DPDPot::comp_pot(Particle *pi, Particle *pj){
    double r1, r2, dist;
    r2=0.0;
    for (int d = 0; d < this->sdim; d++){
        dist = this->ps->comp_rel_position(pi->position, pj->position, d);
        r2 += sqr(dist);
    }
    if (r2 < r_cutoff2  )
    {
        r1 = sqrt(r2);
        *(pi->U) +=  .5 * .5 * k_stiffness * r_cutoff * sqr(1.0-r1/r_cutoff);
    }
}

void DPDPot::comp_laplace(Particle *pi, Particle *pj){
    double r1, r2, dist;
    double dist2[3];
    r2 = 0;
    for (int d = 0; d < this->sdim; d++){
        dist = this->ps->comp_rel_position(pj->position, pi->position, d);
        dist2[d] = sqr(dist);
        r2 += dist2[d];
    }
    if (r2 < r_cutoff2 && r2 > 0 ){
        r1 = sqrt(r2);
        //for (int d = 0; d < this->sdim; d++){
        //    *(pi->laplace) += .5 * k_stiffness * ( 3.0/r_cutoff - 2.0/r1 );
        //}
        for (int d = 0; d < this->sdim; d++){
            *(pi->laplace) += .25 * k_stiffness * ((8.0 * dist2[d])/r_cutoff2-(4.0 * (1.0-r2/r_cutoff))/r_cutoff) * r_cutoff;                                      
        }
    }
}
