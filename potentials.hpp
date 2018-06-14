//
//  potentials.hpp
//  linkedCellMD
//
//  Created by Matthias Sachs on 13/06/2018.
//  Copyright © 2018 Matthias Sachs. All rights reserved.
//

#ifndef potentials_hpp
#define potentials_hpp

#include <stdio.h>
#include "grid.hpp"
class HarmonicPotential : public ExternalPotential {
public:
    double *center;
    double stiffness;
    HarmonicPotential(size_t sdim_a, ParticleSystem *ps_a, double *center_a, double stiffness_a);
    void comp_force(Particle *pi);
    void comp_pot(Particle *pi);
};

class HarmonicPairPot : public PairPotential {
public:
    double stiffness;
    HarmonicPairPot(int sdim_a, ParticleSystem *ps_a, double stiffness) ;
    void comp_force(Particle *pi, Particle *pj);
    void comp_pot(Particle *pi, Particle *pj);
};

class LJPot : public PairPotential {
public:
    double epsilon;
    double sigma;
    LJPot(int sdim_a, ParticleSystem *ps_a, double epsilon_a, double sigma_a) ;
    void comp_force(Particle *pi, Particle *pj);
    void comp_pot(Particle *pi, Particle *pj);
};

class MorsePot : public PairPotential {
public:
    double D;
    double r_e;
    double a;
    MorsePot(int sdim_a, ParticleSystem *ps_a, double D_a, double r_e_a, double a_a);
    void comp_force(Particle *pi, Particle *pj);
    void comp_pot(Particle *pi, Particle *pj);
};

class DPDPot : public PairPotential {
public:
    double k_stiffness;
    double r_cutoff;
    double r_cutoff2;
    DPDPot(int sdim_a, ParticleSystem *ps_a, double k_stiffness_a, double r_cutoff_a);
    void comp_force(Particle *pi, Particle *pj);
    void comp_pot(Particle *pi, Particle *pj);
};


#endif /* potentials_hpp */
