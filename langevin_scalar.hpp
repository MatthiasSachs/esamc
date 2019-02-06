//
//  langevin_scalar.hpp
//  esamc
//
//  Created by Matthias Sachs on 05/02/2019.
//  Copyright © 2019 Matthias Sachs. All rights reserved.
//

#ifndef langevin_scalar_hpp
#define langevin_scalar_hpp

#include <stdio.h>
#include "grid.hpp"
#include "dynamics.hpp"

class Langevin_scalar : public Langevin {
public:
    Langevin_scalar( double stepsize_a, OutputSheduler *outp_a,
             ParticleSystem *ps_a, LcGrid *lcgrid_a, double Tk_B_a, double gamma);
    double gamma;
    double alpha1;
    double zeta1;
    double alpha2;
    double zeta2;
    void O_step(double alpha, double zeta);
    virtual void traverse();
};

class Langevin_scalar_BAOAB : public Langevin_scalar  {
public:
    Langevin_scalar_BAOAB( double stepsize_a, OutputSheduler *outp_a,
                    ParticleSystem *ps_a, LcGrid *lcgrid_a, double Tk_B_a, double gamma);
    void traverse();
};




#endif /* langevin_scalar_hpp */
