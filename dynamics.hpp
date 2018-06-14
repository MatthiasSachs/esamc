//
//  dynamics.hpp
//  linkedCellMD
//
//  Created by Matthias Sachs on 13/06/2018.
//  Copyright © 2018 Matthias Sachs. All rights reserved.
//

#ifndef dynamics_hpp
#define dynamics_hpp

#include <stdio.h>
#include "grid.hpp"

class HamIntegrator{
public:
    double stepsize;
    long int t = 0;
    OutputSheduler *outp;
    ParticleSystem *ps;
    LcGrid *lcgrid;
    void sample();
    virtual void traverse();
    HamIntegrator();
    HamIntegrator(double stepsize_a,
                  OutputSheduler *outp_a,
                  ParticleSystem *ps_a,
                  LcGrid *lcgrid_a);
    void A_step(double stepsize);
    void B_step(double stepsize);
};

class VelocityVerlet : public HamIntegrator{
public:
    VelocityVerlet( double stepsize_a, OutputSheduler *outp_a,
                   ParticleSystem *ps_a,
                   LcGrid *lcgrid_a);
    void traverse();
};



class Thermostat : public HamIntegrator{
public:
    Thermostat( double stepsize_a, OutputSheduler *outp_a,
               ParticleSystem *ps_a,
               LcGrid *lcgrid_a, double Tk_B_a);
    double Tk_B; // target temperature
    virtual void traverse();
};

class Langevin : public Thermostat{
public:
    Langevin( double stepsize_a, OutputSheduler *outp_a,
             ParticleSystem *ps_a,
             LcGrid *lcgrid_a, double Tk_B_a);
    void traverse();
};


#endif /* dynamics_hpp */
