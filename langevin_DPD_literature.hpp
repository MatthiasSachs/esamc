//
//  langevin_DPD_literature.hpp
//  esamc
//
//  Created by Matthias Sachs on 20/06/2018.
//  Copyright © 2018 Matthias Sachs. All rights reserved.
//

#ifndef langevin_DPD_literature_hpp
#define langevin_DPD_literature_hpp

#include <stdio.h>
#include <stdio.h>
#include "grid.hpp"
#include "langevin_DPD.hpp"
#include <gsl/gsl_spblas.h>



class Langevin_DPD_m_vv : public Langevin_DPD_m{
public:
    Langevin_DPD_m_vv( double stepsize_a, OutputSheduler *outp_a,
                          ParticleSystem *ps_a,
                          LcGrid *lcgrid_a, double Tk_B_a, DPD_Tensor *ft, size_t n_subsetps_a);     
    void O_step(double stepsize_factor);
    void O_step(double stepsize_factor, bool need_FD_update);
    void traverse();
};

#endif /* langevin_DPD_literature_hpp */
