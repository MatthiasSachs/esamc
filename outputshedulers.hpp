//
//  outputshedulers.hpp
//  linkedCellMD
//
//  Created by Matthias Sachs on 13/06/2018.
//  Copyright © 2018 Matthias Sachs. All rights reserved.
//

#ifndef outputshedulers_hpp
#define outputshedulers_hpp

#include <stdio.h>
#include "grid.hpp"

class BufferedOutputSheduler : public OutputSheduler{
public:
    std::vector<gsl_matrix*> outputTrajs;
    void feed(long int t);
    void addOutputTask(OutputTask *outputTask);
    BufferedOutputSheduler(long int nsample_a,
                           int modprnt_a,
                           ParticleSystem *ps_a);
    void h5write(std::string outputpath);
    virtual void updateStats();
    
};
class BufferedOutputShedulerU : public BufferedOutputSheduler{
public:
    BufferedOutputShedulerU(long int nsample_a,
                            int modprnt_a,
                            ParticleSystem *ps_a);
    void updateStats();
};

/*

 */
class SingleVarOT : public OutputTask{
public:
    double *variable;
    SingleVarOT(ParticleSystem *ps, std::string variableName); // only used by subclasses
    SingleVarOT(gsl_vector *data, ParticleSystem *ps, std::string variableName);
    SingleVarOT(gsl_matrix *data, ParticleSystem *ps, std::string variableName_a);
    virtual void comp_output(double *outputTraj);
};

class LaplaceOT : public SingleVarOT{
public:
    double *variable;
    LaplaceOT(ParticleSystem *ps);
    void comp_output(double *outputTraj);
};

class ConfigTOutputTask : public OutputTask{
public:
    ParticleSystem *ps;
    ConfigTOutputTask(ParticleSystem *ps);
    void comp_output(double *outputTraj);
};


#endif /* outputshedulers_hpp */
