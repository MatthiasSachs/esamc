//
//  CodeTimer.hpp
//  esamc
//
//  Created by Matthias Sachs on 12/03/2019.
//  Copyright © 2019 Matthias Sachs. All rights reserved.
//

#ifndef CodeTimer_hpp
#define CodeTimer_hpp

#include <stdio.h>
#include <iostream>
#include <hdf5.h>
#include <hdf5_hl.h>


class CodeTimer{
public:
    long int counter;
    long int tstep;
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    
    void start_clock();
    void record_time();
    
    CodeTimer();
    void increment_time();
    virtual void add_time(long int microseconds);
    virtual void h5write(std::string outputpath);
    
};




class CodeTimerTraj : public CodeTimer{
public:
    long int nsample;
    long int *time_vec;
    
    CodeTimerTraj();
    CodeTimerTraj(long int nsample_a);
    void add_time(long int microseconds);
    void h5write(std::string outputpath);
};

class CodeTimerCumulative : public CodeTimer{
public:
    long long cumsum;
    
    CodeTimerCumulative();
    void add_time(long int microseconds);
};

extern CodeTimerTraj force_timer;
extern CodeTimerTraj tensor_timer;
extern CodeTimerTraj dstep_timer;
extern CodeTimerTraj fstep_timer;

#endif /* CodeTimer_hpp */
