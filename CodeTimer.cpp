//
//  CodeTimer.cpp
//  esamc
//
//  Created by Matthias Sachs on 12/03/2019.
//  Copyright © 2019 Matthias Sachs. All rights reserved.
//

#include "CodeTimer.hpp"

#define sqr(x) ((x)*(x))
#define VERIFY(a) do { if((a)<0) { fprintf(stderr,"Failure line in file CodeTimer.cpp %d.\n",__LINE__); exit(-1);}}while(0)
#define GSL_VERIFY(a) do { if((a) != 0) { fprintf(stderr,"Failure line in file CodeTimer.cpp %d.\n",__LINE__); exit(-1);}}while(0)

CodeTimerTraj force_timer;
CodeTimerTraj tensor_timer;
CodeTimerTraj dstep_timer;
CodeTimerTraj fstep_timer;

CodeTimer::CodeTimer(){
    this->counter = 0;
    this->tstep = 0;
};

void CodeTimer::add_time(long int microseconds){};
void CodeTimer::h5write(std::string outputpath){};

void CodeTimer::increment_time(){
    this->tstep+=1;
}
void CodeTimer::start_clock(){
    this->start = std::chrono::high_resolution_clock::now();
}

void CodeTimer::record_time(){
    this->end = std::chrono::high_resolution_clock::now();
    auto elapsed = this->end  - this->start;
    this->add_time(std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count());
}

CodeTimerTraj::CodeTimerTraj() : CodeTimer::CodeTimer(){};
CodeTimerTraj::CodeTimerTraj(long int nsample_a) : CodeTimer::CodeTimer(){
    this->tstep = 0;
    this->nsample=nsample_a;
    this->time_vec= (long int *) calloc(this->nsample, sizeof(*this->time_vec));
}
void CodeTimerTraj::add_time(long int microseconds){
    this->time_vec[this->tstep] += microseconds;
    this->counter+=1;
}

void CodeTimerTraj::h5write(std::string outputpath) {
    hid_t   file_id = H5Fcreate(outputpath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    VERIFY(file_id);
    hid_t group_id = H5Gcreate2(file_id, "/traj", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dim_dataspace[1] = {(hsize_t) this->nsample};
    VERIFY(H5LTmake_dataset_long ( group_id, "values in microseconds", 1, dim_dataspace, this->time_vec));
    VERIFY(H5Gclose(group_id));
    VERIFY(H5Fclose(file_id));
}





CodeTimerCumulative::CodeTimerCumulative() : CodeTimer::CodeTimer(){
    this->cumsum = 0.0;
    this->counter +=1;
};
void CodeTimerCumulative::add_time(long int microseconds){
    this->cumsum+=microseconds;
}


