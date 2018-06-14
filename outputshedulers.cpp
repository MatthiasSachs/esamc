//
//  outputshedulers.cpp
//  linkedCellMD
//
//  Created by Matthias Sachs on 13/06/2018.
//  Copyright © 2018 Matthias Sachs. All rights reserved.
//

#include "outputshedulers.hpp"

#define sqr(x) ((x)*(x))
#define VERIFY(a) do { if((a)<0) { fprintf(stderr,"Failure line in file outputshedulers.cpp %d.\n",__LINE__); exit(-1);}}while(0)
#define GSL_VERIFY(a) do { if((a) != 0) { fprintf(stderr,"Failure line in file outputshedulers.cpp %d.\n",__LINE__); exit(-1);}}while(0)

BufferedOutputSheduler::BufferedOutputSheduler(long int nsample_a,
                                               int modprnt_a,
                                               ParticleSystem *ps_a){
    this->nsample = nsample_a;
    this->modprnt = modprnt_a;
    this->traj_length = (size_t) ceil( ((double)this->nsample)/this->modprnt);
    this->ps = ps_a;
}
void BufferedOutputSheduler::addOutputTask(OutputTask *outputTask){
    this->outputTasks.push_back(outputTask);
    gsl_matrix *output_traj = gsl_matrix_calloc(this->traj_length, outputTask->size);
    this->outputTrajs.push_back(output_traj);
    //insert(vector_size, outputTask);
}
void BufferedOutputSheduler::updateStats(){};

void BufferedOutputSheduler::feed(long int t){
    if (t % this->modprnt == 0) {
        size_t c = t / this->modprnt;
        this->updateStats();
        for (size_t i = 0; i < this->outputTasks.size(); i++) {
            size_t outputsize = this->outputTasks[i]->size;
            for (size_t d=0; d < outputsize; d++) {
                this->outputTrajs[i]->data[outputsize * c + d] = this->outputTasks[i]->variable[d];
            }
            
        }
    }
}
void BufferedOutputSheduler::h5write(std::string outputpath) {
    hid_t   file_id = H5Fcreate(outputpath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    VERIFY(file_id);
    hid_t group_id = H5Gcreate2(file_id, "/traj", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for (size_t i = 0; i < this->outputTrajs.size(); i++) {
        hsize_t dim_dataspace[2] = {this->traj_length, this->outputTasks[i]->size};
        VERIFY(H5LTmake_dataset_double ( group_id, this->outputTasks[i]->variableName.c_str(), 2, dim_dataspace, this->outputTrajs[i]->data));
    }
    VERIFY(H5Gclose(group_id));
    VERIFY(H5Fclose(file_id));
}

BufferedOutputShedulerU::BufferedOutputShedulerU(long int nsample_a,
                                                 int modprnt_a, ParticleSystem *ps_a):BufferedOutputSheduler(nsample_a,modprnt_a, ps_a){};

void BufferedOutputShedulerU::updateStats(){
    this->ps->lcgrid->compPotential();
};
