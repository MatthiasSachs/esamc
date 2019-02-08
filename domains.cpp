//
//  domains.cpp
//  linkedCellMD
//
//  Created by Matthias Sachs on 13/06/2018.
//  Copyright © 2018 Matthias Sachs. All rights reserved.
//

#include "domains.hpp"

RRn::RRn( int sdim_a, double *L_a) : Domain(sdim_a){
    this->L = L_a;
}
void RRn::apply_boundary(){};

void RRn::comp_rel_position(double *rel_pos, double *pos1, double *pos2){
    for (int d = 0; d < this->ps->sdim; d++) {
        rel_pos[d] = pos2[d] - pos1[d];
    }
}
double RRn::comp_rel_position(double *pos1, double *pos2, int d){
    return pos2[d] - pos1[d];
}

Torus::Torus( int sdim_a, double *L_a) : Domain(sdim_a){
    this->L = L_a;
}
void Torus::comp_rel_position(double *rel_pos, double *pos1, double *pos2){
    double cdiff;
    for (int d = 0; d < this->ps->sdim; d++) {
        cdiff = pos2[d] - pos1[d];
        if (cdiff >  this->L[d] * 0.5)
            cdiff = cdiff - this->L[d];
        else if (cdiff <= -this->L[d] * 0.5)
            cdiff = cdiff + this->L[d];
        rel_pos[d] = cdiff;
    }
}
double Torus::comp_rel_position(double *pos1, double *pos2, int d){
    double cdiff = pos2[d] - pos1[d];
    if (cdiff >  this->L[d] * 0.5)
        cdiff = cdiff - this->L[d];
    else if (cdiff <= -this->L[d] * 0.5)
        cdiff = cdiff + this->L[d];
    return cdiff;
}
void Torus::apply_boundary(){
//     double shift_dummy;
    size_t ncol =  this->ps->position->tda;
    for (int ip = 0; ip < this->ps->Np; ip++) {
        for (int d = 0; d < this->ps->sdim; d++) {
//             ps->position->data[ip * ncol  + d]  = modf((this->ps->position->data[ip * ncol  + d]+this->L[d]) /  this->L[d], &shift_dummy) * this->L[d];
            ps->position->data[ip * ncol  + d] -= .5 * this->L[d];
            ps->position->data[ip * ncol  + d] -= this->L[d] * round( ps->position->data[ip * ncol  + d] / this->L[d] );
            ps->position->data[ip * ncol  + d] += .5 * this->L[d];
        }
    }
}
