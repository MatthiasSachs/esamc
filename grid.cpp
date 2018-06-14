//
//  grid.cpp
//  linkedCellMD
//
//  Created by Matthias Sachs on 26/12/2017.
//  Copyright © 2017 Matthias Sachs. All rights reserved.
//

#include "grid.hpp"

#define sqr(x) ((x)*(x))
#define VERIFY(a) do { if((a)<0) { fprintf(stderr,"Failure line in file grid.cpp %d.\n",__LINE__); exit(-1);}}while(0)
#define GSL_VERIFY(a) do { if((a) != 0) { fprintf(stderr,"Failure line in file grid.cpp %d.\n",__LINE__); exit(-1);}}while(0)

inline int positive_modulo(int i, int n) {
    return (i % n + n) % n;
}

ParticleSystem::ParticleSystem(int Np_a, int sdim_a, Domain *domain_a){
    Np     =   Np_a;
    sdim   =   sdim_a;
    dim    =   this->Np * this->sdim;
    domain     =   domain_a; this->domain->ps = this;
    position    = gsl_matrix_calloc(this->Np, this->sdim);
    momentum    = gsl_matrix_calloc(this->Np, this->sdim);
    force       = gsl_matrix_calloc(this->Np, this->sdim);
    mass        = gsl_vector_calloc(this->Np); gsl_vector_set_all(mass, 1.0);
    U           = gsl_vector_calloc(this->Np);
    
    position_as_vec = gsl_vector_view_array(position->data, dim);
    momentum_as_vec = gsl_vector_view_array(momentum->data, dim);
    force_as_vec = gsl_vector_view_array(force->data, dim);
    //this->position   =
    //this->momentum   =
    rel_distance   = gsl_spmatrix_alloc(this->Np, this->Np);
    rel_position   = gsl_spmatrix_alloc(this->dim, this->dim);
    std::cout << "Object is being created" << std::endl;
};

void ParticleSystem::addPotential( PairPotential *potential){
    this->ppotential = potential;
}
void ParticleSystem::addPotential( ExternalPotential *epotential){
    this->epotential = epotential;
}
void ParticleSystem::comp_rel_position(double *rel_pos, double *pos1, double *pos2){
    this->domain->comp_rel_position(rel_pos, pos1, pos2);
}
double ParticleSystem::comp_rel_position(double *pos1, double *pos2, int d){
    return this->domain->comp_rel_position(pos1, pos2, d);
}
void ParticleSystem::apply_boundary(){
    this->domain->apply_boundary();
}
Domain::Domain(int sdim_a){
    this->sdim = sdim_a;
}

void Domain::apply_boundary(){};
void Domain::comp_rel_position(double *rel_pos, double *pos1, double *pos2){};
double Domain::comp_rel_position(double *pos1, double *pos2, int d){ return 0;};




Particle::Particle(double *position_a, double *momentum_a, double *mass_a, double *force_a, double *U_a, int id_a){
    position = position_a;
    momentum = momentum_a;
    mass = mass_a;
    id = id_a;
    force = force_a;
    U = U_a;
    
}
ParticleList::ParticleList(Particle *p_a){
    p = p_a;
    next = NULL;
}
ParticleList::ParticleList(double *position, double *momentum, double *mass, double *force, double *U, int id){
    p = new Particle(position, momentum, mass, force, U, id);
}

LcGrid::LcGrid(ParticleSystem *ps_a,int *nc_a){
    ps = ps_a;
    sdim = this->ps->sdim;
    L = this->ps->domain->L;
    nc_all = 1;
    nc = (int*) malloc(this->sdim * sizeof(int));
    for (int i =0; i<this->sdim; i++) {
        this->nc[i] = nc_a[i] > 3 ? nc_a[i] : 3;
        this->nc_all *=this->nc[i];
    }
    this->cells = (Cell*) malloc(this->nc_all*sizeof(Cell));
    for (int i=0; i < this->nc_all; i++){
        this->cells[i]= NULL;
    }
    stamp = new Stamp(this->sdim);
    
    ps->lcgrid = this;
}

void LcGrid::insertList(ParticleList **root_list, ParticleList *i)
{
    i->next = *root_list;
    *root_list = i;
}
void LcGrid::deleteList(ParticleList **q) {
    *q = (*q)->next;
}

int LcGrid::sub2lin(int *ic){
    switch (this->sdim) {
        case 1:
            return ic[0];
            break;
        case 2:
            return ic[0] + this->nc[0]*ic[1];
            break;
        case 3:
            return ic[0] + this->nc[0]*(ic[1] + this->nc[1]*ic[2]);
            break;
        default:
            return -1;
            break;
    }
};

void LcGrid::lin2sub(int *ic, int lin_index){
    for (int i = 0; i<this->sdim; i++) {
        ic[i] = lin_index  % this->nc[i];
        lin_index/= this->nc[i];
    }
};


void LcGrid::gridIndex(int *ic){
    for (int d=0; d<this->sdim; d++){
        ic[d] = positive_modulo(ic[d], this->nc[d]);
    }
}

void LcGrid::findIndex( double* position, int *index ){
    for (int d=0; d<this->sdim; d++){
        index[d] = ((int) floor(position[d] * this->nc[d] / L[d]) ) % nc[d];
    }
}

void LcGrid::addParticleList(ParticleList *pl){
    int ic[3];
    this->findIndex(pl->p->position, ic);
    this->insertList(&(this->cells[this->sub2lin(ic)]), pl);
}

void LcGrid::sortParticles(){
    int ic[3] = {0,0,0};
    int kc[3] = {0,0,0};
    for (int lin_i = 0; lin_i < this->nc_all; lin_i++){
        this->lin2sub(ic, lin_i);
        ParticleList **q = &(this->cells[this->sub2lin(ic)]);
        ParticleList *i = *q;
        while (NULL != i) {
            for (int d=0; d<this->sdim; d++)
                // pointer to predecessor
                kc[d] = (int) floor(i->p->position[d] * this->nc[d] / this->L[d]);
            if (ic[0]!=kc[0] || ic[1]!=kc[1] || ic[2]!=kc[2] )
            {
                deleteList(q);
                insertList(&(this->cells[this->sub2lin(kc)]), i);
            }
            else{
                q = &i->next;
            }
            i = *q;
        }
    }
}

void LcGrid::compRelPos(gsl_spmatrix *rel_position,gsl_spmatrix *rel_distance, double cutoff){
    int ic[3];
    int kc[3];
    double r2;
    double relpos[3];
    gsl_spmatrix_set_zero(rel_distance);
    gsl_spmatrix_set_zero(rel_position);
    for (int lin_i = 0; lin_i < this->nc_all; lin_i++)
    {
        this->lin2sub(ic, lin_i);
        for (ParticleList *i=this->cells[this->sub2lin(ic)]; NULL!=i; i=i->next)
        {
            for (int li = this->stamp->min_index; li < this->stamp->max_index; li++)
            {
                this->stamp->getNeighbour(kc, ic, li);
                this->gridIndex(kc);
                for (ParticleList *j=this->cells[this->sub2lin(kc)]; NULL!=j; j=j->next)
                {
                    if (i->p->id != j->p->id) {
                        for (int d=0; d<this->sdim; d++){
                            relpos[d] = this->ps->comp_rel_position(i->p->position, j->p->position, d);
                        }
                        r2 = 0;
                        for (int d=0; d<this->sdim; d++){
                            r2+= relpos[d]*relpos[d];
                        }
                        if (r2 < cutoff*cutoff) {
                            gsl_spmatrix_set(rel_distance,
                                             i->p->id,
                                             j->p->id,
                                             sqrt(r2));
                            for (int d=0; d<this->sdim; d++){
                                gsl_spmatrix_set(rel_position,
                                                 this->sdim * i->p->id + d,
                                                 this->sdim * j->p->id + d,
                                                 relpos[d]);
                            }
                        }
                    }
                }
            }
        }
    }
}


void LcGrid::compForce(){
    int ic[3];
    int kc[3];
    gsl_matrix_set_zero(this->ps->force);
    
    if( NULL != this->ps->epotential){
        for (int lin_i = 0; lin_i < this->nc_all; lin_i++){
            this->lin2sub(ic, lin_i);
            for (ParticleList *i=this->cells[this->sub2lin(ic)]; NULL!=i; i=i->next)
                this->ps->epotential->comp_force(i->p);
        }
    }
    
    if( NULL != this->ps->ppotential){
        for (int lin_i = 0; lin_i < this->nc_all; lin_i++)
        {
            this->lin2sub(ic, lin_i);
            for (ParticleList *i=this->cells[this->sub2lin(ic)]; NULL!=i; i=i->next)
            {
                for (int li = this->stamp->min_index; li < this->stamp->max_index; li++)
                {
                    this->stamp->getNeighbour(kc, ic, li);
                    this->gridIndex(kc);
                    for (ParticleList *j=this->cells[this->sub2lin(kc)]; NULL!=j; j=j->next)
                    {
                        if (i->p->id != j->p->id) {
   
                            this->ps->ppotential->comp_force(i->p, j->p);
                        }
                    }
                }
            }
        }
    }
}

void LcGrid::compPotential(){
    int ic[3];
    int kc[3];
    gsl_vector_set_zero(this->ps->U);
    
    if( NULL != this->ps->epotential){
        for (int lin_i = 0; lin_i < this->nc_all; lin_i++){
            this->lin2sub(ic, lin_i);
            for (ParticleList *i=this->cells[this->sub2lin(ic)]; NULL!=i; i=i->next)
                this->ps->epotential->comp_pot(i->p);
        }
    }
    if( NULL != this->ps->ppotential){
        for (int lin_i = 0; lin_i < this->nc_all; lin_i++)
        {
            this->lin2sub(ic, lin_i);
            for (ParticleList *i=this->cells[this->sub2lin(ic)]; NULL!=i; i=i->next)
            {
                for (int li = this->stamp->min_index; li < this->stamp->max_index; li++)
                {
                    this->stamp->getNeighbour(kc, ic, li);
                    this->gridIndex(kc);
                    for (ParticleList *j=this->cells[this->sub2lin(kc)]; NULL!=j; j=j->next)
                    {
                        if (i->p->id != j->p->id)  {
                            this->ps->ppotential->comp_pot(i->p, j->p);
                        }
                    }
                }
            }
            
        }
    }
}
Stamp::Stamp(int sdim_a){
    sdim = sdim_a;
    min_index = 0;
    max_index = (int) pow(3, this->sdim);
}
void Stamp::getNeighbour(int *kc, const int *ic, const int lin_index){
    int dim1[3][1] = {
        {0},{-1},{1}
    };
    int dim2[9][2] = {
        {0,0},{0,1},{0,-1},{1,0},{1,1},{1,-1},{-1,0},{-1,1},{-1,-1}
    };
    int dim3[27][3] = {
                {-1,0,0},{-1,0,1},{-1,0,-1},{-1,1,0},{-1,1,1},{-1,1,-1},{-1,-1,0},{-1,-1,1},{-1,-1,-1},
                {0,0,0},{0,0,1},{0,0,-1},{0,1,0},{0,1,1},{0,1,-1},{0,-1,0},{0,-1,1},{0,-1,-1},
                {1,0,0},{1,0,1},{1,0,-1},{1,1,0},{1,1,1},{1,1,-1},{1,-1,0},{1,-1,1},{1,-1,-1}
    };
    
    switch (this->sdim) {
        case 1:
            for (int d = 0; d < 1; d++) {
                kc[d] = ic[d] + dim1[lin_index][d];
            }
            break;
        case 2:
            for (int d = 0; d < 2; d++) {
                kc[d] = ic[d] + dim2[lin_index][d];
            }
            break;
        case 3:
            for (int d = 0; d < 3; d++) {
                kc[d] = ic[d] + dim3[lin_index][d];
            }
            break;
        default:
            break;
    }
};


ExternalPotential::ExternalPotential(size_t sdim_a, ParticleSystem *ps_a){
    this->sdim = sdim_a;
    this->ps= ps_a;
}
void ExternalPotential::comp_force(Particle *pi){};
void ExternalPotential::comp_pot(Particle *pi){};

PairPotential::PairPotential(size_t sdim_a, ParticleSystem *ps_a){
    this->sdim = sdim_a;
    this->ps= ps_a;
}
void PairPotential::comp_force(Particle *pi, Particle *pj){};
void PairPotential::comp_pot(Particle *pi, Particle *pj){};




OutputTask::OutputTask(gsl_vector *data, std::string variableName_a){
    this->variable = (double *) data->data;
    this->variableName = variableName_a;
    this->size = data->size;
};
OutputTask::OutputTask(gsl_matrix *data, std::string variableName_a){
    this->variable = (double *) data->data;
    this->variableName = variableName_a;
    this->size = data->size1*data->size2;
};
void OutputSheduler::feed(long int t){};




/*
 The following functions are usefull for debuggin purposes
*/

void LcGrid::printState(){
    printf("Grid state: \n");
    int ic[3];
    int kc[3];
    for (int lin_i = 0; lin_i < this->nc_all; lin_i++){
        this->lin2sub(ic, lin_i);
        for (ParticleList *i=this->cells[this->sub2lin(ic)]; NULL!=i; i=i->next)
            printf("CELL: (%d %d), ID %d, position %f, %f \n", ic[0], ic[1], i->p->id, i->p->position[0], i->p->position[1]);
    }
}
void LcGrid::printStateForce(){
    printf("Particles force: \n");
    int ic[3];
    int kc[3];
    for (int lin_i = 0; lin_i < this->nc_all; lin_i++){
        this->lin2sub(ic, lin_i);
        for (ParticleList *i=this->cells[this->sub2lin(ic)]; NULL!=i; i=i->next)
            printf("ID %d, position %f, %f \n", i->p->id, i->p->force[0], i->p->force[1]);
    }
}
void LcGrid::printStateNeighors(){
    printf("Neighbors: \n");
    int ic[3];
    int kc[3];
    for (int lin_i = 0; lin_i < this->nc_all; lin_i++){
        this->lin2sub(ic, lin_i);
        for (ParticleList *i=this->cells[this->sub2lin(ic)]; NULL!=i; i=i->next)
            for (int li = this->stamp->min_index; li < this->stamp->max_index; li++)
            {
                this->stamp->getNeighbour(kc, ic, li);
                this->gridIndex(kc);
                for (ParticleList *j=this->cells[this->sub2lin(kc)]; NULL!=j; j=j->next)
                {
                    printf("ID %d, %d; positions (%f, %f), (%f, %f)\n", i->p->id, j->p->id, i->p->position[0], i->p->position[1],j->p->position[0], j->p->position[1]);
                    
                }
            }
    }
    
}
void LcGrid::printStateNeighorsList(){
    printf("Neighbor List: \n");
    int ic[3];
    int kc[3];
    for (int lin_i = 0; lin_i < this->nc_all; lin_i++){
        this->lin2sub(ic, lin_i);
        for (ParticleList *i=this->cells[this->sub2lin(ic)]; NULL!=i; i=i->next){
            printf("ID %d, {",i->p->id);
            for (int li = this->stamp->min_index; li < this->stamp->max_index; li++)
            {
                this->stamp->getNeighbour(kc, ic, li);
                this->gridIndex(kc);
                for (ParticleList *j=this->cells[this->sub2lin(kc)]; NULL!=j; j=j->next)
                {
                    printf("%d, ", j->p->id);
                    
                }
            }
            printf("} \n");
        }
        
    }
    
}
void LcGrid::printStateNeighorsListDistance(){
    printf("Neighbor List: \n");
    int ic[3];
    int kc[3];
    for (int lin_i = 0; lin_i < this->nc_all; lin_i++){
        this->lin2sub(ic, lin_i);
        for (ParticleList *i=this->cells[this->sub2lin(ic)]; NULL!=i; i=i->next){
            printf("ID %d, {",i->p->id);
            for (int li = this->stamp->min_index; li < this->stamp->max_index; li++)
            {
                this->stamp->getNeighbour(kc, ic, li);
                this->gridIndex(kc);
                for (ParticleList *j=this->cells[this->sub2lin(kc)]; NULL!=j; j=j->next)
                {
                    double dist;
                    double r2 = 0;
                    for (int d=0; d<this->sdim; d++){
                        dist = this->ps->comp_rel_position(i->p->position, j->p->position, d);
                        r2+= dist*dist;
                    }
                    printf("%f, ", sqrt(r2));
                    
                }
            }
            printf("}, \n");
            
        }
        
    }
    
}





