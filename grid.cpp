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


ParticleSystem::ParticleSystem(int Np_a, Domain *domain_a){
    Np     =   Np_a;                   
    sdim   =   domain_a->sdim;
    dim    =   this->Np * this->sdim;
    domain     =   domain_a; this->domain->ps = this;
    position    = gsl_matrix_calloc(this->Np, this->sdim);
    momentum    = gsl_matrix_calloc(this->Np, this->sdim);
    force       = gsl_matrix_calloc(this->Np, this->sdim);
    mass        = gsl_vector_calloc(this->Np); gsl_vector_set_all(mass, 1.0);
    U           = gsl_vector_calloc(this->Np);
    laplace     = gsl_vector_calloc(this->Np);
    
    position_as_vec = gsl_vector_view_array(position->data, dim);
    momentum_as_vec = gsl_vector_view_array(momentum->data, dim);
    force_as_vec = gsl_vector_view_array(force->data, dim);
    //this->position   =
    //this->momentum   =
    rel_distance   = gsl_spmatrix_alloc(this->Np, this->Np);
    rel_position   = gsl_spmatrix_alloc(this->dim, this->dim);
    
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (this->T);
    
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

void ParticleSystem::initMomentum(double beta){
    for (size_t i = 0; i < this->Np; i++) {
        for (size_t j = 0; j < this->sdim; j++) {
            gsl_matrix_set(this->momentum,i,j, gsl_ran_gaussian (this->r, sqrt(this->mass->data[i]/beta)) );
        }
    }
}
void ParticleSystem::initMomentumDPD(double beta){
    // Initialize momentum according to temperature
    this->initMomentum(beta);
    this->centerMomentum();
    this->scaleMomentum(beta);
}

void ParticleSystem::centerMomentum(){
    // Center momentum in each dimension
    double total_momentum;
    for (size_t j = 0; j < this->sdim; j++) {
        total_momentum=0;
        for (size_t i = 0; i < this->Np; i++) {
            total_momentum += this->momentum->data[i * this->momentum->tda + j];
        }
        total_momentum/=this->Np;
        for (size_t i = 0; i < this->Np; i++) {
            this->momentum->data[i * this->momentum->tda + j] += -total_momentum;
        }
    }
}

void ParticleSystem::scaleMomentum(double beta){
    // Scale momentum //MS commment: mass does not seem to be considered here yet
    double total_KT = 0;
    for (size_t i = 0; i < this->Np; i++) {
        for (size_t j = 0; j < this->sdim; j++) {
            total_KT += this->momentum->data[i * this->momentum->tda + j] * this->momentum->data[i * this->momentum->tda + j];
        }
    }
    total_KT /= 3.0*(this->Np-1);
    double KT_scale = sqrt(1.0/beta/total_KT);
    for (size_t i = 0; i < this->Np; i++) {
        for (size_t j = 0; j < this->sdim; j++) {
            this->momentum->data[i * this->momentum->tda + j] *= KT_scale;
        }
    }
}

Domain::Domain(int sdim_a){
    this->sdim = sdim_a;
}

void Domain::apply_boundary(){};
void Domain::comp_rel_position(double *rel_pos, double *pos1, double *pos2){};
double Domain::comp_rel_position(double *pos1, double *pos2, int d){ return 0;};




Particle::Particle(double *position_a, double *momentum_a, double *mass_a, double *force_a, double *U_a, double *laplace_a, int id_a){
    this->position = position_a;
    this->momentum = momentum_a;
    this->mass = mass_a;
    this->id = id_a;
    this->force = force_a;
    this->U = U_a;
    this->laplace = laplace_a;
    
}
ParticleList::ParticleList(Particle *p_a){
    p = p_a;
    next = NULL;
}
ParticleList::ParticleList(double *position, double *momentum, double *mass, double *force, double *U, double *laplace, int id){
    p = new Particle(position, momentum, mass, force, U, laplace, id);
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
    
    force_timer.start_clock();
    
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
    force_timer.record_time();
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

void LcGrid::compLaplacePotential(){
    int ic[3];
    int kc[3];
    gsl_vector_set_zero(this->ps->laplace);
    
    if( NULL != this->ps->epotential){
        for (int lin_i = 0; lin_i < this->nc_all; lin_i++){
            this->lin2sub(ic, lin_i);
            for (ParticleList *i=this->cells[this->sub2lin(ic)]; NULL!=i; i=i->next)
                this->ps->epotential->comp_laplace(i->p);
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
                            this->ps->ppotential->comp_laplace(i->p, j->p);
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
void ExternalPotential::comp_force(Particle *pi){
    printf("Warning: function comp_force(Particle *pi) not implemented");
};
void ExternalPotential::comp_pot(Particle *pi){
    printf("Warning: function comp_pot(Particle *pi) not implemented");
};
void ExternalPotential::comp_laplace(Particle *pi){
    printf("Warning: function comp_laplace(Particle *pi) not implemented");
};

PairPotential::PairPotential(size_t sdim_a, ParticleSystem *ps_a){
    this->sdim = sdim_a;
    this->ps= ps_a;
}
void PairPotential::comp_force(Particle *pi, Particle *pj){
    printf("Warning: function comp_force(Particle *pi, Particle *pj) not implemented");
};
void PairPotential::comp_pot(Particle *pi, Particle *pj){
    printf("Warning: function comp_pot(Particle *pi, Particle *pj) not implemented");
};
void PairPotential::comp_laplace(Particle *pi, Particle *pj){
    printf("Warning: function comp_laplace(Particle *pi, Particle *pj) not implemented");
};


OutputTask::OutputTask(ParticleSystem *ps, std::string variableName){
    this->variableName = variableName;
    this->ps = ps;
};
void OutputTask::comp_output(double *outputTraj){
    printf("Warning: function comp_output(double *outputTraj) not implemented");
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






