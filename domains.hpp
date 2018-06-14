//
//  domains.hpp
//  linkedCellMD
//
//  Created by Matthias Sachs on 13/06/2018.
//  Copyright © 2018 Matthias Sachs. All rights reserved.
//

#ifndef domains_hpp
#define domains_hpp

#include <stdio.h>
#include "grid.hpp"

class Torus : public Domain{
public:
    Torus(int sdim_a, double *L_a);
    void apply_boundary();
    void comp_rel_position(double *rel_pos, double *pos1, double *pos2);
    double comp_rel_position(double *pos1, double *pos2, int d);
};

class RRn : public Domain{
public:
    RRn(int sdim_a, double *L_a);
    void apply_boundary();
    void comp_rel_position(double *rel_pos, double *pos1, double *pos2);
    double comp_rel_position(double *pos1, double *pos2, int d);
};


#endif /* domains_hpp */
