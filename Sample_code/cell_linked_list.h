//
// Created by How Yu on 2022/4/5.
//

#ifndef EC526_CELL_LINKED_LIST_H
#define EC526_CELL_LINKED_LIST_H

#endif //EC526_CELL_LINKED_LIST_H
#include <math.h>
double distance(double ax, double ay, double bx, double by){
    /*
     * This function calculate the distance between 2 particles
     * */
    double lx,ly,d;
    lx=(ax-bx)*(ax-bx);
    ly=(ay-by)*(ay-by);
    d= sqrt(lx+ly);
    return d;
}
void build_cell_list(int num_cell, int num_particle, int dim, double r_cut, double r_skin,double pos[]){
    double head[num_cell];
    double linked_list[8];
    // find the x_range and y_range of each cell
    // find cell for each particle
    // find the neighbor list of each particle
    for(int i=0;i<num_cell;i++){

    }

};
