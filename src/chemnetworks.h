/*************************************************
 * chemnetworks.h                                *
 *                                               *
 * Author: Abdullah Ozkanlar                     *
 *         abdullah.ozkanlar@wsu.edu             *
 *                                               *
 * A. Clark Research Lab, Chemistry Department   *
 * Washington State University, Pullman/WA 99164 *
 *************************************************/

#ifndef CHEMNETWORKS_H
#define CHEMNETWORKS_H

#include <stdio.h>

namespace CN_NS {
  
  class ChemNetwork {
    
  public : 
    ChemNetwork() {};
    virtual ~ChemNetwork() {};
    
    virtual void input_read(int, char * []) {};
    virtual void input_process(int, char * []) {};
    
    virtual void read_xyz_all(char *[]) {};
    virtual void read_xyz_solvent1(char[]) {};
    virtual void load_xyz_solvent1(int, char**, double*) {};
    
    virtual void process_config(char * []) {};

    virtual void deallocate() {};
  };
}
    
#endif // CHEMNETWORKS_H
