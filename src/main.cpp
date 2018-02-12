/*************************************************
 * ChemNetworks.c                                *
 *                                               *
 * Author: Abdullah Ozkanlar                     *
 *         abdullah.ozkanlar@wsu.edu             *
 *                                               *
 * ChemNetworks version 1.0 (6/27/13)            *
 *                                               *
 * A. Clark Research Lab, Chemistry Department   *
 * Washington State University, Pullman/WA 99164 *
 *************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>

#include "chemnetworks.h" 
#include "chemnetworks_orig.h"
#include "chemnetworks_new.h"

using namespace CN_NS;

int main(int argc, char *argv[])
{
  int using_new_algo = 0;
  for(int i=0; i<argc; ++i)
    if(strcmp(argv[i],"-new") == 0) using_new_algo = 1;

  ChemNetwork * cn;
  
  if(using_new_algo) {
    fprintf(stdout,"Launching ChemNetworkNew\n");
    
    cn = new ChemNetworkNew();

    // read command-line and input file
    
    cn->input_read(argc, argv);
    
    // read all solvent+solute configs (single timestep)
    
    cn->read_xyz_all(argv);
    
    // process config
    
    cn->process_config(argv);
    
    cn->deallocate();

  } else {
    cn = new ChemNetworkOrig();

    // read command-line and input file
    
    cn->input_read(argc, argv);
    
    // process config
    
    cn->input_process(argc, argv);
  }

  delete cn;
}
