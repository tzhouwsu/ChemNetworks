// Compile: make

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "chemnetworks_new.h"

#define CN_ELEMENT_SIZE 32

void read_xyz(char * filename, int &n, char ** &t, double * &xyz)
{
  // only allocate memory on first time reading a file
  
  int first_config = 0;
  if(n == 0) first_config = 1;
  
  // we just copy copy from chemnetworks to be lazy

  FILE * finput = NULL;
  if((finput = fopen(filename,"r")) == NULL){
    printf("read_xyz()::Cannot open filename= %s\n", filename);
    exit(1);
  }

  if(fscanf(finput,"%d", &n) != 1){
    printf("Cannot read number of atoms in XYZ file\n");
    exit(1);
  }

  // If first config, then allocate memory
  //  assumes constant number of particles during simulation
  
  if(first_config) {
    xyz = (double*) malloc(3*n*sizeof(double));
  
    t = (char **) malloc(n*sizeof(char *));
    for(int i=0; i < n; ++i) t[i] = (char *) malloc(CN_ELEMENT_SIZE * sizeof(char));
  }
  
  // Skip single header line
  char line[CN_ELEMENT_SIZE];
  fgets(line,CN_ELEMENT_SIZE,finput);

  // read coordinates
  
  for(int i=0; i<n; ++i)
    if(fscanf(finput, "%31s %lf %lf %lf", t[i], &xyz[3*i], &xyz[3*i+1], &xyz[3*i+2]) != 4){
      printf("Solvent1 XYZ file: Cannot read coordinates\n");
      exit(1);
    }

  fclose(finput);
}

int main(int argc, char* argv[]) 
{
  fprintf(stdout,"\nHello, I'm a simple driver\n");

  // create a new instance of ChemNetworks

  fprintf(stdout,"  Driver code is creating instance of ChemNetworks\n");
  
  CN_NS::ChemNetwork * cn = new CN_NS::ChemNetworkNew();

  // pass command-line arguments to ChemNetworks and execute

  fprintf(stdout,"\n  ChemNetworks is reading input file...\n");

  cn->input_read(argc, argv);

  // read single xyz config

  int num_atoms = 0;
  char ** type = NULL;
  double * x = NULL;

  for(int i=0; i<5; ++i) {
  
    fprintf(stdout,"\n  Driver code is reading config from file... i= %i\n",i);
    
    read_xyz(argv[2], num_atoms, type, x);
  
    fprintf(stdout,"\n  Driver code is sending ChemNetworks the coordinates...\n");
  
    cn->load_xyz_solvent1(num_atoms, type, x);
  
    fprintf(stdout,"\n  ChemNetworks is processing...\n");

    cn->process_config(argv);

    fprintf(stdout,"\n  Finished.\n");

  }

  // clean up

  fprintf(stdout,"  Driver code is tearing down ChemNetworks\n");
  
  delete cn;
  
  for(int i=0; i<num_atoms; ++i) free(type[i]);
  free(type);
  free(x);

  fprintf(stdout,"\nBye!!\n");
}
