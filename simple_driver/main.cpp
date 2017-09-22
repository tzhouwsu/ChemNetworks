// Compile: make

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "chemnetworks_new.h"

int main(int argc, char* argv[]) 
{
  fprintf(stdout,"\nHello, I'm a simple driver\n");

  // create a new instance of ChemNetworks

  CN_NS::ChemNetwork * cn = new CN_NS::ChemNetworkNew();

  // pass command-line arguments to ChemNetworks and execute

  fprintf(stdout,"\n  ChemNetworks is reading input file...\n");

  cn->input_read(argc, argv);
  
  fprintf(stdout,"\n  ChemNetworks is processing...\n");

  cn->input_process(argc, argv);

  fprintf(stdout,"\n  Finished.\n");

  // clean up

  delete cn;  

  fprintf(stdout,"\nBye!!\n");
}
