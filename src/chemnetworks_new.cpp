#include <stdlib.h>
#include "chemnetworks_new.h"

using namespace CN_NS;

/* -------------------------------------------------------------------------- */

ChemNetworkNew::ChemNetworkNew()
{

}

/* -------------------------------------------------------------------------- */

void ChemNetworkNew::deallocate()
{
  if(!using_new_interface) return;

  deallocate_solvent(nsolvent1, slvntatm1, nAtomS1, atmTypeS1, atmS1);
}

void ChemNetworkNew::deallocate_solvent(int n1, char ** t1, int n2, char ** t2, double * xyz)
{
  for(int i=0; i<n1; ++i) free(t1[i]);
  free(t1);
  
  for(int i=0; i<n2; ++i) free(t2[i]);
  free(t2);
  
  free(xyz);
}
