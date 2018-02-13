/**************************************************
 * this is for weighted graph in ChemNetwork      *
 * Tiecheng Zhou                                  *
 * A. Clark Research Lab, Department of Chemistry *
 * Washington State University, Pullman, WA 99164 *
 **************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#include<math.h>
#include<time.h>

#include "chemnetworks_orig.h"

using namespace CN_NS;

ChemNetworkOrig::ChemNetwork_Weighted_Graph::ChemNetwork_Weighted_Graph(){   // constructor for weighted graph
  index_weighted_graph = 0; // setting the default values here
  index_wg_st1_cluster = 0;
  for(wgi=0; wgi<NUM_INTER; wgi++)
  {
    cluster_st1_sv_type[wgi] = cluster_st1_atom[wgi] = cluster_st1_sv_atom[wgi] = 0;
    cluster_st1_sv_dist_min[wgi] = 0.0;
    cluster_st1_sv_dist_max[wgi] = 1000.0;
    atom1_wg_st1_sv1[wgi] = atom2_wg_st1_sv1[wgi] = atom1_wg_st1_sv2[wgi] = atom2_wg_st1_sv2[wgi] = 0;
    atom1_wg_sv1_sv1[wgi] = atom2_wg_sv1_sv1[wgi] = atom1_wg_sv2_sv2[wgi] = atom2_wg_sv2_sv2[wgi] = atom1_wg_sv1_sv2[wgi] = atom2_wg_sv1_sv2[wgi] = 0;
  }
  index_wg_st1_sv1 = index_wg_st1_sv2 = 0;
  num_wg_st1_sv1_dist = num_wg_st1_sv2_dist = 0; 
  funct_type_wg_st1_sv1 = funct_type_wg_st1_sv2 = 0;
  funct_par1_wg_st1_sv1 = funct_par2_wg_st1_sv1 = funct_par1_wg_st1_sv2 = funct_par2_wg_st1_sv2 = 0.0;   // parameters in funct, currently 1.0 / (1.0 + exp(n0 * (r - r0)/r0))
  index_wg_sv1_sv1 = index_wg_sv2_sv2 = index_wg_sv1_sv2 = 0;
  num_wg_sv1_sv1_dist = num_wg_sv2_sv2_dist = num_wg_sv1_sv2_dist = 0;
  funct_type_wg_sv1_sv1 = funct_type_wg_sv2_sv2 = funct_type_wg_sv1_sv2;
  funct_par1_wg_sv1_sv1 = funct_par2_wg_sv1_sv1 = funct_par1_wg_sv2_sv2 = funct_par2_wg_sv2_sv2 = funct_par1_wg_sv1_sv2 = funct_par2_wg_sv1_sv2 = 0.0;
};

ChemNetworkOrig::ChemNetwork_Weighted_Graph::~ChemNetwork_Weighted_Graph(){}; // destructor for weighted graph



/* this is a function to create the Adjacency matrix of weighted graph based on the coordinates from selected cluster */
int ChemNetworkOrig::ChemNetwork_Weighted_Graph::create_WG_Adj_from_cluster(FILE *output_weighted_graph, double **WG_Adj, struct Mol_identity *WG_Mol_id,
	double *atom_cluster_st1, int num_mol_cluster_st1, int nsolute1,
	double *atom_cluster_sv1, int num_mol_cluster_sv1, int nsolvent1,
	double *atom_cluster_sv2, int num_mol_cluster_sv2, int nsolvent2,
	int index_wg_st1_sv1, int num_wg_st1_sv1_dist, int atom1_wg_st1_sv1[NUM_INTER], int atom2_wg_st1_sv1[NUM_INTER], int funct_type_wg_st1_sv1, double funct_par1_wg_st1_sv1, double funct_par2_wg_st1_sv1,
	int index_wg_st1_sv2, int num_wg_st1_sv2_dist, int atom1_wg_st1_sv2[NUM_INTER], int atom2_wg_st1_sv2[NUM_INTER], int funct_type_wg_st1_sv2, double funct_par1_wg_st1_sv2, double funct_par2_wg_st1_sv2,
	int index_wg_sv1_sv1, int num_wg_sv1_sv1_dist, int atom1_wg_sv1_sv1[NUM_INTER], int atom2_wg_sv1_sv1[NUM_INTER], int funct_type_wg_sv1_sv1, double funct_par1_wg_sv1_sv1, double funct_par2_wg_sv1_sv1,
	int index_wg_sv2_sv2, int num_wg_sv2_sv2_dist, int atom1_wg_sv2_sv2[NUM_INTER], int atom2_wg_sv2_sv2[NUM_INTER], int funct_type_wg_sv2_sv2, double funct_par1_wg_sv2_sv2, double funct_par2_wg_sv2_sv2,
	int index_wg_sv1_sv2, int num_wg_sv1_sv2_dist, int atom1_wg_sv1_sv2[NUM_INTER], int atom2_wg_sv1_sv2[NUM_INTER], int funct_type_wg_sv1_sv2, double funct_par1_wg_sv1_sv2, double funct_par2_wg_sv1_sv2,
        double xside, double yside, double zside){

  int nodei, nodej, atomi, atomj, typei, typej, indexi, indexj, totalnum, iter;
  double sitei_x,sitei_y, sitei_z, sitej_x, sitej_y, sitej_z;
  double dist_ij, weight_ij;

  totalnum = num_mol_cluster_st1 + num_mol_cluster_sv1 + num_mol_cluster_sv2;  // currently 1 solute and 2 solvents are implemented

  for(nodei = 0; nodei < totalnum; nodei++)
  {
    if(WG_Mol_id[nodei].solute_type == 1 && WG_Mol_id[nodei].solvent_type == 0)
      typei = 1;

    if(WG_Mol_id[nodei].solute_type == 0 && WG_Mol_id[nodei].solvent_type == 1)
      typei = 2;

    if(WG_Mol_id[nodei].solute_type == 0 && WG_Mol_id[nodei].solvent_type == 2)
      typei = 3;

    for(nodej = 0; nodej < totalnum; nodej++)
    {
      if(nodej == nodei)
        WG_Adj[nodei][nodej] = 0.0;
      else
      { 
        if(WG_Mol_id[nodej].solute_type == 1 && WG_Mol_id[nodej].solvent_type == 0)
          typej = 1;

        if(WG_Mol_id[nodej].solute_type == 0 && WG_Mol_id[nodej].solvent_type == 1)
          typej = 2;

        if(WG_Mol_id[nodej].solute_type == 0 && WG_Mol_id[nodej].solvent_type == 2)
          typej = 3;


        /* get the weights for the requested edges*/
        if(typei == 1)
        {
          if(typej == 1)
            WG_Adj[nodei][nodej] = 0.0;   // no solute1 - solute1 interaction is included
          else if(typej == 2)
          {
            if(index_wg_st1_sv1 == 1)  // get edge weights when solute1 - solvent1 interaction is requested
            {
              weight_ij = 0;
              for(iter=0; iter < num_wg_st1_sv1_dist; iter++)  // the weight between nodei and nodej is additive over all the possible interaction sites
              {
                indexi = nodei;
                indexj = nodej - num_mol_cluster_st1;

                dist_ij = wg_site_distance(atom_cluster_st1, indexi, nsolute1, atom1_wg_st1_sv1[iter], atom_cluster_sv1, indexj, nsolvent1, atom2_wg_st1_sv1[iter], xside, yside, zside);
                if(funct_type_wg_st1_sv1 == 1)
                  weight_ij = weight_ij + 1.0 / (1.0 + exp( funct_par1_wg_st1_sv1 * (dist_ij - funct_par2_wg_st1_sv1) / funct_par2_wg_st1_sv1 ) );  // fermi-function, sigmoid-type
                else {                                                                                                      /* add other functional formulisms here */
                  printf("warning in solute1-solvent1 edge: weight-function-type %d is not implemented\n",funct_type_wg_st1_sv1);
                  weight_ij = 0.0;
                }
              }

              WG_Adj[nodei][nodej] = weight_ij;
            }
            else
              WG_Adj[nodei][nodej] = 0.0;

          }
          else if(typej == 3)
          {
            if(index_wg_st1_sv2 == 1)
            {
              weight_ij = 0.0;
              for(iter=0; iter < num_wg_st1_sv2_dist; iter++)
              {
                indexi = nodei;
                indexj = nodej - num_mol_cluster_st1 - num_mol_cluster_sv1;

                dist_ij = wg_site_distance(atom_cluster_st1, indexi, nsolute1, atom1_wg_st1_sv1[iter], atom_cluster_sv2, indexj, nsolvent2, atom2_wg_st1_sv2[iter], xside, yside, zside);
                if(funct_type_wg_st1_sv2 == 1)
                  weight_ij = weight_ij + 1.0 / (1.0 + exp( funct_par1_wg_st1_sv2 * (dist_ij - funct_par2_wg_st1_sv2) / funct_par2_wg_st1_sv2 ) );
                else {
                  printf("warning in solute1-solvent2 edge: weight-function-type %d is not implemented\n",funct_type_wg_st1_sv2);
                  weight_ij = 0.0;
                }
              }
              WG_Adj[nodei][nodej] = weight_ij;
            }
            else
              WG_Adj[nodei][nodej] = 0.0;

          }
        }
        else if (typei == 2)
        {
          if(typej == 1)
          {
            if(index_wg_st1_sv1 == 1)
            {
              weight_ij = 0.0;
              for(iter=0; iter < num_wg_st1_sv1_dist; iter++)
              {
                indexi = nodei - num_mol_cluster_st1;
                indexj = nodej;

                dist_ij = wg_site_distance(atom_cluster_sv1, indexi, nsolvent1, atom2_wg_st1_sv1[iter], atom_cluster_st1, indexj, nsolute1, atom1_wg_st1_sv1[iter], xside, yside, zside);
                if(funct_type_wg_st1_sv1 == 1)
                  weight_ij = weight_ij + 1.0 / (1.0 + exp( funct_par1_wg_st1_sv1 * (dist_ij - funct_par2_wg_st1_sv1) / funct_par2_wg_st1_sv1 ) );
                else {
                  printf("warning in solute1-solvent1 edge: weight-function-type %d is not implemented\n", funct_type_wg_st1_sv1);
                  weight_ij = 0.0;
                }
              }
              WG_Adj[nodei][nodej] = weight_ij;
            }
            else
              WG_Adj[nodei][nodej] = 0.0;
          }
          else if(typej == 2)   // this is solvent1 - solvent1 interaction
          {
            if(index_wg_sv1_sv1 == 1)
            {
              weight_ij = 0.0;
              for(iter=0; iter < num_wg_sv1_sv1_dist; iter++)
              {
                indexi = nodei - num_mol_cluster_st1;
                indexj = nodej - num_mol_cluster_st1;

                dist_ij = wg_site_distance(atom_cluster_sv1, indexi, nsolvent1, atom1_wg_sv1_sv1[iter], atom_cluster_sv1, indexj, nsolvent1, atom2_wg_sv1_sv1[iter], xside, yside, zside);
                if(funct_type_wg_sv1_sv1 == 1)
                  weight_ij = weight_ij + 1.0 / (1.0 + exp( funct_par1_wg_sv1_sv1 * (dist_ij - funct_par2_wg_sv1_sv1) / funct_par2_wg_sv1_sv1 ) );
                else {
                  printf("warning in solvent1-solvent1 edge: weight-function-type %d is not implemented\n", funct_type_wg_sv1_sv1);
                  weight_ij = 0.0;
                }
              }
              WG_Adj[nodei][nodej] = weight_ij;
            }
            else
              WG_Adj[nodei][nodej] = 0.0;
          }
          else if(typej == 3)
          {
            if(index_wg_sv1_sv2 == 1) // this is solvent1 - solvent2 interaction
            {
              weight_ij = 0.0;
              for(iter=0; iter < num_wg_sv1_sv1_dist; iter++)
              {
                indexi = nodei - num_mol_cluster_st1;
                indexj = nodej - num_mol_cluster_st1 - num_mol_cluster_sv1;

                dist_ij = wg_site_distance(atom_cluster_sv1, indexi, nsolvent1, atom1_wg_sv1_sv1[iter], atom_cluster_sv2, indexj, nsolvent2, atom2_wg_sv1_sv1[iter], xside, yside, zside);
                if(funct_type_wg_sv1_sv2 == 1)
                  weight_ij = weight_ij + 1.0 / (1.0 + exp( funct_par1_wg_sv1_sv2 * (dist_ij - funct_par2_wg_sv1_sv2) / funct_par2_wg_sv1_sv2 ) );
                else {
                  printf("warning in solvent1-solvent2 edge: weight-function-type %d is not implemented\n", funct_type_wg_sv1_sv2);
                  weight_ij = 0.0;
                }
              }
              WG_Adj[nodei][nodej] = weight_ij;
            }
            else
              WG_Adj[nodei][nodej] = 0.0;
          }

        }
        else if (typei == 3)
        {
          if(typej == 1)
          {
            if(index_wg_st1_sv2 == 1)
            {
              weight_ij = 0.0;
              for(iter=0; iter < num_wg_st1_sv2_dist; iter++)
              {
                indexi = nodei - num_mol_cluster_st1 - num_mol_cluster_sv1;
                indexj = nodej;

                dist_ij = wg_site_distance(atom_cluster_sv2, indexi, nsolvent2, atom2_wg_st1_sv2[iter], atom_cluster_st1, indexj, nsolute1, atom1_wg_st1_sv2[iter], xside, yside, zside);
                if(funct_type_wg_st1_sv2 == 1)
                  weight_ij = weight_ij + 1.0 / (1.0 + exp( funct_par1_wg_st1_sv2 * (dist_ij - funct_par2_wg_st1_sv2) / funct_par2_wg_st1_sv2 ) );
                else {
                  printf("warning in solute1-solvent2 edge: weight-function-type %d is not implemented\n", funct_type_wg_st1_sv2);
                  weight_ij = 0.0;
                }
              }
              WG_Adj[nodei][nodej] = weight_ij;
            }
            else
              WG_Adj[nodei][nodej] = 0.0;
          }
          else if(typej == 2)
          {
            if(index_wg_sv1_sv2 == 1)
            {
              weight_ij = 0.0;
              for(iter=0; iter < num_wg_sv1_sv2_dist; iter++)
              {
                indexi = nodei - num_mol_cluster_st1 - num_mol_cluster_sv1;
                indexj = nodej - num_mol_cluster_st1;

                dist_ij = wg_site_distance(atom_cluster_sv2, indexi, nsolvent2, atom2_wg_sv1_sv2[iter], atom_cluster_sv1, indexj, nsolvent1, atom1_wg_sv1_sv2[iter], xside, yside, zside);
                if(funct_type_wg_sv1_sv2 == 1)
                  weight_ij = weight_ij + 1.0 / (1.0 + exp( funct_par1_wg_sv1_sv2 * (dist_ij - funct_par2_wg_sv1_sv2) / funct_par2_wg_sv1_sv2 ) );
                else {
                  printf("warning in solvent1-solvent2 edge: weight-function-type %d is not implemented\n", funct_type_wg_sv1_sv2);
                  weight_ij = 0.0;
                }
              }
              WG_Adj[nodei][nodej] = weight_ij;
            }
            else
              WG_Adj[nodei][nodej] = 0.0;
          }
          else if(typej == 3)
          {
            if(index_wg_sv2_sv2 == 1)
            {
              weight_ij = 0.0;
              for(iter=0; iter < num_wg_sv2_sv2_dist; iter++)
              {
                indexi = nodei - num_mol_cluster_st1 - num_mol_cluster_sv1;
                indexj = nodej - num_mol_cluster_st1 - num_mol_cluster_sv1;

                dist_ij = wg_site_distance(atom_cluster_sv2, indexi, nsolvent2, atom1_wg_sv2_sv2[iter], atom_cluster_sv2, indexj, nsolvent2, atom2_wg_sv2_sv2[iter], xside, yside, zside);
                if(funct_type_wg_sv2_sv2 == 1)
                  weight_ij = weight_ij + 1.0 / (1.0 + exp( funct_par1_wg_sv2_sv2 * (dist_ij - funct_par2_wg_sv2_sv2) / funct_par2_wg_sv2_sv2 ) );
                else {
                  printf("warning in solvent2-solvent2 edge: weight-function-type %d is not implemented\n", funct_type_wg_sv2_sv2);
                  weight_ij = 0.0;
                }
              }
              WG_Adj[nodei][nodej] = weight_ij;
            }
            else
              WG_Adj[nodei][nodej] = 0.0;
          }

        }  // end of ' if(typei == 1) else if == 2 else if == 3 '

      } // end of ' if(nodej == nodei) else '

    }  // end of ' for(nodej = 0; nodej < totalnum; nodej++) '
  } // end of ' for(nodei = 0; nodei < totalnum; nodei++) '
  

  /* print out the Adjacency matrix to an output-file */
  for(nodei=0; nodei < totalnum; nodei++)
  {
    for(nodej=0; nodej < totalnum; nodej++)
    {
      if(WG_Adj[nodei][nodej] > 0.0) // only print the non-zero weights
      {
        if(WG_Mol_id[nodei].solute_type == 1)
        {
          if(WG_Mol_id[nodej].solute_type == 1)
            fprintf(output_weighted_graph, "%d ( st%d %d ) - %d ( st%d %d ) edge-weight: %lf \n", nodei, WG_Mol_id[nodei].solute_type, WG_Mol_id[nodei].id, nodej, WG_Mol_id[nodej].solute_type, WG_Mol_id[nodej].id, WG_Adj[nodei][nodej]);
          else if(WG_Mol_id[nodej].solute_type == 0)
            fprintf(output_weighted_graph, "%d ( st%d %d ) - %d ( sv%d %d ) edge-weight: %lf \n", nodei, WG_Mol_id[nodei].solute_type, WG_Mol_id[nodei].id, nodej, WG_Mol_id[nodej].solvent_type, WG_Mol_id[nodej].id, WG_Adj[nodei][nodej]);

        }
        else if(WG_Mol_id[nodei].solute_type == 0)
        {
          if(WG_Mol_id[nodej].solute_type == 1)
            fprintf(output_weighted_graph, "%d ( sv%d %d ) - %d ( st%d %d ) edge-weight: %lf \n", nodei, WG_Mol_id[nodei].solvent_type, WG_Mol_id[nodei].id, nodej, WG_Mol_id[nodej].solute_type, WG_Mol_id[nodej].id, WG_Adj[nodei][nodej]);
          else if(WG_Mol_id[nodej].solute_type == 0)
            fprintf(output_weighted_graph, "%d ( sv%d %d ) - %d ( sv%d %d ) edge-weight: %lf \n", nodei, WG_Mol_id[nodei].solvent_type, WG_Mol_id[nodei].id, nodej, WG_Mol_id[nodej].solvent_type, WG_Mol_id[nodej].id, WG_Adj[nodei][nodej]);

        }
      }

    }
  }


/*   below are used in debug
  fprintf(output_weighted_graph,"%d %d\n",num_mol_cluster_st1,nsolute1);
  fprintf(output_weighted_graph,"%d %d\n",num_mol_cluster_sv1,nsolvent1);
  fprintf(output_weighted_graph,"%d %d\n",num_mol_cluster_sv2,nsolvent2);

  fprintf(output_weighted_graph,"%d %d %d %lf %lf\n",index_wg_st1_sv1,num_wg_st1_sv1,funct_type_wg_st1_sv1,funct_par1_wg_st1_sv1,funct_par2_wg_st1_sv1);
  fprintf(output_weighted_graph,"%d %d %d %lf %lf\n",index_wg_st1_sv2,num_wg_st1_sv2,funct_type_wg_st1_sv2,funct_par1_wg_st1_sv2,funct_par2_wg_st1_sv2);
  fprintf(output_weighted_graph,"%d %d %d %lf %lf\n",index_wg_sv1_sv1,num_wg_sv1_sv1,funct_type_wg_sv1_sv1,funct_par1_wg_sv1_sv1,funct_par2_wg_sv1_sv1);
  fprintf(output_weighted_graph,"%d %d %d %lf %lf\n",index_wg_sv2_sv2,num_wg_sv2_sv2,funct_type_wg_sv2_sv2,funct_par1_wg_sv2_sv2,funct_par2_wg_sv2_sv2);
  fprintf(output_weighted_graph,"%d %d %d %lf %lf\n",index_wg_sv1_sv2,num_wg_sv1_sv2,funct_type_wg_sv1_sv2,funct_par1_wg_sv1_sv2,funct_par2_wg_sv1_sv2);
*/

  return(0);
};





double ChemNetworkOrig::ChemNetwork_Weighted_Graph::wg_site_distance(double *atomM1, int idmolM1, int natmM1, int idatmM1, double *atomM2, int idmolM2, int natmM2, int idatmM2, double xside, double yside, double zside)
{
  double wg_distance = 0.0;
  double site1x, site1y, site1z, site2x, site2y, site2z;
  double minx, miny, minz;

  site1x = atomM1[ (idmolM1 * natmM1 + (idatmM1-1)) * 3 ];
  site1y = atomM1[ (idmolM1 * natmM1 + (idatmM1-1)) * 3 + 1 ];
  site1z = atomM1[ (idmolM1 * natmM1 + (idatmM1-1)) * 3 + 2 ];

  site2x = atomM2[ (idmolM2 * natmM2 + (idatmM2-1)) * 3 ];
  site2y = atomM2[ (idmolM2 * natmM2 + (idatmM2-1)) * 3 + 1 ];
  site2z = atomM2[ (idmolM2 * natmM2 + (idatmM2-1)) * 3 + 2 ];
  
  minx = site1x - site2x;
  miny = site1y - site2y;
  minz = site1z - site2z;

  if(minx > xside*0.5){   // minx is the shortest distance in x-direction, considering pbc
    while(fabs(minx) > xside*0.5){   // if the atom sites are separated by more than 1 pbc-unit, this may happen in CP2K output-file
      minx = minx - xside;
    }
  }
  else if(minx < -xside*0.5){
    while(fabs(minx) > xside*0.5){
      minx = minx + xside;
    }
  }

  if(miny > yside*0.5){
    while(fabs(miny) > yside*0.5){
      miny = miny - yside;
    }
  }
  else if(miny < -yside*0.5){
    while(fabs(miny) > yside*0.5){
      miny = miny + yside;
    }
  }

  if(minz > zside*0.5){
    while(fabs(minz) > zside*0.5){
      minz = minz - zside;
    }
  }
  else if(minz < -zside*0.5){
    while(fabs(minz) > zside*0.5){
      minz = minz + zside;
    }
  }

  wg_distance = sqrt( minx * minx + miny * miny + minz * minz);

  return(wg_distance);
}

int ChemNetworkOrig::ChemNetwork_Weighted_Graph::findsolvent(int id_solvent_type, int number, int array[NUM_INTER]){ // check an solvent-id from an array
  int result = 0;
  int i;

  if(number > NUM_INTER){
    printf("Error: size-overflow in function 'findsolvent'\n");
    result = -1;
  }
  else {
    for(i=0; i<number; i++){
      if(id_solvent_type == array[i]){
        result = 1;
        break;
      }
    }
  }

  return(result);
}



