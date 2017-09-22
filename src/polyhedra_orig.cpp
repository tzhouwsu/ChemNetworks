/*************************************************
 * polyhedra.c                                   *
 * Author:  Marissa Masden                       *
 *          marissa.masden@email.wsu.edu         *
 *                                               *
 * A. Clark Research Lab, Chemistry Department   *
 * Washington State University, Pullman/WA 99164 *
 *************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>

#include "chemnetworks_orig.h"

using namespace CN_NS;

void ChemNetworkOrig::polyhedra_st(double *atmS1, double *atmT1, int nd1, int nd2, int nsolvent1, int nsolute1, int nAtomS1, int s1t1cutoffnum, int *s1t1a, int *s1t1b, double *s1t1cutoff, int pbc, double xside, double yside, double zside, FILE *outputPolys, int maxshellsize, double**edgebds, int whichatmT1, int whichatmdexT1)
{
    
    //TO ADD MORE POLYHEDRA CHANGE THIS VALUE:
    int numofpolies = 51;

    //declare new variables
    int *shellist, i,j;
    double **shelldistmtx;
    int **adjacency_matrix;
    int **temp_adjm;
    int **possible_edgelist;
    int poss_edge_num;
    int shellsize, max_num_of_edges;
    int bestpolynum;
    clock_t starttime, endtime;
    double cpu_time_used;
    double damp=.85;
    double tolerance=10E-7;
    poly **polylist;

    starttime = clock();

    //allocate memory...
    
    //initializing to 0 for easy detection of solvation shell size; node numbering starts at 1.
    shellist=(int*)calloc(maxshellsize,sizeof(double));
    
    //create 2d array for distances
    shelldistmtx=(double**)malloc(maxshellsize*sizeof(double*));
    
    for(i=0; i<maxshellsize;i++)
    {
        shelldistmtx[i]=(double*)malloc(maxshellsize*sizeof(double));
    }
    
    //find shell list and distance matrix; is its own function for later use if needed.
    
    shellist_and_distmtx(shellist, shelldistmtx, atmS1, atmT1, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, pbc, xside, yside, zside, &shellsize, maxshellsize, whichatmT1, whichatmdexT1);
    
    //allocate memory for adjacency matrix (and a placeholder temp) now that we know its size.
    adjacency_matrix = (int**)malloc((shellsize+1)*sizeof(int*));
    
    //the last column/row will represent the solute molecule
    
    temp_adjm = (int**)malloc((shellsize+1)*sizeof(int*));
    for(i=0; i<(shellsize+1); i++)
    {
        adjacency_matrix[i]=(int*)malloc((shellsize+1)*sizeof(int));
        temp_adjm[i] = (int*)malloc((shellsize+1)*sizeof(int));
    }
    
    //allocate memory for list of additional edges (complete graph has n(n-1)/2 edges, but last vertex's edges predefined.)
    
    max_num_of_edges = shellsize*(shellsize-1)/2;
        
    possible_edgelist=(int**)malloc(max_num_of_edges* sizeof(int*));
    for(i=0; i<max_num_of_edges; i++)
    {
        possible_edgelist[i] = (int*)malloc(max_num_of_edges * sizeof(int));
    }
    
    poss_edge_num = 0;
    
    //create "base" adjacency matrix
    
    for(i=0; i<shellsize; i++)
    {
        adjacency_matrix[i][i]=0;
        for(j=i+1; j<shellsize; j++)
        {
            //create adjacency matrix entry if smaller than first threshold. 
            if (shelldistmtx[i][j]<=edgebds[shellsize][0]) {
                
                adjacency_matrix[i][j]= 1;
                adjacency_matrix[j][i]= 1;
            }
            else
            {
                //create "possible edge" entry if smaller than second threshold.
                if(shelldistmtx[i][j] <= edgebds[shellsize][1])
                {
                    possible_edgelist[poss_edge_num][0]=i;
                    possible_edgelist[poss_edge_num][1]=j;
                    poss_edge_num=poss_edge_num + 1;
                    
                }
                adjacency_matrix[i][j] = 0;
                adjacency_matrix[j][i] = 0;
            }
        }
        //add "star graph" (solute) edges.
        adjacency_matrix[i][shellsize] = 1;
        adjacency_matrix[shellsize][i] = 1;
    }
    
    //Output adjacency matrix (for troubleshooting)
    /*
     for(i = 0; i<shellsize+1; i++)
     {
     for(j=0; j<shellsize+1; j++)
     {
     printf("%d ", adjacency_matrix[i][j]);
     }
     printf("\n");
     }
     */
    
    //Output "possible" edge list (for troubleshooting)
    
    /*
     for(i=0; i<poss_edge_num; i++)
     {
     printf("%d, %d \n",shellist[possible_edgelist[i][0]], shellist[possible_edgelist[i][1]]);
     }
     */
    
    //generate list of polyhedra to compare.
    polylist=(poly**)malloc(numofpolies*sizeof(poly*));
    makeallpolys(polylist);
    
    
    //the following function loops through adding other possible edges (all possible combinations).
    
    bestpolynum = polymatch(polylist, adjacency_matrix, shelldistmtx, temp_adjm, shellsize, numofpolies, damp, tolerance, possible_edgelist, poss_edge_num);
    
    endtime = clock();
    
    cpu_time_used = ((double) (endtime - starttime)) / CLOCKS_PER_SEC;

    //print file
    
    if(bestpolynum==-1)
    {
        //printf("No best poly found.");
        fprintf(outputPolys, "%d, NA, %d, %lf \n", shellsize, bestpolynum, cpu_time_used);
    }
    else
    {
    //printf("Best poly: \n %s\n\n",polylist[bestpolynum]->name);
    fprintf(outputPolys, "%d, %s, %d, %lf \n", shellsize, polylist[bestpolynum]->name, bestpolynum, cpu_time_used);
    }
    
    //free memory
    
    free(shellist);
    for (i=0; i<maxshellsize; i++)
    {
        free(shelldistmtx[i]);
    }
    free(shelldistmtx);
    
    for(i=0; i<(shellsize+1); i++)
    {
        free(adjacency_matrix[i]);
        free(temp_adjm[i]);
    }
    free(adjacency_matrix);
    free(temp_adjm);
    
    for(i=0; i<numofpolies; i++)
    {
        poly_destroy(polylist[i]);
    }
    free(polylist);
    
    
    
}

void ChemNetworkOrig::shellist_and_distmtx(int *shellist, double **shelldistmtx, double *atmS1, double *atmT1, int nd1, int nd2, int nsolvent1, int nsolute1, int nAtomS1, int s1t1cutoffnum, int *s1t1a, int *s1t1b, double *s1t1cutoff, int pbc, double xside, double yside, double zside, int *shellsize, int maxshellsize, int whichatmT1, int whichatmdexT1)
{
    int i, j, solute_node, solvent_node, crt;
    int shellist_indx, solvent_node_a, solvent_node_b;
    int solvent_atom_a, solvent_atom_b;
    int crit_atom_holder_a[1];
    int crit_atom_holder_b[1];
    double dist;
    
    int *eachnode_critatm = (int*)malloc(maxshellsize*sizeof(int));
    
    
    //create variables for pbc boxes.
    double *atmS1x, *atmS1y, *atmS1z, *atmS1xy, *atmS1yz, *atmS1zx;
    double *atmS1xminy, *atmS1minyz, *atmS1zminx, *atmS1xyz;
    double *atmS1xyminz, *atmS1minxyz, *atmS1xminyz;
    
    double *atmS1minx, *atmS1miny, *atmS1minz, *atmS1minxminy, *atmS1minyminz, *atmS1minzminx;
    double *atmS1minxy, *atmS1yminz, *atmS1minzx, *atmS1minxminyminz;
    double *atmS1minxminyz, *atmS1xminyminz, *atmS1minxyminz;
    
    //step 1: find list of atoms in solvation shell. Look at only 1 atom in solute.
    
    solute_node=nd2; //ensure solute node has correct numbering
    
    solvent_node=nd1; //start node numbering at 1 though array entries start at 0
    
    shellist_indx=0;
    
    for(i=0; i<nAtomS1; i = i + nsolvent1) //loop through solvent molecules by jumping the number of atoms per molecule
    {
        
        for(crt=0; crt<s1t1cutoffnum; crt++) //loop through diff criteria for s1-t1 interaction
        {
            if(whichatmT1+1==s1t1b[crt]) //ensure this criteria looks at the correct solute atom (in case of polyatomic solute)
            {
            dist=distanceMix(atmT1, atmS1, whichatmdexT1, i, s1t1b, s1t1a, crt); //look at the distance for this criteria
            
            if(shellist_indx<=maxshellsize) //ensure shell list not full already
            {
                if(dist<s1t1cutoff[crt]) //if the solvent molecule is within the cutoff
                {
                    //printf("solvent molecule %d is in the solvation shell.\n", solvent_node);
                    shellist[shellist_indx]=solvent_node; //get node name
                    eachnode_critatm[shellist_indx]=s1t1a[crt]; //get criteria atm for this node. 

                    shellist_indx=shellist_indx+1;
                    
                    
                    if (shellist_indx>maxshellsize) {
                        printf("WARNING: More atoms in solvent shell than given limit. May indicate CUTOFF too large. Vertex search ceased.\n");
                        printf("(Error in Solute molecule number %d. Numbering starts at 1.)\n", solute_node);
                        
                    }
                }
            }
            }
            
        }
        
        solvent_node=solvent_node+1;    //node counting increases by molecule. 
    }
    
    
    
    //create PBC if asked for.
    
    if(pbc !=0)
    {
        // allocate memory
        
        atmS1x = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1y = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1z = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1xy = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1yz = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1zx = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1xminy = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1minyz = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1zminx = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1xyz = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1xyminz = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1minxyz = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1xminyz = (double*)malloc((3*nAtomS1)*sizeof(double));
        
        atmS1minx = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1miny = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1minz = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1minxminy = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1minyminz = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1minzminx = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1minxy = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1yminz = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1minzx = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1minxminyminz = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1minxminyz = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1xminyminz = (double*)malloc((3*nAtomS1)*sizeof(double));
        atmS1minxyminz = (double*)malloc((3*nAtomS1)*sizeof(double));
        
        // get all the PBC Boxes
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1x[3*(i+1)-3] = atmS1[3*(i+1)-3] + xside;
            atmS1x[3*(i+1)-2] = atmS1[3*(i+1)-2];
            atmS1x[3*(i+1)-1] = atmS1[3*(i+1)-1];
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1minx[3*(i+1)-3] = atmS1[3*(i+1)-3] - xside;
            atmS1minx[3*(i+1)-2] = atmS1[3*(i+1)-2];
            atmS1minx[3*(i+1)-1] = atmS1[3*(i+1)-1];
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1y[3*(i+1)-3] = atmS1[3*(i+1)-3];
            atmS1y[3*(i+1)-2] = atmS1[3*(i+1)-2] + yside;
            atmS1y[3*(i+1)-1] = atmS1[3*(i+1)-1];
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1miny[3*(i+1)-3] = atmS1[3*(i+1)-3];
            atmS1miny[3*(i+1)-2] = atmS1[3*(i+1)-2] - yside;
            atmS1miny[3*(i+1)-1] = atmS1[3*(i+1)-1];
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1z[3*(i+1)-3] = atmS1[3*(i+1)-3];
            atmS1z[3*(i+1)-2] = atmS1[3*(i+1)-2];
            atmS1z[3*(i+1)-1] = atmS1[3*(i+1)-1] + zside;
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1minz[3*(i+1)-3] = atmS1[3*(i+1)-3];
            atmS1minz[3*(i+1)-2] = atmS1[3*(i+1)-2];
            atmS1minz[3*(i+1)-1] = atmS1[3*(i+1)-1] - zside;
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1xy[3*(i+1)-3] = atmS1x[3*(i+1)-3];
            atmS1xy[3*(i+1)-2] = atmS1x[3*(i+1)-2] + yside;
            atmS1xy[3*(i+1)-1] = atmS1x[3*(i+1)-1];
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1minxminy[3*(i+1)-3] = atmS1minx[3*(i+1)-3];
            atmS1minxminy[3*(i+1)-2] = atmS1minx[3*(i+1)-2] - yside;
            atmS1minxminy[3*(i+1)-1] = atmS1minx[3*(i+1)-1];
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1yz[3*(i+1)-3] = atmS1y[3*(i+1)-3];
            atmS1yz[3*(i+1)-2] = atmS1y[3*(i+1)-2];
            atmS1yz[3*(i+1)-1] = atmS1y[3*(i+1)-1] + zside;
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1minyminz[3*(i+1)-3] = atmS1miny[3*(i+1)-3];
            atmS1minyminz[3*(i+1)-2] = atmS1miny[3*(i+1)-2];
            atmS1minyminz[3*(i+1)-1] = atmS1miny[3*(i+1)-1] - zside;
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1zx[3*(i+1)-3] = atmS1z[3*(i+1)-3] + xside;
            atmS1zx[3*(i+1)-2] = atmS1z[3*(i+1)-2];
            atmS1zx[3*(i+1)-1] = atmS1z[3*(i+1)-1];
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1minzminx[3*(i+1)-3] = atmS1minz[3*(i+1)-3] - xside;
            atmS1minzminx[3*(i+1)-2] = atmS1minz[3*(i+1)-2];
            atmS1minzminx[3*(i+1)-1] = atmS1minz[3*(i+1)-1];
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1xminy[3*(i+1)-3] = atmS1x[3*(i+1)-3];
            atmS1xminy[3*(i+1)-2] = atmS1x[3*(i+1)-2] - yside;
            atmS1xminy[3*(i+1)-1] = atmS1x[3*(i+1)-1];
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1minxy[3*(i+1)-3] = atmS1minx[3*(i+1)-3];
            atmS1minxy[3*(i+1)-2] = atmS1minx[3*(i+1)-2] + yside;
            atmS1minxy[3*(i+1)-1] = atmS1minx[3*(i+1)-1];
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1minyz[3*(i+1)-3] = atmS1z[3*(i+1)-3];
            atmS1minyz[3*(i+1)-2] = atmS1z[3*(i+1)-2] - yside;
            atmS1minyz[3*(i+1)-1] = atmS1z[3*(i+1)-1];
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1yminz[3*(i+1)-3] = atmS1y[3*(i+1)-3];
            atmS1yminz[3*(i+1)-2] = atmS1y[3*(i+1)-2];
            atmS1yminz[3*(i+1)-1] = atmS1y[3*(i+1)-1] - zside;
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1zminx[3*(i+1)-3] = atmS1z[3*(i+1)-3] - xside;
            atmS1zminx[3*(i+1)-2] = atmS1z[3*(i+1)-2];
            atmS1zminx[3*(i+1)-1] = atmS1z[3*(i+1)-1];
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1minzx[3*(i+1)-3] = atmS1minz[3*(i+1)-3] + xside;
            atmS1minzx[3*(i+1)-2] = atmS1minz[3*(i+1)-2];
            atmS1minzx[3*(i+1)-1] = atmS1minz[3*(i+1)-1];
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1xyz[3*(i+1)-3] = atmS1xy[3*(i+1)-3];
            atmS1xyz[3*(i+1)-2] = atmS1xy[3*(i+1)-2];
            atmS1xyz[3*(i+1)-1] = atmS1xy[3*(i+1)-1] + zside;
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1minxminyminz[3*(i+1)-3] = atmS1minxminy[3*(i+1)-3];
            atmS1minxminyminz[3*(i+1)-2] = atmS1minxminy[3*(i+1)-2];
            atmS1minxminyminz[3*(i+1)-1] = atmS1minxminy[3*(i+1)-1] - zside;
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1xyminz[3*(i+1)-3] = atmS1xy[3*(i+1)-3];
            atmS1xyminz[3*(i+1)-2] = atmS1xy[3*(i+1)-2];
            atmS1xyminz[3*(i+1)-1] = atmS1xy[3*(i+1)-1] - zside;
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1minxminyz[3*(i+1)-3] = atmS1minxminy[3*(i+1)-3];
            atmS1minxminyz[3*(i+1)-2] = atmS1minxminy[3*(i+1)-2];
            atmS1minxminyz[3*(i+1)-1] = atmS1minxminy[3*(i+1)-1] + zside;
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1minxyz[3*(i+1)-3] = atmS1yz[3*(i+1)-3] - xside;
            atmS1minxyz[3*(i+1)-2] = atmS1yz[3*(i+1)-2];
            atmS1minxyz[3*(i+1)-1] = atmS1yz[3*(i+1)-1];
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1xminyminz[3*(i+1)-3] = atmS1minyminz[3*(i+1)-3] + xside;
            atmS1xminyminz[3*(i+1)-2] = atmS1minyminz[3*(i+1)-2];
            atmS1xminyminz[3*(i+1)-1] = atmS1minyminz[3*(i+1)-1];
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1xminyz[3*(i+1)-3] = atmS1minyz[3*(i+1)-3] + xside;
            atmS1xminyz[3*(i+1)-2] = atmS1minyz[3*(i+1)-2];
            atmS1xminyz[3*(i+1)-1] = atmS1minyz[3*(i+1)-1];
        }
        
        for(i=0; i < nAtomS1; i++)
        {
            atmS1minxyminz[3*(i+1)-3] = atmS1minxy[3*(i+1)-3];
            atmS1minxyminz[3*(i+1)-2] = atmS1minxy[3*(i+1)-2];
            atmS1minxyminz[3*(i+1)-1] = atmS1minxy[3*(i+1)-1] - zside;
        }
        
        
        //Go through and look for attachments to the solute in all the pbc boxes.
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1x, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1y, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1z, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1xy, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1yz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1zx, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1xminy, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1minyz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1zminx, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1xyz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1xyminz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1minxyz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1xminyz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1minx, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1miny, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1minz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1minxminy, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1minyminz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1minzminx, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1minxy, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1yminz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1minzx, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1minxminyminz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1minxminyz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1xminyminz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
        
        search_pbc_solute(&shellist_indx, shellist, atmT1, atmS1minxyminz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, eachnode_critatm, whichatmT1, whichatmdexT1);
    }
    
    
    //FOR TROUBLESHOOTING: lists out all molecules in shell (Recall 'node' numbering starts at 1)
    
    /*
     for(i=0; i<15; i=i+1)
     {
     printf("%d, ",shellist[i]);
     }
     */
    
    
    for(i=0;i<15;i=i+1)
    {
        if(shellist[i]==0)
        {
            *shellsize=i;
            break;
        }
    }
    
    //printf("%d\n",*shellsize);
    
    //Next: creating matrix of distances
    
    for(i=0; i<*shellsize; i=i+1)
    {
        solvent_node_a = shellist[i];
        solvent_atom_a = (solvent_node_a - 1)*nsolvent1;
        crit_atom_holder_a[0] = eachnode_critatm[i];
        
        for(j=i; j< *shellsize; j=j+1)
        {
            
            solvent_node_b=shellist[j];
            solvent_atom_b = (solvent_node_b-1)*nsolvent1;
            crit_atom_holder_b[0] = eachnode_critatm[j];
            
            if(pbc!=0)
            {
                dist=solvent_dist(atmS1, atmS1x, atmS1y, atmS1z, atmS1xy, atmS1yz, atmS1zx,
                                  atmS1xminy, atmS1minyz, atmS1zminx, atmS1xyz, atmS1xyminz,
                                  atmS1minxyz, atmS1xminyz, atmS1minx, atmS1miny, atmS1minz,
                                  atmS1minxminy, atmS1minyminz, atmS1minzminx, atmS1minxy,
                                  atmS1yminz, atmS1minzx, atmS1minxminyminz, atmS1minxminyz,
                                  atmS1xminyminz, atmS1minxyminz, *shellsize, shellist, pbc, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b);
            }
            else
            {
                
                dist = distance(atmS1, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b, 0);
                
            }
            shelldistmtx[i][j] = dist;
            shelldistmtx[j][i] = dist;
            
        }
    }
    
    
    //to print the matrix (troubleshooting)
    /*
    for(i=0; i<*shellsize;i=i+1)
    {
        for(j=0; j< *shellsize; j=j+1)
        {
            printf("%lf \t", shelldistmtx[i][j]);
        }
        printf("\n");
    }
    */
    
    
    
    
    
    
    
    
    //free memory
    if(pbc!=0)
    {
        free(atmS1x);
        free(atmS1y);
        free(atmS1z);
        free(atmS1xy);
        free(atmS1yz);
        free(atmS1zx);
        free(atmS1xminy);
        free(atmS1minyz);
        free(atmS1zminx);
        free(atmS1xyz);
        free(atmS1xyminz);
        free(atmS1minxyz);
        free(atmS1xminyz);
        free(atmS1minx);
        free(atmS1miny);
        free(atmS1minz);
        free(atmS1minxminy);
        free(atmS1minyminz);
        free(atmS1minzminx);
        free(atmS1minxy);
        free(atmS1yminz);
        free(atmS1minzx);
        free(atmS1minxminyminz);
        free(atmS1minxminyz);
        free(atmS1xminyminz);
        free(atmS1minxyminz);
    }
    
}

void ChemNetworkOrig::search_pbc_solute(int *shellist_indx, int *shellist, double *atmT1, double *atmS1box, int nd1, int nd2, int nsolvent1, int nsolute1, int nAtomS1, int s1t1cutoffnum, int *s1t1a, int *s1t1b, double *s1t1cutoff, int *eachnode_critatm, int whichatmT1, int whichatmdexT1)
{
    
    
    int i, solute_node, solvent_node, crt;
    double dist;
    
    solute_node=nd2; //ensure nodes have correct numbering
    
    solvent_node=nd1;
    
    for(i=0; i<nAtomS1; i = i + nsolvent1) //loop through solvent molecules by jumping the number of atoms per molecule
    {
        
        for(crt=0; crt<s1t1cutoffnum; crt++) //loop through diff criteria for s1-t1 interaction
        {
            //look at the distance for this criteria
            if(whichatmT1+1==s1t1b[crt]) //ensure this criteria looks at this atom
            {

            dist=distanceMix(atmT1, atmS1box, whichatmdexT1, i, s1t1b, s1t1a, crt);
            
            if(*shellist_indx<=14) //ensure shell list not full already
            {
                if(dist<s1t1cutoff[crt]) //if it is within the cutoff
                {
                    // printf("solvent molecule %d is in the solvation shell.\n", solvent_node);
                    shellist[*shellist_indx]=solvent_node;
                    eachnode_critatm[*shellist_indx]=s1t1a[crt];
                    *shellist_indx=*shellist_indx+1;
                    if (*shellist_indx>14) {
                        printf("WARNING: Too many atoms in solvent shell. May indicate CUTOFF too large. May yield inaccurate results.\n");
                        
                    }
                }
            }
            }
            
        }
        
        solvent_node=solvent_node+1;
    }
    
}

double ChemNetworkOrig::solvent_dist(double *atmS1, double *atmS1x, double *atmS1y, double *atmS1z, double *atmS1xy, double *atmS1yz, double *atmS1zx, double *atmS1xminy, double *atmS1minyz, double *atmS1zminx, double *atmS1xyz, double *atmS1xyminz,
                    double *atmS1minxyz, double *atmS1xminyz, double *atmS1minx, double *atmS1miny, double *atmS1minz,
                    double *atmS1minxminy, double *atmS1minyminz, double *atmS1minzminx, double *atmS1minxy,
                    double *atmS1yminz, double *atmS1minzx, double *atmS1minxminyminz, double *atmS1minxminyz,
                    double *atmS1xminyminz, double *atmS1minxyminz, int shellsize, int *shellist, int pbc, int solvent_atom_a, int solvent_atom_b, int *crit_atom_holder_a, int *crit_atom_holder_b)
{
    
    double truedist;
    
    //go through all the distances and choose the smallest one.
    
    
    truedist = distance(atmS1, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b, 0);
    
    truedist = fmin(truedist, distanceMix(atmS1, atmS1x, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1y, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1z, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1xy, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1yz, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1zx, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1xminy, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1minyz, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1zminx, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1xyz, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1xyminz, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1minxyz, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1xminyz, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1minx, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1miny, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1minz, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1minxminy, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1minyminz, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1minzminx, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1minxy, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1yminz, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1minzx, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1minxminyminz, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1minxminyz, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1xminyminz, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    truedist = fmin(truedist, distanceMix(atmS1, atmS1minxyminz, solvent_atom_a, solvent_atom_b, crit_atom_holder_a, crit_atom_holder_b,   0));
    
    
    return truedist;
    
}

//pagerank function by Tiecheng Zhou
//Needs edit to avoid memory leak!!!
double * ChemNetworkOrig::pagerank(int ** adjM, int nodes, double damp, double tolerance) // regarded as a directed graph
{
	double * a = (double *)malloc(nodes * sizeof(double));
	double * b = (double *)malloc(nodes * sizeof(double));
	int i,j,k,L;
	double diff;
	double sum,factor;
	
    //set initial value
    for(i=0;i<nodes;i++)
    {
        a[i]=1.00/nodes; //if use 1/N then regard as integral the result is 0
        b[i]=1.00/nodes;
    }
    
	diff=1.0;
	while(diff>=tolerance)
    {
        for(i=0;i<nodes;i++)
        {
            b[i]=(1.00-damp)/nodes;
            for(j=0;j<nodes;j++)
            {
                if(adjM[j][i]==1) //j --> i
                {
                    L=0;
                    for(k=0;k<nodes;k++)
                    {
                        if(adjM[j][k]==1)
                            L+=1; //outdegree of node j
                    }
                    b[i]=b[i]+damp*a[j]/L;
                }
            }
        }
        //normalize pagerank to make sure sum=1;
		sum=0.0;
		for(i=0;i<nodes;i++)
		{
			sum+=b[i];
		}
		factor=1/sum;
		for(i=0;i<nodes;i++)
		{
			b[i]*=factor;
		}
        
        diff=0.00;
        for(i=0;i<nodes;i++)
        {
            diff+=fabs(a[i]-b[i]); // absolute value of difference
        }
        for(i=0;i<nodes;i++) // copy the new pagerank information into pgrk[]
        {
            a[i]=b[i];
        }
    }
	return(a);
}

void ChemNetworkOrig::basearray(int n, int length, int *narray, int base)
{
    /*This makes "narray" into an array of values determined by the binary
     (or ternary, or ...) value of n. For base of 2, array[0] is the 2^0
     place, array[1] is the 2^1 place, and so on. If n is too large for the
     given length, any places farther out than 2^length are truncated. That
     is, for base 2 this function will give repeat values modulo 2^length.
     This function is useful any time you wish to get all combinations of values
     of on/off (or possibly a ternary state like red/yellow/green) for a
     certain number of objects. */
    
    int i=0;
    
    if(n<=0) //will return 0 for negative numbers and 0.
    {
        for(i=0; i<length; i++)
        {
            narray[i]=0;
        }
    }
    
    while(n>0){
        narray[i] = n % base;
        n=n/base;
        i=i+1;
        
    }
    
    return;
    
}

void ChemNetworkOrig::makeallpolys(poly **polylist)
{
    
    int sqedgelist[8]={0,1,1,2,2,3,3,0};
    polylist[0]=poly_create((char *)"Square", sqedgelist, 4, 4);
    
    //output adjacency matrix
    /*
     for(i=0; i<sq->nedge+1; i++)
     {
     for(j=0; j<sq->nedge+1; j++)
     {
     printf("%d",sq->adjmat[i][j]);
     }
     printf("\n");
     }
     //output pr value
     printf("%lf \n",sq->char_pr);
     */
    
    int tetedgelist[12] = {0,1, 1,2, 2,0, 0,3, 1,3, 2,3};
    polylist[1]=poly_create((char *)"Tetrahedron", tetedgelist, 4, 6);
    
    int sqpyEL[16] = {0,1, 1,2, 2,3, 3,0, 0,4, 1,4, 2,4, 3,4};
    polylist[2]=poly_create((char *)"Square Pyramid",sqpyEL, 5, 8);
    
    int tbpEL[18] = {0,1, 1,2, 2,3, 3,0, 3,1, 3,4, 0,4, 1,4, 2,4};
    polylist[3]=poly_create((char *)"Triangular Bipyramid",tbpEL,5,9);
    
    int wedgeEL[14] = {0,1, 1,2, 2,3, 3,0, 0,4, 1,4, 3,4};
    polylist[4]=poly_create((char *)"Wedge",wedgeEL,5,7);
    
    int octEL[24] = {0,1, 1,2, 2,3, 3,0, 0,4, 1,4, 2,4, 3,4, 0,5, 1,5, 2,5, 3,5};
    polylist[5]=poly_create((char *)"Octahedron", octEL, 6, 12);
    
    int tprismEL[18] = {0,1, 1,2, 2,0, 3,4, 4,5, 5,3, 1,3, 2,4, 0,5};
    polylist[6] = poly_create((char *)"Trigonal Prism",tprismEL, 6, 9);
    
    int pentpyEL[20] = {0,1, 1,2, 2,3, 3,4, 4,0, 0,5, 1,5, 2,5, 3,5, 4,5};
    polylist[7] = poly_create((char *)"Pentagonal Pyramid", pentpyEL, 6,10);
    
    int csqpyEL[22] = {0,1, 1,2, 2,3, 3,0, 0,4, 1,4, 2,4, 3,4, 0,5, 1,5, 4,5};
    polylist[8] = poly_create((char *)"Capped Square Pyramid", csqpyEL, 6, 11);
    
    int sdppEL[16] = {0,1, 1,2, 2,3, 3,4, 4,0, 0,5, 1,5, 2,5};
    polylist[9] = poly_create((char *)"Side-Deleted Pentagonal Pyramid", sdppEL, 6, 8);
    
    int flapoctEL[20] = {0,1, 1,2, 2,3, 3,0, 0,4, 1,4, 2,4, 3,4, 1,5, 2,5};
    polylist[10] = poly_create((char *)"Flap Octahedron", flapoctEL, 6, 10);
    
    int eltripEL[24] = {0,1, 1,2, 2,0, 0,3, 1,3, 2,3, 4,0, 5,1, 6,2, 4,5, 5,6, 6,4};
    polylist[11] = poly_create((char *)"Elongated Triangular Pyramid", eltripEL, 7, 12);
    
    int pentbpEL[30] = {0,1, 1,2, 2,3, 3,4, 4,0, 0,5, 1,5, 2,5, 3,5, 4,5, 0,6, 1,6, 2,6, 3,6, 4,6};
    polylist[12] = poly_create((char *)"Pentagonal Bipyramid", pentbpEL, 7, 15);
    
    int augtripEL[26] = {0,1, 1,2, 2,0, 3,4, 4,5, 5,3, 1,3, 2,4, 0,5, 0,6, 1,6, 3,6, 5,6};
    polylist[13] = poly_create((char *)"Augmented Trigonal Prism", augtripEL, 7, 13);
    
    int capoctEL[30] = {0,1, 1,2, 2,3, 3,0, 0,4, 1,4, 2,4, 3,4, 0,5, 1,5, 4,5, 0,6, 1,6, 2,6, 3,6};
    polylist[14] = poly_create((char *)"Cap Octahedron", capoctEL, 7, 15);
    
    int ccubeEL[18] = {0,1, 1,2, 2,3, 3,0, 0,4, 3,6, 2,5, 4,6, 5,6};
    polylist[15] = poly_create((char *)"Cut Cube", ccubeEL, 7, 9);
    
    int csdisphEL[26] = {0,1, 1,2, 2,3, 3,4, 4,0, 0,5, 1,5, 2,5, 3,5, 4,5, 0,6, 1,6, 2,6};
    polylist[16] = poly_create((char *)"Cut Snub Disphenoid", csdisphEL, 7, 13);
    
    int cubeEL[24] = {0,1, 1,2, 2,3, 3,0, 0,4, 1,5, 2,6, 3,7, 4,5, 5,6, 6,7, 7,4};
    polylist[17] = poly_create((char *)"Cube", cubeEL, 8, 12);
    
    int sqantiEL[32] = {0,1, 1,2, 2,3, 3,0, 4,5, 5,6, 6,7, 7,4, 0,7, 0,4, 1,4, 1,5, 2,5, 2,6, 3,6, 3,7};
    polylist[18] = poly_create((char *)"Square Antiprism", sqantiEL, 8, 16);
    
    int batpEL[34] = {0,1, 1,2, 2,0, 3,4, 4,5, 5,3, 0,3, 1,4, 2,5, 0,6, 1,6, 3,6, 4,6, 1,7, 2,7, 4,7, 5,7};
    polylist[19] = poly_create((char *)"Bi-augmented Trigonal Prism", batpEL, 8, 17);
    
    int sndisphEL[36] = {0,1, 1,2, 2,3, 3,4, 4,0, 3,5, 5,0, 4,5, 0,6, 1,6, 2,6, 3,6, 4,6, 0,7, 1,7, 2,7, 3,7, 5,7};
    polylist[20] = poly_create((char *)"Snub Disphenoid", sndisphEL, 8, 18);
    
    int gyrobiEL[28] = {0,1, 1,2, 2,0, 2,3, 3,4, 4,0, 1,6, 2,7, 3,7, 4,5, 0,5, 5,6, 6,7, 7,5};
    polylist[21] = poly_create((char *)"Gyrobifastigium", gyrobiEL, 8, 14);
    
    int eltribpEL[30] = {0,1, 1,2, 2,0, 3,4, 4,5, 5,3, 0,3, 1,4, 2,5, 0,6, 1,6, 2,6, 3,7, 4,7, 5,7};
    polylist[22] = poly_create((char *)"Elongated Triangular Bipyramid", eltribpEL, 8, 15);
    
    int triatpEL[42] = {0,1, 1,2, 2,0, 3,4, 4,5, 5,3, 0,3, 1,4, 2,5, 0,6, 1,6, 3,6, 4,6, 1,7, 2,7, 4,7, 5,7, 0,8, 2,8, 3,8, 5,8};
    polylist[23] = poly_create((char *)"Tri-augmented Trigonal Prism", triatpEL, 9, 21);
    
    int mccubeEL[32] = {0,1, 1,2, 2,3, 3,0, 0,4, 1,5, 2,6, 3,7, 4,5, 5,6, 6,7, 7,4, 0,8, 1,8, 2,8, 3,8};
    polylist[24] = poly_create((char *)"Monocapped Cube", mccubeEL, 9, 16);
    
    int mcsqantiEL[40] = {0,1, 1,2, 2,3, 3,0, 4,5, 5,6, 6,7, 7,4, 0,7, 0,4, 1,4, 1,5, 2,5, 2,6, 3,6, 3,7, 0,8, 1,8, 2,8, 3,8};
    polylist[25] = poly_create((char *)"Monocapped Square Antiprism", mcsqantiEL, 9, 20);
    
    int tdimicosEL[30] = {0,1, 1,2, 2,3, 3,4, 4,5, 5,0, 1,3, 3,5, 5,1, 0,6, 2,7, 4,8, 6,7, 7,8, 8,6};
    polylist[26] = poly_create((char *)"Tridiminished Icosahedron", tdimicosEL, 9, 15);
    
    int bcsqantiEL[48] = {0,1, 1,2, 2,3, 3,0, 4,5, 5,6, 6,7, 7,4, 0,7, 0,4, 1,4, 1,5, 2,5, 2,6, 3,6, 3,7, 0,8, 1,8, 2,8, 3,8, 4,9, 5,9, 6,9, 7,9};
    polylist[27] = poly_create((char *)"Bicapped Square Antiprism", bcsqantiEL, 10, 24);
    
    int bccube1EL[40] = {0,1, 1,2, 2,3, 3,0, 0,4, 1,5, 2,6, 3,7, 4,5, 5,6, 6,7, 7,4, 0,8, 1,8, 2,8, 3,8, 4,9, 5,9, 6,9, 7,9};
    polylist[28] = poly_create((char *)"Bicapped Cube 1", bccube1EL, 10, 20);
    //caps are opposite to each other
    
    int bccube2EL[40] = {0,1, 1,2, 2,3, 3,0, 0,4, 1,5, 2,6, 3,7, 4,5, 5,6, 6,7, 7,4, 0,8, 1,8, 2,8, 3,8, 0,9, 1,9, 5,9, 4,9};
    polylist[29] = poly_create((char *)"Bicapped Cube 2", bccube2EL, 10, 20);
    //caps are next to each other
    
    int bidimicosEL[40] = {0,1, 1,2, 2,3, 3,4, 4,5, 5,0, 1,3, 3,5, 5,1, 0,6, 2,7, 4,8, 6,7, 7,8, 8,6, 0,9, 1,9, 2,9, 6,9, 7,9};
    polylist[30] = poly_create((char *)"Bidiminished Icosahedron", bidimicosEL, 10, 20);
    
    int tetatp1EL[48] = {0,1, 1,2, 2,0, 3,4, 4,5, 5,3, 0,3, 1,4, 2,5, 0,6, 1,6, 3,6, 4,6, 1,7, 2,7, 4,7, 5,7, 0,8, 2,8, 3,8, 5,8, 0,9, 1,9, 2,9};
    polylist[31] = poly_create((char *)"Tetra-augmented Trigonal Prism", tetatp1EL, 10, 24);
    
    int sphenoEL[44] = {0,1, 1,2, 2,3, 3,4, 4,0, 0,5, 1,5, 2,5, 3,5, 4,5, 1,6, 2,6, 2,8, 3,8, 3,9, 4,9, 4,7, 0,7, 6,8, 7,9, 6,7, 8,9};
    polylist[32] = poly_create((char *)"Sphenocorona", sphenoEL, 10, 22);
    
    int dctprEL[12] = {0,1, 1,2, 2,0, 3,4, 4,5, 5,3};
    polylist[33] = poly_create((char *)"Doubly Crossed Trigonal Prism", dctprEL, 6, 6);
    
    //polyhedra created following testing:
    
    int cfloc1edgelist[26] = {0,1, 0,2, 0,3, 0,4, 0,5, 1,2, 1,4, 1,6, 2,3, 2,6, 3,4, 3,5, 4,5};
	polylist[34] = poly_create((char *)"Cap Flap Octahedron V1", cfloc1edgelist, 7, 13);
    
    int cfloc2edgelist[26] = {0,1, 0,2, 0,3, 0,4, 0,5, 1,2, 1,4, 1,6, 2,3, 2,6, 3,4, 2,5, 4,5};
	polylist[35] = poly_create((char *)"Cap Flap Octahedron V2", cfloc2edgelist, 7, 13);
    
	int cpentbipyel[36] = {0,1, 1,2, 2,3, 3,4, 4,0, 0,5, 1,5, 2,5, 3,5, 4,5, 0,6, 1,6, 2,6, 3,6, 4,6, 0,7, 1,7, 5,7};
    polylist[36]=poly_create((char *)"Capped Pentagonal Bipyramid", cpentbipyel, 8, 18);
    
    int bisdpentbipy1el[22] = {0,1, 1,2, 2,3, 3,4, 4,0, 0,5, 1,5, 2,5, 0,6, 1,6, 2,6};
    polylist[37] = poly_create((char *)"Double Side Deleted Pentagonal Bipyramid V1", bisdpentbipy1el, 7, 11);
    
    int bisdpentbipy2el[22] = {0,1, 1,2, 2,3, 3,4, 4,0, 0,5, 1,5, 2,5, 1,6, 2,6, 3,6};
    polylist[38] = poly_create((char *)"Double Side Deleted Pentagonal Bipyramid V2", bisdpentbipy2el, 7, 11);
    
    int bisdpentbipy3el[22] = {0,1, 1,2, 2,3, 3,4, 4,0, 0,5, 1,5, 2,5, 2,6, 3,6, 4,6};
    polylist[39] = poly_create((char *)"Double Side Deleted Pentagonal Bipyramid V3", bisdpentbipy3el, 7, 11);
    
    int hexpyel[24] = {0,1, 1,2, 2,3, 3,4, 4,5, 5,0, 0,6, 1,6, 2,6, 3,6, 4,6, 5,6};
    polylist[40] = poly_create((char *)"Hexagonal Pyramid", hexpyel, 7, 12);
    
    int sidedelhexpy[20] = {0,1, 1,2, 2,3, 3,4, 4,5, 5,0, 0,6, 1,6, 2,6, 3,6};
    polylist[41] = poly_create((char *)"Side-Deleted Hexagonal Pyramid",sidedelhexpy, 7, 10);
    
    int hexbipyel[36] = {0,1, 1,2, 2,3, 3,4, 4,5, 5,0, 0,6, 1,6, 2,6, 3,6, 4,6, 5,6, 0,7, 1,7, 2,7, 3,7, 4,7, 5,7};
    polylist[42] = poly_create((char *)"Hexagonal Bipyramid", hexbipyel, 8, 18);
    
    int sdhexbipyel[32] = {0,1, 1,2, 2,3, 3,4, 4,5, 5,0, 0,6, 1,6, 2,6, 3,6, 4,6, 5,6, 0,7, 1,7, 2,7, 3,7};
    polylist[43] = poly_create((char *)"Side-deleted Hexagonal Bipyramid", sdhexbipyel, 8, 16);
    
    int caphexbipyel[42] = {0,1, 1,2, 2,3, 3,4, 4,5, 5,0, 0,6, 1,6, 2,6, 3,6, 4,6, 5,6, 0,7, 1,7, 2,7, 3,7, 4,7, 5,7, 0,8, 1,8, 6,8};
    polylist[44] = poly_create((char *)"Capped Hexagonal Bipyramid", caphexbipyel, 9, 19);
    
    int cpentpyel[26] = {0,1, 1,2, 2,3, 3,4, 4,0, 0,5, 1,5, 2,5, 3,5, 4,5, 0,6, 1,6, 5,6};
    polylist[45] = poly_create((char *)"Capped Pentagonal Pyramid", cpentpyel, 7, 13);
    
    int bicpentpy1el[32] = {0,1, 1,2, 2,3, 3,4, 4,0, 0,5, 1,5, 2,5, 3,5, 4,5, 0,6, 1,6, 5,6, 1,7, 2,7, 5,7};
    polylist[46] = poly_create((char *)"Bicapped Pentagonal Pyramid V1", bicpentpy1el, 8, 16);
    
    int bicpentpy2el[32] = {0,1, 1,2, 2,3, 3,4, 4,0, 0,5, 1,5, 2,5, 3,5, 4,5, 0,6, 1,6, 5,6, 2,7, 3,7, 5,7};
    polylist[47] = poly_create((char *)"Bicapped Pentagonal Pyramid V2", bicpentpy2el, 8, 16);
    
    int openedgecubeel[22] = {0,1, 1,2, 2,3, 3,0, 4,5, 5,6, 6,7, 7,4, 0,4, 1,5, 2,6};
    polylist[48] = poly_create((char *)"Open Edged Cube", openedgecubeel, 8, 11);
    
    int rhombicprismel[28] = {0,1, 1,2, 2,3, 3,0, 4,5, 5,6, 6,7, 7,4, 0,4, 1,5, 2,6, 3,7, 0,2, 4,6};
    polylist[49] = poly_create((char *)"Rhombic Prism", rhombicprismel, 8, 14);
    
    int capatpel[32] = {0,1, 1,2, 2,0, 3,4, 4,5, 5,3, 0,3, 1,4, 2,5, 0,6, 1,6, 2,6, 0,7, 3,7, 1,7, 4,7};
    polylist[50] = poly_create((char *)"Capped Augmented Trigonal Prism", capatpel, 8, 16);
    
    
    return;
    
}

ChemNetworkOrig::poly * ChemNetworkOrig::poly_create(char *name, int *edgelist, int nvtx, int nedge)
{
    int i;
    poly *newpoly;
    double damp=.85;
    double tolerance=10E-7;
    
    newpoly = (poly*)malloc(1*sizeof(poly));
    newpoly->name = name;
    
    newpoly->edgelist = (int*)malloc(2*nedge*sizeof(int));
    for(i=0; i<2*nedge; i++)
    {
        newpoly->edgelist[i] = edgelist[i];
    }
    
    
    newpoly->nvtx = nvtx;
    newpoly->nedge = nedge;
    
    newpoly->adjmat = (int**)malloc((nvtx+1)*sizeof(int*));
    for(i=0; i<nvtx+1; i++)
    {
        newpoly->adjmat[i]=(int*)calloc((nvtx+1),sizeof(int));
    }
    
    
    for(i=0; i<nvtx; i++) //make "star" graph
    {
        newpoly->adjmat[i][nvtx]=1;
        newpoly->adjmat[nvtx][i]=1;
        
    }
    
    for(i=0; i<2*nedge; i+=2) //add other edges (undirected)
    {
        newpoly->adjmat[edgelist[i]][edgelist[i+1]]=1;
        newpoly->adjmat[edgelist[i+1]][edgelist[i]]=1;
    }
    newpoly->prvect=(double*)malloc((nvtx+1)*sizeof(double));
    
    newpoly->prvect = pagerank(newpoly->adjmat, nvtx+1, damp, tolerance);
    newpoly->char_pr = newpoly->prvect[nvtx];
    
    return newpoly;
    
    
}

void ChemNetworkOrig::poly_destroy(poly *polyhedron)
{
    int i;
    free(polyhedron->prvect);
    
    free(polyhedron->edgelist);
    
    for(i=0; i<polyhedron->nvtx+1; i++)
    {
        free(polyhedron->adjmat[i]);
    }
    free(polyhedron->adjmat);
    
    
    free(polyhedron);
    
    return;
}

int ChemNetworkOrig::polymatch(poly **polylist, int **adjacency_matrix, double **shelldistmtx, int **temp_adjm, int shellsize, int numofpolies, double damp, double tolerance, int **possible_edgelist, int poss_edge_num)
{
    int *whichedges;
    double avgdist, vardist, stddist;
    int badexclude;
    int isometrycheck;
    int vtx5;
    int *vtx4 = (int*)calloc(2,sizeof(int));
    int bestpolynum = -1;
    double bestpolystd;
    int i, j, k, l;
    double *adj_PR;
    int ncombos;
    
    
    
    
    
    whichedges = (int*)malloc(poss_edge_num*sizeof(int)); //an array of 1's or 0's since an edge either exists or not
    
    ncombos = pow(2,poss_edge_num); //there are 2^n possible combinations of n edges on/off
    
    for(i=0; i<ncombos; i++)
    {
        //initialize temp_adjm to have same values as adjm.
        
        for(j=0; j<shellsize+1; j++)
        {
            for(k=j; k<shellsize+1; k++)
            {
                temp_adjm[j][k]=adjacency_matrix[j][k];
                temp_adjm[k][j]=adjacency_matrix[k][j];
            }
        }
        
        //determine which edges in the edgelist will be added
        basearray(i,poss_edge_num,whichedges,2);
        
        //add those edges
        
        for(j=0; j<poss_edge_num; j++)
        {
            if(whichedges[j])
            {
                temp_adjm[possible_edgelist[j][0]][possible_edgelist[j][1]] = 1;
                temp_adjm[possible_edgelist[j][1]][possible_edgelist[j][0]] = 1;
            }
        }
        
        //output new adjacency matrix (for troubleshooting)
        /*
         for(j=0; j<shellsize+1; j++)
         {
         for(k=0; k<shellsize+1; k++)
         {
         printf("%d", temp_adjm[j][k]);
         }
         printf("\n");
         }
         printf("\n");
         */
        
        //find page rank vector of new matrix
        adj_PR = pagerank(temp_adjm, shellsize+1, damp, tolerance);
        
        //Output page rank values (for troubleshooting)
        /*
         for(j=0; j<shellsize+1; j++)
         {
         printf("%lf\t",adj_PR[j]);
         if((j+1)%5==0)
         {
         printf("\n");
         }
         }
         printf("\n\n");
         */
        
        //output PR last entry (for troubleshooting)
        /*
         printf("%d: %lf \t| ",i ,adj_PR[shellsize]);
         */
        
        //check for match with each polyhedron
        for(j=0; j<numofpolies; j++)
        {
            if(shellsize == polylist[j]->nvtx) //if it has the right number of vertices then keep checking
            {
                
                if(fabs(polylist[j]->char_pr - adj_PR[shellsize]) < tolerance)
                    //if pageranks are within tolerance of each other they match, since algorithm requires a tolerance-based convergence.
                {
                    
                    // printf("%s \t",polylist[j]->name);
                    
                    //calculate average edge length
                    avgdist = 0;
                    for (k = 0; k<shellsize; k++)
                    {
                        for(l=k; l<shellsize; l++)
                        {
                            avgdist = avgdist + shelldistmtx[k][l] * temp_adjm[l][k]; //add to average entry if temp_adjm != 0
                        }
                    }
                    
                    avgdist = avgdist / polylist[j]->nedge ; //divide by number of edges
                    
                    // printf("%lf \t", avgdist);
                    
                    //calculate variance of edge lenth
                    
                    vardist=0;
                    
                    for(k=0; k<shellsize; k++)
                    {
                        for(l=k+1; l<shellsize; l++)
                        {
                            if(temp_adjm[l][k]!=0)
                            {
                                vardist = vardist + pow((shelldistmtx[k][l]-avgdist),2)/polylist[j]->nedge;
                            }
                        }
                    }
                    // printf("%lf \t", vardist);
                    // printf("%lf \n", vardist/avgdist);
                    
                    stddist=sqrt(vardist);
                    
                    //were distances < mean excluded? What about distances < mean + .2stdev?
                    
                    badexclude = 0;
                    for(k=0; k<shellsize; k++)
                    {
                        for(l=k+1; l<shellsize; l++)
                        {
                            if(temp_adjm[l][k] == 0 && shelldistmtx[k][l]<avgdist+.01)
                            {
                                // printf("Excluded values < mean \n");
                                l=shellsize; //break out of 2 loops
                                k=shellsize;
                                badexclude = 1;
                            }
                            //if you want to check if it excluded values  < mean + 1 stdev
                            else
                            {
                                if(temp_adjm[l][k] ==0 && shelldistmtx[k][l]<(avgdist+.2*stddist+.01))
                                {
                                    // printf("Excluded values < mean+.2 stdev \n");
                                    l=shellsize;
                                    k=shellsize;
                                }
                                else
                                {
                                    if(temp_adjm[l][k]!=0 && shelldistmtx[k][l]>avgdist+stddist)
                                    {
                                        // printf("Included values > mean + stdev \n");
                                        l=shellsize;
                                        k=shellsize;
                                    }
                                }
                            }
                            
                        }
                    }
                    
                    
                    //for bidiminished icosahedron vs. bicapped cube v2
                    if(j==29 || j==30)
                    {
                        isometrycheck =0;
                        
                        for(k=0; k<shellsize; k++)
                        {
                            if(fabs(adj_PR[k]-polylist[29]->prvect[3]) <tolerance) //find 5-connected water.
                            {
                                vtx5=k;
                                break;
                            }
                        }
                        for(k=0; k<shellsize; k++)
                        {
                            //if connnected to 5- and 3-connected (has right pr)
                            if((fabs(adj_PR[k]-polylist[29]->prvect[7]) <tolerance) && temp_adjm[k][vtx5]==1)
                            {
                                if(vtx4[0]!=0)
                                {
                                    vtx4[0]=k;
                                }
                                else
                                {
                                    vtx4[1]=k;
                                    //check if they are connected and should't be, or aren't connected and should be.
                                    if((temp_adjm[vtx4[0]][vtx4[1]]==0 && j==30) || (temp_adjm[vtx4[0]][vtx4[1]]==1 && j==29))
                                    {
                                        badexclude = 1;
                                    }
                                }
                            }
                        }
                    }
                    
                    
                    
                    
                    
                    if(badexclude == 0)
                    {
                        if(bestpolynum == -1)
                        {
                            bestpolynum = j;
                            bestpolystd = stddist;
                        }
                        else
                        {
                            if(polylist[bestpolynum]->nedge < polylist[j]->nedge)
                            {
                                bestpolynum = j;
                                bestpolystd = stddist;
                            }
                            else
                            {
                                if(polylist[bestpolynum]->nedge == polylist[j]->nedge && stddist < bestpolystd)
                                {
                                    bestpolynum = j;
                                    bestpolystd = stddist;
                                    
                                }
                            }
                        }
                        
                        
                        
                        
                        /*
                         //troubleshoot: output name, PR vector etc of polyhedron... just because.
                         
                         printf("%s\n %lf \t %lf \t %d \n\n", polylist[j]->name, avgdist, stddist, polylist[j]->nedge);
                         
                         for(k=0; k<shellsize+1; k++)
                         {
                         printf("\%lf \t ", polylist[j]->prvect[k]);
                         }
                         */
                        
                    }
                    
                    
                }
            }
        }
        
    }
    
    free(vtx4);
    free(whichedges);
    
    return(bestpolynum);
    
    
    
    
}

void ChemNetworkOrig::pdfshell(double **shelldistmtx, int shellsize, FILE *outputshellpdf)
{
    int i, j;
    
    for(i=0; i<shellsize; i++)
    {
        for(j=i+1; j<shellsize; j++)
        {
            fprintf(outputshellpdf, "%d, %lf\n", shellsize, shelldistmtx[i][j]);
        }
    }
    
    return;
}

/*Doesn't seem to be working correctly

void ChemNetworkOrig::varyShellFunction(double varyMin, double varyMax, int varyBreaks, int maxshellsize, FILE *outvariedShellDists, double *atmS1, double *atmT1, int node1Start, int node2Start, int nsolvent1, int nsolute1, int nAtomS1, int *s1t1a, int *s1t1b, double *s1t1cutoff, int pbc, double xside, double yside, double zside, int whichatmT1, int whichatmdexT1, int s1t1cutoffnum)
{
    int i,j;
    int *shellsizes = calloc(varyBreaks+1, sizeof(int));
    int **shellists;
    double **shelldistmtx;
    double test_cutoff[1];
    
    //prepare variables for data entry

    shellists=(int**)malloc((varyBreaks+1)*sizeof(int*));
    for(i=0; i<(varyBreaks+1); i++)
    {
        shellists[i]=calloc(maxshellsize+1, sizeof(int));
    }
    shelldistmtx=(double**)malloc(maxshellsize*sizeof(double));
    for(i=0; i<maxshellsize; i++)
    {
        shelldistmtx[i] = (double*)malloc(maxshellsize*sizeof(double));
    }
    
    
    
    for(i=0; i<varyBreaks+1; i++)
    {

        //get cutoff distance
        test_cutoff[0] = varyMin + i*(varyMax-varyMin)/varyBreaks;
        fprintf(outvariedShellDists, "%lf, ", test_cutoff[0]);
        
        //fill the arrays with values
        shellist_and_distmtx(shellists[i], shelldistmtx, atmS1, atmT1, node1Start, node2Start, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, test_cutoff, pbc, xside, yside, zside, &shellsizes[i], maxshellsize, whichatmT1, whichatmdexT1);
        
        
    }
    
    //output in columnar format for easy read-in to R?
    //Might be more practical to just print j=0 to j=maxshellsize-1
    
    fprintf(outvariedShellDists, "\n");
    
    for(i=0; i<maxshellsize; i++)
    {
        for(j=0; j<varyBreaks+1; j++)
        {
            fprintf(outvariedShellDists, "%d, ", shellists[j][i]);
        }
        fprintf(outvariedShellDists, "\n");
    }
    
    
    for(i=0; i<varyBreaks+1; i++)
    {
        free(shellists[i]);
    }
    free(shellists);
    
    free(shellsizes);
    
}
*/

int * ChemNetworkOrig::findunique(int *list, int numinlist, int *numunique)
{
    int i, j, count;
    int matches;
    *numunique = 0;
    int *newlist;
    
    for(i=0; i<numinlist; i++)
    {
        matches=0;
        for(j=0; j<i; j++)
        {
            if (list[j]==list[i])
            {
                matches=matches+1;
            }
        }
        if(matches==0)
        {
            *numunique=*numunique+1;
        }
    }
    
    newlist=(int*)malloc((*numunique)*sizeof(int));
    
    count=0;
    for(i=0; i<numinlist; i++)
    {
        matches=0;
        for(j=0; j<i; j++)
        {
            if(list[j]==list[i])
            {
                matches=matches+1;
            }
        }
        if(matches==0)
        {
            newlist[count]=list[i];
            count=count+1;
        }
    }

    return newlist;
    
}



