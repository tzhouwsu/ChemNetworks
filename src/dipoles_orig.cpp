/*************************************************
 * dipoles.c                                     *
 *                                               *
 * Author: Abdullah Ozkanlar                     *
 *         abdullah.ozkanlar@wsu.edu             *
 *                                               *
 * A. Clark Research Lab, Chemistry Department   *
 * Washington State University, Pullman/WA 99164 *
 *************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include "chemnetworks_orig.h"

using namespace CN_NS;

// Water - Solvent Graphs and Orientations of Dipole moment vectors of waters wrt solvent molecules.

void ChemNetworkOrig::graph_sAsB_dip(double *atmS1, double *atmS2, int nd1, int nd2, int nsolvent1, int nsolvent2, int nAtomS1, int nAtomS2,
                    int s1s2hbdn, int *s12a, int *s12b, double *s12as12bBD, 
                    int s1s2hban, int *s1s2v1, int *s1s2v2, int *s1s2v3, int *s1s2v4, int *s1s2v5, double *s1s2v6, 
                    int pbc, double xside, double yside, double zside, int opos, int h1pos, int h2pos, int watid, int solid, 
                    FILE *outputfGraphS1S2, FILE *outputfGeodS1S2, FILE *outputfS1S2Dip)
{
  double *atmS2x, *atmS2y, *atmS2z, *atmS2xy, *atmS2yz, *atmS2zx;
  double *atmS2xminy, *atmS2minyz, *atmS2zminx, *atmS2xyz;  
  double *atmS2xyminz, *atmS2minxyz, *atmS2xminyz;

  double *atmS2minx, *atmS2miny, *atmS2minz, *atmS2minxminy, *atmS2minyminz, *atmS2minzminx;
  double *atmS2minxy, *atmS2yminz, *atmS2minzx, *atmS2minxminyminz;
  double *atmS2minxminyz, *atmS2xminyminz, *atmS2minxyminz;

  int nodei, nodej, i, j, crt;
  double dist, dist2, hyptns, ang, dipoleangle; 

  nodei = nd1;

  for(i = 0; i < nAtomS1; i = i + nsolvent1)
  {
    nodej = nd2;
     
    for(j = 0; j < nAtomS2; j = j + nsolvent2)
    {

        for(crt = 0; crt < s1s2hbdn; crt++)
        {

           dist = distanceMix(atmS1, atmS2, i, j, s12a, s12b, crt);

           if(dist < s12as12bBD[crt])
           {
                if(s1s2v1[crt] == 1 && s1s2v2[crt] == 2)
                {
                   dist2 = distance(atmS2, j, j, s1s2v4, s1s2v5, crt);
                   hyptns = distanceMix(atmS1, atmS2, i, j, s1s2v3, s1s2v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(s1s2v1[crt] == 2 && s1s2v2[crt] == 1)
                {
                   dist2 = distance(atmS1, i, i, s1s2v3, s1s2v4, crt);
                   hyptns = distanceMix(atmS1, atmS2, i, j, s1s2v3, s1s2v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(dist < s12as12bBD[crt] && ang > s1s2v6[crt])
                {
                  fprintf(outputfGraphS1S2,"%d\n%d\n",nodei,nodej);
                  fprintf(outputfGeodS1S2,"%d %d 0 0 0\n",nodei,nodej);

                   //Print dipole orienation
                
                 if(solid < watid){
                 dipoleangle = dipole_angle(atmS1, atmS2, i, j, s12a, s12b, crt, opos, h1pos, h2pos);
                 fprintf(outputfS1S2Dip,"Water-Solvent: %d %d  Angle: %.1lf  Distance: %.2lf  WaterAtom-SolventAtom: %d %d\n",nodej,nodei,dipoleangle,dist,s12b[crt],s12a[crt]);
                 }
                 if(solid > watid){
                 dipoleangle = dipole_angle(atmS2, atmS1, j, i, s12b, s12a, crt, opos, h1pos, h2pos);
                 fprintf(outputfS1S2Dip,"Water-Solvent: %d %d  Angle: %.1lf  Distance: %.2lf  WaterAtom-SolventAtom: %d %d\n",nodei,nodej,dipoleangle,dist,s12a[crt],s12b[crt]);
                 }
                }

            }

        } 

    nodej = nodej + 1;
   
    }

  nodei = nodei + 1;

  }

// check PBC and include in the graph if requested

 if(pbc == 1)
 {
   // allocate memory

   atmS2x = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2y = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2z = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2xy = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2yz = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2zx = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2xminy = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2minyz = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2zminx = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2xyz = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2xyminz = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2minxyz = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2xminyz = (double*)malloc((3*nAtomS2)*sizeof(double));   

   atmS2minx = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2miny = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2minz = (double*)malloc((3*nAtomS2)*sizeof(double)); 
   atmS2minxminy = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2minyminz = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2minzminx = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2minxy = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2yminz = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2minzx = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2minxminyminz = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2minxminyz = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2xminyminz = (double*)malloc((3*nAtomS2)*sizeof(double));
   atmS2minxyminz = (double*)malloc((3*nAtomS2)*sizeof(double));

   // get all the PBC Boxes

   for(i=0; i < nAtomS2; i++)
   {
     atmS2x[3*(i+1)-3] = atmS2[3*(i+1)-3] + xside;
     atmS2x[3*(i+1)-2] = atmS2[3*(i+1)-2];
     atmS2x[3*(i+1)-1] = atmS2[3*(i+1)-1];
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2minx[3*(i+1)-3] = atmS2[3*(i+1)-3] - xside;
     atmS2minx[3*(i+1)-2] = atmS2[3*(i+1)-2];
     atmS2minx[3*(i+1)-1] = atmS2[3*(i+1)-1];
   }
   
   for(i=0; i < nAtomS2; i++)
   {
     atmS2y[3*(i+1)-3] = atmS2[3*(i+1)-3];
     atmS2y[3*(i+1)-2] = atmS2[3*(i+1)-2] + yside;
     atmS2y[3*(i+1)-1] = atmS2[3*(i+1)-1];
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2miny[3*(i+1)-3] = atmS2[3*(i+1)-3];
     atmS2miny[3*(i+1)-2] = atmS2[3*(i+1)-2] - yside;
     atmS2miny[3*(i+1)-1] = atmS2[3*(i+1)-1];
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2z[3*(i+1)-3] = atmS2[3*(i+1)-3];
     atmS2z[3*(i+1)-2] = atmS2[3*(i+1)-2];
     atmS2z[3*(i+1)-1] = atmS2[3*(i+1)-1] + zside;
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2minz[3*(i+1)-3] = atmS2[3*(i+1)-3];
     atmS2minz[3*(i+1)-2] = atmS2[3*(i+1)-2];
     atmS2minz[3*(i+1)-1] = atmS2[3*(i+1)-1] - zside;
   }
 
   for(i=0; i < nAtomS2; i++)
   {
     atmS2xy[3*(i+1)-3] = atmS2x[3*(i+1)-3];
     atmS2xy[3*(i+1)-2] = atmS2x[3*(i+1)-2] + yside;
     atmS2xy[3*(i+1)-1] = atmS2x[3*(i+1)-1];
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2minxminy[3*(i+1)-3] = atmS2minx[3*(i+1)-3];
     atmS2minxminy[3*(i+1)-2] = atmS2minx[3*(i+1)-2] - yside;
     atmS2minxminy[3*(i+1)-1] = atmS2minx[3*(i+1)-1];
   }
 
   for(i=0; i < nAtomS2; i++)
   {
     atmS2yz[3*(i+1)-3] = atmS2y[3*(i+1)-3];
     atmS2yz[3*(i+1)-2] = atmS2y[3*(i+1)-2];
     atmS2yz[3*(i+1)-1] = atmS2y[3*(i+1)-1] + zside;
   }   

   for(i=0; i < nAtomS2; i++)
   {
     atmS2minyminz[3*(i+1)-3] = atmS2miny[3*(i+1)-3];
     atmS2minyminz[3*(i+1)-2] = atmS2miny[3*(i+1)-2];
     atmS2minyminz[3*(i+1)-1] = atmS2miny[3*(i+1)-1] - zside;
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2zx[3*(i+1)-3] = atmS2z[3*(i+1)-3] + xside;
     atmS2zx[3*(i+1)-2] = atmS2z[3*(i+1)-2];
     atmS2zx[3*(i+1)-1] = atmS2z[3*(i+1)-1];
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2minzminx[3*(i+1)-3] = atmS2minz[3*(i+1)-3] - xside;
     atmS2minzminx[3*(i+1)-2] = atmS2minz[3*(i+1)-2];
     atmS2minzminx[3*(i+1)-1] = atmS2minz[3*(i+1)-1];
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2xminy[3*(i+1)-3] = atmS2x[3*(i+1)-3];
     atmS2xminy[3*(i+1)-2] = atmS2x[3*(i+1)-2] - yside;
     atmS2xminy[3*(i+1)-1] = atmS2x[3*(i+1)-1];
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2minxy[3*(i+1)-3] = atmS2minx[3*(i+1)-3];
     atmS2minxy[3*(i+1)-2] = atmS2minx[3*(i+1)-2] + yside;
     atmS2minxy[3*(i+1)-1] = atmS2minx[3*(i+1)-1];
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2minyz[3*(i+1)-3] = atmS2z[3*(i+1)-3];
     atmS2minyz[3*(i+1)-2] = atmS2z[3*(i+1)-2] - yside;
     atmS2minyz[3*(i+1)-1] = atmS2z[3*(i+1)-1];
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2yminz[3*(i+1)-3] = atmS2y[3*(i+1)-3];
     atmS2yminz[3*(i+1)-2] = atmS2y[3*(i+1)-2];
     atmS2yminz[3*(i+1)-1] = atmS2y[3*(i+1)-1] - zside;
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2zminx[3*(i+1)-3] = atmS2z[3*(i+1)-3] - xside;
     atmS2zminx[3*(i+1)-2] = atmS2z[3*(i+1)-2];
     atmS2zminx[3*(i+1)-1] = atmS2z[3*(i+1)-1];
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2minzx[3*(i+1)-3] = atmS2minz[3*(i+1)-3] + xside;
     atmS2minzx[3*(i+1)-2] = atmS2minz[3*(i+1)-2];
     atmS2minzx[3*(i+1)-1] = atmS2minz[3*(i+1)-1];
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2xyz[3*(i+1)-3] = atmS2xy[3*(i+1)-3];
     atmS2xyz[3*(i+1)-2] = atmS2xy[3*(i+1)-2];
     atmS2xyz[3*(i+1)-1] = atmS2xy[3*(i+1)-1] + zside;
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2minxminyminz[3*(i+1)-3] = atmS2minxminy[3*(i+1)-3];
     atmS2minxminyminz[3*(i+1)-2] = atmS2minxminy[3*(i+1)-2];
     atmS2minxminyminz[3*(i+1)-1] = atmS2minxminy[3*(i+1)-1] - zside;
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2xyminz[3*(i+1)-3] = atmS2xy[3*(i+1)-3];
     atmS2xyminz[3*(i+1)-2] = atmS2xy[3*(i+1)-2];
     atmS2xyminz[3*(i+1)-1] = atmS2xy[3*(i+1)-1] - zside;
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2minxminyz[3*(i+1)-3] = atmS2minxminy[3*(i+1)-3];
     atmS2minxminyz[3*(i+1)-2] = atmS2minxminy[3*(i+1)-2];
     atmS2minxminyz[3*(i+1)-1] = atmS2minxminy[3*(i+1)-1] + zside;
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2minxyz[3*(i+1)-3] = atmS2yz[3*(i+1)-3] - xside;
     atmS2minxyz[3*(i+1)-2] = atmS2yz[3*(i+1)-2];
     atmS2minxyz[3*(i+1)-1] = atmS2yz[3*(i+1)-1];
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2xminyminz[3*(i+1)-3] = atmS2minyminz[3*(i+1)-3] + xside;
     atmS2xminyminz[3*(i+1)-2] = atmS2minyminz[3*(i+1)-2];
     atmS2xminyminz[3*(i+1)-1] = atmS2minyminz[3*(i+1)-1];
   }
   
   for(i=0; i < nAtomS2; i++)
   {
     atmS2xminyz[3*(i+1)-3] = atmS2minyz[3*(i+1)-3] + xside;
     atmS2xminyz[3*(i+1)-2] = atmS2minyz[3*(i+1)-2];
     atmS2xminyz[3*(i+1)-1] = atmS2minyz[3*(i+1)-1];
   }

   for(i=0; i < nAtomS2; i++)
   {
     atmS2minxyminz[3*(i+1)-3] = atmS2minxy[3*(i+1)-3];
     atmS2minxyminz[3*(i+1)-2] = atmS2minxy[3*(i+1)-2];
     atmS2minxyminz[3*(i+1)-1] = atmS2minxy[3*(i+1)-1] - zside;
   }

   search_pbc_gsAsB_dip(1, atmS1, atmS2x, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(-1, atmS1, atmS2minx, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(2, atmS1, atmS2y, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(-2, atmS1, atmS2miny, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(3, atmS1, atmS2z, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(-3, atmS1, atmS2minz, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(4, atmS1, atmS2xy, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(-4, atmS1, atmS2minxminy, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(5, atmS1, atmS2yz, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(-5, atmS1, atmS2minyminz, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(6, atmS1, atmS2zx, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(-6, atmS1, atmS2minzminx, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(7, atmS1, atmS2xminy, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(-7, atmS1, atmS2minxy, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(8, atmS1, atmS2minyz, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(-8, atmS1, atmS2yminz, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(9, atmS1, atmS2zminx, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(-9, atmS1, atmS2minzx, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(10, atmS1, atmS2xyz, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(-10, atmS1, atmS2minxminyminz, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(11, atmS1, atmS2xyminz, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(-11, atmS1, atmS2minxminyz, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(12, atmS1, atmS2minxyz, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(-12, atmS1, atmS2xminyminz, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(13, atmS1, atmS2xminyz, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                       s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

   search_pbc_gsAsB_dip(-13, atmS1, atmS2minxyminz, nd1, nd2, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
                        s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip);

    // free memory

    free(atmS2x);
    free(atmS2y);
    free(atmS2z);
    free(atmS2xy);
    free(atmS2yz);
    free(atmS2zx);
    free(atmS2xminy);
    free(atmS2minyz);
    free(atmS2zminx);
    free(atmS2xyz);
    free(atmS2xyminz);
    free(atmS2minxyz);
    free(atmS2xminyz);

    free(atmS2minx);
    free(atmS2miny);
    free(atmS2minz);
    free(atmS2minxminy);
    free(atmS2minyminz);
    free(atmS2minzminx);
    free(atmS2minxy);
    free(atmS2yminz);
    free(atmS2minzx);
    free(atmS2minxminyminz);
    free(atmS2minxminyz);
    free(atmS2xminyminz);
    free(atmS2minxyminz);

 } // Close if(pbc == 1)...

}

// search for bonding through periodic boundary conditions (water-solvent)

void ChemNetworkOrig::search_pbc_gsAsB_dip(int boxid, double *atmS1, double *atmS2x, int nd1, int nd2, int nsolvent1, int nsolvent2, int nAtomS1, int nAtomS2, int s1s2hbdn,
                          int *s12a, int *s12b, double *s12as12bBD, int s1s2hban, int *s1s2v1, int *s1s2v2, int *s1s2v3, int *s1s2v4, int *s1s2v5, double *s1s2v6,
                          int opos, int h1pos, int h2pos, int watid, int solid,
                          FILE *outputfGraphS1S2, FILE *outputfGeodS1S2, FILE *outputfS1S2Dip)
{
  int nodei, nodej, i, j, crt;
  double dist, dist2, hyptns, ang, dipoleangle; 

  nodei = nd1;

  for(i = 0; i < nAtomS1; i = i + nsolvent1)
  {
    nodej = nd2;
     
    for(j = 0; j < nAtomS2; j = j + nsolvent2)
    {

        for(crt = 0; crt < s1s2hbdn; crt++)
        {

           dist = distanceMix(atmS1, atmS2x, i, j, s12a, s12b, crt);

           if(dist < s12as12bBD[crt])
           {
                if(s1s2v1[crt] == 1 && s1s2v2[crt] == 2)
                {
                   dist2 = distance(atmS2x, j, j, s1s2v4, s1s2v5, crt);
                   hyptns = distanceMix(atmS1, atmS2x, i, j, s1s2v3, s1s2v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(s1s2v1[crt] == 2 && s1s2v2[crt] == 1)
                {
                   dist2 = distance(atmS1, i, i, s1s2v3, s1s2v4, crt);
                   hyptns = distanceMix(atmS1, atmS2x, i, j, s1s2v3, s1s2v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(dist < s12as12bBD[crt] && ang > s1s2v6[crt])
                {
                   fprintf(outputfGraphS1S2,"%d\n%d\n",nodei,nodej);

                   if(boxid == 1) fprintf(outputfGeodS1S2,"%d %d 1 0 0\n",nodei,nodej);
                   if(boxid == -1) fprintf(outputfGeodS1S2,"%d %d -1 0 0\n",nodei,nodej);
                   if(boxid == 2) fprintf(outputfGeodS1S2,"%d %d 0 1 0\n",nodei,nodej);
                   if(boxid == -2) fprintf(outputfGeodS1S2,"%d %d 0 -1 0\n",nodei,nodej);
                   if(boxid == 3) fprintf(outputfGeodS1S2,"%d %d 0 0 1\n",nodei,nodej);
                   if(boxid == -3) fprintf(outputfGeodS1S2,"%d %d 0 0 -1\n",nodei,nodej);
                   if(boxid == 4) fprintf(outputfGeodS1S2,"%d %d 1 1 0\n",nodei,nodej);
                   if(boxid == -4) fprintf(outputfGeodS1S2,"%d %d -1 -1 0\n",nodei,nodej);
                   if(boxid == 5) fprintf(outputfGeodS1S2,"%d %d 0 1 1\n",nodei,nodej);
                   if(boxid == -5) fprintf(outputfGeodS1S2,"%d %d 0 -1 -1\n",nodei,nodej);
                   if(boxid == 6) fprintf(outputfGeodS1S2,"%d %d 1 0 1\n",nodei,nodej);
                   if(boxid == -6) fprintf(outputfGeodS1S2,"%d %d -1 0 -1\n",nodei,nodej);
                   if(boxid == 7) fprintf(outputfGeodS1S2,"%d %d 1 -1 0\n",nodei,nodej);
                   if(boxid == -7) fprintf(outputfGeodS1S2,"%d %d -1 1 0\n",nodei,nodej);
                   if(boxid == 8) fprintf(outputfGeodS1S2,"%d %d 0 -1 1\n",nodei,nodej);
                   if(boxid == -8) fprintf(outputfGeodS1S2,"%d %d 0 1 -1\n",nodei,nodej);
                   if(boxid == 9) fprintf(outputfGeodS1S2,"%d %d -1 0 1\n",nodei,nodej);
                   if(boxid == -9) fprintf(outputfGeodS1S2,"%d %d 1 0 -1\n",nodei,nodej);
                   if(boxid == 10) fprintf(outputfGeodS1S2,"%d %d 1 1 1\n",nodei,nodej);
                   if(boxid == -10) fprintf(outputfGeodS1S2,"%d %d -1 -1 -1\n",nodei,nodej);
                   if(boxid == 11) fprintf(outputfGeodS1S2,"%d %d 1 1 -1\n",nodei,nodej);
                   if(boxid == -11) fprintf(outputfGeodS1S2,"%d %d -1 -1 1\n",nodei,nodej);
                   if(boxid == 12) fprintf(outputfGeodS1S2,"%d %d -1 1 1\n",nodei,nodej);
                   if(boxid == -12) fprintf(outputfGeodS1S2,"%d %d 1 -1 -1\n",nodei,nodej);
                   if(boxid == 13) fprintf(outputfGeodS1S2,"%d %d 1 -1 1\n",nodei,nodej);
                   if(boxid == -13) fprintf(outputfGeodS1S2,"%d %d -1 1 -1\n",nodei,nodej);
                  
                   //Print dipole orienation
                
                 if(solid < watid){
                 dipoleangle = dipole_angle(atmS1, atmS2x, i, j, s12a, s12b, crt, opos, h1pos, h2pos);
                 fprintf(outputfS1S2Dip,"Water-Solvent: %d %d  Angle: %.1lf  Distance: %.2lf  WaterAtom-SolventAtom: %d %d\n",nodej,nodei,dipoleangle,dist,s12b[crt],s12a[crt]);
                 }
                 if(solid > watid){
                 dipoleangle = dipole_angle(atmS2x, atmS1, j, i, s12b, s12a, crt, opos, h1pos, h2pos);
                 fprintf(outputfS1S2Dip,"Water-Solvent: %d %d  Angle: %.1lf  Distance: %.2lf  WaterAtom-SolventAtom: %d %d\n",nodei,nodej,dipoleangle,dist,s12a[crt],s12b[crt]);
                 }
 
                }

            }

        } 

    nodej = nodej + 1;
   
    }

  nodei = nodei + 1;

  }


}



// Water - Solute Graphs and Orientations of Dipole moment vectors of waters within the cutoff distance.

void ChemNetworkOrig::graph_st_dip(double *atmS1, double *atmT1, int nd1, int nd2, int nsolvent1, int nsolute1, int nAtomS1, int nAtomT1, int s1t1cutoffnum, int *s1t1a, int *s1t1b,
                  double *s1t1cutoff, int pbc, double xside, double yside, double zside, int opos, int h1pos, int h2pos, 
                  FILE *outputfGraphS1T1, FILE *outputfGeodS1T1, FILE *outputfS1T1Dip)
{
  double *atmS1x, *atmS1y, *atmS1z, *atmS1xy, *atmS1yz, *atmS1zx;
  double *atmS1xminy, *atmS1minyz, *atmS1zminx, *atmS1xyz;  
  double *atmS1xyminz, *atmS1minxyz, *atmS1xminyz;

  double *atmS1minx, *atmS1miny, *atmS1minz, *atmS1minxminy, *atmS1minyminz, *atmS1minzminx;
  double *atmS1minxy, *atmS1yminz, *atmS1minzx, *atmS1minxminyminz;
  double *atmS1minxminyz, *atmS1xminyminz, *atmS1minxyminz;

  int nodei, nodej, i, j, crt;
  double dist, dipoleangle; 

  nodei = nd2;

  for(i = 0; i < nAtomT1; i = i + nsolute1)
  {
    nodej = nd1;
     
    for(j = 0; j < nAtomS1; j = j + nsolvent1)
    {

        for(crt = 0; crt < s1t1cutoffnum; crt++)
        {

           dist = distanceMix(atmT1, atmS1, i, j, s1t1b, s1t1a, crt);

           if(dist < s1t1cutoff[crt])
           {    
                  fprintf(outputfGraphS1T1,"%d\n%d\n",nodej,nodei);
                  fprintf(outputfGeodS1T1,"%d %d 0 0 0\n",nodej,nodei);

                  //Print dipole orienation
                  dipoleangle = dipole_angle(atmT1, atmS1, i, j, s1t1b, s1t1a, crt, opos, h1pos, h2pos);
                  fprintf(outputfS1T1Dip,"Water-Solute: %d %d  Angle: %.1lf  Distance: %.2lf  WaterAtom-SoluteAtom: %d %d\n",nodej,nodei,dipoleangle,dist,s1t1a[crt],s1t1b[crt]);

            }

        } 

    nodej = nodej + 1;
   
    }

  nodei = nodei + 1;

  }

// check PBC and include in the graph if requested

 if(pbc == 1)
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

   search_pbc_gst_dip(1, atmT1, atmS1x, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(-1, atmT1, atmS1minx, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(2, atmT1, atmS1y, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(-2, atmT1, atmS1miny, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(3, atmT1, atmS1z, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(-3, atmT1, atmS1minz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(4, atmT1, atmS1xy, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(-4, atmT1, atmS1minxminy, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(5, atmT1, atmS1yz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(-5, atmT1, atmS1minyminz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(6, atmT1, atmS1zx, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(-6, atmT1, atmS1minzminx, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(7, atmT1, atmS1xminy, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(-7, atmT1, atmS1minxy, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(8, atmT1, atmS1minyz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(-8, atmT1, atmS1yminz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(9, atmT1, atmS1zminx, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(-9, atmT1, atmS1minzx, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(10, atmT1, atmS1xyz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(-10, atmT1, atmS1minxminyminz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(11, atmT1, atmS1xyminz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(-11, atmT1, atmS1minxminyz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(12, atmT1, atmS1minxyz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(-12, atmT1, atmS1xminyminz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(13, atmT1, atmS1xminyz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);

   search_pbc_gst_dip(-13, atmT1, atmS1minxyminz, nd1, nd2, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
                      opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);


    // free memory

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

 } // Close if(pbc == 1)...

}


void ChemNetworkOrig::search_pbc_gst_dip(int boxid, double *atmT1, double *atmS1x, int nd1, int nd2, int nsolvent1, int nsolute1, int nAtomS1, int nAtomT1, int s1t1cutoffnum, 
                        int *s1t1a, int *s1t1b, double *s1t1cutoff, int opos, int h1pos, int h2pos, 
                        FILE *outputfGraphS1T1, FILE *outputfGeodS1T1, FILE *outputfS1T1Dip)
{
  int nodei, nodej, i, j, crt;
  double dist, dipoleangle; 

  nodei = nd2;

  for(i = 0; i < nAtomT1; i = i + nsolute1)
  {
    nodej = nd1;
     
    for(j = 0; j < nAtomS1; j = j + nsolvent1)
    {

        for(crt = 0; crt < s1t1cutoffnum; crt++)
        {

           dist = distanceMix(atmT1, atmS1x, i, j, s1t1b, s1t1a, crt);

           if(dist < s1t1cutoff[crt])
           {
                   fprintf(outputfGraphS1T1,"%d\n%d\n",nodej,nodei);

                   if(boxid == 1) fprintf(outputfGeodS1T1,"%d %d 1 0 0\n",nodej,nodei);
                   if(boxid == -1) fprintf(outputfGeodS1T1,"%d %d -1 0 0\n",nodej,nodei);
                   if(boxid == 2) fprintf(outputfGeodS1T1,"%d %d 0 1 0\n",nodej,nodei);
                   if(boxid == -2) fprintf(outputfGeodS1T1,"%d %d 0 -1 0\n",nodej,nodei);
                   if(boxid == 3) fprintf(outputfGeodS1T1,"%d %d 0 0 1\n",nodej,nodei);
                   if(boxid == -3) fprintf(outputfGeodS1T1,"%d %d 0 0 -1\n",nodej,nodei);
                   if(boxid == 4) fprintf(outputfGeodS1T1,"%d %d 1 1 0\n",nodej,nodei);
                   if(boxid == -4) fprintf(outputfGeodS1T1,"%d %d -1 -1 0\n",nodej,nodei);
                   if(boxid == 5) fprintf(outputfGeodS1T1,"%d %d 0 1 1\n",nodej,nodei);
                   if(boxid == -5) fprintf(outputfGeodS1T1,"%d %d 0 -1 -1\n",nodej,nodei);
                   if(boxid == 6) fprintf(outputfGeodS1T1,"%d %d 1 0 1\n",nodej,nodei);
                   if(boxid == -6) fprintf(outputfGeodS1T1,"%d %d -1 0 -1\n",nodej,nodei);
                   if(boxid == 7) fprintf(outputfGeodS1T1,"%d %d 1 -1 0\n",nodej,nodei);
                   if(boxid == -7) fprintf(outputfGeodS1T1,"%d %d -1 1 0\n",nodej,nodei);
                   if(boxid == 8) fprintf(outputfGeodS1T1,"%d %d 0 -1 1\n",nodej,nodei);
                   if(boxid == -8) fprintf(outputfGeodS1T1,"%d %d 0 1 -1\n",nodej,nodei);
                   if(boxid == 9) fprintf(outputfGeodS1T1,"%d %d -1 0 1\n",nodej,nodei);
                   if(boxid == -9) fprintf(outputfGeodS1T1,"%d %d 1 0 -1\n",nodej,nodei);
                   if(boxid == 10) fprintf(outputfGeodS1T1,"%d %d 1 1 1\n",nodej,nodei);
                   if(boxid == -10) fprintf(outputfGeodS1T1,"%d %d -1 -1 -1\n",nodej,nodei);
                   if(boxid == 11) fprintf(outputfGeodS1T1,"%d %d 1 1 -1\n",nodej,nodei);
                   if(boxid == -11) fprintf(outputfGeodS1T1,"%d %d -1 -1 1\n",nodej,nodei);
                   if(boxid == 12) fprintf(outputfGeodS1T1,"%d %d -1 1 1\n",nodej,nodei);
                   if(boxid == -12) fprintf(outputfGeodS1T1,"%d %d 1 -1 -1\n",nodej,nodei);
                   if(boxid == 13) fprintf(outputfGeodS1T1,"%d %d 1 -1 1\n",nodej,nodei);
                   if(boxid == -13) fprintf(outputfGeodS1T1,"%d %d -1 1 -1\n",nodej,nodei);

                  //Print dipole orienation
                  dipoleangle = dipole_angle(atmT1, atmS1x, i, j, s1t1b, s1t1a, crt, opos, h1pos, h2pos);
                  fprintf(outputfS1T1Dip,"Water-Solute: %d %d  Angle: %.1lf  Distance: %.2lf  WaterAtom-SoluteAtom: %d %d\n",nodej,nodei,dipoleangle,dist,s1t1a[crt],s1t1b[crt]);

            }

        } 

    nodej = nodej + 1;
   
    }

  nodei = nodei + 1;

  }


}


