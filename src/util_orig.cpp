/*************************************************
 * util.c                                        *
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

#define PI 3.14159265

using namespace CN_NS;

int ChemNetworkOrig::findf(FILE *fd, int n, ...){
 
  int  s;
  int  ret;
  char buf[1024];
  va_list args;
  char *tem[n];
 
  va_start(args,n);
  for(s=0;s<n;s++)
     tem[s]=va_arg(args,char*);
  va_end(args);
 
  for(s = 0; (ret=fscanf(fd, "%s", buf)) == 1; ){
 
    if( strncmp(buf,tem[s], 1024) == 0 ||
        strcmp("*",tem[s]) == 0) s++; else s = 0;
 
    if( s == n ) break;
  }
 
  return ret;
 }


double ChemNetworkOrig::distance(double *atmS, int i, int j, int *sa, int *sb, int crt)
 {
  double xcomp, ycomp, zcomp, result;

      xcomp = pow((atmS[(i + sa[crt]) * 3 - 3] - atmS[(j + sb[crt]) * 3 - 3]),2);
      ycomp = pow((atmS[(i + sa[crt]) * 3 - 2] - atmS[(j + sb[crt]) * 3 - 2]),2);
      zcomp = pow((atmS[(i + sa[crt]) * 3 - 1] - atmS[(j + sb[crt]) * 3 - 1]),2);

      result = sqrt( xcomp + ycomp + zcomp ); 

  return result;
 }


double ChemNetworkOrig::distanceMix(double *atmS1, double *atmS2, int i, int j, int *sa, int *sb, int crt)
 {
  double xcomp, ycomp, zcomp, result;

      xcomp = pow((atmS1[(i + sa[crt]) * 3 - 3] - atmS2[(j + sb[crt]) * 3 - 3]),2);
      ycomp = pow((atmS1[(i + sa[crt]) * 3 - 2] - atmS2[(j + sb[crt]) * 3 - 2]),2);
      zcomp = pow((atmS1[(i + sa[crt]) * 3 - 1] - atmS2[(j + sb[crt]) * 3 - 1]),2);

      result = sqrt( xcomp + ycomp + zcomp );

  return result;
 }


double ChemNetworkOrig::angle(double dist, double dist2, double hyptns)
 {
//  double PI = 3.14159265;
  double result;

  result = (acos((pow(dist,2) + pow(dist2,2) - pow(hyptns,2)) / (2 * dist * dist2)) / PI) * 180;

  return result;

 }

double ChemNetworkOrig::distanceOxOy(double *atmS, int i, int j, int opos)
 {
  double xcomp, ycomp, zcomp, result;

      xcomp = pow((atmS[(i + opos) * 3 - 3] - atmS[(j + opos) * 3 - 3]),2);
      ycomp = pow((atmS[(i + opos) * 3 - 2] - atmS[(j + opos) * 3 - 2]),2);
      zcomp = pow((atmS[(i + opos) * 3 - 1] - atmS[(j + opos) * 3 - 1]),2);

      result = sqrt( xcomp + ycomp + zcomp );

  return result;
 }

double ChemNetworkOrig::dipole_angle(double *atmT, double *atmS, int i, int j, int *stb, int *sta, int crt, int opos, int h1pos, int h2pos)
{
   double r[4], rmag, p1[4], p2[4], p[4], pmag, dotrp, angle; 
   r[0] = 0; p1[0] = 0; p2[0] = 0; p[0] = 0;

   r[1] = atmS[(j + sta[crt]) * 3 - 3] - atmT[(i + stb[crt]) * 3 - 3];
   r[2] = atmS[(j + sta[crt]) * 3 - 2] - atmT[(i + stb[crt]) * 3 - 2];
   r[3] = atmS[(j + sta[crt]) * 3 - 1] - atmT[(i + stb[crt]) * 3 - 1];
  
   rmag = sqrt((r[1] * r[1]) + (r[2] * r[2]) + (r[3] * r[3]));

   p1[1] = atmS[(j + h1pos) * 3 - 3] - atmS[(j + opos) * 3 - 3]; 
   p1[2] = atmS[(j + h1pos) * 3 - 2] - atmS[(j + opos) * 3 - 2];
   p1[3] = atmS[(j + h1pos) * 3 - 1] - atmS[(j + opos) * 3 - 1];

   p2[1] = atmS[(j + h2pos) * 3 - 3] - atmS[(j + opos) * 3 - 3];     
   p2[2] = atmS[(j + h2pos) * 3 - 2] - atmS[(j + opos) * 3 - 2];
   p2[3] = atmS[(j + h2pos) * 3 - 1] - atmS[(j + opos) * 3 - 1];

   p[1] = p1[1] + p2[1];
   p[2] = p1[2] + p2[2];
   p[3] = p1[3] + p2[3];

   pmag = sqrt((p[1] * p[1]) + (p[2] * p[2]) + (p[3] * p[3])); 

   /* dot product of r and p */

   dotrp = (r[1] * p[1] + r[2] * p[2] + r[3] * p[3]);

   /* dipole angle */

   angle = (acos(dotrp/(rmag*pmag))/PI)*180;
                                   
   if(rmag == 0 || pmag == 0) return 180;
                                   
   else return angle;
 
}

double ChemNetworkOrig::dihedral(double *atmS, int atm1, int atm2, int atm3, int atm4, int opos)
{                                   
  double n1mag, n2mag, dotn1n2, dihAngle;
  double r12[4], r32[4], r23[4], r43[4], n1[4], n2[4];

  r12[0] = 0; r32[0] = 0; r23[0] = 0; r43[0] = 0; n1[0] = 0; n2[0] = 0;
                                   
 /* vectors */
                                   
  r12[1] = atmS[3*(atm1 + opos)-3] - atmS[3*(atm2 + opos)-3];
  r12[2] = atmS[3*(atm1 + opos)-2] - atmS[3*(atm2 + opos)-2];
  r12[3] = atmS[3*(atm1 + opos)-1] - atmS[3*(atm2 + opos)-1];
                                   
  r32[1] = atmS[3*(atm3 + opos)-3] - atmS[3*(atm2 + opos)-3];
  r32[2] = atmS[3*(atm3 + opos)-2] - atmS[3*(atm2 + opos)-2];
  r32[3] = atmS[3*(atm3 + opos)-1] - atmS[3*(atm2 + opos)-1];
                                   
  r23[1] = atmS[3*(atm2 + opos)-3] - atmS[3*(atm3 + opos)-3];
  r23[2] = atmS[3*(atm2 + opos)-2] - atmS[3*(atm3 + opos)-2];
  r23[3] = atmS[3*(atm2 + opos)-1] - atmS[3*(atm3 + opos)-1];
                                   
  r43[1] = atmS[3*(atm4 + opos)-3] - atmS[3*(atm3 + opos)-3];
  r43[2] = atmS[3*(atm4 + opos)-2] - atmS[3*(atm3 + opos)-2];
  r43[3] = atmS[3*(atm4 + opos)-1] - atmS[3*(atm3 + opos)-1];
                                   
  /* normals - cross products */
                                   
  n1[1] = r12[2] * r32[3] - r12[3] * r32[2];
  n1[2] = r12[3] * r32[1] - r12[1] * r32[3];
  n1[3] = r12[1] * r32[2] - r12[2] * r32[1];
                                   
  n2[1] = r23[2] * r43[3] - r23[3] * r43[2];
  n2[2] = r23[3] * r43[1] - r23[1] * r43[3];
  n2[3] = r23[1] * r43[2] - r23[2] * r43[1];
                                   
  n1mag = sqrt((n1[1] * n1[1]) + (n1[2] * n1[2]) + (n1[3] * n1[3]));
  n2mag = sqrt((n2[1] * n2[1]) + (n2[2] * n2[2]) + (n2[3] * n2[3]));
                                   
  /* dot product of n1 and n2 */
                                   
  dotn1n2 = (n1[1] * n2[1] + n1[2] * n2[2] + n1[3] * n2[3]);
                                   
  /* dihedral angle */
                                   
  dihAngle = (acos(dotn1n2/(n1mag*n2mag))/PI)*180;
                                   
  if(n1mag == 0 || n2mag == 0) return 180;
                                   
  else return dihAngle;
                                   
                                   
}

double ChemNetworkOrig::dihedral2(double *atmS, int atm1, int atm2, int atm3, int atm4, int opos)
{                                   
  double n1mag, n2mag, dotn1n2, dihAngle;
  double r12[4], r32[4], r23[4], r43[4], n1[4], n2[4];
  double crossn1n2[4], y, x;

  r12[0] = 0; r32[0] = 0; r23[0] = 0; r43[0] = 0; n1[0] = 0; n2[0] = 0;
                                   
 /* vectors */
                                   
  r12[1] = atmS[3*(atm1 + opos)-3] - atmS[3*(atm2 + opos)-3];
  r12[2] = atmS[3*(atm1 + opos)-2] - atmS[3*(atm2 + opos)-2];
  r12[3] = atmS[3*(atm1 + opos)-1] - atmS[3*(atm2 + opos)-1];
                                   
  r32[1] = atmS[3*(atm3 + opos)-3] - atmS[3*(atm2 + opos)-3];
  r32[2] = atmS[3*(atm3 + opos)-2] - atmS[3*(atm2 + opos)-2];
  r32[3] = atmS[3*(atm3 + opos)-1] - atmS[3*(atm2 + opos)-1];
                                   
  r23[1] = atmS[3*(atm2 + opos)-3] - atmS[3*(atm3 + opos)-3];
  r23[2] = atmS[3*(atm2 + opos)-2] - atmS[3*(atm3 + opos)-2];
  r23[3] = atmS[3*(atm2 + opos)-1] - atmS[3*(atm3 + opos)-1];
                                   
  r43[1] = atmS[3*(atm4 + opos)-3] - atmS[3*(atm3 + opos)-3];
  r43[2] = atmS[3*(atm4 + opos)-2] - atmS[3*(atm3 + opos)-2];
  r43[3] = atmS[3*(atm4 + opos)-1] - atmS[3*(atm3 + opos)-1];
                                   
  /* normals - cross products */
                                   
  n1[1] = r12[2] * r32[3] - r12[3] * r32[2];
  n1[2] = r12[3] * r32[1] - r12[1] * r32[3];
  n1[3] = r12[1] * r32[2] - r12[2] * r32[1];
                                   
  n2[1] = r23[2] * r43[3] - r23[3] * r43[2];
  n2[2] = r23[3] * r43[1] - r23[1] * r43[3];
  n2[3] = r23[1] * r43[2] - r23[2] * r43[1];
                                   
  n1mag = sqrt((n1[1] * n1[1]) + (n1[2] * n1[2]) + (n1[3] * n1[3]));
  n2mag = sqrt((n2[1] * n2[1]) + (n2[2] * n2[2]) + (n2[3] * n2[3]));
                                   
  /* dot product of n1 and n2 */
                                   
  dotn1n2 = (n1[1] * n2[1] + n1[2] * n2[2] + n1[3] * n2[3]);

  /* cross product of n1 and n2 */

  crossn1n2[1] = n1[2] * n2[3] - n1[3] * n2[2];
  crossn1n2[2] = n1[3] * n2[1] - n1[1] * n2[3];
  crossn1n2[3] = n1[1] * n2[2] - n1[2] * n2[1];
                                   
  /* dihedral angle */
                                   
  y = (crossn1n2[1] * r23[1] + crossn1n2[2] * r23[2] + crossn1n2[3] * r23[3]);
                                   
  x = dotn1n2 * sqrt((r23[1] * r23[1]) + (r23[2] * r23[2]) + (r23[3] * r23[3]));
                                   
  dihAngle = (atan(y/x)/PI)*180;

  if(n1mag == 0 || n2mag == 0) return 180;

  else return dihAngle;
                                   
                                   
}


int ChemNetworkOrig::perm6(int v6[], int n, int i, int v2[], int *e){

        int     j,k;

        if(i == n){
               // for (j=0; j<n; j++) printf ("%d ", v[j]); 
               // printf ("\n");

                for(k=0;k<100000;k=k+6){ 
                if(v6[0]==v2[k] && v6[1]==v2[k+1] && v6[2]==v2[k+2] && v6[3]==v2[k+3] && v6[4]==v2[k+4] && v6[5]==v2[k+5]){ *e = *e + 1;}
                }

        } else

                for(j=i; j<n; j++){

                        swap6(v6, i, j);
                        perm6(v6, n, i+1, v2, e);
                        swap6(v6, i, j);
                }
return *e;

}

void ChemNetworkOrig::swap6(int v6[], int i, int j){
    int     t;
    
    t = v6[i];
    v6[i] = v6[j];
    v6[j] = t;
}


int ChemNetworkOrig::perm5(int v5[], int n, int i, int v2[], int *e){

        int     j,k;

        if(i == n){
               // for (j=0; j<n; j++) printf ("%d ", v[j]); 
               // printf ("\n");

                for(k=0;k<100000;k=k+5){ 
                if(v5[0]==v2[k] && v5[1]==v2[k+1] && v5[2]==v2[k+2] && v5[3]==v2[k+3] && v5[4]==v2[k+4]){ *e = *e + 1;}
                }

        } else

                for(j=i; j<n; j++){

                        swap5(v5, i, j);
                        perm5(v5, n, i+1, v2, e);
                        swap5(v5, i, j);
                }
return *e;

}

void ChemNetworkOrig::swap5(int v5[], int i, int j){
    int     t;
    
    t = v5[i];
    v5[i] = v5[j];
    v5[j] = t;
}


int ChemNetworkOrig::perm4(int v4[], int n, int i, int v2[], int *e){

        int     j,k;

        if(i == n){
               // for (j=0; j<n; j++) printf ("%d ", v[j]); 
               // printf ("\n");

                for(k=0;k<100000;k=k+4){ 
                if(v4[0]==v2[k] && v4[1]==v2[k+1] && v4[2]==v2[k+2] && v4[3]==v2[k+3]){ *e = *e + 1;}
                }

        } else

                for(j=i; j<n; j++){

                        swap4(v4, i, j);
                        perm4(v4, n, i+1, v2, e);
                        swap4(v4, i, j);
                }
return *e;

}

void ChemNetworkOrig::swap4(int v4[], int i, int j){
    int     t;
    
    t = v4[i];
    v4[i] = v4[j];
    v4[j] = t;
}

int ChemNetworkOrig::perm3(int v3[], int n, int i, int v2[], int *e){

        int     j,k;

        if(i == n){
               // for (j=0; j<n; j++) printf ("%d ", v[j]); 
               // printf ("\n");

                for(k=0;k<100000;k=k+3){ 
                if(v3[0]==v2[k] && v3[1]==v2[k+1] && v3[2]==v2[k+2]){ *e = *e + 1;}
                }

        } else

                for(j=i; j<n; j++){

                        swap3(v3, i, j);
                        perm3(v3, n, i+1, v2, e);
                        swap3(v3, i, j);
                }
return *e;

}

void ChemNetworkOrig::swap3(int v3[], int i, int j){
    int     t;
    
    t = v3[i];
    v3[i] = v3[j];
    v3[j] = t;
}




