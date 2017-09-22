/*************************************************
 * structures.c                                  *
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

#include "chemnetworks_new.h"

using namespace CN_NS;

// Water Hexamer Ring search

void ChemNetworkNew::hringsearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                 int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                 FILE *outputfHRing, char *foutputHRingIso)
{

   FILE *outputfHRingIso;
   int nodei, i, j, crt;
   double dist, dist2, hyptns, ang;
   int atma, atmb, atmc, atmd, atme, atmf, bonded; 
   int ok[MAX][7], cntr, in, v2[100000], v6[6], en, cek, e;
   double distOaOb, distOaOc, distObOc, distOcOd, distObOd, distOcOe, distOdOe, distOdOf, distOeOf, distOeOa, distOfOa;
   double angleOaObOc, angleObOcOd, angleOcOdOe, angleOdOeOf, angleOeOfOa;

   nodei = nd1;

   if(hiso == 1) outputfHRingIso = fopen(foutputHRingIso,"w"); 
  
   for(i = 0; i < MAX; i++)
       for(j = 0; j < 7; j++)
           ok[i][j] = 0;
    i = 0;          
    cntr = 0; for(in = 0; in < 100000; in++) v2[in] = 0; for(in = 0; in < 6; in++) v6[in] = 0; 


    for(atma = 0; atma <= (nAtomS1 - 3); atma = atma + 3){ i = i + 1; cek = 0;
        for(en = 1; en <= i; en++){
            if(atma != ok[en][1] && atma != ok[en][2] && atma != ok[en][3] && atma != ok[en][4] && atma != ok[en][5] && atma != ok[en][6])
            {
              cek = cek + 1;
            }
        }

    if(cek == i){
                   
    for(atmb = 0; atmb <= (nAtomS1 - 3); atmb = atmb + 3)
    {
        if(atmb != atma){
        
        bonded = 0;
        for(crt = 0; crt < s1s1hbdn; crt++)
        {

           dist = distance(atmS1, atma, atmb, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atma, atma, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        } 
        if(bonded != 0){

    for(atmc = 0; atmc <= (nAtomS1 - 3); atmc = atmc + 3) 
    {
      if(atmc != atma && atmc != atmb){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmb, atmc, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){
        
        distOaOb = distanceOxOy(atmS1, atma, atmb, opos);
        distOaOc = distanceOxOy(atmS1, atma, atmc, opos);
        distObOc = distanceOxOy(atmS1, atmb, atmc, opos);
        angleOaObOc = (acos((pow(distOaOb,2)+pow(distObOc,2)-pow(distOaOc,2))/(2*distOaOb*distObOc))/PI)*180;
        if(angleOaObOc < 130 && angleOaObOc > 110){

    for(atmd = 0; atmd <= (nAtomS1 - 3); atmd = atmd + 3)
    {
      if(atmd != atma && atmd != atmb && atmd != atmc){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmc, atmd, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmd, atmd, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){

        distOcOd = distanceOxOy(atmS1, atmc, atmd, opos);
        distObOd = distanceOxOy(atmS1, atmb, atmd, opos);
        angleObOcOd = (acos((pow(distObOc,2)+pow(distOcOd,2)-pow(distObOd,2))/(2*distObOc*distOcOd))/PI)*180;

        if(angleObOcOd < 130 && angleObOcOd > 110){

        if(fabs(dihedral(atmS1, atma, atmb, atmc, atmd, opos)) <= 30){

        for(atme = 0; atme <= (nAtomS1-3); atme = atme + 3)
        {
          if(atme != atma && atme != atmb && atme != atmc && atme != atmd){

          bonded = 0;
          for(crt = 0; crt < s1s1hbdn; crt++)
          {

              dist = distance(atmS1, atmd, atme, s1a, s1b, crt);

              if(dist < s1as1bBD[crt])
              {
                 if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                 {
                    dist2 = distance(atmS1, atme, atme, s1s1v4, s1s1v5, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                 {
                    dist2 = distance(atmS1, atmd, atmd, s1s1v3, s1s1v4, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                 {
                    bonded = bonded + 1;
                 }

             }

          }
          if(bonded != 0){
 
          distOcOe = distanceOxOy(atmS1, atmc, atme, opos);
          distOdOe = distanceOxOy(atmS1, atmd, atme, opos);
          angleOcOdOe = (acos((pow(distOcOd,2)+pow(distOdOe,2)-pow(distOcOe,2))/(2*distOcOd*distOdOe))/PI)*180;

          if(angleOcOdOe < 130 && angleOcOdOe > 110){

          for(atmf = 0; atmf <= (nAtomS1 - 3); atmf = atmf + 3)
          {
            if(atmf != atma && atmf != atmb && atmf != atmc && atmf != atmd && atmf != atme){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atme, atmf, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atme, atmf, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atme, atme, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atme, atmf, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){
     
            distOdOf = distanceOxOy(atmS1, atmd, atmf, opos);
            distOeOf = distanceOxOy(atmS1, atme, atmf, opos);
            angleOdOeOf = (acos((pow(distOdOe,2)+pow(distOeOf,2)-pow(distOdOf,2))/(2*distOdOe*distOeOf))/PI)*180;

            if(angleOdOeOf < 130 && angleOdOeOf > 110){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmf, atma, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atma, atma, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmf, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmf, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){

            distOeOa = distanceOxOy(atmS1, atme, atma, opos);
            distOfOa = distanceOxOy(atmS1, atmf, atma, opos);
            angleOeOfOa = (acos((pow(distOeOf,2)+pow(distOfOa,2)-pow(distOeOa,2))/(2*distOeOf*distOfOa))/PI)*180;

            if(angleOeOfOa < 130 && angleOeOfOa > 110){

            if(fabs(dihedral(atmS1,atma,atmf,atme,atmd, opos)) <= 30){     

            if(fabs(dihedral(atmS1,atmf,atma,atmd,atmc, opos)) >= 150 && fabs(dihedral(atmS1,atmb,atma,atmd,atme, opos)) >= 150){      

            ok[i][1] = atma; ok[i][2] = atmb; ok[i][3] = atmc; ok[i][4] = atmd; ok[i][5] = atme; ok[i][6] = atmf;

            v6[0] = atma; v6[1] = atmb; v6[2] = atmc; v6[3] = atmd; v6[4] = atme; v6[5] = atmf; 

            v2[cntr] = atma; v2[cntr+1] = atmb; v2[cntr+2] = atmc; v2[cntr+3] = atmd; v2[cntr+4] = atme; v2[cntr+5] = atmf;

            cntr = cntr + 6;

            e = 0;

            if(perm6(v6, 6, 0, v2, &e) == 1){ fprintf(outputfHRing,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);

             
            if(hiso == 1)
            { 
              if(check_isolated_hexamer(atmS1,nAtomS1,s1s1hbdn,s1a,s1b,s1as1bBD,s1s1hban,s1s1v1,s1s1v2,s1s1v3,s1s1v4,s1s1v5,s1s1v6,atma,atmb,atmc,atmd,atme,atmf) == 1)
         
                 fprintf(outputfHRingIso,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);
            }

           }
           }
           }
           }
           }
           }
          }
          }
          }
          }
          }
         }
         }
         }
         }
        }
        }
        }
       }
       }
      }
      }
     }
     }
     }
     }
    }
    
if(hiso == 1) fclose(outputfHRingIso);
}


// Water Hexamer Book search

void ChemNetworkNew::hbooksearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                 int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                 FILE *outputfHBook, char *foutputHBookIso)
{

   FILE *outputfHBookIso;
   int nodei, i, j, crt;
   double dist, dist2, hyptns, ang;
   int atma, atmb, atmc, atmd, atme, atmf, bonded; 
   int ok[MAX][7], cntr, in, v2[100000], v6[6], en, cek, e;
   double distOaOb, distOaOc, distObOc, distOcOd, distObOd, distOcOe, distOdOe, distOdOf, distOeOf, distOeOa, distOfOa;
   double angleOaObOc, angleObOcOd, angleOcOdOe, angleOdOeOf, angleOeOfOa;

   nodei = nd1;

   if(hiso == 1) outputfHBookIso = fopen(foutputHBookIso,"w"); 
  
   for(i = 0; i < MAX; i++)
       for(j = 0; j < 7; j++)
           ok[i][j] = 0;
    i = 0;          
    cntr = 0; for(in = 0; in < 100000; in++) v2[in] = 0; for(in = 0; in < 6; in++) v6[in] = 0; 

    for(atma = 0; atma <= (nAtomS1 - 3); atma = atma + 3){ i = i + 1; cek = 0;
        for(en = 1; en <= i; en++){
            if(atma != ok[en][1] && atma != ok[en][2] && atma != ok[en][3] && atma != ok[en][4] && atma != ok[en][5] && atma != ok[en][6])
            {
              cek = cek + 1;
            }
        }

    if(cek == i){

    for(atmb = 0; atmb <= (nAtomS1 - 3); atmb = atmb + 3)
    {
        if(atmb != atma){
        
        bonded = 0;
        for(crt = 0; crt < s1s1hbdn; crt++)
        {

           dist = distance(atmS1, atma, atmb, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atma, atma, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        } 
        if(bonded != 0){

    for(atmc = 0; atmc <= (nAtomS1 - 3); atmc = atmc + 3) 
    {
      if(atmc != atma && atmc != atmb){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmb, atmc, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){
        
        distOaOb = distanceOxOy(atmS1, atma, atmb, opos);
        distOaOc = distanceOxOy(atmS1, atma, atmc, opos);
        distObOc = distanceOxOy(atmS1, atmb, atmc, opos);
        angleOaObOc = (acos((pow(distOaOb,2)+pow(distObOc,2)-pow(distOaOc,2))/(2*distOaOb*distObOc))/PI)*180;
        if(angleOaObOc < 100 && angleOaObOc > 80){

    for(atmd = 0; atmd <= (nAtomS1 - 3); atmd = atmd + 3)
    {
      if(atmd != atma && atmd != atmb && atmd != atmc){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmc, atmd, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmd, atmd, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){

        distOcOd = distanceOxOy(atmS1, atmc, atmd, opos);
        distObOd = distanceOxOy(atmS1, atmb, atmd, opos);
        angleObOcOd = (acos((pow(distObOc,2)+pow(distOcOd,2)-pow(distObOd,2))/(2*distObOc*distOcOd))/PI)*180;

        if(angleObOcOd < 160 && angleObOcOd > 80){


        for(atme = 0; atme <= (nAtomS1-3); atme = atme + 3)
        {
          if(atme != atma && atme != atmb && atme != atmc && atme != atmd){

          bonded = 0;
          for(crt = 0; crt < s1s1hbdn; crt++)
          {

              dist = distance(atmS1, atmd, atme, s1a, s1b, crt);

              if(dist < s1as1bBD[crt])
              {
                 if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                 {
                    dist2 = distance(atmS1, atme, atme, s1s1v4, s1s1v5, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                 {
                    dist2 = distance(atmS1, atmd, atmd, s1s1v3, s1s1v4, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                 {
                    bonded = bonded + 1;
                 }

             }

          }
          if(bonded != 0){
 
          distOcOe = distanceOxOy(atmS1, atmc, atme, opos);
          distOdOe = distanceOxOy(atmS1, atmd, atme, opos);
          angleOcOdOe = (acos((pow(distOcOd,2)+pow(distOdOe,2)-pow(distOcOe,2))/(2*distOcOd*distOdOe))/PI)*180;

          if(angleOcOdOe < 100 && angleOcOdOe > 80){

          for(atmf = 0; atmf <= (nAtomS1 - 3); atmf = atmf + 3)
          {
            if(atmf != atma && atmf != atmb && atmf != atmc && atmf != atmd && atmf != atme){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atme, atmf, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atme, atmf, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atme, atme, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atme, atmf, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){
     
            distOdOf = distanceOxOy(atmS1, atmd, atmf, opos);
            distOeOf = distanceOxOy(atmS1, atme, atmf, opos);
            angleOdOeOf = (acos((pow(distOdOe,2)+pow(distOeOf,2)-pow(distOdOf,2))/(2*distOdOe*distOeOf))/PI)*180;

            if(angleOdOeOf < 100 && angleOdOeOf > 80){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmf, atma, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atma, atma, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmf, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmf, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){

            distOeOa = distanceOxOy(atmS1, atme, atma, opos);
            distOfOa = distanceOxOy(atmS1, atmf, atma, opos);
            angleOeOfOa = (acos((pow(distOeOf,2)+pow(distOfOa,2)-pow(distOeOa,2))/(2*distOeOf*distOfOa))/PI)*180;

            if(angleOeOfOa < 160 && angleOeOfOa > 80){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmf, atmc, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atmc, atmc, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmf, atmc, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmf, atmc, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){     

            if(fabs(dihedral(atmS1,atma,atmb,atmc,atmf, opos)) <= 20 && fabs(dihedral(atmS1,atmc,atmd,atme,atmf, opos)) <= 20){      

            ok[i][1] = atma; ok[i][2] = atmb; ok[i][3] = atmc; ok[i][4] = atmd; ok[i][5] = atme; ok[i][6] = atmf;

            v6[0] = atma; v6[1] = atmb; v6[2] = atmc; v6[3] = atmd; v6[4] = atme; v6[5] = atmf; 

            v2[cntr] = atma; v2[cntr+1] = atmb; v2[cntr+2] = atmc; v2[cntr+3] = atmd; v2[cntr+4] = atme; v2[cntr+5] = atmf;

            cntr = cntr + 6;

            e = 0;

            if(perm6(v6, 6, 0, v2, &e) == 1){ fprintf(outputfHBook,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);

             
            if(hiso == 1)
            { 
              if(check_isolated_hexamer(atmS1,nAtomS1,s1s1hbdn,s1a,s1b,s1as1bBD,s1s1hban,s1s1v1,s1s1v2,s1s1v3,s1s1v4,s1s1v5,s1s1v6,atma,atmb,atmc,atmd,atme,atmf) == 1)
         
                 fprintf(outputfHBookIso,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);
            }

           }
           }
           }
           }
           }
           }
          }
          }
          }
          }
          }
         }
         }
         }
         }
        }
        }
        }
       }
      }
     }
    }
   }
  }
 }
 }
if(hiso == 1) fclose(outputfHBookIso);
}


// Water Hexamer Prism search

void ChemNetworkNew::hprismsearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                  int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                  FILE *outputfHPrism, char *foutputHPrismIso)
{

   FILE *outputfHPrismIso;
   int nodei, i, j, crt;
   double dist, dist2, hyptns, ang;
   int atma, atmb, atmc, atmd, atme, atmf, bonded; 
   int ok[MAX][7], cntr, in, v2[100000], v6[6], en, cek, e;
   double distOaOb, distOaOc, distObOc, distOcOd, distObOd, distOcOe, distOdOe, distOdOf, distOeOf;
   double angleOaObOc, angleObOcOd, angleOcOdOe, angleOdOeOf, distOeOb, distOfOb, angleOeOfOb;

   nodei = nd1;

   if(hiso == 1) outputfHPrismIso = fopen(foutputHPrismIso,"w"); 
  
   for(i = 0; i < MAX; i++)
       for(j = 0; j < 7; j++)
           ok[i][j] = 0;
    i = 0;          
    cntr = 0; for(in = 0; in < 100000; in++) v2[in] = 0; for(in = 0; in < 6; in++) v6[in] = 0; 

    for(atma = 0; atma <= (nAtomS1 - 3); atma = atma + 3){ i = i + 1; cek = 0;
        for(en = 1; en <= i; en++){
            if(atma != ok[en][1] && atma != ok[en][2] && atma != ok[en][3] && atma != ok[en][4] && atma != ok[en][5] && atma != ok[en][6])
            {
              cek = cek + 1;
            }
        }

    if(cek == i){

    for(atmb = 0; atmb <= (nAtomS1 - 3); atmb = atmb + 3)
    {
        if(atmb != atma){
        
        bonded = 0;
        for(crt = 0; crt < s1s1hbdn; crt++)
        {

           dist = distance(atmS1, atma, atmb, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atma, atma, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        } 
        if(bonded != 0){

    for(atmc = 0; atmc <= (nAtomS1 - 3); atmc = atmc + 3) 
    {
      if(atmc != atma && atmc != atmb){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmb, atmc, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){
        
        distOaOb = distanceOxOy(atmS1, atma, atmb, opos);
        distOaOc = distanceOxOy(atmS1, atma, atmc, opos);
        distObOc = distanceOxOy(atmS1, atmb, atmc, opos);
        angleOaObOc = (acos((pow(distOaOb,2)+pow(distObOc,2)-pow(distOaOc,2))/(2*distOaOb*distObOc))/PI)*180;
        if(angleOaObOc < 70 && angleOaObOc > 50){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmc, atma, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atma, atma, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmc, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmc, atmc, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmc, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){


    for(atmd = 0; atmd <= (nAtomS1 - 3); atmd = atmd + 3)
    {
      if(atmd != atma && atmd != atmb && atmd != atmc){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmc, atmd, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmd, atmd, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){

        distOcOd = distanceOxOy(atmS1, atmc, atmd, opos);
        distObOd = distanceOxOy(atmS1, atmb, atmd, opos);
        angleObOcOd = (acos((pow(distObOc,2)+pow(distOcOd,2)-pow(distObOd,2))/(2*distObOc*distOcOd))/PI)*180;

        if(angleObOcOd < 100 && angleObOcOd > 80){


        for(atme = 0; atme <= (nAtomS1-3); atme = atme + 3)
        {
          if(atme != atma && atme != atmb && atme != atmc && atme != atmd){

          bonded = 0;
          for(crt = 0; crt < s1s1hbdn; crt++)
          {

              dist = distance(atmS1, atmd, atme, s1a, s1b, crt);

              if(dist < s1as1bBD[crt])
              {
                 if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                 {
                    dist2 = distance(atmS1, atme, atme, s1s1v4, s1s1v5, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                 {
                    dist2 = distance(atmS1, atmd, atmd, s1s1v3, s1s1v4, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                 {
                    bonded = bonded + 1;
                 }

             }

          }
          if(bonded != 0){
 
          distOcOe = distanceOxOy(atmS1, atmc, atme, opos);
          distOdOe = distanceOxOy(atmS1, atmd, atme, opos);
          angleOcOdOe = (acos((pow(distOcOd,2)+pow(distOdOe,2)-pow(distOcOe,2))/(2*distOcOd*distOdOe))/PI)*180;

          if(angleOcOdOe < 100 && angleOcOdOe > 80){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atme, atma, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atma, atma, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atme, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atme, atme, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atme, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){


          for(atmf = 0; atmf <= (nAtomS1 - 3); atmf = atmf + 3)
          {
            if(atmf != atma && atmf != atmb && atmf != atmc && atmf != atmd && atmf != atme){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atme, atmf, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atme, atmf, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atme, atme, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atme, atmf, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){
     
            distOdOf = distanceOxOy(atmS1, atmd, atmf, opos);
            distOeOf = distanceOxOy(atmS1, atme, atmf, opos);
            angleOdOeOf = (acos((pow(distOdOe,2)+pow(distOeOf,2)-pow(distOdOf,2))/(2*distOdOe*distOeOf))/PI)*180;

            if(angleOdOeOf < 70 && angleOdOeOf > 50){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmf, atmb, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atmb, atmb, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmf, atmb, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmf, atmb, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){

            distOeOb = distanceOxOy(atmS1, atme, atmb, opos);
            distOfOb = distanceOxOy(atmS1, atmf, atmb, opos);
            angleOeOfOb = (acos((pow(distOeOf,2)+pow(distOfOb,2)-pow(distOeOb,2))/(2*distOeOf*distOfOb))/PI)*180;

            if(angleOeOfOb < 100 && angleOeOfOb > 80){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmf, atmd, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atmd, atmd, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmf, atmd, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmf, atmd, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){     

               

            ok[i][1] = atma; ok[i][2] = atmb; ok[i][3] = atmc; ok[i][4] = atmd; ok[i][5] = atme; ok[i][6] = atmf;

            v6[0] = atma; v6[1] = atmb; v6[2] = atmc; v6[3] = atmd; v6[4] = atme; v6[5] = atmf; 

            v2[cntr] = atma; v2[cntr+1] = atmb; v2[cntr+2] = atmc; v2[cntr+3] = atmd; v2[cntr+4] = atme; v2[cntr+5] = atmf;

            cntr = cntr + 6;

            e = 0;

            if(perm6(v6, 6, 0, v2, &e) == 1){ fprintf(outputfHPrism,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);

             
            if(hiso == 1)
            { 
              if(check_isolated_hexamer(atmS1,nAtomS1,s1s1hbdn,s1a,s1b,s1as1bBD,s1s1hban,s1s1v1,s1s1v2,s1s1v3,s1s1v4,s1s1v5,s1s1v6,atma,atmb,atmc,atmd,atme,atmf) == 1)
         
                 fprintf(outputfHPrismIso,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);
            }

           }
           }
           }
           }
           }
           }
          }
          }
          }
          }
          }
         }
         }
         }
         }
        }
        }
        }
       }
      }
     }
    }
   }
  }
 }
 }
 }
if(hiso == 1) fclose(outputfHPrismIso);
}


// Water Hexamer Cage search

void ChemNetworkNew::hcagesearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                 int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                 FILE *outputfHCage, char *foutputHCageIso)
{

   FILE *outputfHCageIso;
   int nodei, i, j, crt;
   double dist, dist2, hyptns, ang;
   int atma, atmb, atmc, atmd, atme, atmf, bonded; 
   int ok[MAX][7], cntr, in, v2[100000], v6[6], en, cek, e;
   double distOaOb, distOaOc, distObOc, distOcOd, distObOd;
   double angleOaObOc, angleObOcOd;

   nodei = nd1;

   if(hiso == 1) outputfHCageIso = fopen(foutputHCageIso,"w"); 
  
   for(i = 0; i < MAX; i++)
       for(j = 0; j < 7; j++)
           ok[i][j] = 0;
    i = 0;          
    cntr = 0; for(in = 0; in < 100000; in++) v2[in] = 0; for(in = 0; in < 6; in++) v6[in] = 0; 

    for(atma = 0; atma <= (nAtomS1 - 3); atma = atma + 3){ i = i + 1; cek = 0;
        for(en = 1; en <= i; en++){
            if(atma != ok[en][1] && atma != ok[en][2] && atma != ok[en][3] && atma != ok[en][4] && atma != ok[en][5] && atma != ok[en][6])
            {
              cek = cek + 1;
            }
        }

    if(cek == i){

    for(atmb = 0; atmb <= (nAtomS1 - 3); atmb = atmb + 3)
    {
        if(atmb != atma){
        
        bonded = 0;
        for(crt = 0; crt < s1s1hbdn; crt++)
        {

           dist = distance(atmS1, atma, atmb, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atma, atma, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        } 
        if(bonded != 0){

    for(atmc = 0; atmc <= (nAtomS1 - 3); atmc = atmc + 3) 
    {
      if(atmc != atma && atmc != atmb){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmb, atmc, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){
        
        distOaOb = distanceOxOy(atmS1, atma, atmb, opos);
        distOaOc = distanceOxOy(atmS1, atma, atmc, opos);
        distObOc = distanceOxOy(atmS1, atmb, atmc, opos);
        angleOaObOc = (acos((pow(distOaOb,2)+pow(distObOc,2)-pow(distOaOc,2))/(2*distOaOb*distObOc))/PI)*180;
        if(angleOaObOc < 120 && angleOaObOc > 60){

    for(atmd = 0; atmd <= (nAtomS1 - 3); atmd = atmd + 3)
    {
      if(atmd != atma && atmd != atmb && atmd != atmc){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmc, atmd, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmd, atmd, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){

        distOcOd = distanceOxOy(atmS1, atmc, atmd, opos);
        distObOd = distanceOxOy(atmS1, atmb, atmd, opos);
        angleObOcOd = (acos((pow(distObOc,2)+pow(distOcOd,2)-pow(distObOd,2))/(2*distObOc*distOcOd))/PI)*180;

        if(angleObOcOd < 120 && angleObOcOd > 60){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmd, atma, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atma, atma, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmd, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmd, atmd, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmd, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){

            if(fabs(dihedral(atmS1,atma,atmb,atmc,atmd,opos)) < 90 && fabs(dihedral(atmS1,atma,atmb,atmc,atmd,opos)) > 40){   

        for(atme = 0; atme <= (nAtomS1-3); atme = atme + 3)
        {
          if(atme != atma && atme != atmb && atme != atmc && atme != atmd){

          bonded = 0;
          for(crt = 0; crt < s1s1hbdn; crt++)
          {

              dist = distance(atmS1, atmc, atme, s1a, s1b, crt);

              if(dist < s1as1bBD[crt])
              {
                 if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                 {
                    dist2 = distance(atmS1, atme, atme, s1s1v4, s1s1v5, crt);
                    hyptns = distance(atmS1, atmc, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                 {
                    dist2 = distance(atmS1, atmc, atmc, s1s1v3, s1s1v4, crt);
                    hyptns = distance(atmS1, atmc, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                 {
                    bonded = bonded + 1;
                 }

             }

          }
          if(bonded != 0){
 
        
            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atme, atma, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atma, atma, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atme, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atme, atme, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atme, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){
 

          for(atmf = 0; atmf <= (nAtomS1 - 3); atmf = atmf + 3)
          {
            if(atmf != atma && atmf != atmb && atmf != atmc && atmf != atmd && atmf != atme){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmd, atmf, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmd, atmf, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmd, atmd, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmd, atmf, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){
     


            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmf, atmb, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atmb, atmb, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmf, atmb, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmf, atmb, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){


            if(fabs(dihedral(atmS1,atme,atma,atmc,atmf, opos)) <= 180 && fabs(dihedral(atmS1,atme,atma,atmc,atmf, opos)) > 150 && fabs(dihedral(atmS1,atmf,atmb,atmd,atme, opos)) <= 180 && fabs(dihedral(atmS1,atmf,atmb,atmd,atme, opos)) > 150){      

            ok[i][1] = atma; ok[i][2] = atmb; ok[i][3] = atmc; ok[i][4] = atmd; ok[i][5] = atme; ok[i][6] = atmf;

            v6[0] = atma; v6[1] = atmb; v6[2] = atmc; v6[3] = atmd; v6[4] = atme; v6[5] = atmf; 

            v2[cntr] = atma; v2[cntr+1] = atmb; v2[cntr+2] = atmc; v2[cntr+3] = atmd; v2[cntr+4] = atme; v2[cntr+5] = atmf;

            cntr = cntr + 6;

            e = 0;

            if(perm6(v6, 6, 0, v2, &e) == 1){ fprintf(outputfHCage,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);

             
            if(hiso == 1)
            { 
              if(check_isolated_hexamer(atmS1,nAtomS1,s1s1hbdn,s1a,s1b,s1as1bBD,s1s1hban,s1s1v1,s1s1v2,s1s1v3,s1s1v4,s1s1v5,s1s1v6,atma,atmb,atmc,atmd,atme,atmf) == 1)
         
                 fprintf(outputfHCageIso,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);
            }

           }
           }
           }
           }
           }
           }
          }
          }
          }
          }
          }
         }
         }
         }
         }
        }
        }
        }
       }
      }
     }
    }
   }
  }
 }
 
if(hiso == 1) fclose(outputfHCageIso);
}


// Water Hexamer Bag search

void ChemNetworkNew::hbagsearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                FILE *outputfHBag, char *foutputHBagIso)
{

   FILE *outputfHBagIso;
   int nodei, i, j, crt;
   double dist, dist2, hyptns, ang;
   int atma, atmb, atmc, atmd, atme, atmf, bonded; 
   int ok[MAX][7], cntr, in, v2[100000], v6[6], en, cek, e;
   double distOaOb, distOaOc, distObOc, distOcOd, distObOd;
   double angleOaObOc, angleObOcOd;

   nodei = nd1;

   if(hiso == 1) outputfHBagIso = fopen(foutputHBagIso,"w"); 
  
   for(i = 0; i < MAX; i++)
       for(j = 0; j < 7; j++)
           ok[i][j] = 0;
    i = 0;          
    cntr = 0; for(in = 0; in < 100000; in++) v2[in] = 0; for(in = 0; in < 6; in++) v6[in] = 0; 

    for(atma = 0; atma <= (nAtomS1 - 3); atma = atma + 3){ i = i + 1; cek = 0;
        for(en = 1; en <= i; en++){
            if(atma != ok[en][1] && atma != ok[en][2] && atma != ok[en][3] && atma != ok[en][4] && atma != ok[en][5] && atma != ok[en][6])
            {
              cek = cek + 1;
            }
        }

    if(cek == i){

    for(atmb = 0; atmb <= (nAtomS1 - 3); atmb = atmb + 3)
    {
        if(atmb != atma){
        
        bonded = 0;
        for(crt = 0; crt < s1s1hbdn; crt++)
        {

           dist = distance(atmS1, atma, atmb, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atma, atma, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        } 
        if(bonded != 0){

    for(atmc = 0; atmc <= (nAtomS1 - 3); atmc = atmc + 3) 
    {
      if(atmc != atma && atmc != atmb){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmb, atmc, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){
        
        distOaOb = distanceOxOy(atmS1, atma, atmb, opos);
        distOaOc = distanceOxOy(atmS1, atma, atmc, opos);
        distObOc = distanceOxOy(atmS1, atmb, atmc, opos);
        angleOaObOc = (acos((pow(distOaOb,2)+pow(distObOc,2)-pow(distOaOc,2))/(2*distOaOb*distObOc))/PI)*180;
        if(angleOaObOc < 140 && angleOaObOc > 80){

    for(atmd = 0; atmd <= (nAtomS1 - 3); atmd = atmd + 3)
    {
      if(atmd != atma && atmd != atmb && atmd != atmc){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmc, atmd, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmd, atmd, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){

        distOcOd = distanceOxOy(atmS1, atmc, atmd, opos);
        distObOd = distanceOxOy(atmS1, atmb, atmd, opos);
        angleObOcOd = (acos((pow(distObOc,2)+pow(distOcOd,2)-pow(distObOd,2))/(2*distObOc*distOcOd))/PI)*180;

        if(angleObOcOd < 140 && angleObOcOd > 80){

        if(fabs(dihedral(atmS1,atma,atmb,atmc,atmd,opos)) <= 25 && fabs(dihedral(atmS1,atma,atmb,atmc,atmd,opos)) >= 0){

        for(atme = 0; atme <= (nAtomS1-3); atme = atme + 3)
        {
          if(atme != atma && atme != atmb && atme != atmc && atme != atmd){

          bonded = 0;
          for(crt = 0; crt < s1s1hbdn; crt++)
          {

              dist = distance(atmS1, atmd, atme, s1a, s1b, crt);

              if(dist < s1as1bBD[crt])
              {
                 if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                 {
                    dist2 = distance(atmS1, atme, atme, s1s1v4, s1s1v5, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                 {
                    dist2 = distance(atmS1, atmd, atmd, s1s1v3, s1s1v4, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                 {
                    bonded = bonded + 1;
                 }

             }

          }
          if(bonded != 0){
 

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atme, atma, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atma, atma, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atme, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atme, atme, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atme, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){

          for(atmf = 0; atmf <= (nAtomS1 - 3); atmf = atmf + 3)
          {
            if(atmf != atma && atmf != atmb && atmf != atmc && atmf != atmd && atmf != atme){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmd, atmf, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmd, atmf, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmd, atmd, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmd, atmf, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){
     


            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmf, atma, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atma, atma, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmf, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmf, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){

           
            if(fabs(dihedral(atmS1,atme,atma,atmd,atmf,opos)) <= 180 && fabs(dihedral(atmS1,atme,atma,atmd,atmf,opos)) >= 60 && fabs(dihedral(atmS1,atmf,atmd,atma,atme,opos)) <= 180 && fabs(dihedral(atmS1,atmf,atmd,atma,atme,opos)) >= 60){ 

              if(prismbook_check(atmS1, nd1, opos, s1s1hbdn, s1a, s1b, s1as1bBD, 
                                 s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6,
                                 atma, atmb, atmc, atmd, atme, atmf) != 1){  

            ok[i][1] = atma; ok[i][2] = atmb; ok[i][3] = atmc; ok[i][4] = atmd; ok[i][5] = atme; ok[i][6] = atmf;

            v6[0] = atma; v6[1] = atmb; v6[2] = atmc; v6[3] = atmd; v6[4] = atme; v6[5] = atmf; 

            v2[cntr] = atma; v2[cntr+1] = atmb; v2[cntr+2] = atmc; v2[cntr+3] = atmd; v2[cntr+4] = atme; v2[cntr+5] = atmf;

            cntr = cntr + 6;

            e = 0;

            if(perm6(v6, 6, 0, v2, &e) == 1){ fprintf(outputfHBag,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);

             
            if(hiso == 1)
            { 
              if(check_isolated_hexamer(atmS1,nAtomS1,s1s1hbdn,s1a,s1b,s1as1bBD,s1s1hban,s1s1v1,s1s1v2,s1s1v3,s1s1v4,s1s1v5,s1s1v6,atma,atmb,atmc,atmd,atme,atmf) == 1)
         
                 fprintf(outputfHBagIso,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);
            }
           
           }
           }
           }
           }
           }
           }
           }
          }
          }
          }
          }
          }
         }
         }
         }
         }
        }
        }
        }
       }
      }
     }
    }
   }
  }
 
if(hiso == 1) fclose(outputfHBagIso);
}


// Water Hexamer Boat search

void ChemNetworkNew::hboatsearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                 int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                 FILE *outputfHBoat, char *foutputHBoatIso)
{

   FILE *outputfHBoatIso;
   int nodei, i, j, crt;
   double dist, dist2, hyptns, ang;
   int atma, atmb, atmc, atmd, atme, atmf, bonded; 
   int ok[MAX][7], cntr, in, v2[100000], v6[6], en, cek, e;
   double distOaOb, distOaOc, distObOc, distOcOd, distObOd, distOdOe, distOeOf, distOaOe, distOaOd, distOeOb, distOfOd, distOfOb, distOfOc;
   double angleOaObOc, angleOcObOd, angleOaOdOe, angleObOdOe, angleOfOeOd, angleOfOcOb;

   nodei = nd1;

   if(hiso == 1) outputfHBoatIso = fopen(foutputHBoatIso,"w"); 
  
   for(i = 0; i < MAX; i++)
       for(j = 0; j < 7; j++)
           ok[i][j] = 0;
    i = 0;          
    cntr = 0; for(in = 0; in < 100000; in++) v2[in] = 0; for(in = 0; in < 6; in++) v6[in] = 0; 


    for(atma = 0; atma <= (nAtomS1 - 3); atma = atma + 3){ i = i + 1; cek = 0;
        for(en = 1; en <= i; en++){
            if(atma != ok[en][1] && atma != ok[en][2] && atma != ok[en][3] && atma != ok[en][4] && atma != ok[en][5] && atma != ok[en][6])
            {
              cek = cek + 1;
            }
        }

    if(cek == i){
                   
    for(atmb = 0; atmb <= (nAtomS1 - 3); atmb = atmb + 3)
    {
        if(atmb != atma){
        
        bonded = 0;
        for(crt = 0; crt < s1s1hbdn; crt++)
        {

           dist = distance(atmS1, atma, atmb, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atma, atma, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        } 
        if(bonded != 0){

    for(atmc = 0; atmc <= (nAtomS1 - 3); atmc = atmc + 3) 
    {
      if(atmc != atma && atmc != atmb){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmb, atmc, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){
        
        distOaOb = distanceOxOy(atmS1, atma, atmb, opos);
        distOaOc = distanceOxOy(atmS1, atma, atmc, opos);
        distObOc = distanceOxOy(atmS1, atmb, atmc, opos);
        angleOaObOc = (acos((pow(distOaOb,2)+pow(distObOc,2)-pow(distOaOc,2))/(2*distOaOb*distObOc))/PI)*180;
        if(angleOaObOc < 130 && angleOaObOc > 100){

    for(atmd = 0; atmd <= (nAtomS1 - 3); atmd = atmd + 3)
    {
      if(atmd != atma && atmd != atmb && atmd != atmc){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atma, atmd, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmd, atmd, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atma, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atma, atma, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atma, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){

        distOcOd = distanceOxOy(atmS1, atmc, atmd, opos);
        distObOd = distanceOxOy(atmS1, atmb, atmd, opos);
        angleOcObOd = (acos((pow(distObOc,2)+pow(distObOd,2)-pow(distOcOd,2))/(2*distObOc*distObOd))/PI)*180;
       
        if(angleOcObOd < 100 && angleOcObOd > 80){

        if((dihedral2(atmS1,atmc,atmd,atmb,atma,opos) <= 90 && dihedral2(atmS1,atmc,atmd,atmb,atma,opos) >= 10) || (dihedral2(atmS1,atmc,atmd,atmb,atma,opos) <= -10 && dihedral2(atmS1,atmc,atmd,atmb,atma,opos)>= -90) ){
   

        for(atme = 0; atme <= (nAtomS1-3); atme = atme + 3)
        {
          if(atme != atma && atme != atmb && atme != atmc && atme != atmd){

          bonded = 0;
          for(crt = 0; crt < s1s1hbdn; crt++)
          {

              dist = distance(atmS1, atmd, atme, s1a, s1b, crt);

              if(dist < s1as1bBD[crt])
              {
                 if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                 {
                    dist2 = distance(atmS1, atme, atme, s1s1v4, s1s1v5, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                 {
                    dist2 = distance(atmS1, atmd, atmd, s1s1v3, s1s1v4, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                 {
                    bonded = bonded + 1;
                 }

             }

          }
          if(bonded != 0){
 
          distOaOe = distanceOxOy(atmS1, atma, atme, opos);
          distOdOe = distanceOxOy(atmS1, atmd, atme, opos);
          distOaOd = distanceOxOy(atmS1, atma, atmd, opos);
          angleOaOdOe = (acos((pow(distOaOd,2)+pow(distOdOe,2)-pow(distOaOe,2))/(2*distOaOd*distOdOe))/PI)*180;

          if(angleOaOdOe < 130 && angleOaOdOe > 100){

          distObOd = distanceOxOy(atmS1, atmb, atmd, opos);
          distOeOb = distanceOxOy(atmS1, atme, atmb, opos);
          angleObOdOe = (acos((pow(distObOd,2)+pow(distOdOe,2)-pow(distOeOb,2))/(2*distObOd*distOdOe))/PI)*180;

          if(angleObOdOe <= 100 && angleObOdOe >= 80){

          if(fabs(dihedral(atmS1,atmb,atmc,atme,atmd,opos)) <= 15){

          for(atmf = 0; atmf <= (nAtomS1 - 3); atmf = atmf + 3)
          {
            if(atmf != atma && atmf != atmb && atmf != atmc && atmf != atmd && atmf != atme){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atme, atmf, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atme, atmf, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atme, atme, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atme, atmf, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){
     
            distOfOd = distanceOxOy(atmS1, atmf, atmd, opos);
            distOeOf = distanceOxOy(atmS1, atme, atmf, opos);
            angleOfOeOd = (acos((pow(distOeOf,2)+pow(distOdOe,2)-pow(distOfOd,2))/(2*distOeOf*distOdOe))/PI)*180;

            if(angleOfOeOd <= 130 && angleOfOeOd >= 100){ 

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmf, atmc, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atmc, atmc, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmf, atmc, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmf, atmc, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){

            distOfOb = distanceOxOy(atmS1, atmf, atmb, opos);
            distOfOc = distanceOxOy(atmS1, atmf, atmc, opos);
            angleOfOcOb = (acos((pow(distOfOc,2)+pow(distObOc,2)-pow(distOfOb,2))/(2*distOfOc*distObOc))/PI)*180;

            if(angleOfOcOb <= 130 && angleOfOcOb >= 100){

            if(dihedral2(atmS1,atmc,atmd,atmb,atma,opos) <= 90 && dihedral2(atmS1,atmc,atmd,atmb,atma,opos) >= 10 ){     

            if(dihedral2(atmS1,atmd,atmc,atme,atmf,opos) <= 90 && dihedral2(atmS1,atmd,atmc,atme,atmf,opos) >= 10 ){ 

               if(ring_check(atmS1, nd1, opos, s1s1hbdn, s1a, s1b, s1as1bBD, 
                             s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6,
                             atma, atmb, atmc, atmd, atme, atmf) != 1){    

            ok[i][1] = atma; ok[i][2] = atmb; ok[i][3] = atmc; ok[i][4] = atmd; ok[i][5] = atme; ok[i][6] = atmf;

            v6[0] = atma; v6[1] = atmb; v6[2] = atmc; v6[3] = atmd; v6[4] = atme; v6[5] = atmf; 

            v2[cntr] = atma; v2[cntr+1] = atmb; v2[cntr+2] = atmc; v2[cntr+3] = atmd; v2[cntr+4] = atme; v2[cntr+5] = atmf;

            cntr = cntr + 6;

            e = 0;

            if(perm6(v6, 6, 0, v2, &e) == 1){ fprintf(outputfHBoat,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);

             
            if(hiso == 1)
            { 
              if(check_isolated_hexamer(atmS1,nAtomS1,s1s1hbdn,s1a,s1b,s1as1bBD,s1s1hban,s1s1v1,s1s1v2,s1s1v3,s1s1v4,s1s1v5,s1s1v6,atma,atmb,atmc,atmd,atme,atmf) == 1)
         
                 fprintf(outputfHBoatIso,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);
            }

           }
           }
           }
           }

           else if(dihedral2(atmS1,atmc,atmd,atmb,atma,opos) <= -10 && dihedral2(atmS1,atmc,atmd,atmb,atma,opos) >= -90 ){

                if( dihedral2(atmS1,atmd,atmc,atme,atmf,opos) <= -10 && dihedral2(atmS1,atmd,atmc,atme,atmf,opos) >= -90 ){

                if(ring_check(atmS1, nd1, opos, s1s1hbdn, s1a, s1b, s1as1bBD,
                             s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6,
                             atma, atmb, atmc, atmd, atme, atmf) != 1){

            ok[i][1] = atma; ok[i][2] = atmb; ok[i][3] = atmc; ok[i][4] = atmd; ok[i][5] = atme; ok[i][6] = atmf;

            v6[0] = atma; v6[1] = atmb; v6[2] = atmc; v6[3] = atmd; v6[4] = atme; v6[5] = atmf;

            v2[cntr] = atma; v2[cntr+1] = atmb; v2[cntr+2] = atmc; v2[cntr+3] = atmd; v2[cntr+4] = atme; v2[cntr+5] = atmf;

            cntr = cntr + 6;

            e = 0;

            if(perm6(v6, 6, 0, v2, &e) == 1){ fprintf(outputfHBoat,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);


            if(hiso == 1)
            {
              if(check_isolated_hexamer(atmS1,nAtomS1,s1s1hbdn,s1a,s1b,s1as1bBD,s1s1hban,s1s1v1,s1s1v2,s1s1v3,s1s1v4,s1s1v5,s1s1v6,atma,atmb,atmc,atmd,atme,atmf) == 1)

                 fprintf(outputfHBoatIso,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);
            }

           }
           }
           }
           }
          }
          }
          }
          }
          }
         }
         }
         }
         }
        }
        }
        }
       }
       }
      }
      }
     }
     }
     }
     }
    }
   }
   }
  }
  }
 }
    
if(hiso == 1) fclose(outputfHBoatIso);
}

// Water Hexamer Chair search

void ChemNetworkNew::hchairsearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                  int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                  FILE *outputfHChair, char *foutputHChairIso)
{

   FILE *outputfHChairIso;
   int nodei, i, j, crt;
   double dist, dist2, hyptns, ang;
   int atma, atmb, atmc, atmd, atme, atmf, bonded; 
   int ok[MAX][7], cntr, in, v2[100000], v6[6], en, cek, e;
   double distOaOb, distOaOc, distObOc, distOcOd, distObOd, distOdOe, distOeOf, distOaOe, distOaOd, distOeOb, distOfOd, distOfOb, distOfOc;
   double angleOaObOc, angleOcObOd, angleOaOdOe, angleObOdOe, angleOfOeOd, angleOfOcOb;

   nodei = nd1;

   if(hiso == 1) outputfHChairIso = fopen(foutputHChairIso,"w"); 
  
   for(i = 0; i < MAX; i++)
       for(j = 0; j < 7; j++)
           ok[i][j] = 0;
    i = 0;          
    cntr = 0; for(in = 0; in < 100000; in++) v2[in] = 0; for(in = 0; in < 6; in++) v6[in] = 0; 


    for(atma = 0; atma <= (nAtomS1 - 3); atma = atma + 3){ i = i + 1; cek = 0;
        for(en = 1; en <= i; en++){
            if(atma != ok[en][1] && atma != ok[en][2] && atma != ok[en][3] && atma != ok[en][4] && atma != ok[en][5] && atma != ok[en][6])
            {
              cek = cek + 1;
            }
        }

    if(cek == i){
                   
    for(atmb = 0; atmb <= (nAtomS1 - 3); atmb = atmb + 3)
    {
        if(atmb != atma){
        
        bonded = 0;
        for(crt = 0; crt < s1s1hbdn; crt++)
        {

           dist = distance(atmS1, atma, atmb, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atma, atma, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        } 
        if(bonded != 0){

    for(atmc = 0; atmc <= (nAtomS1 - 3); atmc = atmc + 3) 
    {
      if(atmc != atma && atmc != atmb){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmb, atmc, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){
        
        distOaOb = distanceOxOy(atmS1, atma, atmb, opos);
        distOaOc = distanceOxOy(atmS1, atma, atmc, opos);
        distObOc = distanceOxOy(atmS1, atmb, atmc, opos);
        angleOaObOc = (acos((pow(distOaOb,2)+pow(distObOc,2)-pow(distOaOc,2))/(2*distOaOb*distObOc))/PI)*180;
        if(angleOaObOc < 130 && angleOaObOc > 100){

    for(atmd = 0; atmd <= (nAtomS1 - 3); atmd = atmd + 3)
    {
      if(atmd != atma && atmd != atmb && atmd != atmc){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atma, atmd, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmd, atmd, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atma, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atma, atma, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atma, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){

        distOcOd = distanceOxOy(atmS1, atmc, atmd, opos);
        distObOd = distanceOxOy(atmS1, atmb, atmd, opos);
        angleOcObOd = (acos((pow(distObOc,2)+pow(distObOd,2)-pow(distOcOd,2))/(2*distObOc*distObOd))/PI)*180;
       
        if(angleOcObOd < 100 && angleOcObOd > 80){

        if((dihedral2(atmS1,atmc,atmd,atmb,atma,opos) <= 90 && dihedral2(atmS1,atmc,atmd,atmb,atma,opos) >= 10) || (dihedral2(atmS1,atmc,atmd,atmb,atma,opos) <= -10 && dihedral2(atmS1,atmc,atmd,atmb,atma,opos)>= -90) ){
   

        for(atme = 0; atme <= (nAtomS1-3); atme = atme + 3)
        {
          if(atme != atma && atme != atmb && atme != atmc && atme != atmd){

          bonded = 0;
          for(crt = 0; crt < s1s1hbdn; crt++)
          {

              dist = distance(atmS1, atmd, atme, s1a, s1b, crt);

              if(dist < s1as1bBD[crt])
              {
                 if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                 {
                    dist2 = distance(atmS1, atme, atme, s1s1v4, s1s1v5, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                 {
                    dist2 = distance(atmS1, atmd, atmd, s1s1v3, s1s1v4, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                 {
                    bonded = bonded + 1;
                 }

             }

          }
          if(bonded != 0){
 
          distOaOe = distanceOxOy(atmS1, atma, atme, opos);
          distOdOe = distanceOxOy(atmS1, atmd, atme, opos);
          distOaOd = distanceOxOy(atmS1, atma, atmd, opos);
          angleOaOdOe = (acos((pow(distOaOd,2)+pow(distOdOe,2)-pow(distOaOe,2))/(2*distOaOd*distOdOe))/PI)*180;

          if(angleOaOdOe < 130 && angleOaOdOe > 100){

          distObOd = distanceOxOy(atmS1, atmb, atmd, opos);
          distOeOb = distanceOxOy(atmS1, atme, atmb, opos);
          angleObOdOe = (acos((pow(distObOd,2)+pow(distOdOe,2)-pow(distOeOb,2))/(2*distObOd*distOdOe))/PI)*180;

          if(angleObOdOe <= 100 && angleObOdOe >= 80){

          if(fabs(dihedral(atmS1,atmb,atmc,atme,atmd,opos)) <= 15){

          for(atmf = 0; atmf <= (nAtomS1 - 3); atmf = atmf + 3)
          {
            if(atmf != atma && atmf != atmb && atmf != atmc && atmf != atmd && atmf != atme){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atme, atmf, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atme, atmf, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atme, atme, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atme, atmf, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){
     
            distOfOd = distanceOxOy(atmS1, atmf, atmd, opos);
            distOeOf = distanceOxOy(atmS1, atme, atmf, opos);
            angleOfOeOd = (acos((pow(distOeOf,2)+pow(distOdOe,2)-pow(distOfOd,2))/(2*distOeOf*distOdOe))/PI)*180;

            if(angleOfOeOd <= 130 && angleOfOeOd >= 100){ 

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmf, atmc, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atmc, atmc, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmf, atmc, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmf, atmc, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){

            distOfOb = distanceOxOy(atmS1, atmf, atmb, opos);
            distOfOc = distanceOxOy(atmS1, atmf, atmc, opos);
            angleOfOcOb = (acos((pow(distOfOc,2)+pow(distObOc,2)-pow(distOfOb,2))/(2*distOfOc*distObOc))/PI)*180;

            if(angleOfOcOb <= 130 && angleOfOcOb >= 100){

            if(dihedral2(atmS1,atmc,atmd,atmb,atma,opos) <= 90 && dihedral2(atmS1,atmc,atmd,atmb,atma,opos) >= 10 ){     

            if(dihedral2(atmS1,atmd,atmc,atme,atmf,opos) <= -10 && dihedral2(atmS1,atmd,atmc,atme,atmf,opos) >= -90 ){ 

            if(ring_check(atmS1, nd1, opos, s1s1hbdn, s1a, s1b, s1as1bBD,
                          s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6,
                          atma, atmb, atmc, atmd, atme, atmf) != 1){     

            ok[i][1] = atma; ok[i][2] = atmb; ok[i][3] = atmc; ok[i][4] = atmd; ok[i][5] = atme; ok[i][6] = atmf;

            v6[0] = atma; v6[1] = atmb; v6[2] = atmc; v6[3] = atmd; v6[4] = atme; v6[5] = atmf; 

            v2[cntr] = atma; v2[cntr+1] = atmb; v2[cntr+2] = atmc; v2[cntr+3] = atmd; v2[cntr+4] = atme; v2[cntr+5] = atmf;

            cntr = cntr + 6;

            e = 0;

            if(perm6(v6, 6, 0, v2, &e) == 1){ fprintf(outputfHChair,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);

             
            if(hiso == 1)
            { 
              if(check_isolated_hexamer(atmS1,nAtomS1,s1s1hbdn,s1a,s1b,s1as1bBD,s1s1hban,s1s1v1,s1s1v2,s1s1v3,s1s1v4,s1s1v5,s1s1v6,atma,atmb,atmc,atmd,atme,atmf) == 1)
         
                 fprintf(outputfHChairIso,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);
            }

           }
           }
           }
           }

           else if(dihedral2(atmS1,atmc,atmd,atmb,atma,opos) <= -10 && dihedral2(atmS1,atmc,atmd,atmb,atma,opos) >= -90 ){

                if( dihedral2(atmS1,atmd,atmc,atme,atmf,opos) <= 90 && dihedral2(atmS1,atmd,atmc,atme,atmf,opos) >= 10 ){

                if(ring_check(atmS1, nd1, opos, s1s1hbdn, s1a, s1b, s1as1bBD,
                              s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6,
                              atma, atmb, atmc, atmd, atme, atmf) != 1){ 

            ok[i][1] = atma; ok[i][2] = atmb; ok[i][3] = atmc; ok[i][4] = atmd; ok[i][5] = atme; ok[i][6] = atmf;

            v6[0] = atma; v6[1] = atmb; v6[2] = atmc; v6[3] = atmd; v6[4] = atme; v6[5] = atmf;

            v2[cntr] = atma; v2[cntr+1] = atmb; v2[cntr+2] = atmc; v2[cntr+3] = atmd; v2[cntr+4] = atme; v2[cntr+5] = atmf;

            cntr = cntr + 6;

            e = 0;

            if(perm6(v6, 6, 0, v2, &e) == 1){ fprintf(outputfHChair,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);


            if(hiso == 1)
            {
              if(check_isolated_hexamer(atmS1,nAtomS1,s1s1hbdn,s1a,s1b,s1as1bBD,s1s1hban,s1s1v1,s1s1v2,s1s1v3,s1s1v4,s1s1v5,s1s1v6,atma,atmb,atmc,atmd,atme,atmf) == 1)

                 fprintf(outputfHChairIso,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);
            }

           }
           }
           }
           }
          }
          }
          }
          }
          }
         }
         }
         }
         }
        }
        }
        }
       }
       }
      }
      }
     }
     }
     }
     }
    }
   }
   }
  }
  }
 }
    
if(hiso == 1) fclose(outputfHChairIso);
}


// Water Hexamer PrismBook search

void ChemNetworkNew::hprismbooksearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                      int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                      FILE *outputfHPrismBook, char *foutputHPrismBookIso)
{

   FILE *outputfHPrismBookIso;
   int nodei, i, j, crt;
   double dist, dist2, hyptns, ang;
   int atma, atmb, atmc, atmd, atme, atmf, bonded; 
   int ok[MAX][7], cntr, in, v2[100000], v6[6], en, cek, e;
   double distOaOb, distOaOc, distObOc, distOcOd, distObOd, distOcOe, distOdOe, distOdOf, distOeOf;
   double angleOaObOc, angleObOcOd, angleOcOdOe, angleOeOdOf;

   nodei = nd1;

   if(hiso == 1) outputfHPrismBookIso = fopen(foutputHPrismBookIso,"w"); 
  
   for(i = 0; i < MAX; i++)
       for(j = 0; j < 7; j++)
           ok[i][j] = 0;
    i = 0;          
    cntr = 0; for(in = 0; in < 100000; in++) v2[in] = 0; for(in = 0; in < 6; in++) v6[in] = 0; 

    for(atma = 0; atma <= (nAtomS1 - 3); atma = atma + 3){ i = i + 1; cek = 0;
        for(en = 1; en <= i; en++){
            if(atma != ok[en][1] && atma != ok[en][2] && atma != ok[en][3] && atma != ok[en][4] && atma != ok[en][5] && atma != ok[en][6])
            {
              cek = cek + 1;
            }
        }

    if(cek == i){

    for(atmb = 0; atmb <= (nAtomS1 - 3); atmb = atmb + 3)
    {
        if(atmb != atma){
        
        bonded = 0;
        for(crt = 0; crt < s1s1hbdn; crt++)
        {

           dist = distance(atmS1, atma, atmb, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atma, atma, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        } 
        if(bonded != 0){

    for(atmc = 0; atmc <= (nAtomS1 - 3); atmc = atmc + 3) 
    {
      if(atmc != atma && atmc != atmb){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmb, atmc, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){
        
        distOaOb = distanceOxOy(atmS1, atma, atmb, opos);
        distOaOc = distanceOxOy(atmS1, atma, atmc, opos);
        distObOc = distanceOxOy(atmS1, atmb, atmc, opos);
        angleOaObOc = (acos((pow(distOaOb,2)+pow(distObOc,2)-pow(distOaOc,2))/(2*distOaOb*distObOc))/PI)*180;
        if(angleOaObOc < 70 && angleOaObOc > 50){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmc, atma, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atma, atma, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmc, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmc, atmc, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmc, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){


    for(atmd = 0; atmd <= (nAtomS1 - 3); atmd = atmd + 3)
    {
      if(atmd != atma && atmd != atmb && atmd != atmc){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmc, atmd, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmd, atmd, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){

        distOcOd = distanceOxOy(atmS1, atmc, atmd, opos);
        distObOd = distanceOxOy(atmS1, atmb, atmd, opos);
        angleObOcOd = (acos((pow(distObOc,2)+pow(distOcOd,2)-pow(distObOd,2))/(2*distObOc*distOcOd))/PI)*180;

        if(angleObOcOd < 100 && angleObOcOd > 80){


        for(atme = 0; atme <= (nAtomS1-3); atme = atme + 3)
        {
          if(atme != atma && atme != atmb && atme != atmc && atme != atmd){

          bonded = 0;
          for(crt = 0; crt < s1s1hbdn; crt++)
          {

              dist = distance(atmS1, atmd, atme, s1a, s1b, crt);

              if(dist < s1as1bBD[crt])
              {
                 if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                 {
                    dist2 = distance(atmS1, atme, atme, s1s1v4, s1s1v5, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                 {
                    dist2 = distance(atmS1, atmd, atmd, s1s1v3, s1s1v4, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                 {
                    bonded = bonded + 1;
                 }

             }

          }
          if(bonded != 0){
 
          distOcOe = distanceOxOy(atmS1, atmc, atme, opos);
          distOdOe = distanceOxOy(atmS1, atmd, atme, opos);
          angleOcOdOe = (acos((pow(distOcOd,2)+pow(distOdOe,2)-pow(distOcOe,2))/(2*distOcOd*distOdOe))/PI)*180;

          if(angleOcOdOe < 100 && angleOcOdOe > 80){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atme, atma, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atma, atma, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atme, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atme, atme, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atme, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){


          for(atmf = 0; atmf <= (nAtomS1 - 3); atmf = atmf + 3)
          {
            if(atmf != atma && atmf != atmb && atmf != atmc && atmf != atmd && atmf != atme){

            distOdOf = distanceOxOy(atmS1, atmd, atmf, opos);
            distOeOf = distanceOxOy(atmS1, atme, atmf, opos);
            angleOeOdOf = (acos((pow(distOdOe,2)+pow(distOdOf,2)-pow(distOeOf,2))/(2*distOdOe*distOdOf))/PI)*180;

            if(angleOeOdOf > 70){ 

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmf, atmb, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atmb, atmb, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmf, atmb, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmf, atmb, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){
     
          

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmf, atmd, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atmd, atmd, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmf, atmd, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmf, atmd, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){


            ok[i][1] = atma; ok[i][2] = atmb; ok[i][3] = atmc; ok[i][4] = atmd; ok[i][5] = atme; ok[i][6] = atmf;

            v6[0] = atma; v6[1] = atmb; v6[2] = atmc; v6[3] = atmd; v6[4] = atme; v6[5] = atmf; 

            v2[cntr] = atma; v2[cntr+1] = atmb; v2[cntr+2] = atmc; v2[cntr+3] = atmd; v2[cntr+4] = atme; v2[cntr+5] = atmf;

            cntr = cntr + 6;

            e = 0;

            if(perm6(v6, 6, 0, v2, &e) == 1){ fprintf(outputfHPrismBook,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);

             
            if(hiso == 1)
            { 
              if(check_isolated_hexamer(atmS1,nAtomS1,s1s1hbdn,s1a,s1b,s1as1bBD,s1s1hban,s1s1v1,s1s1v2,s1s1v3,s1s1v4,s1s1v5,s1s1v6,atma,atmb,atmc,atmd,atme,atmf) == 1)
         
                 fprintf(outputfHPrismBookIso,"%d %d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3,1+atmf/3);
            }

           }
           }
           }
           }
           }
           }
          }
          }
          }
          }
          }
         }
         }
         }
         }
        }
        }
        }
       }
      }
     }
    }
   }
  }
 }

if(hiso == 1) fclose(outputfHPrismBookIso);
}


// Water Cyclic Pentamer search

void ChemNetworkNew::hpentamersearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                     int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                     FILE *outputfPentamer, char *foutputPentamerIso)
{

   FILE *outputfPentamerIso;
   int nodei, i, j, crt;
   double dist, dist2, hyptns, ang;
   int atma, atmb, atmc, atmd, atme, bonded; 
   int ok[MAX][6], cntr, in, v2[100000], v5[5], en, cek, e;
   double distOaOb, distOaOc, distObOc, distOcOd, distObOd, distOcOe, distOdOe, distOeOa;
   double angleOaObOc, angleObOcOd, angleOcOdOe, distOdOa, angleOdOeOa;

   nodei = nd1;

   if(hiso == 1) outputfPentamerIso = fopen(foutputPentamerIso,"w"); 
  
   for(i = 0; i < MAX; i++)
       for(j = 0; j < 6; j++)
           ok[i][j] = 0;
    i = 0;          
    cntr = 0; for(in = 0; in < 100000; in++) v2[in] = 0; for(in = 0; in < 5; in++) v5[in] = 0; 


    for(atma = 0; atma <= (nAtomS1 - 3); atma = atma + 3){ i = i + 1; cek = 0;
        for(en = 1; en <= i; en++){
            if(atma != ok[en][1] && atma != ok[en][2] && atma != ok[en][3] && atma != ok[en][4] && atma != ok[en][5])
            {
              cek = cek + 1;
            }
        }

    if(cek == i){

    for(atmb = 0; atmb <= (nAtomS1 - 3); atmb = atmb + 3)
    {
        if(atmb != atma){
        
        bonded = 0;
        for(crt = 0; crt < s1s1hbdn; crt++)
        {

           dist = distance(atmS1, atma, atmb, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atma, atma, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        } 
        if(bonded != 0){

    for(atmc = 0; atmc <= (nAtomS1 - 3); atmc = atmc + 3) 
    {
      if(atmc != atma && atmc != atmb){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmb, atmc, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){
        
        distOaOb = distanceOxOy(atmS1, atma, atmb, opos);
        distOaOc = distanceOxOy(atmS1, atma, atmc, opos);
        distObOc = distanceOxOy(atmS1, atmb, atmc, opos);
        angleOaObOc = (acos((pow(distOaOb,2)+pow(distObOc,2)-pow(distOaOc,2))/(2*distOaOb*distObOc))/PI)*180;
        if(angleOaObOc < 118 && angleOaObOc > 98){

    for(atmd = 0; atmd <= (nAtomS1 - 3); atmd = atmd + 3)
    {
      if(atmd != atma && atmd != atmb && atmd != atmc){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmc, atmd, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmd, atmd, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){

        distOcOd = distanceOxOy(atmS1, atmc, atmd, opos);
        distObOd = distanceOxOy(atmS1, atmb, atmd, opos);
        angleObOcOd = (acos((pow(distObOc,2)+pow(distOcOd,2)-pow(distObOd,2))/(2*distObOc*distOcOd))/PI)*180;

        if(angleObOcOd < 118 && angleObOcOd > 98){

        if(fabs(dihedral(atmS1,atma,atmb,atmc,atmd,opos)) <= 20){

        for(atme = 0; atme <= (nAtomS1-3); atme = atme + 3)
        {
          if(atme != atma && atme != atmb && atme != atmc && atme != atmd){

          bonded = 0;
          for(crt = 0; crt < s1s1hbdn; crt++)
          {

              dist = distance(atmS1, atmd, atme, s1a, s1b, crt);

              if(dist < s1as1bBD[crt])
              {
                 if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                 {
                    dist2 = distance(atmS1, atme, atme, s1s1v4, s1s1v5, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                 {
                    dist2 = distance(atmS1, atmd, atmd, s1s1v3, s1s1v4, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                 {
                    bonded = bonded + 1;
                 }

             }

          }
          if(bonded != 0){

          distOcOe = distanceOxOy(atmS1, atmc, atme, opos);
          distOdOe = distanceOxOy(atmS1, atmd, atme, opos);
          angleOcOdOe = (acos((pow(distOcOd,2)+pow(distOdOe,2)-pow(distOcOe,2))/(2*distOcOd*distOdOe))/PI)*180;

          if(angleOcOdOe < 118 && angleOcOdOe > 98){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atme, atma, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atma, atma, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atme, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atme, atme, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atme, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){

           distOeOa = distanceOxOy(atmS1, atme, atma, opos);
           distOdOa = distanceOxOy(atmS1, atmd, atma, opos);
           angleOdOeOa = (acos((pow(distOdOe,2)+pow(distOeOa,2)-pow(distOdOa,2))/(2*distOdOe*distOeOa))/PI)*180;

           if(angleOdOeOa < 118 && angleOdOeOa > 98){

           if(fabs(dihedral(atmS1,atma,atme,atmd,atmc,opos)) <= 20){

            ok[i][1] = atma; ok[i][2] = atmb; ok[i][3] = atmc; ok[i][4] = atmd; ok[i][5] = atme;

            v5[0] = atma; v5[1] = atmb; v5[2] = atmc; v5[3] = atmd; v5[4] = atme;

            v2[cntr] = atma; v2[cntr+1] = atmb; v2[cntr+2] = atmc; v2[cntr+3] = atmd; v2[cntr+4] = atme;

            cntr = cntr + 5;

            e = 0;

            if(perm5(v5, 5, 0, v2, &e) == 1){ fprintf(outputfPentamer,"%d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3);

             
            if(hiso == 1)
            { 
              if(check_isolated_pentamer(atmS1,nAtomS1,s1s1hbdn,s1a,s1b,s1as1bBD,s1s1hban,s1s1v1,s1s1v2,s1s1v3,s1s1v4,s1s1v5,s1s1v6,atma,atmb,atmc,atmd,atme) == 1)
         
                 fprintf(outputfPentamerIso,"%d %d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3,1+atme/3);
            }
           }
          }
         }
        }
        }
        }
       }
      }
     }
     }
     }
     }
    }
    }
    }
   }
   }
   }
   }
  }
  }
  }

if(hiso == 1) fclose(outputfPentamerIso);
}


// Water Cyclic Tetramer search

void ChemNetworkNew::htetramersearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                     int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                     FILE *outputfTetramer, char *foutputTetramerIso)
{

   FILE *outputfTetramerIso;
   int nodei, i, j, crt;
   double dist, dist2, hyptns, ang;
   int atma, atmb, atmc, atmd, bonded; 
   int ok[MAX][5], cntr, in, v2[100000], v4[4], en, cek, e;
   double distOaOb, distOaOc, distObOc, distOcOd, distObOd;
   double angleOaObOc, angleObOcOd, distOdOa, distOcOa, angleOcOdOa, distOdOb, angleOdOaOb;

   nodei = nd1;

   if(hiso == 1) outputfTetramerIso = fopen(foutputTetramerIso,"w"); 
  
   for(i = 0; i < MAX; i++)
       for(j = 0; j < 5; j++)
           ok[i][j] = 0;
    i = 0;          
    cntr = 0; for(in = 0; in < 100000; in++) v2[in] = 0; for(in = 0; in < 4; in++) v4[in] = 0; 


    for(atma = 0; atma <= (nAtomS1 - 3); atma = atma + 3){ i = i + 1; cek = 0;
        for(en = 1; en <= i; en++){
            if(atma != ok[en][1] && atma != ok[en][2] && atma != ok[en][3] && atma != ok[en][4])
            {
              cek = cek + 1;
            }
        }

    if(cek == i){

    for(atmb = 0; atmb <= (nAtomS1 - 3); atmb = atmb + 3)
    {
        if(atmb != atma){
        
        bonded = 0;
        for(crt = 0; crt < s1s1hbdn; crt++)
        {

           dist = distance(atmS1, atma, atmb, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atma, atma, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        } 
        if(bonded != 0){

    for(atmc = 0; atmc <= (nAtomS1 - 3); atmc = atmc + 3) 
    {
      if(atmc != atma && atmc != atmb){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmb, atmc, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){
        
        distOaOb = distanceOxOy(atmS1, atma, atmb, opos);
        distOaOc = distanceOxOy(atmS1, atma, atmc, opos);
        distObOc = distanceOxOy(atmS1, atmb, atmc, opos);
        angleOaObOc = (acos((pow(distOaOb,2)+pow(distObOc,2)-pow(distOaOc,2))/(2*distOaOb*distObOc))/PI)*180;
        if(angleOaObOc < 100 && angleOaObOc > 80){

    for(atmd = 0; atmd <= (nAtomS1 - 3); atmd = atmd + 3)
    {
      if(atmd != atma && atmd != atmb && atmd != atmc){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmc, atmd, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmd, atmd, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){

        distOcOd = distanceOxOy(atmS1, atmc, atmd, opos);
        distObOd = distanceOxOy(atmS1, atmb, atmd, opos);
        angleObOcOd = (acos((pow(distObOc,2)+pow(distOcOd,2)-pow(distObOd,2))/(2*distObOc*distOcOd))/PI)*180;

        if(angleObOcOd < 100 && angleObOcOd > 80){

        if(fabs(dihedral(atmS1,atma,atmb,atmc,atmd,opos)) <= 20){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmd, atma, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atma, atma, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmd, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmd, atmd, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmd, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){

            distOdOa = distanceOxOy(atmS1, atmd, atma, opos);
            distOcOa = distanceOxOy(atmS1, atmc, atma, opos);
            angleOcOdOa = (acos((pow(distOcOd,2)+pow(distOdOa,2)-pow(distOcOa,2))/(2*distOcOd*distOdOa))/PI)*180;

            if(angleOcOdOa < 100 && angleOcOdOa > 80){

            distOdOb = distanceOxOy(atmS1, atmd, atmb, opos);
            angleOdOaOb = (acos((pow(distOdOa,2)+pow(distOaOb,2)-pow(distOdOb,2))/(2*distOdOa*distOaOb))/PI)*180;

            if(angleOdOaOb < 100 && angleOdOaOb > 80){

            ok[i][1] = atma; ok[i][2] = atmb; ok[i][3] = atmc; ok[i][4] = atmd;

            v4[0]=atma; v4[1]=atmb; v4[2]=atmc; v4[3]=atmd;

            v2[cntr]=atma; v2[cntr+1]=atmb; v2[cntr+2]=atmc; v2[cntr+3]=atmd;

            cntr = cntr + 4;

            e = 0;

            if(perm4(v4, 4, 0, v2, &e)==1){ fprintf(outputfTetramer,"%d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3);

            if(hiso == 1)
            { 
              if(check_isolated_tetramer(atmS1,nAtomS1,s1s1hbdn,s1a,s1b,s1as1bBD,s1s1hban,s1s1v1,s1s1v2,s1s1v3,s1s1v4,s1s1v5,s1s1v6,atma,atmb,atmc,atmd) == 1)
         
                 fprintf(outputfTetramerIso,"%d %d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3,1+atmd/3);
            }
           }
          }
          }
          }
         }
         }
         }
        }
        }
        }
       }
       }
      }
      }
     }
    }
   }
  }

if(hiso == 1) fclose(outputfTetramerIso);
}

// Water Cyclic Trimer search

void ChemNetworkNew::htrimersearch(double *atmS1, int nd1, int opos, int nsolvent1, int nAtomS1, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, 
                   int s1s1hban, int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, int hiso,
                   FILE *outputfTrimer, char *foutputTrimerIso)
{

   FILE *outputfTrimerIso;
   int nodei, i, j, crt;
   double dist, dist2, hyptns, ang;
   int atma, atmb, atmc, bonded; 
   int ok[MAX][4], cntr, in, v2[100000], v3[3], en, cek, e;
   double distOaOb, distOaOc, distObOc;
   double angleOaObOc, distOcOa, angleObOcOa, angleOcOaOb;

   nodei = nd1;

   if(hiso == 1) outputfTrimerIso = fopen(foutputTrimerIso,"w"); 
  
   for(i = 0; i < MAX; i++)
       for(j = 0; j < 4; j++)
           ok[i][j] = 0;
    i = 0;          
    cntr = 0; for(in = 0; in < 100000; in++) v2[in] = 0; for(in = 0; in < 3; in++) v3[in] = 0; 


    for(atma = 0; atma <= (nAtomS1 - 3); atma = atma + 3){ i = i + 1; cek = 0;
        for(en = 1; en <= i; en++){
            if(atma != ok[en][1] && atma != ok[en][2] && atma != ok[en][3])
            {
              cek = cek + 1;
            }
        }

    if(cek == i){

    for(atmb = 0; atmb <= (nAtomS1 - 3); atmb = atmb + 3)
    {
        if(atmb != atma){
        
        bonded = 0;
        for(crt = 0; crt < s1s1hbdn; crt++)
        {

           dist = distance(atmS1, atma, atmb, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atma, atma, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        } 
        if(bonded != 0){

    for(atmc = 0; atmc <= (nAtomS1 - 3); atmc = atmc + 3) 
    {
      if(atmc != atma && atmc != atmb){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmb, atmc, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){
        
        distOaOb = distanceOxOy(atmS1, atma, atmb, opos);
        distOaOc = distanceOxOy(atmS1, atma, atmc, opos);
        distObOc = distanceOxOy(atmS1, atmb, atmc, opos);
        angleOaObOc = (acos((pow(distOaOb,2)+pow(distObOc,2)-pow(distOaOc,2))/(2*distOaOb*distObOc))/PI)*180;
        if(angleOaObOc < 70 && angleOaObOc > 50){


      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmc, atma, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atma, atma, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmc, atma, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmc, atma, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){

        distOcOa = distanceOxOy(atmS1, atmc, atma, opos);
        angleObOcOa = (acos((pow(distObOc,2)+pow(distOcOa,2)-pow(distOaOb,2))/(2*distObOc*distOcOa))/PI)*180;

        if(angleObOcOa < 70 && angleObOcOa > 50){

        
        angleOcOaOb = (acos((pow(distOcOa,2)+pow(distOaOb,2)-pow(distObOc,2))/(2*distOcOa*distOaOb))/PI)*180;

        if(angleOcOaOb < 70 && angleOcOaOb > 50){

       

            ok[i][1] = atma; ok[i][2] = atmb; ok[i][3] = atmc;

            v3[0]=atma; v3[1]=atmb; v3[2]=atmc;

            v2[cntr]=atma; v2[cntr+1]=atmb; v2[cntr+2]=atmc;

            cntr = cntr + 3;

            e = 0;

            if(perm3(v3, 3, 0, v2, &e)==1){ fprintf(outputfTrimer,"%d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3);

            if(hiso == 1)
            { 
              if(check_isolated_trimer(atmS1,nAtomS1,s1s1hbdn,s1a,s1b,s1as1bBD,s1s1hban,s1s1v1,s1s1v2,s1s1v3,s1s1v4,s1s1v5,s1s1v6,atma,atmb,atmc) == 1)
         
                 fprintf(outputfTrimerIso,"%d %d %d\n",1+atma/3,1+atmb/3,1+atmc/3);
            }
           }
          }
          }
          }
         }
         }
         }
        }
        }
        }
       }
       }
      }


if(hiso == 1) fclose(outputfTrimerIso);
}



// Isolated Hexamer search

int ChemNetworkNew::check_isolated_hexamer(double *atmS, int nAtomS, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, int s1s1hban, 
                           int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, 
                           int atma, int atmb, int atmc, int atmd, int atme, int atmf)
{
 /* Returns 1 if Hexamer is Isolated */

 int atm, test, bonded, crt;
 double dist, dist2, hyptns, ang;
       
 test = 0;

 /*** Atm a ***/

 for(atm = 0; atm <= (nAtomS-3); atm = atm + 3){

    if(atm != atma && atm != atmb && atm != atmc && atm != atmd && atm != atme && atm != atmf){

    bonded = 0;
    for(crt = 0; crt < s1s1hbdn; crt++)
    {

       dist = distance(atmS, atm, atma, s1a, s1b, crt);

       if(dist < s1as1bBD[crt])
       {
            if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
            {
               dist2 = distance(atmS, atma, atma, s1s1v4, s1s1v5, crt);
               hyptns = distance(atmS, atm, atma, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns); 
            }
            if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
            {
               dist2 = distance(atmS, atm, atm, s1s1v3, s1s1v4, crt);
               hyptns = distance(atmS, atm, atma, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns); 
            }
            if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
            {
               bonded = bonded + 1;
            }

        }

    } 
    if(bonded != 0)  test = test + 1;
  }
 }
 
 /*** Atm a ***/

 /*** Atm b ***/

 for(atm = 0; atm <= (nAtomS - 3); atm = atm + 3){

   if(atm != atma && atm != atmb && atm != atmc && atm != atmd && atm != atme && atm != atmf){

    bonded = 0;
    for(crt = 0; crt < s1s1hbdn; crt++)
    {

       dist = distance(atmS, atm, atmb, s1a, s1b, crt);

       if(dist < s1as1bBD[crt])
       {
            if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
            {
               dist2 = distance(atmS, atmb, atmb, s1s1v4, s1s1v5, crt);
               hyptns = distance(atmS, atm, atmb, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
            {
               dist2 = distance(atmS, atm, atm, s1s1v3, s1s1v4, crt);
               hyptns = distance(atmS, atm, atmb, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
            {
               bonded = bonded + 1;
            }

        }

    }
    if(bonded != 0)  test = test + 1;
  }
 }

 /*** Atm b ***/

 /*** Atm c ***/

 for(atm = 0; atm <= (nAtomS - 3); atm = atm + 3){

   if(atm != atma && atm != atmb && atm != atmc && atm != atmd && atm != atme && atm != atmf){

    bonded = 0;
    for(crt = 0; crt < s1s1hbdn; crt++)
    {

       dist = distance(atmS, atm, atmc, s1a, s1b, crt);

       if(dist < s1as1bBD[crt])
       {
            if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
            {
               dist2 = distance(atmS, atmc, atmc, s1s1v4, s1s1v5, crt);
               hyptns = distance(atmS, atm, atmc, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
            {
               dist2 = distance(atmS, atm, atm, s1s1v3, s1s1v4, crt);
               hyptns = distance(atmS, atm, atmc, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
            {
               bonded = bonded + 1;
            }

        }

    }
    if(bonded != 0)  test = test + 1;
  }
 }

 /*** Atm c ***/

 /*** Atm d ***/

 for(atm = 0; atm <= (nAtomS - 3); atm = atm + 3){

   if(atm != atma && atm != atmb && atm != atmc && atm != atmd && atm != atme && atm != atmf){

    bonded = 0;
    for(crt = 0; crt < s1s1hbdn; crt++)
    {

       dist = distance(atmS, atm, atmd, s1a, s1b, crt);

       if(dist < s1as1bBD[crt])
       {
            if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
            {
               dist2 = distance(atmS, atmd, atmd, s1s1v4, s1s1v5, crt);
               hyptns = distance(atmS, atm, atmd, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
            {
               dist2 = distance(atmS, atm, atm, s1s1v3, s1s1v4, crt);
               hyptns = distance(atmS, atm, atmd, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
            {
               bonded = bonded + 1;
            }

        }

    }
    if(bonded != 0)  test = test + 1;
  }
 }

 /*** Atm d ***/

 /*** Atm e ***/

 for(atm = 0; atm <= (nAtomS - 3); atm = atm + 3){

   if(atm != atma && atm != atmb && atm != atmc && atm != atmd && atm != atme && atm != atmf){

    bonded = 0;
    for(crt = 0; crt < s1s1hbdn; crt++)
    {

       dist = distance(atmS, atm, atme, s1a, s1b, crt);

       if(dist < s1as1bBD[crt])
       {
            if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
            {
               dist2 = distance(atmS, atme, atme, s1s1v4, s1s1v5, crt);
               hyptns = distance(atmS, atm, atme, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
            {
               dist2 = distance(atmS, atm, atm, s1s1v3, s1s1v4, crt);
               hyptns = distance(atmS, atm, atme, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
            {
               bonded = bonded + 1;
            }

        }

    }
    if(bonded != 0)  test = test + 1;
  }
 }

 /*** Atm e ***/

 /*** Atm f ***/

 for(atm = 0; atm <= (nAtomS - 3); atm = atm + 3){

   if(atm != atma && atm != atmb && atm != atmc && atm != atmd && atm != atme && atm != atmf){

    bonded = 0;
    for(crt = 0; crt < s1s1hbdn; crt++)
    {

       dist = distance(atmS, atm, atmf, s1a, s1b, crt);

       if(dist < s1as1bBD[crt])
       {
            if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
            {
               dist2 = distance(atmS, atmf, atmf, s1s1v4, s1s1v5, crt);
               hyptns = distance(atmS, atm, atmf, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
            {
               dist2 = distance(atmS, atm, atm, s1s1v3, s1s1v4, crt);
               hyptns = distance(atmS, atm, atmf, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
            {
               bonded = bonded + 1;
            }

        }

    }
    if(bonded != 0)  test = test + 1;
  }
 }

 /*** Atm f ***/

if(test == 0) return 1;
else return 0;

}


// Isolated Pentamer search

int ChemNetworkNew::check_isolated_pentamer(double *atmS, int nAtomS, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, int s1s1hban, 
                            int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, 
                            int atma, int atmb, int atmc, int atmd, int atme)
{
 /* Returns 1 if Pentamer is Isolated */

 int atm, test, bonded, crt;
 double dist, dist2, hyptns, ang;
       
 test = 0;

 /*** Atm a ***/

 for(atm = 0; atm <= (nAtomS-3); atm = atm + 3){

    if(atm != atma && atm != atmb && atm != atmc && atm != atmd && atm != atme){

    bonded = 0;
    for(crt = 0; crt < s1s1hbdn; crt++)
    {

       dist = distance(atmS, atm, atma, s1a, s1b, crt);

       if(dist < s1as1bBD[crt])
       {
            if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
            {
               dist2 = distance(atmS, atma, atma, s1s1v4, s1s1v5, crt);
               hyptns = distance(atmS, atm, atma, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns); 
            }
            if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
            {
               dist2 = distance(atmS, atm, atm, s1s1v3, s1s1v4, crt);
               hyptns = distance(atmS, atm, atma, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns); 
            }
            if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
            {
               bonded = bonded + 1;
            }

        }

    } 
    if(bonded != 0)  test = test + 1;
  }
 }
 
 /*** Atm a ***/

 /*** Atm b ***/

 for(atm = 0; atm <= (nAtomS - 3); atm = atm + 3){

   if(atm != atma && atm != atmb && atm != atmc && atm != atmd && atm != atme){

    bonded = 0;
    for(crt = 0; crt < s1s1hbdn; crt++)
    {

       dist = distance(atmS, atm, atmb, s1a, s1b, crt);

       if(dist < s1as1bBD[crt])
       {
            if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
            {
               dist2 = distance(atmS, atmb, atmb, s1s1v4, s1s1v5, crt);
               hyptns = distance(atmS, atm, atmb, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
            {
               dist2 = distance(atmS, atm, atm, s1s1v3, s1s1v4, crt);
               hyptns = distance(atmS, atm, atmb, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
            {
               bonded = bonded + 1;
            }

        }

    }
    if(bonded != 0)  test = test + 1;
  }
 }

 /*** Atm b ***/

 /*** Atm c ***/

 for(atm = 0; atm <= (nAtomS - 3); atm = atm + 3){

   if(atm != atma && atm != atmb && atm != atmc && atm != atmd && atm != atme){

    bonded = 0;
    for(crt = 0; crt < s1s1hbdn; crt++)
    {

       dist = distance(atmS, atm, atmc, s1a, s1b, crt);

       if(dist < s1as1bBD[crt])
       {
            if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
            {
               dist2 = distance(atmS, atmc, atmc, s1s1v4, s1s1v5, crt);
               hyptns = distance(atmS, atm, atmc, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
            {
               dist2 = distance(atmS, atm, atm, s1s1v3, s1s1v4, crt);
               hyptns = distance(atmS, atm, atmc, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
            {
               bonded = bonded + 1;
            }

        }

    }
    if(bonded != 0)  test = test + 1;
  }
 }

 /*** Atm c ***/

 /*** Atm d ***/

 for(atm = 0; atm <= (nAtomS - 3); atm = atm + 3){

   if(atm != atma && atm != atmb && atm != atmc && atm != atmd && atm != atme){

    bonded = 0;
    for(crt = 0; crt < s1s1hbdn; crt++)
    {

       dist = distance(atmS, atm, atmd, s1a, s1b, crt);

       if(dist < s1as1bBD[crt])
       {
            if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
            {
               dist2 = distance(atmS, atmd, atmd, s1s1v4, s1s1v5, crt);
               hyptns = distance(atmS, atm, atmd, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
            {
               dist2 = distance(atmS, atm, atm, s1s1v3, s1s1v4, crt);
               hyptns = distance(atmS, atm, atmd, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
            {
               bonded = bonded + 1;
            }

        }

    }
    if(bonded != 0)  test = test + 1;
  }
 }

 /*** Atm d ***/

 /*** Atm e ***/

 for(atm = 0; atm <= (nAtomS - 3); atm = atm + 3){

   if(atm != atma && atm != atmb && atm != atmc && atm != atmd && atm != atme){

    bonded = 0;
    for(crt = 0; crt < s1s1hbdn; crt++)
    {

       dist = distance(atmS, atm, atme, s1a, s1b, crt);

       if(dist < s1as1bBD[crt])
       {
            if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
            {
               dist2 = distance(atmS, atme, atme, s1s1v4, s1s1v5, crt);
               hyptns = distance(atmS, atm, atme, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
            {
               dist2 = distance(atmS, atm, atm, s1s1v3, s1s1v4, crt);
               hyptns = distance(atmS, atm, atme, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
            {
               bonded = bonded + 1;
            }

        }

    }
    if(bonded != 0)  test = test + 1;
  }
 }

 /*** Atm e ***/


if(test == 0) return 1;
else return 0;

}


// Isolated Tetramer search

int ChemNetworkNew::check_isolated_tetramer(double *atmS, int nAtomS, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, int s1s1hban, 
                            int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, 
                            int atma, int atmb, int atmc, int atmd)
{
 /* Returns 1 if Tetramer is Isolated */

 int atm, test, bonded, crt;
 double dist, dist2, hyptns, ang;
       
 test = 0;

 /*** Atm a ***/

 for(atm = 0; atm <= (nAtomS-3); atm = atm + 3){

    if(atm != atma && atm != atmb && atm != atmc && atm != atmd){

    bonded = 0;
    for(crt = 0; crt < s1s1hbdn; crt++)
    {

       dist = distance(atmS, atm, atma, s1a, s1b, crt);

       if(dist < s1as1bBD[crt])
       {
            if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
            {
               dist2 = distance(atmS, atma, atma, s1s1v4, s1s1v5, crt);
               hyptns = distance(atmS, atm, atma, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns); 
            }
            if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
            {
               dist2 = distance(atmS, atm, atm, s1s1v3, s1s1v4, crt);
               hyptns = distance(atmS, atm, atma, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns); 
            }
            if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
            {
               bonded = bonded + 1;
            }

        }

    } 
    if(bonded != 0)  test = test + 1;
  }
 }
 
 /*** Atm a ***/

 /*** Atm b ***/

 for(atm = 0; atm <= (nAtomS - 3); atm = atm + 3){

   if(atm != atma && atm != atmb && atm != atmc && atm != atmd){

    bonded = 0;
    for(crt = 0; crt < s1s1hbdn; crt++)
    {

       dist = distance(atmS, atm, atmb, s1a, s1b, crt);

       if(dist < s1as1bBD[crt])
       {
            if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
            {
               dist2 = distance(atmS, atmb, atmb, s1s1v4, s1s1v5, crt);
               hyptns = distance(atmS, atm, atmb, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
            {
               dist2 = distance(atmS, atm, atm, s1s1v3, s1s1v4, crt);
               hyptns = distance(atmS, atm, atmb, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
            {
               bonded = bonded + 1;
            }

        }

    }
    if(bonded != 0)  test = test + 1;
  }
 }

 /*** Atm b ***/

 /*** Atm c ***/

 for(atm = 0; atm <= (nAtomS - 3); atm = atm + 3){

   if(atm != atma && atm != atmb && atm != atmc && atm != atmd){

    bonded = 0;
    for(crt = 0; crt < s1s1hbdn; crt++)
    {

       dist = distance(atmS, atm, atmc, s1a, s1b, crt);

       if(dist < s1as1bBD[crt])
       {
            if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
            {
               dist2 = distance(atmS, atmc, atmc, s1s1v4, s1s1v5, crt);
               hyptns = distance(atmS, atm, atmc, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
            {
               dist2 = distance(atmS, atm, atm, s1s1v3, s1s1v4, crt);
               hyptns = distance(atmS, atm, atmc, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
            {
               bonded = bonded + 1;
            }

        }

    }
    if(bonded != 0)  test = test + 1;
  }
 }

 /*** Atm c ***/

 /*** Atm d ***/

 for(atm = 0; atm <= (nAtomS - 3); atm = atm + 3){

   if(atm != atma && atm != atmb && atm != atmc && atm != atmd){

    bonded = 0;
    for(crt = 0; crt < s1s1hbdn; crt++)
    {

       dist = distance(atmS, atm, atmd, s1a, s1b, crt);

       if(dist < s1as1bBD[crt])
       {
            if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
            {
               dist2 = distance(atmS, atmd, atmd, s1s1v4, s1s1v5, crt);
               hyptns = distance(atmS, atm, atmd, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
            {
               dist2 = distance(atmS, atm, atm, s1s1v3, s1s1v4, crt);
               hyptns = distance(atmS, atm, atmd, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
            {
               bonded = bonded + 1;
            }

        }

    }
    if(bonded != 0)  test = test + 1;
  }
 }

 /*** Atm d ***/


if(test == 0) return 1;
else return 0;

}

// Isolated Trimer search

int ChemNetworkNew::check_isolated_trimer(double *atmS, int nAtomS, int s1s1hbdn, int *s1a, int *s1b, double *s1as1bBD, int s1s1hban, 
                          int *s1s1v1, int *s1s1v2, int *s1s1v3, int *s1s1v4, int *s1s1v5, double *s1s1v6, 
                          int atma, int atmb, int atmc)
{
 /* Returns 1 if Trimer is Isolated */

 int atm, test, bonded, crt;
 double dist, dist2, hyptns, ang;
       
 test = 0;

 /*** Atm a ***/

 for(atm = 0; atm <= (nAtomS-3); atm = atm + 3){

    if(atm != atma && atm != atmb && atm != atmc){

    bonded = 0;
    for(crt = 0; crt < s1s1hbdn; crt++)
    {

       dist = distance(atmS, atm, atma, s1a, s1b, crt);

       if(dist < s1as1bBD[crt])
       {
            if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
            {
               dist2 = distance(atmS, atma, atma, s1s1v4, s1s1v5, crt);
               hyptns = distance(atmS, atm, atma, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns); 
            }
            if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
            {
               dist2 = distance(atmS, atm, atm, s1s1v3, s1s1v4, crt);
               hyptns = distance(atmS, atm, atma, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns); 
            }
            if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
            {
               bonded = bonded + 1;
            }

        }

    } 
    if(bonded != 0)  test = test + 1;
  }
 }
 
 /*** Atm a ***/

 /*** Atm b ***/

 for(atm = 0; atm <= (nAtomS - 3); atm = atm + 3){

   if(atm != atma && atm != atmb && atm != atmc){

    bonded = 0;
    for(crt = 0; crt < s1s1hbdn; crt++)
    {

       dist = distance(atmS, atm, atmb, s1a, s1b, crt);

       if(dist < s1as1bBD[crt])
       {
            if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
            {
               dist2 = distance(atmS, atmb, atmb, s1s1v4, s1s1v5, crt);
               hyptns = distance(atmS, atm, atmb, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
            {
               dist2 = distance(atmS, atm, atm, s1s1v3, s1s1v4, crt);
               hyptns = distance(atmS, atm, atmb, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
            {
               bonded = bonded + 1;
            }

        }

    }
    if(bonded != 0)  test = test + 1;
  }
 }

 /*** Atm b ***/

 /*** Atm c ***/

 for(atm = 0; atm <= (nAtomS - 3); atm = atm + 3){

   if(atm != atma && atm != atmb && atm != atmc){

    bonded = 0;
    for(crt = 0; crt < s1s1hbdn; crt++)
    {

       dist = distance(atmS, atm, atmc, s1a, s1b, crt);

       if(dist < s1as1bBD[crt])
       {
            if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
            {
               dist2 = distance(atmS, atmc, atmc, s1s1v4, s1s1v5, crt);
               hyptns = distance(atmS, atm, atmc, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
            {
               dist2 = distance(atmS, atm, atm, s1s1v3, s1s1v4, crt);
               hyptns = distance(atmS, atm, atmc, s1s1v3, s1s1v5, crt);
               ang = angle(dist, dist2, hyptns);
            }
            if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
            {
               bonded = bonded + 1;
            }

        }

    }
    if(bonded != 0)  test = test + 1;
  }
 }

 /*** Atm c ***/


if(test == 0) return 1;
else return 0;

}


int ChemNetworkNew::prismbook_check(double *atmSX, int ndX, int opos, int sXsXhbdn, int *sXa, int *sXb, double *sXasXbBD, 
                    int sXsXhban, int *sXsXv1, int *sXsXv2, int *sXsXv3, int *sXsXv4, int *sXsXv5, double *sXsXv6,
                    int atmaX, int atmbX, int atmcX, int atmdX, int atmeX, int atmfX)
{
   /* Returns 1, if atoms "a b c d e f" form a PRISM/BOOK hexamer with the same criteria set for the PRISM/BOOK HEXAMER */

   int nodei, i, crt;
   double dist, dist2, hyptns, ang;
   int nAtomS1, atma, atmb, atmc, atmd, atme, atmf, bonded, test; 
   double distOaOb, distOaOc, distObOc, distOcOd, distObOd, distOcOe, distOdOe, distOdOf, distOeOf;
   double angleOaObOc, angleObOcOd, angleOcOdOe, angleOeOdOf;
   int s1s1hbdn, s1s1hban;
   int s1a[50], s1b[50], s1s1v1[50], s1s1v2[50], s1s1v3[50], s1s1v4[50], s1s1v5[50];
   double s1as1bBD[50], s1s1v6[50], *atmS1;

   nodei = ndX;
   s1s1hbdn = sXsXhbdn;
   s1s1hban = sXsXhban;
   nAtomS1 = 18;
 
   test = 0; 
 
for(i=0; i < s1s1hbdn; i++)
{

   s1a[i] = sXa[i];  s1b[i] = sXb[i];  s1as1bBD[i] = sXasXbBD[i];

}


for(i=0; i < s1s1hban; i++)
{

   s1s1v1[i] = sXsXv1[i]; s1s1v2[i] = sXsXv2[i]; s1s1v3[i] = sXsXv3[i]; 
   s1s1v4[i] = sXsXv4[i]; s1s1v5[i] = sXsXv5[i]; s1s1v6[i] = sXsXv6[i];

}


 
   atmS1 = (double*)malloc((3*nAtomS1)*sizeof(double));

if(opos == 1)
{
   atmS1[0] = atmSX[3*(atmaX + opos + 0) - 3];    atmS1[9] = atmSX[3*(atmbX + opos + 0) - 3];   atmS1[18] = atmSX[3*(atmcX + opos + 0) - 3];
   atmS1[1] = atmSX[3*(atmaX + opos + 0) - 2];   atmS1[10] = atmSX[3*(atmbX + opos + 0) - 2];   atmS1[19] = atmSX[3*(atmcX + opos + 0) - 2];
   atmS1[2] = atmSX[3*(atmaX + opos + 0) - 1];   atmS1[11] = atmSX[3*(atmbX + opos + 0) - 1];   atmS1[20] = atmSX[3*(atmcX + opos + 0) - 1];
   atmS1[3] = atmSX[3*(atmaX + opos + 1) - 3];   atmS1[12] = atmSX[3*(atmbX + opos + 1) - 3];   atmS1[21] = atmSX[3*(atmcX + opos + 1) - 3];
   atmS1[4] = atmSX[3*(atmaX + opos + 1) - 2];   atmS1[13] = atmSX[3*(atmbX + opos + 1) - 2];   atmS1[22] = atmSX[3*(atmcX + opos + 1) - 2];
   atmS1[5] = atmSX[3*(atmaX + opos + 1) - 1];   atmS1[14] = atmSX[3*(atmbX + opos + 1) - 1];   atmS1[23] = atmSX[3*(atmcX + opos + 1) - 1];
   atmS1[6] = atmSX[3*(atmaX + opos + 2) - 3];   atmS1[15] = atmSX[3*(atmbX + opos + 2) - 3];   atmS1[24] = atmSX[3*(atmcX + opos + 2) - 3];
   atmS1[7] = atmSX[3*(atmaX + opos + 2) - 2];   atmS1[16] = atmSX[3*(atmbX + opos + 2) - 2];   atmS1[25] = atmSX[3*(atmcX + opos + 2) - 2];
   atmS1[8] = atmSX[3*(atmaX + opos + 2) - 1];   atmS1[17] = atmSX[3*(atmbX + opos + 2) - 1];   atmS1[26] = atmSX[3*(atmcX + opos + 2) - 1];

   atmS1[27] = atmSX[3*(atmdX + opos + 0) - 3];   atmS1[36] = atmSX[3*(atmeX + opos + 0) - 3];   atmS1[45] = atmSX[3*(atmfX + opos + 0) - 3];
   atmS1[28] = atmSX[3*(atmdX + opos + 0) - 2];   atmS1[37] = atmSX[3*(atmeX + opos + 0) - 2];   atmS1[46] = atmSX[3*(atmfX + opos + 0) - 2];
   atmS1[29] = atmSX[3*(atmdX + opos + 0) - 1];   atmS1[38] = atmSX[3*(atmeX + opos + 0) - 1];   atmS1[47] = atmSX[3*(atmfX + opos + 0) - 1];
   atmS1[30] = atmSX[3*(atmdX + opos + 1) - 3];   atmS1[39] = atmSX[3*(atmeX + opos + 1) - 3];   atmS1[48] = atmSX[3*(atmfX + opos + 1) - 3];
   atmS1[31] = atmSX[3*(atmdX + opos + 1) - 2];   atmS1[40] = atmSX[3*(atmeX + opos + 1) - 2];   atmS1[49] = atmSX[3*(atmfX + opos + 1) - 2];
   atmS1[32] = atmSX[3*(atmdX + opos + 1) - 1];   atmS1[41] = atmSX[3*(atmeX + opos + 1) - 1];   atmS1[50] = atmSX[3*(atmfX + opos + 1) - 1];
   atmS1[33] = atmSX[3*(atmdX + opos + 2) - 3];   atmS1[42] = atmSX[3*(atmeX + opos + 2) - 3];   atmS1[51] = atmSX[3*(atmfX + opos + 2) - 3];
   atmS1[34] = atmSX[3*(atmdX + opos + 2) - 2];   atmS1[43] = atmSX[3*(atmeX + opos + 2) - 2];   atmS1[52] = atmSX[3*(atmfX + opos + 2) - 2];
   atmS1[35] = atmSX[3*(atmdX + opos + 2) - 1];   atmS1[44] = atmSX[3*(atmeX + opos + 2) - 1];   atmS1[53] = atmSX[3*(atmfX + opos + 2) - 1];

}

if(opos == 2)
{
   atmS1[0] = atmSX[3*(atmaX + opos - 1) - 3];    atmS1[9] = atmSX[3*(atmbX + opos - 1) - 3];   atmS1[18] = atmSX[3*(atmcX + opos - 1) - 3];
   atmS1[1] = atmSX[3*(atmaX + opos - 1) - 2];   atmS1[10] = atmSX[3*(atmbX + opos - 1) - 2];   atmS1[19] = atmSX[3*(atmcX + opos - 1) - 2];
   atmS1[2] = atmSX[3*(atmaX + opos - 1) - 1];   atmS1[11] = atmSX[3*(atmbX + opos - 1) - 1];   atmS1[20] = atmSX[3*(atmcX + opos - 1) - 1];
   atmS1[3] = atmSX[3*(atmaX + opos + 0) - 3];   atmS1[12] = atmSX[3*(atmbX + opos + 0) - 3];   atmS1[21] = atmSX[3*(atmcX + opos + 0) - 3];
   atmS1[4] = atmSX[3*(atmaX + opos + 0) - 2];   atmS1[13] = atmSX[3*(atmbX + opos + 0) - 2];   atmS1[22] = atmSX[3*(atmcX + opos + 0) - 2];
   atmS1[5] = atmSX[3*(atmaX + opos + 0) - 1];   atmS1[14] = atmSX[3*(atmbX + opos + 0) - 1];   atmS1[23] = atmSX[3*(atmcX + opos + 0) - 1];
   atmS1[6] = atmSX[3*(atmaX + opos + 1) - 3];   atmS1[15] = atmSX[3*(atmbX + opos + 1) - 3];   atmS1[24] = atmSX[3*(atmcX + opos + 1) - 3];
   atmS1[7] = atmSX[3*(atmaX + opos + 1) - 2];   atmS1[16] = atmSX[3*(atmbX + opos + 1) - 2];   atmS1[25] = atmSX[3*(atmcX + opos + 1) - 2];
   atmS1[8] = atmSX[3*(atmaX + opos + 1) - 1];   atmS1[17] = atmSX[3*(atmbX + opos + 1) - 1];   atmS1[26] = atmSX[3*(atmcX + opos + 1) - 1];

   atmS1[27] = atmSX[3*(atmdX + opos - 1) - 3];   atmS1[36] = atmSX[3*(atmeX + opos - 1) - 3];   atmS1[45] = atmSX[3*(atmfX + opos - 1) - 3];
   atmS1[28] = atmSX[3*(atmdX + opos - 1) - 2];   atmS1[37] = atmSX[3*(atmeX + opos - 1) - 2];   atmS1[46] = atmSX[3*(atmfX + opos - 1) - 2];
   atmS1[29] = atmSX[3*(atmdX + opos - 1) - 1];   atmS1[38] = atmSX[3*(atmeX + opos - 1) - 1];   atmS1[47] = atmSX[3*(atmfX + opos - 1) - 1];
   atmS1[30] = atmSX[3*(atmdX + opos + 0) - 3];   atmS1[39] = atmSX[3*(atmeX + opos + 0) - 3];   atmS1[48] = atmSX[3*(atmfX + opos + 0) - 3];
   atmS1[31] = atmSX[3*(atmdX + opos + 0) - 2];   atmS1[40] = atmSX[3*(atmeX + opos + 0) - 2];   atmS1[49] = atmSX[3*(atmfX + opos + 0) - 2];
   atmS1[32] = atmSX[3*(atmdX + opos + 0) - 1];   atmS1[41] = atmSX[3*(atmeX + opos + 0) - 1];   atmS1[50] = atmSX[3*(atmfX + opos + 0) - 1];
   atmS1[33] = atmSX[3*(atmdX + opos + 1) - 3];   atmS1[42] = atmSX[3*(atmeX + opos + 1) - 3];   atmS1[51] = atmSX[3*(atmfX + opos + 1) - 3];
   atmS1[34] = atmSX[3*(atmdX + opos + 1) - 2];   atmS1[43] = atmSX[3*(atmeX + opos + 1) - 2];   atmS1[52] = atmSX[3*(atmfX + opos + 1) - 2];
   atmS1[35] = atmSX[3*(atmdX + opos + 1) - 1];   atmS1[44] = atmSX[3*(atmeX + opos + 1) - 1];   atmS1[53] = atmSX[3*(atmfX + opos + 1) - 1];

}

if(opos == 3)
{
   atmS1[0] = atmSX[3*(atmaX + opos - 2) - 3];    atmS1[9] = atmSX[3*(atmbX + opos - 2) - 3];   atmS1[18] = atmSX[3*(atmcX + opos - 2) - 3];
   atmS1[1] = atmSX[3*(atmaX + opos - 2) - 2];   atmS1[10] = atmSX[3*(atmbX + opos - 2) - 2];   atmS1[19] = atmSX[3*(atmcX + opos - 2) - 2];
   atmS1[2] = atmSX[3*(atmaX + opos - 2) - 1];   atmS1[11] = atmSX[3*(atmbX + opos - 2) - 1];   atmS1[20] = atmSX[3*(atmcX + opos - 2) - 1];
   atmS1[3] = atmSX[3*(atmaX + opos - 1) - 3];   atmS1[12] = atmSX[3*(atmbX + opos - 1) - 3];   atmS1[21] = atmSX[3*(atmcX + opos - 1) - 3];
   atmS1[4] = atmSX[3*(atmaX + opos - 1) - 2];   atmS1[13] = atmSX[3*(atmbX + opos - 1) - 2];   atmS1[22] = atmSX[3*(atmcX + opos - 1) - 2];
   atmS1[5] = atmSX[3*(atmaX + opos - 1) - 1];   atmS1[14] = atmSX[3*(atmbX + opos - 1) - 1];   atmS1[23] = atmSX[3*(atmcX + opos - 1) - 1];
   atmS1[6] = atmSX[3*(atmaX + opos + 0) - 3];   atmS1[15] = atmSX[3*(atmbX + opos + 0) - 3];   atmS1[24] = atmSX[3*(atmcX + opos + 0) - 3];
   atmS1[7] = atmSX[3*(atmaX + opos + 0) - 2];   atmS1[16] = atmSX[3*(atmbX + opos + 0) - 2];   atmS1[25] = atmSX[3*(atmcX + opos + 0) - 2];
   atmS1[8] = atmSX[3*(atmaX + opos + 0) - 1];   atmS1[17] = atmSX[3*(atmbX + opos + 0) - 1];   atmS1[26] = atmSX[3*(atmcX + opos + 0) - 1];

   atmS1[27] = atmSX[3*(atmdX + opos - 2) - 3];   atmS1[36] = atmSX[3*(atmeX + opos - 2) - 3];   atmS1[45] = atmSX[3*(atmfX + opos - 2) - 3];
   atmS1[28] = atmSX[3*(atmdX + opos - 2) - 2];   atmS1[37] = atmSX[3*(atmeX + opos - 2) - 2];   atmS1[46] = atmSX[3*(atmfX + opos - 2) - 2];
   atmS1[29] = atmSX[3*(atmdX + opos - 2) - 1];   atmS1[38] = atmSX[3*(atmeX + opos - 2) - 1];   atmS1[47] = atmSX[3*(atmfX + opos - 2) - 1];
   atmS1[30] = atmSX[3*(atmdX + opos - 1) - 3];   atmS1[39] = atmSX[3*(atmeX + opos - 1) - 3];   atmS1[48] = atmSX[3*(atmfX + opos - 1) - 3];
   atmS1[31] = atmSX[3*(atmdX + opos - 1) - 2];   atmS1[40] = atmSX[3*(atmeX + opos - 1) - 2];   atmS1[49] = atmSX[3*(atmfX + opos - 1) - 2];
   atmS1[32] = atmSX[3*(atmdX + opos - 1) - 1];   atmS1[41] = atmSX[3*(atmeX + opos - 1) - 1];   atmS1[50] = atmSX[3*(atmfX + opos - 1) - 1];
   atmS1[33] = atmSX[3*(atmdX + opos + 0) - 3];   atmS1[42] = atmSX[3*(atmeX + opos + 0) - 3];   atmS1[51] = atmSX[3*(atmfX + opos + 0) - 3];
   atmS1[34] = atmSX[3*(atmdX + opos + 0) - 2];   atmS1[43] = atmSX[3*(atmeX + opos + 0) - 2];   atmS1[52] = atmSX[3*(atmfX + opos + 0) - 2];
   atmS1[35] = atmSX[3*(atmdX + opos + 0) - 1];   atmS1[44] = atmSX[3*(atmeX + opos + 0) - 1];   atmS1[53] = atmSX[3*(atmfX + opos + 0) - 1];

}


    for(atma = 0; atma <= (nAtomS1 - 3); atma = atma + 3){ 


    for(atmb = 0; atmb <= (nAtomS1 - 3); atmb = atmb + 3)
    {
        if(atmb != atma){
        
        bonded = 0;
        for(crt = 0; crt < s1s1hbdn; crt++)
        {

           dist = distance(atmS1, atma, atmb, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atma, atma, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        } 
        if(bonded != 0){

    for(atmc = 0; atmc <= (nAtomS1 - 3); atmc = atmc + 3) 
    {
      if(atmc != atma && atmc != atmb){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmb, atmc, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){
        
        distOaOb = distanceOxOy(atmS1, atma, atmb, opos);
        distOaOc = distanceOxOy(atmS1, atma, atmc, opos);
        distObOc = distanceOxOy(atmS1, atmb, atmc, opos);
        angleOaObOc = (acos((pow(distOaOb,2)+pow(distObOc,2)-pow(distOaOc,2))/(2*distOaOb*distObOc))/PI)*180;
        if(angleOaObOc < 70 && angleOaObOc > 50){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmc, atma, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atma, atma, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmc, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmc, atmc, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmc, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){


    for(atmd = 0; atmd <= (nAtomS1 - 3); atmd = atmd + 3)
    {
      if(atmd != atma && atmd != atmb && atmd != atmc){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmc, atmd, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmd, atmd, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){

        distOcOd = distanceOxOy(atmS1, atmc, atmd, opos);
        distObOd = distanceOxOy(atmS1, atmb, atmd, opos);
        angleObOcOd = (acos((pow(distObOc,2)+pow(distOcOd,2)-pow(distObOd,2))/(2*distObOc*distOcOd))/PI)*180;

        if(angleObOcOd < 100 && angleObOcOd > 80){


        for(atme = 0; atme <= (nAtomS1-3); atme = atme + 3)
        {
          if(atme != atma && atme != atmb && atme != atmc && atme != atmd){

          bonded = 0;
          for(crt = 0; crt < s1s1hbdn; crt++)
          {

              dist = distance(atmS1, atmd, atme, s1a, s1b, crt);

              if(dist < s1as1bBD[crt])
              {
                 if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                 {
                    dist2 = distance(atmS1, atme, atme, s1s1v4, s1s1v5, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                 {
                    dist2 = distance(atmS1, atmd, atmd, s1s1v3, s1s1v4, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                 {
                    bonded = bonded + 1;
                 }

             }

          }
          if(bonded != 0){
 
          distOcOe = distanceOxOy(atmS1, atmc, atme, opos);
          distOdOe = distanceOxOy(atmS1, atmd, atme, opos);
          angleOcOdOe = (acos((pow(distOcOd,2)+pow(distOdOe,2)-pow(distOcOe,2))/(2*distOcOd*distOdOe))/PI)*180;

          if(angleOcOdOe < 100 && angleOcOdOe > 80){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atme, atma, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atma, atma, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atme, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atme, atme, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atme, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){


          for(atmf = 0; atmf <= (nAtomS1 - 3); atmf = atmf + 3)
          {
            if(atmf != atma && atmf != atmb && atmf != atmc && atmf != atmd && atmf != atme){

            distOdOf = distanceOxOy(atmS1, atmd, atmf, opos);
            distOeOf = distanceOxOy(atmS1, atme, atmf, opos);
            angleOeOdOf = (acos((pow(distOdOe,2)+pow(distOdOf,2)-pow(distOeOf,2))/(2*distOdOe*distOdOf))/PI)*180;

            if(angleOeOdOf > 70){ 

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmf, atmb, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atmb, atmb, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmf, atmb, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmf, atmb, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){
     
          

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmf, atmd, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atmd, atmd, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmf, atmd, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmf, atmd, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){

              test = test + 1;
            
           
           }
           }
           }
           }
           }
          }
          }
          }
          }
          }
         }
         }
         }
         }
        }
        }
        }
       }
      }
     }
    }
   }
  }

    if(test != 0) return 1;
    else return 0;
    
    free(atmS1);
}


int ChemNetworkNew::ring_check(double *atmSX, int ndX, int opos, int sXsXhbdn, int *sXa, int *sXb, double *sXasXbBD, 
               int sXsXhban, int *sXsXv1, int *sXsXv2, int *sXsXv3, int *sXsXv4, int *sXsXv5, double *sXsXv6,
               int atmaX, int atmbX, int atmcX, int atmdX, int atmeX, int atmfX)
{
   /* Returns 1, if atoms "a b c d e f" form a RING hexamer with the same criteria set for the RING HEXAMER. */

   int nodei, i, crt;
   double dist, dist2, hyptns, ang;
   int nAtomS1, atma, atmb, atmc, atmd, atme, atmf, bonded, test; 
   double distOaOb, distOaOc, distObOc, distOcOd, distObOd, distOcOe, distOdOe, distOdOf, distOeOf;
   double angleOaObOc, angleObOcOd, angleOcOdOe, angleOdOeOf, distOeOa, distOfOa, angleOeOfOa;
   int s1s1hbdn, s1s1hban;
   int s1a[50], s1b[50], s1s1v1[50], s1s1v2[50], s1s1v3[50], s1s1v4[50], s1s1v5[50];
   double s1as1bBD[50], s1s1v6[50], *atmS1;

   nodei = ndX;
   s1s1hbdn = sXsXhbdn;
   s1s1hban = sXsXhban;
   nAtomS1 = 18;
 
   test = 0; 
 
for(i=0; i < s1s1hbdn; i++)
{

   s1a[i] = sXa[i];  s1b[i] = sXb[i];  s1as1bBD[i] = sXasXbBD[i];

}


for(i=0; i < s1s1hban; i++)
{

   s1s1v1[i] = sXsXv1[i]; s1s1v2[i] = sXsXv2[i]; s1s1v3[i] = sXsXv3[i]; 
   s1s1v4[i] = sXsXv4[i]; s1s1v5[i] = sXsXv5[i]; s1s1v6[i] = sXsXv6[i];

}


 
   atmS1 = (double*)malloc((3*nAtomS1)*sizeof(double));

if(opos == 1)
{
   atmS1[0] = atmSX[3*(atmaX + opos + 0) - 3];    atmS1[9] = atmSX[3*(atmbX + opos + 0) - 3];   atmS1[18] = atmSX[3*(atmcX + opos + 0) - 3];
   atmS1[1] = atmSX[3*(atmaX + opos + 0) - 2];   atmS1[10] = atmSX[3*(atmbX + opos + 0) - 2];   atmS1[19] = atmSX[3*(atmcX + opos + 0) - 2];
   atmS1[2] = atmSX[3*(atmaX + opos + 0) - 1];   atmS1[11] = atmSX[3*(atmbX + opos + 0) - 1];   atmS1[20] = atmSX[3*(atmcX + opos + 0) - 1];
   atmS1[3] = atmSX[3*(atmaX + opos + 1) - 3];   atmS1[12] = atmSX[3*(atmbX + opos + 1) - 3];   atmS1[21] = atmSX[3*(atmcX + opos + 1) - 3];
   atmS1[4] = atmSX[3*(atmaX + opos + 1) - 2];   atmS1[13] = atmSX[3*(atmbX + opos + 1) - 2];   atmS1[22] = atmSX[3*(atmcX + opos + 1) - 2];
   atmS1[5] = atmSX[3*(atmaX + opos + 1) - 1];   atmS1[14] = atmSX[3*(atmbX + opos + 1) - 1];   atmS1[23] = atmSX[3*(atmcX + opos + 1) - 1];
   atmS1[6] = atmSX[3*(atmaX + opos + 2) - 3];   atmS1[15] = atmSX[3*(atmbX + opos + 2) - 3];   atmS1[24] = atmSX[3*(atmcX + opos + 2) - 3];
   atmS1[7] = atmSX[3*(atmaX + opos + 2) - 2];   atmS1[16] = atmSX[3*(atmbX + opos + 2) - 2];   atmS1[25] = atmSX[3*(atmcX + opos + 2) - 2];
   atmS1[8] = atmSX[3*(atmaX + opos + 2) - 1];   atmS1[17] = atmSX[3*(atmbX + opos + 2) - 1];   atmS1[26] = atmSX[3*(atmcX + opos + 2) - 1];

   atmS1[27] = atmSX[3*(atmdX + opos + 0) - 3];   atmS1[36] = atmSX[3*(atmeX + opos + 0) - 3];   atmS1[45] = atmSX[3*(atmfX + opos + 0) - 3];
   atmS1[28] = atmSX[3*(atmdX + opos + 0) - 2];   atmS1[37] = atmSX[3*(atmeX + opos + 0) - 2];   atmS1[46] = atmSX[3*(atmfX + opos + 0) - 2];
   atmS1[29] = atmSX[3*(atmdX + opos + 0) - 1];   atmS1[38] = atmSX[3*(atmeX + opos + 0) - 1];   atmS1[47] = atmSX[3*(atmfX + opos + 0) - 1];
   atmS1[30] = atmSX[3*(atmdX + opos + 1) - 3];   atmS1[39] = atmSX[3*(atmeX + opos + 1) - 3];   atmS1[48] = atmSX[3*(atmfX + opos + 1) - 3];
   atmS1[31] = atmSX[3*(atmdX + opos + 1) - 2];   atmS1[40] = atmSX[3*(atmeX + opos + 1) - 2];   atmS1[49] = atmSX[3*(atmfX + opos + 1) - 2];
   atmS1[32] = atmSX[3*(atmdX + opos + 1) - 1];   atmS1[41] = atmSX[3*(atmeX + opos + 1) - 1];   atmS1[50] = atmSX[3*(atmfX + opos + 1) - 1];
   atmS1[33] = atmSX[3*(atmdX + opos + 2) - 3];   atmS1[42] = atmSX[3*(atmeX + opos + 2) - 3];   atmS1[51] = atmSX[3*(atmfX + opos + 2) - 3];
   atmS1[34] = atmSX[3*(atmdX + opos + 2) - 2];   atmS1[43] = atmSX[3*(atmeX + opos + 2) - 2];   atmS1[52] = atmSX[3*(atmfX + opos + 2) - 2];
   atmS1[35] = atmSX[3*(atmdX + opos + 2) - 1];   atmS1[44] = atmSX[3*(atmeX + opos + 2) - 1];   atmS1[53] = atmSX[3*(atmfX + opos + 2) - 1];

}

if(opos == 2)
{
   atmS1[0] = atmSX[3*(atmaX + opos - 1) - 3];    atmS1[9] = atmSX[3*(atmbX + opos - 1) - 3];   atmS1[18] = atmSX[3*(atmcX + opos - 1) - 3];
   atmS1[1] = atmSX[3*(atmaX + opos - 1) - 2];   atmS1[10] = atmSX[3*(atmbX + opos - 1) - 2];   atmS1[19] = atmSX[3*(atmcX + opos - 1) - 2];
   atmS1[2] = atmSX[3*(atmaX + opos - 1) - 1];   atmS1[11] = atmSX[3*(atmbX + opos - 1) - 1];   atmS1[20] = atmSX[3*(atmcX + opos - 1) - 1];
   atmS1[3] = atmSX[3*(atmaX + opos + 0) - 3];   atmS1[12] = atmSX[3*(atmbX + opos + 0) - 3];   atmS1[21] = atmSX[3*(atmcX + opos + 0) - 3];
   atmS1[4] = atmSX[3*(atmaX + opos + 0) - 2];   atmS1[13] = atmSX[3*(atmbX + opos + 0) - 2];   atmS1[22] = atmSX[3*(atmcX + opos + 0) - 2];
   atmS1[5] = atmSX[3*(atmaX + opos + 0) - 1];   atmS1[14] = atmSX[3*(atmbX + opos + 0) - 1];   atmS1[23] = atmSX[3*(atmcX + opos + 0) - 1];
   atmS1[6] = atmSX[3*(atmaX + opos + 1) - 3];   atmS1[15] = atmSX[3*(atmbX + opos + 1) - 3];   atmS1[24] = atmSX[3*(atmcX + opos + 1) - 3];
   atmS1[7] = atmSX[3*(atmaX + opos + 1) - 2];   atmS1[16] = atmSX[3*(atmbX + opos + 1) - 2];   atmS1[25] = atmSX[3*(atmcX + opos + 1) - 2];
   atmS1[8] = atmSX[3*(atmaX + opos + 1) - 1];   atmS1[17] = atmSX[3*(atmbX + opos + 1) - 1];   atmS1[26] = atmSX[3*(atmcX + opos + 1) - 1];

   atmS1[27] = atmSX[3*(atmdX + opos - 1) - 3];   atmS1[36] = atmSX[3*(atmeX + opos - 1) - 3];   atmS1[45] = atmSX[3*(atmfX + opos - 1) - 3];
   atmS1[28] = atmSX[3*(atmdX + opos - 1) - 2];   atmS1[37] = atmSX[3*(atmeX + opos - 1) - 2];   atmS1[46] = atmSX[3*(atmfX + opos - 1) - 2];
   atmS1[29] = atmSX[3*(atmdX + opos - 1) - 1];   atmS1[38] = atmSX[3*(atmeX + opos - 1) - 1];   atmS1[47] = atmSX[3*(atmfX + opos - 1) - 1];
   atmS1[30] = atmSX[3*(atmdX + opos + 0) - 3];   atmS1[39] = atmSX[3*(atmeX + opos + 0) - 3];   atmS1[48] = atmSX[3*(atmfX + opos + 0) - 3];
   atmS1[31] = atmSX[3*(atmdX + opos + 0) - 2];   atmS1[40] = atmSX[3*(atmeX + opos + 0) - 2];   atmS1[49] = atmSX[3*(atmfX + opos + 0) - 2];
   atmS1[32] = atmSX[3*(atmdX + opos + 0) - 1];   atmS1[41] = atmSX[3*(atmeX + opos + 0) - 1];   atmS1[50] = atmSX[3*(atmfX + opos + 0) - 1];
   atmS1[33] = atmSX[3*(atmdX + opos + 1) - 3];   atmS1[42] = atmSX[3*(atmeX + opos + 1) - 3];   atmS1[51] = atmSX[3*(atmfX + opos + 1) - 3];
   atmS1[34] = atmSX[3*(atmdX + opos + 1) - 2];   atmS1[43] = atmSX[3*(atmeX + opos + 1) - 2];   atmS1[52] = atmSX[3*(atmfX + opos + 1) - 2];
   atmS1[35] = atmSX[3*(atmdX + opos + 1) - 1];   atmS1[44] = atmSX[3*(atmeX + opos + 1) - 1];   atmS1[53] = atmSX[3*(atmfX + opos + 1) - 1];

}

if(opos == 3)
{
   atmS1[0] = atmSX[3*(atmaX + opos - 2) - 3];    atmS1[9] = atmSX[3*(atmbX + opos - 2) - 3];   atmS1[18] = atmSX[3*(atmcX + opos - 2) - 3];
   atmS1[1] = atmSX[3*(atmaX + opos - 2) - 2];   atmS1[10] = atmSX[3*(atmbX + opos - 2) - 2];   atmS1[19] = atmSX[3*(atmcX + opos - 2) - 2];
   atmS1[2] = atmSX[3*(atmaX + opos - 2) - 1];   atmS1[11] = atmSX[3*(atmbX + opos - 2) - 1];   atmS1[20] = atmSX[3*(atmcX + opos - 2) - 1];
   atmS1[3] = atmSX[3*(atmaX + opos - 1) - 3];   atmS1[12] = atmSX[3*(atmbX + opos - 1) - 3];   atmS1[21] = atmSX[3*(atmcX + opos - 1) - 3];
   atmS1[4] = atmSX[3*(atmaX + opos - 1) - 2];   atmS1[13] = atmSX[3*(atmbX + opos - 1) - 2];   atmS1[22] = atmSX[3*(atmcX + opos - 1) - 2];
   atmS1[5] = atmSX[3*(atmaX + opos - 1) - 1];   atmS1[14] = atmSX[3*(atmbX + opos - 1) - 1];   atmS1[23] = atmSX[3*(atmcX + opos - 1) - 1];
   atmS1[6] = atmSX[3*(atmaX + opos + 0) - 3];   atmS1[15] = atmSX[3*(atmbX + opos + 0) - 3];   atmS1[24] = atmSX[3*(atmcX + opos + 0) - 3];
   atmS1[7] = atmSX[3*(atmaX + opos + 0) - 2];   atmS1[16] = atmSX[3*(atmbX + opos + 0) - 2];   atmS1[25] = atmSX[3*(atmcX + opos + 0) - 2];
   atmS1[8] = atmSX[3*(atmaX + opos + 0) - 1];   atmS1[17] = atmSX[3*(atmbX + opos + 0) - 1];   atmS1[26] = atmSX[3*(atmcX + opos + 0) - 1];

   atmS1[27] = atmSX[3*(atmdX + opos - 2) - 3];   atmS1[36] = atmSX[3*(atmeX + opos - 2) - 3];   atmS1[45] = atmSX[3*(atmfX + opos - 2) - 3];
   atmS1[28] = atmSX[3*(atmdX + opos - 2) - 2];   atmS1[37] = atmSX[3*(atmeX + opos - 2) - 2];   atmS1[46] = atmSX[3*(atmfX + opos - 2) - 2];
   atmS1[29] = atmSX[3*(atmdX + opos - 2) - 1];   atmS1[38] = atmSX[3*(atmeX + opos - 2) - 1];   atmS1[47] = atmSX[3*(atmfX + opos - 2) - 1];
   atmS1[30] = atmSX[3*(atmdX + opos - 1) - 3];   atmS1[39] = atmSX[3*(atmeX + opos - 1) - 3];   atmS1[48] = atmSX[3*(atmfX + opos - 1) - 3];
   atmS1[31] = atmSX[3*(atmdX + opos - 1) - 2];   atmS1[40] = atmSX[3*(atmeX + opos - 1) - 2];   atmS1[49] = atmSX[3*(atmfX + opos - 1) - 2];
   atmS1[32] = atmSX[3*(atmdX + opos - 1) - 1];   atmS1[41] = atmSX[3*(atmeX + opos - 1) - 1];   atmS1[50] = atmSX[3*(atmfX + opos - 1) - 1];
   atmS1[33] = atmSX[3*(atmdX + opos + 0) - 3];   atmS1[42] = atmSX[3*(atmeX + opos + 0) - 3];   atmS1[51] = atmSX[3*(atmfX + opos + 0) - 3];
   atmS1[34] = atmSX[3*(atmdX + opos + 0) - 2];   atmS1[43] = atmSX[3*(atmeX + opos + 0) - 2];   atmS1[52] = atmSX[3*(atmfX + opos + 0) - 2];
   atmS1[35] = atmSX[3*(atmdX + opos + 0) - 1];   atmS1[44] = atmSX[3*(atmeX + opos + 0) - 1];   atmS1[53] = atmSX[3*(atmfX + opos + 0) - 1];

}


    for(atma = 0; atma <= (nAtomS1 - 3); atma = atma + 3){ 
                   
    for(atmb = 0; atmb <= (nAtomS1 - 3); atmb = atmb + 3)
    {
        if(atmb != atma){
        
        bonded = 0;
        for(crt = 0; crt < s1s1hbdn; crt++)
        {

           dist = distance(atmS1, atma, atmb, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atma, atma, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atma, atmb, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns); 
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        } 
        if(bonded != 0){

    for(atmc = 0; atmc <= (nAtomS1 - 3); atmc = atmc + 3) 
    {
      if(atmc != atma && atmc != atmb){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmb, atmc, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmb, atmb, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmb, atmc, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){
        
        distOaOb = distanceOxOy(atmS1, atma, atmb, opos);
        distOaOc = distanceOxOy(atmS1, atma, atmc, opos);
        distObOc = distanceOxOy(atmS1, atmb, atmc, opos);
        angleOaObOc = (acos((pow(distOaOb,2)+pow(distObOc,2)-pow(distOaOc,2))/(2*distOaOb*distObOc))/PI)*180;
        if(angleOaObOc < 130 && angleOaObOc > 110){

    for(atmd = 0; atmd <= (nAtomS1 - 3); atmd = atmd + 3)
    {
      if(atmd != atma && atmd != atmb && atmd != atmc){

      bonded = 0;
      for(crt = 0; crt < s1s1hbdn; crt++)
      {

           dist = distance(atmS1, atmc, atmd, s1a, s1b, crt);

           if(dist < s1as1bBD[crt])
           {
                if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                {
                   dist2 = distance(atmS1, atmd, atmd, s1s1v4, s1s1v5, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                {
                   dist2 = distance(atmS1, atmc, atmc, s1s1v3, s1s1v4, crt);
                   hyptns = distance(atmS1, atmc, atmd, s1s1v3, s1s1v5, crt);
                   ang = angle(dist, dist2, hyptns);
                }
                if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                {
                   bonded = bonded + 1;
                }

            }

        }
        if(bonded != 0){

        distOcOd = distanceOxOy(atmS1, atmc, atmd, opos);
        distObOd = distanceOxOy(atmS1, atmb, atmd, opos);
        angleObOcOd = (acos((pow(distObOc,2)+pow(distOcOd,2)-pow(distObOd,2))/(2*distObOc*distOcOd))/PI)*180;

        if(angleObOcOd < 130 && angleObOcOd > 110){

        if(fabs(dihedral(atmS1, atma, atmb, atmc, atmd, opos)) <= 30){

        for(atme = 0; atme <= (nAtomS1-3); atme = atme + 3)
        {
          if(atme != atma && atme != atmb && atme != atmc && atme != atmd){

          bonded = 0;
          for(crt = 0; crt < s1s1hbdn; crt++)
          {

              dist = distance(atmS1, atmd, atme, s1a, s1b, crt);

              if(dist < s1as1bBD[crt])
              {
                 if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                 {
                    dist2 = distance(atmS1, atme, atme, s1s1v4, s1s1v5, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                 {
                    dist2 = distance(atmS1, atmd, atmd, s1s1v3, s1s1v4, crt);
                    hyptns = distance(atmS1, atmd, atme, s1s1v3, s1s1v5, crt);
                    ang = angle(dist, dist2, hyptns);
                 }
                 if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                 {
                    bonded = bonded + 1;
                 }

             }

          }
          if(bonded != 0){
 
          distOcOe = distanceOxOy(atmS1, atmc, atme, opos);
          distOdOe = distanceOxOy(atmS1, atmd, atme, opos);
          angleOcOdOe = (acos((pow(distOcOd,2)+pow(distOdOe,2)-pow(distOcOe,2))/(2*distOcOd*distOdOe))/PI)*180;

          if(angleOcOdOe < 130 && angleOcOdOe > 110){

          for(atmf = 0; atmf <= (nAtomS1 - 3); atmf = atmf + 3)
          {
            if(atmf != atma && atmf != atmb && atmf != atmc && atmf != atmd && atmf != atme){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atme, atmf, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atme, atmf, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atme, atme, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atme, atmf, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){
     
            distOdOf = distanceOxOy(atmS1, atmd, atmf, opos);
            distOeOf = distanceOxOy(atmS1, atme, atmf, opos);
            angleOdOeOf = (acos((pow(distOdOe,2)+pow(distOeOf,2)-pow(distOdOf,2))/(2*distOdOe*distOeOf))/PI)*180;

            if(angleOdOeOf < 130 && angleOdOeOf > 110){

            bonded = 0;
            for(crt = 0; crt < s1s1hbdn; crt++)
            {

                dist = distance(atmS1, atmf, atma, s1a, s1b, crt);

                if(dist < s1as1bBD[crt])
                {
                   if(s1s1v1[crt] == 1 && s1s1v2[crt] == 2)
                   {
                      dist2 = distance(atmS1, atma, atma, s1s1v4, s1s1v5, crt);
                      hyptns = distance(atmS1, atmf, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(s1s1v1[crt] == 2 && s1s1v2[crt] == 1)
                   {
                      dist2 = distance(atmS1, atmf, atmf, s1s1v3, s1s1v4, crt);
                      hyptns = distance(atmS1, atmf, atma, s1s1v3, s1s1v5, crt);
                      ang = angle(dist, dist2, hyptns);
                   }
                   if(dist < s1as1bBD[crt] && ang > s1s1v6[crt])
                   {
                      bonded = bonded + 1;
                   }

               }

            }
            if(bonded != 0){

            distOeOa = distanceOxOy(atmS1, atme, atma, opos);
            distOfOa = distanceOxOy(atmS1, atmf, atma, opos);
            angleOeOfOa = (acos((pow(distOeOf,2)+pow(distOfOa,2)-pow(distOeOa,2))/(2*distOeOf*distOfOa))/PI)*180;

            if(angleOeOfOa < 130 && angleOeOfOa > 110){

            if(fabs(dihedral(atmS1,atma,atmf,atme,atmd, opos)) <= 30){     

            if(fabs(dihedral(atmS1,atmf,atma,atmd,atmc, opos)) >= 150 && fabs(dihedral(atmS1,atmb,atma,atmd,atme, opos)) >= 150){      

               test = test + 1;            
           
           
           }
           }
           }
           }
          }
          }
          }
          }
          }
         }
         }
         }
         }
        }
        }
        }
       }
       }
      }
      }
     }
     }
     }
     }
    }
    
    if(test != 0) return 1;
    else return 0;

    free(atmS1);
}




