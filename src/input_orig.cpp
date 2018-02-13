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

#include "chemnetworks_orig.h"

#define FLN 256

using namespace CN_NS;

void ChemNetworkOrig::input_read(int argc, char *argv[]){

// 2016.July : here I want to initilize the cutoffs so that the code can handle the case without cutoffs, for example discard the angle cutoffs in solvent-solvent graphs
    for(i=0;i<NUM_INTER;i++)  // NUM_INTER is defined in common.h, and the array is defined in chemnetworks_orig.h, Tiecheng
    {
        s1s1v6[i] = s2s2v6[i] = s3s3v6[i] = s1s2v6[i] = s1s3v6[i] = s2s3v6[i] = -1.0; // default lower cutoff in angle of solvent-solvent, (real angle is in the range of 0.0 to 180.0)
        s1s1v7[i] = s2s2v7[i] = s3s3v7[i] = s1s2v7[i] = s1s3v7[i] = s2s3v7[i] = 181.0; // default upper cutoff in angle of solvent-solvent, (real angle is in the range of 0.0 to 180.0)
        s1as1bBDmin[i] = s2as2bBDmin[i] = s3as3bBDmin[i] = s12as12bBDmin[i] = s13as13bBDmin[i] = s23as23bBDmin[i] = -0.1; // defalut lower value in distance of solvent-solvent, (real distance is in the range of 0.0 to 100.0)
        s1as1bBDmax[i] = s2as2bBDmax[i] = s3as3bBDmax[i] = s12as12bBDmax[i] = s13as13bBDmax[i] = s23as23bBDmax[i] = 100.0; // defalut upper value in distance of solvent-solvent, (real distance is in the range of 0.0 to 100.0)
        s1t1cutoffmin[i] = s1t2cutoffmin[i] = s2t1cutoffmin[i] = s2t2cutoffmin[i] = s3t1cutoffmin[i] = s3t2cutoffmin[i] = -1.0; // defalut lower value in distance of solvent-solute, (real distance is in the range of 0.0 to 100.0)
        s1t1cutoffmax[i] = s1t2cutoffmax[i] = s2t1cutoffmax[i] = s2t2cutoffmax[i] = s3t1cutoffmax[i] = s3t2cutoffmax[i] = 100.0; // defalut upper value in distance of solvent-solute, (real distance is in the range of 0.0 to 100.0)
    }


// Start reading Input Files
    {
  if(argc < 3){
    printf("Usuage: %s inputfile coordinate-file(s)\n"
           "inputfile is assumed to have extension .input\n"
           "coordinate-file(s) is(are) assumed to have extension .xyz\n", argv[0]);

    exit(-1);
  }

  sprintf(finput, "%s", argv[1]);

  if((fd = fopen(finput,"r")) == NULL){
    printf("Cannot open file %s\n", finput);
    exit(-1);
  }

// number of solvent types keyword
  if(findf(fd, 4, "[NUMBER", "OF", "SOLVENT", "TYPES]")==EOF){
    printf("Cannot find [NUMBER OF SOLVENT TYPES] keyword\n");
    exit(-1);
  }

// number of solvent types value
  if(fscanf(fd,"%d", &nsolventtype)!=1){
    printf("Cannot read number of solvent types\n");
    exit(-1);
  }

  rewind(fd);

// number of solute types keyword
  if(findf(fd, 4, "[NUMBER", "OF", "SOLUTE", "TYPES]")==EOF){
    printf("Cannot find [NUMBER OF SOLUTE TYPES] keyword\n");
    exit(-1);
  }

// number of solute types value
  if(fscanf(fd,"%d", &nsolutetype)!=1){
    printf("Cannot read number of solute types\n");
    exit(-1);
  }

  rewind(fd);

  if(argc != (nsolventtype + nsolutetype + 2)){
    printf("Usuage: %s inputfile coordinate-file(s)\n"
           "inputfile is assumed to have extension .input\n"
           "coordinate-file(s) is(are) assumed to have extension .xyz\n", argv[0]);

    exit(-1);
  }

// Continue reading input file

  // number of atoms in solvent 1 keyword
  if(findf(fd, 5, "[NUMBER", "OF", "ATOMS", "IN", "SOLVENT1]")==EOF){
    printf("Cannot find [NUMBER OF ATOMS IN SOLVENT1] keyword\n");
    exit(-1);
  }

  // number of atoms in solvent 1 value
  if(fscanf(fd,"%d", &nsolvent1)!=1){
    printf("Cannot read number of atoms in solvent 1\n");
    exit(-1);
  }
     
    slvntatm1 = (char **)malloc(sizeof(char *)*nsolvent1);
    for(i=0; i < nsolvent1; i++) slvntatm1[i] = (char *)malloc(sizeof(char)*32);

    for(i=0; i < nsolvent1; i++){
      if(fscanf(fd, "%31s %d", slvntatm1[i], &slvntatmnum1[i]) != 2){
        printf("Solvent1: Cannot read atom types and their order in xyz file\n");
        exit(-1);
      }
    }

  rewind(fd);


// If there are two solvent types 

if(nsolventtype == 2 || nsolventtype == 3)
{ 
    // number of atoms in solvent 2 keyword
    if(findf(fd, 5, "[NUMBER", "OF", "ATOMS", "IN", "SOLVENT2]")==EOF){
      printf("Cannot find [NUMBER OF ATOMS IN SOLVENT2] keyword\n");
      exit(-1);
    }

    // number of atoms in solvent 2 value
    if(fscanf(fd,"%d", &nsolvent2)!=1){
      printf("Cannot read number of atoms in solvent 2\n");
      exit(-1);
    }

    slvntatm2 = (char **)malloc(sizeof(char *)*nsolvent2);
    for(i=0; i < nsolvent2; i++) slvntatm2[i] = (char *)malloc(sizeof(char)*32);

     for(i=0; i < nsolvent2; i++){
       if(fscanf(fd, "%31s %d", slvntatm2[i], &slvntatmnum2[i]) != 2){
         printf("Solvent2: Cannot read atom types and their order in xyz file\n");
         exit(-1);
       }
     }

    rewind(fd);

}

// If there are three solvent types

if(nsolventtype == 3)
{
    // number of atoms in solvent 3 keyword
    if(findf(fd, 5, "[NUMBER", "OF", "ATOMS", "IN", "SOLVENT3]")==EOF){
      printf("Cannot find [NUMBER OF ATOMS IN SOLVENT3] keyword\n");
      exit(-1);
    }

    // number of atoms in solvent 3 value
    if(fscanf(fd,"%d", &nsolvent3)!=1){
      printf("Cannot read number of atoms in solvent 3\n");
      exit(-1);
    }

    slvntatm3 = (char **)malloc(sizeof(char *)*nsolvent3);
    for(i=0; i < nsolvent3; i++) slvntatm3[i] = (char *)malloc(sizeof(char)*32);

     for(i=0; i < nsolvent3; i++){
       if(fscanf(fd, "%31s %d", slvntatm3[i], &slvntatmnum3[i]) != 2){
         printf("Solvent3: Cannot read atom types and their order in xyz file\n");
         exit(-1);
       }
     }

     rewind(fd);

}


// Read solute information

if(nsolutetype != 0)
{
  // number of atoms in solute 1 keyword
  if(findf(fd, 5, "[NUMBER", "OF", "ATOMS", "IN", "SOLUTE1]")==EOF){
    printf("Cannot find [NUMBER OF ATOMS IN SOLUTE1] keyword\n");
    exit(-1);
  }

  // number of atoms in solute 1 value
  if(fscanf(fd,"%d", &nsolute1)!=1){
    printf("Cannot read number of atoms in solute 1\n");
    exit(-1);
  }

    sltatm1 = (char **)malloc(sizeof(char *)*nsolute1);
    for(i=0; i < nsolute1; i++) sltatm1[i] = (char *)malloc(sizeof(char)*32);

     for(i=0; i < nsolute1; i++){
       if(fscanf(fd, "%31s %d", sltatm1[i], &sltatmnum1[i]) != 2){
         printf("Solute1: Cannot read atom types and their order in xyz file\n");
         exit(-1);
       }
     }

  rewind(fd);

}

// If there are two solute types

if(nsolutetype == 2)
{
  // number of atoms in solute 2 keyword
  if(findf(fd, 5, "[NUMBER", "OF", "ATOMS", "IN", "SOLUTE2]")==EOF){
    printf("Cannot find [NUMBER OF ATOMS IN SOLUTE2] keyword\n");
    exit(-1);
  }

  // number of atoms in solute 2 value
  if(fscanf(fd,"%d", &nsolute2)!=1){
    printf("Cannot read number of atoms in solute 2\n");
    exit(-1);
  }

    sltatm2 = (char **)malloc(sizeof(char *)*nsolute2);
    for(i=0; i < nsolute2; i++) sltatm2[i] = (char *)malloc(sizeof(char)*32);

     for(i=0; i < nsolute2; i++){
       if(fscanf(fd, "%31s %d", sltatm2[i], &sltatmnum2[i]) != 2){
         printf("Solute2: Cannot read atom types and their order in xyz file\n");
         exit(-1);
       }
     }

  rewind(fd);

}

// Read Periodic Boundary Conditions

  if(findf(fd, 3, "[PERIODIC", "BOUNDARY", "CONDITIONS]")==EOF){
    printf("Cannot find [PERIODIC BOUNDARY CONDITIONS] keyword\n");
    exit(-1);
  }

  if(fscanf(fd,"%d", &pbc)!=1){
    printf("Cannot read YES (1) or NO (0) for PERIODIC BOUNDARY CONDITIONS\n");
    exit(-1);
  }

  rewind(fd);

  // Read Box dimensions

      if(findf(fd, 2, "[BOX", "XSIDE]")==EOF){
        printf("Cannot find [BOX XSIDE] keyword\n");
        exit(-1);
      }
      if(fscanf(fd,"%lf", &xside)!=1){
        printf("Cannot read the x dimension of box\n");
        exit(-1);
      }
      rewind(fd);

      if(findf(fd, 2, "[BOX", "YSIDE]")==EOF){
        printf("Cannot find [BOX YSIDE] keyword\n");
        exit(-1);
      }
      if(fscanf(fd,"%lf", &yside)!=1){
        printf("Cannot read the y dimension of box\n");
        exit(-1);
      }
      rewind(fd);

      if(findf(fd, 2, "[BOX", "ZSIDE]")==EOF){
        printf("Cannot find [BOX ZSIDE] keyword\n");
        exit(-1);
      }
      if(fscanf(fd,"%lf", &zside)!=1){
        printf("Cannot read the z dimension of box\n");
        exit(-1);
      }
      rewind(fd);


// Read weighted graph keyword, added, Feb-2018

  if(findf(fd, 2, "[WEIGHTED", "GRAPH]") != EOF){   // find this keyword, read the other informations
    if(fscanf(fd, "%d", &cnwg.index_weighted_graph) != 1){
      printf("Cannot read the value YES(1) or NO(0) for WEIGHTED GRAPH keyword\n");
      exit(-1);
    }
    
    rewind(fd);

    if(cnwg.index_weighted_graph == 1)   // the weighted graph is requested, only works for { 1 solute + 1 or 2 solvent }, Feb-2018
    {
      /* one way to create the weighted graph is based on cluster of solute1 */
      if(findf(fd, 5, "[WEIGHTED", "GRAPH", "BY", "SOLUTE1", "CLUSTER]") == EOF){
        printf("Cannot find the keyword [WEIGHTED GRAPH BY SOLUTE1 CLUSTER]\n");
        exit(-1);
      }

      if(fscanf(fd, "%d", &cnwg.index_wg_st1_cluster) != 1){
        printf("Cannot find the value of keyword [WEIGHTED GRAPH BY SOLUTE1 CLUSTER]\n");
        exit(-1);
      }

      if(cnwg.index_wg_st1_cluster < 0)
      {
        printf("Wrong value at keyword [WEIGHTED GRAPH BY SOLUTE1 CLUSTER]\n");
        exit(-1);
      }
      else if(cnwg.index_wg_st1_cluster > NUM_INTER)
      {
        printf("Size-overflow at keyword [WEIGHTED GRAPH BY SOLUTE1 CLUSTER], default maximum %d\n",NUM_INTER);
        exit(-1);
      }
      else
      {
        for(i=0; i < cnwg.index_wg_st1_cluster; i++){
          if(fscanf(fd, "%d %d %d %d %lf %lf", &cnwg.wg_solute_type, &cnwg.cluster_st1_sv_type[i], &cnwg.cluster_st1_atom[i], &cnwg.cluster_st1_sv_atom[i], 
                    &cnwg.cluster_st1_sv_dist_min[i], &cnwg.cluster_st1_sv_dist_max[i]) != 6){
             printf("Cannot read parameters in keyword [WEIGHTED GRAPH BY SOLUTE1 CLUSTER]\n");
             exit(-1);
          }

          if(cnwg.wg_solute_type != 1){  // it only works for 1 solute, the solute-id must be 1 currently, Feb-2018
            printf("Wrong value in keyword [WEIGHTED GRAPH BY SOLUTE1 CLUSTER]\n Only solute1 is implimented, wait for the next version\n");
            exit(-1);
          }

          if(cnwg.cluster_st1_sv_type[i] < 0  || cnwg.cluster_st1_sv_type[i] >= 3){  // it only works for 1 or 2 solvent, the solvent-id must be 1 or 2, Feb-2018
            printf("Wrong value in keyword [WEIGHTED GRAPH BY SOLUTE1 CLUSTER]\n Only 1 or 2 solvent is implemented, wait for the next version\n");
            exit(-1);
          }
        }

        rewind(fd);

      }  // this is the end of ' if(cnwg.index_wg_st1_cluster < 0) else if.. '

      rewind(fd);

      /* create weighted graphs using other ways, for example, by solvent-solvent, solute-solvent type */
      // .....

      rewind(fd);

      /* below are the edge definition, weighting function, parameters..  */
      if(cnwg.findsolvent(1, cnwg.index_wg_st1_cluster, cnwg.cluster_st1_sv_type) == 1){

        if(findf(fd, 5, "[WEIGHTED", "GRAPH", "EDGE", "SOLUTE1", "SOLVENT1]") == EOF){
          printf("Cannot find the keyword [WEIGHTED GRAPH EDGE SOLUTE1 SOLVENT1]\n");
          exit(-1);
         }

        if(fscanf(fd, "%d", &cnwg.index_wg_st1_sv1) != 1){
          printf("Cannot read the value of keyword [WEIGHTED GRAPH EDGE SOLUTE1 SOLVENT1]\n");
          exit(-1);
        }

        rewind(fd);

        if(cnwg.index_wg_st1_sv1 == 1){   // weighted egde of solute1-solvent1 requested
          if(findf(fd, 6, "[WEIGHTED", "GRAPH", "EDGE", "SOLUTE1", "SOLVENT1", "DISTANCE]") == EOF){
            printf("Cannot find the keyword [WEIGHTED GRAPH EDGE SOLUTE1 SOLVENT1 DISTANCE]\n");
            exit(-1);
          }

          if(fscanf(fd, "%d", &cnwg.num_wg_st1_sv1_dist) != 1){
            printf("Cannot read the value in keyword [WEIGHTED GRAPH EDGE SOLUTE1 SOLVENT1 DISTANCE]\n");
            exit(-1);
          }

          for(i=0; i < cnwg.num_wg_st1_sv1_dist; i++){
            if(fscanf(fd, "%d %d", &cnwg.atom1_wg_st1_sv1[i], &cnwg.atom2_wg_st1_sv1[i]) != 2){
              printf("Cannot read the atom index in [WEIGHTED GRAPH EDGE SOLUTE1 SOLVENT1 DISTANCE]\n");
              exit(-1);
            }
          }

          if(!(nsolventtype >= 1 && nsolutetype >= 1)){
            printf("There is not enough solvent, solute type for [WEIGHTED GRAPH EDGE SOLUTE1 SOLVENT1]\n");
            exit(-1);
          }

          rewind(fd);

          if(findf(fd, 7, "[WEIGHTED", "GRAPH", "EDGE", "SOLUTE1", "SOLVENT1", "DISTANCE", "FUNCTION]") == EOF){
            printf("Cannot find keyword [WEIGHTED GRAPH EDGE SOLUTE1 SOLVENT1 DISTANCE FUNCTION]\n");
            exit(-1);
          }

          if(fscanf(fd, "%d", &cnwg.funct_type_wg_st1_sv1) != 1){
            printf("Cannot read function-type in keyword [WEIGHTED GRAPH EDGE SOLUTE1 SOLVENT1 DISTANCE FUNCTION]\n");
            exit(-1);
          }

          if(cnwg.funct_type_wg_st1_sv1 == 1){  // function-type = 1, which is 1.0 / (1.0 + exp(n0 * (r - r0)/r0)), there are two parameters
            if(fscanf(fd, "%lf %lf", &cnwg.funct_par1_wg_st1_sv1, &cnwg.funct_par2_wg_st1_sv1) != 2){
              printf("Cannot read the parameters in function [WEIGHTED GRAPH EDGE SOLUTE1 SOLVENT1 DISTANCE FUNCTION]\n");
              exit(-1);
            }
          }
          else {  // other function formalism to be implemented here
            printf("function formalism other than 1 is not implemented\n");
            exit(-1);
          }

          rewind(fd);
        } // this is the end of ' if(cnwg.index_wg_st1_sv1 == 1) '

        rewind(fd);
      }

      rewind(fd);

      if(cnwg.findsolvent(2, cnwg.index_wg_st1_cluster, cnwg.cluster_st1_sv_type) == 1) {

        if(findf(fd, 5, "[WEIGHTED", "GRAPH", "EDGE", "SOLUTE1", "SOLVENT2]") == EOF){
          printf("Cannot find the keyword [WEIGHTED GRAPH EDGE SOLUTE1 SOLVENT2]\n");
          exit(-1);
         }

        if(fscanf(fd, "%d", &cnwg.index_wg_st1_sv2) != 1){
          printf("Cannot read the value of keyword [WEIGHTED GRAPH EDGE SOLUTE1 SOLVENT2]\n");
          exit(-1);
        }

        if(cnwg.index_wg_st1_sv2 == 1){   // weighted egde of solute1-solvent2 requested
          if(findf(fd, 6, "[WEIGHTED", "GRAPH", "EDGE", "SOLUTE1", "SOLVENT2", "DISTANCE]") == EOF){
            printf("Cannot find the keyword [WEIGHTED GRAPH EDGE SOLUTE1 SOLVENT2 DISTANCE]\n");
            exit(-1);
          }

          if(fscanf(fd, "%d", &cnwg.num_wg_st1_sv2_dist) != 1){
            printf("Cannot read the value in keyword [WEIGHTED GRAPH EDGE SOLUTE1 SOLVENT2 DISTANCE]\n");
            exit(-1);
          }

          for(i=0; i < cnwg.num_wg_st1_sv2_dist; i++){
            if(fscanf(fd, "%d %d", &cnwg.atom1_wg_st1_sv2[i], &cnwg.atom2_wg_st1_sv2[i]) != 2){
              printf("Cannot read the atom index in [WEIGHTED GRAPH EDGE SOLUTE1 SOLVENT2 DISTANCE]\n");
              exit(-1);
            }
          }

          if(!(nsolventtype >= 2 && nsolutetype >= 1)){
            printf("There is not enough solvent, solute type for [WEIGHTED GRAPH EDGE SOLUTE1 SOLVENT2]\n");
            exit(-1);
          }

          rewind(fd);

          if(findf(fd, 7, "[WEIGHTED", "GRAPH", "EDGE", "SOLUTE1", "SOLVENT2", "DISTANCE", "FUNCTION]") == EOF){
            printf("Cannot find keyword [WEIGHTED GRAPH EDGE SOLUTE1 SOLVENT2 DISTANCE FUNCTION]\n");
            exit(-1);
          }

          if(fscanf(fd, "%d", &cnwg.funct_type_wg_st1_sv2) != 1){
            printf("Cannot read function-type in keyword [WEIGHTED GRAPH EDGE SOLUTE1 SOLVENT2 DISTANCE FUNCTION]\n");
            exit(-1);
          }

          if(cnwg.funct_type_wg_st1_sv2 == 1){  // function-type = 1, which is 1.0 / (1.0 + exp(n0 * (r - r0)/r0)), there are two parameters
            if(fscanf(fd, "%lf %lf", &cnwg.funct_par1_wg_st1_sv2, &cnwg.funct_par2_wg_st1_sv2) != 2){
              printf("Cannot read the parameters in function [WEIGHTED GRAPH EDGE SOLUTE1 SOLVENT2 DISTANCE FUNCTION]\n");
              exit(-1);
            }
          }
         else { // other weight function formalism to be implemented here
            printf("function formalism other than 1 is not implemented\n");
            exit(-1);
          }

          rewind(fd);
        }  // this is the end of ' if(cnwg.index_wg_st1_sv2 == 1) ' 

        rewind(fd);
      }

      rewind(fd);

      if(cnwg.findsolvent(1, cnwg.index_wg_st1_cluster, cnwg.cluster_st1_sv_type) == 1){

        if(findf(fd, 5, "[WEIGHTED", "GRAPH", "EDGE", "SOLVENT1", "SOLVENT1]") == EOF){  // weighted graph for solvent1-solvent1
          printf("Cannot find keyword [WEIGHTED GRAPH EDGE SOLVENT1 SOLVENT1]\n");
          exit(-1);
        }

        if(fscanf(fd, "%d", &cnwg.index_wg_sv1_sv1) != 1){
          printf("Cannot read the value of keyword [WEIGHTED GRAPH EDGE SOLVENT1 SOLVENT1]\n");
          exit(-1);
        }

        rewind(fd);

        if(cnwg.index_wg_sv1_sv1 == 1){  // solvent1-solvent1 weighted graph requested
          if(findf(fd, 6, "[WEIGHTED", "GRAPH", "EDGE", "SOLVENT1", "SOLVENT1", "DISTANCE]") == EOF){
            printf("Cannot find keyword [WEIGHTED GRAPH EDGE SOLVENT1 SOLVENT1 DISTANCE]\n");
            exit(-1);
          }

          if(fscanf(fd, "%d", &cnwg.num_wg_sv1_sv1_dist) != 1){
            printf("Cannot read value in keyword [WEIGHTED GRAPH EDGE SOLVENT1 SOLVENT1 DISTANCE]\n");
            exit(-1);
          }

          for(i=0; i < cnwg.num_wg_sv1_sv1_dist; i++){
            if(fscanf(fd, "%d %d", &cnwg.atom1_wg_sv1_sv1[i], &cnwg.atom2_wg_sv1_sv1[i]) != 2){
              printf("Cannot read atom index in keyword [WEIGHTED GRAPH EDGE SOLVENT1 SOLVENT1 DISTANCE]\n");
              exit(-1);
            }
          }

          rewind(fd);

          if(findf(fd, 7, "[WEIGHTED", "GRAPH", "EDGE", "SOLVENT1", "SOLVENT1", "DISTANCE", "FUNCTION]") == EOF){
            printf("Cannot find keyword [WEIGHTED GRAPH EDGE SOLVENT1 SOLVENT1 DISTANCE FUNCTION]\n");
            exit(-1);
          }

          if(fscanf(fd, "%d", &cnwg.funct_type_wg_sv1_sv1) != 1){
            printf("Cannot read the value in keyword [WEIGHTED GRAPH EDGE SOLVENT1 SOLVENT1 DISTANCE FUNCTION]\n");
            exit(-1);
          }

          if(cnwg.funct_type_wg_sv1_sv1 == 1) {
            if(fscanf(fd, "%lf %lf", &cnwg.funct_par1_wg_sv1_sv1, &cnwg.funct_par2_wg_sv1_sv1) != 2){
              printf("Cannot read the parameters in keyword [WEIGHTED GRAPH EDGE SOLVENT1 SOLVENT1 DISTANCE FUNCTION]\n");
              exit(-1);
            }
          }
          else {
            printf("function formalism other than 1 is not implemented\n");
            exit(-1);
          }

          rewind(fd);
        }  // this is the end of ' if(cnwg.index_wg_sv1_sv1==1) '

        rewind(fd);
      }

      rewind(fd);

      if(cnwg.findsolvent(2, cnwg.index_wg_st1_cluster, cnwg.cluster_st1_sv_type) == 1) {

        if(findf(fd, 5, "[WEIGHTED", "GRAPH", "EDGE", "SOLVENT2", "SOLVENT2]") == EOF){  // weighted graph for solvent2-solvent2
          printf("Cannot find keyword [WEIGHTED GRAPH EDGE SOLVENT2 SOLVENT2]\n");
          exit(-1);
        }

        if(fscanf(fd, "%d", &cnwg.index_wg_sv2_sv2) != 1){
          printf("Cannot read the value of keyword [WEIGHTED GRAPH EDGE SOLVENT2 SOLVENT2]\n");
          exit(-1);
        }

        rewind(fd);

        if(cnwg.index_wg_sv2_sv2 == 1){  // solvent2-solvent2 weighted graph requested
          if(findf(fd, 6, "[WEIGHTED", "GRAPH", "EDGE", "SOLVENT2", "SOLVENT2", "DISTANCE]") == EOF){
            printf("Cannot find keyword [WEIGHTED GRAPH EDGE SOLVENT2 SOLVENT2 DISTANCE]\n");
            exit(-1);
          }

          if(fscanf(fd, "%d", &cnwg.num_wg_sv2_sv2_dist) != 1){
            printf("Cannot read value in keyword [WEIGHTED GRAPH EDGE SOLVENT2 SOLVENT2 DISTANCE]\n");
            exit(-1);
          }

          for(i=0; i < cnwg.num_wg_sv2_sv2_dist; i++){
            if(fscanf(fd, "%d %d", &cnwg.atom1_wg_sv2_sv2[i], &cnwg.atom2_wg_sv2_sv2[i]) != 2){
              printf("Cannot read atom index in keyword [WEIGHTED GRAPH EDGE SOLVENT2 SOLVENT2 DISTANCE]\n");
              exit(-1);
            }
          }

          rewind(fd);

          if(findf(fd, 7, "[WEIGHTED", "GRAPH", "EDGE", "SOLVENT2", "SOLVENT2", "DISTANCE", "FUNCTION]") == EOF){
            printf("Cannot find keyword [WEIGHTED GRAPH EDGE SOLVENT2 SOLVENT2 DISTANCE FUNCTION]\n");
            exit(-1);
          }

          if(fscanf(fd, "%d", &cnwg.funct_type_wg_sv2_sv2) != 1){
            printf("Cannot read the value in keyword [WEIGHTED GRAPH EDGE SOLVENT2 SOLVENT2 DISTANCE FUNCTION]\n");
            exit(-1);
          }

          if(cnwg.funct_type_wg_sv2_sv2 == 1) {
            if(fscanf(fd, "%lf %lf", &cnwg.funct_par1_wg_sv2_sv2, &cnwg.funct_par2_wg_sv2_sv2) != 2){
              printf("Cannot read the parameters in keyword [WEIGHTED GRAPH EDGE SOLVENT2 SOLVENT2 DISTANCE FUNCTION]\n");
              exit(-1);
            }
          }
          else {
            printf("function formalism other than 1 is not implemented\n");
            exit(-1);
          }

          rewind(fd);
        }  // this is the end of ' if(cnwg.index_wg_sv2_sv2==1) '

        rewind(fd);
      }

      rewind(fd);

      if(cnwg.findsolvent(1, cnwg.index_wg_st1_cluster, cnwg.cluster_st1_sv_type) == 1 && cnwg.findsolvent(2, cnwg.index_wg_st1_cluster, cnwg.cluster_st1_sv_type) == 1){

        if(findf(fd, 5, "[WEIGHTED", "GRAPH", "EDGE", "SOLVENT1", "SOLVENT2]") == EOF){  // weighted graph for solvent1-solvent2
          printf("Cannot find keyword [WEIGHTED GRAPH EDGE SOLVENT1 SOLVENT2]\n");
          exit(-1);
        }

        if(fscanf(fd, "%d", &cnwg.index_wg_sv1_sv2) != 1){
          printf("Cannot read the value of keyword [WEIGHTED GRAPH EDGE SOLVENT1 SOLVENT2]\n");
          exit(-1);
        }

        rewind(fd);

        if(cnwg.index_wg_sv1_sv2 == 1){  // solvent1-solvent2 weighted graph requested
          if(findf(fd, 6, "[WEIGHTED", "GRAPH", "EDGE", "SOLVENT1", "SOLVENT2", "DISTANCE]") == EOF){
            printf("Cannot find keyword [WEIGHTED GRAPH EDGE SOLVENT1 SOLVENT2 DISTANCE]\n");
            exit(-1);
          }

          if(fscanf(fd, "%d", &cnwg.num_wg_sv1_sv2_dist) != 1){
            printf("Cannot read value in keyword [WEIGHTED GRAPH EDGE SOLVENT1 SOLVENT2 DISTANCE]\n");
            exit(-1);
          }

          for(i=0; i < cnwg.num_wg_sv1_sv2_dist; i++){
            if(fscanf(fd, "%d %d", &cnwg.atom1_wg_sv1_sv2[i], &cnwg.atom2_wg_sv1_sv2[i]) != 2){
              printf("Cannot read atom index in keyword [WEIGHTED GRAPH EDGE SOLVENT1 SOLVENT2 DISTANCE]\n");
              exit(-1);
            }
          }

          rewind(fd);

          if(findf(fd, 7, "[WEIGHTED", "GRAPH", "EDGE", "SOLVENT1", "SOLVENT2", "DISTANCE", "FUNCTION]") == EOF){
            printf("Cannot find keyword [WEIGHTED GRAPH EDGE SOLVENT1 SOLVENT2 DISTANCE FUNCTION]\n");
            exit(-1);
          }

          if(fscanf(fd, "%d", &cnwg.funct_type_wg_sv1_sv2) != 1){
            printf("Cannot read the value in keyword [WEIGHTED GRAPH EDGE SOLVENT1 SOLVENT2 DISTANCE FUNCTION]\n");
            exit(-1);
          }

          if(cnwg.funct_type_wg_sv1_sv2 == 1) {
            if(fscanf(fd, "%lf %lf", &cnwg.funct_par1_wg_sv1_sv2, &cnwg.funct_par2_wg_sv1_sv2) != 2){
              printf("Cannot read the parameters in keyword [WEIGHTED GRAPH EDGE SOLVENT1 SOLVENT2 DISTANCE FUNCTION]\n");
              exit(-1);
            }
          }
          else {
            printf("function formalism other than 1 is not implemented\n");
            exit(-1);
          }

          rewind(fd);
        }  // this is the end of ' if(cnwg.index_wg_sv1_sv2==1) '

        rewind(fd);
      }

      rewind(fd);
    } // this is the end of ' if(cnwg.index_weighted_graph == 1) ' 

    rewind(fd);

  }  // this is the end of ' if(findf(fd, 2, "[WEIGHTED", "GRAPH]") != EOF) '

  rewind(fd);



// Read Graph Types 

 // GRAPH SOLVENT1 SOLVENT1 keyword 
  if(findf(fd, 3, "[GRAPH", "SOLVENT1", "SOLVENT1]")==EOF){
    printf("Cannot find [GRAPH SOLVENT1 SOLVENT1] keyword\n");
    exit(-1);
  }

 // GRAPH SOLVENT1 SOLVENT1 value  
  if(fscanf(fd,"%d", &gs1s1)!=1){
    printf("Cannot read YES (1) or NO (0) for GRAPH SOLVENT1 SOLVENT1 search\n");
    exit(-1);
  }

  rewind(fd);

 // If GRAPH SOLVENT1 SOLVENT1 requested

  if(gs1s1 == 1)
  {
    if(findf(fd, 4, "[SOLVENT1", "SOLVENT1", "HBOND", "DISTANCE]")==EOF){
      printf("Cannot find [SOLVENT1 SOLVENT1 HBOND DISTANCE] keyword\n");
      exit(-1);
    }

    if(fscanf(fd,"%d", &s1s1hbdn)!=1){
      printf("Cannot read [SOLVENT1 SOLVENT1 HBOND DISTANCE] value\n");
      exit(-1);
    }

    for(i=0; i < s1s1hbdn; i++){
      if(fscanf(fd, "%d %d %lf %lf", &s1a[i], &s1b[i], &s1as1bBDmin[i], &s1as1bBDmax[i]) != 4){   // 2015.12.14 , I changed here, Tiecheng
        printf("Solvent1-Solvent1: Cannot read HBOND DISTANCE pairs and values\n");
        exit(-1);
      }
    }

    rewind(fd);

    if(findf(fd, 4, "[SOLVENT1", "SOLVENT1", "HBOND", "ANGLE]")==EOF){
      printf("Cannot find [SOLVENT1 SOLVENT1 HBOND ANGLE] keyword\n");
      exit(-1);
    }

    if(fscanf(fd,"%d", &s1s1hban)!=1){
      printf("Cannot read [SOLVENT1 SOLVENT1 HBOND ANGLE] value\n");
      exit(-1);
    }

    for(i=0; i < s1s1hban; i++){
      if(fscanf(fd, "%d %d %d %d %d %lf %lf", &s1s1v1[i], &s1s1v2[i], &s1s1v3[i], &s1s1v4[i], &s1s1v5[i], &s1s1v6[i], &s1s1v7[i]) != 7){   // 2015.12.14 , I changed here, Tiecheng
        printf("Solvent1-Solvent1: Cannot read HBOND ANGLE pairs and values\n");
        exit(-1);
      }
    }

    rewind(fd);

    // here I'm adding keywords for energetic definition for water-water, Tiecheng
    E_s1s1_num = 0; E_s1s1_charge_num = 0; E_s1s1_LJ_num = 0;    // initialize
    for(i=0; i < NUM_INTER; i++)
    {
       E_s1s1v1[i] = E_s1s1v2[i] = 1;
       E_s1s1_min[i] = -100.0; E_s1s1_max[i] = 100.0;
       E_s1s1_charge_index[i] = 1;
       E_s1s1_charge_value[i] = 0.0;
       E_s1s1_LJ_index_a[i] = E_s1s1_LJ_index_b[i] = 1;
       E_s1s1_LJ_value_sigma[i] = E_s1s1_LJ_value_epsilon[i] = 0.0;
    }

    if(findf(fd, 4, "[SOLVENT1", "SOLVENT1", "HBOND", "ENERGY]")==EOF){  // this keyword should be optional
      printf("Warning cannot find [SOLVENT1 SOLVENT1 HBOND ENERGY] keyword, no pair-energy definition is used\n");
    }
    else   // pair energy calculation (only for SPC, SPC/E water currently) is requested
    {
      if(fscanf(fd,"%d",&E_s1s1_num)!=1){
        printf("Cannot read [SOLVENT1 SOLVENT1 HBOND ENERGY] value\n");
        exit(-1);
      }

      for(i=0; i < E_s1s1_num; i++){
        if(fscanf(fd, "%d %d %lf %lf", &E_s1s1v1[i], &E_s1s1v2[i], &E_s1s1_min[i], &E_s1s1_max[i]) != 4) {
          printf("Solvent1-Solvent1: Cannot read HBOND ENERGY pairs and values\n");
          exit(-1);
        }
      }

      rewind(fd);
 
      if(findf(fd, 5, "[CHARGES", "OF", "ATOMS", "IN", "SOLVENT1]")==EOF){
        printf("Cannot find [CHARGES OF ATOMS IN SOLVENT1] keyword\n");
        exit(-1);
      }

      if(fscanf(fd,"%d",&E_s1s1_charge_num)!=1){
        printf("Cannot read [CHARGES OF ATOMS IN SOLVENT1] value\n");
        exit(-1);
      }

      for(i=0; i<E_s1s1_charge_num; i++){
        if(fscanf(fd, "%d %lf", &E_s1s1_charge_index[i], &E_s1s1_charge_value[i]) != 2){
           printf("Solvent1-Solvent1: Cannot read [CHARGES OF ATOMS IN SOLVENT1] values\n");
           exit(-1);
        }
      }
 
      rewind(fd);

      if(findf(fd, 6, "[LJ", "PARAMETERS", "OF", "ATOMS", "IN", "SOLVENT1]")==EOF){
        printf("Cannot find [LJ PARAMETERS OF ATOMS IN SOLVENT1] keyword\n");
        exit(-1);
      }

      if(fscanf(fd,"%d",&E_s1s1_LJ_num)!=1){
        printf("Cannot find [LJ PARAMETERS OF ATOMS IN SOLVENT1] value\n");
        exit(-1);
      }

      for(i=0; i<E_s1s1_LJ_num; i++){
        if(fscanf(fd, "%d %d %lf %lf",&E_s1s1_LJ_index_a[i],&E_s1s1_LJ_index_b[i],&E_s1s1_LJ_value_sigma[i],&E_s1s1_LJ_value_epsilon[i]) != 4){
          printf("Solvent1-Solven1: Cannot read [LJ PARAMETERS OF ATOMS IN SOLVENT1] values\n");
          exit(-1);
        }
      }

      rewind(fd);
   }  // this is the end of if-else statement for ENERGY keyword

   rewind(fd);
  }  // this is the end of if(gs1s1==1)


 // GRAPH SOLVENT2 SOLVENT2 keyword 
  if(findf(fd, 3, "[GRAPH", "SOLVENT2", "SOLVENT2]")==EOF){
    printf("Cannot find [GRAPH SOLVENT2 SOLVENT2] keyword\n");
    exit(-1);
  }

 // GRAPH SOLVENT2 SOLVENT2 value   
  if(fscanf(fd,"%d", &gs2s2)!=1){
    printf("Cannot read YES (1) or NO (0) for GRAPH SOLVENT2 SOLVENT2 search\n");
    exit(-1);
  }

  rewind(fd);

 // If GRAPH SOLVENT2 SOLVENT2 requested

  if(gs2s2 == 1)
  {
    if(findf(fd, 4, "[SOLVENT2", "SOLVENT2", "HBOND", "DISTANCE]")==EOF){
      printf("Cannot find [SOLVENT2 SOLVENT2 HBOND DISTANCE] keyword\n");
      exit(-1);
    }

    if(fscanf(fd,"%d", &s2s2hbdn)!=1){
      printf("Cannot read [SOLVENT2 SOLVENT2 HBOND DISTANCE] value\n");
      exit(-1);
    }

    for(i=0; i < s2s2hbdn; i++){
      if(fscanf(fd, "%d %d %lf %lf", &s2a[i], &s2b[i], &s2as2bBDmin[i], &s2as2bBDmax[i]) != 4){  // 2015.12.14 , I changed here, Tiecheng
        printf("Solvent2-Solvent2: Cannot read HBOND DISTANCE pairs and values\n");
        exit(-1);
      }
    }

    rewind(fd);

    if(findf(fd, 4, "[SOLVENT2", "SOLVENT2", "HBOND", "ANGLE]")==EOF){
      printf("Cannot find [SOLVENT2 SOLVENT2 HBOND ANGLE] keyword\n");
      exit(-1);
    }

    if(fscanf(fd,"%d", &s2s2hban)!=1){
      printf("Cannot read [SOLVENT2 SOLVENT2 HBOND ANGLE] value\n");
      exit(-1);
    }

    for(i=0; i < s2s2hban; i++){
      if(fscanf(fd, "%d %d %d %d %d %lf %lf", &s2s2v1[i], &s2s2v2[i], &s2s2v3[i], &s2s2v4[i], &s2s2v5[i], &s2s2v6[i], &s2s2v7[i]) != 7){    // 2015.12.14 , I changed here, Tiecheng
        printf("Solvent2-Solvent2: Cannot read HBOND ANGLE pairs and values\n");
        exit(-1);
      }
    }

    rewind(fd);

  }

 // GRAPH SOLVENT3 SOLVENT3 keyword
  if(findf(fd, 3, "[GRAPH", "SOLVENT3", "SOLVENT3]")==EOF){
    printf("Cannot find [GRAPH SOLVENT3 SOLVENT3] keyword\n");
    exit(-1);
  }

 // GRAPH SOLVENT3 SOLVENT3 value
  if(fscanf(fd,"%d", &gs3s3)!=1){
    printf("Cannot read YES (1) or NO (0) for GRAPH SOLVENT3 SOLVENT3 search\n");
    exit(-1);
  }

  rewind(fd);

 // If GRAPH SOLVENT3 SOLVENT3 requested

  if(gs3s3 == 1)
  {
    if(findf(fd, 4, "[SOLVENT3", "SOLVENT3", "HBOND", "DISTANCE]")==EOF){
      printf("Cannot find [SOLVENT3 SOLVENT3 HBOND DISTANCE] keyword\n");
      exit(-1);
    }

    if(fscanf(fd,"%d", &s3s3hbdn)!=1){
      printf("Cannot read [SOLVENT3 SOLVENT3 HBOND DISTANCE] value\n");
      exit(-1);
    }

    for(i=0; i < s3s3hbdn; i++){
      if(fscanf(fd, "%d %d %lf %lf", &s3a[i], &s3b[i], &s3as3bBDmin[i], &s3as3bBDmax[i]) != 4){   // 2015.12.14 , I changed here, Tiecheng
        printf("Solvent3-Solvent3: Cannot read HBOND DISTANCE pairs and values\n");
        exit(-1);
      }
    }

    rewind(fd);

    if(findf(fd, 4, "[SOLVENT3", "SOLVENT3", "HBOND", "ANGLE]")==EOF){
      printf("Cannot find [SOLVENT3 SOLVENT3 HBOND ANGLE] keyword\n");
      exit(-1);
    }

    if(fscanf(fd,"%d", &s3s3hban)!=1){
      printf("Cannot read [SOLVENT3 SOLVENT3 HBOND ANGLE] value\n");
      exit(-1);
    }

    for(i=0; i < s3s3hban; i++){
      if(fscanf(fd, "%d %d %d %d %d %lf %lf", &s3s3v1[i], &s3s3v2[i], &s3s3v3[i], &s3s3v4[i], &s3s3v5[i], &s3s3v6[i], &s3s3v7[i]) != 7){   // 2015.12.14 , I changed here, Tiecheng
        printf("Solvent3-Solvent3: Cannot read HBOND ANGLE pairs and values\n");
        exit(-1);
      }
    }

    rewind(fd);

  }


 // GRAPH SOLVENT1 SOLVENT2 keyword
  if(findf(fd, 3, "[GRAPH", "SOLVENT1", "SOLVENT2]")==EOF){
    printf("Cannot find [GRAPH SOLVENT1 SOLVENT2] keyword\n");
    exit(-1);
  }

 // GRAPH SOLVENT1 SOLVENT2 value
  if(fscanf(fd,"%d", &gs1s2)!=1){
    printf("Cannot read YES (1) or NO (0) for GRAPH SOLVENT1 SOLVENT2 search\n");
    exit(-1);
  }

  rewind(fd);

 // If GRAPH SOLVENT1 SOLVENT2 requested

  if(gs1s2 == 1)
  {
    if(findf(fd, 4, "[SOLVENT1", "SOLVENT2", "HBOND", "DISTANCE]")==EOF){
      printf("Cannot find [SOLVENT1 SOLVENT2 HBOND DISTANCE] keyword\n");
      exit(-1);
    }

    if(fscanf(fd,"%d", &s1s2hbdn)!=1){
      printf("Cannot read [SOLVENT1 SOLVENT2 HBOND DISTANCE] value\n");
      exit(-1);
    }

    for(i=0; i < s1s2hbdn; i++){
      if(fscanf(fd, "%d %d %lf %lf", &s12a[i], &s12b[i], &s12as12bBDmin[i], &s12as12bBDmax[i]) != 4){   // 2015.12.14 , I changed here, Tiecheng
        printf("Solvent1-Solvent2: Cannot read HBOND DISTANCE pairs and values\n");
        exit(-1);
      }
    }

    rewind(fd);

    if(findf(fd, 4, "[SOLVENT1", "SOLVENT2", "HBOND", "ANGLE]")==EOF){
      printf("Cannot find [SOLVENT1 SOLVENT2 HBOND ANGLE] keyword\n");
      exit(-1);
    }

    if(fscanf(fd,"%d", &s1s2hban)!=1){
      printf("Cannot read [SOLVENT1 SOLVENT2 HBOND ANGLE] value\n");
      exit(-1);
    }

    for(i=0; i < s1s2hban; i++){
      if(fscanf(fd, "%d %d %d %d %d %lf %lf", &s1s2v1[i], &s1s2v2[i], &s1s2v3[i], &s1s2v4[i], &s1s2v5[i], &s1s2v6[i], &s1s2v7[i]) != 7){  // 2015.12.14 , I changed here, Tiecheng
        printf("Solvent1-Solvent2: Cannot read HBOND ANGLE pairs and values\n");
        exit(-1);
      }
    }

    rewind(fd);

  }


 // GRAPH SOLVENT1 SOLVENT3 keyword
  if(findf(fd, 3, "[GRAPH", "SOLVENT1", "SOLVENT3]")==EOF){
    printf("Cannot find [GRAPH SOLVENT1 SOLVENT3] keyword\n");
    exit(-1);
  }

 // GRAPH SOLVENT1 SOLVENT3 value
  if(fscanf(fd,"%d", &gs1s3)!=1){
    printf("Cannot read YES (1) or NO (0) for GRAPH SOLVENT1 SOLVENT3 search\n");
    exit(-1);
  }

  rewind(fd);

 // If GRAPH SOLVENT1 SOLVENT3 requested

  if(gs1s3 == 1)
  {
    if(findf(fd, 4, "[SOLVENT1", "SOLVENT3", "HBOND", "DISTANCE]")==EOF){
      printf("Cannot find [SOLVENT1 SOLVENT3 HBOND DISTANCE] keyword\n");
      exit(-1);
    }

    if(fscanf(fd,"%d", &s1s3hbdn)!=1){
      printf("Cannot read [SOLVENT1 SOLVENT3 HBOND DISTANCE] value\n");
      exit(-1);
    }

    for(i=0; i < s1s3hbdn; i++){
      if(fscanf(fd, "%d %d %lf %lf", &s13a[i], &s13b[i], &s13as13bBDmin[i], &s13as13bBDmax[i]) != 4){   // 2015.12.14 , I changed here, Tiecheng
        printf("Solvent1-Solvent3: Cannot read HBOND DISTANCE pairs and values\n");
        exit(-1);
      }
    }

    rewind(fd);

    if(findf(fd, 4, "[SOLVENT1", "SOLVENT3", "HBOND", "ANGLE]")==EOF){
      printf("Cannot find [SOLVENT1 SOLVENT3 HBOND ANGLE] keyword\n");
      exit(-1);
    }

    if(fscanf(fd,"%d", &s1s3hban)!=1){
      printf("Cannot read [SOLVENT1 SOLVENT3 HBOND ANGLE] value\n");
      exit(-1);
    }

    for(i=0; i < s1s3hban; i++){
      if(fscanf(fd, "%d %d %d %d %d %lf %lf", &s1s3v1[i], &s1s3v2[i], &s1s3v3[i], &s1s3v4[i], &s1s3v5[i], &s1s3v6[i], &s1s3v7[i]) != 7){   // 2015.12.14 , I changed here, Tiecheng
        printf("Solvent1-Solvent3: Cannot read HBOND ANGLE pairs and values\n");
        exit(-1);
      }
    }

    rewind(fd);

  }


 // GRAPH SOLVENT1 SOLUTE1 keyword
  if(findf(fd, 3, "[GRAPH", "SOLVENT1", "SOLUTE1]")==EOF){
    printf("Cannot find [GRAPH SOLVENT1 SOLUTE1] keyword\n");
    exit(-1);
  }

 // GRAPH SOLVENT1 SOLUTE1 value
  if(fscanf(fd,"%d", &gs1t1)!=1){
    printf("Cannot read YES (1) or NO (0) for GRAPH SOLVENT1 SOLUTE1 search\n");
    exit(-1);
  }

  rewind(fd);

// If GRAPH SOLVENT1 SOLUTE1 requested

  if(gs1t1 == 1)
  {
    if(findf(fd, 3, "[SOLVENT1", "SOLUTE1", "CUTOFF]")==EOF){
      printf("Cannot find [SOLVENT1 SOLUTE1 CUTOFF] keyword\n");
      exit(-1);
    }

    if(fscanf(fd,"%d", &s1t1cutoffnum)!=1){
      printf("Cannot read [SOLVENT1 SOLUTE1 CUTOFF] value\n");
      exit(-1);
    }    
    
    for(i=0; i < s1t1cutoffnum; i++){
      if(fscanf(fd, "%d %d %lf %lf", &s1t1a[i], &s1t1b[i], &s1t1cutoffmin[i], &s1t1cutoffmax[i]) != 4){   // 2015.12.14 , I changed here, Tiecheng
        printf("Solvent1-Solute1: Cannot read atom pairs and cutoff distances\n");
        exit(-1);
      }
    } 

    rewind(fd);
  }

 // GRAPH SOLVENT1 SOLUTE2 keyword
  if(findf(fd, 3, "[GRAPH", "SOLVENT1", "SOLUTE2]")==EOF){
    printf("Cannot find [GRAPH SOLVENT1 SOLUTE2] keyword\n");
    exit(-1);
  }

 // GRAPH SOLVENT1 SOLUTE2 value
  if(fscanf(fd,"%d", &gs1t2)!=1){
    printf("Cannot read YES (1) or NO (0) for GRAPH SOLVENT1 SOLUTE2 search\n");
    exit(-1);
  }

  rewind(fd);

// If GRAPH SOLVENT1 SOLUTE2 requested

  if(gs1t2 == 1)
  {
    if(findf(fd, 3, "[SOLVENT1", "SOLUTE2", "CUTOFF]")==EOF){
      printf("Cannot find [SOLVENT1 SOLUTE2 CUTOFF] keyword\n");
      exit(-1);
    }
 
    if(fscanf(fd,"%d", &s1t2cutoffnum)!=1){
      printf("Cannot read [SOLVENT1 SOLUTE2 CUTOFF] value\n");
      exit(-1);
    }

    for(i=0; i < s1t2cutoffnum; i++){
      if(fscanf(fd, "%d %d %lf %lf", &s1t2a[i], &s1t2b[i], &s1t2cutoffmin[i], &s1t2cutoffmax[i]) != 4){   // 2015.12.14 , I changed here, Tiecheng
        printf("Solvent1-Solute2: Cannot read atom pairs and cutoff distances\n");
        exit(-1);
      }
    }

    rewind(fd);
  }


 // GRAPH SOLVENT2 SOLVENT3 keyword
  if(findf(fd, 3, "[GRAPH", "SOLVENT2", "SOLVENT3]")==EOF){
    printf("Cannot find [GRAPH SOLVENT2 SOLVENT3] keyword\n");
    exit(-1);
  }

 // GRAPH SOLVENT2 SOLVENT3 value
  if(fscanf(fd,"%d", &gs2s3)!=1){
    printf("Cannot read YES (1) or NO (0) for GRAPH SOLVENT2 SOLVENT3 search\n");
    exit(-1);
  }

  rewind(fd);

 // If GRAPH SOLVENT2 SOLVENT3 requested

  if(gs2s3 == 1)
  {
    if(findf(fd, 4, "[SOLVENT2", "SOLVENT3", "HBOND", "DISTANCE]")==EOF){
      printf("Cannot find [SOLVENT2 SOLVENT3 HBOND DISTANCE] keyword\n");
      exit(-1);
    }

    if(fscanf(fd,"%d", &s2s3hbdn)!=1){
      printf("Cannot read [SOLVENT2 SOLVENT3 HBOND DISTANCE] value\n");
      exit(-1);
    }

    for(i=0; i < s2s3hbdn; i++){
      if(fscanf(fd, "%d %d %lf %lf", &s23a[i], &s23b[i], &s23as23bBDmin[i], &s23as23bBDmax[i]) != 4){   // 2015.12.14 , I changed here, Tiecheng
        printf("Solvent2-Solvent3: Cannot read HBOND DISTANCE pairs and values\n");
        exit(-1);
      }
    }

    rewind(fd);

    if(findf(fd, 4, "[SOLVENT2", "SOLVENT3", "HBOND", "ANGLE]")==EOF){
      printf("Cannot find [SOLVENT2 SOLVENT3 HBOND ANGLE] keyword\n");
      exit(-1);
    }

    if(fscanf(fd,"%d", &s2s3hban)!=1){
      printf("Cannot read [SOLVENT2 SOLVENT3 HBOND ANGLE] value\n");
      exit(-1);
    }

    for(i=0; i < s2s3hban; i++){
      if(fscanf(fd, "%d %d %d %d %d %lf %lf", &s2s3v1[i], &s2s3v2[i], &s2s3v3[i], &s2s3v4[i], &s2s3v5[i], &s2s3v6[i], &s2s3v7[i]) != 7){   // 2015.12.14 , I changed here, Tiecheng
        printf("Solvent2-Solvent3: Cannot read HBOND ANGLE pairs and values\n");
        exit(-1);
      }
    }

    rewind(fd);

  }


 // GRAPH SOLVENT2 SOLUTE1 keyword
  if(findf(fd, 3, "[GRAPH", "SOLVENT2", "SOLUTE1]")==EOF){
    printf("Cannot find [GRAPH SOLVENT2 SOLUTE1] keyword\n");
    exit(-1);
  }

 // GRAPH SOLVENT2 SOLUTE1 value
  if(fscanf(fd,"%d", &gs2t1)!=1){
    printf("Cannot read YES (1) or NO (0) for GRAPH SOLVENT2 SOLUTE1 search\n");
    exit(-1);
  }

  rewind(fd);

// If GRAPH SOLVENT2 SOLUTE1 requested

  if(gs2t1 == 1)
  {
    if(findf(fd, 3, "[SOLVENT2", "SOLUTE1", "CUTOFF]")==EOF){
      printf("Cannot find [SOLVENT2 SOLUTE1 CUTOFF] keyword\n");
      exit(-1);
    }
 
    if(fscanf(fd,"%d", &s2t1cutoffnum)!=1){
      printf("Cannot read [SOLVENT2 SOLUTE1 CUTOFF] value\n");
      exit(-1);
    }

    for(i=0; i < s2t1cutoffnum; i++){
      if(fscanf(fd, "%d %d %lf %lf", &s2t1a[i], &s2t1b[i], &s2t1cutoffmin[i], &s2t1cutoffmax[i]) != 4){  // 2015.12.14 , I changed here, Tiecheng
        printf("Solvent2-Solute1: Cannot read atom pairs and cutoff distances\n");
        exit(-1);
      }
    }

    rewind(fd);
  }


 // GRAPH SOLVENT2 SOLUTE2 keyword
  if(findf(fd, 3, "[GRAPH", "SOLVENT2", "SOLUTE2]")==EOF){
    printf("Cannot find [GRAPH SOLVENT2 SOLUTE2] keyword\n");
    exit(-1);
  }

 // GRAPH SOLVENT2 SOLUTE2 value
  if(fscanf(fd,"%d", &gs2t2)!=1){
    printf("Cannot read YES (1) or NO (0) for GRAPH SOLVENT2 SOLUTE2 search\n");
    exit(-1);
  }

  rewind(fd);

// If GRAPH SOLVENT2 SOLUTE2 requested

  if(gs2t2 == 1)
  {
    if(findf(fd, 3, "[SOLVENT2", "SOLUTE2", "CUTOFF]")==EOF){
      printf("Cannot find [SOLVENT2 SOLUTE2 CUTOFF] keyword\n");
      exit(-1);
    }
 
    if(fscanf(fd,"%d", &s2t2cutoffnum)!=1){
      printf("Cannot read [SOLVENT2 SOLUTE2 CUTOFF] value\n");
      exit(-1);
    }

    for(i=0; i < s2t2cutoffnum; i++){
      if(fscanf(fd, "%d %d %lf %lf", &s2t2a[i], &s2t2b[i], &s2t2cutoffmin[i], &s2t2cutoffmax[i]) != 4){  // 2015.12.14 , I changed here, Tiecheng
        printf("Solvent2-Solute2: Cannot read atom pairs and cutoff distances\n");
        exit(-1);
      }
    }

    rewind(fd);
  }


 // GRAPH SOLVENT3 SOLUTE1 keyword
  if(findf(fd, 3, "[GRAPH", "SOLVENT3", "SOLUTE1]")==EOF){
    printf("Cannot find [GRAPH SOLVENT3 SOLUTE1] keyword\n");
    exit(-1);
  }

 // GRAPH SOLVENT3 SOLUTE1 value
  if(fscanf(fd,"%d", &gs3t1)!=1){
    printf("Cannot read YES (1) or NO (0) for GRAPH SOLVENT3 SOLUTE1 search\n");
    exit(-1);
  }

  rewind(fd);

// If GRAPH SOLVENT3 SOLUTE1 requested

  if(gs3t1 == 1)
  {
    if(findf(fd, 3, "[SOLVENT3", "SOLUTE1", "CUTOFF]")==EOF){
      printf("Cannot find [SOLVENT3 SOLUTE1 CUTOFF] keyword\n");
      exit(-1);
    }
 
    if(fscanf(fd,"%d", &s3t1cutoffnum)!=1){
      printf("Cannot read [SOLVENT3 SOLUTE1 CUTOFF] value\n");
      exit(-1);
    }

    for(i=0; i < s3t1cutoffnum; i++){
      if(fscanf(fd, "%d %d %lf %lf", &s3t1a[i], &s3t1b[i], &s3t1cutoffmin[i], &s3t1cutoffmax[i]) != 4){  // 2015.12.14 , I changed here, Tiecheng
        printf("Solvent3-Solute1: Cannot read atom pairs and cutoff distances\n");
        exit(-1);
      }
    }

    rewind(fd);
  }


 // GRAPH SOLVENT3 SOLUTE2 keyword
  if(findf(fd, 3, "[GRAPH", "SOLVENT3", "SOLUTE2]")==EOF){
    printf("Cannot find [GRAPH SOLVENT3 SOLUTE2] keyword\n");
    exit(-1);
  }

 // GRAPH SOLVENT3 SOLUTE2 value
  if(fscanf(fd,"%d", &gs3t2)!=1){
    printf("Cannot read YES (1) or NO (0) for GRAPH SOLVENT3 SOLUTE2 search\n");
    exit(-1);
  }

  rewind(fd);

// If GRAPH SOLVENT3 SOLUTE2 requested

  if(gs3t2 == 1)
  {
    if(findf(fd, 3, "[SOLVENT3", "SOLUTE2", "CUTOFF]")==EOF){
      printf("Cannot find [SOLVENT3 SOLUTE2 CUTOFF] keyword\n");
      exit(-1);
    }
 
    if(fscanf(fd,"%d", &s3t2cutoffnum)!=1){
      printf("Cannot read [SOLVENT3 SOLUTE2 CUTOFF] value\n");
      exit(-1);
    }

    for(i=0; i < s3t2cutoffnum; i++){
      if(fscanf(fd, "%d %d %lf %lf", &s3t2a[i], &s3t2b[i], &s3t2cutoffmin[i], &s3t2cutoffmax[i]) != 4){  // 2015.12.14 , I changed here, Tiecheng
        printf("Solvent3-Solute2: Cannot read atom pairs and cutoff distances\n");
        exit(-1);
      }
    }

    rewind(fd);
  }

 // GRAPH SOLVENT1 SOLVENT2 SOLVENT3 keyword
  if(findf(fd, 4, "[GRAPH", "SOLVENT1", "SOLVENT2", "SOLVENT3]")==EOF){
    printf("Cannot find [GRAPH SOLVENT1 SOLVENT2 SOLVENT3] keyword\n");
    exit(-1);
  }

 // GRAPH SOLVENT1 SOLVENT2 SOLVENT3 value
  if(fscanf(fd,"%d", &gs1s2s3)!=1){
    printf("Cannot read YES (1) or NO (0) for GRAPH SOLVENT1 SOLVENT2 SOLVENT3 search\n");
    exit(-1);
  }

// 2016.July the SOLVENT1 SOLVENT2 SOLVENT3 part has not been included



  rewind(fd);


// Read if print number of nodes requested

 // PRINT NUMBER OF NODES keyword
  if(findf(fd, 4, "[PRINT", "NUMBER", "OF", "NODES]")==EOF){
    printf("Cannot find [PRINT NUMBER OF NODES] keyword\n");
    exit(-1);
  }

 // PRINT NUMBER OF NODES value
  if(fscanf(fd,"%d", &pnumnodes)!=1){
    printf("Cannot read YES (1) or NO (0) for PRINT NUMBER OF NODES keyword\n");
    exit(-1);
  }

  rewind(fd);

// Read GEODESICS keywords

// GEODESICS GD keyword

  if(findf(fd, 2, "[GEODESICS", "GD]")==EOF){
    printf("Cannot find [GEODESICS GD] keyword\n");
    exit(-1);
  }

// GEODESICS GD value

  if(fscanf(fd,"%d", &gdist)!=1){
    printf("Cannot read YES (1) or NO (0) for GEODESICS GD keyword\n");
    exit(-1);
  }

  rewind(fd);

if(gdist == 1)
{

   if(findf(fd, 2, "[GD", "SOLVENT1]")==EOF){
    printf("Cannot find [GD SOLVENT1] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &gdistS1)!=1){
    printf("Cannot read YES (1) or NO (0) for GD SOLVENT1 keyword\n");
    exit(-1);
  } rewind(fd);

  if(gdistS1 == 1){ if(gs1s1 != 1){ printf("GD SOLVENT1: Geodesics require the Graph constructed\n"); exit(-1); } }

   if(findf(fd, 4, "[GD", "SOLVENT1", "EUCLIDEAN", "DISTANCE]")==EOF){
    printf("Cannot find [GD SOLVENT1 EUCLIDEAN DISTANCE] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &EucDistS1)!=1){
    printf("Cannot read YES (1) or NO (0) for GD SOLVENT1 EUCLIDEAN DISTANCE keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 4, "[GD", "SOLVENT1", "EUCLIDEAN", "REFERENCE]")==EOF){
    printf("Cannot find [GD SOLVENT1 EUCLIDEAN REFERENCE] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &EucRefS1)!=1){
    printf("Cannot read parameter for GD SOLVENT1 EUCLIDEAN REFERENCE keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 2, "[GD", "SOLVENT2]")==EOF){
    printf("Cannot find [GD SOLVENT2] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &gdistS2)!=1){
    printf("Cannot read YES (1) or NO (0) for GD SOLVENT2 keyword\n");
    exit(-1);
  } rewind(fd);

  if(gdistS2 == 1){ if(gs2s2 != 1){ printf("GD SOLVENT2: Geodesics require the Graph constructed\n"); exit(-1); } }

   if(findf(fd, 4, "[GD", "SOLVENT2", "EUCLIDEAN", "DISTANCE]")==EOF){
    printf("Cannot find [GD SOLVENT2 EUCLIDEAN DISTANCE] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &EucDistS2)!=1){
    printf("Cannot read YES (1) or NO (0) for GD SOLVENT2 EUCLIDEAN DISTANCE keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 4, "[GD", "SOLVENT2", "EUCLIDEAN", "REFERENCE]")==EOF){
    printf("Cannot find [GD SOLVENT2 EUCLIDEAN REFERENCE] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &EucRefS2)!=1){
    printf("Cannot read parameter for GD SOLVENT2 EUCLIDEAN REFERENCE keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 2, "[GD", "SOLVENT3]")==EOF){
    printf("Cannot find [GD SOLVENT3] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &gdistS3)!=1){
    printf("Cannot read YES (1) or NO (0) for GD SOLVENT3 keyword\n");
    exit(-1);
  } rewind(fd);

  if(gdistS3 == 1){ if(gs3s3 != 1){ printf("GD SOLVENT3: Geodesics require the Graph constructed\n"); exit(-1); } }

   if(findf(fd, 4, "[GD", "SOLVENT3", "EUCLIDEAN", "DISTANCE]")==EOF){
    printf("Cannot find [GD SOLVENT3 EUCLIDEAN DISTANCE] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &EucDistS3)!=1){
    printf("Cannot read YES (1) or NO (0) for GD SOLVENT3 EUCLIDEAN DISTANCE keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 4, "[GD", "SOLVENT3", "EUCLIDEAN", "REFERENCE]")==EOF){
    printf("Cannot find [GD SOLVENT3 EUCLIDEAN REFERENCE] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &EucRefS3)!=1){
    printf("Cannot read parameter for GD SOLVENT3 EUCLIDEAN REFERENCE keyword\n");
    exit(-1);
  } rewind(fd);

// binaries

   if(findf(fd, 3, "[GD", "SOLVENT1", "SOLVENT2]")==EOF){
    printf("Cannot find [GD SOLVENT1 SOLVENT2] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &gdistS1S2)!=1){
    printf("Cannot read YES (1) or NO (0) for GD SOLVENT1 SOLVENT2 keyword\n");
    exit(-1);
  } rewind(fd);

  if(gdistS1S2 == 1){ if(gs1s2 != 1){ printf("GD SOLVENT1 SOLVENT2: Geodesics require the Graph constructed\n"); exit(-1); } }

   if(findf(fd, 5, "[GD", "SOLVENT1", "SOLVENT2", "EUCLIDEAN", "DISTANCE]")==EOF){
    printf("Cannot find [GD SOLVENT1 SOLVENT2 EUCLIDEAN DISTANCE] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &EucDistS1S2)!=1){
    printf("Cannot read YES (1) or NO (0) for GD SOLVENT1 SOLVENT2 EUCLIDEAN DISTANCE keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 5, "[GD", "SOLVENT1", "SOLVENT2", "EUCLIDEAN", "REFERENCE]")==EOF){
    printf("Cannot find [GD SOLVENT1 SOLVENT2 EUCLIDEAN REFERENCE] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d %d", &EucRefS1S2s1, &EucRefS1S2s2)!=2){
    printf("Cannot read parameter for GD SOLVENT1 SOLVENT2 EUCLIDEAN REFERENCE keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 3, "[GD", "SOLVENT1", "SOLVENT3]")==EOF){
    printf("Cannot find [GD SOLVENT1 SOLVENT3] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &gdistS1S3)!=1){
    printf("Cannot read YES (1) or NO (0) for GD SOLVENT1 SOLVENT3 keyword\n");
    exit(-1);
  } rewind(fd);

  if(gdistS1S3 == 1){ if(gs1s3 != 1){ printf("GD SOLVENT1 SOLVENT3: Geodesics require the Graph constructed\n"); exit(-1); } }
  
   if(findf(fd, 5, "[GD", "SOLVENT1", "SOLVENT3", "EUCLIDEAN", "DISTANCE]")==EOF){
    printf("Cannot find [GD SOLVENT1 SOLVENT3 EUCLIDEAN DISTANCE] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &EucDistS1S3)!=1){
    printf("Cannot read YES (1) or NO (0) for GD SOLVENT1 SOLVENT3 EUCLIDEAN DISTANCE keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 5, "[GD", "SOLVENT1", "SOLVENT3", "EUCLIDEAN", "REFERENCE]")==EOF){
    printf("Cannot find [GD SOLVENT1 SOLVENT3 EUCLIDEAN REFERENCE] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d %d", &EucRefS1S3s1, &EucRefS1S3s3)!=2){
    printf("Cannot read parameter for GD SOLVENT1 SOLVENT3 EUCLIDEAN REFERENCE keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 3, "[GD", "SOLVENT2", "SOLVENT3]")==EOF){
    printf("Cannot find [GD SOLVENT2 SOLVENT3] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &gdistS2S3)!=1){
    printf("Cannot read YES (1) or NO (0) for GD SOLVENT2 SOLVENT3 keyword\n");
    exit(-1);
  } rewind(fd);

  if(gdistS2S3 == 1){ if(gs2s3 != 1){ printf("GD SOLVENT2 SOLVENT3: Geodesics require the Graph constructed\n"); exit(-1); } }

   if(findf(fd, 5, "[GD", "SOLVENT2", "SOLVENT3", "EUCLIDEAN", "DISTANCE]")==EOF){
    printf("Cannot find [GD SOLVENT2 SOLVENT3 EUCLIDEAN DISTANCE] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &EucDistS2S3)!=1){
    printf("Cannot read YES (1) or NO (0) for GD SOLVENT2 SOLVENT3 EUCLIDEAN DISTANCE keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 5, "[GD", "SOLVENT2", "SOLVENT3", "EUCLIDEAN", "REFERENCE]")==EOF){
    printf("Cannot find [GD SOLVENT2 SOLVENT3 EUCLIDEAN REFERENCE] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d %d", &EucRefS2S3s2, &EucRefS2S3s3)!=2){
    printf("Cannot read parameter for GD SOLVENT2 SOLVENT3 EUCLIDEAN REFERENCE keyword\n");
    exit(-1);
  } rewind(fd);

// ternary

   if(findf(fd, 4, "[GD", "SOLVENT1", "SOLVENT2", "SOLVENT3]")==EOF){
    printf("Cannot find [GD SOLVENT1 SOLVENT2 SOLVENT3] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &gdistS1S2S3)!=1){
    printf("Cannot read YES (1) or NO (0) for GD SOLVENT1 SOLVENT2 SOLVENT3 keyword\n");
    exit(-1);
  } rewind(fd);

  if(gdistS1S2S3 == 1){ if(gs1s2s3 != 1){ printf("GD SOLVENT1 SOLVENT2 SOLVENT3: Geodesics require the Graph constructed\n"); exit(-1); } }

   if(findf(fd, 6, "[GD", "SOLVENT1", "SOLVENT2", "SOLVENT3", "EUCLIDEAN", "DISTANCE]")==EOF){
    printf("Cannot find [GD SOLVENT1 SOLVENT2 SOLVENT3 EUCLIDEAN DISTANCE] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &EucDistS1S2S3)!=1){
    printf("Cannot read YES (1) or NO (0) for GD SOLVENT1 SOLVENT2 SOLVENT3 EUCLIDEAN DISTANCE keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 6, "[GD", "SOLVENT1", "SOLVENT2", "SOLVENT3", "EUCLIDEAN", "REFERENCE]")==EOF){
    printf("Cannot find [GD SOLVENT1 SOLVENT2 SOLVENT3 EUCLIDEAN REFERENCE] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d %d %d", &EucRefS1S2S3s1, &EucRefS1S2S3s2, &EucRefS1S2S3s3)!=3){
    printf("Cannot read parameter for GD SOLVENT1 SOLVENT2 SOLVENT3 EUCLIDEAN REFERENCE keyword\n");
    exit(-1);
  } rewind(fd);


}

// Read if SOLVENT WATER DIPOLE ORIENTATIONS requested

 // SOLVENT WATER DIPOLE ORIENTATIONS keyword

  if(findf(fd, 4, "[SOLVENT", "WATER", "DIPOLE", "ORIENTATIONS]")==EOF){
    printf("Cannot find [SOLVENT WATER DIPOLE ORIENTATIONS] keyword\n");
    exit(-1);
  }

 // SOLVENT WATER DIPOLE ORIENTATIONS value
  if(fscanf(fd,"%d", &swatdip)!=1){
    printf("Cannot read YES (1) or NO (0) for SOLVENT WATER DIPOLE ORIENTATIONS keyword\n");
    exit(-1);
  }

  if(swatdip == 1){
     if(fscanf(fd,"%d %d", &solid, &watid)!=2){
       printf("Cannot read the order of SOLVENT and WATER within solvent(s) for SOLVENT WATER DIPOLE ORIENTATIONS keyword\n");
       exit(-1);
     }
  }

  rewind(fd);

 if((solid == 2 && watid == 1 && swatdip == 1) || (solid == 1 && watid == 2 && swatdip == 1))
    if(gs1s2 != 1){ printf("DIPOLE: [GRAPH SOLVENT1 SOLVENT2] must be turned on\n"); exit(-1); }

 if((solid == 3 && watid == 1 && swatdip == 1) || (solid == 1 && watid == 3 && swatdip == 1))
    if(gs1s3 != 1){ printf("DIPOLE: [GRAPH SOLVENT1 SOLVENT3] must be turned on\n"); exit(-1); }

 if((solid == 2 && watid == 3 && swatdip == 1) || (solid == 3 && watid == 2 && swatdip == 1))
    if(gs2s3 != 1){ printf("DIPOLE: [GRAPH SOLVENT2 SOLVENT3] must be turned on\n"); exit(-1); }


// Read if SOLUTE WATER DIPOLE ORIENTATIONS requested

 // SOLUTE1 WATER DIPOLE ORIENTATIONS keyword
  if(findf(fd, 4, "[SOLUTE1", "WATER", "DIPOLE", "ORIENTATIONS]")==EOF){
    printf("Cannot find [SOLUTE1 WATER DIPOLE ORIENTATIONS] keyword\n");
    exit(-1);
  }

 // SOLUTE1 WATER DIPOLE ORIENTATIONS value
  if(fscanf(fd,"%d", &t1watdip)!=1){
    printf("Cannot read YES (1) or NO (0) for SOLUTE1 WATER DIPOLE ORIENTATIONS keyword\n");
    exit(-1);
  }

  if(t1watdip == 1){

     if(fscanf(fd,"%d", &watid)!=1){
       printf("Cannot read the order of WATER within solvent(s) for SOLUTE1 WATER DIPOLE ORIENTATIONS keyword\n");
       exit(-1);
     }

 }

  rewind(fd);

 if(watid == 1 && t1watdip == 1)
    if(gs1t1 != 1){ printf("DIPOLE: [GRAPH SOLVENT1 SOLUTE1] must be turned on\n"); exit(-1); }

 if(watid == 2 && t1watdip == 1)
    if(gs2t1 != 1){ printf("DIPOLE: [GRAPH SOLVENT2 SOLUTE1] must be turned on\n"); exit(-1); }

 if(watid == 3 && t1watdip == 1)
    if(gs3t1 != 1){ printf("DIPOLE: [GRAPH SOLVENT3 SOLUTE1] must be turned on\n"); exit(-1); }

 // SOLUTE2 WATER DIPOLE ORIENTATIONS keyword
  if(findf(fd, 4, "[SOLUTE2", "WATER", "DIPOLE", "ORIENTATIONS]")==EOF){
    printf("Cannot find [SOLUTE2 WATER DIPOLE ORIENTATIONS] keyword\n");
    exit(-1);
  }

 // SOLUTE2 WATER DIPOLE ORIENTATIONS value
  if(fscanf(fd,"%d", &t2watdip)!=1){
    printf("Cannot read YES (1) or NO (0) for SOLUTE2 WATER DIPOLE ORIENTATIONS keyword\n");
    exit(-1);
  }

  if(t2watdip == 1){

     if(fscanf(fd,"%d", &watid)!=1){
       printf("Cannot read the order of WATER within solvent(s) for SOLUTE2 WATER DIPOLE ORIENTATIONS keyword\n");
       exit(-1);
     }

 }

  rewind(fd);

 if(watid == 1 && t2watdip == 1)
    if(gs1t2 != 1){ printf("DIPOLE: [GRAPH SOLVENT1 SOLUTE2] must be turned on\n"); exit(-1); }

 if(watid == 2 && t2watdip == 1)
    if(gs2t2 != 1){ printf("DIPOLE: [GRAPH SOLVENT2 SOLUTE2] must be turned on\n"); exit(-1); }

 if(watid == 3 && t2watdip == 1)
    if(gs3t2 != 1){ printf("DIPOLE: [GRAPH SOLVENT3 SOLUTE2] must be turned on\n"); exit(-1); }

// Read if Water Structures requested

 // WATER STRUCTURES keyword
  if(findf(fd, 2, "[WATER", "STRUCTURES]")==EOF){
    printf("Cannot find [WATER STRUCTURES] keyword\n");
    exit(-1);
  }

 // WATER STRUCTURES value
  if(fscanf(fd,"%d", &ws)!=1){
    printf("Cannot read YES (1) or NO (0) for WATER STRUCTURES keyword\n");
    exit(-1);
  }

  if(ws == 1){

     if(fscanf(fd,"%d", &watid)!=1){
       printf("Cannot read the order of WATER within solvent(s) for WATER STRUCTURES keyword\n");
       exit(-1);
     }

 }

  rewind(fd);

// Check if order of WATER is correct - start

if((ws == 1 && watid == 1) || (t1watdip == 1 && watid == 1) || (t2watdip == 1 && watid == 1)){

  if(findf(fd, 5, "[NUMBER", "OF", "ATOMS", "IN", "SOLVENT1]")==EOF){
     printf("SOLUTE1(2) WATER DIPOLE ORIENTATIONS/WATER STRUCTURES: Solvent 1 cannot be WATER as it does not exist\n");
     exit(-1);
  }

  fscanf(fd,"%d", &ntest);     
     
    slvntatmtest = (char **)malloc(sizeof(char *)*ntest);
    for(i=0; i < ntest; i++) slvntatmtest[i] = (char *)malloc(sizeof(char)*32);

    for(i=0; i < ntest; i++){
      fscanf(fd, "%31s %d", slvntatmtest[i], &scr[i]);
    }

  rewind(fd);      
  cks = 0;
  if(ntest != 3) { printf("SOLUTE1(2) WATER DIPOLE ORIENTATIONS/WATER STRUCTURES: Solvent 1 is NOT WATER\n"); exit(-1); }
  
  else{

     if(strcmp(slvntatmtest[0], "O") == 0 || strcmp(slvntatmtest[0], "H") == 0 || strcmp(slvntatmtest[0], "OW") == 0 || strcmp(slvntatmtest[0], "HW") == 0) cks = cks + 1; 
     if(strcmp(slvntatmtest[1], "O") == 0 || strcmp(slvntatmtest[1], "H") == 0 || strcmp(slvntatmtest[1], "OW") == 0 || strcmp(slvntatmtest[1], "HW") == 0) cks = cks + 1; 
     if(strcmp(slvntatmtest[2], "O") == 0 || strcmp(slvntatmtest[2], "H") == 0 || strcmp(slvntatmtest[2], "OW") == 0 || strcmp(slvntatmtest[2], "HW") == 0) cks = cks + 1; 

  } 
  if(cks != 3) { printf("SOLUTE1(2) WATER DIPOLE ORIENTATIONS/WATER STRUCTURES: Solvent 1 is NOT WATER\n"); exit(-1); }

}

if((ws == 1 && watid == 2) || (t1watdip == 1 && watid == 2) || (t2watdip == 1 && watid == 2)){

    if(findf(fd, 5, "[NUMBER", "OF", "ATOMS", "IN", "SOLVENT2]")==EOF){
       printf("SOLUTE1(2) WATER DIPOLE ORIENTATIONS/WATER STRUCTURES: Solvent 2 cannot be WATER as it does not exist\n");
       exit(-1);
    }
  
  fscanf(fd,"%d", &ntest);
     
    slvntatmtest = (char **)malloc(sizeof(char *)*ntest);
    for(i=0; i < ntest; i++) slvntatmtest[i] = (char *)malloc(sizeof(char)*32);
  
    for(i=0; i < ntest; i++){
      fscanf(fd, "%31s %d", slvntatmtest[i], &scr[i]);
    }
    
  rewind(fd);
  cks = 0;
  if(ntest != 3) { printf("SOLUTE1(2) WATER DIPOLE ORIENTATIONS/WATER STRUCTURES: Solvent 2 is NOT WATER\n"); exit(-1); }
  
  else{

     if(strcmp(slvntatmtest[0], "O") == 0 || strcmp(slvntatmtest[0], "H") == 0 || strcmp(slvntatmtest[0], "OW") == 0 || strcmp(slvntatmtest[0], "HW") == 0) cks = cks + 1;
     if(strcmp(slvntatmtest[1], "O") == 0 || strcmp(slvntatmtest[1], "H") == 0 || strcmp(slvntatmtest[1], "OW") == 0 || strcmp(slvntatmtest[1], "HW") == 0) cks = cks + 1;
     if(strcmp(slvntatmtest[2], "O") == 0 || strcmp(slvntatmtest[2], "H") == 0 || strcmp(slvntatmtest[2], "OW") == 0 || strcmp(slvntatmtest[2], "HW") == 0) cks = cks + 1;      

  }
  if(cks != 3) { printf("SOLUTE1(2) WATER DIPOLE ORIENTATIONS/WATER STRUCTURES: Solvent 2 is NOT WATER\n"); exit(-1); }

}

if((ws == 1 && watid == 3) || (t1watdip == 1 && watid == 3) || (t2watdip == 1 && watid == 3)){

    if(findf(fd, 5, "[NUMBER", "OF", "ATOMS", "IN", "SOLVENT3]")==EOF){
       printf("SOLUTE1(2) WATER DIPOLE ORIENTATIONS/WATER STRUCTURES: Solvent 3 cannot be WATER as it does not exist\n");
       exit(-1);
    }
  
  fscanf(fd,"%d", &ntest);
     
    slvntatmtest = (char **)malloc(sizeof(char *)*ntest);
    for(i=0; i < ntest; i++) slvntatmtest[i] = (char *)malloc(sizeof(char)*32);
  
    for(i=0; i < ntest; i++){
      fscanf(fd, "%31s %d", slvntatmtest[i], &scr[i]);
    }
    
  rewind(fd);
  cks = 0;
  if(ntest != 3) { printf("SOLUTE1(2) WATER DIPOLE ORIENTATIONS/WATER STRUCTURES: Solvent 3 is NOT WATER\n"); exit(-1); }
  
  else{

     if(strcmp(slvntatmtest[0], "O") == 0 || strcmp(slvntatmtest[0], "H") == 0 || strcmp(slvntatmtest[0], "OW") == 0 || strcmp(slvntatmtest[0], "HW") == 0) cks = cks + 1;
     if(strcmp(slvntatmtest[1], "O") == 0 || strcmp(slvntatmtest[1], "H") == 0 || strcmp(slvntatmtest[1], "OW") == 0 || strcmp(slvntatmtest[1], "HW") == 0) cks = cks + 1;
     if(strcmp(slvntatmtest[2], "O") == 0 || strcmp(slvntatmtest[2], "H") == 0 || strcmp(slvntatmtest[2], "OW") == 0 || strcmp(slvntatmtest[2], "HW") == 0) cks = cks + 1;     

  }
  if(cks != 3) { printf("SOLUTE1(2) WATER DIPOLE ORIENTATIONS/WATER STRUCTURES: Solvent 3 is NOT WATER\n"); exit(-1); } 

}

if(!((ws == 1 && watid == 1) || (ws == 1 && watid == 2) || (ws == 1 && watid == 3) || ws == 0))
{ printf("WATER STRUCTURES: Invalid input parameters\n"); exit(-1); }

if(!((t1watdip == 1 && watid == 1) || (t1watdip == 1 && watid == 2) || (t1watdip == 1 && watid == 3) || t1watdip == 0))
{ printf("SOLUTE1 WATER DIPOLE ORIENTATIONS: Invalid input parameters\n"); exit(-1); }

if(!((t2watdip == 1 && watid == 1) || (t2watdip == 1 && watid == 2) || (t2watdip == 1 && watid == 3) || t2watdip == 0))
{ printf("SOLUTE2 WATER DIPOLE ORIENTATIONS: Invalid input parameters\n"); exit(-1); }

// Check if order of WATER is correct - end



// IF WATER STRUCTURES is requested

if(ws == 1)
{

   if(findf(fd, 3, "[WATER", "HEXAMER", "RING]")==EOF){
    printf("Cannot find [WATER HEXAMER RING] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &hring)!=1){
    printf("Cannot read YES (1) or NO (0) for WATER HEXAMER RING keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 3, "[WATER", "HEXAMER", "BOOK]")==EOF){
    printf("Cannot find [WATER HEXAMER BOOK] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &hbook)!=1){
    printf("Cannot read YES (1) or NO (0) for WATER HEXAMER BOOK keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 3, "[WATER", "HEXAMER", "PRISM]")==EOF){
    printf("Cannot find [WATER HEXAMER PRISM keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &hprism)!=1){
    printf("Cannot read YES (1) or NO (0) for WATER HEXAMER PRISM keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 3, "[WATER", "HEXAMER", "CAGE]")==EOF){
    printf("Cannot find [WATER HEXAMER CAGE] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &hcage)!=1){
    printf("Cannot read YES (1) or NO (0) for WATER HEXAMER CAGE keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 3, "[WATER", "HEXAMER", "BAG]")==EOF){
    printf("Cannot find [WATER HEXAMER BAG] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &hbag)!=1){
    printf("Cannot read YES (1) or NO (0) for WATER HEXAMER BAG keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 3, "[WATER", "HEXAMER", "BOAT]")==EOF){
    printf("Cannot find [WATER HEXAMER BOAT] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &hboat)!=1){
    printf("Cannot read YES (1) or NO (0) for WATER HEXAMER BOAT keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 3, "[WATER", "HEXAMER", "CHAIR]")==EOF){
    printf("Cannot find [WATER HEXAMER CHAIR] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &hchair)!=1){
    printf("Cannot read YES (1) or NO (0) for WATER HEXAMER CHAIR keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 3, "[WATER", "HEXAMER", "PRISMBOOK]")==EOF){
    printf("Cannot find [WATER HEXAMER PRISMBOOK] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &hprismbook)!=1){
    printf("Cannot read YES (1) or NO (0) for WATER HEXAMER PRISMBOOK keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 3, "[WATER", "PENTAMER", "SEARCH]")==EOF){
    printf("Cannot find [WATER PENTAMER SEARCH] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &hpentamer)!=1){
    printf("Cannot read YES (1) or NO (0) for WATER PENTAMER SEARCH keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 3, "[WATER", "TETRAMER", "SEARCH]")==EOF){
    printf("Cannot find [WATER TETRAMER SEARCH] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &htetramer)!=1){
    printf("Cannot read YES (1) or NO (0) for WATER TETRAMER SEARCH keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 3, "[WATER", "TRIMER", "SEARCH]")==EOF){
    printf("Cannot find [WATER TRIMER SEARCH] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &htrimer)!=1){
    printf("Cannot read YES (1) or NO (0) for WATER TRIMER SEARCH keyword\n");
    exit(-1);
  } rewind(fd);

   if(findf(fd, 3, "[WATER", "ISOLATED", "STRUCTURES]")==EOF){
    printf("Cannot find [WATER ISOLATED STRUCTURES] keyword\n");
    exit(-1);
  }
  if(fscanf(fd,"%d", &hiso)!=1){
    printf("Cannot read YES (1) or NO (0) for WATER ISOLATED STRUCTURES keyword\n");
    exit(-1);
  } rewind(fd);

}

//POLYHEDRA Keyword
        if(findf(fd,1,"[POLYHEDRA]")==EOF){
            printf("Cannot find [POLYHEDRA] keyword\n");
            exit(-1);
        }
// POLYHEDRA value
        if(fscanf(fd,"%d",&findpolys)!=1){
            printf("Cannot read YES (1) or NO(0) for POLYHEDRA search\n");
            exit(-1);
        }
        rewind(fd);
        
// If POLYHEDRA requested
        if(findpolys==1)
        {
            if(nsolutetype!=1 || nsolventtype!=1)
            {
                printf("Polyhedra Search only valid with 1 Solvent and 1 Solute");
                exit(-1);
            }
            
            if(findf(fd, 3, "[MAX", "SHELL","SIZE]")==EOF){
                printf("Cannot find [MAX SHELL SIZE] keyword.\n");
                exit(-1);
            }
            if(fscanf(fd,"%d", &maxshellsize)!=1){
                printf("Cannot read parameter for MAX SHELL SIZE keyword.\n");
                exit(-1);
            }
            rewind(fd);
            
            edgebds=(double**)malloc((1+maxshellsize)*sizeof(double*));
            for(i=0; i<maxshellsize+1; i++)
            {
                edgebds[i]=(double*)malloc(2*sizeof(double));
            }
            
            //need better keyword ...
            if(findf(fd, 4, "[POLYHEDRA","EDGE","CERTAINTY", "BOUNDS]")==EOF){
                printf("Cannot find [POLYHEDRA EDGE CERTAINTY BOUNDS] keyword.\n");
                exit(-1);
            }
            if(fscanf(fd, "%d", &n_edg_bd)!=1)
            {
                printf("Cannot read parameter for POLYHEDRA EDGE CERTAINTY BOUNDS keyword.\n");
                exit(-1);
            }
            
                edgbd_idx=(int*)malloc(n_edg_bd*sizeof(int));
                
                for(i=0; i<n_edg_bd; i++)
                {
                    if(fscanf(fd, "%d", &edgbd_idx[i])!=1)
                    {
                        printf("Edge Certainty bounds: Reading bound shellsize %d failed", i);
                        exit(-1);
                    }
                    if(fscanf(fd, "%lf %lf",&(edgebds[edgbd_idx[i]][0]), &(edgebds[edgbd_idx[i]][1]))!=2)
                    {
                        printf("Edge Certainty bounds: Reading bounds %d failed",i);
                        exit(-1);
                    }
                }
                //fill in rest of edge bounds.
                
                for(i=0; i<edgbd_idx[0]; i++)
                {
                    edgebds[i][0]=edgebds[edgbd_idx[0]][0];
                    edgebds[i][1]=edgebds[edgbd_idx[0]][1];
                }
                for(i=0; i<n_edg_bd-1; i++)
                {
                    
                    for(j=edgbd_idx[i]; j<edgbd_idx[i+1]; j++) // Linear Extrapolation Between Given Edge Certainty Bounds
                    {
                        edgebds[j][0]=edgebds[edgbd_idx[i]][0] + ((float)(j-edgbd_idx[i]))/(edgbd_idx[i+1]-edgbd_idx[i])*(edgebds[edgbd_idx[i+1]][0]-edgebds[edgbd_idx[i]][0]);
                        edgebds[j][1]=edgebds[edgbd_idx[i]][1] + ((float)(j-edgbd_idx[i]))/(edgbd_idx[i+1]-edgbd_idx[i])*(edgebds[edgbd_idx[i+1]][1]-edgebds[edgbd_idx[i]][1]);
                        
                    }
                    
                }
                for(i=edgbd_idx[n_edg_bd-1]; i<maxshellsize+1; i++)
                {
                    edgebds[i][0]=edgebds[edgbd_idx[n_edg_bd-1]][0];
                    edgebds[i][1]=edgebds[edgbd_idx[n_edg_bd-1]][1];
                }
            
             // print for troubleshooting
            /*
                for(i=0; i<maxshellsize+1; i++)
                {
                    printf("%d \t %lf \t %lf \t \n", i, edgebds[i][0], edgebds[i][1]);
                }
               */ 
                
                    
            
            rewind(fd);
        
        }
        
        
        //PAIRWISE DISTS Keyword
        if(findf(fd, 3, "[PAIRED", "SHELL", "DISTANCES]")==EOF){
            printf("Cannot find [PAIRED SHELL DISTANCES] keyword\n");
            exit(-1);
        }
        //PAIRWISE DISTS values
        
        if (fscanf(fd, "%d", &findpdfshell)!=1){
            printf("Cannot read YES (1) or NO (0) for PAIRED SHELL DISTANCES\n");
            exit(-1);
        }
        rewind(fd);
        
        if(findpdfshell==1){
            if(nsolutetype!=1 || nsolventtype!=1)
            {
                printf("shell PDF calculation only valid with 1 Solvent");
                exit(-1);
            }
            
            if(findf(fd, 3, "[MAX", "SHELL","SIZE]")==EOF){
                printf("Cannot find [MAX SHELL SIZE] keyword.\n");
                exit(-1);
            }
            if(fscanf(fd,"%d", &maxshellsize)!=1){
                printf("Cannot read parameter for MAX SHELL SIZE keyword.\n");
                exit(-1);
            }
            rewind(fd);
        }
        
        //VARY SHELL DISTANCES Keyword (For ease of finding MCD)
        /*
        if (findf(fd, 3, "[VARY", "SHELL", "DISTANCES]")==EOF){
            printf("Cannot find [VARY SHELL DISTANCES] keyword.\n");
            exit(-1);
        }
        
        if (fscanf(fd, "%d", &varyshellsz)!=1){
            printf("Cannot read YES(1) or NO(0) for VARY SHELL DISTANCES.\n");
            exit(-1);
        }
        
        if(varyshellsz==1){
            
            if(fscanf(fd, "%lf %lf %d", &varyMin, &varyMax, &varyBreaks)!=3){
                printf("Cannot read parameters for VARY SHELL DISTANCES keyword.\n");
                exit(-1);
            }
        }*/

        

        
        

// Closing the input file

fclose(fd);

}

}

void ChemNetworkOrig::input_process(int argc, char *argv[])
{

// Read Coordinates from the XYZ files
    {
// Solvent

if((nsolventtype == 1 && nsolutetype == 0) || (nsolventtype == 1 && nsolutetype == 1) || nsolventtype == 2 || nsolventtype == 3)
{
   sprintf(finput, "%s", argv[2]);  // the first xyz file is the solvent1

   if((fsolvent1 = fopen(finput,"r")) == NULL){
     printf("Cannot open file %s\n", finput);
     exit(-1);
   }

   rewind(fsolvent1);

   if(fscanf(fsolvent1,"%d", &nAtomS1)!=1){
      printf("Cannot read number of atoms in solvent1 XYZ file\n");
      exit(-1);
   }

   // allocate memory

   atmS1 = (double*)malloc((3*nAtomS1)*sizeof(double));
  
   atmTypeS1 = (char **)malloc(sizeof(char *)*nAtomS1);
     for(i=0; i < nAtomS1; i++) atmTypeS1[i] = (char *)malloc(sizeof(char)*32);

   // Skip 2 lines from top of file
      rewind(fsolvent1); 
      fgets(line,256,fsolvent1);
      fgets(line,256,fsolvent1);

   for(i=0; i < nAtomS1; i++){
      if(fscanf(fsolvent1, "%31s %lf %lf %lf", atmTypeS1[i], &atmS1[3*(i+1)-3], &atmS1[3*(i+1)-2], &atmS1[3*(i+1)-1]) != 4){
        printf("Solvent1 XYZ file: Cannot read coordinates\n");
        exit(-1);
      }
    } 
   for(i = 0; i < nsolvent1; i++)
   {
     if(strcmp(slvntatm1[i], atmTypeS1[i]) != 0) {printf("Solvent1 XYZ file: Atom types and/or their order do not match with the input file\n"); exit(-1);}
   }
   fclose(fsolvent1);
}

// Solvent Solute

if((nsolventtype == 1 && nsolutetype == 1) || (nsolventtype == 1 && nsolutetype == 2))
{
   sprintf(finput, "%s", argv[3]);   // 1 solvent with solute, the second xyz file is solute

   if((fsolute1 = fopen(finput,"r")) == NULL){
     printf("Cannot open file %s\n", finput);
     exit(-1);
   }

   rewind(fsolute1);

   if(fscanf(fsolute1,"%d", &nAtomT1)!=1){
      printf("Cannot read number of atoms in solute1 XYZ file\n");
      exit(-1);
   }

   // allocate memory

   atmT1 = (double*)malloc((3*nAtomT1)*sizeof(double));

   atmTypeT1 = (char **)malloc(sizeof(char *)*nAtomT1);
     for(i=0; i < nAtomT1; i++) atmTypeT1[i] = (char *)malloc(sizeof(char)*32);

   // Skip 2 lines from top of file
      rewind(fsolute1);
      fgets(line,256,fsolute1);
      fgets(line,256,fsolute1);

   for(i=0; i < nAtomT1; i++){
      if(fscanf(fsolute1, "%31s %lf %lf %lf", atmTypeT1[i], &atmT1[3*(i+1)-3], &atmT1[3*(i+1)-2], &atmT1[3*(i+1)-1]) != 4){
        printf("Solute1 XYZ file: Cannot read coordinates\n");
        exit(-1);
      }
    }
   for(i = 0; i < nsolute1; i++)
   {
     if(strcmp(sltatm1[i], atmTypeT1[i]) != 0) {printf("Solute1 XYZ file: Atom types and/or their order do not match with the input file\n"); exit(-1);}
   }
   fclose(fsolute1);
}

// Solvent Solute Solute

if(nsolventtype == 1 && nsolutetype == 2)
{
   sprintf(finput, "%s", argv[4]);

   if((fsolute2 = fopen(finput,"r")) == NULL){
     printf("Cannot open file %s\n", finput);
     exit(-1);
   }

   rewind(fsolute2);

   if(fscanf(fsolute2,"%d", &nAtomT2)!=1){
      printf("Cannot read number of atoms in solute2 XYZ file\n");
      exit(-1);
   }

   // allocate memory

   atmT2 = (double*)malloc((3*nAtomT2)*sizeof(double));

   atmTypeT2 = (char **)malloc(sizeof(char *)*nAtomT2);
     for(i=0; i < nAtomT2; i++) atmTypeT2[i] = (char *)malloc(sizeof(char)*32);

   // Skip 2 lines from top of file
      rewind(fsolute2);
      fgets(line,256,fsolute2);
      fgets(line,256,fsolute2);

   for(i=0; i < nAtomT2; i++){
      if(fscanf(fsolute2, "%31s %lf %lf %lf", atmTypeT2[i], &atmT2[3*(i+1)-3], &atmT2[3*(i+1)-2], &atmT2[3*(i+1)-1]) != 4){
        printf("Solute2 XYZ file: Cannot read coordinates\n");
        exit(-1);
      }
    }
   for(i = 0; i < nsolute2; i++)
   {
     if(strcmp(sltatm2[i], atmTypeT2[i]) != 0) {printf("Solute2 XYZ file: Atom types and/or their order do not match with the input file\n"); exit(-1);}
   }
   fclose(fsolute2);
}

// Solvent Solvent

if((nsolventtype == 2 && nsolutetype == 0) || (nsolventtype == 2 && nsolutetype == 1) || (nsolventtype == 2 && nsolutetype == 2) || nsolventtype == 3)
{
   sprintf(finput, "%s", argv[3]);

   if((fsolvent2 = fopen(finput,"r")) == NULL){
     printf("Cannot open file %s\n", finput);
     exit(-1);
   }

   rewind(fsolvent2);

   if(fscanf(fsolvent2,"%d", &nAtomS2)!=1){
      printf("Cannot read number of atoms in solvent2 XYZ file\n");
      exit(-1);
   }

   // allocate memory

   atmS2 = (double*)malloc((3*nAtomS2)*sizeof(double));

   atmTypeS2 = (char **)malloc(sizeof(char *)*nAtomS2);
     for(i=0; i < nAtomS2; i++) atmTypeS2[i] = (char *)malloc(sizeof(char)*32);

   // Skip 2 lines from top of file
      rewind(fsolvent2);
      fgets(line,256,fsolvent2);
      fgets(line,256,fsolvent2);

   for(i=0; i < nAtomS2; i++){
      if(fscanf(fsolvent2, "%31s %lf %lf %lf", atmTypeS2[i], &atmS2[3*(i+1)-3], &atmS2[3*(i+1)-2], &atmS2[3*(i+1)-1]) != 4){
        printf("Solvent2 XYZ file: Cannot read coordinates\n");
        exit(-1);
      }
    }
   for(i = 0; i < nsolvent2; i++)
   {
     if(strcmp(slvntatm2[i], atmTypeS2[i]) != 0) {printf("Solvent2 XYZ file: Atom types and/or their order do not match with the input file\n"); exit(-1);}
   }
   fclose(fsolvent2);
}

// Solvent Solvent Solute

if((nsolventtype == 2 && nsolutetype == 1) || (nsolventtype == 2 && nsolutetype == 2))
{

   sprintf(finput, "%s", argv[4]);

   if((fsolute1 = fopen(finput,"r")) == NULL){
     printf("Cannot open file %s\n", finput);
     exit(-1);
   }

   rewind(fsolute1);

   if(fscanf(fsolute1,"%d", &nAtomT1)!=1){
      printf("Cannot read number of atoms in solute1 XYZ file\n");
      exit(-1);
   }

   // allocate memory

   atmT1 = (double*)malloc((3*nAtomT1)*sizeof(double));

   atmTypeT1 = (char **)malloc(sizeof(char *)*nAtomT1);
     for(i=0; i < nAtomT1; i++) atmTypeT1[i] = (char *)malloc(sizeof(char)*32);

   // Skip 2 lines from top of file
      rewind(fsolute1);
      fgets(line,256,fsolute1);
      fgets(line,256,fsolute1);

   for(i=0; i < nAtomT1; i++){
      if(fscanf(fsolute1, "%31s %lf %lf %lf", atmTypeT1[i], &atmT1[3*(i+1)-3], &atmT1[3*(i+1)-2], &atmT1[3*(i+1)-1]) != 4){
        printf("Solute1 XYZ file: Cannot read coordinates\n");
        exit(-1);
      }
    }
   for(i = 0; i < nsolute1; i++)
   {
     if(strcmp(sltatm1[i], atmTypeT1[i]) != 0) {printf("Solute1 XYZ file: Atom types and/or their order do not match with the input file\n"); exit(-1);}
   }
   fclose(fsolute1);

}

// Solvent Solvent Solute Solute

if(nsolventtype == 2 && nsolutetype == 2)
{

   sprintf(finput, "%s", argv[5]);

   if((fsolute2 = fopen(finput,"r")) == NULL){
     printf("Cannot open file %s\n", finput);
     exit(-1);
   }

   rewind(fsolute2);

   if(fscanf(fsolute2,"%d", &nAtomT2)!=1){
      printf("Cannot read number of atoms in solute2 XYZ file\n");
      exit(-1);
   }

   // allocate memory

   atmT2 = (double*)malloc((3*nAtomT2)*sizeof(double));

   atmTypeT2 = (char **)malloc(sizeof(char *)*nAtomT2);
     for(i=0; i < nAtomT2; i++) atmTypeT2[i] = (char *)malloc(sizeof(char)*32);

   // Skip 2 lines from top of file
      rewind(fsolute2);
      fgets(line,256,fsolute2);
      fgets(line,256,fsolute2);

   for(i=0; i < nAtomT2; i++){
      if(fscanf(fsolute2, "%31s %lf %lf %lf", atmTypeT2[i], &atmT2[3*(i+1)-3], &atmT2[3*(i+1)-2], &atmT2[3*(i+1)-1]) != 4){
        printf("Solute2 XYZ file: Cannot read coordinates\n");
        exit(-1);
      }
    }
   for(i = 0; i < nsolute2; i++)
   {
     if(strcmp(sltatm2[i], atmTypeT2[i]) != 0) {printf("Solute2 XYZ file: Atom types and/or their order do not match with the input file\n"); exit(-1);}
   }
   fclose(fsolute2);

}

// Solvent Solvent Solvent

if(nsolventtype == 3)
{

   sprintf(finput, "%s", argv[4]);

   if((fsolvent3 = fopen(finput,"r")) == NULL){
     printf("Cannot open file %s\n", finput);
     exit(-1);
   }

   rewind(fsolvent3);

   if(fscanf(fsolvent3,"%d", &nAtomS3)!=1){
      printf("Cannot read number of atoms in solvent3 XYZ file\n");
      exit(-1);
   }

   // allocate memory

   atmS3 = (double*)malloc((3*nAtomS3)*sizeof(double));

   atmTypeS3 = (char **)malloc(sizeof(char *)*nAtomS3);
     for(i=0; i < nAtomS3; i++) atmTypeS3[i] = (char *)malloc(sizeof(char)*32);

   // Skip 2 lines from top of file
      rewind(fsolvent3);
      fgets(line,256,fsolvent3);
      fgets(line,256,fsolvent3);

   for(i=0; i < nAtomS3; i++){
      if(fscanf(fsolvent3, "%31s %lf %lf %lf", atmTypeS3[i], &atmS3[3*(i+1)-3], &atmS3[3*(i+1)-2], &atmS3[3*(i+1)-1]) != 4){
        printf("Solvent3 XYZ file: Cannot read coordinates\n");
        exit(-1);
      }
    }
   for(i = 0; i < nsolvent3; i++)
   {
     if(strcmp(slvntatm3[i], atmTypeS3[i]) != 0) {printf("Solvent3 XYZ file: Atom types and/or their order do not match with the input file\n"); exit(-1);}
   }
   fclose(fsolvent3);

}

// Solvent Solvent Solvent Solute

if((nsolventtype == 3 && nsolutetype == 1) || (nsolventtype == 3 && nsolutetype == 2))
{

   sprintf(finput, "%s", argv[5]);

   if((fsolute1 = fopen(finput,"r")) == NULL){
     printf("Cannot open file %s\n", finput);
     exit(-1);
   }

   rewind(fsolute1);

   if(fscanf(fsolute1,"%d", &nAtomT1)!=1){
      printf("Cannot read number of atoms in solute1 XYZ file\n");
      exit(-1);
   }

   // allocate memory

   atmT1 = (double*)malloc((3*nAtomT1)*sizeof(double));

   atmTypeT1 = (char **)malloc(sizeof(char *)*nAtomT1);
     for(i=0; i < nAtomT1; i++) atmTypeT1[i] = (char *)malloc(sizeof(char)*32);

   // Skip 2 lines from top of file
      rewind(fsolute1);
      fgets(line,256,fsolute1);
      fgets(line,256,fsolute1);

   for(i=0; i < nAtomT1; i++){
      if(fscanf(fsolute1, "%31s %lf %lf %lf", atmTypeT1[i], &atmT1[3*(i+1)-3], &atmT1[3*(i+1)-2], &atmT1[3*(i+1)-1]) != 4){
        printf("Solute1 XYZ file: Cannot read coordinates\n");
        exit(-1);
      }
    }
   for(i = 0; i < nsolute1; i++)
   {
     if(strcmp(sltatm1[i], atmTypeT1[i]) != 0) {printf("Solute1 XYZ file: Atom types and/or their order do not match with the input file\n"); exit(-1);}
   }
   fclose(fsolute1);

}

// Solvent Solvent Solvent Solute Solute

if(nsolventtype == 3 && nsolutetype == 2)
{

   sprintf(finput, "%s", argv[6]);

   if((fsolute2 = fopen(finput,"r")) == NULL){
     printf("Cannot open file %s\n", finput);
     exit(-1);
   }

   rewind(fsolute2);

   if(fscanf(fsolute2,"%d", &nAtomT2)!=1){
      printf("Cannot read number of atoms in solute2 XYZ file\n");
      exit(-1);
   }

   // allocate memory

   atmT2 = (double*)malloc((3*nAtomT2)*sizeof(double));

   atmTypeT2 = (char **)malloc(sizeof(char *)*nAtomT2);
     for(i=0; i < nAtomT2; i++) atmTypeT2[i] = (char *)malloc(sizeof(char)*32);

   // Skip 2 lines from top of file
      rewind(fsolute2);
      fgets(line,256,fsolute2);
      fgets(line,256,fsolute2);

   for(i=0; i < nAtomT2; i++){
      if(fscanf(fsolute2, "%31s %lf %lf %lf", atmTypeT2[i], &atmT2[3*(i+1)-3], &atmT2[3*(i+1)-2], &atmT2[3*(i+1)-1]) != 4){
        printf("Solute2 XYZ file: Cannot read coordinates\n");
        exit(-1);
      }
    }
   for(i = 0; i < nsolute2; i++)
   {
     if(strcmp(sltatm2[i], atmTypeT2[i]) != 0) {printf("Solute2 XYZ file: Atom types and/or their order do not match with the input file\n"); exit(-1);}
   }
   fclose(fsolute2);

}
    
//ensure polyhedral search valid
        
    }
// End reading Input Files


// Output Test
//
// sprintf(foutput, "%s.output.XYZ.xyz", argv[2]);
// outputf=fopen(foutput,"w");
//
//       for(i=0;i<nAtomS1;i++)
//       {   
//           fprintf(outputf,"%s  %lf %lf %lf\n",atmTypeS1[i],atmS1[3*(i+1)-3],atmS1[3*(i+1)-2],atmS1[3*(i+1)-1]);
//       } 
//

// General Output to LOG File

if((nsolventtype + nsolutetype) == 1)
   sprintf(foutput, "%s.%s.log", argv[1], argv[2]);

if((nsolventtype + nsolutetype) == 2)
   sprintf(foutput, "%s.%s.%s.log", argv[1], argv[2], argv[3]);

if((nsolventtype + nsolutetype) == 3)
   sprintf(foutput, "%s.%s.%s.%s.log", argv[1], argv[2], argv[3], argv[4]);

if((nsolventtype + nsolutetype) == 4)
   sprintf(foutput, "%s.%s.%s.%s.%s.log", argv[1], argv[2], argv[3], argv[4], argv[5]);

if((nsolventtype + nsolutetype) == 5)
   sprintf(foutput, "%s.%s.%s.%s.%s.%s.log", argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);

/*

 outputf=fopen(foutput,"w");

 fprintf(outputf,
  "*------------------------------------------------------*\n" 
  "| Chemical Networks Program (ChemNetworks)             |\n"
  "| Version 1.0  June 29, 2013                           |\n"
  "| Dr. Abdullah Ozkanlar                                |\n"
  "| A. Clark Research Lab, Department of Chemistry       |\n"
  "| Washington State University, Pullman/WA 99164        |\n"
  "*------------------------------------------------------*\n\n");

time(&now);
fprintf(outputf,"Start time: %s", ctime(&now));

if((nsolventtype + nsolutetype) == 1) fprintf(outputf,"\nCoordinates File 1: %s\n", argv[2]);

if((nsolventtype + nsolutetype) == 2) 
{
  fprintf(outputf,"\nCoordinates File 1: %s\n", argv[2]); 
  fprintf(outputf,"\nCoordinates File 2: %s\n", argv[3]); 
}

if((nsolventtype + nsolutetype) == 3)
{
  fprintf(outputf,"\nCoordinates File 1: %s\n", argv[2]);
  fprintf(outputf,"\nCoordinates File 2: %s\n", argv[3]);
  fprintf(outputf,"\nCoordinates File 3: %s\n", argv[4]);
} 
if((nsolventtype + nsolutetype) == 4)
{
  fprintf(outputf,"\nCoordinates File 1: %s\n", argv[2]);
  fprintf(outputf,"\nCoordinates File 2: %s\n", argv[3]);
  fprintf(outputf,"\nCoordinates File 3: %s\n", argv[4]);
  fprintf(outputf,"\nCoordinates File 4: %s\n", argv[5]);
} 
if((nsolventtype + nsolutetype) == 5) 
{
  fprintf(outputf,"\nCoordinates File 1: %s\n", argv[2]);
  fprintf(outputf,"\nCoordinates File 2: %s\n", argv[3]);
  fprintf(outputf,"\nCoordinates File 3: %s\n", argv[4]);
  fprintf(outputf,"\nCoordinates File 4: %s\n", argv[5]);
  fprintf(outputf,"\nCoordinates File 5: %s\n", argv[6]);
}

*/

    
// Get the positions (orders) of water atoms
opos = 0; h1pos = 0; h2pos = 0;
trs = 0;

for(i = 0; i < 3; i++)    // 2015.12.14, O and OW refer to water Oxygen
{
  if(watid == 1)   
  { 
    if(strcmp(slvntatm1[i], "O") == 0 || strcmp(slvntatm1[i], "OW") == 0) opos = slvntatmnum1[i];
    if(strcmp(slvntatm1[i], "H") == 0 || strcmp(slvntatm1[i], "HW") == 0){ if(trs == 0){ h1pos = slvntatmnum1[i]; trs = trs + 1; } else h2pos = slvntatmnum1[i]; }
  }
  if(watid == 2)
  { 
    if(strcmp(slvntatm2[i], "O") == 0 || strcmp(slvntatm2[i], "OW") == 0) opos = slvntatmnum2[i];
    if(strcmp(slvntatm2[i], "H") == 0 || strcmp(slvntatm2[i], "HW") == 0){ if(trs == 0){ h1pos = slvntatmnum2[i]; trs = trs + 1; } else h2pos = slvntatmnum2[i]; } 
  }
  if(watid == 3)
  { 
    if(strcmp(slvntatm3[i], "O") == 0 || strcmp(slvntatm3[i], "OW") == 0) opos = slvntatmnum3[i]; 
    if(strcmp(slvntatm3[i], "H") == 0 || strcmp(slvntatm3[i], "HW") == 0){ if(trs == 0){ h1pos = slvntatmnum3[i]; trs = trs + 1; } else h2pos = slvntatmnum3[i]; }
  }
}




// construct the requested Weighted Graph, Feb-2018
if(cnwg.index_weighted_graph == 1)   // when the weighted graph is requested
{
  sprintf(cnwg.foutput_weighted_graph, "Weighted.Graph.%s",argv[1]);
  cnwg.output_weighted_graph = fopen(cnwg.foutput_weighted_graph,"w");

  /* select the molecules/atoms within the cluster calculation */
  if(cnwg.index_wg_st1_cluster >= 1)
  {
     if(nAtomT1 % nsolute1 != 0){
       printf("Error in the coordinate file of solute1\n");
       exit(-1);
     }

     cnwg.num_mol_cluster_st1 = nAtomT1 / nsolute1;

     cnwg.WG_Mol_id = (struct Mol_identity *)malloc( cnwg.num_mol_cluster_st1 * sizeof(struct Mol_identity) );
     for(i=0; i < cnwg.num_mol_cluster_st1; i++){  // keep record of the solute1 id
       cnwg.WG_Mol_id[i].solute_type = 1;
       cnwg.WG_Mol_id[i].solvent_type = 0;
       cnwg.WG_Mol_id[i].id = i+1;
     }

     cnwg.atom_cluster_st1 = (double *)malloc( (cnwg.num_mol_cluster_st1 * nsolute1 * 3) * sizeof(double));
     for(i=0; i < cnwg.num_mol_cluster_st1 * nsolute1 * 3; i++){
       cnwg.atom_cluster_st1[i] = atmT1[i];
     }

     if(nsolventtype >= 1 && cnwg.index_wg_st1_sv1 == 1 && cnwg.findsolvent(1, cnwg.index_wg_st1_cluster, cnwg.cluster_st1_sv_type) == 1)
     {
       if(nAtomS1 % nsolvent1 != 0){
         printf("Error in coordinate file of solvent1\n");
         exit(-1);
       }

       cnwg.id_sv1 = 0;
       for(i=0; i < cnwg.index_wg_st1_cluster; i++){  
         if(cnwg.cluster_st1_sv_type[i] == 1){
           cnwg.id_sv1 = i;      // solvent1 is in line 'id_sv1' in [WEIGHTED GRAPH BY SOLUTE1 CLUSTER] parameters
           break;
         }
       }

       cnwg.num_mol_cluster_sv1 = 0;

       for(i=0; i < nAtomS1/nsolvent1; i++){
         cnwg.index_select_mol = 0;

         for(j=0; j < cnwg.num_mol_cluster_st1; j++){
           cnwg.temp_wg_site_dist = cnwg.wg_site_distance(atmT1, j, nsolute1, cnwg.cluster_st1_atom[cnwg.id_sv1], atmS1, i, nsolvent1, cnwg.cluster_st1_sv_atom[cnwg.id_sv1], xside, yside, zside);
           if(cnwg.temp_wg_site_dist >= cnwg.cluster_st1_sv_dist_min[cnwg.id_sv1] && cnwg.temp_wg_site_dist <= cnwg.cluster_st1_sv_dist_max[cnwg.id_sv1]){ // select solvent1 into graph cluster, based on distance to solute1
             cnwg.index_select_mol = 1;
             break;
           };
         }

         if(cnwg.index_select_mol == 1){    // copy the coordinate of solvent1 molecule i to cluster
           if(cnwg.num_mol_cluster_sv1 == 0) {
             cnwg.atom_cluster_sv1 = (double *)malloc( (1 * nsolvent1 * 3) * sizeof(double));  // size of one solvent molecule
           }
           else {
             cnwg.atom_cluster_temp = (double *)realloc( cnwg.atom_cluster_sv1 , ((cnwg.num_mol_cluster_sv1 + 1) * nsolvent1 * 3) * sizeof(double)); // adjust the size by increasing one solvent molecule
             if(cnwg.atom_cluster_temp == NULL){
               printf("Cannot allocate memory for coordinates of atom cluster in solvent1\n");
               exit(-1);
             }
             else
               cnwg.atom_cluster_sv1 = cnwg.atom_cluster_temp; 
          }

           for(j = i*nsolvent1*3; j < (i+1)*nsolvent1*3; j++){
             cnwg.atom_cluster_sv1[ j - (i-cnwg.num_mol_cluster_sv1)*nsolvent1*3 ] = atmS1[j];
           }

           cnwg.num_mol_cluster_sv1++;

           cnwg.temp_Mol_id = (struct Mol_identity *)realloc( cnwg.WG_Mol_id, (cnwg.num_mol_cluster_st1 + cnwg.num_mol_cluster_sv1) * sizeof(struct Mol_identity) );  // adjust the size of WG_Mol_id
           if(cnwg.temp_Mol_id == NULL){
             printf("Cannot allocate memory for molecule identity at solvent1\n");
             exit(-1);
           }
           else {    // keep record of the solute1 id
             cnwg.WG_Mol_id = cnwg.temp_Mol_id;    
             cnwg.WG_Mol_id[cnwg.num_mol_cluster_st1 + cnwg.num_mol_cluster_sv1 - 1].solute_type = 0;
             cnwg.WG_Mol_id[cnwg.num_mol_cluster_st1 + cnwg.num_mol_cluster_sv1 - 1].solvent_type = 1;
             cnwg.WG_Mol_id[cnwg.num_mol_cluster_st1 + cnwg.num_mol_cluster_sv1 - 1].id = i+1;     // id of that molecule is 1 to num_of_that_molecule
           }

         } // this is the end of ' if(cnwg.index_select_mol == 1) '

       }  // this is the end of ' for(i=0; i < nAtomS1/nsolvent1; i++) '

     } // this is the end of ' if(nsolventtype >= 1 && cnwg.index_wg_st1_sv1 == 1 && cnwg.findsolvent(1, cnwg.index_wg_st1_cluster, cnwg.cluster_st1_sv_type) == 1) 


     if(nsolventtype >= 2 && cnwg.index_wg_st1_sv2 == 1 && cnwg.findsolvent(2, cnwg.index_wg_st1_cluster, cnwg.cluster_st1_sv_type) == 1)  // read the cluster for solvent2-solute1
     {
       if(nAtomS2 % nsolvent2 != 0){
         printf("Error in coordinate file of solvent2\n");
         exit(-1);
       }

       cnwg.id_sv2 = 0;
       for(i=0; i < cnwg.index_wg_st1_cluster; i++){  
         if(cnwg.cluster_st1_sv_type[i] == 2){
           cnwg.id_sv2 = i;      // solvent2 is in line 'id_sv2' in [WEIGHTED GRAPH BY SOLUTE1 CLUSTER] parameters
           break;
         }
       }

       cnwg.num_mol_cluster_sv2 = 0;

       for(i=0; i < nAtomS2/nsolvent2; i++){
         cnwg.index_select_mol = 0;

         for(j=0; j < cnwg.num_mol_cluster_st1; j++){
           cnwg.temp_wg_site_dist = cnwg.wg_site_distance(atmT1, j, nsolute1, cnwg.cluster_st1_atom[cnwg.id_sv2], atmS2, i, nsolvent2, cnwg.cluster_st1_sv_atom[cnwg.id_sv2], xside, yside, zside);
           if(cnwg.temp_wg_site_dist >= cnwg.cluster_st1_sv_dist_min[cnwg.id_sv2] && cnwg.temp_wg_site_dist <= cnwg.cluster_st1_sv_dist_max[cnwg.id_sv2]){ // select solvent2 into graph cluster, based on distance to solute1
             cnwg.index_select_mol = 1;
             break;
           };
         }

         if(cnwg.index_select_mol == 1){    // copy the coordinate of solvent2 molecule i to cluster
           if(cnwg.num_mol_cluster_sv2 == 0) {
             cnwg.atom_cluster_sv2 = (double *)malloc( (1 * nsolvent2 * 3) * sizeof(double));  // size of one solvent molecule
           }
           else {
             cnwg.atom_cluster_temp = (double *)realloc( cnwg.atom_cluster_sv2 , ((cnwg.num_mol_cluster_sv2 + 1) * nsolvent2 * 3) * sizeof(double)); // adjust the size by increasing one solvent molecule
             if(cnwg.atom_cluster_temp == NULL){
               printf("Cannot allocate memory for coordinates of atom cluster in solvent2\n");
               exit(-1);
             }
             else
               cnwg.atom_cluster_sv2 = cnwg.atom_cluster_temp; 
          }

           for(j = i*nsolvent2*3; j < (i+1)*nsolvent2*3; j++){
             cnwg.atom_cluster_sv2[ j - (i-cnwg.num_mol_cluster_sv2)*nsolvent2*3 ] = atmS2[j];
           }

           cnwg.num_mol_cluster_sv2++;

           cnwg.temp_Mol_id = (struct Mol_identity *)realloc( cnwg.WG_Mol_id, (cnwg.num_mol_cluster_st1 + cnwg.num_mol_cluster_sv1 + cnwg.num_mol_cluster_sv2) * sizeof(struct Mol_identity) );  // adjust the size of WG_Mol_id
           if(cnwg.temp_Mol_id == NULL){
             printf("Cannot allocate memory for molecule identity at solvent2\n");
             exit(-1);
           }
           else {    // keep record of the solvent2 id, which is after solute1 & solvent1
             cnwg.WG_Mol_id = cnwg.temp_Mol_id;    
             cnwg.WG_Mol_id[cnwg.num_mol_cluster_st1 + cnwg.num_mol_cluster_sv1 + cnwg.num_mol_cluster_sv2 - 1].solute_type = 0;
             cnwg.WG_Mol_id[cnwg.num_mol_cluster_st1 + cnwg.num_mol_cluster_sv1 + cnwg.num_mol_cluster_sv2 - 1].solvent_type = 2;  // it is solvent 2
             cnwg.WG_Mol_id[cnwg.num_mol_cluster_st1 + cnwg.num_mol_cluster_sv1 + cnwg.num_mol_cluster_sv2 - 1].id = i+1;
           }

         } // this is the end of ' if(cnwg.index_select_mol == 1) '

       }  // this is the end of ' for(i=0; i < nAtomS2/nsolvent2; i++) '

     } // this is the end of ' if(nsolventtype >= 2 && cnwg.index_wg_st1_sv2 == 1 && cnwg.findsolvent(2, cnwg.index_wg_st1_cluster, cnwg.cluster_st1_sv_type) == 1) 


    /* allocate the Adjacency matrix for weighted graph */
    cnwg.total_num_nodes = cnwg.num_mol_cluster_st1 + cnwg.num_mol_cluster_sv1 + cnwg.num_mol_cluster_sv2;
    cnwg.WG_Adj = (double **)malloc( cnwg.total_num_nodes * sizeof(double *) ); 
    for(i=0; i< cnwg.total_num_nodes; i++) {
      cnwg.WG_Adj[i] = (double *)calloc( cnwg.total_num_nodes, sizeof(double) );
    }

    /* create the weighted graph, based on the selected cluster from solute1, solvent1, solvent2 coordinates */
    cnwg.index_create_Adj = cnwg.create_WG_Adj_from_cluster(cnwg.output_weighted_graph, cnwg.WG_Adj, cnwg.WG_Mol_id, 
                              cnwg.atom_cluster_st1, cnwg.num_mol_cluster_st1, nsolute1, cnwg.atom_cluster_sv1, cnwg.num_mol_cluster_sv1, nsolvent1, cnwg.atom_cluster_sv2, cnwg.num_mol_cluster_sv2, nsolvent2, 
                              cnwg.index_wg_st1_sv1, cnwg.num_wg_st1_sv1_dist, cnwg.atom1_wg_st1_sv1, cnwg.atom2_wg_st1_sv1, cnwg.funct_type_wg_st1_sv1, cnwg.funct_par1_wg_st1_sv1, cnwg.funct_par2_wg_st1_sv1,
                              cnwg.index_wg_st1_sv2, cnwg.num_wg_st1_sv2_dist, cnwg.atom1_wg_st1_sv2, cnwg.atom2_wg_st1_sv2, cnwg.funct_type_wg_st1_sv2, cnwg.funct_par1_wg_st1_sv2, cnwg.funct_par2_wg_st1_sv2, 
                              cnwg.index_wg_sv1_sv1, cnwg.num_wg_sv1_sv1_dist, cnwg.atom1_wg_sv1_sv1, cnwg.atom2_wg_sv1_sv1, cnwg.funct_type_wg_sv1_sv1, cnwg.funct_par1_wg_sv1_sv1, cnwg.funct_par2_wg_sv1_sv1,
                              cnwg.index_wg_sv2_sv2, cnwg.num_wg_sv2_sv2_dist, cnwg.atom1_wg_sv2_sv2, cnwg.atom2_wg_sv2_sv2, cnwg.funct_type_wg_sv2_sv2, cnwg.funct_par1_wg_sv2_sv2, cnwg.funct_par2_wg_sv2_sv2,
                              cnwg.index_wg_sv1_sv2, cnwg.num_wg_sv1_sv2_dist, cnwg.atom1_wg_sv1_sv2, cnwg.atom2_wg_sv1_sv2, cnwg.funct_type_wg_sv1_sv2, cnwg.funct_par1_wg_sv1_sv2, cnwg.funct_par2_wg_sv1_sv2,
                              xside, yside, zside);

    
    
  } // this is the end of ' if(index_wg_st1_cluster == 1) '

 
  fclose(cnwg.output_weighted_graph);

}





// Construct the requested Graph Types

if(gs1s1 == 1)   // 2015.12.14, Output file name format here, Tiecheng
{
 sprintf(foutputGraphS1S1, "%s.%s.%s.Graph", argv[1], argv[2], argv[2]);
 outputfGraphS1S1 = fopen(foutputGraphS1S1,"w");

 sprintf(foutputGeodS1S1, "%s.%s.%s.GraphGeod", argv[1], argv[2], argv[2]);
 outputfGeodS1S1 = fopen(foutputGeodS1S1,"w");

 if(pnumnodes==1){
 sprintf(foutputNumnodesS1S1, "%s.%s.%s.GraphNumnodes", argv[1], argv[2], argv[2]);
 outputfNumnodesS1S1 = fopen(foutputNumnodesS1S1,"w"); }

 // Construct the solvent1-solvent1 graph

  // add enegetic calculation is requested for water
  if(E_s1s1_num > 0 && E_s1s1_charge_num > 0 && E_s1s1_LJ_num > 0) 
  {
    if( strcmp(slvntatm1[0],"O")==0 && strcmp(slvntatm1[1],"H")==0 && strcmp(slvntatm1[1],"H")==0 )  // energetics is implemented for water as solvent1 only
    {
      graph_ss_E(atmS1, 1, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBDmin, s1as1bBDmax, s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6, s1s1v7,
                 pbc, xside, yside, zside, outputfGraphS1S1, outputfGeodS1S1, E_s1s1_num, E_s1s1v1, E_s1s1v2, E_s1s1_min, E_s1s1_max, E_s1s1_charge_num,
                 E_s1s1_charge_value, E_s1s1_LJ_num, E_s1s1_LJ_index_a, E_s1s1_LJ_index_b, E_s1s1_LJ_value_sigma, E_s1s1_LJ_value_epsilon);    
    }
    else
    {
      printf("Warning: energetic calculation is only implimented for water\n");

// graph_ss(atmS1, 1, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBD, s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6,
//          pbc, xside, yside, zside, outputfGraphS1S1, outputfGeodS1S1);
      graph_ss(atmS1, 1, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBDmin, s1as1bBDmax, s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6, s1s1v7,
               pbc, xside, yside, zside, outputfGraphS1S1, outputfGeodS1S1);    // 2015.12.15, I changed here, Tiecheng
    }
  }
  else  // the pair energy definition is not requested, use the previous function
  {
// graph_ss(atmS1, 1, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBD, s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6,
//          pbc, xside, yside, zside, outputfGraphS1S1, outputfGeodS1S1);
  graph_ss(atmS1, 1, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBDmin, s1as1bBDmax, s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6, s1s1v7,
           pbc, xside, yside, zside, outputfGraphS1S1, outputfGeodS1S1);    // 2015.12.15, I changed here, Tiecheng
  }


if(pnumnodes==1) fprintf(outputfNumnodesS1S1,"%d\n",(nAtomS1 / nsolvent1));
// fprintf(outputf,"\n[GRAPH SOLVENT1 SOLVENT1]: Done...\n");
}

if(gs2s2 == 1)
{
 sprintf(foutputGraphS2S2, "%s.%s.%s.Graph", argv[1], argv[3], argv[3]);
 outputfGraphS2S2 = fopen(foutputGraphS2S2,"w");

 sprintf(foutputGeodS2S2, "%s.%s.%s.GraphGeod", argv[1], argv[3], argv[3]);
 outputfGeodS2S2 = fopen(foutputGeodS2S2,"w");

 if(pnumnodes==1){
 sprintf(foutputNumnodesS2S2, "%s.%s.%s.GraphNumnodes", argv[1], argv[3], argv[3]);
 outputfNumnodesS2S2 = fopen(foutputNumnodesS2S2,"w"); }

 // Construct the solvent2-solvent2 graph

// graph_ss(atmS2, 1, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBD, s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6,
//          pbc, xside, yside, zside, outputfGraphS2S2, outputfGeodS2S2);
 graph_ss(atmS2, 1, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBDmin, s2as2bBDmax, s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6, s2s2v7,
          pbc, xside, yside, zside, outputfGraphS2S2, outputfGeodS2S2);    // 2015.12.15, I changed here, Tiecheng

if(pnumnodes==1) fprintf(outputfNumnodesS2S2,"%d\n",(nAtomS2 / nsolvent2));
// fprintf(outputf,"\n[GRAPH SOLVENT2 SOLVENT2]: Done...\n");
}

if(gs3s3 == 1)
{
 sprintf(foutputGraphS3S3, "%s.%s.%s.Graph", argv[1], argv[4], argv[4]);
 outputfGraphS3S3 = fopen(foutputGraphS3S3,"w");

 sprintf(foutputGeodS3S3, "%s.%s.%s.GraphGeod", argv[1], argv[4], argv[4]);
 outputfGeodS3S3 = fopen(foutputGeodS3S3,"w");

 if(pnumnodes==1){
 sprintf(foutputNumnodesS3S3, "%s.%s.%s.GraphNumnodes", argv[1], argv[4], argv[4]);
 outputfNumnodesS3S3 = fopen(foutputNumnodesS3S3,"w"); }

 // Construct the solvent3-solvent3 graph

// graph_ss(atmS3, 1, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBD, s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6,
//          pbc, xside, yside, zside, outputfGraphS3S3, outputfGeodS3S3);
 graph_ss(atmS3, 1, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBDmin, s3as3bBDmax, s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6, s3s3v7,
          pbc, xside, yside, zside, outputfGraphS3S3, outputfGeodS3S3);  // 2015.12.15, I changed here, Tiecheng

 if(pnumnodes==1) fprintf(outputfNumnodesS3S3,"%d\n",(nAtomS3 / nsolvent3));
// fprintf(outputf,"\n[GRAPH SOLVENT3 SOLVENT3]: Done...\n");
}

if(gs1s2 == 1)
{
 sprintf(foutputGraphS1S2, "%s.%s.%s.Graph", argv[1], argv[2], argv[3]);
 outputfGraphS1S2 = fopen(foutputGraphS1S2,"w");

 sprintf(foutputGeodS1S2, "%s.%s.%s.GraphGeod", argv[1], argv[2], argv[3]);
 outputfGeodS1S2 = fopen(foutputGeodS1S2,"w");

 if(pnumnodes==1){
 sprintf(foutputNumnodesS1S2, "%s.%s.%s.GraphNumnodes", argv[1], argv[2], argv[3]);
 outputfNumnodesS1S2 = fopen(foutputNumnodesS1S2,"w"); }

 // Construct the solvent1-solvent2 graph
  
// graph_ss(atmS1, 1, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBD, s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6,
//          pbc, xside, yside, zside, outputfGraphS1S2, outputfGeodS1S2);
 graph_ss(atmS1, 1, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBDmin, s1as1bBDmax, s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6, s1s1v7,
          pbc, xside, yside, zside, outputfGraphS1S2, outputfGeodS1S2);  // 2015.12.15, I changed here, Tiecheng

 node1Start = 1 + (nAtomS1 / nsolvent1);     // 2015.12.15, the starting node id for solvent2 is right after solvent1, Tiecheng commented
// graph_ss(atmS2, node1Start, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBD, s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6,
//          pbc, xside, yside, zside, outputfGraphS1S2, outputfGeodS1S2);
 graph_ss(atmS2, node1Start, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBDmin, s2as2bBDmax, s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6, s2s2v7,
          pbc, xside, yside, zside, outputfGraphS1S2, outputfGeodS1S2);  // 2015.12.15, I changed here, Tiecheng

 node1Start = 1;
 node2Start = 1 + (nAtomS1 / nsolvent1);   // 2015.12.15, the starting node id for solvent2 is right after solvent1, Tiecheng commented

if((swatdip ==1 && watid == 1 && solid == 2) || (swatdip ==1 && watid == 2 && solid == 1))    // 2015.12.15, the dipole part, Tiecheng
{
  sprintf(foutputS1S2Dip, "%s.Dipole.angles", foutputGraphS1S2);
  outputfS1S2Dip = fopen(foutputS1S2Dip,"w"); 
 // graph_sAsB_dip(atmS1, atmS2, node1Start, node2Start, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
  graph_sAsB_dip(atmS1, atmS2, node1Start, node2Start, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBDmax,  // 2015.12.15, I only change the name here, the function is not updated, Tiecheng
                 s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6,
                 pbc, xside, yside, zside, opos, h1pos, h2pos, watid, solid, outputfGraphS1S2, outputfGeodS1S2, outputfS1S2Dip); 
 }
 else
// graph_sAsB(atmS1, atmS2, node1Start, node2Start, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
//            s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6,
//            pbc, xside, yside, zside, outputfGraphS1S2, outputfGeodS1S2);  
 graph_sAsB(atmS1, atmS2, node1Start, node2Start, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBDmin, s12as12bBDmax, // 2015.12.15, I changed here, Tiecheng
            s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, s1s2v7,
            pbc, xside, yside, zside, outputfGraphS1S2, outputfGeodS1S2);  

 if(pnumnodes==1) fprintf(outputfNumnodesS1S2,"%d\n",((nAtomS1 / nsolvent1) + (nAtomS2 / nsolvent2)));
// fprintf(outputf,"\n[GRAPH SOLVENT1 SOLVENT2]: Done...\n");
}

if(gs1s3 == 1)
{
 sprintf(foutputGraphS1S3, "%s.%s.%s.Graph", argv[1], argv[2], argv[4]);
 outputfGraphS1S3 = fopen(foutputGraphS1S3,"w");

 sprintf(foutputGeodS1S3, "%s.%s.%s.GraphGeod", argv[1], argv[2], argv[4]);
 outputfGeodS1S3 = fopen(foutputGeodS1S3,"w");

 if(pnumnodes==1){
 sprintf(foutputNumnodesS1S3, "%s.%s.%s.GraphNumnodes", argv[1], argv[2], argv[4]);
 outputfNumnodesS1S3 = fopen(foutputNumnodesS1S3,"w"); }

 // Construct the solvent1-solvent3 graph

// graph_ss(atmS1, 1, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBD, s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6,
//          pbc, xside, yside, zside, outputfGraphS1S3, outputfGeodS1S3);
 graph_ss(atmS1, 1, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBDmin, s1as1bBDmax, s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6, s1s1v7, // 2015.12.15, I changed here, Tiecheng
          pbc, xside, yside, zside, outputfGraphS1S3, outputfGeodS1S3);

 node1Start = 1 + (nAtomS1 / nsolvent1);
// graph_ss(atmS3, node1Start, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBD, s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6,
//          pbc, xside, yside, zside, outputfGraphS1S3, outputfGeodS1S3);
 graph_ss(atmS3, node1Start, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBDmin, s3as3bBDmax, s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6, s3s3v7, // 2015.12.15, I changed here, Tiecheng
          pbc, xside, yside, zside, outputfGraphS1S3, outputfGeodS1S3);

 node1Start = 1;
 node2Start = 1 + (nAtomS1 / nsolvent1);

if((swatdip ==1 && watid == 1 && solid == 3) || (swatdip ==1 && watid == 3 && solid == 1))
{
  sprintf(foutputS1S3Dip, "%s.Dipole.angles", foutputGraphS1S3);
  outputfS1S3Dip = fopen(foutputS1S3Dip,"w");
  graph_sAsB_dip(atmS1, atmS3, node1Start, node2Start, nsolvent1, nsolvent3, nAtomS1, nAtomS3, s1s3hbdn, s13a, s13b, s13as13bBDmax,  // 2015.12.15, I only change the name here, Tiecheng
                 s1s3hban, s1s3v1, s1s3v2, s1s3v3, s1s3v4, s1s3v5, s1s3v6,
                 pbc, xside, yside, zside, opos, h1pos, h2pos, watid, solid, outputfGraphS1S3, outputfGeodS1S3, outputfS1S3Dip);
 }
 else
// graph_sAsB(atmS1, atmS3, node1Start, node2Start, nsolvent1, nsolvent3, nAtomS1, nAtomS3, s1s3hbdn, s13a, s13b, s13as13bBD,
//            s1s3hban, s1s3v1, s1s3v2, s1s3v3, s1s3v4, s1s3v5, s1s3v6,
//            pbc, xside, yside, zside, outputfGraphS1S3, outputfGeodS1S3);
 graph_sAsB(atmS1, atmS3, node1Start, node2Start, nsolvent1, nsolvent3, nAtomS1, nAtomS3, s1s3hbdn, s13a, s13b, s13as13bBDmin, s13as13bBDmax,  // 2015.12.15, I changed here, Tiecheng
            s1s3hban, s1s3v1, s1s3v2, s1s3v3, s1s3v4, s1s3v5, s1s3v6, s1s3v7,
            pbc, xside, yside, zside, outputfGraphS1S3, outputfGeodS1S3);

 if(pnumnodes==1) fprintf(outputfNumnodesS1S3,"%d\n",((nAtomS1 / nsolvent1) + (nAtomS3 / nsolvent3)));
// fprintf(outputf,"\n[GRAPH SOLVENT1 SOLVENT3]: Done...\n");
}

if(gs2s3 == 1)
{
 sprintf(foutputGraphS2S3, "%s.%s.%s.Graph", argv[1], argv[3], argv[4]);
 outputfGraphS2S3 = fopen(foutputGraphS2S3,"w");

 sprintf(foutputGeodS2S3, "%s.%s.%s.GraphGeod", argv[1], argv[3], argv[4]);
 outputfGeodS2S3 = fopen(foutputGeodS2S3,"w");

 if(pnumnodes==1){
 sprintf(foutputNumnodesS2S3, "%s.%s.%s.GraphNumnodes", argv[1], argv[3], argv[4]);
 outputfNumnodesS2S3 = fopen(foutputNumnodesS2S3,"w"); }

 // Construct the solvent2-solvent3 graph

// graph_ss(atmS2, 1, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBD, s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6,
//          pbc, xside, yside, zside, outputfGraphS2S3, outputfGeodS2S3);
 graph_ss(atmS2, 1, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBDmin, s2as2bBDmax, s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6, s2s2v7,  // 2015.12.15, I changed here, Tiecheng
          pbc, xside, yside, zside, outputfGraphS2S3, outputfGeodS2S3);

 node1Start = 1 + (nAtomS2 / nsolvent2);
// graph_ss(atmS3, node1Start, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBD, s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6,
//         pbc, xside, yside, zside, outputfGraphS2S3, outputfGeodS2S3);
 graph_ss(atmS3, node1Start, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBDmin, s3as3bBDmax, s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6, s3s3v7, // 2015.12.15, I changed here, Tiecheng
          pbc, xside, yside, zside, outputfGraphS2S3, outputfGeodS2S3);

 node1Start = 1;
 node2Start = 1 + (nAtomS2 / nsolvent2);

if((swatdip ==1 && watid == 2 && solid == 3) || (swatdip ==1 && watid == 3 && solid == 2))
{
  sprintf(foutputS2S3Dip, "%s.Dipole.angles", foutputGraphS2S3);
  outputfS2S3Dip = fopen(foutputS2S3Dip,"w");
  graph_sAsB_dip(atmS2, atmS3, node1Start, node2Start, nsolvent2, nsolvent3, nAtomS2, nAtomS3, s2s3hbdn, s23a, s23b, s23as23bBDmax,   // 2015.12.15, I only change the name here
                 s2s3hban, s2s3v1, s2s3v2, s2s3v3, s2s3v4, s2s3v5, s2s3v6,
                 pbc, xside, yside, zside, opos, h1pos, h2pos, watid, solid, outputfGraphS2S3, outputfGeodS2S3, outputfS2S3Dip);
 }
 else
// graph_sAsB(atmS2, atmS3, node1Start, node2Start, nsolvent2, nsolvent3, nAtomS2, nAtomS3, s2s3hbdn, s23a, s23b, s23as23bBD,
//            s2s3hban, s2s3v1, s2s3v2, s2s3v3, s2s3v4, s2s3v5, s2s3v6,
//            pbc, xside, yside, zside, outputfGraphS2S3, outputfGeodS2S3);
 graph_sAsB(atmS2, atmS3, node1Start, node2Start, nsolvent2, nsolvent3, nAtomS2, nAtomS3, s2s3hbdn, s23a, s23b, s23as23bBDmin, s23as23bBDmax,  // 2015.12.15, I changed here, Tiecheng
            s2s3hban, s2s3v1, s2s3v2, s2s3v3, s2s3v4, s2s3v5, s2s3v6, s2s3v7,
            pbc, xside, yside, zside, outputfGraphS2S3, outputfGeodS2S3);

 if(pnumnodes==1) fprintf(outputfNumnodesS2S3,"%d\n",((nAtomS2 / nsolvent2) + (nAtomS3 / nsolvent3)));
// fprintf(outputf,"\n[GRAPH SOLVENT2 SOLVENT3]: Done...\n");
}


if(gs1t1 == 1)
{
 if(nsolventtype == 3 && (nsolutetype == 1 || nsolutetype == 2))
 {
    sprintf(foutputGraphS1T1, "%s.%s.%s.Graph", argv[1], argv[2], argv[5]);
    sprintf(foutputGeodS1T1, "%s.%s.%s.GraphGeod", argv[1], argv[2], argv[5]);
    sprintf(foutputNumnodesS1T1, "%s.%s.%s.GraphNumnodes", argv[1], argv[2], argv[5]);
 }
 if(nsolventtype == 2 && (nsolutetype == 1 || nsolutetype == 2))
 {
    sprintf(foutputGraphS1T1, "%s.%s.%s.Graph", argv[1], argv[2], argv[4]);
    sprintf(foutputGeodS1T1, "%s.%s.%s.GraphGeod", argv[1], argv[2], argv[4]);
    sprintf(foutputNumnodesS1T1, "%s.%s.%s.GraphNumnodes", argv[1], argv[2], argv[4]);
 }
 if(nsolventtype == 1 && (nsolutetype == 1 || nsolutetype == 2))
 {
    sprintf(foutputGraphS1T1, "%s.%s.%s.Graph", argv[1], argv[2], argv[3]);
    sprintf(foutputGeodS1T1, "%s.%s.%s.GraphGeod", argv[1], argv[2], argv[3]);
    sprintf(foutputNumnodesS1T1, "%s.%s.%s.GraphNumnodes", argv[1], argv[2], argv[3]);
 }

 outputfGraphS1T1 = fopen(foutputGraphS1T1,"w");
 outputfGeodS1T1 = fopen(foutputGeodS1T1,"w");
 if(pnumnodes==1) outputfNumnodesS1T1 = fopen(foutputNumnodesS1T1,"w");

 // Construct the solvent1-solute1 graph
  
// graph_ss(atmS1, 1, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBD, s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6,
//          pbc, xside, yside, zside, outputfGraphS1T1, outputfGeodS1T1);
 graph_ss(atmS1, 1, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBDmin, s1as1bBDmax, s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6, s1s1v7,  // 2015.12.15, I changed here, Tiecheng
          pbc, xside, yside, zside, outputfGraphS1T1, outputfGeodS1T1);

 node1Start = 1;
 node2Start = 1 + (nAtomS1 / nsolvent1);

 if(t1watdip == 1 && watid == 1)
 {
    sprintf(foutputS1T1Dip, "%s.Dipole.angles", foutputGraphS1T1);
    outputfS1T1Dip = fopen(foutputS1T1Dip,"w");
    graph_st_dip(atmS1, atmT1, node1Start, node2Start, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoffmax,   // 2015.12.15, I only change the name here, Tiecheng
                 pbc, xside, yside, zside, opos, h1pos, h2pos, outputfGraphS1T1, outputfGeodS1T1, outputfS1T1Dip);
 }
 else
//    graph_st(atmS1, atmT1, node1Start, node2Start, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff,
//            pbc, xside, yside, zside, outputfGraphS1T1, outputfGeodS1T1);  
    graph_st(atmS1, atmT1, node1Start, node2Start, nsolvent1, nsolute1, nAtomS1, nAtomT1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoffmin, s1t1cutoffmax, // 2015.12.15, I changed here, Tiecheng
             pbc, xside, yside, zside, outputfGraphS1T1, outputfGeodS1T1);  

 if(pnumnodes==1) fprintf(outputfNumnodesS1T1,"%d\n",((nAtomS1 / nsolvent1) + (nAtomT1 / nsolute1)));
// fprintf(outputf,"\n[GRAPH SOLVENT1 SOLUTE1]: Done...\n");
}

if(gs1t2 == 1)
{
 if(nsolventtype == 3)
 {
    sprintf(foutputGraphS1T2, "%s.%s.%s.Graph", argv[1], argv[2], argv[6]);
    sprintf(foutputGeodS1T2, "%s.%s.%s.GraphGeod", argv[1], argv[2], argv[6]);
    sprintf(foutputNumnodesS1T2, "%s.%s.%s.GraphNumnodes", argv[1], argv[2], argv[6]);
 }
 if(nsolventtype == 2)
 {
    sprintf(foutputGraphS1T2, "%s.%s.%s.Graph", argv[1], argv[2], argv[5]);
    sprintf(foutputGeodS1T2, "%s.%s.%s.GraphGeod", argv[1], argv[2], argv[5]);
    sprintf(foutputNumnodesS1T2, "%s.%s.%s.GraphNumnodes", argv[1], argv[2], argv[5]);
 }
 if(nsolventtype == 1)
 {
    sprintf(foutputGraphS1T2, "%s.%s.%s.Graph", argv[1], argv[2], argv[4]);
    sprintf(foutputGeodS1T2, "%s.%s.%s.GraphGeod", argv[1], argv[2], argv[4]);
    sprintf(foutputNumnodesS1T2, "%s.%s.%s.GraphNumnodes", argv[1], argv[2], argv[4]);
 }

 outputfGraphS1T2 = fopen(foutputGraphS1T2,"w");
 outputfGeodS1T2 = fopen(foutputGeodS1T2,"w");
 if(pnumnodes==1) outputfNumnodesS1T2 = fopen(foutputNumnodesS1T2,"w");

 // Construct the solvent1-solute2 graph
  
// graph_ss(atmS1, 1, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBD, s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6,
//          pbc, xside, yside, zside, outputfGraphS1T2, outputfGeodS1T2);
 graph_ss(atmS1, 1, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBDmin, s1as1bBDmax, s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6, s1s1v7,  // 2015.12.15, I changed here, Tiecheng
          pbc, xside, yside, zside, outputfGraphS1T2, outputfGeodS1T2);

 node1Start = 1;
 node2Start = 1 + (nAtomS1 / nsolvent1);

 if(t2watdip == 1 && watid == 1)
 {
    sprintf(foutputS1T2Dip, "%s.Dipole.angles", foutputGraphS1T2);
    outputfS1T2Dip = fopen(foutputS1T2Dip,"w");
    graph_st_dip(atmS1, atmT2, node1Start, node2Start, nsolvent1, nsolute2, nAtomS1, nAtomT2, s1t2cutoffnum, s1t2a, s1t2b, s1t2cutoffmax,  // 2015.12.15, I only changed the name here, Tiecheng
                 pbc, xside, yside, zside, opos, h1pos, h2pos, outputfGraphS1T2, outputfGeodS1T2, outputfS1T2Dip);

 }
 else
//    graph_st(atmS1, atmT2, node1Start, node2Start, nsolvent1, nsolute2, nAtomS1, nAtomT2, s1t2cutoffnum, s1t2a, s1t2b, s1t2cutoff,
//             pbc, xside, yside, zside, outputfGraphS1T2, outputfGeodS1T2);  
    graph_st(atmS1, atmT2, node1Start, node2Start, nsolvent1, nsolute2, nAtomS1, nAtomT2, s1t2cutoffnum, s1t2a, s1t2b, s1t2cutoffmin, s1t2cutoffmax,  // 2015.12.15, I changed here, Tiecheng
             pbc, xside, yside, zside, outputfGraphS1T2, outputfGeodS1T2);  

 if(pnumnodes==1) fprintf(outputfNumnodesS1T2,"%d\n",((nAtomS1 / nsolvent1) + (nAtomT2 / nsolute2)));
// fprintf(outputf,"\n[GRAPH SOLVENT1 SOLUTE2]: Done...\n");
}

if(gs2t1 == 1)
{
 if(nsolventtype == 3 && (nsolutetype == 1 || nsolutetype == 2))
 {
    sprintf(foutputGraphS2T1, "%s.%s.%s.Graph", argv[1], argv[3], argv[5]);
    sprintf(foutputGeodS2T1, "%s.%s.%s.GraphGeod", argv[1], argv[3], argv[5]);
    sprintf(foutputNumnodesS2T1, "%s.%s.%s.GraphNumnodes", argv[1], argv[3], argv[5]);
 }
 if(nsolventtype == 2 && (nsolutetype == 1 || nsolutetype == 2))
 {
    sprintf(foutputGraphS2T1, "%s.%s.%s.Graph", argv[1], argv[3], argv[4]);
    sprintf(foutputGeodS2T1, "%s.%s.%s.GraphGeod", argv[1], argv[3], argv[4]);
    sprintf(foutputNumnodesS2T1, "%s.%s.%s.GraphNumnodes", argv[1], argv[3], argv[4]);
 }

 outputfGraphS2T1 = fopen(foutputGraphS2T1,"w");
 outputfGeodS2T1 = fopen(foutputGeodS2T1,"w");
 if(pnumnodes==1) outputfNumnodesS2T1 = fopen(foutputNumnodesS2T1,"w");

 // Construct the solvent2-solute1 graph
  
// graph_ss(atmS2, 1, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBD, s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6,
//          pbc, xside, yside, zside, outputfGraphS2T1, outputfGeodS2T1);
 graph_ss(atmS2, 1, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBDmin, s2as2bBDmax, s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6, s2s2v7,  // 2015.12.15, I changed here, Tiecheng
          pbc, xside, yside, zside, outputfGraphS2T1, outputfGeodS2T1);

 node1Start = 1;
 node2Start = 1 + (nAtomS2 / nsolvent2);

 if(t1watdip == 1 && watid == 2)
 {
    sprintf(foutputS2T1Dip, "%s.Dipole.angles", foutputGraphS2T1);
    outputfS2T1Dip = fopen(foutputS2T1Dip,"w");
    graph_st_dip(atmS2, atmT1, node1Start, node2Start, nsolvent2, nsolute1, nAtomS2, nAtomT1, s2t1cutoffnum, s2t1a, s2t1b, s2t1cutoffmax,  // 2015.12.15, I only changed the name here, Tiecheng
                 pbc, xside, yside, zside, opos, h1pos, h2pos, outputfGraphS2T1, outputfGeodS2T1, outputfS2T1Dip);
 }
 else
//    graph_st(atmS2, atmT1, node1Start, node2Start, nsolvent2, nsolute1, nAtomS2, nAtomT1, s2t1cutoffnum, s2t1a, s2t1b, s2t1cutoff,
//             pbc, xside, yside, zside, outputfGraphS2T1, outputfGeodS2T1);  
    graph_st(atmS2, atmT1, node1Start, node2Start, nsolvent2, nsolute1, nAtomS2, nAtomT1, s2t1cutoffnum, s2t1a, s2t1b, s2t1cutoffmin, s2t1cutoffmax,  // 2015.12.15, I changed here, Tiecheng
             pbc, xside, yside, zside, outputfGraphS2T1, outputfGeodS2T1);  

 if(pnumnodes==1) fprintf(outputfNumnodesS2T1,"%d\n",((nAtomS2 / nsolvent2) + (nAtomT1 / nsolute1)));
// fprintf(outputf,"\n[GRAPH SOLVENT2 SOLUTE1]: Done...\n");
}

if(gs2t2 == 1)
{
 if(nsolventtype == 3)
 {
    sprintf(foutputGraphS2T2, "%s.%s.%s.Graph", argv[1], argv[3], argv[6]);
    sprintf(foutputGeodS2T2, "%s.%s.%s.GraphGeod", argv[1], argv[3], argv[6]);
    sprintf(foutputNumnodesS2T2, "%s.%s.%s.GraphNumnodes", argv[1], argv[3], argv[6]);
 }
 if(nsolventtype == 2)
 {
    sprintf(foutputGraphS2T2, "%s.%s.%s.Graph", argv[1], argv[3], argv[5]);
    sprintf(foutputGeodS2T2, "%s.%s.%s.GraphGeod", argv[1], argv[3], argv[5]);
    sprintf(foutputNumnodesS2T2, "%s.%s.%s.GraphNumnodes", argv[1], argv[3], argv[5]);
 }

 outputfGraphS2T2 = fopen(foutputGraphS2T2,"w");
 outputfGeodS2T2 = fopen(foutputGeodS2T2,"w");
 if(pnumnodes==1) outputfNumnodesS2T2 = fopen(foutputNumnodesS2T2,"w");

 // Construct the solvent2-solute2 graph
  
// graph_ss(atmS2, 1, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBD, s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6,
//          pbc, xside, yside, zside, outputfGraphS2T2, outputfGeodS2T2);
 graph_ss(atmS2, 1, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBDmin, s2as2bBDmax, s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6, s2s2v7,   // 2015.12.15, I changed here, Tiecheng
          pbc, xside, yside, zside, outputfGraphS2T2, outputfGeodS2T2);

 node1Start = 1;
 node2Start = 1 + (nAtomS2 / nsolvent2);

 if(t2watdip == 1 && watid == 2)
 {
    sprintf(foutputS2T2Dip, "%s.Dipole.angles", foutputGraphS2T2);
    outputfS2T2Dip = fopen(foutputS2T2Dip,"w");
    graph_st_dip(atmS2, atmT2, node1Start, node2Start, nsolvent2, nsolute2, nAtomS2, nAtomT2, s2t2cutoffnum, s2t2a, s2t2b, s2t2cutoffmax,  // 2015.12.15, I only changed the name, Tiecheng
                 pbc, xside, yside, zside, opos, h1pos, h2pos, outputfGraphS2T2, outputfGeodS2T2, outputfS2T2Dip);
 }
 else 
//    graph_st(atmS2, atmT2, node1Start, node2Start, nsolvent2, nsolute2, nAtomS2, nAtomT2, s2t2cutoffnum, s2t2a, s2t2b, s2t2cutoff,
//             pbc, xside, yside, zside, outputfGraphS2T2, outputfGeodS2T2);  
    graph_st(atmS2, atmT2, node1Start, node2Start, nsolvent2, nsolute2, nAtomS2, nAtomT2, s2t2cutoffnum, s2t2a, s2t2b, s2t2cutoffmin, s2t2cutoffmax,  // 2015.12.15, I changed here, Tiecheng
             pbc, xside, yside, zside, outputfGraphS2T2, outputfGeodS2T2);  

 if(pnumnodes==1) fprintf(outputfNumnodesS2T2,"%d\n",((nAtomS2 / nsolvent2) + (nAtomT2 / nsolute2)));
// fprintf(outputf,"\n[GRAPH SOLVENT2 SOLUTE2]: Done...\n");
}

if(gs3t1 == 1)
{
 if(nsolutetype == 1 || nsolutetype == 2)
 {
    sprintf(foutputGraphS3T1, "%s.%s.%s.Graph", argv[1], argv[4], argv[5]);
    sprintf(foutputGeodS3T1, "%s.%s.%s.GraphGeod", argv[1], argv[4], argv[5]);
    sprintf(foutputNumnodesS3T1, "%s.%s.%s.GraphNumnodes", argv[1], argv[4], argv[5]);
 }

 outputfGraphS3T1 = fopen(foutputGraphS3T1,"w");
 outputfGeodS3T1 = fopen(foutputGeodS3T1,"w");
 if(pnumnodes==1) outputfNumnodesS3T1 = fopen(foutputNumnodesS3T1,"w");

 // Construct the solvent3-solute1 graph
  
// graph_ss(atmS3, 1, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBD, s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6,
//          pbc, xside, yside, zside, outputfGraphS3T1, outputfGeodS3T1);
 graph_ss(atmS3, 1, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBDmin, s3as3bBDmax, s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6, s3s3v7,  // 2015.12.15, I changed here, Tiecheng
          pbc, xside, yside, zside, outputfGraphS3T1, outputfGeodS3T1);

 node1Start = 1;
 node2Start = 1 + (nAtomS3 / nsolvent3);

 if(t1watdip == 1 && watid == 3)
 {
    sprintf(foutputS3T1Dip, "%s.Dipole.angles", foutputGraphS3T1);
    outputfS3T1Dip = fopen(foutputS3T1Dip,"w");
    graph_st_dip(atmS3, atmT1, node1Start, node2Start, nsolvent3, nsolute1, nAtomS3, nAtomT1, s3t1cutoffnum, s3t1a, s3t1b, s3t1cutoffmax,  // 2015.12.15, I only changed the name here, Tiecheng
                 pbc, xside, yside, zside, opos, h1pos, h2pos, outputfGraphS3T1, outputfGeodS3T1, outputfS3T1Dip);
 }
 else
//    graph_st(atmS3, atmT1, node1Start, node2Start, nsolvent3, nsolute1, nAtomS3, nAtomT1, s3t1cutoffnum, s3t1a, s3t1b, s3t1cutoff,
//             pbc, xside, yside, zside, outputfGraphS3T1, outputfGeodS3T1);  
    graph_st(atmS3, atmT1, node1Start, node2Start, nsolvent3, nsolute1, nAtomS3, nAtomT1, s3t1cutoffnum, s3t1a, s3t1b, s3t1cutoffmin, s3t1cutoffmax,  // 2015.12.15, I changed here, Tiecheng
             pbc, xside, yside, zside, outputfGraphS3T1, outputfGeodS3T1);  

 if(pnumnodes==1) fprintf(outputfNumnodesS3T1,"%d\n",((nAtomS3 / nsolvent3) + (nAtomT1 / nsolute1)));
// fprintf(outputf,"\n[GRAPH SOLVENT3 SOLUTE1]: Done...\n");
}

if(gs3t2 == 1)
{
  sprintf(foutputGraphS3T2, "%s.%s.%s.Graph", argv[1], argv[4], argv[6]);
  sprintf(foutputGeodS3T2, "%s.%s.%s.GraphGeod", argv[1], argv[4], argv[6]);
  sprintf(foutputNumnodesS3T2, "%s.%s.%s.GraphNumnodes", argv[1], argv[4], argv[6]);

  outputfGraphS3T2 = fopen(foutputGraphS3T2,"w");
  outputfGeodS3T2 = fopen(foutputGeodS3T2,"w");
 if(pnumnodes==1)  outputfNumnodesS3T2 = fopen(foutputNumnodesS3T2,"w");

 // Construct the solvent3-solute2 graph

// graph_ss(atmS3, 1, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBD, s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6,
//          pbc, xside, yside, zside, outputfGraphS3T2, outputfGeodS3T2);
 graph_ss(atmS3, 1, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBDmin, s3as3bBDmax, s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6, s3s3v7,  // 2015.12.15, I changed here, Tiecheng
          pbc, xside, yside, zside, outputfGraphS3T2, outputfGeodS3T2);

 node1Start = 1;
 node2Start = 1 + (nAtomS3 / nsolvent3);

 if(t2watdip == 1 && watid == 3)
 {
    sprintf(foutputS3T2Dip, "%s.Dipole.angles", foutputGraphS3T2);
    outputfS3T2Dip = fopen(foutputS3T2Dip,"w");
    graph_st_dip(atmS3, atmT2, node1Start, node2Start, nsolvent3, nsolute2, nAtomS3, nAtomT2, s3t2cutoffnum, s3t2a, s3t2b, s3t2cutoffmax,   // 2015.12.15, I only changed the name here, Tiecheng
                 pbc, xside, yside, zside, opos, h1pos, h2pos, outputfGraphS3T2, outputfGeodS3T2, outputfS3T2Dip);
 }
 else
//    graph_st(atmS3, atmT2, node1Start, node2Start, nsolvent3, nsolute2, nAtomS3, nAtomT2, s3t2cutoffnum, s3t2a, s3t2b, s3t2cutoff,
//             pbc, xside, yside, zside, outputfGraphS3T2, outputfGeodS3T2);
    graph_st(atmS3, atmT2, node1Start, node2Start, nsolvent3, nsolute2, nAtomS3, nAtomT2, s3t2cutoffnum, s3t2a, s3t2b, s3t2cutoffmin, s3t2cutoffmax,  //2015.12.15, I changed here Tiecheng
             pbc, xside, yside, zside, outputfGraphS3T2, outputfGeodS3T2);

 if(pnumnodes==1) fprintf(outputfNumnodesS3T2,"%d\n",((nAtomS3 / nsolvent3) + (nAtomT2 / nsolute2)));
// fprintf(outputf,"\n[GRAPH SOLVENT3 SOLUTE2]: Done...\n");
}


if(gs1s2s3 == 1)
{
 sprintf(foutputGraphS1S2S3, "%s.%s.%s.%s.Graph", argv[1], argv[2], argv[3], argv[4]);
 outputfGraphS1S2S3 = fopen(foutputGraphS1S2S3,"w");

 sprintf(foutputGeodS1S2S3, "%s.%s.%s.%s.GraphGeod", argv[1], argv[2], argv[3], argv[4]);
 outputfGeodS1S2S3 = fopen(foutputGeodS1S2S3,"w");

 if(pnumnodes==1){
 sprintf(foutputNumnodesS1S2S3, "%s.%s.%s.%s.GraphNumnodes", argv[1], argv[2], argv[3], argv[4]);
 outputfNumnodesS1S2S3 = fopen(foutputNumnodesS1S2S3,"w"); }

 // Construct the solvent1-solvent2-solvent3 graph
  
// graph_ss(atmS1, 1, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBD, s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6,
//          pbc, xside, yside, zside, outputfGraphS1S2S3, outputfGeodS1S2S3);
 graph_ss(atmS1, 1, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBDmin, s1as1bBDmax, s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6, s1s1v7,  //2015.12.15, I changed here, Tiecheng
          pbc, xside, yside, zside, outputfGraphS1S2S3, outputfGeodS1S2S3);

 node1Start = 1 + (nAtomS1 / nsolvent1);
// graph_ss(atmS2, node1Start, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBD, s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6,
//          pbc, xside, yside, zside, outputfGraphS1S2S3, outputfGeodS1S2S3);
 graph_ss(atmS2, node1Start, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBDmin, s2as2bBDmax, s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6, s2s2v7,   //2015.12.15, I changed here, Tiecheng
          pbc, xside, yside, zside, outputfGraphS1S2S3, outputfGeodS1S2S3);

 node1Start = 1;
 node2Start = 1 + (nAtomS1 / nsolvent1); 
// graph_sAsB(atmS1, atmS2, node1Start, node2Start, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBD,
//            s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6,
//            pbc, xside, yside, zside, outputfGraphS1S2S3, outputfGeodS1S2S3);
 graph_sAsB(atmS1, atmS2, node1Start, node2Start, nsolvent1, nsolvent2, nAtomS1, nAtomS2, s1s2hbdn, s12a, s12b, s12as12bBDmin, s12as12bBDmax,    // 2015.12.15, I changed here, Tiecheng
            s1s2hban, s1s2v1, s1s2v2, s1s2v3, s1s2v4, s1s2v5, s1s2v6, s1s2v7, 
            pbc, xside, yside, zside, outputfGraphS1S2S3, outputfGeodS1S2S3);

 node1Start = 1 + (nAtomS1 / nsolvent1) + (nAtomS2 / nsolvent2);
// graph_ss(atmS3, node1Start, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBD, s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6,
//          pbc, xside, yside, zside, outputfGraphS1S2S3, outputfGeodS1S2S3);
 graph_ss(atmS3, node1Start, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBDmin, s3as3bBDmax, s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6, s3s3v7,  // 2015.12.15, I changed here, Tiecheng
          pbc, xside, yside, zside, outputfGraphS1S2S3, outputfGeodS1S2S3);

 node1Start = 1;
 node2Start = 1 + (nAtomS1 / nsolvent1) + (nAtomS2 / nsolvent2);
// graph_sAsB(atmS1, atmS3, node1Start, node2Start, nsolvent1, nsolvent3, nAtomS1, nAtomS3, s1s3hbdn, s13a, s13b, s13as13bBD,
//            s1s3hban, s1s3v1, s1s3v2, s1s3v3, s1s3v4, s1s3v5, s1s3v6,
//            pbc, xside, yside, zside, outputfGraphS1S2S3, outputfGeodS1S2S3);
 graph_sAsB(atmS1, atmS3, node1Start, node2Start, nsolvent1, nsolvent3, nAtomS1, nAtomS3, s1s3hbdn, s13a, s13b, s13as13bBDmin, s13as13bBDmax,   // 2015.12.15, I changed here, Tiecheng
            s1s3hban, s1s3v1, s1s3v2, s1s3v3, s1s3v4, s1s3v5, s1s3v6, s1s3v7,
            pbc, xside, yside, zside, outputfGraphS1S2S3, outputfGeodS1S2S3);

 node1Start = 1 + (nAtomS1 / nsolvent1);
 node2Start = 1 + (nAtomS1 / nsolvent1) + (nAtomS2 / nsolvent2);
// graph_sAsB(atmS2, atmS3, node1Start, node2Start, nsolvent2, nsolvent3, nAtomS2, nAtomS3, s2s3hbdn, s23a, s23b, s23as23bBD,
//            s2s3hban, s2s3v1, s2s3v2, s2s3v3, s2s3v4, s2s3v5, s2s3v6,
//            pbc, xside, yside, zside, outputfGraphS1S2S3, outputfGeodS1S2S3);
 graph_sAsB(atmS2, atmS3, node1Start, node2Start, nsolvent2, nsolvent3, nAtomS2, nAtomS3, s2s3hbdn, s23a, s23b, s23as23bBDmin, s23as23bBDmax,   // 2015.12.15, I changed here, Tiecheng
            s2s3hban, s2s3v1, s2s3v2, s2s3v3, s2s3v4, s2s3v5, s2s3v6, s2s3v7,
            pbc, xside, yside, zside, outputfGraphS1S2S3, outputfGeodS1S2S3);
  
 if(pnumnodes==1) fprintf(outputfNumnodesS1S2S3,"%d\n",((nAtomS1 / nsolvent1) + (nAtomS2 / nsolvent2) + (nAtomS3 / nsolvent3)));
// fprintf(outputf,"\n[GRAPH SOLVENT1 SOLVENT2 SOLVENT3]: Done...\n");
}




// 2015.12.15, I did not change the water structures part, what I did is just replace distance terms with an upper boundary max

// Search for water structures hexamers - pentamers - tetramers - trimers


opos = 0;

for(i = 0; i < 3; i++)
{
  if(watid == 1){ if(strcmp(slvntatm1[i], "O") == 0 || strcmp(slvntatm1[i], "OW") == 0) opos = slvntatmnum1[i]; }
  if(watid == 2){ if(strcmp(slvntatm2[i], "O") == 0 || strcmp(slvntatm2[i], "OW") == 0) opos = slvntatmnum2[i]; }
  if(watid == 3){ if(strcmp(slvntatm3[i], "O") == 0 || strcmp(slvntatm3[i], "OW") == 0) opos = slvntatmnum3[i]; }
}

// Hexamers

if(ws == 1 && hring == 1)
{
 if(opos == 0) {printf("Water structures: Position of water oxygen ( O ) cannot be read\n"); exit(-1);} 
 if(watid == 1) sprintf(foutputHRing, "%s.%s.Water.Hexamer.Rings", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputHRing, "%s.%s.Water.Hexamer.Rings", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputHRing, "%s.%s.Water.Hexamer.Rings", argv[1], argv[4]);
 outputfHRing = fopen(foutputHRing,"w");
 if(watid == 1) sprintf(foutputHRingIso, "%s.%s.Water.Hexamer.Rings.isolated", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputHRingIso, "%s.%s.Water.Hexamer.Rings.isolated", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputHRingIso, "%s.%s.Water.Hexamer.Rings.isolated", argv[1], argv[4]);

 if(watid == 1) 
// hringsearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBD, 
 hringsearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBDmax,   // 2015.12.15, I changed the name, Tiecheng
             s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6, hiso, 
             outputfHRing, foutputHRingIso);

 if(watid == 2)
// hringsearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBD,
 hringsearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBDmax,  // 2015.12.15, I changed the name, Tiecheng
             s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6, hiso,
             outputfHRing, foutputHRingIso);

 if(watid == 3)
 //hringsearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBD,
 hringsearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBDmax,  // 2015.12.15, I changed the name, Tiecheng
             s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6, hiso,
             outputfHRing, foutputHRingIso);

}

if(ws == 1 && hbook == 1)
{
 if(opos == 0) {printf("Water structures: Position of water oxygen ( O ) cannot be read\n"); exit(-1);}
 if(watid == 1) sprintf(foutputHBook, "%s.%s.Water.Hexamer.Books", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputHBook, "%s.%s.Water.Hexamer.Books", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputHBook, "%s.%s.Water.Hexamer.Books", argv[1], argv[4]);
 outputfHBook = fopen(foutputHBook,"w");
 if(watid == 1) sprintf(foutputHBookIso, "%s.%s.Water.Hexamer.Books.isolated", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputHBookIso, "%s.%s.Water.Hexamer.Books.isolated", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputHBookIso, "%s.%s.Water.Hexamer.Books.isolated", argv[1], argv[4]);

 if(watid == 1)
 //hbooksearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBD,
 hbooksearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBDmax,   // 2015.12.15, I changed the name here, Tiecheng
             s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6, hiso,
             outputfHBook, foutputHBookIso);

 if(watid == 2)
// hbooksearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBD,
 hbooksearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBDmax,  // 2015.12.15, I changed the name here, Tiecheng
             s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6, hiso,
             outputfHBook, foutputHBookIso);

 if(watid == 3)
 //hbooksearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBD,
 hbooksearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBDmax,   // 2015.12.15, I changed the name here, Tiecheng
             s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6, hiso,
             outputfHBook, foutputHBookIso);

}

if(ws == 1 && hprism == 1)
{
 if(opos == 0) {printf("Water structures: Position of water oxygen ( O ) cannot be read\n"); exit(-1);}
 if(watid == 1) sprintf(foutputHPrism, "%s.%s.Water.Hexamer.Prisms", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputHPrism, "%s.%s.Water.Hexamer.Prisms", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputHPrism, "%s.%s.Water.Hexamer.Prisms", argv[1], argv[4]);
 outputfHPrism = fopen(foutputHPrism,"w");
 if(watid == 1) sprintf(foutputHPrismIso, "%s.%s.Water.Hexamer.Prisms.isolated", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputHPrismIso, "%s.%s.Water.Hexamer.Prisms.isolated", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputHPrismIso, "%s.%s.Water.Hexamer.Prisms.isolated", argv[1], argv[4]);

 if(watid == 1)
 //hprismsearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBD,
 hprismsearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBDmax,   // 2015.12.15, I changed the name here, Tiecheng
              s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6, hiso,
              outputfHPrism, foutputHPrismIso);

 if(watid == 2)
 //hprismsearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBD,
 hprismsearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBDmax,   // 2015.12.15,  I changed the name here, Tiecheng
              s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6, hiso,
              outputfHPrism, foutputHPrismIso);

 if(watid == 3)
// hprismsearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBD,
//             s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6, hiso,
//              outputfHPrism, foutputHPrismIso);
 hprismsearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBDmax,    // 2015.12.15, I changed the name here, Tiecheng
              s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6, hiso,
              outputfHPrism, foutputHPrismIso);

}

if(ws == 1 && hcage == 1)
{
 if(opos == 0) {printf("Water structures: Position of water oxygen ( O ) cannot be read\n"); exit(-1);}
 if(watid == 1) sprintf(foutputHCage, "%s.%s.Water.Hexamer.Cages", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputHCage, "%s.%s.Water.Hexamer.Cages", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputHCage, "%s.%s.Water.Hexamer.Cages", argv[1], argv[4]);
 outputfHCage = fopen(foutputHCage,"w");
 if(watid == 1) sprintf(foutputHCageIso, "%s.%s.Water.Hexamer.Cages.isolated", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputHCageIso, "%s.%s.Water.Hexamer.Cages.isolated", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputHCageIso, "%s.%s.Water.Hexamer.Cages.isolated", argv[1], argv[4]);

 if(watid == 1)
//hcagesearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBD,
 hcagesearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBDmax,   // 2015.12.16, I changed the name here, Tiecheng
             s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6, hiso,
             outputfHCage, foutputHCageIso);

 if(watid == 2)
// hcagesearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBD,
 hcagesearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBDmax,    // 2015.12.16, I changed the name here, Tiecheng
             s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6, hiso,
             outputfHCage, foutputHCageIso);

 if(watid == 3)
// hcagesearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBD,
 hcagesearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBDmax,   // 2015.12.16, I changed the name here, Tiecheng
             s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6, hiso,
             outputfHCage, foutputHCageIso);

}

if(ws == 1 && hbag == 1)
{
 if(opos == 0) {printf("Water structures: Position of water oxygen ( O ) cannot be read\n"); exit(-1);}
 if(watid == 1) sprintf(foutputHBag, "%s.%s.Water.Hexamer.Bags", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputHBag, "%s.%s.Water.Hexamer.Bags", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputHBag, "%s.%s.Water.Hexamer.Bags", argv[1], argv[4]);
 outputfHBag = fopen(foutputHBag,"w");
 if(watid == 1) sprintf(foutputHBagIso, "%s.%s.Water.Hexamer.Bags.isolated", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputHBagIso, "%s.%s.Water.Hexamer.Bags.isolated", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputHBagIso, "%s.%s.Water.Hexamer.Bags.isolated", argv[1], argv[4]);

 if(watid == 1)
// hbagsearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBD,
 hbagsearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBDmax,    // 2015.12.16, I changed the name here, Tiecheng
            s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6, hiso,
            outputfHBag, foutputHBagIso);

 if(watid == 2)
// hbagsearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBD,
 hbagsearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBDmax,    // 2015.12.16, I changed the name here, Tiecheng
            s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6, hiso,
            outputfHBag, foutputHBagIso);

 if(watid == 3)
// hbagsearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBD,
 hbagsearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBDmax,    // 2015.12.16, I changed the name here, Tiecheng
            s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6, hiso,
            outputfHBag, foutputHBagIso);

}

if(ws == 1 && hboat == 1)
{
 if(opos == 0) {printf("Water structures: Position of water oxygen ( O ) cannot be read\n"); exit(-1);}
 if(watid == 1) sprintf(foutputHBoat, "%s.%s.Water.Hexamer.Boats", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputHBoat, "%s.%s.Water.Hexamer.Boats", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputHBoat, "%s.%s.Water.Hexamer.Boats", argv[1], argv[4]);
 outputfHBoat = fopen(foutputHBoat,"w");
 if(watid == 1) sprintf(foutputHBoatIso, "%s.%s.Water.Hexamer.Boats.isolated", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputHBoatIso, "%s.%s.Water.Hexamer.Boats.isolated", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputHBoatIso, "%s.%s.Water.Hexamer.Boats.isolated", argv[1], argv[4]);

 if(watid == 1)
// hboatsearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBD,
 hboatsearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBDmax,    // 2015.12.16, I changed the name here, Tiecheng
             s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6, hiso,
             outputfHBoat, foutputHBoatIso);

 if(watid == 2)
// hboatsearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBD,
 hboatsearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBDmax,    // 2015.12.16, I changed the name here, Tiecheng
             s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6, hiso,
             outputfHBoat, foutputHBoatIso);

 if(watid == 3)
// hboatsearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBD,
 hboatsearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBDmax,     // 2015.12.16, I changed the name here, Tiecheng
             s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6, hiso,
             outputfHBoat, foutputHBoatIso);

}

if(ws == 1 && hchair == 1)
{
 if(opos == 0) {printf("Water structures: Position of water oxygen ( O ) cannot be read\n"); exit(-1);}
 if(watid == 1) sprintf(foutputHChair, "%s.%s.Water.Hexamer.Chairs", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputHChair, "%s.%s.Water.Hexamer.Chairs", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputHChair, "%s.%s.Water.Hexamer.Chairs", argv[1], argv[4]);
 outputfHChair = fopen(foutputHChair,"w");
 if(watid == 1) sprintf(foutputHChairIso, "%s.%s.Water.Hexamer.Chairs.isolated", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputHChairIso, "%s.%s.Water.Hexamer.Chairs.isolated", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputHChairIso, "%s.%s.Water.Hexamer.Chairs.isolated", argv[1], argv[4]);

 if(watid == 1)
// hchairsearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBD,
 hchairsearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBDmax,     // 2015.12.16, I changed the name here, Tiecheng
              s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6, hiso,
              outputfHChair, foutputHChairIso);

 if(watid == 2)
// hchairsearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBD,
 hchairsearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBDmax,    // 2015.12.16, I changed the name here, Tiecheng
              s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6, hiso,
              outputfHChair, foutputHChairIso);

 if(watid == 3)
// hchairsearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBD,
 hchairsearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBDmax,    // 2015.12.16, I changed the name here, Tiecheng
              s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6, hiso,
              outputfHChair, foutputHChairIso);

}

if(ws == 1 && hprismbook == 1)
{
 if(opos == 0) {printf("Water structures: Position of water oxygen ( O ) cannot be read\n"); exit(-1);}
 if(watid == 1) sprintf(foutputHPrismBook, "%s.%s.Water.Hexamer.PrismBooks", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputHPrismBook, "%s.%s.Water.Hexamer.PrismBooks", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputHPrismBook, "%s.%s.Water.Hexamer.PrismBooks", argv[1], argv[4]);
 outputfHPrismBook = fopen(foutputHPrismBook,"w");
 if(watid == 1) sprintf(foutputHPrismBookIso, "%s.%s.Water.Hexamer.PrismBooks.isolated", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputHPrismBookIso, "%s.%s.Water.Hexamer.PrismBooks.isolated", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputHPrismBookIso, "%s.%s.Water.Hexamer.PrismBooks.isolated", argv[1], argv[4]);

 if(watid == 1)
// hprismbooksearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBD,
 hprismbooksearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBDmax,    // 2015.12.16, I changed the name, Tiecheng
                  s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6, hiso,
                  outputfHPrismBook, foutputHPrismBookIso);

 if(watid == 2)
// hprismbooksearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBD,
 hprismbooksearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBDmax,    // 2015.12.16, I changed the name, Tiecheng
                  s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6, hiso,
                  outputfHPrismBook, foutputHPrismBookIso);

 if(watid == 3)
// hprismbooksearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBD,
 hprismbooksearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBDmax,     // 2015.12.16, I changed the name, Tiecheng
                  s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6, hiso,
                  outputfHPrismBook, foutputHPrismBookIso);

}

// Water Cyclic Pentmer search

if(ws == 1 && hpentamer == 1)
{
 if(opos == 0) {printf("Water structures: Position of water oxygen ( O ) cannot be read\n"); exit(-1);}
 if(watid == 1) sprintf(foutputPentamer, "%s.%s.Water.Cyclic.Pentamers", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputPentamer, "%s.%s.Water.Cyclic.Pentamers", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputPentamer, "%s.%s.Water.Cyclic.Pentamers", argv[1], argv[4]);
 outputfPentamer = fopen(foutputPentamer,"w");
 if(watid == 1) sprintf(foutputPentamerIso, "%s.%s.Water.Cyclic.Pentamers.isolated", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputPentamerIso, "%s.%s.Water.Cyclic.Pentamers.isolated", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputPentamerIso, "%s.%s.Water.Cyclic.Pentamers.isolated", argv[1], argv[4]);

 if(watid == 1)
// hpentamersearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBD,
 hpentamersearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBDmax,     // 2015.12.16, I changed the name, Tiecheng
                 s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6, hiso,
                 outputfPentamer, foutputPentamerIso);

 if(watid == 2)
// hpentamersearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBD,
 hpentamersearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBDmax,     // 2015.12.16, I changed the name, Tiecheng
                 s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6, hiso,
                 outputfPentamer, foutputPentamerIso);

 if(watid == 3)
// hpentamersearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBD,
 hpentamersearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBDmax,     // 2015.12.16, I changed the name, Tiecheng
                 s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6, hiso,
                 outputfPentamer, foutputPentamerIso);

}

// Water Cyclic Tetramer search

if(ws == 1 && htetramer == 1)
{
 if(opos == 0) {printf("Water structures: Position of water oxygen ( O ) cannot be read\n"); exit(-1);}
 if(watid == 1) sprintf(foutputTetramer, "%s.%s.Water.Cyclic.Tetramers", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputTetramer, "%s.%s.Water.Cyclic.Tetramers", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputTetramer, "%s.%s.Water.Cyclic.Tetramers", argv[1], argv[4]);
 outputfTetramer = fopen(foutputTetramer,"w");
 if(watid == 1) sprintf(foutputTetramerIso, "%s.%s.Water.Cyclic.Tetramers.isolated", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputTetramerIso, "%s.%s.Water.Cyclic.Tetramers.isolated", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputTetramerIso, "%s.%s.Water.Cyclic.Tetramers.isolated", argv[1], argv[4]);

 if(watid == 1)
// htetramersearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBD,
 htetramersearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBDmax,     // 2015.12.16, I changed the name, Tiecheng
                 s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6, hiso,
                 outputfTetramer, foutputTetramerIso);

 if(watid == 2)
// htetramersearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBD,
 htetramersearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBDmax,     // 2015.12.16, I changed the name, Tiecheng
                 s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6, hiso,
                 outputfTetramer, foutputTetramerIso);

 if(watid == 3)
// htetramersearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBD,
 htetramersearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBDmax,      // 2015.12.16, I changed the name, Tiecheng
                 s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6, hiso,
                 outputfTetramer, foutputTetramerIso);

}

// Water Cyclic Trimer search

if(ws == 1 && htrimer == 1)
{
 if(opos == 0) {printf("Water structures: Position of water oxygen ( O ) cannot be read\n"); exit(-1);}
 if(watid == 1) sprintf(foutputTrimer, "%s.%s.Water.Cyclic.Trimers", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputTrimer, "%s.%s.Water.Cyclic.Trimers", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputTrimer, "%s.%s.Water.Cyclic.Trimers", argv[1], argv[4]);
 outputfTrimer = fopen(foutputTrimer,"w");
 if(watid == 1) sprintf(foutputTrimerIso, "%s.%s.Water.Cyclic.Trimers.isolated", argv[1], argv[2]);
 if(watid == 2) sprintf(foutputTrimerIso, "%s.%s.Water.Cyclic.Trimers.isolated", argv[1], argv[3]);
 if(watid == 3) sprintf(foutputTrimerIso, "%s.%s.Water.Cyclic.Trimers.isolated", argv[1], argv[4]);

 if(watid == 1)
// htrimersearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBD,
 htrimersearch(atmS1, 1, opos, nsolvent1, nAtomS1, s1s1hbdn, s1a, s1b, s1as1bBDmax,     // 2015.12.16, I changed the name, Tiecheng
               s1s1hban, s1s1v1, s1s1v2, s1s1v3, s1s1v4, s1s1v5, s1s1v6, hiso,
               outputfTrimer, foutputTrimerIso);

 if(watid == 2)
// htrimersearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBD,
 htrimersearch(atmS2, 1, opos, nsolvent2, nAtomS2, s2s2hbdn, s2a, s2b, s2as2bBDmax,     // 2015.12.16, I changed the name, Tiecheng
               s2s2hban, s2s2v1, s2s2v2, s2s2v3, s2s2v4, s2s2v5, s2s2v6, hiso,
               outputfTrimer, foutputTrimerIso);

 if(watid == 3)
// htrimersearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBD,
 htrimersearch(atmS3, 1, opos, nsolvent3, nAtomS3, s3s3hbdn, s3a, s3b, s3as3bBDmax,     // 2015.12.16, I changed the name, Tiecheng
               s3s3hban, s3s3v1, s3s3v2, s3s3v3, s3s3v4, s3s3v5, s3s3v6, hiso,
               outputfTrimer, foutputTrimerIso);

}

    
//POLYHEDRA FUNCTIONS
//POLYHEDRA
    if(findpolys==1)
    {
        sprintf(foutputPolys,"%s.Polys",argv[1]);
        outputPolys = fopen(foutputPolys,"a");
        
        node1Start = 1;
        nT1molecules = nAtomT1/nsolute1;
        
        for(i=0; i<nT1molecules; i++)
        {
        
            critlist=findunique(s1t1b, s1t1cutoffnum, &numunique);
                        
            for(j=0; j<numunique; j++)
            {
            
            node2Start=nAtomS1/nsolvent1+i;
            whichatmT1=critlist[j]-1; //this indexing starts at 0, one less than expected in input file. Says "which atom" is being looked at in the molecule. 
                
            whichatmdexT1=i*nsolute1+whichatmT1; //this says the index of the atom in the array atomT1 given the molecule number.
            
            fprintf(outputPolys, "%s, %d, %d, ", argv[2], i, critlist[j]); //input file name, solute molecule number, atom number (whichatmT1+1)

            //start the polyhedra search on this solute atom.
                
//            polyhedra_st(atmS1, atmT1, node1Start, node2Start, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum,  s1t1a, s1t1b, s1t1cutoff, pbc, xside, yside, zside,outputPolys, maxshellsize, edgebds, whichatmT1, whichatmdexT1);
            polyhedra_st(atmS1, atmT1, node1Start, node2Start, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum,  s1t1a, s1t1b, s1t1cutoffmax, pbc, xside, yside, zside,outputPolys, maxshellsize, edgebds, whichatmT1, whichatmdexT1);     // 2015.12.16, I changed the name's1t1cutoffmax' here, Tiecheng
            }
            free(critlist);
            }
    }
    
//PAIRED SHELL DISTANCES
    if(findpdfshell==1)
    {
        sprintf(foutputPdfShell, "%s.pdfs.csv", argv[1]);
        outputPdfShell = fopen(foutputPdfShell, "a");
        
        shellist=(int*)calloc(maxshellsize, sizeof(int));
        shelldistmtx=(double**)malloc(maxshellsize*sizeof(double*));
        
        for(i=0; i<maxshellsize; i++)
        {
            shelldistmtx[i]=(double*)malloc(maxshellsize*sizeof(double));
        }
        
        
        
        node1Start = 1;
        nT1molecules = nAtomT1/nsolute1;

        for(i=0; i<nT1molecules; i++)
        {
            critlist=findunique(s1t1b, s1t1cutoffnum, &numunique);

            
            for(j=0; j<numunique; j++)
            {
                
                node2Start=nAtomS1/nsolvent1+i;
                whichatmT1=critlist[j]-1; //this indexing starts at 0, one less than expected in input file. Says "which atom" is being looked at.
                
                whichatmdexT1=i*nsolute1+whichatmT1; //this says the index of the atom in atomT1 given the molecule number.
                
                shellsize=0;
                free(shellist);
                shellist=(int*)calloc(maxshellsize,sizeof(int));
        
//        shellist_and_distmtx(shellist, shelldistmtx, atmS1, atmT1, node1Start, node2Start, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoff, pbc, xside, yside, zside, &shellsize, maxshellsize, whichatmT1, whichatmdexT1);
        shellist_and_distmtx(shellist, shelldistmtx, atmS1, atmT1, node1Start, node2Start, nsolvent1, nsolute1, nAtomS1, s1t1cutoffnum, s1t1a, s1t1b, s1t1cutoffmax, pbc, xside, yside, zside, &shellsize, maxshellsize, whichatmT1, whichatmdexT1);   // 2015.12.16, I changed the name's1t1cutoffmax' here, Tiecheng

        pdfshell(shelldistmtx, shellsize, outputPdfShell);
                
            }
            free(critlist);

        }

        
        //free values in case.
        /*
        free(shellist);
        for(i=0; i<maxshellsize; i++)
        {
            free(shelldistmtx[i]);
        }
        free(shelldistmtx);
        */

    }
    
    
    //VARY SHELL DISTANCES

    /* Doesn't seem to be working correctly
if(varyshellsz==1)
{
    sprintf(foutputVShellDists, "%s.%s.%s.VarShellDist", argv[1],argv[2],argv[3]);
    outvariedShellDists = fopen(foutputVShellDists, "w+");
    
    
    node1Start = 1;
    nT1molecules = nAtomT1/nsolute1;
    
    for(i=0; i<nT1molecules; i++)
    {
        
        critlist=findunique(s1t1b, s1t1cutoffnum, &numunique);
        
        for(j=0; j<numunique; j++)
        {
            
            node2Start=nAtomS1/nsolvent1+i;
            whichatmT1=critlist[j]-1; //this indexing starts at 0, one less than expected in input file. Says "which atom" is being looked at.
            
            whichatmdexT1=i*nsolute1+whichatmT1; //this says the index of the atom in atomT1 given the molecule number.
            
            
            printf("%d %d \n", whichatmT1, whichatmdexT1);

    
    varyShellFunction(varyMin, varyMax, varyBreaks, maxshellsize, outvariedShellDists, atmS1, atmT1, node1Start, node2Start, nsolvent1, nsolute1, nAtomS1, s1t1a, s1t1b, s1t1cutoff, pbc, xside, yside, zside, whichatmT1, whichatmdexT1, s1t1cutoffnum);
            
        }
        free(critlist);
    }
    
}*/


// Ending............

// time(&now);
// fprintf(outputf,"\nEnd time: %s\n", ctime(&now));

// fclose(outputf);

// closing all output files

if(gs1s1 == 1){ fclose(outputfGraphS1S1); fclose(outputfGeodS1S1); if(pnumnodes==1) fclose(outputfNumnodesS1S1); }
if(gs2s2 == 1){ fclose(outputfGraphS2S2); fclose(outputfGeodS2S2); if(pnumnodes==1) fclose(outputfNumnodesS2S2); }
if(gs3s3 == 1){ fclose(outputfGraphS3S3); fclose(outputfGeodS3S3); if(pnumnodes==1) fclose(outputfNumnodesS3S3); }
if(gs1s2 == 1){ fclose(outputfGraphS1S2); fclose(outputfGeodS1S2); if(pnumnodes==1) fclose(outputfNumnodesS1S2); }
if(gs1s3 == 1){ fclose(outputfGraphS1S3); fclose(outputfGeodS1S3); if(pnumnodes==1) fclose(outputfNumnodesS1S3); }
if(gs2s3 == 1){ fclose(outputfGraphS2S3); fclose(outputfGeodS2S3); if(pnumnodes==1) fclose(outputfNumnodesS2S3); }
if(gs1t1 == 1){ fclose(outputfGraphS1T1); fclose(outputfGeodS1T1); if(pnumnodes==1) fclose(outputfNumnodesS1T1); }
if(gs1t2 == 1){ fclose(outputfGraphS1T2); fclose(outputfGeodS1T2); if(pnumnodes==1) fclose(outputfNumnodesS1T2); }
if(gs2t1 == 1){ fclose(outputfGraphS2T1); fclose(outputfGeodS2T1); if(pnumnodes==1) fclose(outputfNumnodesS2T1); }
if(gs2t2 == 1){ fclose(outputfGraphS2T2); fclose(outputfGeodS2T2); if(pnumnodes==1) fclose(outputfNumnodesS2T2); }
if(gs3t1 == 1){ fclose(outputfGraphS3T1); fclose(outputfGeodS3T1); if(pnumnodes==1) fclose(outputfNumnodesS3T1); }
if(gs3t2 == 1){ fclose(outputfGraphS3T2); fclose(outputfGeodS3T2); if(pnumnodes==1) fclose(outputfNumnodesS3T2); }
if(gs1s2s3 == 1){ fclose(outputfGraphS1S2S3); fclose(outputfGeodS1S2S3); if(pnumnodes==1) fclose(outputfNumnodesS1S2S3); }
if(ws == 1 && hring ==1) fclose(outputfHRing);
if(ws == 1 && hbook ==1) fclose(outputfHBook);
if(ws == 1 && hprism ==1) fclose(outputfHPrism);
if(ws == 1 && hcage ==1) fclose(outputfHCage);
if(ws == 1 && hbag ==1) fclose(outputfHBag);
if(ws == 1 && hboat ==1) fclose(outputfHBoat);
if(ws == 1 && hchair ==1) fclose(outputfHChair);
if(ws == 1 && hprismbook ==1) fclose(outputfHPrismBook);
if(ws == 1 && hpentamer ==1) fclose(outputfPentamer);
if(ws == 1 && htetramer ==1) fclose(outputfTetramer);
if(ws == 1 && htrimer ==1) fclose(outputfTrimer);

    
    if(findpdfshell==1) fclose(outputPdfShell);
    if(findpolys==1) fclose(outputPolys);
 
// free memory

 for(i=0; i < nsolvent1; i++) 
 {
   free(slvntatm1[i]);
 }
 free(slvntatm1);

 for(i=0; i < nAtomS1; i++) 
 {
   free(atmTypeS1[i]);
 }
 free(atmTypeS1);
 
 free(atmS1);

if(nsolventtype == 2 || nsolventtype == 3)
{
 for(i=0; i < nsolvent2; i++)
 {
   free(slvntatm2[i]);
 }
 free(slvntatm2);

 for(i=0; i < nAtomS2; i++) 
 {
   free(atmTypeS2[i]);
 }
 free(atmTypeS2);

 free(atmS2);
}

if(nsolventtype == 3)
{
 for(i=0; i < nsolvent3; i++)
 {
   free(slvntatm3[i]);
 }
 free(slvntatm3);

 for(i=0; i < nAtomS3; i++) 
 {
   free(atmTypeS3[i]);
 }
 free(atmTypeS3);

 free(atmS3);
}

if(nsolutetype != 0)
{
 for(i=0; i < nsolute1; i++)
 {
   free(sltatm1[i]);
 }
 free(sltatm1);

 for(i=0; i < nAtomT1; i++)
 {
   free(atmTypeT1[i]);
 }
 free(atmTypeT1);

 free(atmT1);
}

if(nsolutetype == 2)
{
 for(i=0; i < nsolute2; i++)
 {
   free(sltatm2[i]);
 }
 free(sltatm2);

 for(i=0; i < nAtomT2; i++)
 {
   free(atmTypeT2[i]);
 }
 free(atmTypeT2);

 free(atmT2);
}
    
    if(findpdfshell==1)
    {
        for(i=0; i<maxshellsize; i++)
        {
            free(shelldistmtx[i]);
        }
        free(shelldistmtx);
    }



// Construct the GEODESIC DISTANCE - GD - Matrix

if(gdist == 1 && gdistS1 == 1)
{
 // name of solvent1.xyz file
 sprintf(finput, "%s", argv[2]);

 geodesics_ss(nsolvent1, nAtomS1, EucDistS1, EucRefS1, xside, yside, zside, foutputGeodS1S1, finput);

}

if(gdist == 1 && gdistS2 == 1)
{ 
 // name of solvent2.xyz file
 sprintf(finput, "%s", argv[3]);

 geodesics_ss(nsolvent2, nAtomS2, EucDistS2, EucRefS2, xside, yside, zside, foutputGeodS2S2, finput);

}

if(gdist == 1 && gdistS3 == 1)
{
 // name of solvent3.xyz file
 sprintf(finput, "%s", argv[4]);

 geodesics_ss(nsolvent3, nAtomS3, EucDistS3, EucRefS3, xside, yside, zside, foutputGeodS3S3, finput);

}

// binaries

if(gdist == 1 && gdistS1S2 == 1)
{
 // names of solvent1.xyz and solvent2.xyz files
 sprintf(finput, "%s", argv[2]);
 sprintf(finput2, "%s", argv[3]);

 geodesics_sAsB(nsolvent1, nAtomS1, nsolvent2, nAtomS2, EucDistS1S2, EucRefS1S2s1, EucRefS1S2s2, xside, yside, zside, foutputGeodS1S2, finput, finput2);

}

if(gdist == 1 && gdistS1S3 == 1)
{
 // names of solvent1.xyz and solvent3.xyz files
 sprintf(finput, "%s", argv[2]);
 sprintf(finput3, "%s", argv[4]);

 geodesics_sAsB(nsolvent1, nAtomS1, nsolvent3, nAtomS3, EucDistS1S3, EucRefS1S3s1, EucRefS1S3s3, xside, yside, zside, foutputGeodS1S3, finput, finput3);

}

if(gdist == 1 && gdistS2S3 == 1)
{
 // names of solvent2.xyz and solvent3.xyz files
 sprintf(finput2, "%s", argv[3]);
 sprintf(finput3, "%s", argv[4]);

 geodesics_sAsB(nsolvent2, nAtomS2, nsolvent3, nAtomS3, EucDistS2S3, EucRefS2S3s2, EucRefS2S3s3, xside, yside, zside, foutputGeodS2S3, finput2, finput3);

}

// ternary

if(gdist == 1 && gdistS1S2S3 == 1)
{
 // names of solvent1.xyz solvent2.xyz and solvent3.xyz files
 sprintf(finput, "%s", argv[2]);
 sprintf(finput2, "%s", argv[3]);
 sprintf(finput3, "%s", argv[4]);

 geodesics_sAsBsC(nsolvent1, nAtomS1, nsolvent2, nAtomS2, nsolvent3, nAtomS3, EucDistS1S2S3, EucRefS1S2S3s1, EucRefS1S2S3s2, EucRefS1S2S3s3, xside, yside, zside, 
                  foutputGeodS1S2S3, finput, finput2, finput3);

}

}
