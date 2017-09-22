/*************************************************
 * geodesics.c                                   *
 *                                               *
 * Author: Abdullah Ozkanlar                     *
 *         abdullah.ozkanlar@wsu.edu             *
 *                                               *
 * Based on CK's implementation for water box    *
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

void ChemNetworkNew::geodesics_ss(int nsolvent, int nAtomS, int EucDistS, int EucRefS, double xside, double yside, double zside, char *foutputGeodSS, char *finput)
{

  FILE *fGeodSS; // GraphGeod   input file
  FILE *fxyzS;  //  Solvent.xyz input file
  FILE *outputfPath;
  FILE *outputfCard;
  FILE *fin;

  char foutputPath[256];
  char foutputCard[256];

  char *cptr, buffer[256], line[256];
  int changeval[3];
  int nodea, nodeb, i, j, k, numnodes;
  int numedgelist = 0;
  struct edge *edgelist;
  int **adjmatrix, **gdmatrix, **next;
  double **coord;
  char **aTyp;
  char *str;
  int startnode, endnode;
  int prevnode, curnode, totalchange[3], changeindex;
  double dist;

  numnodes = nAtomS / nsolvent;

  if((fGeodSS = fopen(foutputGeodSS,"r")) == NULL){
      printf("Cannot open file %s\n", foutputGeodSS);
      exit(-1);
  } 

        cptr = fgets(buffer, 256, fGeodSS);
        while(cptr != NULL)
        {
            sscanf(buffer, "%d %d %d %d %d", &nodea, &nodeb, &changeval[0], &changeval[1], &changeval[2]);
            
            if (numedgelist == 0)
                edgelist = (struct edge *) malloc(sizeof(struct edge));
            else
                edgelist = (struct edge *) realloc(edgelist, sizeof(struct edge) * (numedgelist + 1));
            edgelist[numedgelist].nodea = nodea-1;
            edgelist[numedgelist].nodeb = nodeb-1;
            edgelist[numedgelist].change[0] = changeval[0];
            edgelist[numedgelist].change[1] = changeval[1];
            edgelist[numedgelist].change[2] = changeval[2];
            numedgelist++;
            cptr = fgets(buffer, 256, fGeodSS);
            
        }
        fclose(fGeodSS);

        // compute the geodesic path
        // build the adjaceny matrix
        adjmatrix = (int **)malloc(sizeof(int *) * numnodes);
        for (i=0; i < numnodes; i++)
            adjmatrix[i] = (int *)calloc(numnodes, sizeof(int));
        
        for (i=0; i < numedgelist; i++)
        {
            adjmatrix[edgelist[i].nodea][edgelist[i].nodeb] = 1;
            adjmatrix[edgelist[i].nodeb][edgelist[i].nodea] = 1;
        }
        
        // compute the geodesic distance - gd - between the nodes
        gdmatrix = (int **)malloc(sizeof(float *) * numnodes);
        for (i=0; i < numnodes; i++)
            gdmatrix[i] = (int *)calloc(numnodes, sizeof(int));
        
        next = (int **)malloc(sizeof(int *) * numnodes);
        for (i=0; i<numnodes; i++)
            next[i] = (int *)calloc(numnodes, sizeof(int));

        // initialize the geodesic distance matrix with values in the adjacency matrix
        // gd(i, i) = 0;
        // gd(i, j) = adjmatrix(i, j) if there is an edge between i and j
        // gd(i, j) = LARGE if there is no edge between i and j
        for (i=0; i < numnodes; i++)
            for (j=0; j < numnodes; j++)
            {
                if (adjmatrix[i][j] != 0)
                    gdmatrix[i][j] = adjmatrix[i][j];
                else if (i!=j)
                    gdmatrix[i][j] = LARGE;
                
                next[i][j] = LARGE;
            }
        
        // floyd warshall algorithm for computing the shortest path between the nodes
        for (k=0; k < numnodes; k++)
            for (i=0; i < numnodes; i++)
                for (j=0; j < numnodes; j++)
                    if (gdmatrix[i][j] > (gdmatrix[i][k] + gdmatrix[k][j]))
                    {
                        gdmatrix[i][j] = gdmatrix[i][k] + gdmatrix[k][j];
                        next[i][j] = k;
                    }

        // print the geodesic distance between the nodes that are connected to each other.
        sprintf(foutputPath, "%s.geopath", finput);
        outputfPath = fopen(foutputPath, "w");
        for (i=0; i < numnodes; i++)
            for( j=i+1; j < numnodes; j++)
                if (gdmatrix[i][j] < LARGE)
                {
                    if (gdmatrix[i][j] > 1)
                    {
                        fprintf(outputfPath, "%d %d %d path: ", i+1, j+1, gdmatrix[i][j]);
                        GetPath(i, j, next, gdmatrix, outputfPath);
                        fprintf(outputfPath, "%d", j+1);
                    }
                    else
                        fprintf(outputfPath, "%d %d %d path: %d %d", i+1, j+1, gdmatrix[i][j], i+1, j+1);
                    fprintf(outputfPath, "\n");
                }
        fclose(outputfPath);



// If Euclidean distance is requested

if(EucDistS == 1)
{

  if((fxyzS = fopen(finput,"r")) == NULL){
      printf("Cannot open file %s\n", finput);
      exit(-1);
  } 

   coord = (double**) malloc(numnodes*sizeof(double*));
   for (i = 0; i < numnodes; i++)
      coord[i] = (double*) malloc(3*sizeof(double));

   aTyp = (char **)malloc(sizeof(char *)*numnodes);
     for(i = 0; i < numnodes; i++) aTyp[i] = (char *)malloc(sizeof(char)*32);


      // Skip 1 line from top of .xyz file
      rewind(fxyzS);
      fgets(line,256,fxyzS);
      // Go to line of first reference node atom - EucRef
      for(i = 1; i <= EucRefS; i++)
         fgets(line,256,fxyzS);

     // Read coordinates of reference node atom for euclidean dist. calc. - EucRef

    for(i = 0; i < numnodes; i++)
    {

       fscanf(fxyzS, "%31s %lf %lf %lf", aTyp[i], &coord[i][0], &coord[i][1], &coord[i][2]); 
       
        // go to next ref node 
        for(k = 1; k <= nsolvent; k++)
           fgets(line,256,fxyzS);

    } 

    fclose(fxyzS);


   fin = fopen(foutputPath, "r");
   if (fin == NULL)
   {
       printf("error reading the path file: %s\n", foutputPath);
       exit(-1);
   }

   sprintf(foutputCard, "%s.geocard", finput);
   outputfCard = fopen(foutputCard, "w");

        // start reading the paths one by one and 
        // compute the distance in the cartesian coordinates
        // the first line will have snapshot in it.
        cptr = fgets(buffer, 256, fin);
        // start node
        str = strtok(cptr, " \n");
        fprintf(outputfCard, "%s ", str);
        // end node
        str = strtok(NULL, " \n");
        fprintf(outputfCard, "%s ", str);
        // length of the paths
        str = strtok(NULL, " \n");
        fprintf(outputfCard, "%s ", str);
        // read path:
        str = strtok(NULL, " \n");
        fprintf(outputfCard, "%s ", str);
        // obtain the start node
        str = strtok(NULL, " \n");
        prevnode = atoi(str)-1;
        startnode = prevnode;
        fprintf(outputfCard, "%d ", startnode+1);

        while (cptr != NULL)
        {
            str = strtok(NULL, " \n");
            totalchange[0] = 0;
            totalchange[1] = 0;
            totalchange[2] = 0;
            while (str != NULL)
            {
                // determine the successor node
                curnode = atoi(str)-1;
                // find the change
                changeindex = FindChange(prevnode, curnode, edgelist, numedgelist);
                // add it to the total change
                if (changeindex < 0)
                {
                    totalchange[0] += edgelist[-1*changeindex].change[0]*-1;
                    totalchange[1] += edgelist[-1*changeindex].change[1]*-1;
                    totalchange[2] += edgelist[-1*changeindex].change[2]*-1;
                }
                else
                {
                    totalchange[0] += edgelist[changeindex].change[0];
                    totalchange[1] += edgelist[changeindex].change[1];
                    totalchange[2] += edgelist[changeindex].change[2];
                }
                // assign the current node as the start node and
                // move to the next node
                prevnode = curnode;
                fprintf(outputfCard, "%d ", curnode+1);
                str = strtok(NULL, " \n");
            }
            endnode = prevnode;
            dist = 0;
            
            dist += (double)pow((coord[startnode][0] - (coord[endnode][0]+xside*totalchange[0])), 2);
            
            dist += (double)pow((coord[startnode][1] - (coord[endnode][1]+yside*totalchange[1])), 2);
            
            dist += (double)pow(((coord[startnode][2] - (coord[endnode][2]+zside*totalchange[2]))), 2);
            
            fprintf(outputfCard, "distance: %lf\n", pow(dist, 0.5));
            cptr = fgets(buffer, 256, fin);
            if (cptr != NULL)
            {
                // start node
                str = strtok(cptr, " \n");
                fprintf(outputfCard, "%s ", str);
                // end node
                str = strtok(NULL, " \n");
                fprintf(outputfCard, "%s ", str);
                // length of the paths
                str = strtok(NULL, " \n");
                fprintf(outputfCard, "%s ", str);
                // read path:
                str = strtok(NULL, " \n");
                fprintf(outputfCard, "%s ", str);
                // obtain the start node
                str = strtok(NULL, " \n");
                prevnode = atoi(str)-1;
                startnode = prevnode;
                fprintf(outputfCard, "%d ", startnode+1);
            }

        }
        
        fclose(fin);
        fclose(outputfCard);

}



        // free memory
        for (i=0; i < numnodes; i++)
        {
            free(coord[i]);    
            free(gdmatrix[i]);
            free(adjmatrix[i]);
            free(next[i]);
        }
        free(coord);           
        free(edgelist);
        free(gdmatrix);
        free(adjmatrix);
        free(next);

}

// Geodesics for binaries

void ChemNetworkNew::geodesics_sAsB(int nsolventA, int nAtomSA, int nsolventB, int nAtomSB, int EucDistSASB, int EucRefSASBsa, int EucRefSASBsb, double xside, double yside, double zside, 
                    char *foutputGeodSASB, char *finputA, char *finputB)
{

  FILE *fGeodSS; // GraphGeod   input file
  FILE *fxyzSA;  //  SolventA.xyz input file
  FILE *fxyzSB;  //  SolventB.xyz input file
  FILE *outputfPath;
  FILE *outputfCard;
  FILE *fin;

  char foutputPath[256];
  char foutputCard[256];

  char *cptr, buffer[256], line[256];
  int changeval[3];
  int nodea, nodeb, i, j, k, numnodes, numnodesA, numnodesB;
  int numedgelist = 0;
  struct edge *edgelist;
  int **adjmatrix, **gdmatrix, **next;
  double **coord;
  char **aTyp;
  char *str;
  int startnode, endnode;
  int prevnode, curnode, totalchange[3], changeindex;
  double dist;

  numnodesA = nAtomSA / nsolventA;
  numnodesB = nAtomSB / nsolventB;

  numnodes = numnodesA + numnodesB;

  if((fGeodSS = fopen(foutputGeodSASB,"r")) == NULL){
      printf("Cannot open file %s\n", foutputGeodSASB);
      exit(-1);
  } 

        cptr = fgets(buffer, 256, fGeodSS);
        while(cptr != NULL)
        {
            sscanf(buffer, "%d %d %d %d %d", &nodea, &nodeb, &changeval[0], &changeval[1], &changeval[2]);
            
            if (numedgelist == 0)
                edgelist = (struct edge *) malloc(sizeof(struct edge));
            else
                edgelist = (struct edge *) realloc(edgelist, sizeof(struct edge) * (numedgelist + 1));
            edgelist[numedgelist].nodea = nodea-1;
            edgelist[numedgelist].nodeb = nodeb-1;
            edgelist[numedgelist].change[0] = changeval[0];
            edgelist[numedgelist].change[1] = changeval[1];
            edgelist[numedgelist].change[2] = changeval[2];
            numedgelist++;
            cptr = fgets(buffer, 256, fGeodSS);
            
        }
        fclose(fGeodSS);

        // compute the geodesic path
        // build the adjaceny matrix
        adjmatrix = (int **)malloc(sizeof(int *) * numnodes);
        for (i=0; i < numnodes; i++)
            adjmatrix[i] = (int *)calloc(numnodes, sizeof(int));
        
        for (i=0; i < numedgelist; i++)
        {
            adjmatrix[edgelist[i].nodea][edgelist[i].nodeb] = 1;
            adjmatrix[edgelist[i].nodeb][edgelist[i].nodea] = 1;
        }
        
        // compute the geodesic distance - gd - between the nodes
        gdmatrix = (int **)malloc(sizeof(float *) * numnodes);
        for (i=0; i < numnodes; i++)
            gdmatrix[i] = (int *)calloc(numnodes, sizeof(int));
        
        next = (int **)malloc(sizeof(int *) * numnodes);
        for (i=0; i<numnodes; i++)
            next[i] = (int *)calloc(numnodes, sizeof(int));

        // initialize the geodesic distance matrix with values in the adjacency matrix
        // gd(i, i) = 0;
        // gd(i, j) = adjmatrix(i, j) if there is an edge between i and j
        // gd(i, j) = LARGE if there is no edge between i and j
        for (i=0; i < numnodes; i++)
            for (j=0; j < numnodes; j++)
            {
                if (adjmatrix[i][j] != 0)
                    gdmatrix[i][j] = adjmatrix[i][j];
                else if (i!=j)
                    gdmatrix[i][j] = LARGE;
                
                next[i][j] = LARGE;
            }
        
        // floyd warshall algorithm for computing the shortest path between the nodes
        for (k=0; k < numnodes; k++)
            for (i=0; i < numnodes; i++)
                for (j=0; j < numnodes; j++)
                    if (gdmatrix[i][j] > (gdmatrix[i][k] + gdmatrix[k][j]))
                    {
                        gdmatrix[i][j] = gdmatrix[i][k] + gdmatrix[k][j];
                        next[i][j] = k;
                    }

        // print the geodesic distance between the nodes that are connected to each other.
        sprintf(foutputPath, "%s.%s.geopath", finputA, finputB);
        outputfPath = fopen(foutputPath, "w");
        for (i=0; i < numnodes; i++)
            for( j=i+1; j < numnodes; j++)
                if (gdmatrix[i][j] < LARGE)
                {
                    if (gdmatrix[i][j] > 1)
                    {
                        fprintf(outputfPath, "%d %d %d path: ", i+1, j+1, gdmatrix[i][j]);
                        GetPath(i, j, next, gdmatrix, outputfPath);
                        fprintf(outputfPath, "%d", j+1);
                    }
                    else
                        fprintf(outputfPath, "%d %d %d path: %d %d", i+1, j+1, gdmatrix[i][j], i+1, j+1);
                    fprintf(outputfPath, "\n");
                }
        fclose(outputfPath);



// If Euclidean distance is requested

if(EucDistSASB == 1)
{

   coord = (double**) malloc(numnodes*sizeof(double*));
   for (i = 0; i < numnodes; i++)
      coord[i] = (double*) malloc(3*sizeof(double));

   aTyp = (char **)malloc(sizeof(char *)*numnodes);
     for(i = 0; i < numnodes; i++) aTyp[i] = (char *)malloc(sizeof(char)*32);

  // Read coordinates of solvent A
  if((fxyzSA = fopen(finputA,"r")) == NULL){
      printf("Cannot open file %s\n", finputA);
      exit(-1);
  }

      // Skip 1 line from top of .xyz file
      rewind(fxyzSA);
      fgets(line,256,fxyzSA);
      // Go to line of first reference node atom - EucRef
      for(i = 1; i <= EucRefSASBsa; i++)
         fgets(line,256,fxyzSA);

     // Read coordinates of reference node atom for euclidean dist. calc. - EucRef

    for(i = 0; i < numnodesA; i++)
    {

       fscanf(fxyzSA, "%31s %lf %lf %lf", aTyp[i], &coord[i][0], &coord[i][1], &coord[i][2]); 
       
        // go to next ref node 
        for(k = 1; k <= nsolventA; k++)
           fgets(line,256,fxyzSA);

    } 

    fclose(fxyzSA);

  // Read coordinates of solvent B
  if((fxyzSB = fopen(finputB,"r")) == NULL){
      printf("Cannot open file %s\n", finputB);
      exit(-1);
  }

      // Skip 1 line from top of .xyz file
      rewind(fxyzSB);
      fgets(line,256,fxyzSB);
      // Go to line of first reference node atom - EucRef
      for(i = 1; i <= EucRefSASBsb; i++)
         fgets(line,256,fxyzSB);

     // Read coordinates of reference node atom for euclidean dist. calc. - EucRef

    for(i = numnodesA; i < numnodes; i++)
    {

       fscanf(fxyzSB, "%31s %lf %lf %lf", aTyp[i], &coord[i][0], &coord[i][1], &coord[i][2]);

        // go to next ref node 
        for(k = 1; k <= nsolventB; k++)
           fgets(line,256,fxyzSB);

    }

    fclose(fxyzSB);




   fin = fopen(foutputPath, "r");
   if (fin == NULL)
   {
       printf("error reading the path file: %s\n", foutputPath);
       exit(-1);
   }

   sprintf(foutputCard, "%s.%s.geocard", finputA, finputB);
   outputfCard = fopen(foutputCard, "w");

        // start reading the paths one by one and 
        // compute the distance in the cartesian coordinates
        // the first line will have snapshot in it.
        cptr = fgets(buffer, 256, fin);
        // start node
        str = strtok(cptr, " \n");
        fprintf(outputfCard, "%s ", str);
        // end node
        str = strtok(NULL, " \n");
        fprintf(outputfCard, "%s ", str);
        // length of the paths
        str = strtok(NULL, " \n");
        fprintf(outputfCard, "%s ", str);
        // read path:
        str = strtok(NULL, " \n");
        fprintf(outputfCard, "%s ", str);
        // obtain the start node
        str = strtok(NULL, " \n");
        prevnode = atoi(str)-1;
        startnode = prevnode;
        fprintf(outputfCard, "%d ", startnode+1);

        while (cptr != NULL)
        {
            str = strtok(NULL, " \n");
            totalchange[0] = 0;
            totalchange[1] = 0;
            totalchange[2] = 0;
            while (str != NULL)
            {
                // determine the successor node
                curnode = atoi(str)-1;
                // find the change
                changeindex = FindChange(prevnode, curnode, edgelist, numedgelist);
                // add it to the total change
                if (changeindex < 0)
                {
                    totalchange[0] += edgelist[-1*changeindex].change[0]*-1;
                    totalchange[1] += edgelist[-1*changeindex].change[1]*-1;
                    totalchange[2] += edgelist[-1*changeindex].change[2]*-1;
                }
                else
                {
                    totalchange[0] += edgelist[changeindex].change[0];
                    totalchange[1] += edgelist[changeindex].change[1];
                    totalchange[2] += edgelist[changeindex].change[2];
                }
                // assign the current node as the start node and
                // move to the next node
                prevnode = curnode;
                fprintf(outputfCard, "%d ", curnode+1);
                str = strtok(NULL, " \n");
            }
            endnode = prevnode;
            dist = 0;
            
            dist += (double)pow((coord[startnode][0] - (coord[endnode][0]+xside*totalchange[0])), 2);
            
            dist += (double)pow((coord[startnode][1] - (coord[endnode][1]+yside*totalchange[1])), 2);
            
            dist += (double)pow(((coord[startnode][2] - (coord[endnode][2]+zside*totalchange[2]))), 2);
            
            fprintf(outputfCard, "distance: %lf\n", pow(dist, 0.5));
            cptr = fgets(buffer, 256, fin);
            if (cptr != NULL)
            {
                // start node
                str = strtok(cptr, " \n");
                fprintf(outputfCard, "%s ", str);
                // end node
                str = strtok(NULL, " \n");
                fprintf(outputfCard, "%s ", str);
                // length of the paths
                str = strtok(NULL, " \n");
                fprintf(outputfCard, "%s ", str);
                // read path:
                str = strtok(NULL, " \n");
                fprintf(outputfCard, "%s ", str);
                // obtain the start node
                str = strtok(NULL, " \n");
                prevnode = atoi(str)-1;
                startnode = prevnode;
                fprintf(outputfCard, "%d ", startnode+1);
            }

        }
        
        fclose(fin);
        fclose(outputfCard);

}



        // free memory
        for (i=0; i < numnodes; i++)
        {
            free(coord[i]);    
            free(gdmatrix[i]);
            free(adjmatrix[i]);
            free(next[i]);
        }
        free(coord);           
        free(edgelist);
        free(gdmatrix);
        free(adjmatrix);
        free(next);

}

// ternary

void ChemNetworkNew::geodesics_sAsBsC(int nsolventA, int nAtomSA, int nsolventB, int nAtomSB, int nsolventC, int nAtomSC, int EucDistSASBSC, int EucRefSASBSCsa, int EucRefSASBSCsb, int EucRefSASBSCsc, 
                      double xside, double yside, double zside, char *foutputGeodSASBSC, char *finputA, char *finputB, char *finputC)
{

  FILE *fGeodSS; // GraphGeod   input file
  FILE *fxyzSA;  //  SolventA.xyz input file
  FILE *fxyzSB;  //  SolventB.xyz input file
  FILE *fxyzSC;  //  SolventC.xyz input file
  FILE *outputfPath;
  FILE *outputfCard;
  FILE *fin;

  char foutputPath[256];
  char foutputCard[256];

  char *cptr, buffer[256], line[256];
  int changeval[3];
  int nodea, nodeb, i, j, k, numnodes, numnodesA, numnodesB, numnodesC;
  int numedgelist = 0;
  struct edge *edgelist;
  int **adjmatrix, **gdmatrix, **next;
  double **coord;
  char **aTyp;
  char *str;
  int startnode, endnode;
  int prevnode, curnode, totalchange[3], changeindex;
  double dist;

  numnodesA = nAtomSA / nsolventA;
  numnodesB = nAtomSB / nsolventB;
  numnodesC = nAtomSC / nsolventC;

  numnodes = numnodesA + numnodesB + numnodesC;

  if((fGeodSS = fopen(foutputGeodSASBSC,"r")) == NULL){
      printf("Cannot open file %s\n", foutputGeodSASBSC);
      exit(-1);
  } 

        cptr = fgets(buffer, 256, fGeodSS);
        while(cptr != NULL)
        {
            sscanf(buffer, "%d %d %d %d %d", &nodea, &nodeb, &changeval[0], &changeval[1], &changeval[2]);
            
            if (numedgelist == 0)
                edgelist = (struct edge *) malloc(sizeof(struct edge));
            else
                edgelist = (struct edge *) realloc(edgelist, sizeof(struct edge) * (numedgelist + 1));
            edgelist[numedgelist].nodea = nodea-1;
            edgelist[numedgelist].nodeb = nodeb-1;
            edgelist[numedgelist].change[0] = changeval[0];
            edgelist[numedgelist].change[1] = changeval[1];
            edgelist[numedgelist].change[2] = changeval[2];
            numedgelist++;
            cptr = fgets(buffer, 256, fGeodSS);
            
        }
        fclose(fGeodSS);

        // compute the geodesic path
        // build the adjaceny matrix
        adjmatrix = (int **)malloc(sizeof(int *) * numnodes);
        for (i=0; i < numnodes; i++)
            adjmatrix[i] = (int *)calloc(numnodes, sizeof(int));
        
        for (i=0; i < numedgelist; i++)
        {
            adjmatrix[edgelist[i].nodea][edgelist[i].nodeb] = 1;
            adjmatrix[edgelist[i].nodeb][edgelist[i].nodea] = 1;
        }
        
        // compute the geodesic distance - gd - between the nodes
        gdmatrix = (int **)malloc(sizeof(float *) * numnodes);
        for (i=0; i < numnodes; i++)
            gdmatrix[i] = (int *)calloc(numnodes, sizeof(int));
        
        next = (int **)malloc(sizeof(int *) * numnodes);
        for (i=0; i<numnodes; i++)
            next[i] = (int *)calloc(numnodes, sizeof(int));

        // initialize the geodesic distance matrix with values in the adjacency matrix
        // gd(i, i) = 0;
        // gd(i, j) = adjmatrix(i, j) if there is an edge between i and j
        // gd(i, j) = LARGE if there is no edge between i and j
        for (i=0; i < numnodes; i++)
            for (j=0; j < numnodes; j++)
            {
                if (adjmatrix[i][j] != 0)
                    gdmatrix[i][j] = adjmatrix[i][j];
                else if (i!=j)
                    gdmatrix[i][j] = LARGE;
                
                next[i][j] = LARGE;
            }
        
        // floyd warshall algorithm for computing the shortest path between the nodes
        for (k=0; k < numnodes; k++)
            for (i=0; i < numnodes; i++)
                for (j=0; j < numnodes; j++)
                    if (gdmatrix[i][j] > (gdmatrix[i][k] + gdmatrix[k][j]))
                    {
                        gdmatrix[i][j] = gdmatrix[i][k] + gdmatrix[k][j];
                        next[i][j] = k;
                    }

        // print the geodesic distance between the nodes that are connected to each other.
        sprintf(foutputPath, "%s.%s.%s.geopath", finputA, finputB, finputC);
        outputfPath = fopen(foutputPath, "w");
        for (i=0; i < numnodes; i++)
            for( j=i+1; j < numnodes; j++)
                if (gdmatrix[i][j] < LARGE)
                {
                    if (gdmatrix[i][j] > 1)
                    {
                        fprintf(outputfPath, "%d %d %d path: ", i+1, j+1, gdmatrix[i][j]);
                        GetPath(i, j, next, gdmatrix, outputfPath);
                        fprintf(outputfPath, "%d", j+1);
                    }
                    else
                        fprintf(outputfPath, "%d %d %d path: %d %d", i+1, j+1, gdmatrix[i][j], i+1, j+1);
                    fprintf(outputfPath, "\n");
                }
        fclose(outputfPath);



// If Euclidean distance is requested

if(EucDistSASBSC == 1)
{

   coord = (double**) malloc(numnodes*sizeof(double*));
   for (i = 0; i < numnodes; i++)
      coord[i] = (double*) malloc(3*sizeof(double));

   aTyp = (char **)malloc(sizeof(char *)*numnodes);
     for(i = 0; i < numnodes; i++) aTyp[i] = (char *)malloc(sizeof(char)*32);

  // Read coordinates of solvent A
  if((fxyzSA = fopen(finputA,"r")) == NULL){
      printf("Cannot open file %s\n", finputA);
      exit(-1);
  }

      // Skip 1 line from top of .xyz file
      rewind(fxyzSA);
      fgets(line,256,fxyzSA);
      // Go to line of first reference node atom - EucRef
      for(i = 1; i <= EucRefSASBSCsa; i++)
         fgets(line,256,fxyzSA);

     // Read coordinates of reference node atom for euclidean dist. calc. - EucRef

    for(i = 0; i < numnodesA; i++)
    {

       fscanf(fxyzSA, "%31s %lf %lf %lf", aTyp[i], &coord[i][0], &coord[i][1], &coord[i][2]); 
       
        // go to next ref node 
        for(k = 1; k <= nsolventA; k++)
           fgets(line,256,fxyzSA);

    } 

    fclose(fxyzSA);

  // Read coordinates of solvent B
  if((fxyzSB = fopen(finputB,"r")) == NULL){
      printf("Cannot open file %s\n", finputB);
      exit(-1);
  }

      // Skip 1 line from top of .xyz file
      rewind(fxyzSB);
      fgets(line,256,fxyzSB);
      // Go to line of first reference node atom - EucRef
      for(i = 1; i <= EucRefSASBSCsb; i++)
         fgets(line,256,fxyzSB);

     // Read coordinates of reference node atom for euclidean dist. calc. - EucRef

    for(i = numnodesA; i < (numnodesA + numnodesB); i++)
    {

       fscanf(fxyzSB, "%31s %lf %lf %lf", aTyp[i], &coord[i][0], &coord[i][1], &coord[i][2]);

        // go to next ref node 
        for(k = 1; k <= nsolventB; k++)
           fgets(line,256,fxyzSB);

    }

    fclose(fxyzSB);

  // Read coordinates of solvent C
  if((fxyzSC = fopen(finputC,"r")) == NULL){
      printf("Cannot open file %s\n", finputC);
      exit(-1);
  }

      // Skip 1 line from top of .xyz file
      rewind(fxyzSC);
      fgets(line,256,fxyzSC);
      // Go to line of first reference node atom - EucRef
      for(i = 1; i <= EucRefSASBSCsc; i++)
         fgets(line,256,fxyzSC);

     // Read coordinates of reference node atom for euclidean dist. calc. - EucRef

    for(i = (numnodesA + numnodesB); i < numnodes; i++)
    {

       fscanf(fxyzSC, "%31s %lf %lf %lf", aTyp[i], &coord[i][0], &coord[i][1], &coord[i][2]);

        // go to next ref node 
        for(k = 1; k <= nsolventC; k++)
           fgets(line,256,fxyzSC);

    }

    fclose(fxyzSC);



   fin = fopen(foutputPath, "r");
   if (fin == NULL)
   {
       printf("error reading the path file: %s\n", foutputPath);
       exit(-1);
   }

   sprintf(foutputCard, "%s.%s.%s.geocard", finputA, finputB, finputC);
   outputfCard = fopen(foutputCard, "w");

        // start reading the paths one by one and 
        // compute the distance in the cartesian coordinates
        // the first line will have snapshot in it.
        cptr = fgets(buffer, 256, fin);
        // start node
        str = strtok(cptr, " \n");
        fprintf(outputfCard, "%s ", str);
        // end node
        str = strtok(NULL, " \n");
        fprintf(outputfCard, "%s ", str);
        // length of the paths
        str = strtok(NULL, " \n");
        fprintf(outputfCard, "%s ", str);
        // read path:
        str = strtok(NULL, " \n");
        fprintf(outputfCard, "%s ", str);
        // obtain the start node
        str = strtok(NULL, " \n");
        prevnode = atoi(str)-1;
        startnode = prevnode;
        fprintf(outputfCard, "%d ", startnode+1);

        while (cptr != NULL)
        {
            str = strtok(NULL, " \n");
            totalchange[0] = 0;
            totalchange[1] = 0;
            totalchange[2] = 0;
            while (str != NULL)
            {
                // determine the successor node
                curnode = atoi(str)-1;
                // find the change
                changeindex = FindChange(prevnode, curnode, edgelist, numedgelist);
                // add it to the total change
                if (changeindex < 0)
                {
                    totalchange[0] += edgelist[-1*changeindex].change[0]*-1;
                    totalchange[1] += edgelist[-1*changeindex].change[1]*-1;
                    totalchange[2] += edgelist[-1*changeindex].change[2]*-1;
                }
                else
                {
                    totalchange[0] += edgelist[changeindex].change[0];
                    totalchange[1] += edgelist[changeindex].change[1];
                    totalchange[2] += edgelist[changeindex].change[2];
                }
                // assign the current node as the start node and
                // move to the next node
                prevnode = curnode;
                fprintf(outputfCard, "%d ", curnode+1);
                str = strtok(NULL, " \n");
            }
            endnode = prevnode;
            dist = 0;
            
            dist += (double)pow((coord[startnode][0] - (coord[endnode][0]+xside*totalchange[0])), 2);
            
            dist += (double)pow((coord[startnode][1] - (coord[endnode][1]+yside*totalchange[1])), 2);
            
            dist += (double)pow(((coord[startnode][2] - (coord[endnode][2]+zside*totalchange[2]))), 2);
            
            fprintf(outputfCard, "distance: %lf\n", pow(dist, 0.5));
            cptr = fgets(buffer, 256, fin);
            if (cptr != NULL)
            {
                // start node
                str = strtok(cptr, " \n");
                fprintf(outputfCard, "%s ", str);
                // end node
                str = strtok(NULL, " \n");
                fprintf(outputfCard, "%s ", str);
                // length of the paths
                str = strtok(NULL, " \n");
                fprintf(outputfCard, "%s ", str);
                // read path:
                str = strtok(NULL, " \n");
                fprintf(outputfCard, "%s ", str);
                // obtain the start node
                str = strtok(NULL, " \n");
                prevnode = atoi(str)-1;
                startnode = prevnode;
                fprintf(outputfCard, "%d ", startnode+1);
            }

        }
        
        fclose(fin);
        fclose(outputfCard);

}



        // free memory
        for (i=0; i < numnodes; i++)
        {
            free(coord[i]);    
            free(gdmatrix[i]);
            free(adjmatrix[i]);
            free(next[i]);
        }
        free(coord);           
        free(edgelist);
        free(gdmatrix);
        free(adjmatrix);
        free(next);

}



void ChemNetworkNew::GetPath(int i, int j, int **next, int **gdmatrix, FILE *fout)
{        
    if (gdmatrix[i][j] == 1)
    {
        fprintf(fout, "%d ", i+1);
        return;
    }
    else
    {
        GetPath(i, next[i][j], next, gdmatrix, fout);
        GetPath(next[i][j], j, next, gdmatrix, fout);
        return;
    } 
}

int ChemNetworkNew::FindChange(int startnode, int endnode, struct edge* edgelist, int numedgelist)
{
    
    int i;
    
    for (i=0; i<numedgelist; i++)
    {
        if (edgelist[i].nodea == startnode)
            if (edgelist[i].nodeb == endnode)
                return(i);
        if (edgelist[i].nodeb == startnode)
            if (edgelist[i].nodea == endnode)
                return(-i);
    }
    
}





