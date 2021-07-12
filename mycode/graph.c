/*
graph.c

Set of vertices and edges implementation.

Implementations for helper functions for graph construction and manipulation.
 
Skeleton written by Grady Fitzpatrick for COMP20007 Assignment 1 2021
*/
#include<stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include "graph.h"
#include "utils.h"
#include "pq.h"

#define INITIALEDGES 32


typedef struct item{
  int cost;
  int prev;
}items;

/* Definition of a graph. */
struct graph {
  int numVertices;
  int numEdges;
  int allocedEdges;
  struct edge **edgeList;
};

/* Definition of an edge. */
struct edge {
  int start;
  int end;
  int cost;
};
int prims(struct graph* g, int numHouse, items **lst, int table[(g->numVertices)][(g->numVertices)], 
int* marked, enum problemPart part, int antennaCost);
struct graph *newGraph(int numVertices){
  struct graph *g = (struct graph *) malloc(sizeof(struct graph));
  assert(g);
  /* Initialise edges. */
  g->numVertices = numVertices;
  g->numEdges = 0;
  g->allocedEdges = 0;
  g->edgeList = NULL;
  return g;
}

/* Adds an edge to the given graph. */
void addEdge(struct graph *g, int start, int end, int cost){
  assert(g);
  struct edge *newEdge = NULL;
  /* Check we have enough space for the new edge. */
  if((g->numEdges + 1) > g->allocedEdges){
    if(g->allocedEdges == 0){
      g->allocedEdges = INITIALEDGES;
    } else {
      (g->allocedEdges) *= 2;
    }
    g->edgeList = (struct edge **) realloc(g->edgeList,
      sizeof(struct edge *) * g->allocedEdges);
    assert(g->edgeList);
  }

  /* Create the edge */
  newEdge = (struct edge *) malloc(sizeof(struct edge));
  assert(newEdge);
  newEdge->start = start;
  newEdge->end = end;
  newEdge->cost = cost;

  /* Add the edge to the list of edges. */
  g->edgeList[g->numEdges] = newEdge;
  (g->numEdges)++;
}

/* Frees all memory used by graph. */
void freeGraph(struct graph *g){
  int i;
  for(i = 0; i < g->numEdges; i++){
    free((g->edgeList)[i]);
  }
  if(g->edgeList){
    free(g->edgeList);
  }
  free(g);
}

struct solution *graphSolve(struct graph *g, enum problemPart part,
  int antennaCost, int numHouses){

  struct solution *solution = (struct solution *)
    malloc(sizeof(struct solution));
  assert(solution);
  struct edge** edgeList;
    edgeList = g->edgeList;
    int *marked = malloc((g->numVertices) * sizeof(int));
    items *lst [g->numVertices];
    
    
     
    
    lst[0] = malloc(sizeof(items));
    lst[0]-> prev = 0;
    lst[0] ->cost = 0;
   
    marked[0] = 0;
    
    
   
    int table [g->numVertices][g->numVertices];
    for(int i = 0; i < g->numVertices; i++){
      for(int j = 0; j < g->numVertices; j++){
        table[i][j] = 0;

      }
     
    }
   
    for(int i = 1; i < g->numVertices; i++){
      
      lst[i] = malloc(sizeof(items));
      lst[i]->prev = i;
     
      lst[i]->cost =INT_MAX;
      
      marked[i] = 0;
    }
  
    for(int i = 0; i < g->numEdges; i++){

      table[edgeList[i]->start][edgeList[i]->start] = 0;
      table[edgeList[i]->start][edgeList[i]->end] = edgeList[i]->cost;
      table[edgeList[i]->end][edgeList[i]->start] = edgeList[i]->cost;
      
    

      
    }
  if(part == PART_A){
    /* IMPLEMENT 2A SOLUTION HERE */
    solution->antennaTotal = antennaCost * numHouses;
   
    solution->cableTotal = prims(g, numHouses, lst, table,  marked, part, antennaCost);
  } else {
    /* IMPLEMENT 2C SOLUTION HERE */
    solution->mixedTotal =prims(g, numHouses, lst, table,  marked, part, antennaCost);
  }
  return solution;
}

int prims(struct graph* g, int numHouse, items **lst, int table[(g->numVertices)][(g->numVertices)], 
int* marked, enum problemPart part, int antennaCost){
    
    items* minNode;
   
    struct pq* pq = newPQ();
  
    enqueue(pq, lst[0], lst[0]->cost);

    while(!empty(pq)){
      
      minNode = deletemin(pq);
      
      marked[minNode->prev] = 1;
     
      for(int i = 0; i < (g ->numVertices); i++){
       
        if(table[minNode->prev][i] == 0)  continue;
       
        if ( !marked[i] && lst[i]->cost > table[minNode->prev][i]){
          lst[i]->cost = table[minNode->prev][i];
          if(part == PART_C && lst[i]->cost > antennaCost)
            lst[i] ->cost = antennaCost;
          enqueue(pq,  lst[i], lst[i]->cost);
          
          
        }
      }
      

    }
    int tot = 0;
    for(int i = 0; i < g->numVertices; i++){
      tot += lst[i]->cost;
        
    }
    
    return tot;
}


