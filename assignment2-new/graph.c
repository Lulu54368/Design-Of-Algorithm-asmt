/*
graph.c

Set of vertices and edges implementation.

Implementations for helper functions for graph construction and manipulation.

Skeleton written by Grady Fitzpatrick for COMP20007 Assignment 1 2021
*/
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include "graph.h"
#include "utils.h"
#include "pq.h"
#include"list.h"
#define INITIALEDGES 32

struct edge;

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
};
typedef struct vertice{

  int *neighbours;
  int numNeighbour;
}vertices;
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
struct node{
  int parent;
  int value;
  int marked;
};

int dfs(struct graph *g, int table[g->numVertices][g->numVertices], int * group, int * marked,  int *outages, int numOutages);
int dfsExplorer(int v, struct graph *g, int table[g->numVertices][g->numVertices], int *group, int *currgroup, int *count, int *marked);
int countNum(struct graph * g, int *group);
int getMaxGroup(struct graph * g, int *group, int *list, int currgroup);
int find_min(struct graph *g, int *marked, int group, int table[g->numVertices][g->numVertices], 
struct node *track[g->numVertices], int v);
void traceback(struct graph * g, struct node* track[g->numVertices],int * trace,int num,int  start);
void modifytable(struct graph * g, int table[g->numVertices][g->numVertices], int *outages, int numOutages);
int isancester(int end, int start, int *parent);
int findBackEdge(int *count, int origin, struct graph *g, vertices** verticeList,
int *parent, int * marked, int * critical);
void addEdge(struct graph *g, int start, int end){
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

/* Finds:
  - Number of connected subnetworks (before outage) (Task 2)
  - Number of servers in largest subnetwork (before outage) (Task 3)
  - SIDs of servers in largest subnetwork (before outage) (Task 3)
  - Diameter of largest subnetworks (after outage) (Task 4)
  - Number of servers in path with largest diameter - should be one more than
    Diameter if a path exists (after outage) (Task 4)
  - SIDs in path with largest diameter (after outage) (Task 4)
  - Number of critical servers (before outage) (Task 7)
  - SIDs of critical servers (before outage) (Task 7)
*/
struct solution *graphSolve(struct graph *g, enum problemPart part,
  int numServers, int numOutages, int *outages){
  struct solution *solution = (struct solution *)
    malloc(sizeof(struct solution));
  assert(solution);
  /* Initialise solution values */
  initaliseSolution(solution);
  
 

  
 
  int marked[g->numVertices];
  int group[g->numVertices];
  int i, j;
  for(i = 0; i < g->numVertices; i++){
    marked[i] = 0;
    group[i] = 0;
  }

  
  
  if(part == TASK_2){
    /* IMPLEMENT TASK 2 SOLUTION HERE */
    int table [g->numVertices][g->numVertices];
    for(i = 0; i < g->numVertices; i++){
      for(j = 0; j < g->numVertices; j++){
        table[i][j] = 0;
      }
    }
 
    for(i = 0; i < g->numEdges; i++){
      table[g->edgeList[i]->start][g->edgeList[i]->end] = 1;
      table[g->edgeList[i]->end][g->edgeList[i]->start] = 1;
    
    }
    int groupNum = dfs(g, table,group, marked, outages, numOutages);
    solution->connectedSubnets = groupNum;
  } else if(part == TASK_3) {
    int table [g->numVertices][g->numVertices];
    for(i = 0; i < g->numVertices; i++){
      for(j = 0; j < g->numVertices; j++){
        table[i][j] = 0;
      }
    }
    int *destroyList = NULL;
    int destroyNum = 0;
    for( i = 0; i < g->numEdges; i++){
      table[g->edgeList[i]->start][g->edgeList[i]->end] = 1;
      table[g->edgeList[i]->end][g->edgeList[i]->start] = 1;
    
    }
    int groupNum = dfs(g, table,group, marked, destroyList, destroyNum);
    int curr = 0, max = 0;
    int *currList, *maxList;
    maxList = malloc(g->numVertices *sizeof(int));
   
    
    for(i = 1; i <= groupNum; i ++){
      currList = malloc(g->numVertices * sizeof(int));
      curr = getMaxGroup( g, group, currList, i);
      if(max < curr){
        max = curr;
        maxList = currList;
      }
    }
   
    solution->largestSubnet = max;
    solution->largestSubnetSIDs = maxList;
  } else if(part == TASK_4) {
  
    int table [g->numVertices][g->numVertices];
    for(i = 0; i < g->numVertices; i++){
      for(j = 0; j < g->numVertices; j++){
        table[i][j] = 0;
      }
    }
    
    for( i = 0; i < g->numEdges; i++){
      table[g->edgeList[i]->start][g->edgeList[i]->end] = 1;
      table[g->edgeList[i]->end][g->edgeList[i]->start] = 1;
    
    }
   
    modifytable(g, table, outages, numOutages);
    
    for( i = 0; i < g->numVertices; i++){
      marked[i] = 0;
      group[i] = 0;
    }
    int numgroup;
    numgroup = dfs(g,table ,group,marked, outages, numOutages);

    
  
    struct node **alltrack[numgroup];
    
    
    for( i = 0; i < numgroup; i++){
      alltrack[i] = malloc(g->numVertices * sizeof(*alltrack));  
      for( j = 0; j < g->numVertices; j ++){
        alltrack[i][j] = malloc(sizeof(struct node) );
        alltrack[i][j] ->value = 0;
        alltrack[i][j] -> parent = i;
        alltrack[i][j]->marked = 0;
      }

    }
    
    
    int outageInd = 0;
    int currgroup = 1;
    for(i = 0; i < g->numVertices; i++){
      if (currgroup > numgroup)  break;
      else if (i == outages[outageInd]){
        if(outageInd +1 < numOutages) outageInd++;
        continue;
      }
      
    
      
      if(group[i] == currgroup){
        
        alltrack[currgroup-1][i]->marked = 1;
        int t;
        for(t = 0; t < g->numVertices ; t++){
          if(group[t] == currgroup)
            alltrack[currgroup-1][t]->value = g->numVertices;
        }
        alltrack[currgroup-1][i]-> value = 0;
        
        find_min(g, group, currgroup, table, alltrack[currgroup-1] , i); 
        
       
        currgroup ++;
      }
     
       
    }
    
    int maxGroup , maxstart, currmax = 0;
    for(i = 0; i < numgroup; i++){
      for(j = 0; j < g->numVertices; j++){
        if(alltrack[i][j]->value > currmax){
          maxGroup = i;
          maxstart = j;
          currmax = alltrack[i][j]->value;
        }
      }
    }
    int *trace;
    trace = malloc(sizeof(int) * (currmax+1));
    traceback(g, alltrack[maxGroup], trace, currmax,  maxstart);
   
    solution->postOutageDiameter =  currmax;
    solution->postOutageDiameterCount = currmax+1;
    solution->postOutageDiameterSIDs = trace;
  } else if(part == TASK_7) {

    vertices **verticeList;
    verticeList = malloc(sizeof(*verticeList) * g->numVertices);
   
    int edgeNum = g->numEdges;
    for(i = 0; i < g->numVertices; i++){
      
      verticeList[i] = malloc(sizeof(vertices));
      
      verticeList[i]->neighbours = malloc(sizeof(int));
      verticeList[i]->numNeighbour = 0;
    }
  
    for(i = 0; i < edgeNum; i++){
    
      int start = g->edgeList[i]->start;
      int end = g->edgeList[i]->end;
      
      verticeList[start]->neighbours[verticeList[start]->numNeighbour++] = end;
   
      verticeList[start]->neighbours = realloc( verticeList[start]->neighbours , (verticeList[start]->numNeighbour + 1)*sizeof(int));
      
      verticeList[end]->neighbours[verticeList[end]->numNeighbour++] = start;

      verticeList[end]->neighbours = realloc(verticeList[end]->neighbours , (verticeList[end]->numNeighbour + 1)*sizeof(int));
      
    }
  

  
  
    int parent[g->numVertices], critical[g->numVertices];
    for(i = 0; i < g->numVertices; i ++){
      parent[i] = 0;
      critical[i] = 0;
      marked[i] = 0;
    }
    int v;
    for(v = 0; v < g->numVertices; v ++){
      
      if(!marked[v]){
   
        int count = 0;
      
        findBackEdge(&count, v, g, verticeList, parent, marked, critical);
        if(critical[v] && verticeList[v]->numNeighbour <= 2){
        
          critical[v] = 0;
        }
        if(!critical[v] && verticeList[v]->numNeighbour > 2){
        
          critical[v] = 1;
        }
      }
      
        
        
          
      
    }
    int *final;
    int tot = 0;
    final = malloc(sizeof(int) * g->numVertices);
    for(i = 0; i < g->numVertices; i++){
      if(critical[i]){
        final[tot++] = i;
      }
    }
    solution->criticalServerCount = tot;
    solution->criticalServerSIDs = final;

  }
  return solution;
}
  

int dfs(struct graph *g, int table[g->numVertices][g->numVertices],int * group,int * marked, int *outages, int numOutages){
  int count = 1, currgroup = 0;
  int ind = 0, curr = -1;
  if(numOutages) curr = outages[ind];
  int v;
  for(v = 0; v < g->numVertices; v++){
    if(!marked[v] && (!numOutages || v != curr)){
      currgroup++;
    
      dfsExplorer(v, g, table, group, &currgroup,&count, marked);
    }
    else if(numOutages && v == curr){
    
      marked[v] = 0;
      group[v] = 0;
      ind++;
      curr = outages[ind];
      
    }
    else continue;
  }
  return currgroup;
}

int dfsExplorer(int v, struct graph *g, int table[g->numVertices][g->numVertices], int *group, int *currgroup, int *count, int *marked){
  
  marked[v] = *count;
  (*count)++;
  group[v] = *currgroup;
 
  int i;
  for(i = 0; i < g->numVertices; i++){
    if(table[v][i] && !marked[i]){
      
      dfsExplorer(i, g, table, group, currgroup, count, marked);
      
    }
  }
  return 0;
  
}

int getMaxGroup(struct graph * g, int *group, int *list, int currgroup){
  int Ind = 0;
  int i;
  for(i = 0; i < g->numVertices; i ++){
    if (group[i] == currgroup){
      list[Ind++] = i;
      
    }
  }
  return Ind;
}

int find_min(struct graph *g, int *marked, int group, int table[g->numVertices][g->numVertices], 
struct node *track[g->numVertices], int v){
  track[v]->marked = 1;
  int i;
 
 
  int num[g->numVertices];
  int order = 0;



  for(i = 0; i < g->numVertices; i++){
  
    if(table[v][i] + track[v]->value < track[i]->value  && table[v][i])
    {
      
        if(track[i]->marked)    {
         
         
            num[order] = i;
            order++;
        }
        else{
            track[i] ->parent = v;
            track[i]->value = table[v][i] + track[v]->value;
            
            find_min(g, marked, group, table, track, i);
        }
         for(i = 0; i < order; i++){
          
            track[num[i]] ->parent = v;
            track[num[i]]->value = table[v][num[i]] + track[v]->value;
           
            
        }
    }
      
    }
    
      
   
  
  return group;
}




void modifytable(struct graph * g, int table[g->numVertices][g->numVertices], int *outages, int numOutages){
  int Ind = 0;
  int i, j;
  while(Ind < numOutages){
    for( i = 0; i < numOutages; i++){
      table[outages[Ind]][i] = 0;
    }
    Ind++;
    
  }
  for(i = 0; i < g->numVertices; i++){
    for(j =0; j < numOutages; j ++)
      table[i][outages[j]] = 0;
  }
}

void traceback(struct graph * g, struct node* track[g->numVertices],int * trace,int num,int  start){
 
  int curr = start;
  struct list * backtrace;
  int i = 0;
  backtrace = newlist(track[start]);

  while(i < num-1){
    curr = track[curr]->parent;
    backtrace = prependList(backtrace, track[curr]);
 
    i++;
  }

  for(i = 0; i < num; i++){
  
    trace[i] = ((struct node*)deleteHead(&backtrace))->parent;
 
  }
  trace[num] = start;
}





int isancester(int end, int start, int *parent){
  if(parent[start]== end) return 0;
  int before;
  before = start;
  
  while(parent[before]!= before){
   
    if(before == end) {
      
      return 1;
    }
    else  before = parent[before];
  }
  if(before == end) {
      
      return 1;
  }
  return 0;
}

int findBackEdge(int *count, int origin, struct graph *g, vertices** verticeList,
int *parent, int * marked, int * critical){
 int min = 0;
  (*count)++;
  marked[origin] = *count;

  struct pq * pq = newPQ();
  int *node; 
  int prev, i;
  int totmin = g->numVertices;
  for(i = 0; i <verticeList[origin]->numNeighbour; i++){
      int w = verticeList[origin]->neighbours[i];
      if(marked[w] == 0){
        parent[w] = origin;
        prev = findBackEdge(count, w, g, verticeList, parent, marked, critical);
        
        if(prev < totmin) totmin = prev;
      }
      else if(isancester(w, origin, parent)){
        node = malloc(sizeof(int));
        *node = w;
        enqueue(pq, node, marked[w]);
      
      }
    
    
   
  }
  if(!empty(pq)) {
    min = *(int*)deletemin(pq);

    if(!empty(pq)|| prev > min){
      critical[min] = 1;
      
    
    }
    
  }
  else{
    if(totmin >= origin && totmin != g->numVertices){
      critical[origin] = 1;
      min = totmin;
    } 
    else if (totmin == g->numVertices)
      min = origin;
    else  min = totmin;
  }
  
  return min;
}