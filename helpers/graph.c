#include "graph.h"
#include <stdlib.h>
#include <stdio.h>
/**
* Type concret représentant un graphe
*/

struct graph_s{
  int n; /* Nombre de noeuds */
  int** links; /* Matrice d'adjacence */
};

graph build_graph_from_matrix(int n, int** links){
  graph g = (graph) malloc(sizeof(struct graph_s));
  g->n = n;
  
  /* Attention à ne pas réécrire la matrice */
  g->links = links;
  
  /* Version safe
  int** linkscp = (int**) malloc(sizeof(int*) * n);
  int i,j;
  for(i = 0; i < n; i++){
    linkscp[i] = (int*) malloc(sizeof(int) * n);
    for(j = 0; j < n; j++){
      linkscp[i][j] = links[i][j];
    }
  }
  g->links = linkscp;
  */
  return g;
}

int nbnodes(graph g){
  return g->n;
}

int** build_matrix_from_graph(graph g){

  return g->links;
   
  /* Version safe
  int n = nbnodes(g);
 
  int** linkscp = (int**) malloc(sizeof(int*) * n);
  int i,j;
  for(i = 0; i < n; i++){
    linkscp[i] = (int*) malloc(sizeof(int) * n);
    for(j = 0; j < n; j++){
      linkscp[i][j] = g->links[i][j];
    }
  }
  return linkscp;
  */
}

void destroy(graph g){

  /* Si links a été créé lors de la création du graphe
  int n = nbnodes(g);
  int i;
  for(i = 0; i < n; i++){
    free(g->links[i]);
  } 
  free(g->links);
  */
  free(g);
}

int are_neighbors(graph g, int u, int v){
  return g->links[u][v];
}

