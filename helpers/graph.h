
#ifndef GRAPH_H

#define GRAPH_H


/**

* Bibliotheque d'aide pour la manipulation de graphes

* @author Dimitri Watel

*/

/**
* Type abstrait représentant un graphe non orienté
*/

typedef struct graph_s* graph;


/**
* Crée un graphe à partir de la matrice d'adjacence* @param n : nombre de noeuds du graphe 
* @param links : matrice d'ajacence du graphe 
* @return un graphe dont la matrice d'adjacence est links et le nombre de noeuds est n
* @require links est de taille nxn et ne contient que des 1 et des 0
*/
graph build_graph_from_matrix(int n, int** links);

/**
* Renvoie le nombre de noeuds du graphe g
* @param g : un graphe
* @return le nombre de noeud du graphe g
*/
int nbnodes(graph g);

/**
* Renvoie la matrice d'adjacence du graphe g
* @param g : un graphe
* @return la matrice d'adjacence du graphe g
*/
int** build_matrix_from_graph(graph g);

/**
* Destructeur du graphe g
* @param g : un graphe
*/
void destroy(graph g);

/**
* Vérifie si u et v sont voisins
* @param g : un graphe
* @param u : un noeud de g
* @param v : un noeud de g
* @return 1 si u et v sont reliés dans g et 0 sinon
* @require 0 <= u,v < n
*/ 
int are_neighbors(graph g, int u, int v);

#endif
