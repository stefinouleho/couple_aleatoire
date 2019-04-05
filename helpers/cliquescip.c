#include "proglin_helper.h"
#include "clique.h"
#include "graph.h"
#include <stdlib.h>
#include <stdio.h>

int* clique_max(graph g, long timeout){
  /* Clique calculée avec un programme linéaire */

  
  /* Créer le programme  linéaire, n variables, une par noeud */ 

  int n = nbnodes(g);

  abpl pl = create_lp(n);
  
  int i,j;


  for(i = 0; i < n; i++){
    char name[12];
    sprintf(name, "x%d", i); /* en espérant qu'il y ait moins de 10^11 noeuds dans le graphe */
    add_variable_binary(pl, name);
  }

  /* Contraintes, une par non-arête : si deux noeuds u et v ne sont pas reliés, alors xu + xv <= 1 */


  for(i = 0; i < n; i++){
    for(j = i+1; j < n; j++){
      if(are_neighbors(g, i, j))
        continue;
      int colno[2] = {i, j};
      double row[2] = {1, 1};
      add_le_constraint(pl, 2, colno, row, 1); /* xi + xj <= 1*/
    }
  }



  /* Fonction objectif : maximiser le nombre de noeuds */

  int* colnoobj = (int*) malloc(sizeof(int)*n);
  double* rowobj = (double*) malloc(sizeof(double)*n);
  for(i = 0; i < n; i++){
    colnoobj[i] = i;
    rowobj[i] = 1;
  }

  set_maximization(pl);

  set_objective(pl, n, colnoobj, rowobj); /* max sum(x_u, u noeud de g) */

    

  /* Quantité de débug affiché à l'écran, ici on affiche rien */

  set_silent_verbose(pl);
  /* on indique le timeout en seconde */
  if(timeout != 0)
    set_max_timeout(pl, timeout);

  /* Résolution, ret est le code de retour, qui explique comment s'est terminé la résolution */

  int ret = solve_pl(pl);


  int* clique = (int*) malloc(sizeof(int)*n);

  /* Si ret indique qu'on a trouvé une solution optimale */
  if(ret == ABLP_OPTIMAL || ret == ABLP_SUBOPTIMAL){
    for(i = 0; i < n; i++){
      double x = get_variable_value(pl, i);
      if(x > 0.5)
        clique[i] = 1;
      else
        clique[i] = 0;
    }
  }
  else{
    for(i = 0; i < n; i++)
      clique[i] = -1;
  }

  /* Destruction du programme linéaire */

  destroy_lp(pl);

  return clique;

}

