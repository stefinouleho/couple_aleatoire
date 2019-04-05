/*

   NE PAS MODIFIER CE FICHIER

*/



#include "scip/scip.h"

#include "scip/scipdefplugins.h"

#include "proglin_helper.h"



/**

* Type concret représentant un programme linéaire

*/

struct prog_lin{

  SCIP* scip;

  int nbv; /** < Nombre de variables du programme linéaire, sert à vérifier les indices des variables en entrée */

  int current_var_index;

  int current_cons_index;

  SCIP_VAR** variables;

  SCIP_SOL* solution;

};





/****************************  CREATE DESTROY  ******************************/



/**

* Constructeur de programme linéaire.

* @param nb_variables Nombre de variables dans le programme linéaire.

* @return un programme linéaire sans contrainte sans fonction objectif avec nb_variables variables

* @require nb_variables > 0

*

* Il faut appeler cette méthode pour créer une variable

* de type abpl et pouvoir l'utiliser. 

*/

abpl create_lp(int nb_variables){

  SCIP* scip;

  SCIPcreate(&scip);

  SCIPincludeDefaultPlugins(scip);

  SCIPcreateProb(scip, "Programme linéaire Puzzle", NULL, NULL,

NULL, NULL, NULL, NULL, NULL);

  

  abpl pl = (abpl) malloc(sizeof(struct prog_lin));

  pl->scip = scip;

  pl->variables = (SCIP_VAR**) malloc(nb_variables * sizeof(SCIP_VAR*));

  pl->nbv = nb_variables;

  pl->current_var_index = 0;

  pl->solution = NULL;

  set_debug_verbose(pl);

  return pl;

}



/**

* Destructeur de programme linéaire.

* @param pl Un programme linéaire déjà initialisé

* @require pl a été créé avec la fonction create_lp

*/

void destroy_lp(abpl pl){

  SCIPfree(&(pl->scip));

  free(pl->variables);

  free(pl);

}





/****************************  VARIABLES  ******************************/





/**

* Renvoie le nom d'une variables

* @param pl Un programme linéaire déjà initialisé

* @param variable_index Numéro d'une variable du programme

* @require variable_index > 0 et variable_index <= nombre de variables du programme

* @require pl a été créé avec la fonction create_lp

* @require la fonction set_variable_name a déjà été appelée au moins une fois avec les paramêtres pl et variable_index

*

* Renvoie le nom associé à la variable dont l'indice est variable_index dans le programme pl par la fonction set_variable_name.

*/

char* get_variable_name(abpl pl, int variable_index){

  return (char*)SCIPvarGetName(pl->variables[variable_index]);

}





int add_variable_type(abpl pl, double lower_bound, double upper_bound, char* variable_name, SCIP_VARTYPE type){

  SCIP_VAR* variable = pl->variables[pl->current_var_index];

  SCIPcreateVarBasic(pl->scip, &variable, variable_name, lower_bound, upper_bound, 0.0, type);

  SCIPaddVar(pl->scip, variable);

  pl->variables[pl->current_var_index] = variable;

  pl->current_var_index++;

  return pl->current_var_index - 1;

}





/**

* Ajoute une variable binaire et renvoie son indice.

* @param pl Un programme linéaire déjà initialisé

* @param variable_name Nom de la variable

* @return L'indice de la variable

* @require pl a été créé avec la fonction create_lp

* @require variable_name n'est pas le nom d'une autre variable existance

*

* Crée une nouvelle variable binaire (entier compris entre 0 et 1) dont le nom est variable_name. Renvoie l'indice de cette variable qui servira a identifier les variables pour créer les contraintes et la fonction objectif. La première variable est d'indice 0, la seconde d'indice 1, ...

*/

int add_variable_binary(abpl pl, char* variable_name){

  return add_variable_type(pl, 0.0, 1.0, variable_name, SCIP_VARTYPE_BINARY);

}





/**

* Ajoute une variable entière, définit ses bornes et renvoie son indice.

* @param pl Un programme linéaire déjà initialisé

* @param lower_bound Borne inférieure

* @param upper_bound Borne supérieure

* @param variable_name Nom de la variable

* @return L'indice de la variable

* @require pl a été créé avec la fonction create_lp

* @require lower_bound <= upper_bound

* @require variable_name n'est pas le nom d'une autre variable existance

*

* Crée une nouvelle variable entière comprise entre lower_bound et upper_bound dont le nom est variable_name. Renvoie l'indice de cette variable qui servira a identifier les variables pour créer les contraintes et la fonction objectif. La première variable est d'indice 0, la seconde d'indice 1, ...

*/

int add_variable_integer(abpl pl, int lower_bound, int upper_bound, char* variable_name){

  return add_variable_type(pl, lower_bound, upper_bound, variable_name, SCIP_VARTYPE_INTEGER);

}





/**

* Ajoute une variable réelle, définit ses bornes et renvoie son indice.

* @param pl Un programme linéaire déjà initialisé

* @param lower_bound Borne inférieure

* @param upper_bound Borne supérieure

* @param variable_name Nom de la variable

* @return L'indice de la variable

* @require pl a été créé avec la fonction create_lp

* @require lower_bound <= upper_bound

* @require variable_name n'est pas le nom d'une autre variable existance

*

* Crée une nouvelle variable réelle comprise entre lower_bound et upper_bound dont le nom est variable_name. Renvoie l'indice de cette variable qui servira a identifier les variables pour créer les contraintes et la fonction objectif. La première variable est d'indice 0, la seconde d'indice 1, ...

*/

int add_variable_real(abpl pl, double lower_bound, double upper_bound, char* variable_name){

  return add_variable_type(pl, lower_bound, upper_bound, variable_name, SCIP_VARTYPE_CONTINUOUS);

}







/****************************  CONTRAINTES  ******************************/





void add_constraint_to_pl(abpl pl, int nb_variables, int* variables_indexes, SCIP_Real* variables_coeffs, SCIP_Real left_member, SCIP_Real right_member){

  SCIP_CONS* cons;

  SCIP_VAR ** variables = (SCIP_VAR**) malloc(nb_variables * sizeof(SCIP_VAR*));

  int i;

  for(i = 0; i<nb_variables; i++){

    variables[i] = pl->variables[variables_indexes[i]];

  }

  

  char name[30];

  sprintf(name, "C%d", pl->current_cons_index);

  pl->current_cons_index++;

  

  SCIPcreateConsLinear(pl->scip, &cons, name, nb_variables, variables, variables_coeffs, left_member, right_member, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);

  SCIPaddCons(pl->scip, cons);

  free(variables);

}





/**

* Définit une contrainte de type "<="

* @param pl Un programme linéaire déjà initialisé

* @param nb_variables Nombre de variables avec un coefficient non nul dans la contrainte

* @param variables_indexes Indices des variables avec un coefficient non nul dans la contrainte

* @param variables_coeffs Coefficients des variables dont les indices sont définits dans variables_indexes

* @param right_member Membre de droite de la contrainte

* @require pl a été créé avec la fonction create_lp

* @require La fonction set_objective n'a pas été appelée avec le paramètre pl

* @require variables_indexes et variables_coeffs ont la même taille égale à nb_variables

* @require Pour tout index de variables_indexes, index > 0 et index <= nombre de variables du programme

*

* Spécifie une contrainte de type "<=" du programme pl. Seules les variables dont le coefficient est non nul doivent être définies (cela améliore les performances). Le tableau variables_indexes contient les indices de ces variables et le tableau variables_coeffs les coefficients de ces variables dans le même ordre. right_member définit le nombre à droite de "<=".

*

* Par exemple, si vous voulez créer la contrainte 10x + 15z <= 12 dans un programme avec 3 variables x, y, z, alors, il vous faudra appeler la fonction avec ces paramètres :

add_le_constraint(pl, 2, [1,3], [10,15], 12);

*

* ATTENTION, dans un soucis de performance, il ne faut pas avoir définit la fonction objective avec set_objective avant d'avoir définit une contrainte.

*

*/

void add_le_constraint(abpl pl, int nb_variables, int* variables_indexes, double* variables_coeffs, double right_member){

  add_constraint_to_pl(pl, nb_variables, variables_indexes, variables_coeffs, -1 * SCIPinfinity(pl->scip), right_member);

}





/**

* Définit une contrainte de type ">="

* @param pl Un programme linéaire déjà initialisé

* @param nb_variables Nombre de variables avec un coefficient non nul dans la contrainte

* @param variables_indexes Indices des variables avec un coefficient non nul dans la contrainte

* @param variables_coeffs Coefficients des variables dont les indices sont définits dans variables_indexes

* @param right_member Membre de droite de la contrainte

* @require pl a été créé avec la fonction create_lp

* @require La fonction set_objective n'a pas été appelée avec le paramètre pl

* @require variables_indexes et variables_coeffs ont la même taille égale à nb_variables

* @require Pour tout index de variables_indexes, index > 0 et index <= nombre de variables du programme

*

* Spécifie une contrainte de type ">=" du programme pl. Seules les variables dont le coefficient est non nul doivent être définies (cela améliore les performances). Le tableau variables_indexes contient les indices de ces variables et le tableau variables_coeffs les coefficients de ces variables dans le même ordre. right_member définit le nombre à droite de ">=".

*

* Par exemple, si vous voulez créer la contrainte 10x + 15z >= 12 dans un programme avec 3 variables x, y, z, alors, il vous faudra appeler la fonction avec ces paramètres :

add_ge_constraint(pl, 2, [1,3], [10,15], 12);

*

* ATTENTION, dans un soucis de performance, il ne faut pas avoir définit la fonction objective avec set_objective avant d'avoir définit une contrainte.

*

*/

void add_ge_constraint(abpl pl, int nb_variables, int* variables_indexes, double* variables_coeffs, double right_member){

  add_constraint_to_pl(pl, nb_variables, variables_indexes, variables_coeffs, right_member, SCIPinfinity(pl->scip));

}





/**

* Définit une contrainte de type "="

* @param pl Un programme linéaire déjà initialisé

* @param nb_variables Nombre de variables avec un coefficient non nul dans la contrainte

* @param variables_indexes Indices des variables avec un coefficient non nul dans la contrainte

* @param variables_coeffs Coefficients des variables dont les indices sont définits dans variables_indexes

* @param right_member Membre de droite de la contrainte

* @require pl a été créé avec la fonction create_lp

* @require La fonction set_objective n'a pas été appelée avec le paramètre pl

* @require variables_indexes et variables_coeffs ont la même taille égale à nb_variables

* @require Pour tout index de variables_indexes, index > 0 et index <= nombre de variables du programme

*

* Spécifie une contrainte de type "=" du programme pl. Seules les variables dont le coefficient est non nul doivent être définies (cela améliore les performances). Le tableau variables_indexes contient les indices de ces variables et le tableau variables_coeffs les coefficients de ces variables dans le même ordre. right_member définit le nombre à droite de "=".

*

* Par exemple, si vous voulez créer la contrainte 10x + 15z = 12 dans un programme avec 3 variables x, y, z, alors, il vous faudra appeler la fonction avec ces paramètres :

add_eq_constraint(pl, 2, [1,3], [10,15], 12);

*

* ATTENTION, dans un soucis de performance, il ne faut pas avoir définit la fonction objective avec set_objective avant d'avoir définit une contrainte.

*

*/

void add_eq_constraint(abpl pl, int nb_variables, int* variables_indexes, double* variables_coeffs, double right_member){

  add_constraint_to_pl(pl, nb_variables, variables_indexes, variables_coeffs, right_member, right_member);

}





/****************************  FONCTION OBJECTIF  ******************************/





/**

* Définit la fonction objective

* @param pl Un programme linéaire déjà initialisé

* @param nb_variables Nombre de variables avec un coefficient non nul dans la contrainte

* @param variables_indexes Indices des variables avec un coefficient non nul dans la contrainte

* @param variables_coeffs Coefficients des variables dont les indices sont définits dans variables_indexes

* @require pl a été créé avec la fonction create_lp

* @require Toutes les contraintes ont déjà été créées avec les fonctions add_le_constraint, add_ge_constraint et add_eq_constraint

* @require variables_indexes et variables_coeffs ont la même taille égale à nb_variables

* @require Pour tout index de variables_indexes, index > 0 et index <= nombre de variables du programme

*

* Spécifie la fonction objective du programme pl. Seules les variables dont le coefficient est non nul doivent être définies (cela améliore les performances). Le tableau variables_indexes contient les indices de ces variables et le tableau variables_coeffs les coefficients de ces variables dans le même ordre.

*

* Par exemple, si vous voulez créer la fonction objective 10x + 15z dans un programme avec 3 variables x, y, z, alors, il vous faudra appeler la fonction avec ces paramètres :

set_objective(pl, 2, [1,3], [10,15]);

*

*ATTENTION, dans un soucis de performance, il ne faut pas avoir définit toutes les contraintes avec les fonctions add_le_constraint, add_ge_constraint, add_eq_constraint

*

*/

void set_objective(abpl pl, int nb_variables, int* variables_indexes, double* variables_coeffs){

  int i;

  for(i = 0; i < nb_variables; i++){

    int index = variables_indexes[i];

    double coeff = variables_coeffs[i];

    SCIP_VAR* variable = pl->variables[index];

    SCIPchgVarObj(pl->scip, variable, coeff);

  }

}





/**

* Définit le programme linéaire comme étant un problème de maximisation.

* @param pl Un programme linéaire déjà initialisé

* @require pl a été créé avec la fonction create_lp

* @require La fonction objective a déjà été définie avec la fonction set_objective

*/

void set_maximization(abpl pl){

  SCIPsetObjsense(pl->scip, SCIP_OBJSENSE_MAXIMIZE);

}





/**

* Définit le programme linéaire comme étant un problème de minimisation.

* @param pl Un programme linéaire déjà initialisé

* @require pl a été créé avec la fonction create_lp

* @require La fonction objective a déjà été définie avec la fonction set_objective

*/

void set_minimization(abpl pl){

  SCIPsetObjsense(pl->scip, SCIP_OBJSENSE_MINIMIZE);

}





/****************************  VERBOSITE DU PROGRAMME  ******************************/





/**

* Définit le mode de verbosité du programme comme étant maximal.

* @param pl Un programme linéaire déjà initialisé

* @require pl a été créé avec la fonction create_lp

* 

* Lors de la résolution du programme, tous les messages seront affichés. Il est conseillé d'appeler cette fonction uniquement dans un but de débogage car les performances sont réduites. Par défaut, le mode de verbosité est le mode normal.

*/

void set_debug_verbose(abpl pl){

  SCIPenableDebugSol(pl->scip);
  SCIPsetIntParam(pl->scip, "display/verblevel", 4);

}





/**

* Définit le mode de verbosité du programme comme étant minimal.

* @param pl Un programme linéaire déjà initialisé

* @require pl a été créé avec la fonction create_lp

* 

* Lors de la résolution du programme, aucun message ne sera affiché dans la sortie standard. Quelques messages seront affichés en sortie d'erreur. Par défaut, le mode de verbosité est le mode normal.

*/

void set_silent_verbose(abpl pl){

  SCIPdisableDebugSol(pl->scip);
  SCIPsetIntParam(pl->scip, "display/verblevel", 0);
}





/**

* Affiche le programme linéaire en console

* @param pl Un programme linéaire déjà initialisé

* @require pl a été créé avec la fonction create_lp

* 

* Affiche tout le programme linéaire en console. Il est conseillé de n'utiliser cette fonction que pour du debug.

*/

void print_proglin(abpl pl){

  SCIPprintOrigProblem(pl->scip, stdout, NULL, FALSE);

}





/****************************  RESOLUTION DU PROGRAMME  ******************************/





/**

* Définit le timeout du programme linéaire

* @param pl Un programme linéaire déjà initialisé

* @param timeout Durée maximum en secondes après lesquelles le programme s'arête même s'il n'a pas fini de résoudre le programme. Si égal 0, aucun timeout n'est définit.

* @require pl a été créé avec la fonction create_lp

*

* Si timeout est non nul, définit une durée maximum en secondes après lesquelles le programme s'arrête même s'il n'a pas fini de résoudre le programme. Il renvoie alors la meilleure solution qu'il a trouvé. Si timeout est nul, le programme ne s'arrête que lorsqu'il a trouvé une solution ou décrété que le programme ne peut être résolu.

*/

void set_max_timeout(abpl pl, long timeout){

  SCIPsetRealParam(pl->scip, "limits/time", timeout);

}







/**

* Résoud le programme linéaire

* @param pl Un programme linéaire déjà initialisé

* @return un code de retour en fonction de ce qui est renvoyé par le programme 

* @require pl a été créé avec la fonction create_lp

* @require Toutes les contraintes et la fonction objective de pl ont été définies

* 

* Cette fonction résoud le programme linéaire et renvoie un code de retour.

*    ABLP_OPTIMAL si une solution optimale est renvoyée, 

*    ABLP_SUBOPTIMAL en cas de timeout si une solution non optimale n'est pas renvoyée, 

*    ABLP_UNSOLVABLE si le problème n'a pas de solution optimale (non borné ou non réalisable), 

*    ABLP_TIMEOUT en cas de timeout si aucune solution réalisable n'a été trouvée, 

*    ABLP_ERROR sinon

*/

int solve_pl(abpl pl){

  SCIPsolve(pl->scip);

  SCIP_STATUS status = SCIPgetStatus(pl->scip);

  pl->solution = SCIPgetBestSol(pl->scip);

  switch(status) {

    case SCIP_STATUS_OPTIMAL:

      return ABLP_OPTIMAL;

    case SCIP_STATUS_INFEASIBLE:

      return ABLP_UNSOLVABLE;

    case SCIP_STATUS_UNBOUNDED:

      return ABLP_UNSOLVABLE;

    case SCIP_STATUS_INFORUNBD:

      return ABLP_UNSOLVABLE;

    case SCIP_STATUS_TIMELIMIT:

      if(pl->solution == NULL)

        return ABLP_TIMEOUT;

      else

        return ABLP_SUBOPTIMAL;

    default: 

      return ABLP_ERROR;

  }

}





/**

* Renvoie la meilleure solution trouvée pour le programme linéaire

* @param pl Un programme linéaire déjà initialisé

* @require pl a été créé avec la fonction create_lp

* @require pl a été résolu avec la fonction solve et un code de retour ABLP_OPTIMAL ou ABLP_SUBOPTIMAL a été renvoyé

*

* Renvoie la meilleure solution trouvée pour le programme linéaire pl par la fonction solve. Renvoie la solution optimale si ABLP_OPTIMAL a été renvoyé par la fonctin solve, ou la meilleure solution réalisable si ABLP_SUBOPTIMAL a été renvoyé. 

*/

double get_optimal_solution(abpl pl){

  return SCIPgetSolOrigObj(pl->scip, pl->solution);

}





/**

* Renvoie la valeur d'une variable de la meilleure solution trouvée pour le programme linéaire

* @param pl Un programme linéaire déjà initialisé

* @param variable_index Numéro d'une variable du programme

* @require variable_index > 0 et variable_index <= nombre de variables du programme

* @require pl a été créé avec la fonction create_lp

* @require pl a été résolu avec la fonction solve et un code de retour ABLP_OPTIMAL ou ABLP_SUBOPTIMAL a été renvoyé

*

* Renvoie la valeur d'une variable d'indice variable_index dans la meilleure solution trouvée pour le programme linéaire pl par la fonction solve. Renvoie la valeur de la variable dans la solution optimale si ABLP_OPTIMAL a été renvoyé par la fonctin solve, ou la valeur dans la meilleure solution réalisable si ABLP_SUBOPTIMAL a été renvoyé. 

*/

double get_variable_value(abpl pl, int variable_index){

  SCIP_VAR* variable = pl->variables[variable_index];

  return SCIPgetSolVal(pl->scip, pl->solution, variable);

}





/******************* SORTIES ***************************/





/**

* Renvoie le nombre d'itérations du programme lors de sa résolution.

* @param pl Un programme linéaire déjà initialisé

* @require pl a été créé avec la fonction create_lp

* @require pl a été résolu avec la fonction solve et un code de retour ABLP_OPTIMAL ou ABLP_SUBOPTIMAL a été renvoyé

*

* Renvoie le nombre d'itérations du programme.

*/

int get_nb_iterations(abpl pl){

  return SCIPgetNLPIterations(pl->scip);

}





/**

* Renvoie le nombre de noeuds de Branch and Bounds visités par le programme lors de sa résolution.

* @param pl Un programme linéaire déjà initialisé

* @require pl a été créé avec la fonction create_lp

* @require pl a été résolu avec la fonction solve et un code de retour ABLP_OPTIMAL ou ABLP_SUBOPTIMAL a été renvoyé

*

* Renvoie le nombre de noeuds de Branch and Bounds.

*/

int get_nb_nodes(abpl pl){

  return SCIPgetNStrongbranchs(pl->scip);

}
