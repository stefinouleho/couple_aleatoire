/**

* Bibliotheque d'aide pour la création, la destruction et la manipulation de programmes linéaire

* Cette interface vous évite de savoir avec quel langage de programmation linéaire vous travaillez.

*

* @author Dimitri Watel

*/


#ifndef PROGLIN_HELPER_H

#define PROGLIN_HELPER_H





/**

* Type abstrait représentant un programme linéaire

*/

typedef struct prog_lin* abpl;







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

abpl create_lp(int nb_variables);





/**

* Destructeur de programme linéaire.

* @param pl Un programme linéaire déjà initialisé

* @require pl a été créé avec la fonction create_lp

*/

void destroy_lp(abpl pl);







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

char* get_variable_name(abpl pl, int variable_index);





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

int add_variable_binary(abpl pl, char* variable_name);





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

int add_variable_integer(abpl pl, int lower_bound, int upper_bound, char* variable_name);





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

int add_variable_real(abpl pl, double lower_bound, double upper_bound, char* variable_name);





/****************************  CONTRAINTES  ******************************/





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

*ATTENTION, dans un soucis de performance, il ne faut pas avoir définit la fonction objective avec set_objective avant d'avoir définit une contrainte.

*

*/

void add_le_constraint(abpl pl, int nb_variables, int* variables_indexes, double* variables_coeffs, double right_member);





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

*ATTENTION, dans un soucis de performance, il ne faut pas avoir définit la fonction objective avec set_objective avant d'avoir définit une contrainte.

*

*/

void add_ge_constraint(abpl pl, int nb_variables, int* variables_indexes, double* variables_coeffs, double right_member);





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

*ATTENTION, dans un soucis de performance, il ne faut pas avoir définit la fonction objective avec set_objective avant d'avoir définit une contrainte.

*

*/

void add_eq_constraint(abpl pl, int nb_variables, int* variables_indexes, double* variables_coeffs, double right_member);



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

*ATTENTION, dans un soucis de performance, il faut avoir définit toutes les contraintes avec les fonctions add_le_constraint, add_ge_constraint, add_eq_constraint

*

*/

void set_objective(abpl pl, int nb_variables, int* variables_indexes, double* variables_coeffs);





/**

* Définit le programme linéaire comme étant un problème de maximisation.

* @param pl Un programme linéaire déjà initialisé

* @require pl a été créé avec la fonction create_lp

* @require La fonction objective a déjà été définie avec la fonction set_objective

*/

void set_maximization(abpl pl);





/**

* Définit le programme linéaire comme étant un problème de minimisation.

* @param pl Un programme linéaire déjà initialisé

* @require pl a été créé avec la fonction create_lp

* @require La fonction objective a déjà été définie avec la fonction set_objective

*/

void set_minimization(abpl pl);





/****************************  VERBOSITE DU PROGRAMME  ******************************/



/**

* Définit le mode de verbosité du programme comme étant maximal.

* @param pl Un programme linéaire déjà initialisé

* @require pl a été créé avec la fonction create_lp

* 

* Lors de la résolution du programme, tous les messages seront affichés. Il est conseillé d'appeler cette fonction uniquement dans un but de débogage car les performances sont réduites. Par défaut, c'est le mode choisi.

*/

void set_debug_verbose(abpl pl);





/**

* Définit le mode de verbosité du programme comme étant minimal.

* @param pl Un programme linéaire déjà initialisé

* @require pl a été créé avec la fonction create_lp

* 

* Lors de la résolution du programme, aucun message ne sera affiché dans la sortie standard. Quelques messages seront affichés en sortie d'erreur. Par défaut, le mode de verbosité est le mode debug

*/

void set_silent_verbose(abpl pl);





/**

* Affiche le programme linéaire en console

* @param pl Un programme linéaire déjà initialisé

* @require pl a été créé avec la fonction create_lp

* 

* Affiche tout le programme linéaire en console. Il est conseillé de n'utiliser cette fonction que pour du debug.

*/

void print_proglin(abpl pl);





/****************************  RESOLUTION DU PROGRAMME  ******************************/





#define ABLP_OPTIMAL 0 /** < Code de retour de résolution renvoyé quand la solution optimale a été calculée */

#define ABLP_SUBOPTIMAL 1 /** < Code de retour de résolution renvoyé quand la solution optimale a été calculée */

#define ABLP_UNSOLVABLE 2 /** < Code de retour de résolution renvoyé quand une solution réalisable mais pas nécessairement optimale a été calculée, après un temps de résolution maximum */

#define ABLP_TIMEOUT 3 /** < Code de retour de résolution renvoyé quand aucune solution n'a été renvoyée après un temps de résolution maximum */

#define ABLP_ERROR -1 /** < Code de retour de résolution renvoyé quand une erreur s'est produite lors de la résolution */





/**

* Définit le timeout du programme linéaire

* @param pl Un programme linéaire déjà initialisé

* @param timeout Durée maximum en secondes après lesquelles le programme s'arête même s'il n'a pas fini de résoudre le programme. Si égal 0, aucun timeout n'est définit.

* @require pl a été créé avec la fonction create_lp

*

* Si timeout est non nul, définit une durée maximum en secondes après lesquelles le programme s'arrête même s'il n'a pas fini de résoudre le programme. Il renvoie alors la meilleure solution qu'il a trouvé. Si timeout est nul, le programme ne s'arrête que lorsqu'il a trouvé une solution ou décrété que le programme ne peut être résolu.

*/

void set_max_timeout(abpl pl, long timeout);





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

int solve_pl(abpl pl);





/******************* SORTIES ***************************/





/**

* Renvoie la meilleure solution trouvée pour le programme linéaire

* @param pl Un programme linéaire déjà initialisé

* @require pl a été créé avec la fonction create_lp

* @require pl a été résolu avec la fonction solve et un code de retour ABLP_OPTIMAL ou ABLP_SUBOPTIMAL a été renvoyé

*

* Renvoie la meilleure solution trouvée pour le programme linéaire pl par la fonction solve. Renvoie la solution optimale si ABLP_OPTIMAL a été renvoyé par la fonctin solve, ou la meilleure solution réalisable si ABLP_SUBOPTIMAL a été renvoyé. 

*/

double get_optimal_solution(abpl pl);





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

double get_variable_value(abpl pl, int variable_index);












/**

* Renvoie le nombre d'itérations du programme lors de sa résolution.

* @param pl Un programme linéaire déjà initialisé

* @require pl a été créé avec la fonction create_lp

* @require pl a été résolu avec la fonction solve et un code de retour ABLP_OPTIMAL ou ABLP_SUBOPTIMAL a été renvoyé

*

* Renvoie le nombre d'itérations du programme.

*/

int get_nb_iterations(abpl pl);





/**

* Renvoie le nombre de noeuds de Branch and Bounds visités par le programme lors de sa résolution.

* @param pl Un programme linéaire déjà initialisé

* @require pl a été créé avec la fonction create_lp

* @require pl a été résolu avec la fonction solve et un code de retour ABLP_OPTIMAL ou ABLP_SUBOPTIMAL a été renvoyé

*

* Renvoie le nombre de noeuds de Branch and Bounds.

*/

int get_nb_nodes(abpl pl);


#endif
