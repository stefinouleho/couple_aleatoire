#include "fonctions_mces.h"

double last_chrono;
int main(int argc, char *argv[])
{
	
	
	if( argc != 5){
		fprintf(stdout,"Missing arguments (num chebi1 , num chebi 2 ( 0  for all) , date , lenght limit )\n");
		exit(20);
	}

	//lecture des molecules dans le fichier chebi_lite.sdf
	struct molecule *M = lecture_fichier_chebi();
	
	//recuperation des parametres
	int num_chebi1, num_chebi2,taille_limite;
	float date;
	num_chebi1 = atoi(argv[1]);
	num_chebi2 = atoi(argv[2]);
	date 	   = atof(argv[3]); //temps max de calcul en secondes
	taille_limite = atoi(argv[4]); // la taille max du graphe produit
	
	fprintf(stdout," donnees recus :\n molecule de chebi 1 = %d et chebi 2 = %d  avec date limite de  %f s et taille limite de %d sommets\n" , num_chebi1,num_chebi2,date,taille_limite);
	
	float r;
	if( num_chebi2 != 0)
	{
		last_chrono = chrono();
		r = mesure_similarite( num_chebi1, num_chebi2,M, date, taille_limite);
		fprintf(stdout,"similarite mesure entre %d et %d est de : %f\n",num_chebi1,num_chebi2,r);
	}
	
	else
	{
		similarite_all(num_chebi1, M, date, taille_limite);
	}	
	
	int nb_mol;
	printf("3. Libération de la mémoire : %.3lf s\n",chrono());

	for(nb_mol=0 ; nb_mol < NB_MOLECULES ; nb_mol++) 
		liberer_molecule(M[nb_mol]);
	free(M);

	exit(0);
}
