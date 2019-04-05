#include "fonctions_mces.h"
#include "helpers/graph.h"
#include "helpers/clique.h"

#define AUCUNE_LIAISON (-1024)

int liaison_max = 0;
double last_chrono;

double chrono() 
{
	struct timeval tv;
	static double date_deb = 0.0;
	double date_courante;
	gettimeofday(&tv,NULL);
	date_courante = tv.tv_sec + ((double)tv.tv_usec)*1e-6;
	if (date_deb == 0.0) date_deb = date_courante;
	return date_courante-date_deb;
}

int position_M( int g1_chebi,struct molecule *M)
{ // trouve la poition d'une molecule de chebi g1_chebi dans M
	int i;
	
	for(i=0;i< NB_MOLECULES ;i++)
	{
		if(M[i].chebi_id == g1_chebi)
			return i;
	}
	if(i == NB_MOLECULES)
	{ 
		fprintf(stderr,"numero de chebi non present dans la base \n");
		exit(5); 
	}
	
	return 0;
}



struct molecule construction_matrice_mol(struct molecule m)
{//construction de la matrice de liaison d'une molecule
	
	int i,j;
	if(m.matrice_liaisons == NULL)
	{
		m.matrice_liaisons =  malloc(m.nb_atomes * sizeof(int *));
		
		for(i=0;i< m.nb_atomes;i++) m.matrice_liaisons[i] =  malloc(m.nb_atomes * sizeof(int));
		
		for(i=0;i< m.nb_atomes;i++)
		{
			for(j=0;j< m.nb_atomes;j++)
			 m.matrice_liaisons[i][j] = AUCUNE_LIAISON;
		}

		for(i =0; i< m.nb_liaisons;i++)
		{
			m.matrice_liaisons[m.liste_liaisons[i].A1-1][m.liste_liaisons[i].A2-1]= m.liste_liaisons[i].l_type;
			m.matrice_liaisons[m.liste_liaisons[i].A2-1][m.liste_liaisons[i].A1-1]= m.liste_liaisons[i].l_type;
		}
		

	}

	return m;
	
}


void affiche_matrice(struct molecule m)
{//affcihe la matrice d'une molecule m
	int i,j;
	printf("Affichage de la matrice m \n");
	
	if( m.matrice_liaisons == 	NULL)
	{
		printf("La matrice de cette molecule n'a pas encore eté défini \n");
		return;
	}
	for(i=0;i< m.nb_atomes;i++)
	{
		for(j=0;j< m.nb_atomes;j++)
			printf("%d ",m.matrice_liaisons[i][j] );
		
		printf("\n");
	}
	
}

struct couple *construction_couples(struct molecule *M,int pos1, int pos2,int taille)
{//Construction des couples d'atomes compatibles
	
	int n = 0,i,j;
	
	struct couple *couple_at;
	couple_at = malloc(taille * sizeof(couple));

	for(i= 0; i < M[pos1].nb_atomes;i++)
	{ 
		for(j= 0; j < M[pos2].nb_atomes;j++)
		{
			if( M[pos1].liste_atomes[i] == M[pos2].liste_atomes[j])
			{
				couple_at[n].a1 = i;
				couple_at[n].a2 = j;
				n++;
			}
		}
	}
	return couple_at;

}



graph graphe_produit(int g1_chebi,int g2_chebi,struct molecule *M)
{ //prend en entrée les chebi id de deux molecules  et contruit le graphe produit
		
		
	//trouve la position des molecules g1 et g2
	int taille= 0,pos1,pos2;
	pos1 =position_M(g1_chebi,M);
	pos2 =position_M(g2_chebi,M);
	
	//calcul de la taille du graphe produit
	int i,j;
	for(i= 0; i < M[pos1].nb_atomes; i++) 
		for(j= 0; j < M[pos2].nb_atomes;j++)
			if( M[pos1].liste_atomes[i] == M[pos2].liste_atomes[j]) 
        taille ++;
	
	//couple de liaisons entre les nouveaux sommets
	struct couple * couple_atome = construction_couples(M,pos1,pos2,taille);
	
	//construction de la matrice de liaison d'une molecule
	M[pos1]= construction_matrice_mol(M[pos1]);
	M[pos2]=construction_matrice_mol(M[pos2]);
	
    

	//initialisation de g12
	/*if( taille > 0)
	{
		g12.liste_atomes = malloc(g12.nb_atomes * sizeof(int));
		if(g12.liste_atomes == NULL)
		{
			fprintf(stdout, "Cannot allocate memory \n" );
			exit(33);
		}
	}*/
	
  int** matrice_liaisons;

	if( taille > 0)
	{	
		matrice_liaisons =  malloc(taille * sizeof(int *));	
		for(i = 0 ;i < taille;i++) 
      matrice_liaisons[i] =  malloc(taille * sizeof(int));
	}
	
  //remplissage des liaisons
	int i1,i2,j1,j2;
	for(i= 0; i < taille ;i++)
		for(j= 0; j < taille ;j++)
			matrice_liaisons[i][j] = 0;

	for(i= 0; i < taille ;i++)
		matrice_liaisons[i][i] = 1;

	for(i= 0; i < taille ;i++){ 
		i1 = couple_atome[i].a1;
		i2 = couple_atome[i].a2;
		for(j= i + 1 ; j < taille ;j++){
			j1=couple_atome[j].a1;
			j2=couple_atome[j].a2;
			
			if( M[pos1].matrice_liaisons[i1][j1] == M[pos2].matrice_liaisons[i2][j2] ){
				if(((i1 == j1) && (M[pos2].matrice_liaisons[i2][j2] != AUCUNE_LIAISON)) || ((i2 == j2) && (M[pos1].matrice_liaisons[i1][j1] != AUCUNE_LIAISON)) ||( (i1 != j1) && (i2!=j2))||( (i1 ==j1) && (i2==j2)) )
				{
					matrice_liaisons[i][j] = 1;
					matrice_liaisons[j][i] = 1;					
				}
			}
		}
	}
	free(couple_atome);
	//affiche_matrice(g12);
	return build_graph_from_matrix(taille, matrice_liaisons);
}

void liberer_molecule(struct molecule g) 
{ //liberation de l'espace memoire d'une molecule
	if (g.liste_atomes != NULL)
	{
		free(g.liste_atomes);
	}
	if (g.liste_liaisons != NULL ) free(g.liste_liaisons);
	if (g.matrice_liaisons)
	{	
		int i;
		for (i=0 ; i<g.nb_atomes ; i++) free(g.matrice_liaisons[i]);
		free(g.matrice_liaisons);
	}
} 


int*  graphe_g12(graph g12, int* clique_max,  struct molecule *M, int g1_chebi, int g2_chebi)
{ //contruction du graphe commun

	int* taille_graphe_commun = (int*) malloc(sizeof(int)*2);
  int nb_at =0,nb_liaisons=0,i,j,i1,j1;
	int pos1,pos2;
	pos1 =position_M(g1_chebi,M);
	pos2 =position_M(g2_chebi,M);
	
  int taille = nbnodes(g12);
	struct couple *couple_atome = construction_couples(M,pos1,pos2,taille);

	
	int tab[taille];
	for(i=0;i < taille ; i++)
	{
		tab[i] = clique_max[i];
	}
	
	for(i=0;i < taille - 1; i++)
	{
		if(tab[i] == 1)
		{
			for(j=i+1;j < taille ; j++)
			{
				if(tab[j]==1 && (couple_atome[i].a1 == couple_atome[j].a1))
						tab[j] = 0;
			}
		}
	}

	for(i=0;i < taille ; i++)
	{
		if(tab[i] ==1)
			nb_at++;
	}

	for(i=0;i < taille - 1; i++)
	{
		if(tab[i] == 1)
		{
			for(j = i+1;j < taille; j++)
			{
				if(tab[j] == 1)
				{
					
					i1 = couple_atome[i].a1;
					j1 = couple_atome[j].a1;
					if(M[pos1].matrice_liaisons[i1][j1] != AUCUNE_LIAISON)
						nb_liaisons ++;
				}
			}
		}
	}
	free(couple_atome);
	taille_graphe_commun[0] = nb_at;
	taille_graphe_commun[1] = nb_liaisons;
	return taille_graphe_commun;
	
}




struct molecule * lecture_fichier_chebi()
{//lecture du fichier chebi.sdf
	
	FILE *F;
	F = fopen("ChEBI_lite.sdf","r");
	
	if ( F == NULL ) 
	{
		 fprintf(stderr,"Cannot open ChEBI_lite.pdf file\n"); 
		 exit(1); 
	}
	init_atom_num();
	int nb_mol, DEB = 0, FIN = NB_MOLECULES;
	struct molecule *M = malloc(NB_MOLECULES*sizeof(struct molecule));
	
	if (M == NULL)
	{
		fprintf(stderr,"Not enough memory for M\n"); 
		exit(3); 
	}
	struct molecule m;
	printf("1. Lecture des molecules : %.3lf s\n",chrono());

	for(nb_mol = DEB ; nb_mol < FIN ; nb_mol++)
	{
		if (nb_mol % 1000 == 0) 
		{ 
			fprintf(stdout,"\r%5d / %d",nb_mol,FIN);
			fflush(stdout); 
		}
		m = lire_molecule(F);
		M[nb_mol] = m;
	}
	
	fclose(F);
	fprintf(stdout,"\r%5d / %d\n",nb_mol,FIN); 
	printf("Fin de la Lecture des molecules : %.3lf s\n",chrono());
	return M;
	
}

float mesure_similarite (int g1_chebi,int g2_chebi,struct molecule *M,double date,int taille_limite)
{//calcul du degré de similarité
	
	float similarite= -14;
	int pos1,pos2;
	pos1 = position_M(g1_chebi,M);
	pos2 = position_M(g2_chebi,M);
	
	graph g12 = graphe_produit(g1_chebi,g2_chebi,M);
  int taille = nbnodes(g12);
  int** liaisons = build_matrix_from_graph(g12);	
  
	if( taille_limite != 0 && ( taille > taille_limite))
	{
		similarite = -2;
	}
	else
	{
    
		int* clique = clique_max(g12, (long)date);

		int* taille_graphe_commun = graphe_g12(g12,clique,M,g1_chebi,g2_chebi);
    free(clique);
	  int nb_atomes_communs = taille_graphe_commun[0];	
	  int nb_laisons_communs = taille_graphe_commun[1];	
		
		if(date == 0 || (chrono() - last_chrono <= date))
		{
			float num = (float)((nb_atomes_communs + nb_laisons_communs)*(nb_atomes_communs + nb_laisons_communs));
			float denum = (float)((M[pos1].nb_atomes + M[pos1].nb_liaisons)*(M[pos2].nb_atomes + M[pos2].nb_liaisons));
			similarite = num/denum;
		}
		else
		{
			similarite = -1;
		}

    free(taille_graphe_commun);
	}
	
	if(liaisons != NULL)
	{
		int i;
		for (i = 0; i < taille;i++)
			free(liaisons[i]);
	}
	free(liaisons);
  destroy(g12);
  
	return similarite;
}

void similarite_all(int g1_chebi,struct molecule *M,double date,int taille_limite)
{
	int i;
	char nom[64];
	char nom2[64];
	sprintf(nom2,"%d_%d_%d",g1_chebi,(int)date,taille_limite);
	strcpy(nom, "resultats/similarite_");
	strcat(nom,nom2);
	strcat(nom,"_all.data");
	FILE *F;
	F = fopen( nom, "w");
	if(F == NULL){
		fprintf(stdout,"Impossible de creer le fichier de resultat\n");
		exit(45);
	}
	
	float r;
	
	for ( i = 0;  i < NB_MOLECULES;  i++)
	{
		last_chrono = chrono();
		if (i % 1 == 0) 
		{ 
			fprintf(stdout,"\r%5d / %d (%3d atomes) %.3lf ",i,NB_MOLECULES,M[i].nb_atomes,last_chrono);
			fflush(stdout); 
		}
		r = mesure_similarite( g1_chebi, M[i].chebi_id,M, date, taille_limite);
		fprintf(F, "%d %f\n",M[i].chebi_id,r);
		fflush(F);
	}
	
	fclose(F);
}
