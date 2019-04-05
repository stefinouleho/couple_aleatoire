#include "proglin_helper.h"
#include "clique.h"
#include "graph.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>


long chrono_clique() 
{
	struct timeval tv;
	long date_courante;
	gettimeofday(&tv,NULL);
	date_courante = tv.tv_sec;
	return date_courante;
}

void calcul_cl(graph g, int* dans_clique_max, int* taille_clique_max, int *dans_clique,int taille_clique,int *candidat,int taille_candidat, long maxdate)
{	//calcul de la clique max recursif
	
  int taille = nbnodes(g);
  
	if (maxdate != 0)
	{	
		if(chrono_clique()  > maxdate) 
			return;
	}
	int i,j;
	if( taille_candidat == 0)
	{
		if(taille_clique > *taille_clique_max){
			*taille_clique_max = taille_clique;
			for (i = 0 ;  i < taille ; i ++)
				dans_clique_max[i] = dans_clique[i];
		}
		return;
	}
	
	if (taille_candidat + taille_clique <= *taille_clique_max)
	{

		return;
	}
	
	//else 
	int taille_candidat_temp;
	int *candidat_temp;
	candidat_temp = malloc( taille * sizeof(int));



	for (i = 0 ;  i < taille ; i ++)
	{
		if ( candidat[i] == 1)
		{
			candidat[i] = 0;
			dans_clique[i] = 1 ;
			taille_candidat_temp = taille_candidat;
			
			for (j = 0 ;  j < taille ; j ++)
			{
				candidat_temp[j] = candidat[j]; 
				if ((candidat[j] == 1) && (are_neighbors(g, i, j) == 0))
				{
					candidat_temp[j] = 0;
					taille_candidat_temp--;	
				}	
			}
			
			taille_candidat_temp --;
			calcul_cl(g,dans_clique_max,taille_clique_max,dans_clique,taille_clique + 1,candidat_temp,taille_candidat_temp,maxdate);
			dans_clique[i] = 0;
			
		}	
		
	}
	free(candidat_temp);
		 
}
int* clique_max(graph g, long timeout)
{	//Debut calcul de la clique -- Initialisation
	int i;
	int *candidat;
	int *dans_clique;
  int taille_clique_max;
	int *dans_clique_max;

  int taille = nbnodes(g);
	
  taille_clique_max = 0;
  dans_clique_max = malloc( taille *sizeof(int));
	if (!dans_clique_max) { fprintf(stderr,"cannot malloc dans_clique_max %d\n",taille); exit(41); }
  dans_clique = malloc( taille *sizeof(int));
	if (!dans_clique) { fprintf(stderr,"cannot malloc dans_clique %d\n",taille); exit(42); }
	candidat = malloc( taille *sizeof(int));
	if (!candidat) { fprintf(stderr,"cannot malloc candidat %d\n",taille); exit(43); }
	
	//initialisation 
	for(i = 0; i < taille ; i++ )
	{
		candidat[i] 	= 1;
		dans_clique[i]	= 0;
		dans_clique_max[i] = 0;
	}
	
	
	calcul_cl(g,dans_clique_max,&taille_clique_max,dans_clique,0,candidat,taille,chrono_clique() + timeout); // 0 taille de la clique initial  et m.nb_atome = nb sommets candidats

   
	free(dans_clique);
	free(candidat);
  return dans_clique_max;
}
