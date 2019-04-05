CFLAGS=-g -Wall
#CFLAGS=-O2 -Wall

CC = gcc -Wall -Wextra
CXX = g++ -Wall -Wextra

all:distribution_taille_graphe_cycles

run: distribution_taille_graphe_cycles
	./distribution_taille_graphe_cycles 165 598 20 1600

val: distribution_taille_graphe_cycles
	valgrind --leak-check=full --show-leak-kinds=all ./distribution_taille_graphe_cycles 165 598 20 1600	

distribution_taille_graphe_cycles: distribution_taille_graphe_cycles.o lecture_molecule_sdf.o fonctions_mces.o helpers/graph.o helpers/cliquerecursif.o
	gcc ${CFLAGS} -o $@ $^

distribution_taille_graphe_cycles_scip : distribution_taille_graphe_cycles.o lecture_molecule_sdf.o fonctions_mces.o helpers/graph.o helpers/proglin_helper_scip.o helpers/sciplib.a helpers/cliquescip.o 
	$(CXX) -I helpers/scip -o $@ $^ -lpopt -lgmp -lm -lz -lreadline -lncurses 	

distribution_taille_graphe_cycles.o: distribution_taille_graphe_cycles.c
	gcc ${CFLAGS} -c distribution_taille_graphe_cycles.c

fonctions_mces.o: fonctions_mces.c fonctions_mces.h 
	gcc ${CFLAGS} -c fonctions_mces.c
	
lecture_molecule_sdf.o: lecture_molecule_sdf.c lecture_molecule_sdf.h
	gcc ${CFLAGS} -c lecture_molecule_sdf.c

helpers/proglin_helper_scip.o : helpers/proglin_helper_scip.c helpers/proglin_helper.h helpers/sciplib.a
	$(CC) -I helpers/scip -o $@ -c $<


	
clean: 
	rm -f distribution_taille_graphe_cycles
	rm -f *.o
	rm -f helpers/*.o

