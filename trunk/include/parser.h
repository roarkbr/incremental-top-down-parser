#ifndef PARSHDR
#define PARSHDR

#include "pstruct.h"

void train_parser(TDTrPtr TC);
void run_parser(char *ifile, TDParPtr PC);
void train_percep(char *ifile, TDParPtr PC, float step);
void show_copytree(FILE *fp, TreePtr otree, int tabtree, int both);

#endif /* PARSHDR */
