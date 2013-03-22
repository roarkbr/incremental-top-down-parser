#ifndef PIOHDR
#define PIOHDR

#include "pstruct.h"

FILE *gopenfile(char *file, char *mode);
void write_parser(ParPtr par, char *ofile);
ParPtr read_parser(FILE *fp, int lrn);
void read_norms(char *readtree, TDTrPtr TC, LexPtr oldlex);
void read_terms(char *readtree, TDTrPtr TC);
SLPtr get_string(char *readstr, TDParPtr PC);
SLPtr get_postring(char *readstr, TDParPtr PC);
int quotes(char *readterm);

#endif /* PIOHDR */
