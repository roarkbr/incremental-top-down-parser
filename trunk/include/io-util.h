// io-util.h 
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Copyright 2004-2013 Brian Roark
// Authors: roarkbr@gmail.com  (Brian Roark)

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
