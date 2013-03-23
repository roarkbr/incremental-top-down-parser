// parser.h
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

#ifndef PARSHDR
#define PARSHDR

#include "pstruct.h"

void train_parser(TDTrPtr TC);
void run_parser(char *ifile, TDParPtr PC);
void train_percep(char *ifile, TDParPtr PC, float step);
void show_copytree(FILE *fp, TreePtr otree, int tabtree, int both);

#endif /* PARSHDR */
