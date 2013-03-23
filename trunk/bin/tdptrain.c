// tdptrain.c
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pstruct.h"
#include "putil.h"
#include "parser.h"

#define USAGE "Usage: %s [-opts] treebank                                  \n\
  Note that the treebank cannot be piped, as it may be required repeatedly.\n\
                                                                           \n\
Options:                                                                   \n\
 -F file         output file                                               \n\
 -p              selective LC transform for all left-recursive productions \n\
 -n              selective LC transform for all non-POS left-children      \n\
 -s \"(X(Y\"       selective LC transform for 'X -> Y a' productions       \n\
 -l file         selective LC file                                         \n\
 -f \"a:b:c...\"   feature set {a,ab,abc,...} for the models               \n\
 -m file         feature set file                                          \n\
 -t thresh       threshold on min -log backoff probability                 \n\
 -x x            number of cross-validation sets for EM backoff (def:10)   \n\
 -v              verbose                                                   \n\
 -?              info/options                                              \n"

int main(int ac, char *av[])
{
  int c, err=0;
  extern char *optarg;
  extern int optind;
  TDTrPtr TrainConf = TDTrainDefault();

  while ((c = getopt(ac, av, "F:s:f:m:t:l:x:npv?")) != -1)
    switch (c) {
    case 'v':
      TrainConf->verbose = 1;
      break;
    case 'p':
      TrainConf->slcprods = 1;
      break;
    case 'n':
      TrainConf->slcnts = 1;
      break;
    case 'x':
      TrainConf->xvals = atof(optarg);
      break;
    case 't':
      TrainConf->thresh = atof(optarg);
      break;
    case 's':
      add_slcprods(TrainConf,optarg);
      break;
    case 'l':
      add_slcfile(TrainConf,optarg);
      break;
    case 'f':
      add_featschemas(TrainConf,optarg);
      break;
    case 'm':
      add_featfile(TrainConf,optarg);
      break;
    case 'F':
      TrainConf->ofile = optarg;
      break;
    case '?':
    default:
      err++;
    }

  if (err || ac != optind+1) {
    fprintf(stderr, USAGE, av[0]);
    exit(1);
  }
  TrainConf->tfile = av[optind];
  train_parser(TrainConf);
}
