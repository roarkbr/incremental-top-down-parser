// pstruct.h
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

#ifndef PSTRHDR
#define PSTRHDR

#define MAXLABLEN 8092
#define WHASHM  2015179
#define FHASHM  12195257
#define SCHEMM  5800079
#define LEVELM  2750159
#define STATEM  73
#define LABELM  113
#define INITWD  100000       /* initial Wordlist allocation */
#define INITNT  200          /* initial NT allocation */
#define BINS    20           /* smoothing bins for parameter tying */
#define XVALS   10           /* cross validation for smoothing */
#define MAXPOS  8
#define MAXVALS 8000
#define MNXVALS 2000

struct WordList;
typedef struct WordList *WLPtr;

struct SLCPairs;
typedef struct SLCPairs *SLCPtr;

struct SLCPairs
{
  WLPtr LChild;
  WLPtr NewCat;
  int orig;
  SLCPtr next;
};

struct WordList 
{ 
  char *label;
  char *hlabel;
  SLCPtr slc;
  int windex;
  int children;        /* 0 STOP or T              1 NT                */
  int terminal;        /* 0 STOP or non-PT NT      1 if PT or terminal */
  int sym;
  int cclass;
  int ccount;
  int caps;
  int root;
  int oov;
  float tot;
  float *NTs;
  WLPtr unkcat;
};

struct StringList;
typedef struct StringList *SLPtr;

struct StringList
{
  WLPtr word;
  char *outlabel;
  char *leftsym;
  char *rightsym;
  SLPtr next;
};

struct Tree;
typedef struct Tree *TreePtr;

struct GoldList;
typedef struct GoldList *GLPtr;

struct Tree
{
  WLPtr NT;
  SLPtr LookWord;     /* look-ahead word from StringList */
  int LookPT;         /* look-ahead PT number from StringList entry */
  float Score;        /* parse score */
  float WScore;       /* POS -> word scores only */
  float PScore;       /* percep score */
  float *LookScore;   /* look-ahead scores */
  SLPtr LookUpd;      /* look-ahead for scores */
  TreePtr LSib;       /* left-sibling */
  TreePtr Par;        /* parent */
  TreePtr HChild;     /* head constituent of parent constituent */
  TreePtr CCnode;     /* left c-commanding node in tree */
  GLPtr gold;         /* pointer to gold standard derivation, if there */
  int *featarr;       /* feature array */
  int points;         /* number of trees pointing to this node */
  TreePtr prev;
};

struct GoldList
{
  WLPtr move;
  TreePtr best;
  TreePtr gtree;
  GLPtr next;
};

struct TreeHeap;
typedef struct TreeHeap *TreeHPtr;

struct TreeHeap
{
  TreePtr tree;
  int under;
  int preterm;
  float score;
  int rank;
  TreeHPtr left;
  TreeHPtr right;
  TreeHPtr next;
};

struct WordHash;
typedef struct WordHash *WHashPtr;

struct WordHash
{
  WLPtr WD;
  WHashPtr next;
};

typedef struct Lexicon
{ 
  WLPtr *WDs;
  WHashPtr *WHash;
  int numWDs;
  int numWH;
  int root;
  int allocWDs;
  int numNTs;
  int numPTs;
  float *nts;
  WLPtr allword;
} *LexPtr;

struct TreeWalk;
typedef struct TreeWalk *TWPtr;

struct TreeWalk
{
  int findex;
  int numoves;
  char *moves;
  int numqueries;
  WLPtr *WQueries;
  int SQuery;
  TWPtr suffmove;
};

struct TWSchema;
typedef struct TWSchema *TWSchPtr;

struct TWSchema
{
  int numfuncts;
  TWPtr *functs;
};

struct LambdaStruct;
typedef struct LambdaStruct *LamPtr;

struct LambData;
typedef struct LambData *LamDaPtr;

struct FeatStruct;
typedef struct FeatStruct *FeatPtr;

struct FeatData;
typedef struct FeatData *FDataPtr;

struct FeatHash;
typedef struct FeatHash *FHashPtr;

struct CFeatData;
typedef struct CFeatData *CFDataPtr;

struct LambData
{
  FeatPtr feat;
  CFDataPtr cfeat;
  float RootExp;
  float AddExp;
  float phat;
  float psmoo;
  LamDaPtr next;
  LamDaPtr anext;
  LamDaPtr backoff;
};

struct LambdaStruct
{
  float *olambdas;
  float *lambdas;
  int batch;
  LamDaPtr ldata;
  LamDaPtr alldata;
};

struct FDataLam;
typedef struct FDataLam *FLamPtr;

struct FDataLam
{
  int bucket;
  int ind;
  FDataPtr next;
};

struct FeatData
{
  float score;
  float bscore;
  FeatPtr prefix;
  FLamPtr flam;
  int numc;
  CFDataPtr *cfeat;
};

struct FeatStruct
{
  int state;
  int label;
  int level;
  int schema;
  int dstate;
  FDataPtr data;
};

struct CDataLam;
typedef struct CDataLam *CLamPtr;

struct CDataLam
{
  int sbatch;
  LamDaPtr lambda;
  float *pscore;
};

struct CFeatData
{
  int label;
  int idx;
};

struct FeatHash
{
  FeatPtr feat;
  FHashPtr next;
};

typedef struct FltScore
{
  float *rsc;
} *FScorePtr;

struct BORev;
typedef struct BORev *BORPtr;

struct BORev
{
  int state;
  BORPtr next;
};

struct FeatUpd;
typedef struct FeatUpd *FUDPtr;

struct FeatUpd
{
  int idx;
  int lastidx;
  float bcost;
  float count;
  FUDPtr next;
};

typedef struct LearnScore
{
  FScorePtr *fscores;
  float *scores;
  float *ascores;
  float *pscores;
  int *lastupd;
  float *ppsc;
  float *apsc;
  int *lastpsc;
  int iter;
  int xvals;
  int bsize;
  BORPtr *backs;
} *LearnPtr;

typedef struct ScoreStruct
{
  int numind;
  int numalloc;
  CLamPtr *clams;
  float *scores;
  float *pscores;
  float *ppsc;
  int numpsc;
  int alcpsc;
  LearnPtr lrn;
} *ScorePtr;

typedef struct Features
{
  int numschema;
  TWSchPtr *schema;
  int numfeats;
  TWPtr *featfuncts;
  int fstates;
  FHashPtr *FHash;
  int numFH;
  ScorePtr scores;
} *FtrPtr;

typedef struct ParserStats
{ 
  int owds;
  int cwds;
  int lupd;
  float surprisal;
  float nonlex;
  float lex;
  float ambig;
  float open;
  float rernk;
  float toprr;
  float dsteps;
  float cl_surprisal;
  float cl_nonlex;
  float cl_lex;
  float cl_ambig;
  float cl_open;
  float cl_rernk;
  float cl_toprr;
  float cl_dsteps;
} *PStats;

typedef struct Parser 
{ 
  LexPtr lexicon;
  FtrPtr features;
  TreePtr ntree;
  TreePtr stree;
  float norm;
  float snorm;
  float cond;
  float lcond;
  float step;
  float *exps;
  float *pexps;
  int *fd;
  int sch;
  int wsch;
  char *hlabel;
  PStats pstat;
} *ParPtr; 

typedef struct TDTrainer
{ 
  ParPtr parser;
  int verbose;
  float thresh;
  int slcprods;
  int slcnts;
  int xvals;
  int xlast;
  char *tfile;
  char *hfile;
  char *ofile;
} *TDTrPtr; 

typedef struct TDParser
{ 
  ParPtr parser;
  int verbose;
  int nbest;
  int tabtree;
  int prefix;
  int partial;
  int allpossible;
  int maxpos;
  float thresh;
  double *wscores;
  double *sscores;
  char *pfile;
  char *sfile;
  char *ofile;
} *TDParPtr; 

#endif /* PSTRHDR */
