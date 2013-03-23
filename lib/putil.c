// putil.c
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
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "pstruct.h"
#include "putil.h"
#include "io-util.h"
#include "parser.h"

void failproc(char *msg)
{
  if (msg != NULL) fprintf(stderr,"%s\n",msg);
  exit(1);
}

void verbose_count(int verbose, int k, int mod)
{
  if (verbose && k%mod == 0) fprintf(stderr," %d ",k);
  else if (verbose && 10*k%mod == 0) fprintf(stderr,".");
}

void wordtoken(char *readtree, char *readterm, int *i, int len, int STOPs)
{
  int j=0;
  while (i[0] < len && ws(readtree[i[0]])) i[0]++;
  while (i[0] < len && !ws(readtree[i[0]])) {
    if (readtree[i[0]] == ')') {
      if (j > 0 && STOPs) break;
      else if (STOPs) readterm[j++]=readtree[i[0]++];
      else i[0]++;
    }
    else readterm[j++]=readtree[i[0]++];
  }
  readterm[j]=0;
}

void make_stop(LexPtr lex)
{
  int hash;
  WLPtr word;
  char STOP[2];
  STOP[0]=')';
  STOP[1]=0;
  word=find_word(&STOP[0],lex,&hash);
  if (word == NULL) word=insert_word(&STOP[0],lex,hash,0,0);
}

void slcentry(WLPtr pword, WLPtr lword, WLPtr word, int orig)
{
  SLCPtr slc=malloc(sizeof(struct SLCPairs));
  slc->LChild=lword;
  slc->NewCat=word;
  slc->orig=orig;
  slc->next=pword->slc;
  pword->slc=slc;
}

void make_slcentry(WLPtr pword, WLPtr lword, WLPtr word)
{
  slcentry(pword,lword,word,1);
  slcentry(word,lword,pword,0);
}

void add_slcprods(TDTrPtr TC,char *label)
{
  int i=0, j, len=strlen(label), hash;
  LexPtr lex=TC->parser->lexicon;
  WLPtr word=NULL, pword, lword;
  char str[MAXLABLEN];
  j=0;
  str[j++]=label[i++];
  while (i < len && label[i] != '(') str[j++]=label[i++];
  str[j]=0;
  pword=find_word(&str[0],lex,&hash);
  if (pword == NULL) pword=insert_word(&str[0],lex,hash,1,0);
  j=0;
  str[j++]=label[i];
  label[i++] = '/';
  while (i < len) str[j++]=label[i++];
  str[j]=0;
  lword=find_word(&str[0],lex,&hash);
  if (lword == NULL) lword=insert_word(&str[0],lex,hash,1,0);
  word=find_word(label,lex,&hash);
  if (word == NULL) {
    word=insert_word(label,lex,hash,1,0);
    make_slcentry(pword,lword,word);
  }
}

TWSchPtr init_fschema()
{
  TWSchPtr fschema;
  fschema=malloc(sizeof(struct TWSchema));
  fschema->numfuncts=0;
  fschema->functs=NULL;
  return fschema;
}

TWPtr init_twalk(int numoves, int numqs, int SQuery)
{
  TWPtr twalk=malloc(sizeof(struct TreeWalk));
  twalk->numoves=numoves;
  twalk->moves=malloc(numoves*sizeof(char));
  twalk->numqueries=numqs;
  twalk->WQueries=malloc(numqs*sizeof(WLPtr));
  twalk->SQuery=SQuery;
  twalk->suffmove=NULL;
  twalk->findex=-1;
  return twalk;
}

int twop(char c)
{
  if (c == 'u' || c == 'U') return 1;
  if (c == 'l' || c == 'L') return 1;
  if (c == 'f' || c == 'F') return 1;
  if (c == 'r' || c == 'R') return 1;
  if (c == 'd' || c == 'D') return 1;
  if (c == 'c' || c == 'C') return 1;
  if (c == 'h' || c == 'H') return 1;
  if (c == '*') return 1;
  if (c == 's') return 1;
  return 0;
}

int valop(char c)
{
  if (c == 'N') return 1;
  if (c == 'W') return 1;
  if (c == 'P') return 1;
  if (c == 'T') return 1;
  return 0;
}

int getvalop(TreePtr tree, char c)
{
  SLPtr LW;
  if (tree == NULL) return -1;
  if (c == 'N') return tree->NT->windex;
  LW = tree->LookWord;
  if (LW == NULL) return -1;
  if (c == 'T' || (c == 'W' && tree->NT->terminal)) return LW->word->windex;
  if (c == 'P') return tree->LookPT;
  return -1; 
}

TreePtr askquery(TreePtr ftree, int s, int q, WLPtr *Q, TreePtr tree)
{
  if (s==q && Q[s]==NULL) return NULL;
  if (ftree == NULL && Q[s] == NULL) return tree;
  if (ftree == NULL) return ftree;
  if (ftree->NT == Q[s]) return ftree;
  return NULL;
}

TreePtr getstar(TreePtr ftree, TreePtr tree, char lc)
{
  TreePtr rtree = gettwop(ftree, tree, lc, lc), rftree=NULL;
  while (rtree != NULL && rtree != rftree) {
    rftree = rtree;
    rtree = gettwop(rtree, tree, lc, lc);
  }
  return rftree;
}

TreePtr getright(TreePtr ftree, TreePtr tree, char c)
{ 
  TreePtr rtree=NULL, stree=tree;
  while (stree != NULL && stree->Par != ftree->Par) stree = stree->prev;
  if (stree == NULL || stree == ftree) return rtree;
  rtree = stree;
  while (stree != NULL && stree->LSib != ftree) {
    stree = stree->LSib;
    if (stree != NULL && !stree->NT->sym) rtree=stree;
  }
  if (c == 'r') return stree;
  return rtree;
}

TreePtr getdown(TreePtr ftree, TreePtr tree, char c)
{
  TreePtr stree = tree;
  if (ftree->NT->terminal) return NULL;
  while (stree != NULL && stree->Par != ftree) stree = stree->prev;
  if (stree == NULL || c == 'd') return stree;
  while (stree != NULL && stree->NT->sym) stree = stree->LSib;
  return stree;
}

TreePtr gethchild(TreePtr ftree, TreePtr tree, char c)
{
  TreePtr stree = tree, rftree=NULL;
  if (!ftree->NT->children) return NULL;
  if (ftree->NT->terminal) return ftree;
  while (stree != NULL && stree->Par != ftree) stree = stree->prev;
  rftree=stree->HChild;
  if (rftree == NULL) rftree=stree->LSib;
  if (c == 'h') return rftree;
  while (rftree != NULL && rftree->NT->sym) rftree = rftree->LSib;
  return rftree;    
}

TreePtr nosym(TreePtr tree)
{
  while (tree != NULL && tree->NT->sym) tree = tree->LSib;
  return tree;    
}

TreePtr oldleft(TreePtr ftree, int nsym)
{
  TreePtr rftree=ftree->LSib;
  if (nsym) rftree=nosym(rftree);
  if (rftree == NULL && ftree->Par != NULL && ftree->Par->NT->slc != NULL &&
      !ftree->Par->NT->slc->orig) rftree=ftree->Par->Par;
  return rftree;
}

TreePtr gettwop(TreePtr ftree, TreePtr tree, char c, char lc)
{
  TreePtr rftree=NULL, rtree, stree;
  if (ftree == NULL) return ftree;
  switch (c) {
  case '*': return getstar(ftree,tree,lc);
  case 'u': return ftree->Par;
  case 'U': return ftree->Par;
  case 'l': return ftree->LSib;
  case 'L': return nosym(ftree->LSib);
  case 'f': return oldleft(ftree,0);
  case 'F': return oldleft(ftree,1);
  case 'r': return getright(ftree,tree,c);
  case 'R': return getright(ftree,tree,c);
  case 'd': return getdown(ftree,tree,c);
  case 'D': return getdown(ftree,tree,c);
  case 'c': return ftree->CCnode;
  case 'C': return nosym(ftree->CCnode);
  case 'h': return gethchild(ftree,tree,c);
  case 'H': return gethchild(ftree,tree,c);
  default:
    fprintf(stderr,"code %c not found\n",c);
    exit(1);
  }
}

int trimoves(char *moves, char *nmoves, int numoves, int *newnumoves)
{
  int i, last=numoves-2, plast=numoves-3;
  char cl;
  if (last < 0) return 0;
  cl=moves[last];
  if (cl == 'q') return 0;
  if (!twop(cl)) return 0;
  if (cl == '*') {
    newnumoves[0]=numoves-2;
    for (i=0; i < last-1; i++) nmoves[i]=moves[i];
    nmoves[last-1]=moves[last+1];
    return 1;
  }
  if ((plast >= 0 && moves[plast] == cl) || 
      (cl != 'l' && cl != 'L' && cl != 'R' && cl != 'r')) {
    newnumoves[0]=numoves-1;
    for (i=0; i < last; i++) nmoves[i]=moves[i];
    nmoves[last]=moves[last+1];
    return 1;
  }
  if (cl != 'l' || cl != 'L' || cl != 'R' || cl != 'r') {
    newnumoves[0]=numoves;
    for (i=0; i <= last+1; i++) nmoves[i]=moves[i];
    nmoves[last]='u';
    return 1;
  }
  return 0;
}

TWPtr find_twalk(ParPtr parser,char *moves,int numoves, WLPtr *queries,int numqs,int SQuery)
{
  int i, j;
  TWPtr tfunct;
  for (i=0; i < parser->features->numfeats; i++) {
    tfunct=parser->features->featfuncts[i];
    if (tfunct->numoves != numoves || tfunct->numqueries != numqs) continue;
    for (j=0; j < numoves; j++) if (moves[j] != tfunct->moves[j]) break;
    if (j < numoves) continue;
    for (j=0; j < numqs; j++) if (queries[j] != tfunct->WQueries[j]) break;
    if (j < numqs) continue;
    if (tfunct->SQuery != SQuery) continue;
    return tfunct;
  }
  parser->features->numfeats++;
  if (parser->features->numfeats == 1)
    parser->features->featfuncts=malloc(sizeof(TWPtr));
  else parser->features->featfuncts=realloc(parser->features->featfuncts,
				  parser->features->numfeats*sizeof(TWPtr));
  tfunct=parser->features->featfuncts[parser->features->numfeats-1]=
    init_twalk(numoves,numqs,SQuery);
  for (j=0; j < numoves; j++) tfunct->moves[j]=moves[j];
  for (j=0; j < numqs; j++) tfunct->WQueries[j]=queries[j];
  return tfunct;
}

void update_twalk(ParPtr parser,TWPtr *functs, TWSchPtr fschema, char *moves, int numoves, WLPtr *queries, int numqs, int SQuery, int incl)
{
  int i, newnumoves;
  TWPtr funct;
  char newmoves[1000];
  if (incl) {
    if (trimoves(moves,&newmoves[0],numoves,&newnumoves))
      update_twalk(parser,functs,fschema,&newmoves[0],newnumoves,
		   queries,numqs,SQuery,incl);
  }
  funct=find_twalk(parser,moves,numoves,queries,numqs,SQuery);
  for (i=0; i < fschema->numfuncts; i++)
    if (functs[i] == funct) break;
  if (i==fschema->numfuncts) functs[fschema->numfuncts++]=funct;
}

void update_fschema(ParPtr parser,TWPtr *functs, TWSchPtr fschema)
{
  int i;
  fschema->functs=malloc(fschema->numfuncts*sizeof(TWPtr));
  for (i=0; i < fschema->numfuncts; i++)
    fschema->functs[i]=functs[i];
}

int getqueries(char *label, char c, int i, int len, char *moves, int *numov, WLPtr *queries, int *numqs, int *SQuery, LexPtr lex)
{
  int j=0, hash;
  char readtoken[MAXLABLEN];
  if (c=='q') {
    if (label[i] != '\'') readtoken[j++]=0;
    else {
      i++;
      while (i < len && label[i] != '\'') readtoken[j++]=label[i++];
      if (i >= len) {
	sprintf(&moves[0],"Mis-parse of label query: %s\n",label);
	failproc(&moves[0]);
      }
      i++;
    }
  }
  else if (c == '0') moves[numov[0]++]='l';
  else if (c == 'k') {
    j=3;
    sprintf(&readtoken[0],"(CC");
  }
  if (j == 0) queries[numqs[0]++]=NULL;
  else if (j==1 && readtoken[0]==0) {
    SQuery[0]=numqs[0];
    queries[numqs[0]++]=NULL;
  }
  else {
    readtoken[j]=0;
    queries[numqs[0]]=find_word(&readtoken[0],lex,&hash);
    if (queries[numqs[0]++] == NULL) 
      queries[numqs[0]-1]=insert_word(&readtoken[0],lex,hash,1,0);
  }
  moves[numov[0]++]='q';
  return i;
}

void parse_fschema(char *label, ParPtr parser, TWSchPtr fschema, TWPtr *functs)
{
  int i=0, len=strlen(label), numov, numqs, incl, hash, SQuery;
  char c, moves[MAXLABLEN];
  LexPtr lex=parser->lexicon;
  WLPtr queries[1000],word;
  while (i<len) {
    while (i<len && label[i] == ':') i++;
    if (i >= len) return;
    numov=numqs=incl=0;
    SQuery=-1;
    while (i<len && label[i] != ':') {
      c=label[i++];
      if (c == 'i') incl=1;
      else if (twop(c)) moves[numov++]=c;
      else if (valop(c)) {
	if (i < len && label[i] != ':') {
	  sprintf(&moves[0],"Value returning operation %c must be last: %s\n",
		  c,label);
	  failproc(&moves[0]);
	}
	moves[numov++]=c;
	update_twalk(parser,functs,fschema,&moves[0],numov,
		     &queries[0],numqs,SQuery,incl);
      }
      else if (c == 'j') {
	moves[numov++]='d';
	moves[numov++]='l';
	moves[numov++]='*';
      }
      else if (c=='q' || c == 'k' || c == '0') 
	i=getqueries(label,c,i,len,&moves[0],&numov,&queries[0],&numqs,&SQuery,lex);
      else {
	sprintf(&moves[0],"Unrecognized move %c: %s\n",c,label);
	failproc(&moves[0]);    
      }
    }
  }
}

void add_featschemas(TDTrPtr TC,char *label)
{
  ParPtr parser=TC->parser;
  FtrPtr feats=parser->features;
  TWPtr funct, functs[1000];
  TWSchPtr fschema;
  feats->numschema++;
  if (feats->schema == NULL) feats->schema=malloc(sizeof(TWSchPtr));
  else feats->schema=realloc(feats->schema,feats->numschema*sizeof(TWSchPtr));
  fschema=feats->schema[feats->numschema-1]=init_fschema();
  parse_fschema(label,parser,fschema,&functs[0]);
  update_fschema(parser,&functs[0],fschema);
}

void add_featfile(TDTrPtr TC, char *filename)
{
  int i;
  FILE *fp=fopen(filename,"r");
  char readtree[MAXLABLEN];
  while (fgets(readtree,MAXLABLEN,fp) != NULL) {
    i=0;
    while (i<MAXLABLEN && readtree[i] != '\n' && readtree[i] != EOF) i++;
    readtree[i]=0;
    add_featschemas(TC,&readtree[0]);
  }
  fclose(fp);
}

void compile_feats(TDTrPtr TC, char *filename)
{
  int i, j, k;
  ParPtr parser=TC->parser;
  TWPtr tfunct, hfunct;
  for (i=0; i < parser->features->numfeats; i++) {
    tfunct=parser->features->featfuncts[i];
    tfunct->findex = i;
    if (tfunct->numqueries != 0) continue;
    for (j=0; j < parser->features->numfeats; j++) {
      hfunct=parser->features->featfuncts[j];
      if (j==i || hfunct->numqueries != 0 || 
	  hfunct->numoves != tfunct->numoves-1) continue;
      for (k=0; k < hfunct->numoves; k++)
	if (tfunct->moves[k+1] != hfunct->moves[k]) break;
      if (k < hfunct->numoves) continue;
      tfunct->suffmove = hfunct;
    }
  }
}

void add_slcfile(TDTrPtr TC, char *filename)
{
  int i;
  FILE *fp=fopen(filename,"r");
  char readtree[MAXLABLEN];
  while (fgets(readtree,MAXLABLEN,fp) != NULL) {
    i=0;
    while (i<MAXLABLEN && readtree[i] != '\n' && readtree[i] != EOF) i++;
    readtree[i]=0;
    add_slcprods(TC,&readtree[0]);
  }
  fclose(fp);
}

LexPtr LexDefault(int numWH, int numWD, int stop)
{
  int i;
  LexPtr lexicon=malloc(sizeof(struct Lexicon));
  lexicon->numWDs=lexicon->numNTs=lexicon->numPTs=0;
  lexicon->allword=NULL;
  lexicon->nts=NULL;
  lexicon->allocWDs=numWD;
  lexicon->WDs=malloc(lexicon->allocWDs*sizeof(WLPtr));
  for (i=0; i < lexicon->allocWDs; i++) lexicon->WDs[i]=NULL;
  lexicon->numWH=numWH;
  lexicon->WHash=malloc(numWH*sizeof(WHashPtr));
  for (i=0; i < numWH; i++) lexicon->WHash[i]=NULL;
  if (stop) make_stop(lexicon);
  return lexicon;
}

void WDFree(WLPtr word)
{
  free(word->label);
  free(word);
}

void WHashFree(WHashPtr *whash, int numWH)
{
  int i;
  WHashPtr thash;
  for (i=0; i < numWH; i++)
    while (whash[i] != NULL) {
      thash=whash[i];
      whash[i]=thash->next;
      free(thash);
    }
  free(whash);
}

void LexFree(LexPtr lex)
{
  int i;
  for (i=0; i < lex->numWDs; i++) WDFree(lex->WDs[i]);
  free(lex->WDs);
  WHashFree(lex->WHash,lex->numWH);
  free(lex);
}

void FeatFree(FHashPtr fhash)
{
  FHashPtr thash;
  while (fhash != NULL) {
    thash=fhash;
    fhash=fhash->next;
    free(thash->feat->data);
    free(thash->feat);
    free(thash);
  }
}

void FHashFree(FtrPtr features)
{
  int i;
  for (i=0; i < features->numFH; i++) FeatFree(features->FHash[i]); 
  free(features->FHash);
}

ScorePtr ScrDefault()
{
  ScorePtr scr=malloc(sizeof(struct ScoreStruct));
  scr->numind=scr->numalloc=scr->numpsc=scr->alcpsc=0;
  scr->clams=NULL;
  scr->scores=scr->pscores=scr->ppsc=NULL;
  scr->lrn=NULL;
  return scr;
}

FtrPtr FtrDefault(int numFH)
{
  int i;
  FtrPtr features=malloc(sizeof(struct Features));
  features->numschema=0;
  features->schema=NULL;
  features->numfeats=0;
  features->featfuncts=NULL;
  features->fstates=1;
  features->numFH=numFH;
  features->FHash=malloc(numFH*sizeof(FHashPtr));
  for (i=0; i < numFH; i++) features->FHash[i]=NULL;
  features->scores=ScrDefault();
  return features;
}

void reset_pstat(PStats pstat)
{
  pstat->owds=pstat->cwds=pstat->lupd=0;
  pstat->surprisal=pstat->cl_surprisal=0;
  pstat->nonlex=pstat->cl_nonlex=0;
  pstat->lex=pstat->cl_lex=0;
  pstat->ambig=pstat->cl_ambig=0;
  pstat->open=pstat->cl_open=0;
  pstat->rernk=pstat->cl_rernk=0;
  pstat->toprr=pstat->cl_toprr=0;
  pstat->dsteps=pstat->cl_dsteps=0;
}

PStats init_pstat()
{
  PStats pstat=malloc(sizeof(struct ParserStats));
  reset_pstat(pstat);
  return pstat;
}

ParPtr ParserDefault(int numWH, int numFH, int numWD, int stop)
{
  ParPtr parser=malloc(sizeof(struct Parser));
  parser->lexicon=LexDefault(numWH,numWD,stop);
  parser->features=FtrDefault(numFH);
  parser->ntree=NULL;
  parser->stree=NULL;
  parser->step=1.0;
  parser->pstat=init_pstat();
  return parser;
}

TDTrPtr TDTrainDefault()
{
  TDTrPtr TC=malloc(sizeof(struct TDTrainer));
  TC->parser=ParserDefault(WHASHM,FHASHM,INITWD,1);
  TC->verbose=0;
  TC->thresh=0;
  TC->slcprods=0;
  TC->slcnts=0;
  TC->xvals=XVALS;
  TC->xlast=1;
  TC->ofile=TC->tfile=TC->hfile=NULL;
  return TC;
}

TDParPtr TDParseDefault()
{
  TDParPtr TC=malloc(sizeof(struct TDParser));
  TC->parser=NULL;
  TC->verbose=TC->prefix=TC->partial=TC->allpossible=0;
  TC->maxpos=MAXPOS;
  TC->nbest=1;
  TC->tabtree=0;
  TC->thresh=12;
  TC->ofile=TC->pfile=TC->sfile=NULL;
  return TC;
}

int ws(char n)
{
  if (n == ' ') return 1;
  if (n == '\t') return 1;
  if (n == '\n') return 1;
  return 0;
}

int ucase(char c)
{
  if (c > 64 && c < 91) return 1;
  return 0;
}

int lcase(char c)
{
  if (c > 96 && c < 123) return 1;
  return 0;
}

int digit(char c)
{
  if (c > 47 && c < 58) return 1;
  return 0;
}

int alpha(char *str)
{
  int i, len=strlen(str);
  for (i=0; i < len; i++)
    if (ucase(str[i]) || lcase(str[i])) return 1;
  return 0;
}

int nonalpha(char *str)
{
  int i, len=strlen(str);
  for (i=0; i < len; i++)
    if (!(ucase(str[i]) || lcase(str[i]) || digit(str[i]) || str[i] == 39))
      return 1;
  return 0;
}

int mpunc(char c)
{
  if (c==39 || c==35 || c==36) return 1;
  return 0;
}

int noalpha(char *str)
{
  int i, len=strlen(str);
  if (len==2 && str[0]=='\'' && str[1]=='\'') return 1;
  for (i=0; i < len; i++)
    if (ucase(str[i]) || lcase(str[i]) || digit(str[i]) || mpunc(str[i])) return 0;
  return 1;
}

int numeral(char *str)
{
  int i, len=strlen(str);
  for (i=0; i < len; i++)
    if (ucase(str[i]) || lcase(str[i])) return 0;
  return 1;
}

int numst(char *str)
{
  int i, len=strlen(str);
  for (i=0; i < len; i++) if (!digit(str[i])) break;
  if (i != len-2) return 0;
  if (str[i] == 't' && str[i+1] == 'h') return 1;
  if (str[i] == 'r' && str[i+1] == 'd') return 1;
  if (str[i] == 'n' && str[i+1] == 'd') return 1;
  if (str[i] == 's' && str[i+1] == 't') return 1;
  return 0;
}

int suffclass(char *str)
{
  if (strcmp(str,"ble") == 0) return 1;
  if (strcmp(str,"ful") == 0) return 1;
  if (strcmp(str,"ous") == 0) return 1;
  if (strcmp(str,"ing") == 0) return 2;
  if (strcmp(&str[1],"al") == 0) return 1;
  if (strcmp(&str[1],"ly") == 0) return 3;
  if (strcmp(&str[1],"ss") == 0) return 4;
  if (strcmp(&str[1],"ll") == 0) return 5;
  if (strcmp(&str[1],"ar") == 0) return 5;
  if (strcmp(&str[1],"ry") == 0) return 5;
  if (strcmp(&str[1],"ed") == 0) return 6;
  if (str[2] == 's') return 7;
  return 0;
}

int prime(int i)
{
  int primes[300]={2,3,5,7,11,13,17,19,23,29,
                     31,37,41,43,47,53,59,61,67,71,
                     73,79,83,89,97,101,103,107,109,113,
                     127,131,137,139,149,151,157,163,167,173,
                     179,181,191,193,197,199,211,223,227,229,
                     233,239,241,251,257,263,269,271,277,281,
                     283,293,307,311,313,317,331,337,347,349,
                     353,359,367,373,379,383,389,397,401,409,
                     419,421,431,433,439,443,449,457,461,463,
                     467,479,487,491,499,503,509,521,523,541,
                     547,557,563,569,571,577,587,593,599,601,
                     607,613,617,619,631,641,643,647,653,659,
                     661,673,677,683,691,701,709,719,727,733,
                     739,743,751,757,761,769,773,787,797,809,
                     811,821,823,827,829,839,853,857,859,863,
                     877,881,883,887,907,911,919,929,937,941,
                     947,953,967,971,977,983,991,997,1009,1013,
                     1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,
                     1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,
                     1153,1163,1171,1181,1187,1193,1201,1213,1217,1223,
                     1229,1231,1237,1249,1259,1277,1279,1283,1289,1291,
                     1297,1301,1303,1307,1319,1321,1327,1361,1367,1373,
                     1381,1399,1409,1423,1427,1429,1433,1439,1447,1451,
                     1453,1459,1471,1481,1483,1487,1489,1493,1499,1511,
                     1523,1531,1543,1549,1553,1559,1567,1571,1579,1583,
                     1597,1601,1607,1609,1613,1619,1621,1627,1637,1657,
                     1663,1667,1669,1693,1697,1699,1709,1721,1723,1733,
                     1741,1747,1753,1759,1777,1783,1787,1789,1801,1811,
                     1823,1831,1847,1861,1867,1871,1873,1877,1879,1889,
                     1901,1907,1913,1931,1933,1949,1951,1973,1979,1987};
  return primes[i%300];
}

int get_whash(char *word, LexPtr lex)
{
  int i, hash, len=strlen(word);
  long lhash=0;
  for (i=0; i < len; i++) lhash += word[i]*word[i]*prime(i);
  hash=lhash%lex->numWH;
  return hash;
}

WLPtr find_word(char *word, LexPtr lex, int *hash)
{
  WHashPtr WHash;
  hash[0]=get_whash(word,lex);
  WHash=lex->WHash[hash[0]];
  while (WHash != NULL) {
    if (strcmp(word,WHash->WD->label) == 0) return WHash->WD;
    WHash=WHash->next;
  }
  return NULL;
}

WLPtr must_find_word(char *word, LexPtr lex, int *hash)
{
  char msg[MAXLABLEN];
  WHashPtr WHash;
  hash[0]=get_whash(word,lex);
  WHash=lex->WHash[hash[0]];
  while (WHash != NULL) {
    if (strcmp(word,WHash->WD->label) == 0) return WHash->WD;
    WHash=WHash->next;
  }
  sprintf(&msg[0],"%s not in lexicon!\n",word);
  failproc(&msg[0]);
}

int makelower(char *word)
{
  int i, len=strlen(word), ch=0;
  for (i=0; i < len; i++)
    if (word[i] > 64 && word[i] < 91) {
      word[i] += 32; 
      ch=1;
    }
  return ch;
}

void fill_ufeats(int *ufeats, char *llabel, char *ulabel)
{
  int i, len=strlen(llabel);
  char *suffix=NULL;
  if (lcase(ulabel[0])) ufeats[11]=1;
  if (ucase(ulabel[0])) ufeats[12]=1;
  if (llabel[0] == '-') ufeats[14]=1;
  for (i=1; i < len; i++) {
    if (!ufeats[1] && lcase(ulabel[i])) ufeats[1]=1;        /* some lower */
    if (!ufeats[2] && ucase(ulabel[i])) ufeats[2]=1;        /* some upper */
    if (!ufeats[3] && digit(llabel[i])) ufeats[3]=1;        /* some digit */
    if (!ufeats[4] && llabel[i] == '-') ufeats[4]=1;        /* some hyphen */
  }
  if (ufeats[13]) ufeats[5]=numeral(llabel);               /* CD likely */
  if (len < 3) ufeats[6]=1;                                /* short word */
  else {
    suffix=&llabel[len-3];
    if (ufeats[13] && !ufeats[5]) ufeats[7]=numst(llabel);  /* ordinal */
    ufeats[8]=suffclass(suffix);                           /* 3 letter suffix */ 
  }
}

void make_unkcat(char *word, char *nword, char *fword, int *ufeats)
{
  int i, j;
  fill_ufeats(ufeats,nword,fword);
  if (!ufeats[0] && ufeats[12]) {
    if (ufeats[8] == 7) sprintf(word,"oovNNPS");
    else sprintf(word,"oovNNP");
    return;
  }
  if (ufeats[1]+ufeats[2]+ufeats[3]+ufeats[11]+ufeats[12]+ufeats[13]==0) {
    sprintf(word,"oovSYM");
    return;
  }
  if (ufeats[5]) {
    sprintf(word,"oovCD");
    return;    
  }
  if (ufeats[4]) {
    sprintf(word,"oovHYP");
    return;
  }
  if (ufeats[7]) {
    sprintf(word,"oovORD");
    return;
  }
  if (ufeats[6]) sprintf(word,"oovSHORT");
  else sprintf(word,"oovSUF%d",ufeats[8]);
}

void find_defunk(char *word, char *outlabel, int num)
{
  if (num < 8) {
    sprintf(word,"oovSUF%d",num);
    return;
  }
  switch (num) {
  case 8:
    sprintf(word,"oovNNPS");
    return;
  case 9:
    sprintf(word,"oovNNP");
    return;
  case 10:
    sprintf(word,"oovSYM");
    return;
  case 11:
    sprintf(word,"oovCD");
    return;
  case 12:
    sprintf(word,"oovHYP");
    return;
  case 13:
    sprintf(word,"oovORD");
    return;
  case 14:
    sprintf(word,"oovSHORT");
    return;
  default:
    fprintf(stderr,"failed to find word or unk: %s\n",outlabel);
    exit(1);
  }
}

WLPtr find_nword(char *word, LexPtr lex, int *hash, WLPtr fword, char *flabel, WLPtr lword, int first, int climit, int oov)
{
  int i, ufeats[20];
  char *nlabel=flabel;
  WLPtr nword=NULL;
  WHashPtr WHash;
  for (i=1; i < 20; i++) ufeats[i]=0;
  makelower(word);
  hash[0]=get_whash(word,lex);
  WHash=lex->WHash[hash[0]];
  while (WHash != NULL) {
    if (strcmp(word,WHash->WD->label) == 0) {
      nword=WHash->WD;
      break;
      }
    WHash=WHash->next;
  }
  if (nword == NULL) nword=fword;
  else if (lword->caps && fword != NULL && 
	   (fword->tot > 3 || fword->tot >= nword->tot)) nword=fword;
  if (!climit && nword != NULL) return nword;
  if (nword != NULL) nlabel=nword->label;
  ufeats[0]=first;
  if (digit(word[0])) ufeats[13]=1;
  else if (!lword->caps && nword != NULL && nword->tot > 1) return nword;
  else if (nword != NULL && nword->tot > 3) return nword;
  make_unkcat(word,nlabel,flabel,&ufeats[0]);
  nword=find_word(word,lex,hash);
  if (nword!=NULL||!oov) return nword;
  sprintf(word,"oovSUF%d",ufeats[8]);
  return find_word(word,lex,hash);
}

WLPtr find_zword(char *word, LexPtr lex, int *hash, WLPtr fword, char *flabel, WLPtr lword)
{
  int i;
  char *nlabel=flabel;
  WLPtr nword=NULL;
  WHashPtr WHash;
  makelower(word);
  hash[0]=get_whash(word,lex);
  WHash=lex->WHash[hash[0]];
  while (WHash != NULL) {
    if (strcmp(word,WHash->WD->label) == 0) {
      nword=WHash->WD;
      break;
      }
    WHash=WHash->next;
  }
  if (nword == NULL) nword=fword;
  else if (lword->caps && fword != NULL && nword->tot < fword->tot) nword=fword;
  return nword;
}

WLPtr find_nsword(char *word, LexPtr lex, int *hash, char *flabel, int first)
{
  int i, ufeats[20];
  char *nlabel=flabel;
  WLPtr nword=NULL;
  WHashPtr WHash;
  for (i=1; i < 20; i++) ufeats[i]=0;
  makelower(word);
  hash[0]=get_whash(word,lex);
  WHash=lex->WHash[hash[0]];
  while (WHash != NULL) {
    if (strcmp(word,WHash->WD->label) == 0) {
      nword=WHash->WD;
      break;
      }
    WHash=WHash->next;
  }
  if (nword != NULL) return nword;
  ufeats[0]=first;
  if (digit(word[0])) ufeats[13]=1;
  make_unkcat(word,nlabel,flabel,&ufeats[0]);
  nword=find_word(word,lex,hash);
  if (nword != NULL) return nword;
  i=0;
  while (1) {
    find_defunk(word,flabel,i++);
    nword=find_word(word,lex,hash);
    if (nword != NULL) return nword;
  }
}

void wrealloc(LexPtr lex)
{
  int i, init=lex->allocWDs;
  lex->allocWDs += INITWD;
  lex->WDs=realloc(lex->WDs,lex->allocWDs*sizeof(WLPtr));
  for (i=init; i < lex->allocWDs; i++) lex->WDs[i]=NULL;
}

WLPtr insert_word(char *word, LexPtr lex, int hash, int children, int terminal)
{
  int len=strlen(word)+1;
  WHashPtr WHash=lex->WHash[hash], nhash;
  WLPtr wentry=malloc(sizeof(struct WordList));
  wentry->label=malloc(len*sizeof(char));
  strcpy(wentry->label,word);
  wentry->hlabel=wentry->label;
  wentry->slc=NULL;
  wentry->windex=lex->numWDs++;
  if (lex->numWDs > lex->allocWDs) wrealloc(lex);
  wentry->children=children;
  wentry->terminal=terminal;
  wentry->caps=wentry->sym=wentry->oov=wentry->root=0;
  wentry->cclass=wentry->ccount=0;
  wentry->NTs=NULL;
  wentry->unkcat=NULL;
  wentry->tot=0;
  lex->WDs[wentry->windex]=wentry;
  nhash=malloc(sizeof(struct WordHash));
  nhash->next=WHash;
  nhash->WD=wentry;
  lex->WHash[hash]=nhash;
  return wentry;
}

TreePtr get_ptree(TreePtr tree)
{
  if (tree==NULL || tree->prev == NULL) return NULL;
  if (tree->prev->NT->children && !tree->prev->NT->terminal) return tree->prev;
  if (tree->prev->NT->children) return tree->prev->Par;
  return tree->prev->Par->Par;
}

void default_link(WLPtr nt, TreePtr lastree, TreePtr tree)
{
  tree->NT=nt;
  tree->points=0;
  tree->LookWord=tree->LookUpd=NULL;
  tree->LookPT=0;
  tree->Score=tree->PScore=tree->WScore=0;
  tree->prev=lastree;
  tree->LSib=tree->HChild=tree->CCnode=NULL;
  tree->Par=get_ptree(tree);
  if (tree->Par != NULL && tree->Par != lastree) {
    if (!lastree->NT->children) tree->LSib=lastree->Par;      /* STOP */
    else tree->LSib=lastree;                                  /* NT */
  }
}

TreePtr default_tree(WLPtr nt, TreePtr lastree)
{
  TreePtr tree=malloc(sizeof(struct Tree));
  tree->NT=nt;
  tree->points=0;
  tree->LookWord=tree->LookUpd=NULL;
  tree->LookPT=0;
  tree->LookScore=NULL;
  tree->Score=tree->PScore=tree->WScore=0;
  tree->prev=lastree;
  if (lastree != NULL) lastree->points++;
  tree->LSib=tree->HChild=tree->CCnode=NULL;
  tree->featarr=NULL;
  tree->gold=NULL;
  tree->Par=get_ptree(tree);
  if (tree->Par != NULL && tree->Par != lastree) {
    if (!lastree->NT->children) tree->LSib=lastree->Par;      /* STOP */
    else tree->LSib=lastree;                                  /* NT */
  }
  return tree;
}

TreePtr copy_tree(WLPtr nt, TreePtr ntree, int nf)
{
  int i;
  TreePtr tree=default_tree(nt,ntree->prev);
  tree->LookWord=tree->LookUpd=ntree->LookWord;
  tree->LookPT=ntree->LookPT;
  tree->HChild=ntree->HChild;
  tree->CCnode=ntree->CCnode;
  if (ntree->featarr != NULL) {
    tree->featarr = malloc(nf*sizeof(int));
    for (i=0; i < nf; i++) tree->featarr[i] = ntree->featarr[i];
  }
  return tree;
}

SLPtr default_string(char *label, SLPtr lastring)
{
  SLPtr str=malloc(sizeof(struct StringList));
  str->word=NULL;
  str->outlabel=malloc((strlen(label)+1)*sizeof(char));
  strcpy(str->outlabel,label);
  str->leftsym=NULL;
  str->rightsym=NULL;
  str->next=NULL;
  if (lastring != NULL) lastring->next=str;
  return str;
}

SLCPtr lc_node(TreePtr ptree, TreePtr tree)
{
  SLCPtr slc;
  if (ptree == NULL || tree->LSib != NULL || tree->NT->terminal) return NULL; 
  slc=ptree->NT->slc;
  while (slc != NULL) {
    if (slc->LChild == tree->NT && slc->orig) break;
    slc=slc->next;
  }
  return slc;
}

int make_lctran(TreePtr ptree, TreePtr ltree, SLCPtr slc, TreePtr rtree)
{
  TreePtr pstop=rtree, psib=NULL, lstop, lsib, ttree;
  while (pstop != NULL && pstop->Par != ptree) {
    psib=pstop;
    pstop=pstop->prev;
  }
  if (pstop==NULL) failproc("no stop for parent of lc trans");
  lsib=pstop;
  while (lsib->LSib != NULL && lsib->LSib != ltree) lsib = lsib->LSib;
  if (lsib->LSib != ltree) failproc("ltree not lsib of pstop!");
  if (pstop == lsib && ptree->NT != ltree->NT) return 0;
  lstop=lsib->prev;
  if (lstop->Par != ltree) failproc("lstop not child of ltree!");
  if (psib != NULL) psib->LSib=ltree;
  ltree->Par = ptree->Par;
  ltree->LSib = ptree->LSib;
  ltree->NT = ptree->NT;
  ltree->prev = ptree->prev;
  if (pstop == lsib) {      /* remove unary transforms (NT (NT ...)) */
    ttree = ptree;
    ptree = lstop->LSib;
    free_tree(ttree);
    ttree = lstop;
    lstop = lstop->prev;
    free_tree(ttree);
  }
  else {
    ptree->Par = ltree;
    ptree->LSib = lstop->LSib;
    ptree->NT = slc->NewCat;
    ptree->prev = lstop->prev;
    lstop->Par = ptree;
    lstop->LSib = pstop->LSib;
    lstop->prev = pstop->prev;
    lsib->prev = ptree;
    lsib->LSib = NULL;
  }
  pstop->Par = ltree;
  pstop->LSib = ptree;
  pstop->prev = lstop;
  return 1;
}

SLCPtr delc_node(TreePtr tree)
{
  SLCPtr slc;
  if (tree == NULL || tree->LSib == NULL || tree->Par == NULL) return NULL;
  slc=tree->LSib->NT->slc;
  if (slc != NULL && !slc->orig) return slc;
  return NULL;
}

void lc_transform(TreePtr tree)
{
  TreePtr ltree=tree;
  SLCPtr slc;
  while (ltree != NULL) {
    slc=lc_node(ltree->Par,ltree);
    if (slc == NULL || !make_lctran(ltree->Par,ltree,slc,tree))
      ltree=ltree->prev;
  }
  ltree=tree;
  while (ltree != NULL) {
    slc=delc_node(ltree);
    if (slc != NULL) ltree->Par->NT=slc->LChild;
    ltree=ltree->prev;
  }
}

void unmake_lctran(TreePtr pstop, TreePtr psib, SLCPtr slc, TreePtr rtree)
{
  TreePtr ptree=pstop->LSib,ltree=pstop->Par,lstop=pstop->prev,lsib=lstop->LSib;
  while (lsib->LSib != NULL) lsib = lsib->LSib;
  if (psib != NULL) psib->LSib=ptree;
  lsib->LSib = ltree;
  lsib->prev = lstop;
  pstop->Par = ptree;
  pstop->LSib = lstop->LSib;
  pstop->prev = lstop->prev;
  lstop->Par = ltree;
  lstop->LSib = ptree->LSib;
  lstop->prev = ptree->prev;
  ptree->Par = ltree->Par;
  ptree->LSib = ltree->LSib;
  ptree->prev = ltree->prev;
  ptree->NT = slc->NewCat;
  ltree->Par = ptree;
  ltree->LSib = NULL;
  ltree->prev = ptree;
}

void lc_detransform(TreePtr tree)
{
  TreePtr ltree=tree,lastree=NULL;
  SLCPtr slc;
  while (ltree != NULL) {
    slc=delc_node(ltree);
    if (slc != NULL) unmake_lctran(ltree,lastree,slc,tree);
    else {
      lastree=ltree;
      ltree=ltree->prev;
    }
  }
}

TreePtr best_hcat(TreePtr rchild)
{
  TreePtr lchild=rchild, lnt=NULL, hcat=NULL, psame=NULL, pns=NULL, nchild;
  WLPtr par=rchild->Par->NT;
  while (lchild != NULL) {
    nchild = lchild->LSib;
    while (nchild != NULL && nchild->NT->sym) nchild = nchild->LSib;
    if (!lchild->NT->terminal && lchild->NT->children && nchild != NULL &&
	nchild->NT->terminal) lnt=lchild;
    lchild = lchild->LSib;
  }
  if (lnt != NULL && lnt->LSib != NULL) {
    lnt = lnt->LSib;
    while (lnt != NULL) {
      if (lnt->NT->terminal && !lnt->NT->sym) {
	if (lnt->NT->label[1] == par->label[1] && psame == NULL) psame=lnt;
	if (lnt->NT->label[1] != par->label[1] && pns == NULL) pns=lnt;	
      }
      lnt = lnt->LSib;
    }
    if (psame != NULL) return psame;
    if (pns != NULL) return pns;
  }
  lchild = rchild;
  while (lchild != NULL) {
    if (!lchild->NT->sym) return lchild;
    lchild = lchild->LSib;
  }
  return rchild;
}

TreePtr get_hnode(TreePtr tree)
{
  TreePtr hcand, xcand;
  if (tree==NULL || tree->LSib == NULL) return NULL;
  if (!tree->NT->children && tree->LSib->LSib == NULL) return tree->LSib;
  if (tree->LSib->LSib == NULL) return NULL;
  if (tree->LSib->HChild != NULL) return tree->LSib->HChild;
  xcand=tree->LSib;
  while (xcand != NULL && xcand->NT->sym) xcand=xcand->LSib;
  if (xcand == NULL || xcand->LSib == NULL) {
    if (!tree->NT->children && xcand==NULL) return tree->LSib;
    if (!tree->NT->children) return xcand;
    return NULL;
  }
  hcand = best_hcat(xcand);
  if (tree->NT->children && hcand != xcand->LSib) return NULL;
  return hcand;
}

void fill_hnodes(TreePtr tree)
{
  if (tree==NULL) return;
  fill_hnodes(tree->prev);
  tree->HChild=get_hnode(tree);
}

void fill_ccnodes(TreePtr tree,int all)
{
  if (tree == NULL) return;
  if (all) fill_ccnodes(tree->prev,all);
  tree->CCnode = NULL;
  if (tree->LSib != NULL && tree->HChild != NULL) tree->CCnode = tree->HChild;
  else if (tree->LSib != NULL) tree->CCnode = tree->LSib;
  else if (tree->Par != NULL) tree->CCnode = tree->Par->CCnode;
}

void fill_farrays(TreePtr tree, TreePtr otree, FtrPtr feats, int all)
{
  int q, i, j, k, gind;
  char c, lc;
  TreePtr ftree;
  TWPtr tfunct, hfunct;
  if (tree == NULL) return;
  if (all) fill_farrays(tree->prev,otree,feats, all);
  if (tree->featarr==NULL) tree->featarr = malloc(feats->numfeats*sizeof(int));
  for (i=0; i < feats->numfeats; i++) {
    tfunct=feats->featfuncts[i];
    if (tfunct->suffmove == NULL) gind = -1;
    else gind = tfunct->suffmove->findex;
    ftree = tree;
    q=c=0;
    if (tfunct->SQuery >= 0) {
      if (tfunct->SQuery >= tfunct->numqueries) 
	failproc("not enough queries for storage");
      tfunct->WQueries[tfunct->SQuery]=NULL;
    }
    for (j=0; j < tfunct->numoves; j++) {
      lc=c;
      c=tfunct->moves[j];
      if (valop(c)) tree->featarr[i] = getvalop(ftree,c);
      else if (c == 'q') 
	ftree = askquery(ftree,q++,tfunct->SQuery,tfunct->WQueries,tree);
      else if (c == 's') {
	if (ftree != NULL && tfunct->SQuery >= 0) 
	  tfunct->WQueries[tfunct->SQuery]=ftree->NT;
      }
      else ftree = gettwop(ftree,tree,c,lc);
      if (j==0 && gind >= 0) {
	if (ftree == NULL) tree->featarr[i] = -1;
	else if (ftree->featarr == NULL) continue;
	else tree->featarr[i] = ftree->featarr[gind];
	break;
      }
    }
  }
}

void fill_lclfarrays(TreePtr tree, TreePtr otree, FtrPtr feats)
{
  int q, i, j;
  char c, lc;
  TreePtr ftree;
  TWPtr tfunct, hfunct;
  if (tree == NULL) return;
  for (i=0; i < feats->numfeats; i++) {
    tfunct=feats->featfuncts[i];
    if (tfunct->suffmove != NULL) continue;
    c=tfunct->moves[0];
    if (valop(c)) tree->featarr[i] = getvalop(tree,c);
  }
}

TreePtr punctmove(TreePtr tree, TreePtr last, TreePtr plast, TreePtr orig)
{
  TreePtr ftree;
  if (tree == NULL) return tree;
  if (tree->LSib != NULL && (last->NT->terminal || last->NT->children))
    return tree;
  if (tree->Par == NULL || tree->Par->Par == orig->Par) return tree;
  if (tree->LSib == NULL) {
    tree->prev=tree->Par->prev;
    tree->LSib=tree->Par->LSib;
    tree->Par->prev=tree;
    tree->Par->LSib=tree;
    tree->Par=tree->Par->Par;
    if (!last->NT->terminal && !last->NT->children) {
      plast->prev=plast->LSib=tree;
      ftree=last->Par;
      free_tree(last);
      free_tree(ftree);
    }
    else {
      last->LSib=NULL;
      last->prev=last->Par;
    }    
    return orig;
  }
  last->LSib=tree->LSib;
  last->prev=tree->prev;
  tree->prev=last;
  plast->prev=plast->LSib=tree;
  tree->LSib=tree->Par;
  tree->Par=tree->Par->Par;
  return orig;
}

TreePtr getprevterm(TreePtr tree)
{
  if (tree==NULL) return NULL;
  if (tree->NT->terminal) return tree;
  return getprevterm(tree->prev);
}

TreePtr punctnorm(TreePtr tree)
{
  TreePtr ttree=tree, last=NULL, plast=NULL, tlast=NULL, nterm, punc=NULL;
  while (ttree != NULL) {
    if (ttree->NT->sym) {
      nterm=getprevterm(ttree->prev);
      if (nterm != NULL || last==NULL || unary(ttree->prev,tree)) 
	ttree=punctmove(ttree,last,plast,tree);
      else {
	punc=ttree;
	last->prev=ttree->prev;
	last->LSib=ttree->LSib;
	ttree=last;
	last=plast;
      }
    }
    plast=last;
    last=ttree;
    ttree=ttree->prev;
  }
  return punc;
}

TreePtr read_tree(char *readtree, ParPtr par, int oov)
{
  int i=0, hash, len=strlen(readtree), firstw=1;
  char readtoken[MAXLABLEN], holdtoken[MAXLABLEN];
  WLPtr nt, word;
  TreePtr tree=NULL,lastree;
  SLPtr lastring=NULL;
  while (i < len) {
    lastree=tree;
    wordtoken(readtree,&readtoken[0],&i,len,1);
    if (i >= len) continue;
    nt=must_find_word(&readtoken[0],par->lexicon,&hash);
    tree=default_tree(nt,lastree);
    if (tree->NT->children && tree->NT->terminal) {
      wordtoken(readtree,&readtoken[0],&i,len,1);
      if (i >= len || readtree[i] != ')') {
 	sprintf(&readtoken[0],"misparse: %s%s\n",readtree,&readtree[i]);
	failproc(&readtoken[0]);
      }
      i++;
      lastring=tree->LookWord=default_string(&readtoken[0],lastring);
      nt=find_word(&readtoken[0],par->lexicon,&hash);
      if (nt == NULL) {
	nt=find_nword(&readtoken[0],par->lexicon,&hash,nt,
		      tree->LookWord->outlabel,tree->NT,firstw,0,1);
	if (nt == NULL || (nt->oov && !oov)) {
	  strcpy(holdtoken,tree->LookWord->outlabel);
	  nt=find_zword(&holdtoken[0],par->lexicon,&hash,NULL,
			tree->LookWord->outlabel,tree->NT);
	}
	if (nt == NULL) {
	  fprintf(stderr,"%s\n",tree->LookWord->outlabel);
	  failproc("failed to insert word!");
	}
      }
      if (firstw && !nt->sym) firstw=0;
      tree->LookWord->word=nt;
    }
  }
  return tree;
}

void fill_tree(TreePtr tree, ParPtr par)
{
  int j;
  TreePtr lastree=tree;
  SLPtr lastring=NULL;
  fill_hnodes(tree);
  fill_ccnodes(tree,1);
  while (lastree != NULL) {
    if (lastree->NT->terminal) {
      j=lastree->NT->windex;
      lastring=lastree->LookWord;
    }
    else {
      lastree->LookWord = lastring;
      lastree->LookPT = j;
    }
    lastree = lastree->prev;
  }
  fill_farrays(tree,tree,par->features,1);
}

void putpuncback(TreePtr tree, TreePtr punc)
{
  TreePtr ttree=tree->prev;
  while (ttree->LSib != NULL) ttree=ttree->LSib;
  punc->Par=ttree->Par;
  punc->LSib=ttree->LSib;
  punc->prev=ttree->prev;
  ttree->prev=punc;
  ttree->LSib=punc;
}

void remove_dblq(TreePtr tree)
{
  TreePtr ttree=tree, plast=NULL, last=NULL;
  while (ttree != NULL) {
    if (ttree->NT->children && ttree->NT->terminal) {
      if (quotes(ttree->LookWord->outlabel)) {
	last->prev=ttree->prev;
	last->LSib=ttree->LSib;
	ttree=last;
	last=plast;
      }
    }
    plast=last;
    last=ttree;
    ttree=ttree->prev;
  }
}

TreePtr build_tree(char *readtree, ParPtr par, int oov)
{
  TreePtr tree=NULL, punc=NULL;
  tree=read_tree(readtree,par,oov);
  remove_dblq(tree);
  punc=punctnorm(tree);
  lc_transform(tree);
  if (punc!=NULL) putpuncback(tree,punc);
  fill_tree(tree,par);
  return tree;
}

SLPtr get_tstring(ParPtr par, TreePtr tree)
{
  SLPtr string=NULL;
  while (tree != NULL) {
    if (tree->LookWord != NULL) string=tree->LookWord;
    tree=tree->prev;
  }
  return string;
}

GLPtr get_gold(TreePtr tree, GLPtr pgold)
{
  GLPtr gold;
  if (tree->prev==NULL) return pgold;
  gold=malloc(sizeof(struct GoldList));
  gold->move=tree->NT;
  gold->best=NULL;
  gold->gtree=tree;
  gold->next=pgold;
  return get_gold(tree->prev,gold);
}

void free_gold(GLPtr gold)
{
  GLPtr tmp;
  while (gold!=NULL) {
    tmp=gold;
    gold=gold->next;
    free(tmp);
  }
}

GLPtr goldthere(TreeHPtr ntrees,int last)
{
  GLPtr gold=NULL;
  if (ntrees==NULL) return gold;
  if (!last && ntrees->tree->gold!=NULL) return ntrees->tree->gold;
  if (last && ntrees->tree->prev->gold!=NULL) return ntrees->tree->prev->gold;
  gold=goldthere(ntrees->next,last);
  if (gold != NULL) return gold;
  gold=goldthere(ntrees->right,last);
  if (gold != NULL) return gold;
  return goldthere(ntrees->left,last);
}

TreePtr goldcand(TreeHPtr ntrees,int last)
{
  TreePtr gold=NULL;
  if (ntrees==NULL) return gold;
  if (!last && ntrees->tree->gold!=NULL) return ntrees->tree;
  if (last && ntrees->tree->prev->gold!=NULL) return ntrees->tree;
  gold=goldcand(ntrees->next,last);
  if (gold != NULL) return gold;
  gold=goldcand(ntrees->right,last);
  if (gold != NULL) return gold;
  return goldcand(ntrees->left,last);
}

GLPtr adv_move(GLPtr gold)
{
  while (gold->next != NULL && !gold->move->terminal) gold=gold->next;
  return gold;
}

GLPtr get_rgold(TreeHPtr ntrees,GLPtr gold,int last)
{
  GLPtr ngold;
  if (gold->best != NULL || gold->next == NULL) return gold;
  ngold=goldthere(ntrees,last);
  if (ngold!=NULL) return ngold;
  ngold=adv_move(gold);
  if (ntrees==NULL) return ngold;
  if (ntrees->tree->LookWord!=ngold->gtree->LookWord) {
    fprintf(stderr,"oops, rgold mismatch: %s %s %s %s\n",ntrees->tree->NT->label, ntrees->tree->LookWord->word->label,ngold->move->label,ngold->gtree->LookWord->word->label);
    exit(1);
  }
  ngold->best=ntrees->tree;
  ngold->best->points++;
  return ngold;
}

TreePtr build_copytree(TreePtr btree)
{
  TreePtr tree, prev=NULL;
  if (btree==NULL) return prev;
  prev=build_copytree(btree->prev);
  tree=default_tree(btree->NT,prev);
  tree->LookWord=tree->LookUpd=btree->LookWord;
  tree->Score=btree->Score;
  tree->WScore=btree->WScore;
  tree->PScore=btree->PScore;
  return tree;
}

void show_tree(FILE *fp, TreePtr tree, int first)
{
  if (tree == NULL) return;
  show_tree(fp,tree->prev,0);
  if (tree->NT->children) {
    if (tree->LSib != NULL) fprintf(fp," ");
    if (tree->NT->terminal && tree->LookWord->leftsym!= NULL) {
      if (tree->LookWord->leftsym[0]!='"')
	fprintf(fp,"(%s %s) ",tree->LookWord->leftsym,
		tree->LookWord->leftsym);
      else fprintf(fp,"('' %s) ",tree->LookWord->leftsym);
    }
    fprintf(fp,"%s ",tree->NT->label);
    if (tree->NT->terminal) {
      fprintf(fp,"%s)",tree->LookWord->outlabel);
      if (tree->LookWord->rightsym!= NULL)
	fprintf(fp," ('' %s) ",tree->LookWord->rightsym);      
    }
  }
  else fprintf(fp,")");
  if (first) fprintf(fp,"\n");
}

TreePtr firstchild(TreePtr tree, TreePtr orig)
{
  while (orig != NULL && orig->Par != tree) orig=orig->prev;
  if (orig==NULL) {
    fprintf(stderr,"oops, no children\n");
    exit(1);
  }
  return orig;
}

int unary(TreePtr tree, TreePtr orig)
{
  orig=firstchild(tree,orig);
  if (orig->LSib == NULL || orig->LSib->LSib == NULL) return 1;
  return 0;
}

int nppos(TreePtr tree, TreePtr orig)
{
  orig=firstchild(tree,orig);
  if (orig->LSib==NULL) return 0;
  if (strcmp(orig->LSib->NT->label,"(POS")==0) return 1;
  return 0;
}

int npbase(TreePtr tree, TreePtr orig)
{
  orig=firstchild(tree,orig);
  orig=orig->LSib;
  while (orig != NULL) {
    if (!orig->NT->terminal) return 0;
    orig=orig->LSib;
  }
  return 1;
}

int gapess(TreePtr tree, TreePtr orig)
{
  orig=firstchild(tree,orig);
  if (orig->LSib == NULL || orig->LSib->LSib != NULL) return 0;
  if (strcmp(orig->LSib->NT->label,"(VP")==0) return 1; 
  return 0;
}

char *verbhead(TreePtr tree, TreePtr orig)
{
  orig=firstchild(tree,orig);
  while (orig->LSib != NULL) orig=orig->LSib;
  if (orig != NULL) { 
    if (strncmp(orig->NT->label,"(VB",3)==0) 
      return orig->NT->label; 
    if (strncmp(orig->NT->label,"(AU",3)==0) 
      return orig->NT->label; 
    if (strncmp(orig->NT->label,"(TO",3)==0) 
      return orig->NT->label; 
    if (strncmp(orig->NT->label,"(MD",3)==0) 
      return orig->NT->label; 
  }
  return NULL;
}

void show_kmod(FILE *fp, TreePtr tree, TreePtr orig, int first)
{
  char *hd;
  if (tree == NULL) return;
  show_kmod(fp,tree->prev,orig,0);
  if (tree->NT->children) {
    if (tree->LSib != NULL) fprintf(fp," ");
    fprintf(fp,"%s",tree->NT->label);
    if (tree->NT->terminal) {
      if (strcmp(tree->NT->label,"(DT")==0 ||
	  strcmp(tree->NT->label,"(RB")==0)
	if (unary(tree->Par,orig)) fprintf(fp,"~U");
      fprintf(fp," ");
      fprintf(fp,"%s)",tree->LookWord->outlabel);
    }
    else {
      if (strcmp(tree->NT->label,"(NP")==0) {
	if (nppos(tree,orig)) fprintf(fp,"~P");
	if (npbase(tree,orig)) fprintf(fp,"~B");
      }
      else if (strcmp(tree->NT->label,"(S")==0) {
	if (gapess(tree,orig)) fprintf(fp,"~G");
      }
      else if (strcmp(tree->NT->label,"(VP")==0) {
	hd=verbhead(tree,orig);
	if (hd!=NULL) fprintf(fp,"~%s",&hd[1]);
      }
      if (unary(tree,orig)) fprintf(fp,"~U");
      fprintf(fp," ");
    }
  }
  else fprintf(fp,")");
  if (first) fprintf(fp,"\n");
}

int show_ttree(FILE *fp, TreePtr tree, int first)
{
  int i, sp;
  TreePtr lsib;
  if (tree == NULL) return -1;
  sp=show_ttree(fp,tree->prev,0);
  if (tree->NT->children) {
    lsib=tree->LSib;
    if (lsib==NULL) sp+=1;
    if (tree->NT->terminal && (lsib==NULL || (lsib->NT->terminal && !lsib->NT->sym)))
      fprintf(fp," ");
    else {
      if (sp > 0) fprintf(fp,"\n");
      for (i=0; i < sp; i++) {
	if (i%2==0) fprintf(fp,"  ");
	else fprintf(fp,"  ");
      }
    }
    fprintf(fp,"%s",tree->NT->label);
    if (tree->NT->terminal) fprintf(fp," %s)",tree->LookWord->outlabel);
  }
  else {
    sp-=1;
    fprintf(fp,")");
  }
  if (first) fprintf(fp,"\n");
  return sp;
}

void show_defstring(FILE *fp, SLPtr string, LexPtr lex)
{
  int i, bpos;
  float best;
  fprintf(fp,"(TOP");
  while (string != NULL) {
    bpos=-1;
    for (i=1; i < lex->numPTs; i++) 
      if (bpos < 0 || string->word->NTs[i] < best) {
	bpos=i;
	best=string->word->NTs[i];
      }
    fprintf(fp," %s %s)",lex->WDs[bpos]->label,string->outlabel);
    string=string->next;
  }
  fprintf(fp,")\n");
}

void free_slist(SLPtr slist)
{
  if (slist == NULL) return;
  free_slist(slist->next);
  free(slist->outlabel);
  if (slist->leftsym!=NULL) free(slist->leftsym);
  if (slist->rightsym!=NULL) free(slist->rightsym);
  free(slist);
}

void free_tree(TreePtr tree)
{
  if (tree->LookScore != NULL) free(tree->LookScore);
  if (tree->featarr != NULL) free(tree->featarr);
  free(tree);
}

void free_trees(TreePtr tree, SLPtr tstring, int dostr)
{
  SLPtr lstring=tstring;
  if (tree == NULL) {
    if (tstring != NULL) free_slist(tstring);
    return;
  }
  if (!dostr) lstring=NULL;
  else if (tree->NT->terminal) lstring=tree->LookWord;
  if (tree->points==0) {
    if (tree->prev != NULL) tree->prev->points--;
    free_trees(tree->prev,lstring,dostr);
    free_tree(tree);
  }
}

void free_theap(TreeHPtr trees)
{
  if (trees==NULL) return;
  if (trees->tree != NULL) trees->tree->points--;
  free_trees(trees->tree,NULL,0);
  free(trees);
}

void free_theaps(TreeHPtr trees)
{
  if (trees==NULL) return;
  free_theaps(trees->next);
  free_theaps(trees->left);
  free_theaps(trees->right);
  free_theap(trees);
}

TreeHPtr get_lunder(TreeHPtr htree)
{
  TreeHPtr ret;
  if (htree->left==NULL) {
    ret=htree->right;
    htree->right=NULL;
  }
  else if (htree->right==NULL) {
    ret=htree->left;
    htree->left=NULL;
  }
  else if (htree->right->under <= htree->left->under) {
    ret=htree->right;
    htree->right=NULL;
  }
  else {
    ret=htree->left;
    htree->left=NULL;
  }
  if (ret == NULL) htree->under=1;
  else htree->under-=ret->under;
  fill_tnull(htree);
  return ret;
}

void fill_tnull(TreeHPtr tree)
{
  TreeHPtr htree;
  if (tree->right != NULL && tree->left != NULL) return;
  if (tree->right == NULL && tree->left == NULL) return;
  if (tree->right==NULL) {
    htree=tree->left;
    tree->right=get_lunder(htree);
  }
  else {
    htree=tree->right;
    tree->left=get_lunder(htree);
  }
}

TreeHPtr null_thinsert(TreeHPtr trees, TreeHPtr tree)
{
  trees->right = tree->right;
  trees->left = tree->left;
  tree->right = tree->left = NULL;
  tree->under = 1;
  if (trees->right==NULL) trees->right=tree;
  else if (trees->left==NULL) trees->left=tree;
  else {
    if (trees->right->under <= trees->left->under) {
      tree->right=trees->right;
      trees->right=tree;
    }
    else {
      tree->right=trees->left;
      trees->left=tree;
    }
    tree->under+=tree->right->under;	
    fill_tnull(tree);
  }
  return trees;
}

TreeHPtr insert_theap(TreeHPtr tree, TreeHPtr trees)
{
  TreeHPtr whichild;
  if (trees==NULL) return tree;
  if (tree==NULL) return trees;
  if (tree->score == trees->score && tree->right==NULL && 
      tree->left==NULL && tree->next == NULL) {
    tree->next=trees->next;
    trees->next=tree;
    return trees;
  }
  if (tree->score < trees->score) return insert_theap(trees, tree);
  trees->under+=tree->under;
  if (trees->right == NULL && trees->left == NULL) 
    return null_thinsert(trees,tree);
  if (trees->right==NULL) trees->right=tree;
  else if (trees->left==NULL) trees->left=tree;
  else if (trees->right->under <= trees->left->under)
    trees->right=insert_theap(tree,trees->right);
  else trees->left=insert_theap(tree,trees->left);
  return trees;
}

TreeHPtr perc_heap(TreeHPtr ctrees)
{
  TreeHPtr ret=NULL;
  if (ctrees==NULL) return ret;
  if (ctrees->next != NULL) {
    ret=ctrees->next;
    ret->left=ctrees->left;
    ret->right=ctrees->right;
    ret->under=ctrees->under;
  }
  else if (ctrees->left==NULL) ret=ctrees->right;
  else if (ctrees->right==NULL) ret=ctrees->left;
  else ret=insert_theap(ctrees->left,ctrees->right);
  ctrees->left=ctrees->right=ctrees->next=NULL;
  ctrees->under=1;
  return ret;
}

TreeHPtr insert_theapns(TreeHPtr tree, TreeHPtr trees)
{
  if (trees==NULL) return tree;
  if (tree==NULL) return trees;
  if (tree->score < trees->score) return insert_theapns(trees, tree);
  if (tree->right != NULL && trees->right != NULL) 
    failproc("invalid no-sort insert");
  if (trees->right != NULL) tree->right = trees->right;
  trees->right = tree;
  return trees;
}

TreeHPtr default_theap(TreePtr tree, float score, int pt, int rank)
{
  TreeHPtr theap=malloc(sizeof(struct TreeHeap));
  theap->tree=tree;
  if (tree != NULL) tree->points++;
  theap->score=score;
  theap->left=theap->right=theap->next=NULL;
  theap->under=1;
  theap->rank=rank;
  theap->preterm=pt;
  return theap;
}

TreeHPtr init_parse(TDParPtr PC, GLPtr gold)
{
  int root=PC->parser->lexicon->root;
  float score=0;
  PC->parser->norm=PC->parser->snorm=PC->parser->cond=PC->parser->lcond=score;
  TreePtr tree=default_tree(PC->parser->lexicon->WDs[root],NULL);
  TreeHPtr theap=default_theap(tree,score,-1,1);
  tree->gold=gold;
  return theap;
}

float *defscores(float *scores, int numind, int numalloc)
{
  int i;
  float *nscores=scores;
  if (nscores==NULL) nscores=malloc(numind*sizeof(float));
  else nscores=realloc(nscores,numind*sizeof(float));
  for (i=numalloc; i < numind; i++) nscores[i]=0;
  return nscores;
}

int *deflupd(int *lupd, int numind, int numalloc)
{
  int i, *nlupd=lupd;
  if (nlupd==NULL) nlupd=malloc(numind*sizeof(float));
  else nlupd=realloc(nlupd,numind*sizeof(float));
  for (i=numalloc; i < numind; i++) nlupd[i]=0;
  return nlupd;
}

CFDataPtr def_cfeat(ScorePtr scr, int clbl, int lam)
{
  int i, ind;
  CFDataPtr cf=malloc(sizeof(struct CFeatData));
  cf->label=clbl;
  cf->idx=scr->numind++;
  if (!lam) return cf;
  if (scr->numind >= scr->numalloc) {
    ind=scr->numind+100;
    scr->scores=defscores(scr->scores,ind,scr->numalloc);
    if (scr->clams==NULL) scr->clams=malloc(ind*sizeof(CLamPtr));
    else scr->clams=realloc(scr->clams,ind*sizeof(CLamPtr));
    for (i=scr->numalloc; i < ind; i++) scr->clams[i]=NULL;
    scr->numalloc=ind;
  }
  scr->scores[cf->idx]=0;
  scr->clams[cf->idx]=malloc(sizeof(struct CDataLam));
  scr->clams[cf->idx]->sbatch=0;
  scr->clams[cf->idx]->lambda=NULL;
  scr->clams[cf->idx]->pscore=NULL;
  return cf;
}

void insert_cfeat(FtrPtr feats, FDataPtr fdata, int clbl, int lam)
{
  int i, alloc;
  alloc=fdata->numc+1;
  if (fdata->cfeat == NULL) fdata->cfeat = malloc(alloc*sizeof(CFDataPtr));
  else fdata->cfeat = realloc(fdata->cfeat,alloc*sizeof(CFDataPtr));
  fdata->cfeat[fdata->numc++]=def_cfeat(feats->scores,clbl,lam);
}

int get_fhash(int state, int label, int level, int schema, FtrPtr features)
{
  int i, hash;
  long lhash=0;
  lhash += fabs(state+2)*STATEM;
  lhash += fabs(label+2)*LABELM;
  lhash += fabs(level+2)*LEVELM;
  lhash += fabs(schema+2)*SCHEMM;
  hash = lhash%features->numFH;
  return hash;
}

FHashPtr find_feat(int state, int label, int level, int schema, FtrPtr features, int *hash)
{
  FHashPtr fhash;
  hash[0]=get_fhash(state,label,level,schema,features);
  fhash=features->FHash[hash[0]];
  while (fhash != NULL) {
    if (fhash->feat->state == state && fhash->feat->label == label && 
	fhash->feat->level == level && fhash->feat->schema == schema) 
      return fhash;
    fhash=fhash->next;
  }
  return NULL;
}

FHashPtr insert_feat(int state, int label, int level, int schema, FtrPtr features, int hash, int afeat)
{
  FHashPtr fhash=malloc(sizeof(struct FeatHash));
  FeatPtr feat=fhash->feat=malloc(sizeof(struct FeatStruct));
  feat->state=state;
  feat->label=label;
  feat->level=level;
  feat->schema=schema;
  if (afeat) {
    feat->dstate=features->fstates++;
    feat->data=malloc(sizeof(struct FeatData));
    feat->data->score=0;
    feat->data->bscore=-log(0.5);
    feat->data->flam=NULL;
    feat->data->prefix=NULL;
    feat->data->numc=0;
    feat->data->cfeat=NULL;
  }
  else feat->data=NULL;
  fhash->next = features->FHash[hash];
  features->FHash[hash] = fhash;
  return fhash;
}

void presum_fsc(TDTrPtr TC, CFDataPtr cfeat, CFDataPtr bfeat, double s0, double s1)
{
  int i;
  double sum0, sum1;
  ScorePtr scr=TC->parser->features->scores;
  for (i=0; i < TC->xvals; i++) {
    sum0=s0;
    sum1=s1;
    sum0*=exp(-scr->lrn->fscores[i]->rsc[bfeat->idx]);
    sum1*=exp(-scr->lrn->fscores[i]->rsc[cfeat->idx]);
    scr->lrn->fscores[i]->rsc[cfeat->idx] = -log(sum0+sum1);
  }
}

int presum_fhash(FtrPtr feats, FHashPtr fhash, int k, int m, TDTrPtr TC)
{
  int j, bstate, hash, lbl;
  double norm, sum0, sum1;
  FHashPtr bhash;
  ScorePtr scr=feats->scores;
  FDataPtr data;
  FeatPtr bfeat;
  CFDataPtr cfeat;
  while (fhash != NULL) {
    data=fhash->feat->data;
    if (data != NULL && fhash->feat->level == m && data->prefix != NULL) {
      bfeat=data->prefix;
      bstate=bfeat->dstate;
      for (j=0; j < data->numc; j++) {
	lbl=data->cfeat[j]->label;
	sum0=exp(-data->bscore);
	sum1=1-sum0;
	bhash=find_feat(bfeat->dstate,lbl,bfeat->level+1,bfeat->schema+1,feats,&hash);
	if (bhash == NULL) failproc("no backoff arc");
	presum_fsc(TC,data->cfeat[j],bfeat->data->cfeat[bhash->feat->dstate],
		   sum0,sum1);
	sum0*=exp(-scr->scores[bfeat->data->cfeat[bhash->feat->dstate]->idx]);
	sum1*=exp(-scr->scores[data->cfeat[j]->idx]);
	scr->scores[data->cfeat[j]->idx] = -log(sum0+sum1);
      }
      verbose_count(TC->verbose,++k,1000000);
    }
    fhash=fhash->next;
  }
  return k;
}

int presum_feats(TDTrPtr TC)
{
  int max=0, i, k=0, m;
  FtrPtr feats=TC->parser->features;
  if (TC->verbose) fprintf(stderr,"Grammar: prsum "); 
  for (i=0; i < feats->numschema; i++) 
    if (feats->schema[i]->numfuncts > max) max=feats->schema[i]->numfuncts;
  for (m=0; m < max; m++)
    for (i=0; i < feats->numFH; i++) 
      k=presum_fhash(feats,feats->FHash[i],k,m,TC);
  if (TC->verbose) fprintf(stderr," %d\n",k);
  return max;
}

int set_bscore(TDTrPtr TC)
{
  int i, k=0, numFH=TC->parser->features->numFH, maxlvl=0;
  FHashPtr fhash;
  FeatPtr feat;
  FDataPtr data;
  for (i=0; i < numFH; i++) {
    fhash=TC->parser->features->FHash[i];
    while (fhash != NULL) {
      feat=fhash->feat;
      data=feat->data;
      fhash=fhash->next;
      if (data == NULL) continue;
      if (feat->level > maxlvl) maxlvl=feat->level;
      data->flam=malloc(sizeof(struct FDataLam));
      data->flam->next=NULL;
      data->flam->bucket=data->flam->ind=0;
      data->bscore=data->score;
      data->bscore/=data->numc;
      verbose_count(TC->verbose,++k,1000000);
    }
  }
  if (TC->verbose) fprintf(stderr," %d   ",k);  
  return maxlvl;
}

int binscore(FDataPtr data, FeatPtr feat, int numsch, int bins)
{
  int score;
  if (data->bscore >= MAXLABLEN) score=MAXLABLEN-1;
  else score=floor(data->bscore);
  if (score < 0) score=0;
  score += feat->level*MAXLABLEN;
  data->flam->ind=feat->level*numsch*bins+feat->schema*bins;
  return score;
}

FDataPtr bindatasort(FDataPtr *datas, FDataPtr *lowdatas, int maxlvl)
{
  int k=(maxlvl+1)*MAXLABLEN-1;
  FDataPtr data=NULL,ldata;
  while (k >= 0) {
    if (k >=0 && datas[k] != NULL) {
      if (data==NULL) data=datas[k];
      else ldata->flam->next=datas[k];
      ldata=lowdatas[k];
    }
    k--;
  }
  return data;
}

FDataPtr order_bins(TDTrPtr TC, FDataPtr *datas, FDataPtr *lowdatas, float *counts, int bins, int numsch, int maxlvl)
{
  int i, k=0, numFH=TC->parser->features->numFH, score;
  FHashPtr fhash;
  FeatPtr feat;
  FDataPtr data, ldata;
  for (i=0; i < numFH; i++) {
    fhash=TC->parser->features->FHash[i];
    while (fhash != NULL) {
      feat=fhash->feat;
      data=feat->data;
      fhash=fhash->next;
      if (data == NULL || data->bscore < 2) continue;
      score=binscore(data,feat, numsch, bins);
      counts[data->flam->ind]+=data->score;
      verbose_count(TC->verbose,++k,1000000);
      if (datas[score]==NULL || datas[score]->bscore <= data->bscore) {
	data->flam->next=datas[score];
	datas[score]=data;
	if (lowdatas[score]==NULL) lowdatas[score]=data;
	continue;
      }
      if (data->bscore <= lowdatas[score]->bscore) {
	lowdatas[score]->flam->next=data;
	lowdatas[score]=data;
	continue;
      }
      ldata=datas[score];
      while (ldata->flam->next != NULL && ldata->flam->next->bscore > data->bscore) 
	ldata = ldata->flam->next;
      data->flam->next=ldata->flam->next;
      ldata->flam->next=data;
    }
  }
  if (TC->verbose) fprintf(stderr," %d\n",k);
  return bindatasort(datas,lowdatas,maxlvl);
}

FDataPtr first_data(int j, FDataPtr *datas)
{
  int ind=(j+1)*MAXLABLEN-1, min=j*MAXLABLEN-1;
  while (ind > min && datas[ind] == NULL) ind--;
  if (ind <= min) return NULL;
  return datas[ind];
}

FDataPtr updata(FDataPtr tdata, FDataPtr ldata, int bins, int thismin, int *thisbin, float *cands, float *thiscands, float *totcands, float *lscore)
{
  if (ldata != NULL) ldata->flam->next=tdata->flam->next;
  if (thisbin[0]<bins-1 && totcands[0]>=thismin && tdata->bscore!=lscore[0] &&
      (cands[0]==0 || thiscands[0] > cands[0])) {
    thisbin[0]++;
    thiscands[0]=0;
  }
  cands[0]=totcands[0]/(bins-thisbin[0]);
  if (cands[0] < thismin) cands[0]=thismin;
  totcands[0]-=tdata->score;
  thiscands[0]+=tdata->score;
  tdata->flam->bucket=thisbin[0];
  lscore[0]=tdata->bscore;
  return tdata->flam->next;
}

void assign_bins(FDataPtr data, FDataPtr *datas, FDataPtr *lowdatas, float *counts, int bins, int numsch, int maxlvl)
{
  int i, j, k=0, tot, score, ind, thisbin;
  float cands, totcands, minobs=200, thiscands, thismin, lscore;
  FDataPtr ldata, tdata, vldata;
  for (i=0; i < numsch; i++) {
    vldata=NULL;
    for (j=1; j <= maxlvl; j++) {
      tdata=data=first_data(j,datas);
      if (data == NULL) continue;
      ind=j*numsch*bins+i*bins;
      thisbin=thismin=cands=0;
      totcands=counts[ind];
      thismin=totcands/(2*bins);
      if (thismin < minobs) thismin=minobs;
      lscore=-1;
      ldata=NULL;
      while (tdata != vldata) {
	while (tdata != vldata && tdata->flam->ind!=ind) {
	  ldata=tdata;
	  tdata=tdata->flam->next;
	}
	if (tdata == vldata) continue;
	tdata=updata(tdata,ldata,bins,thismin,&thisbin,
		     &cands,&thiscands,&totcands,&lscore);
      }
      vldata=data;
    }
  }
}

int bucket_counts(TDTrPtr TC)
{
  int i, maxlvl, bins=BINS, numsch=TC->parser->features->numschema, ind;
  float *counts;
  FHashPtr fhash;
  FeatPtr feat;
  FDataPtr data, *datas, *lowdatas;
  if (TC->verbose) fprintf(stderr,"Grammar: bins  "); 
  maxlvl=set_bscore(TC);
  ind=(maxlvl+1)*numsch*bins;
  counts=malloc(ind*sizeof(float));
  datas=malloc((maxlvl+1)*MAXLABLEN*sizeof(FDataPtr));
  lowdatas=malloc((maxlvl+1)*MAXLABLEN*sizeof(FDataPtr));
  for (i=0; i < ind; i++) counts[i]=0.0;
  for (i=0; i < (maxlvl+1)*MAXLABLEN; i++) datas[i]=lowdatas[i]=NULL;
  data=order_bins(TC,datas,lowdatas,counts,bins,numsch,maxlvl);
  assign_bins(data,datas,lowdatas,counts,bins,numsch,maxlvl);
  free(datas);
  free(lowdatas);
  free(counts);
  return ind;
}

LamPtr init_lambda(int xvals)
{
  int i;
  LamPtr lambda=malloc(sizeof(struct LambdaStruct));
  lambda->olambdas=malloc(xvals*sizeof(float));
  lambda->lambdas=malloc(xvals*sizeof(float));
  lambda->batch=0;
  lambda->ldata=NULL;
  lambda->alldata=NULL;
  for (i=0; i < xvals; i++) lambda->lambdas[i]=lambda->olambdas[i]=0.49;
  return lambda;
}

void addnewlambda(FtrPtr feats, LamPtr *lambdas, int ind, FeatPtr feat, FDataPtr data, int j)
{
  LamDaPtr lambda=malloc(sizeof(struct LambData));
  lambda->feat=feat;
  lambda->cfeat=data->cfeat[j];
  lambda->RootExp=lambda->AddExp=lambda->phat=lambda->psmoo=0;
  lambda->next=lambda->backoff=NULL;
  lambda->anext=lambdas[ind]->alldata;
  lambdas[ind]->alldata=lambda;
  feats->scores->clams[lambda->cfeat->idx]->lambda=lambda;
}

void init_lambdas(TDTrPtr TC, LamPtr *lambdas, int tot)
{
  FtrPtr feats=TC->parser->features;
  int i, j, numFH=feats->numFH, ind, bins=BINS, numsch=feats->numschema, hash;
  FHashPtr fhash;
  FHashPtr currhash;
  FeatPtr feat, bf;
  FDataPtr data;
  LamDaPtr lambda;
  for (i=0; i < numFH; i++) {
    fhash=TC->parser->features->FHash[i];
    while (fhash != NULL) {
      feat=fhash->feat;
      data=feat->data;
      fhash=fhash->next;
      if (data == NULL) continue;
      ind=feat->level*numsch*bins+feat->schema*bins+data->flam->bucket;
      if (lambdas[ind]==NULL) lambdas[ind]=init_lambda(TC->xvals);
      for (j=0; j < data->numc; j++) 
	addnewlambda(feats,lambdas,ind,feat,data,j);
     }
  }
  for (i=0; i < tot; i++) {
    if (lambdas[i] == NULL) continue;
    lambda=lambdas[i]->alldata;
    while (lambda != NULL) {
      bf=lambda->feat->data->prefix;
      if (bf != NULL) {
	currhash=find_feat(bf->dstate,lambda->cfeat->label,bf->level+1,
			   bf->schema+1,feats,&hash);
	if (currhash == NULL) failproc("no backoff arc!");
	lambda->backoff=feats->scores->clams[bf->data->cfeat[currhash->feat->dstate]->idx]->lambda;
      }
      lambda=lambda->anext;
    }
  }
}

void free_clams(ScorePtr scr)
{
  int i;
  for (i=0; i < scr->numalloc; i++) {
    if (scr->clams[i] == NULL) continue;
    if (scr->clams[i]->pscore != NULL) free(scr->clams[i]->pscore);
    free(scr->clams[i]);
  }
  free(scr->clams);
  scr->clams=NULL;
}

void free_lambdas(LamPtr *lambdas, FtrPtr feats, int ind)
{
  int i;
  LamDaPtr lambda, tlambda;
  for (i=0; i < ind; i++) {
    if (lambdas[i] == NULL) continue;
    free(lambdas[i]->lambdas);
    free(lambdas[i]->olambdas);
    lambda=lambdas[i]->alldata;
    while (lambda != NULL) {
      tlambda=lambda;
      lambda=lambda->next;
      free(tlambda);
    }
    free(lambdas[i]);
  }
  free(lambdas);
  free_clams(feats->scores);
}

float batchc(ScorePtr scr, CFDataPtr cfeat, int batch)
{
  int idx=cfeat->idx;
  float cnt=scr->scores[idx];
  if (scr->clams[idx]->pscore != NULL) cnt-=scr->clams[idx]->pscore[batch];
  else if (scr->clams[idx]->sbatch==batch) cnt=0;
  return cnt;
}

float get_expcnt(LamDaPtr lambda, float lamval)
{
  float ret=1-lamval;
  lambda->AddExp+=lambda->RootExp;
  ret*=lambda->AddExp;
  return ret;
}

float get_smprob(LamDaPtr lambda, float lamval)
{
  float olam=1-lamval, pr;
  if (lambda->backoff==NULL) return lambda->phat;
  pr=lamval*lambda->phat;
  olam *= lambda->backoff->psmoo;
  pr+=olam;
  return pr;
}

void set_psmooth(LamPtr *lambdas, int tot, int batch)
{
  int i, j, cnt=0, cnt0=0;
  float lam, ev, pr;
  LamDaPtr lambda, llambda, flambda;
  FDataPtr data;
  for (i=0; i < tot; i++) {
    if (lambdas[i] == NULL) continue;
    lambda=lambdas[i]->ldata;
    lam=lambdas[i]->lambdas[batch];
    flambda=llambda=NULL;
    while (lambda != NULL) {
      lambda->psmoo=get_smprob(lambda,lam);
      if (lambda->psmoo>0) {
	lambda->backoff->RootExp-=lambda->RootExp;
	if (flambda==NULL) flambda=lambda;
	llambda=lambda;
      }
      else if (llambda!=NULL) llambda->next=lambda->next;
      lambda=lambda->next;
    }
    lambdas[i]->ldata=flambda;
  }
}

void init_lcountX(TDTrPtr TC, int batch)
{
  FtrPtr feats=TC->parser->features;
  int i, j, numFH=feats->numFH;
  FHashPtr fhash;
  FDataPtr data;
  for (i=0; i < numFH; i++) {
    fhash=feats->FHash[i];
    while (fhash != NULL) {
      data=fhash->feat->data;
      fhash=fhash->next;
      if (data == NULL) continue;
      data->bscore=0.0;
      for (j=0; j<data->numc; j++) 
	data->bscore+=batchc(feats->scores,data->cfeat[j],batch);
    }
  }
}

void set_lambdaX(LamPtr *lambdas, int tot, int batch)
{
  int i;
  float lam, ev;
  LamDaPtr lambda;
  for (i=tot-1; i >= 0; i--) {
    if (lambdas[i] == NULL) continue;
    lambda=lambdas[i]->alldata;
    lambdas[i]->ldata=NULL;
    lam=lambdas[i]->lambdas[batch];
    while (lambda != NULL) {
      if (lambda->feat->data->bscore > 0) {
	ev=get_expcnt(lambda,lam);
	if (ev > 0 && lambda->backoff != NULL) {
	  lambda->backoff->AddExp+=ev;
	  if (lambdas[i]->ldata!=NULL) lambda->next=lambdas[i]->ldata;
	  lambdas[i]->ldata=lambda;
	}
      }
      lambda=lambda->anext;
    }
  }
  set_psmooth(lambdas,tot,batch);
}

void init_lambdaX(ScorePtr scr, LamPtr *lambdas, int tot, int batch)
{
  int i, j;
  LamDaPtr lambda;
  FDataPtr data;
  for (i=0; i < tot; i++) {
    if (lambdas[i] == NULL) continue;
    lambda=lambdas[i]->alldata;
    lambdas[i]->batch=batch;
    while (lambda != NULL) {
      data=lambda->feat->data;
      lambda->AddExp=lambda->psmoo=lambda->RootExp=lambda->phat=0;
      if (data->bscore > 0) {
	lambda->RootExp=scr->scores[lambda->cfeat->idx];
	lambda->phat=batchc(scr,lambda->cfeat,batch);
	lambda->RootExp-=lambda->phat;
	lambda->phat/=data->bscore;
	lambda->psmoo=lambda->phat;
      }
      lambda->next=NULL;
      lambda=lambda->anext;
    }
  }
  set_lambdaX(lambdas,tot,batch);
}

float update_exps(LamPtr *lambdas, int tot, int batch)
{
  int i;
  float lam, llike=0, tlike;
  LamDaPtr lambda;
  for (i=0; i < tot; i++) {
    if (lambdas[i] == NULL) continue;
    lambda=lambdas[i]->alldata;
    lam=lambdas[i]->olambdas[batch]=lambdas[i]->lambdas[batch];
    while (lambda != NULL) {
      lambda->AddExp=0;
      lambda->psmoo=get_smprob(lambda,lam);
      lambda=lambda->anext;
    }
  }
  for (i=tot-1; i >= 0; i--) {
    if (lambdas[i] == NULL) continue;
    lambda=lambdas[i]->ldata;
    lam=lambdas[i]->lambdas[batch];
    while (lambda != NULL) {
      lambda->backoff->AddExp+=get_expcnt(lambda,lam);
      if (lambda->RootExp > 0) {
	tlike=-log(lambda->psmoo);
	tlike*=lambda->RootExp;
	llike+=tlike;
      }
      lambda=lambda->next;
    }
  }
  return llike;
}

void max_likes(LamPtr *lambdas, int tot, int batch)
{
  int i;
  float lam, olam, den, eps=0.001, psm;
  LamDaPtr lambda;
  for (i=0; i < tot; i++) {
    if (lambdas[i] == NULL) continue;
    lambda=lambdas[i]->ldata;
    lam=lambdas[i]->lambdas[batch];
    olam=1-lam;
    lambdas[i]->lambdas[batch]=den=0;
    while (lambda != NULL) {
      psm=lam*lambda->phat;
      psm+=olam*lambda->backoff->psmoo;
      lambdas[i]->lambdas[batch]+=lambda->AddExp*lam*lambda->phat/psm;
      den+=lambda->AddExp*olam*lambda->backoff->psmoo/psm;
      lambda=lambda->next;
    }
    if (lambdas[i]->lambdas[batch] <= 0) lambdas[i]->lambdas[batch]=eps;
    den+=lambdas[i]->lambdas[batch];
    if (lambdas[i]->lambdas[batch]>=den) den=lambdas[i]->lambdas[batch]+eps;
    lambdas[i]->lambdas[batch]/=den;
  }
}

void max_lambdaX(LamPtr *lambdas, int tot, int batch)
{
  int i, j, iter=0;
  float llike=0, olike=0, odelta=1.0001, zlike;
  while (iter < 1000 && (llike*odelta < olike || iter < 2)) {
    olike=llike;
    max_likes(lambdas,tot,batch);
    llike=update_exps(lambdas,tot,batch);
    if (iter==0) zlike=llike;
    iter++;
  }
}

void avg_lambdas(LamPtr *lambdas, int xvals, int tot)
{
  int i, batch;
  float lam;
  LamDaPtr lambda;
  for (i=0; i < tot; i++) {
    if (lambdas[i] == NULL) continue;
    lam=0;
    for (batch=0; batch < xvals; batch++) lam+=lambdas[i]->lambdas[batch];
    lambdas[i]->lambdas[0]=lam;
    lambdas[i]->lambdas[0]/=xvals;
  }
}

void set_lambdas(TDTrPtr TC, LamPtr *lambdas)
{
  FtrPtr feats=TC->parser->features;
  int i, numFH=feats->numFH, ind, bins=BINS, numsch=feats->numschema;
  float def=0.5, score, eps=0.00001;
  FHashPtr fhash;
  FeatPtr feat;
  FDataPtr data;
  for (i=0; i < numFH; i++) {
    fhash=TC->parser->features->FHash[i];
    while (fhash != NULL) {
      feat=fhash->feat;
      data=feat->data;
      fhash=fhash->next;
      if (data == NULL) continue;
      ind=feat->level*numsch*bins+feat->schema*bins+data->flam->bucket;
      score=1;
      if (lambdas[ind]==NULL) score=def;
      else if (lambdas[ind]->lambdas[0] <= 0) score-=eps;
      else score-=lambdas[ind]->lambdas[0];
      if (score <= 0) score=eps;
      data->bscore=-log(score);
    }
  }
}

void em_smooth(TDTrPtr TC, int ind)
{
  int i, xvals=TC->xvals;
  LamPtr *lambdas=malloc(ind*sizeof(LamPtr));
  if (!TC->xlast) xvals--;
  if (TC->verbose) fprintf(stderr,"Grammar: EM %d-way X-val: ",xvals); 
  for (i=0; i < ind; i++) lambdas[i]=NULL;
  init_lambdas(TC,lambdas,ind);
  for (i=0; i < xvals; i++) {
    if (TC->verbose) fprintf(stderr," %d",i);
    init_lcountX(TC,i);
    init_lambdaX(TC->parser->features->scores,lambdas,ind,i);
    max_lambdaX(lambdas,ind,i);
  }
  if (TC->verbose) fprintf(stderr,"\n");
  avg_lambdas(lambdas,xvals,ind);
  set_lambdas(TC,lambdas);
}

void init_lrn(TDTrPtr TC, int bsize)
{
  int i, j;
  ScorePtr scr=TC->parser->features->scores;
  scr->lrn=malloc(sizeof(struct LearnScore));
  scr->lrn->iter=0;
  scr->lrn->xvals=TC->xvals;
  scr->lrn->bsize=bsize;
  scr->lrn->fscores=malloc(TC->xvals*sizeof(FScorePtr));
  scr->lrn->ascores=malloc(scr->numalloc*sizeof(float));
  scr->lrn->lastupd=malloc(scr->numalloc*sizeof(int));
  scr->lrn->pscores=scr->pscores;
  scr->lrn->scores=scr->scores;
  scr->lrn->ppsc=scr->ppsc;
  scr->lrn->apsc=malloc(scr->alcpsc*sizeof(float));
  scr->lrn->lastpsc=malloc(scr->alcpsc*sizeof(int));
  for (i=0; i < TC->xvals; i++) {
    scr->lrn->fscores[i]=malloc(sizeof(struct FltScore));
    scr->lrn->fscores[i]->rsc=malloc(scr->numalloc*sizeof(float));
  }
  for (i=0; i < scr->alcpsc; i++) {
    scr->lrn->apsc[i]=0.0;
    scr->lrn->lastpsc[i]=0;
  }
  for (i=0; i < scr->numalloc; i++) {
    scr->lrn->ascores[i]=scr->pscores[i]=0.0;
    scr->lrn->lastupd[i]=0;
    for (j=0; j < TC->xvals; j++) scr->lrn->fscores[j]->rsc[i]=0.0;
  }
}

float *init_percep(TDTrPtr TC, int bsize)
{
  int i, j;
  float *den=malloc(TC->xvals*sizeof(float));
  ScorePtr scr=TC->parser->features->scores;
  scr->pscores=malloc(scr->numalloc*sizeof(float));
  scr->alcpsc=scr->numalloc;
  scr->ppsc=malloc(scr->numalloc*sizeof(float));
  for (i=0; i < scr->numalloc; i++) scr->ppsc[i]=0;
  init_lrn(TC,bsize);
  return den;
}

void fill_den(float *den, FDataPtr data, TDTrPtr TC)
{
  int i,j;
  ScorePtr scr=TC->parser->features->scores;
  CLamPtr clam;
  for (i=0; i < TC->xvals; i++) den[i]=data->score;
  for (j=0; j < data->numc; j++) {
    clam=scr->clams[data->cfeat[j]->idx];
    if (clam->pscore==NULL) 
      den[clam->sbatch]-=scr->scores[data->cfeat[j]->idx];
    else for (i=0; i < TC->xvals; i++) 
      if (clam->pscore[i]>0) den[i]-=clam->pscore[i];
  }
}

void fill_fsc(float *den, CFDataPtr cfeat, TDTrPtr TC)
{
  int i;
  float num;
  ScorePtr scr=TC->parser->features->scores;
  CLamPtr clam;
  for (i=0; i < TC->xvals; i++) {
    num=scr->scores[cfeat->idx];
    clam=scr->clams[cfeat->idx];
    if (clam->pscore==NULL) {
      if (clam->sbatch==i) num=0.0;
    }
    else num-=clam->pscore[i];
    if (num > 0) num/=den[i];
    scr->lrn->fscores[i]->rsc[cfeat->idx] = -log(num);
  }
}

int norm_feats(TDTrPtr TC, int ind, int sents)
{
  int i, j, k=0, m, numFH=TC->parser->features->numFH;
  float num, *den=NULL;
  ScorePtr scr=TC->parser->features->scores;
  FHashPtr fhash;
  FDataPtr data;
  em_smooth(TC,ind);
  den=init_percep(TC,sents);
  if (TC->verbose) fprintf(stderr,"Grammar: norm. "); 
  for (i=0; i < numFH; i++) {
    fhash=TC->parser->features->FHash[i];
    while (fhash != NULL) {
      data=fhash->feat->data;
      if (data != NULL) {
	if (den != NULL) fill_den(den,data,TC);
	for (j=0; j < data->numc; j++) {
	  if (den!=NULL) fill_fsc(den,data->cfeat[j],TC);
	  scr->scores[data->cfeat[j]->idx] = 
	    -log(scr->scores[data->cfeat[j]->idx]/data->score);
	}
	verbose_count(TC->verbose,++k,1000000);
      }
      fhash=fhash->next;
    }
  }
  if (TC->verbose) fprintf(stderr," %d\n",k);
  if (den != NULL) free(den);
  return presum_feats(TC);
}

void count_backoffs(TDTrPtr TC, CFDataPtr *cfeats)
{
  FtrPtr feats=TC->parser->features;
  int i, j, numFH=feats->numFH;
  FHashPtr fhash;
  FeatPtr feat;
  FDataPtr data;
  for (i=0; i < numFH; i++) {
    fhash=TC->parser->features->FHash[i];
    while (fhash != NULL) {
      feat=fhash->feat;
      data=feat->data;
      fhash=fhash->next;
      if (data != NULL) {
	data->flam->ind=0;
	for (j=0; j < data->numc; j++) 
	  cfeats[data->cfeat[j]->idx]=data->cfeat[j];
      }
    }
  }
  for (i=0; i < numFH; i++) {
    fhash=TC->parser->features->FHash[i];
    while (fhash != NULL) {
      feat=fhash->feat;
      data=feat->data;
      fhash=fhash->next;
      if (data == NULL || data->prefix == NULL) continue;
      data->prefix->data->flam->ind++;
    }
  }  
}

void remscores(CFDataPtr cfeat, ScorePtr scr, CFDataPtr *cfeats, int xvals)
{
  int i;
  CFDataPtr redir=cfeats[scr->numind-1];
  CLamPtr clam;
  scr->numind--;
  if (redir==cfeat) return;
  if (scr->clams != NULL) {
    clam=scr->clams[cfeat->idx];  
    scr->clams[cfeat->idx]=scr->clams[redir->idx];
    scr->clams[redir->idx]=clam;
  }
  scr->scores[cfeat->idx]=scr->scores[redir->idx];
  if (scr->pscores != NULL) scr->pscores[cfeat->idx]=scr->pscores[redir->idx];
  if (scr->lrn != NULL) {
    scr->lrn->scores[cfeat->idx]=scr->lrn->scores[redir->idx];
    scr->lrn->ascores[cfeat->idx]=scr->lrn->ascores[redir->idx];
    scr->lrn->pscores[cfeat->idx]=scr->lrn->pscores[redir->idx];
    scr->lrn->lastupd[cfeat->idx]=scr->lrn->lastupd[redir->idx];
    for (i=0; i < xvals; i++) 
      scr->lrn->fscores[i]->rsc[cfeat->idx]=
	scr->lrn->fscores[i]->rsc[redir->idx];
  }
  redir->idx=cfeat->idx;
  cfeats[scr->numind]=NULL;
  cfeats[redir->idx]=redir;
}

void zapcfeat(FtrPtr feats,int j,FeatPtr feat,CFDataPtr *cfeats,int xvals)
{
  int hash, lb, st=feat->dstate, lv=feat->level+1, sch=feat->schema+1;
  FHashPtr fhash, lhash;
  remscores(feat->data->cfeat[j],feats->scores,cfeats,xvals);
  lb=feat->data->cfeat[j]->label;
  hash=get_fhash(st,lb,lv,sch,feats);
  fhash=feats->FHash[hash];
  lhash=NULL;
  while (fhash != NULL) {
    if (fhash->feat->state == st && fhash->feat->label == lb &&
	fhash->feat->level == lv && fhash->feat->schema == sch) break;
    lhash=fhash;
    fhash=fhash->next;
  }
  if (fhash==NULL) {
    fprintf(stderr,"couldn't find entry %d %d %d %d",hash,j,lb,st);
    exit(1);
  }
  zapentry(feats,hash,lhash,fhash,cfeats,xvals);
}

void updzcfeat(FtrPtr feats,int j,FeatPtr feat,CFDataPtr *cfeats)
{
  int hash, lb, st=feat->dstate, lv=feat->level+1, sch=feat->schema+1;
  FHashPtr fhash, lhash;
  lb=feat->data->cfeat[j]->label;
  hash=get_fhash(st,lb,lv,sch,feats);
  fhash=feats->FHash[hash];
  while (fhash != NULL) {
    if (fhash->feat->state == st && fhash->feat->label == lb &&
	fhash->feat->level == lv && fhash->feat->schema == sch) break;
    fhash=fhash->next;
  }
  if (fhash==NULL) {
    fprintf(stderr,"couldn't find entry %d %d %d %d",hash,j,lb,st);
    exit(1);
  }
  fhash->feat->dstate=j;
}

void zapentry(FtrPtr feats, int i, FHashPtr last, FHashPtr this, CFDataPtr *cfeats, int xvals)
{
  int j, numc;
  FeatPtr feat=this->feat;
  if (last==NULL) feats->FHash[i] = this->next;
  else last->next = this->next;
  if (feat->data == NULL) return;
  numc=feat->data->numc;
  for (j=0; j < numc; j++) zapcfeat(feats,j,feat,cfeats,xvals);
}

void prune_zeros(TDTrPtr TC, CFDataPtr *cfeats, int k)
{
  int i, j;
  FtrPtr feats=TC->parser->features;
  FHashPtr fhash, last;
  FDataPtr data;
  if (TC->verbose) fprintf(stderr," %d  Pruning: ",k);
  for (i=0; i < feats->numFH; i++) {
    last=NULL;
    fhash=TC->parser->features->FHash[i];
    while (fhash != NULL) {
      data=fhash->feat->data;
      if (data != NULL && data->score == 0.0) {
	zapentry(TC->parser->features,i,last,fhash,cfeats,TC->xvals);
	verbose_count(TC->verbose,--k,1000000);
      }
      else last=fhash;
      fhash=fhash->next;
    }
  }
  if (TC->verbose) fprintf(stderr," %d\n",k);
}

float getdiff(FtrPtr feats, FeatPtr bfeat, float bcost, CFDataPtr cfeat, int st,int lv, int sch)
{
  int hash;
  float diff=exp(-feats->scores->scores[cfeat->idx]),bsc;
  FHashPtr bhash=find_feat(st,cfeat->label,lv,sch,feats,&hash);
  if (bhash == NULL) return 0.0;
  bsc=getcost(feats->scores,bfeat,bhash,bcost,NULL);
  diff-=exp(-bsc);
  diff/=exp(-bsc);
  return diff;
}

void prune_lowdiff(TDTrPtr TC, CFDataPtr *cfeats, int m)
{
  int i, j, st, lv, sch, hash;
  float bsc, diff;
  FtrPtr feats=TC->parser->features;
  ScorePtr scr=feats->scores;
  FDataPtr data;
  FeatPtr feat;
  FHashPtr fhash, last, bhash;
  for (i=0; i < feats->numFH; i++) {
    fhash=TC->parser->features->FHash[i];
    while (fhash != NULL) {
      feat=fhash->feat;
      fhash=fhash->next;
      data=feat->data;
      if (data == NULL || data->prefix==NULL || feat->level != m) continue;
      st=data->prefix->dstate;
      lv=data->prefix->level+1;
      sch=data->prefix->schema+1;
      for (j=0; j < data->numc; j++) {
	diff=getdiff(feats,data->prefix,data->bscore,data->cfeat[j],st,lv,sch);
	if (diff >= TC->thresh) continue;
	zapcfeat(feats,j,feat,cfeats,TC->xvals);
	data->numc--;
	if (j == data->numc) continue;
	data->cfeat[j]=data->cfeat[data->numc];
	updzcfeat(feats,j,feat,cfeats);
	j--;
      }
      if (data->numc==0) {
	free(data->cfeat);
	data->cfeat=NULL;
      }
    }
  }
}

int prunefeat(FDataPtr data, TDTrPtr TC)
{
  FeatPtr prefix=data->prefix;
  if (data->bscore < TC->thresh) return 1;
  if (data->numc==1 && data->numc==prefix->data->numc) return 1;
  return 0;
}

void prune_feats(TDTrPtr TC, int max)
{
  int i, j, k=0, m;
  FtrPtr feats=TC->parser->features;
  CFDataPtr *cfeats=malloc(feats->scores->numind*sizeof(CFDataPtr));
  FHashPtr fhash, last;
  FDataPtr data;
  for (i=0; i < feats->scores->numind; i++) cfeats[i]=NULL;
  count_backoffs(TC,cfeats);
  if (TC->verbose) fprintf(stderr,"Grammar: prune "); 
  for (m=max-1; m >=0; m--) {
    for (i=0; i < feats->numFH; i++) {
      fhash=TC->parser->features->FHash[i];
      while (fhash != NULL) {
	data=fhash->feat->data;
	if (data != NULL && fhash->feat->level == m && data->prefix != NULL && 
	    data->flam->ind == 0 && prunefeat(data,TC)) {
	  data->prefix->data->flam->ind--;
	  data->score=0.0;
	  verbose_count(TC->verbose,++k,1000000);
	}
	fhash=fhash->next;
      }
    }
  }
  prune_zeros(TC,cfeats,k);
  for (m=1; m < max; m++) prune_lowdiff(TC,cfeats,m);
  free(cfeats);
}

void update_wpos(ParPtr par)
{
  int i, j, k, hash;
  LexPtr lex=par->lexicon;
  FtrPtr feats=par->features;
  ScorePtr scr=feats->scores;
  FHashPtr fhash;
  FDataPtr fdata;
  for (i=0; i < feats->numschema; i++) {
    if (feats->schema[i]->functs[0]->moves[0] != 'W') continue;
    for (j=1; j < lex->numPTs; j++) {
      fhash=find_feat(0,j,1,i,feats,&hash);
      if (fhash == NULL) continue;
      fdata=fhash->feat->data;
      for (k=0; k < fdata->numc; k++) {
	lex->WDs[fdata->cfeat[k]->label]->NTs[j]=
	  scr->scores[fdata->cfeat[k]->idx];
	if (scr->pscores!=NULL) lex->WDs[fdata->cfeat[k]->label]->NTs[j]+=
				  scr->pscores[fdata->cfeat[k]->idx];
      }
    }
    break;
  }
  if (i==feats->numschema) failproc("no schema predicting words!");
}

void update_ntpos(ParPtr par)
{
  int i, j;
  double sum;
  LexPtr lex=par->lexicon;
  for (i=lex->numPTs; i < lex->numNTs; i++) {
    if (lex->WDs[i]->slc != NULL && !lex->WDs[i]->slc->orig) continue;
    sum=0;
    for (j=1; j < lex->numPTs; j++) sum += lex->WDs[i]->NTs[j];
    for (j=1; j < lex->numPTs; j++)
      lex->WDs[i]->NTs[j] = -log(lex->WDs[i]->NTs[j]/sum);
  }
}

void upPTs(LexPtr lex, double *PTs, int lbl, float prob)
{
  int m;
  float lz=-log(0);
  for (m=1; m < lex->numPTs; m++) 
    if (lex->WDs[lbl]->NTs[m] < lz) PTs[m] += prob*exp(-lex->WDs[lbl]->NTs[m]);
}

void update_lcpos(ParPtr par)
{
  int i, j, k, m, hash;
  double *PTs, prob;
  LexPtr lex=par->lexicon;
  FtrPtr feats=par->features;
  ScorePtr scr=feats->scores;
  FHashPtr fhash;
  FDataPtr fdata;
  PTs = malloc(lex->numPTs*sizeof(double));
  for (i=0; i < feats->numschema; i++) {
    if (feats->schema[i]->functs[0]->moves[0] != 'N') continue;
    for (j=lex->numPTs; j < lex->numNTs; j++) {
      if (lex->WDs[j]->slc == NULL || lex->WDs[j]->slc->orig) continue;
      for (k=1; k < lex->numPTs; k++) PTs[k]=0;
      fhash=find_feat(0,j,1,i,feats,&hash);
      if (fhash != NULL) {
	fhash=find_feat(fhash->feat->dstate,-1,2,i,feats,&hash);
	if (fhash == NULL) failproc("no observation of leftmost from slc");
	fdata=fhash->feat->data;
	for (k=0; k < fdata->numc; k++) {
	  if (fdata->cfeat[k]->label == 0) continue;
	  prob=exp(-scr->scores[fdata->cfeat[k]->idx]);
	  if (fdata->cfeat[k]->label < lex->numPTs) 
	    PTs[fdata->cfeat[k]->label] += prob;
	  else upPTs(lex,PTs,fdata->cfeat[k]->label,prob);
	}
      }
      for (k=1; k < lex->numPTs; k++) lex->WDs[j]->NTs[k] = -log(PTs[k]);
    }
    break;
  }
  if (i==feats->numschema) failproc("no schema predicting non-terms!");
  free(PTs);
}

void update_pos(ParPtr par)
{
  update_wpos(par);
  update_ntpos(par);
  update_lcpos(par);
}

float getcost(ScorePtr scr, FeatPtr feat, FHashPtr fhash, float bcost, float *p)
{
  int idx=feat->data->cfeat[fhash->feat->dstate]->idx;
  float cost=bcost+scr->scores[idx];
  if (p!=NULL && scr->pscores!=NULL) p[0]=scr->pscores[idx];
  return cost;
}

float getboff(FtrPtr feats, FeatPtr feat, int lbl, float bcost, int bottom, float *p)
{
  int hash,idx;
  FeatPtr bfeat;
  ScorePtr scr=feats->scores;
  FHashPtr fhash;
  float z=0, cost=-log(z), nbcost=bcost, boff;
  if (feat == NULL) return cost;
  nbcost+=feat->data->bscore;
  bfeat=feat->data->prefix;
  if (bfeat==NULL && !bottom) return cost;
  fhash=find_feat(feat->dstate,lbl,feat->level+1,feat->schema+1,feats,&hash);
  if (fhash == NULL) return getboff(feats,bfeat, lbl, nbcost, bottom, p);
  return getcost(scr,feat,fhash,bcost,p);
}

void get_ntprobs(int lb, int siblb, TDParPtr PC)
{
  int i, j, hash;
  float score;
  FeatPtr bf;
  LexPtr lex=PC->parser->lexicon;
  FtrPtr feats=PC->parser->features;
  ScorePtr scr=feats->scores;
  FHashPtr fhash;
  FDataPtr fdata;
  for (i=0; i < lex->numNTs; i++) lex->nts[i]=0;
  for (i=0; i < feats->numschema; i++) {
    if (feats->schema[i]->functs[0]->moves[0] != 'N') continue;
    fhash=find_feat(0,lb,1,i,feats,&hash);
    if (fhash == NULL) break;
    fhash=find_feat(fhash->feat->dstate,siblb,2,i,feats,&hash);
    if (fhash == NULL) break;
    fdata=fhash->feat->data;
    bf=fdata->prefix;
    for (j=0; j < fdata->numc; j++) 
      lex->nts[fdata->cfeat[j]->label]=exp(-scr->scores[fdata->cfeat[j]->idx]);
    break;
  }
}

void fill_parser(ParPtr parser)
{
  int i;
  parser->lexicon->nts=malloc(parser->lexicon->numNTs*sizeof(double));
  parser->ntree=default_tree(parser->lexicon->WDs[1],NULL);
  parser->stree=default_tree(parser->lexicon->WDs[0],NULL);
  parser->fd=malloc(parser->lexicon->numNTs*sizeof(int));
  parser->exps=malloc(parser->lexicon->numNTs*sizeof(float));
  parser->pexps=malloc(parser->lexicon->numNTs*sizeof(float));
  parser->sch=parser->wsch=-1;
  for (i=0; i < parser->features->numschema; i++) {
    if (parser->sch < 0 && parser->features->schema[i]->functs[0]->moves[0] == 'N') 
      parser->sch=i;
    if (parser->wsch < 0 && parser->features->schema[i]->functs[0]->moves[0] == 'W') 
      parser->wsch=i;
  }
  for (i = 0; i < parser->lexicon->numWDs; i++)
    if (parser->lexicon->WDs[i]->root) {
      parser->lexicon->root=i;
      break;
    }
  if (parser->sch<0) failproc("no N emitting schema");
  if (parser->wsch<0) failproc("no W emitting schema");
}

ParPtr init_parser(FILE *fp, int lrn)
{
  ParPtr parser=read_parser(fp,lrn);
  fill_parser(parser);
  return parser;
}

void fill_borev(FtrPtr feats)
{
  int i, j, st, lv, sch, hash, idx;
  FHashPtr fhash;
  FeatPtr feat, bfeat;
  BORPtr bor;
  for (i=0; i < feats->numFH; i++) {
    fhash=feats->FHash[i];
    while (fhash != NULL) {
      feat=fhash->feat;
      fhash=fhash->next;
      if (feat->data==NULL) continue;
      bfeat=feat->data->prefix;
      if (bfeat==NULL) continue;
      st=bfeat->dstate;
      lv=bfeat->level+1;
      sch=bfeat->schema+1;
      for (j=0; j < feat->data->numc; j++) {
	fhash=find_feat(st,feat->data->cfeat[j]->label,lv,sch,feats,&hash);
	if (fhash==NULL) {
	  for (j=0; j < feat->data->numc; j++)
	    fprintf(stderr,"%d ",feat->data->cfeat[j]->label);
	  fprintf(stderr," %d\n",feat->data->numc);
	  for (j=0; j < bfeat->data->numc; j++)
	    fprintf(stderr,"%d ",bfeat->data->cfeat[j]->label);
	  fprintf(stderr," %d\n",bfeat->data->numc);
	  failproc("no backoff label");
	}
	idx=bfeat->data->cfeat[fhash->feat->dstate]->idx;
	bor=malloc(sizeof(struct BORev));
	bor->next=feats->scores->lrn->backs[idx];
	bor->state=feat->data->cfeat[j]->idx;
	feats->scores->lrn->backs[idx]=bor;
      }
    }
  }
}
