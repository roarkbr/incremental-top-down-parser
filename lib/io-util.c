// io-util.c
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
#include "parser.h"
#include "io-util.h"

FILE *gopenfile(char *file, char *mode)
{
  FILE *fp;
  char filename[MAXLABLEN];
  if ( mode[0] == 'r' ) {
    if (file != NULL) sprintf(filename,"gzip -cd %s",file);
    else sprintf(filename,"gzip -cd");
  }
  else {
    if (file != NULL) sprintf(filename,"gzip -c >%s",file);
    else sprintf(filename,"gzip -c");
  }
  if ((fp = popen(filename,mode)) == NULL) {
    fprintf(stderr,"could not open: %s\n",filename);
    exit(1);
  }
  return fp;
}

void byteswp(void *out, void *in, size_t l)
{
  size_t i;
  char c, *inp = (char *) in, *outp = (char *) out;
  for(i=0; i<l/2; i++) {
    c = inp[i];
    outp[i] = inp[l-i-1];
    outp[l-i-1] = c;
  }
}

int rdsw_int(FILE *fp, int swp)
{
  int i, j;
  if (fread(&i, sizeof(int), 1, fp) < 1) {
    fprintf(stderr,"integer read failed");
    exit(1);
  }
  if (swp) byteswp(&j,&i,sizeof(int));
  else j = i;
  return j;
}

float rdsw_float(FILE *fp, int swp)
{
  float i, j;
  if (fread(&i, sizeof(float), 1, fp) < 1) {
    fprintf(stderr,"float read failed");
    exit(1);
  }
  if (swp) byteswp(&j,&i,sizeof(float));
  else j = i;
  return j;
}

double rdsw_dbl(FILE *fp, int swp)
{
  double i, j;
  if (fread(&i, sizeof(double), 1, fp) < 1) {
    fprintf(stderr,"double read failed");
    exit(1);
  }
  if (swp) byteswp(&j,&i,sizeof(double));
  else j = i;
  return j;
}

void write_bslc(FILE *fp,SLCPtr slc)
{
  int o=1, z=0;
  while (slc != NULL) {
    fwrite(&(o), sizeof(int), 1, fp);
    fwrite(&(slc->LChild->windex), sizeof(int), 1, fp);
    fwrite(&(slc->NewCat->windex), sizeof(int), 1, fp);
    fwrite(&(slc->orig), sizeof(int), 1, fp);
    slc=slc->next;
  }
  fwrite(&(z), sizeof(int), 1, fp);  
}

SLCPtr read_bslc(FILE *fp,LexPtr lex,int swp)
{
  int ind, next;
  SLCPtr slc=NULL, last;
  next=rdsw_int(fp,swp);
  while (next) {
    last=slc;
    slc=malloc(sizeof(struct SLCPairs));
    slc->next=last;
    ind=rdsw_int(fp,swp);
    slc->LChild = lex->WDs[ind];
    ind=rdsw_int(fp,swp);
    slc->NewCat = lex->WDs[ind];
    slc->orig=rdsw_int(fp,swp);
    next=rdsw_int(fp,swp);
  }
  return slc;
}

void write_blexicon(FILE *fp, LexPtr lex)
{
  int i, j, len, unk;
  fwrite(&(lex->numNTs), sizeof(int), 1, fp);
  fwrite(&(lex->numPTs), sizeof(int), 1, fp);
  for (i = 0; i < lex->numWDs; i++) {
    len=strlen(lex->WDs[i]->label);
    fwrite(&i, sizeof(int), 1, fp);
    fwrite(&len, sizeof(int), 1, fp);
    fwrite(lex->WDs[i]->label, sizeof(char), len, fp);
    fwrite(&(lex->WDs[i]->children), sizeof(int), 1, fp);
    fwrite(&(lex->WDs[i]->terminal), sizeof(int), 1, fp);
    fwrite(&(lex->WDs[i]->sym), sizeof(int), 1, fp);
    fwrite(&(lex->WDs[i]->caps), sizeof(int), 1, fp);
    fwrite(&(lex->WDs[i]->root), sizeof(int), 1, fp);
    fwrite(&(lex->WDs[i]->oov), sizeof(int), 1, fp);
    fwrite(&(lex->WDs[i]->tot), sizeof(float), 1, fp);
    if (lex->WDs[i]->unkcat==NULL) unk=-1;
    else unk=lex->WDs[i]->unkcat->windex;
    fwrite(&(unk), sizeof(int), 1, fp);    
    if (i >= lex->numPTs)
      for (j=1; j < lex->numPTs; j++)
	fwrite(&(lex->WDs[i]->NTs[j]), sizeof(float), 1, fp);
  }
  for (i = lex->numPTs; i < lex->numNTs; i++) write_bslc(fp,lex->WDs[i]->slc);
}

void read_blexicon(FILE *fp,LexPtr lex,int swp)
{
  int i, j, len, hash, ch, t, unk;
  WLPtr word;
  char readstr[MAXLABLEN];
  lex->numNTs=rdsw_int(fp,swp);
  lex->numPTs=rdsw_int(fp,swp);
  for (i = 0; i < lex->allocWDs; i++) {
    j=rdsw_int(fp,swp);
    if (i != j) failproc("file format incorrect: lex");
    len=rdsw_int(fp,swp);
    fread(&readstr[0], sizeof(char), len, fp);
    readstr[len]=0;
    hash=get_whash(&readstr[0],lex);
    ch=rdsw_int(fp,swp);
    t=rdsw_int(fp,swp);
    word=insert_word(&readstr[0],lex,hash,ch,t);
    if (word->windex != i) failproc("wrong index!");
    word->sym=rdsw_int(fp,swp);
    word->caps=rdsw_int(fp,swp);
    word->root=rdsw_int(fp,swp);
    word->oov=rdsw_int(fp,swp);
    word->tot=rdsw_float(fp,swp);
    unk=rdsw_int(fp,swp);
    if (unk>=0) word->unkcat=lex->WDs[unk];
    if (i >= lex->numPTs) {
      word->NTs = malloc(lex->numPTs*sizeof(float));
      for (j=1; j < lex->numPTs; j++) word->NTs[j]=rdsw_float(fp,swp);
    }
  }
  for (i=lex->numPTs; i < lex->numNTs; i++) 
    lex->WDs[i]->slc=read_bslc(fp,lex,swp);
}

void write_bfsch(FILE *fp,FtrPtr feats)
{
  int i, j, no=-1;
  TWPtr ffunct;
  for (i=0; i < feats->numfeats; i++) {
    ffunct=feats->featfuncts[i];
    fwrite(&(i), sizeof(int), 1, fp);
    fwrite(&(ffunct->numoves), sizeof(int), 1, fp);
    fwrite(ffunct->moves, sizeof(char), ffunct->numoves, fp);
    fwrite(&(ffunct->numqueries), sizeof(int), 1, fp);
    fwrite(&(ffunct->SQuery), sizeof(int), 1, fp);
    for (j=0; j < ffunct->numqueries; j++) {
      if (ffunct->WQueries[j]==NULL) fwrite(&(no), sizeof(int), 1, fp);
      else fwrite(&(ffunct->WQueries[j]->windex), sizeof(int), 1, fp);
    }
    if (ffunct->suffmove == NULL) fwrite(&(no), sizeof(int), 1, fp);
    else fwrite(&(ffunct->suffmove->findex), sizeof(int), 1, fp);
   }
  for (i=0; i < feats->numschema; i++) {
    fwrite(&(i), sizeof(int), 1, fp);
    fwrite(&(feats->schema[i]->numfuncts), sizeof(int), 1, fp);
    for (j=0; j < feats->schema[i]->numfuncts; j++)
      fwrite(&(feats->schema[i]->functs[j]->findex), sizeof(int), 1, fp);      
  }
}

void write_bfback(FILE *fp, FtrPtr feats, int pcount)
{
  int i, j=0, no=-1, o=1, z=0;
  FHashPtr fhash;
  FeatPtr feat;
  fwrite(&(pcount), sizeof(int), 1, fp);     
  for (i=0; i < feats->numFH; i++) {
    fhash=feats->FHash[i];
    while (fhash != NULL) {
      feat=fhash->feat;
      fhash=fhash->next;
      if (feat->data != NULL && feat->data->prefix != NULL) {
	fwrite(&(feat->state), sizeof(int), 1, fp);
	fwrite(&(feat->label), sizeof(int), 1, fp);
	fwrite(&(feat->level), sizeof(int), 1, fp);
	fwrite(&(feat->schema), sizeof(int), 1, fp);
	if (feat->level != feat->data->prefix->level+1) 
	  failproc("levels don't match");
	if (feat->schema != feat->data->prefix->schema) 
	  failproc("schemas don't match");
	/*
	if (feat->data->prefix->data->score==0) 
	  failproc("writing pruned backoff");
	*/
	fwrite(&(feat->data->prefix->state), sizeof(int), 1, fp);	
	fwrite(&(feat->data->prefix->label), sizeof(int), 1, fp);	
	j++;
      }
    }
  }
  if (j != pcount) failproc("no match in pcount");
  fwrite(&(pcount), sizeof(int), 1, fp);
}

int write_bfmain(FILE *fp, FtrPtr feats)
{
  int i, j, no=-1, o=1, z=0, pcount=0, idx;
  FHashPtr fhash;
  FeatPtr feat;
  for (i=0; i < feats->numFH; i++) {
    fhash=feats->FHash[i];
    while (fhash != NULL) {
      fwrite(&(o), sizeof(int), 1, fp);     
      feat=fhash->feat;
      fhash=fhash->next;
      fwrite(&(feat->state), sizeof(int), 1, fp);
      fwrite(&(feat->label), sizeof(int), 1, fp);
      fwrite(&(feat->level), sizeof(int), 1, fp);
      fwrite(&(feat->schema), sizeof(int), 1, fp);
      fwrite(&(feat->dstate), sizeof(int), 1, fp);
      if (feat->data == NULL) {
	fwrite(&(z), sizeof(int), 1, fp);
	continue;
      }
      fwrite(&(o), sizeof(int), 1, fp);
      if (feat->data->prefix != NULL) pcount++;
      fwrite(&(feat->data->score), sizeof(float), 1, fp);
      fwrite(&(feat->data->bscore), sizeof(float), 1, fp);
      fwrite(&(feat->data->numc), sizeof(int), 1, fp);
      for (j=0; j < feat->data->numc; j++) {
	fwrite(&(feat->data->cfeat[j]->label), sizeof(int), 1, fp);
	fwrite(&(feat->data->cfeat[j]->idx), sizeof(int), 1, fp);
      }
    }
  }
  fwrite(&(z), sizeof(int), 1, fp);     
  return pcount;
}

void write_bfscores(FILE *fp, ScorePtr scores)
{
  int p=0, i;
  fwrite(&(scores->numind), sizeof(int), 1, fp);
  if (scores->pscores != NULL) p=1;
  fwrite(&(p), sizeof(int), 1, fp);  
  for (i=0; i < scores->numind; i++) {
    fwrite(&(scores->scores[i]), sizeof(float), 1, fp);    
    if (p) fwrite(&(scores->pscores[i]), sizeof(float), 1, fp);
  }
  fwrite(&(scores->numind), sizeof(int), 1, fp);
  fwrite(&(scores->alcpsc), sizeof(int), 1, fp);
  fwrite(&(scores->numpsc), sizeof(int), 1, fp);
  for (i=0; i < scores->numpsc; i++)
    fwrite(&(scores->ppsc[i]), sizeof(float), 1, fp);
}

void write_bfeatures(FILE *fp,FtrPtr feats)
{
  int i, j, no=-1, o=1, z=0, pcount=0;
  FHashPtr fhash;
  FeatPtr feat;
  fwrite(&(feats->numschema), sizeof(int), 1, fp);
  fwrite(&(feats->numfeats), sizeof(int), 1, fp);
  fwrite(&(feats->fstates), sizeof(int), 1, fp);
  write_bfsch(fp,feats);
  pcount=write_bfmain(fp,feats);
  write_bfback(fp,feats,pcount);
  write_bfscores(fp,feats->scores);
}

void write_blrn(FILE *fp,ScorePtr scr)
{
  int i, j, no=-1, o=1, z=0, pcount=0;
  fwrite(&(scr->lrn->xvals), sizeof(int), 1, fp);
  fwrite(&(scr->lrn->bsize), sizeof(int), 1, fp);
  fwrite(&(scr->lrn->iter), sizeof(int), 1, fp);
  for (i=0; i < scr->numind; i++) {
    fwrite(&(scr->lrn->pscores[i]), sizeof(float), 1, fp);    
    fwrite(&(scr->lrn->lastupd[i]), sizeof(int), 1, fp);    
    for (j=0; j < scr->lrn->xvals; j++)
      fwrite(&(scr->lrn->fscores[j]->rsc[i]), sizeof(float), 1, fp);    
  }
  fwrite(&(scr->numind), sizeof(int), 1, fp);
  for (i=0; i < scr->numpsc; i++) {
    fwrite(&(scr->lrn->pscores[i]), sizeof(float), 1, fp);    
    fwrite(&(scr->lrn->lastupd[i]), sizeof(int), 1, fp);    
  }
}

void read_bfsch(FILE *fp, FtrPtr feats, LexPtr lex, int swp)
{
  int i, j, ind, *suff=malloc(feats->numfeats*sizeof(int));
  TWPtr ffunct;
  for (i=0; i < feats->numfeats; i++) {
    ffunct=feats->featfuncts[i] = malloc(sizeof(struct TreeWalk));
    j=ffunct->findex=rdsw_int(fp,swp);
    if (i != j) failproc("file format incorrect: featfuncts");
    ffunct->numoves=rdsw_int(fp,swp);
    ffunct->moves=malloc(ffunct->numoves*sizeof(char));
    fread(ffunct->moves, sizeof(char), ffunct->numoves, fp);
    ffunct->numqueries=rdsw_int(fp,swp);
    ffunct->SQuery=rdsw_int(fp,swp);
    if (ffunct->numqueries==0) ffunct->WQueries=NULL;
    else ffunct->WQueries=malloc(ffunct->numqueries*sizeof(WLPtr));
    for (j=0; j < ffunct->numqueries; j++) {
      ind=rdsw_int(fp,swp);
      if (ind < 0) ffunct->WQueries[j]=NULL;
      else ffunct->WQueries[j]=lex->WDs[ind];
    }    
    suff[i]=rdsw_int(fp,swp);
    ffunct->suffmove=NULL;
  }
  for (i=0; i < feats->numfeats; i++)
    if (suff[i]>=0) feats->featfuncts[i]->suffmove=feats->featfuncts[suff[i]];
  for (i=0; i < feats->numschema; i++) {
    j=rdsw_int(fp,swp);
    if (i != j) failproc("file format incorrect: featschema");
    feats->schema[i]=malloc(sizeof(struct TWSchema));
    feats->schema[i]->numfuncts=rdsw_int(fp,swp);
    feats->schema[i]->functs=malloc(feats->schema[i]->numfuncts*sizeof(TWPtr));
    for (j=0; j < feats->schema[i]->numfuncts; j++) {
      ind=rdsw_int(fp,swp);
      feats->schema[i]->functs[j]=feats->featfuncts[ind];
    }
  }
  free(suff);
}

void read_bfmain(FILE *fp, FtrPtr feats, int swp)
{
  int i, j, ind, st, lb, lv, sc, ds, next, hash;
  FHashPtr fhash;
  FeatPtr feat;
  next=rdsw_int(fp,swp);
  while (next) {
    st=rdsw_int(fp,swp);
    lb=rdsw_int(fp,swp);
    lv=rdsw_int(fp,swp);
    sc=rdsw_int(fp,swp);
    ds=rdsw_int(fp,swp);
    hash=get_fhash(st,lb,lv,sc,feats);
    ind=rdsw_int(fp,swp);
    fhash=insert_feat(st,lb,lv,sc,feats,hash,ind);
    feat=fhash->feat;
    feat->dstate=ds;
    if (ind) {
      feat->data->score=rdsw_float(fp,swp);
      feat->data->bscore=rdsw_float(fp,swp);
      feat->data->numc=rdsw_int(fp,swp);
      feat->data->cfeat=malloc(feat->data->numc*sizeof(CFDataPtr));
      for (j=0; j < feat->data->numc; j++) {
	feat->data->cfeat[j]=malloc(sizeof(struct CFeatData));
	feat->data->cfeat[j]->label=rdsw_int(fp,swp);
	feat->data->cfeat[j]->idx=rdsw_int(fp,swp);
      }
    }
    next=rdsw_int(fp,swp);
  }
}

void read_bfback(FILE *fp, FtrPtr feats, int swp)
{
  int i, j, ind, st, lb, lv, sc, ds, hash, pcount=rdsw_int(fp,swp);
  FHashPtr fhash, bhash;
  FeatPtr feat;
  for (i=0; i < pcount; i++) {
    st=rdsw_int(fp,swp);
    lb=rdsw_int(fp,swp);
    lv=rdsw_int(fp,swp);
    sc=rdsw_int(fp,swp);
    ds=rdsw_int(fp,swp);
    ind=rdsw_int(fp,swp);
    fhash=find_feat(st,lb,lv,sc,feats,&hash);
    if (fhash == NULL) {
      fprintf(stderr,"%d %d %d %d %d %d\n",i,pcount,st,lb,lv,sc);
      failproc("hash not found");
    }
    bhash=find_feat(ds,ind,lv-1,sc,feats,&hash);
    if (bhash == NULL) failproc("bhash not found");
    fhash->feat->data->prefix=bhash->feat;
  }
  ind=rdsw_int(fp,swp);
  if (ind != pcount) failproc("indices don't match");
}

void read_bfscores(FILE *fp, FtrPtr feats, int swp)
{
  int i, j, p, ind;
  float *tmp;
  ScorePtr scr=feats->scores;
  scr->numind=scr->numalloc=rdsw_int(fp,swp);
  p=rdsw_int(fp,swp);
  scr->scores=malloc(scr->numind*sizeof(float));
  if (p) scr->pscores=malloc(scr->numind*sizeof(float));
  for (i=0; i < scr->numind; i++) {
    scr->scores[i]=rdsw_float(fp,swp);
    if (p) scr->pscores[i]=rdsw_float(fp,swp);
  }
  ind=rdsw_int(fp,swp);
  if (ind != scr->numind) failproc("indices don't match");
  scr->alcpsc=rdsw_int(fp,swp);
  scr->numpsc=rdsw_int(fp,swp);
  if (scr->alcpsc>0) scr->ppsc=malloc(scr->alcpsc*sizeof(float));
  for (i=0; i < scr->alcpsc; i++) scr->ppsc=0;
  for (i=0; i < scr->numpsc; i++) scr->ppsc[i]=rdsw_float(fp,swp);
}

void read_bfeatures(FILE *fp, ParPtr par, int swp)
{
  FtrPtr feats=par->features;
  LexPtr lex=par->lexicon;
  feats->numschema=rdsw_int(fp,swp);
  feats->numfeats=rdsw_int(fp,swp);
  feats->fstates=rdsw_int(fp,swp);
  feats->schema=malloc(feats->numschema*sizeof(TWSchPtr));
  feats->featfuncts=malloc(feats->numfeats*sizeof(TWPtr));
  read_bfsch(fp,feats,lex,swp);
  read_bfmain(fp,feats,swp);
  read_bfback(fp,feats,swp);
  read_bfscores(fp,feats,swp);
}

void read_blrn(FILE *fp,ScorePtr scr,int swp)
{
  int i, j, ind;
  scr->lrn=malloc(sizeof(struct LearnScore));
  scr->lrn->xvals=rdsw_int(fp,swp);
  scr->lrn->bsize=rdsw_int(fp,swp);
  scr->lrn->iter=rdsw_int(fp,swp);
  scr->lrn->ascores=scr->pscores;
  scr->lrn->scores=scr->scores;
  scr->lrn->apsc=scr->ppsc;
  scr->lrn->pscores=malloc(scr->numind*sizeof(float));
  scr->lrn->lastupd=malloc(scr->numind*sizeof(int));
  scr->lrn->ppsc=malloc(scr->alcpsc*sizeof(float));
  scr->lrn->lastpsc=malloc(scr->alcpsc*sizeof(int));
  scr->lrn->backs=malloc(scr->numind*sizeof(BORPtr));
  scr->lrn->fscores=malloc(scr->lrn->xvals*sizeof(FScorePtr));
  for (j=0; j < scr->lrn->xvals; j++) {
    scr->lrn->fscores[j]=malloc(sizeof(struct FltScore));
    scr->lrn->fscores[j]->rsc=malloc(scr->numind*sizeof(float));
  }
  for (i=0; i < scr->alcpsc; i++) {
    scr->lrn->ppsc[i]=0;
    scr->lrn->lastpsc[i]=0;
  }
  for (i=0; i < scr->numind; i++) {
    scr->lrn->backs[i]=NULL;
    scr->lrn->pscores[i]=rdsw_float(fp,swp);
    scr->lrn->lastupd[i]=rdsw_int(fp,swp);
    for (j=0; j < scr->lrn->xvals; j++)
      scr->lrn->fscores[j]->rsc[i]=rdsw_float(fp,swp);
  }
  ind=rdsw_int(fp,swp);
  if (ind != scr->numind) failproc("input mis-parse");
  for (i=0; i < scr->numpsc; i++) {
    scr->lrn->ppsc[i]=rdsw_float(fp,swp);
    scr->lrn->lastpsc[i]=rdsw_int(fp,swp);
  }
}

void write_parser(ParPtr par, char *ofile)
{
  int ind=1;                        /* for byte swp detection */
  FILE *fp;
  fp=gopenfile(ofile,"w");
  fwrite(&ind, sizeof(int), 1, fp);
  fwrite(&(par->lexicon->numWH), sizeof(int), 1, fp);
  fwrite(&(par->features->numFH), sizeof(int), 1, fp);
  fwrite(&(par->lexicon->numWDs), sizeof(int), 1, fp);
  write_blexicon(fp,par->lexicon);
  write_bfeatures(fp,par->features);
  if (par->features->scores->lrn != NULL)
    write_blrn(fp,par->features->scores);
  pclose(fp);
}

ParPtr read_parser(FILE *fp, int lrn)
{
  int ind=1, swp=0, numWH, numFH, numWD;
  ParPtr par;
  fread(&ind, sizeof(int), 1, fp);
  if (ind != 1) swp = 1;
  numWH = rdsw_int(fp,swp);
  numFH = rdsw_int(fp,swp);
  numWD = rdsw_int(fp,swp);
  par=ParserDefault(numWH,numFH,numWD,0);
  read_blexicon(fp,par->lexicon,swp);
  read_bfeatures(fp,par,swp);
  if (lrn) read_blrn(fp,par->features->scores,swp);
  return par;
}

WLPtr addcword(LexPtr lex, char *str, int *nhash)
{
  int j;
  WLPtr word=find_word(str,lex,nhash);
  if (word == NULL) {
    word=insert_word(str,lex,nhash[0],0,1);
    word->NTs = malloc(lex->numPTs*sizeof(float));
    for (j=0; j < lex->numPTs; j++) word->NTs[j]=-log(0.0);
  }
  return word;
}

void read_norms(char *readtree, TDTrPtr TC, LexPtr oldlex)
{
  int i=0, j, len = strlen(readtree), child=0, hash, nhash, firstw=1, oov;
  float tot;
  char readterm[MAXLABLEN], *str;
  WLPtr word=NULL, lword, nword, oword, uword;
  LexPtr lex = TC->parser->lexicon;
  while (i < len) {
    if (word != NULL && word->children==0 && !word->sym) firstw=0;
    lword = word;
    wordtoken(readtree,&readterm[0],&i,len,0);
    if (i >= len) continue;
    if (readterm[0] == '(') {
      word=find_word(&readterm[0],lex,&hash);
      word->tot++;
      continue;
    }
    oword=word=find_word(&readterm[0],oldlex,&hash);
    if (word == NULL) failproc("word not found in old lexicon!");
    if (lword == NULL) failproc("no POS tag for word!");
    nhash = hash;
    nword=find_nword(&readterm[0],oldlex,&nhash,word,word->label,lword,firstw,1,0);
    if (nword == NULL) str = &readterm[0];
    else str = nword->label;
    uword=word=addcword(lex,str,&nhash);
    if (nword == NULL) word->oov = 1;
    word->tot++;
    if (!word->oov) continue;
    strcpy(readterm,oword->label);
    word=find_zword(&readterm[0],oldlex,&nhash,oword,oword->label,lword);
    if (word == NULL) str=&readterm[0];
    else str = word->label;
    word=addcword(lex,str,&nhash);
    word->unkcat=uword;
  }
}

void read_terms(char *readtree, TDTrPtr TC)
{
  int i=0, j, len = strlen(readtree), child=1, hash, firstw=1, firstnt=1;
  char readterm[MAXLABLEN];
  WLPtr word=NULL, lword=NULL, llword;
  LexPtr lex = TC->parser->lexicon;
  while (i < len) {
    if (child==0) firstw=0;
    if ((TC->slcnts && llword != NULL && !llword->root && llword->children && 
	 !llword->terminal && lword != NULL && lword->children && !lword->terminal) || 
	(TC->slcprods && lword != NULL && llword == lword)) {
      sprintf(readterm,"%s%s",llword->label,lword->label);
      add_slcprods(TC,&readterm[0]);
    }
    llword = lword;
    lword = word;
    wordtoken(readtree,&readterm[0],&i,len,0);
    if (i >= len) continue;
    if (readterm[0] == '(') child=1;
    else child=0;
    word=find_word(&readterm[0],lex,&hash);
    if (word == NULL) {
      word=insert_word(&readterm[0],lex,hash,child,1-child);
      if (firstnt) word->root=1;
    }
    if (lword != NULL && !child) {
      if (readterm[0] > 64 && readterm[0] < 91 && !firstw) lword->caps++;
      if (noalpha(&readterm[0])) lword->sym++;
      lword->terminal=1;
    }
    word->tot++;
    if (firstnt) firstnt=0;
  }
}

int quotes(char *readterm)
{
  if (readterm[0]=='"' && readterm[1]==0) return 1;
  if (readterm[0]=='\'' && readterm[1]=='\'' && readterm[2]==0) return 1;
  if (readterm[0]=='`' && readterm[1]=='`' && readterm[2]==0) return 1;
  return 0;
}

SLPtr get_vstring(char *readstr, TDParPtr PC, int pos)
{
  int i=0, j, len = strlen(readstr), child=0, hash, nhash, firstw=1, oov;
  float tot;
  char readterm[MAXLABLEN], *str, sym[MAXLABLEN];
  WLPtr word=NULL, nword;
  LexPtr lex = PC->parser->lexicon;
  SLPtr string=NULL, fstring=NULL;
  sym[0]=0;
  while (i < len) {
    if (word != NULL && !word->sym) firstw=0;
    j=0;
    while (i < len && ws(readstr[i])) i++;
    if (i >= len) continue;
    while (i < len && !ws(readstr[i])) readterm[j++]=readstr[i++];
    readterm[j]=0;
    if (j==0) continue;
    if (quotes(&readterm[0])) {
      if (sym[0]!=0 && string !=NULL) {
	string->rightsym=malloc(5*sizeof(char));
	strcpy(string->rightsym,sym);
	sym[0]=0;
      }
      strcpy(sym,readterm);
      continue;
    }
    string=default_string(&readterm[0],string);
    string->word=find_word(&readterm[0],lex,&hash);
    if (string->word == NULL) 
      string->word=find_nsword(&readterm[0],lex,&hash,string->outlabel,firstw);
    if (string->word == NULL) {
      fprintf(stderr,"failed to find word: %s",string->outlabel);
      exit(1);
    }
    if (fstring == NULL) fstring=string;
    word=string->word;
    if (sym[0]!=0) {
      string->leftsym=malloc(5*sizeof(char));
      strcpy(string->leftsym,sym);
      sym[0]=0;
    }
    if (pos) {
      j=0;
      while (i < len && ws(readstr[i])) i++;
      if (i >= len) continue;
      while (i < len && !ws(readstr[i])) i++;
    }
  }
  if (sym[0]!=0) {
    string->rightsym=malloc(5*sizeof(char));
    strcpy(string->rightsym,sym);
    sym[0]=0;
  }
  return fstring;
}

SLPtr get_string(char *readstr, TDParPtr PC) {
  return get_vstring(readstr,PC,0);
}

SLPtr get_postring(char *readstr, TDParPtr PC) {
  return get_vstring(readstr,PC,1);
}
