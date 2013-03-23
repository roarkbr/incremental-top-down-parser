// parser.c
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

int yes_over_thresh(int num, int den, float thresh)
{
  float score=num;
  score/=den;
  if (score > thresh) return 1;
  return 0;
}

void nonterm_lex(LexPtr lex, LexPtr oldlex, int child, int term, int arr)
{
  int i, j, hash;
  float sym, caps, thresh=0.9;
  WLPtr pword, lword, word, oword;
  char *label;
  for (i = 1; i < oldlex->numWDs; i++) {
    oword=oldlex->WDs[i];
    if (oword->children!=child || oword->terminal!=term) continue;
    label=oword->label;
    word=find_word(label,lex,&hash);
    word=insert_word(label,lex,hash,child,term);
    if (arr > 0) {
      word->NTs = malloc(arr*sizeof(float));
      for (j=0; j < arr; j++) word->NTs[j]=0;
    }
    if (oword->slc != NULL && !oword->slc->orig) {
      pword=find_word(oword->slc->NewCat->label,lex,&hash);
      lword=find_word(oword->slc->LChild->label,lex,&hash);
      if (pword == NULL || lword == NULL) 
	failproc("lchild and parent not added in nonterm_lex");
      make_slcentry(pword,lword,word);
    }
    word->root = oword->root;
    oword->unkcat=word;
    if (term) {
      word->sym=yes_over_thresh(oword->sym,oword->tot,thresh);
      word->caps=yes_over_thresh(oword->caps,oword->tot,thresh);
    }
  }
}

int compress_lexicon(TDTrPtr TC,char *readtree)
{
  int i, j, k=0, ret=0, PT;
  FILE *fp;
  FtrPtr feats=TC->parser->features;
  TWPtr tfunct;
  LexPtr olex = TC->parser->lexicon, lex = LexDefault(olex->numWH, INITWD, 1);
  if (TC->verbose) fprintf(stderr,"Lexicon: comp. ");
  TC->parser->lexicon = lex;
  nonterm_lex(lex,olex,1,1,0);
  lex->numPTs=lex->numWDs;
  nonterm_lex(lex,olex,1,0,lex->numWDs);
  lex->numNTs=lex->numWDs;
  ret=lex->numWDs;
  for (j=0; j < feats->numfeats; j++) {
    tfunct=feats->featfuncts[j];
    if (tfunct->numqueries <= 0) continue;
    for (i=0; i < tfunct->numqueries; i++) {
      if (tfunct->WQueries[i]==NULL) continue;
      tfunct->WQueries[i]=tfunct->WQueries[i]->unkcat;
    }
  }
  fp = fopen(TC->tfile,"r");
  while (fgets(readtree,MAXLABLEN,fp) != NULL) {
    read_norms(readtree,TC,olex);
    verbose_count(TC->verbose,++k,10000);
  }
  fclose(fp);
  if (TC->verbose) fprintf(stderr," %d\n",k);
  return ret;
}

int acquire_lexicon(TDTrPtr TC,char *readtree)
{
  int k=0;
  FILE *fp = fopen(TC->tfile,"r");
  if (TC->verbose) fprintf(stderr,"Lexicon: build ");
  while (fgets(readtree,MAXLABLEN,fp) != NULL) {
    read_terms(readtree,TC);
    verbose_count(TC->verbose,++k,10000);
  }
  fclose(fp);
  if (TC->verbose) fprintf(stderr," %d\n",k);
  return k;
}

int divvy_xvals(int ind, TDTrPtr TC)
{
  int sents=ind;
  sents/=TC->xvals;
  if (sents > MAXVALS) {
    sents=MAXVALS-1;
    TC->xvals++;
    TC->xlast=0;
  }
  else if (sents < MNXVALS) {
    sents=MNXVALS;
    TC->xvals=ind;
    TC->xvals/=sents;
  }
  if (TC->xvals < 1) {
    sents=2;
    TC->xvals=ind;
    TC->xvals/=sents;    
  }
  return sents;
}

int get_lexicon(TDTrPtr TC,char *readtree)
{
  int sents=acquire_lexicon(TC,readtree);
  compress_lexicon(TC,readtree);
  return divvy_xvals(sents,TC);
}

void add_lookarr(TreePtr tree, int *NTs, TDTrPtr TC)
{
  int i;
  TreePtr ltree=tree;
  while (ltree != NULL && ltree->LSib == NULL) {
    ltree = ltree->Par;
    if (ltree != NULL && (ltree->NT->slc == NULL || ltree->NT->slc->orig)) {
      if (ltree->NT->windex >= TC->parser->lexicon->numNTs ||
	  ltree->NT->windex < TC->parser->lexicon->numPTs) {
	fprintf(stderr,"label out of range %d\n",ltree->NT->windex);
	exit(1);
      }
      NTs[ltree->NT->windex]=1;
    }
  }
  for (i=TC->parser->lexicon->numPTs; i < TC->parser->lexicon->numNTs; i++) {
    if (!NTs[i]) continue;
    TC->parser->lexicon->WDs[i]->NTs[tree->NT->windex]++;
    NTs[i]=0;
  }
}

void update_schema(FtrPtr feats, int xvals, int batch, int sch, int clbl, TreePtr tree, int unk)
{
  int j, k, st, lb, hash;
  TWSchPtr fschema=feats->schema[sch];
  FeatPtr currfeat=NULL, prefix;
  ScorePtr scr=feats->scores;
  FHashPtr currhash, holdhash;
  CFDataPtr cfeat;
  for (j=1; j < fschema->numfuncts; j++) {
    prefix=currfeat;
    k=fschema->functs[j]->findex;
    if (prefix==NULL) st=0;
    else st=prefix->dstate;
    lb=tree->featarr[k];
    holdhash=currhash=find_feat(st,lb,j,sch,feats,&hash);
    if (currhash == NULL) holdhash=currhash=insert_feat(st,lb,j,sch,feats,hash,1);
    currfeat=currhash->feat;
    currfeat->data->score++;
    currfeat->data->prefix=prefix;
    currhash=find_feat(currfeat->dstate,clbl,j+1,sch+1,feats,&hash);
    if (currhash == NULL) {
      currhash=insert_feat(currfeat->dstate,clbl,j+1,sch+1,feats,hash,0);
      currhash->feat->dstate=holdhash->feat->data->numc;
      insert_cfeat(feats,holdhash->feat->data,clbl,1);
    }
    cfeat=holdhash->feat->data->cfeat[currhash->feat->dstate];
    if (scr->scores[cfeat->idx]>0 && scr->clams[cfeat->idx]->pscore==NULL && 
	scr->clams[cfeat->idx]->sbatch!=batch) {
      scr->clams[cfeat->idx]->pscore=malloc(xvals*sizeof(float));
      for (k=0; k < xvals; k++) scr->clams[cfeat->idx]->pscore[k]=0;
      scr->clams[cfeat->idx]->pscore[scr->clams[cfeat->idx]->sbatch]=
	scr->scores[cfeat->idx];
    }
    if (scr->clams[cfeat->idx]->pscore == NULL) 
      scr->clams[cfeat->idx]->sbatch=batch;
    else scr->clams[cfeat->idx]->pscore[batch]++;
    scr->scores[cfeat->idx]++;
    if (unk) break;
  }
}

void update_feats(TreePtr tree, TDTrPtr TC, int batch)
{
  int i, clbl;
  FtrPtr feats = TC->parser->features;
  LexPtr lex=TC->parser->lexicon;
  TWSchPtr fschema;
  for (i=0; i < feats->numschema; i++) {
    fschema = feats->schema[i];
    clbl=tree->featarr[fschema->functs[0]->findex];
    if (clbl < 0) continue;
    update_schema(feats, TC->xvals, batch, i, clbl, tree, 0);
    if (lex->WDs[clbl]->unkcat != NULL) 
      update_schema(feats,TC->xvals,batch,i,lex->WDs[clbl]->unkcat->windex ,tree,1);
  }
}
 
void count_feats(TreePtr tree, int *NTs, TDTrPtr TC, int batch)
{
  while (tree != NULL) {
    if (tree->NT->terminal && tree->LSib == NULL) add_lookarr(tree,NTs,TC);
    update_feats(tree, TC, batch);
    tree = tree->prev;
  }
}

void collect_counts(TDTrPtr TC, char *readtree, int sents)
{
  int i, j, k=0, *NTs, batch=0;
  FILE *fp = fopen(TC->tfile,"r");
  TreePtr tree;
  NTs = malloc(TC->parser->lexicon->numNTs*sizeof(int));
  for (i=0; i < TC->parser->lexicon->numNTs; i++) NTs[i]=0;
  if (TC->verbose) fprintf(stderr,"Grammar: count ");
  while (fgets(readtree,MAXLABLEN,fp) != NULL) {
    tree=build_tree(readtree,TC->parser,0);
    count_feats(tree,NTs,TC,batch);
    free_trees(tree,NULL,1);
    verbose_count(TC->verbose,++k,10000);
    if (k%sents == 0 && batch < TC->xvals-1) batch++;
  }
  fclose(fp);
  free(NTs);
  if (TC->verbose) fprintf(stderr," %d\n",k);
}

float getschprob(ParPtr parser, int sch, TreePtr tree, int bottom, float *p)
{
  int i, j, lb, hash, ds=0, idx;
  float z=0, lz=-log(z);
  FtrPtr feats=parser->features;
  ScorePtr scr=feats->scores;
  FHashPtr fhash;
  FeatPtr feat=NULL, bfeat;
  TWSchPtr fsch=feats->schema[sch];
  for (j=1; j < fsch->numfuncts; j++) {
    lb=tree->featarr[fsch->functs[j]->findex];
    fhash=find_feat(ds,lb,j,sch,feats,&hash);
    if (fhash == NULL) break;
    feat=fhash->feat;
    ds=feat->dstate;
  }
  if (feat == NULL) return lz;
  bfeat=feat->data->prefix;
  lb=tree->featarr[fsch->functs[0]->findex];
  fhash=find_feat(feat->dstate,lb,feat->level+1,feat->schema+1,feats,&hash);
  if (fhash == NULL) return getboff(feats,bfeat,lb,feat->data->bscore,bottom,p);
  idx=feat->data->cfeat[fhash->feat->dstate]->idx;
  if (p!=NULL && scr->pscores!=NULL) p[0]=scr->pscores[idx];
  return scr->scores[idx];
}

void fill_thisarray(TreePtr tree, FtrPtr feats)
{
  tree->HChild=get_hnode(tree);
  fill_ccnodes(tree,0);
  fill_farrays(tree,tree,feats,0);     
}

void eval_gold(TreePtr tree)
{
  if (tree==NULL || tree->prev==NULL || tree->prev->gold==NULL) return;
  if (tree->prev->gold->move==tree->NT) tree->gold=tree->prev->gold->next;
  else if (tree->prev->gold->move->terminal && tree->NT->terminal) {
    if (tree->prev->gold->gtree->LookWord != tree->LookWord) {
      fprintf(stderr,"misalignment in gold: %s %s\n",tree->prev->gold->gtree->LookWord->word->label,tree->LookWord->word->label);
      exit(1);
    }
    tree->gold=tree->prev->gold->next;
  }
}

TreeHPtr stop_trees(TreeHPtr ctrees, TDParPtr PC)
{
  TreePtr tree, lsib, ntree;
  float p=0;
  if (ctrees==NULL) return NULL;
  lsib=tree=ctrees->tree;
  if (!tree->NT->children) lsib=tree->Par;
  if (lsib != NULL && lsib->Par != NULL) {
    tree->points--;
    tree=ctrees->tree=default_tree(PC->parser->lexicon->WDs[0],tree);
    eval_gold(tree);
    tree->points++;
    fill_thisarray(tree,PC->parser->features);
    tree->Score=tree->prev->Score+
      getschprob(PC->parser,PC->parser->sch,tree,1,&p);
    /*tree->Score+=p;*/            /* add pscore every move */
    tree->PScore=tree->prev->PScore+p;
    return stop_trees(ctrees,PC);
  }
  ctrees->next=stop_trees(ctrees->next,PC);
  ctrees->right=stop_trees(ctrees->right,PC);
  ctrees->left=stop_trees(ctrees->left,PC);
  return ctrees;
}

TreePtr finishtree(TreePtr tree, LexPtr lex)
{
  int i, bpos, cnt=0, addstop=0;
  float best;
  TreePtr lastree=tree;
  SLPtr string=tree->LookUpd;
  if (tree->LookWord == NULL) addstop=1;
  else string=tree->LookWord->next;
  while (string != NULL) {
    if (cnt++ >= 1000) failproc("oops, 1000 0 adds");
    bpos=-1;
    for (i=1; i < lex->numPTs; i++) {
      if (bpos>=0 && string->word->NTs[i] >= best) continue;
      bpos=i;
      best=string->word->NTs[i];
    }
    tree=default_tree(lex->WDs[bpos],tree);
    tree->LookWord=string;
    string=string->next;
  }
  while (addstop || (cnt < 1000 && tree->Par->NT->windex != lex->root)) {
    tree=default_tree(lex->WDs[0],tree);
    if (addstop) addstop=0;
    if (cnt++ >= 1000) failproc("oops, 1000 1 adds");
  }
  return tree;
}

TreeHPtr failheap(TreePtr tree, LexPtr lex)
{
  float score=tree->Score;
  tree=finishtree(tree,lex);
  tree->Score=score;
  return default_theap(tree,score,-1,0);
}

double renorm(TreeHPtr ctrees,double norm)
{
  double amb=0;
  if (ctrees==NULL) return amb;
  ctrees->tree->Score-=norm;
  ctrees->score=ctrees->tree->Score;
  amb=exp(-ctrees->score);
  amb*=ctrees->score;
  amb+=renorm(ctrees->next,norm);
  amb+=renorm(ctrees->right,norm);
  amb+=renorm(ctrees->left,norm);
  return amb;
}

TreeHPtr stop_heap(TreeHPtr ctrees, TDParPtr PC, TreePtr holdtree, int *fails, double *plp, double *normp, double *slp, double *snormp)
{
  TreeHPtr ntrees=NULL, thish;
  double lp, pr, cpr=0;
  ctrees=stop_trees(ctrees,PC);
  while (ctrees != NULL) {
    thish=ctrees;
    ctrees=perc_heap(ctrees);
    lp=thish->tree->Score;
    pr=exp(-lp);
    plp[0]+=pr*lp;
    normp[0]+=pr;
    cpr+=pr;
    slp[0]+=pr*lp;
    snormp[0]+=pr;
    thish->tree->Score+=thish->tree->PScore;
    thish->tree->PScore=0;
    thish->score=thish->tree->Score;
    ntrees=insert_theap(thish,ntrees);
  }
  if (PC->allpossible) PC->wscores[0]=PC->sscores[0]=cpr;
  if (holdtree!=NULL) holdtree->points--;
  if (ntrees != NULL && holdtree!=NULL) free_trees(holdtree,NULL,0);
  if (ntrees != NULL) return ntrees;
  fails[0]++;
  if (holdtree==NULL) return ntrees;
  return failheap(holdtree,PC->parser->lexicon);
}

void calc_stop_heap(TreeHPtr ctrees, TDParPtr PC, double *plp, double *normp, double *slp, double *snormp)
{
  TreeHPtr thish;
  double lp, pr, cpr=0;
  ctrees=stop_trees(ctrees,PC);
  while (ctrees != NULL) {
    thish=ctrees;
    ctrees=perc_heap(ctrees);
    lp=thish->tree->Score;
    pr=exp(-lp);
    cpr+=pr;
    plp[0]+=pr*lp;
    normp[0]+=pr;
    slp[0]+=pr*lp;
    snormp[0]+=pr;
    free_theap(thish);
  }
  if (PC->allpossible) PC->wscores[0]=PC->sscores[0]=cpr;
}

float *init_lscore(TreePtr tree, int numPT)
{
  int i;
  float *scores=tree->LookScore, def=0;
  if (scores == NULL) scores=malloc(numPT*sizeof(float));
  for (i=0; i < numPT; i++) scores[i]=-log(def);
  return scores;
}

void upnew_looks(TreePtr tree, TreePtr ctree, SLPtr string, TDParPtr PC)
{
  LexPtr lex=PC->parser->lexicon;
  int i,j, numPT=lex->numPTs, numNT=lex->numNTs, lb, clbl=-1;
  float stop=0, lz=-log(stop);
  double score;
  if (tree == NULL || tree->LookUpd == string) return;
  lb=tree->NT->windex;
  upnew_looks(tree->Par, tree, string, PC);
  tree->LookScore = init_lscore(tree,numPT);
  tree->LookUpd=string;
  if (ctree==NULL && lb != lex->root) {
    for (i=1; i < numPT; i++) tree->LookScore[i]=tree->Par->LookScore[i];
    return;
  }
  if (ctree != NULL) clbl=ctree->NT->windex;
  get_ntprobs(lb,clbl,PC);
  if (lex->nts[0] < lz) stop=lex->nts[0];
  for (i=1; i < numPT; i++) { 
    tree->LookScore[i]=lex->WDs[string->word->windex]->NTs[i];
    if (tree->LookScore[i] >= lz) continue;
    score=stop;
    if (tree->Par != NULL && tree->Par->LookScore[i]<lz) 
      score *= exp(-tree->Par->LookScore[i]);
    else score=0.0;
    if (lex->nts[i] > 0.0) score+=lex->nts[i];
    for (j=numPT; j < numNT; j++) {
      if (lex->nts[j] <= 0.0 || lex->WDs[j]->NTs[i] == lz) continue;
      score+=lex->nts[j]*exp(-lex->WDs[j]->NTs[i]);
    }
    tree->LookScore[i]-=log(score);
  }
}

double logsum(double A, double B)
{
  double ret=A, comb=1;
  if (A > B) return logsum(B,A);
  comb += exp(A-B);
  ret -= log(comb);
  return ret;
}

float calc_wscore(TreePtr tree)
{
  float wscore=tree->WScore;
  if (tree->prev != NULL) wscore+=calc_wscore(tree->prev);
  return wscore;
}

double calc_norm(TreeHPtr ctrees, TDParPtr PC, int wn)
{
  double tscore=0, score=-log(tscore);
  if (ctrees == NULL) return score;
  if (wn) {
    score=ctrees->score+PC->parser->norm;
    score-=ctrees->tree->WScore;
  }
  else {
    ctrees->tree->Score+=ctrees->tree->PScore;
    ctrees->tree->PScore=0;
    score=ctrees->tree->Score;
  }
  tscore=calc_norm(ctrees->next,PC,wn);
  score=logsum(score,tscore);
  tscore=calc_norm(ctrees->right,PC,wn);
  score=logsum(score,tscore);
  tscore=calc_norm(ctrees->left,PC,wn);
  score=logsum(score,tscore);
  return score;
}

double calc_onenorm(TreeHPtr ctrees, TDParPtr PC)
{
  double tscore=0, zs=-log(tscore), score=zs;
  if (ctrees == NULL) return score;
  if (ctrees->rank==1) score=ctrees->score;
  tscore=calc_onenorm(ctrees->next,PC);
  if (tscore != zs) {
    if (score==zs) score=tscore;
    else score=logsum(score,tscore);
  }
  tscore=calc_onenorm(ctrees->right,PC);
  if (tscore != zs) {
    if (score==zs) score=tscore;
    else score=logsum(score,tscore);
  }
  tscore=calc_onenorm(ctrees->left,PC);
  if (tscore != zs) {
    if (score==zs) score=tscore;
    else score=logsum(score,tscore);
  }
  return score;
}

int calc_openbr(TreePtr tree)
{
  TreePtr par=tree->Par;
  int openbr=-1;
  while (par != NULL) {
    par=par->Par;
    openbr++;
  }
  return openbr;
}

double weighted_openb(TreeHPtr ctrees)
{
  double amb=0;
  if (ctrees==NULL) return amb;
  amb=calc_openbr(ctrees->tree);
  amb*=exp(-ctrees->score);
  amb+=weighted_openb(ctrees->next);
  amb+=weighted_openb(ctrees->right);
  amb+=weighted_openb(ctrees->left);
  return amb;
}

double upnorm(TreeHPtr ctrees,TDParPtr PC, int wn)
{
  double norm=calc_norm(ctrees,PC,wn), amb=0;
  if (wn) {
    PC->parser->snorm=norm;
    amb=weighted_openb(ctrees);
  }
  else {
    amb=renorm(ctrees,norm);
    PC->parser->norm+=norm;
    PC->parser->cond=calc_onenorm(ctrees,PC);
  }
  return amb;
}

double upnorms(TreeHPtr ctrees, TDParPtr PC, double *openb)
{
  if (ctrees==NULL) return 0;
  double amb=upnorm(ctrees,PC,0);
  openb[0]=upnorm(ctrees,PC,1);
  return amb;
}

double new_looks(TreeHPtr ctrees, SLPtr string, TDParPtr PC)
{
  double tscore=0, score=-log(tscore);
  if (ctrees == NULL) return score;
  ctrees->tree->Score+=ctrees->tree->PScore;
  ctrees->tree->PScore=0;
  score=ctrees->tree->Score;
  upnew_looks(ctrees->tree, NULL, string, PC);
  tscore=new_looks(ctrees->next, string, PC);
  score=logsum(score,tscore);
  tscore=new_looks(ctrees->right, string, PC);
  score=logsum(score,tscore);
  tscore=new_looks(ctrees->left, string, PC);
  score=logsum(score,tscore);
  return score;
}

float calc_derv(TreeHPtr ctrees, int stopterm)
{
  float score=0, cnt=0;
  TreePtr tree;
  if (ctrees == NULL) return score;
  tree=ctrees->tree;
  score=exp(-tree->Score);
  tree=tree->prev;
  cnt++;
  while (tree != NULL && (!stopterm || !tree->NT->terminal)) {
    tree=tree->prev;
    cnt++;
  }
  score*=cnt;
  score+=calc_derv(ctrees->next,stopterm);
  score+=calc_derv(ctrees->right,stopterm);
  score+=calc_derv(ctrees->left,stopterm);
  return score;
}

float calc_thresh(float best, int ins, float bthresh)
{
  float thresh=best+bthresh, ifact=log10(ins), lc=3*ifact;
  thresh -= lc;
  return thresh;
}

void pos_expt(TreePtr tree, float *nts, int pts, float norm)
{
  int i;
  float score=0, lz=-log(score);
  if (tree==NULL) return;
  tree->Score-=norm;
  if (tree->LookScore == NULL) return;
  for (i=1; i < pts; i++) {
    if (tree->LookScore[i]>=lz) continue;
    if (nts[i]==lz) score=0;
    else score=exp(-nts[i]);
    score+=exp(-tree->Score-tree->LookScore[i]);
    nts[i]=-log(score);
  }
}

void pos_expts(TreeHPtr ctrees, float *nts, int pts, float norm)
{
  if (ctrees == NULL) return;
  pos_expt(ctrees->tree, nts, pts, norm);
  pos_expts(ctrees->next, nts, pts, norm);
  pos_expts(ctrees->right, nts, pts, norm);
  pos_expts(ctrees->left, nts, pts, norm);
}

int best_pos(TreeHPtr ctrees, float *nts, int pts, float *bs, int *best, int cands, float norm)
{
  float lz=-log(0);
  int i, j, k, ret=0;
  pos_expts(ctrees,nts,pts,norm);
  for (i=1; i < pts; i++) {
    if (nts[i]>=lz) continue;
    for (j=0; j < cands; j++) {
      if (best[j] != 0 && bs[j] <= nts[i]) continue;
      for (k=cands-1; k > j; k--) {
	bs[k]=bs[k-1];
	best[k]=best[k-1];
      }
      bs[j]=nts[i];
      best[j]=i;
      if (ret < cands) ret++;
      break;
    }
  }
  return ret;
}

TreeHPtr update_look(TreeHPtr ctrees, SLPtr string, TDParPtr PC)
{
  TreePtr tree;
  TreeHPtr ntrees=NULL, thish, thath;
  int i, pos, cands=PC->maxpos, pts=PC->parser->lexicon->numPTs, 
    bestpos[PC->maxpos];
  float minsc=0, thresh=-log(minsc), *nts=PC->parser->lexicon->nts, 
    bs[PC->maxpos], norm, best=thresh;
  norm=new_looks(ctrees,string,PC);
  for (i=0; i < pts; i++) nts[i]=thresh;
  for (i=0; i < cands; i++) bestpos[i]=0;
  cands=best_pos(ctrees,nts,pts,&bs[0],&bestpos[0],cands,norm);
  while (ctrees != NULL) {
    thish=ctrees;
    ctrees=perc_heap(ctrees);
    tree=thish->tree;
    for (i=0; i < cands; i++) {
      pos=bestpos[i];
      if (tree->LookScore[pos] >= thresh) continue;
      thath=default_theap(tree,tree->Score+tree->LookScore[pos],pos,thish->rank);
      ntrees=insert_theap(thath,ntrees);
      if (thath->score < best) {
	best=thath->score;
	thresh=calc_thresh(best,1,PC->thresh);
      }
    }
    free_theap(thish);
  }
  return ntrees;
}

FeatPtr getmaxc(ParPtr parser, FeatPtr firfeat, TreePtr tree, int min)
{
  int j, lb, hash;
  FtrPtr feats=parser->features;
  FeatPtr feat=firfeat;
  FHashPtr fhash;
  TWSchPtr fsch=feats->schema[parser->sch];
  for (j=min; j < fsch->numfuncts; j++) {
    lb=tree->featarr[fsch->functs[j]->findex];
    fhash=find_feat(feat->dstate,lb,j,parser->sch,feats,&hash);
    if (fhash == NULL) break;
    feat=fhash->feat;
  }
  return feat;
}

FeatPtr getminc(ParPtr parser, TreePtr tree, int preterm, int min)
{
  FtrPtr feats=parser->features;
  LexPtr lex=parser->lexicon;
  int i, st=0, lb, hash;
  float z=0, lz=-log(z);
  FDataPtr fdata;
  FHashPtr fhash;
  TWSchPtr fsch=feats->schema[parser->sch];
  for (i=1; i < min; i++) {
    lb=tree->featarr[fsch->functs[i]->findex];
    fhash=find_feat(st,lb,i,parser->sch,feats,&hash);
    if (fhash == NULL) return NULL;
    st=fhash->feat->dstate;
  }
  for (i=0; i < lex->numNTs; i++) {
    parser->exps[i]=lz;
    parser->pexps[i]=0.0;
    parser->fd[i]=1;
  }
  fdata=fhash->feat->data;
  for (i=0; i < fdata->numc; i++) parser->fd[fdata->cfeat[i]->label]=0;
  if (tree->LSib == NULL) parser->fd[0]=1;
  for (i=1; i < lex->numNTs; i++) {
    if (parser->fd[i]) continue;
    if (i >= lex->numPTs) {
      if (lex->WDs[i]->NTs[preterm] == lz) parser->fd[i]=1;
    }
    else if (i != preterm || tree->Par->LookScore[i] == lz) parser->fd[i]=1;
  }
  return fhash->feat;
}

void getexps(FtrPtr feats, LexPtr lex, FeatPtr feat, FeatPtr ifeat, float *exps, float *pexps, int *fd, int numNT)
{
  int i,numc=feat->data->numc,inumc=ifeat->data->numc,same=0;
  float bcost=feat->data->bscore, ibcost=ifeat->data->bscore;
  ScorePtr scr=feats->scores;
  FeatPtr bfeat=feat->data->prefix, ibfeat=ifeat->data->prefix;
  CFDataPtr *cfeat=feat->data->cfeat, *icfeat=ifeat->data->cfeat;
  if (ifeat==feat) same=1;
  for (i=0; i < numc; i++) {
    if (fd[cfeat[i]->label] || (!same && cfeat[i]->label == 0)) continue;
    exps[cfeat[i]->label]=scr->scores[cfeat[i]->idx];
    if (scr->pscores!=NULL) pexps[cfeat[i]->label]=scr->pscores[cfeat[i]->idx];
    fd[cfeat[i]->label]=1;
  }
  if (!same) {
    while (i < inumc && icfeat[i]->label != 0) i++;
    if (i < inumc) {
      exps[icfeat[i]->label]=scr->scores[icfeat[i]->idx];
      if (scr->pscores!=NULL) 
	pexps[icfeat[i]->label]=scr->pscores[icfeat[i]->idx];
    }
    else exps[0]=getboff(feats,ibfeat,0,ibcost,0,&(pexps[0]));
    fd[0];
  }
  for (i=0; i < numNT; i++) {
    if (fd[i]) continue;
    exps[i]=getboff(feats,bfeat,i,bcost,0,&(pexps[i]));
  }
}

void fill_looksc(TreePtr tree, LexPtr lex, int preterm, int wd, float loox)
{
  int i;
  if (tree->LookScore==NULL) tree->LookScore=malloc(lex->numNTs*sizeof(float));
  for (i=0; i < lex->numNTs; i++) tree->LookScore[i]=-log(0);
  tree->LookScore[preterm]=loox;
}

float new_loox(WLPtr NT, TreePtr ptree, LexPtr lex, int preterm, int wd)
{
  if (NT->terminal) return lex->WDs[wd]->NTs[preterm];
  if (NT->children) return NT->NTs[preterm]+lex->WDs[wd]->NTs[preterm];
  if (ptree->NT->windex != lex->root) return ptree->Par->LookScore[preterm];
  return -log(0);
}

TreePtr copy_ntree(TreePtr ntree, int pt, int nt, LexPtr lex, float *exps, float *pexps, float loox, FtrPtr feats)
{
  TreePtr xtree=ntree, prev=ntree->prev;
  xtree=copy_tree(lex->WDs[nt],ntree,feats->numfeats);
  xtree->Score=prev->Score+exps[nt];  /*+pexps[nt];  add pscore every move */
  xtree->PScore=prev->PScore+pexps[nt];
  fill_lclfarrays(xtree,xtree,feats);
  fill_looksc(xtree,lex,xtree->LookPT,xtree->LookWord->word->windex,loox);
  return xtree;
}

TreePtr newcand_tree(LexPtr lex, FtrPtr feats, TreePtr tree, SLPtr string, int ntcat, int ptcat, TreePtr ntree)
{
  default_link(lex->WDs[ntcat],tree,ntree);
  ntree->LookPT=ptcat;
  ntree->LookUpd=ntree->LookWord=string;
  fill_thisarray(ntree,feats);
  return ntree;
}

void fillexps(TDParPtr PC, FeatPtr ifeat, TreePtr ntree, TreePtr stree, int min)
{
  LexPtr lex=PC->parser->lexicon;
  FtrPtr feats=PC->parser->features;
  FeatPtr feat, sfeat;
  sfeat=feat=getmaxc(PC->parser,ifeat,ntree,min);
  if (stree != NULL) sfeat=getmaxc(PC->parser,ifeat,stree,min);
  getexps(feats,lex,feat,sfeat,PC->parser->exps,PC->parser->pexps,
	  PC->parser->fd,lex->numNTs);
}

TreeHPtr make_ctrees(TreeHPtr ctrees, TreeHPtr ptermh, TreeHPtr thish)
{
  free_theap(thish);
  if (ptermh == NULL) return ctrees;
  ptermh->right=ctrees;
  return ptermh;
}

TreeHPtr make_moves(TreeHPtr ctrees, TreeHPtr thish, TreePtr ntree, TreePtr stree, WLPtr sword, float thresh, TDParPtr PC)
{
  int i, pt=thish->preterm, *fd=PC->parser->fd;
  float z=0, lz=-log(z), loox, wds, *exps=PC->parser->exps, bascore, p,
    *pexps=PC->parser->pexps;
  LexPtr lex=PC->parser->lexicon;
  FtrPtr feats=PC->parser->features;
  TreeHPtr treeh, ptermh=NULL;
  TreePtr xtree, ptree=get_ptree(ntree);
  bascore=ntree->prev->Score;
  if (exps[pt] < lz) {
    p=0;
    if (PC->allpossible && sword->windex == PC->parser->lexicon->numWDs-1)
      wds=0;
    else wds=getschprob(PC->parser,PC->parser->wsch,ntree,1,&p);
    if (sword->oov) wds=new_loox(lex->WDs[pt],ptree,lex,pt,sword->windex);
    if (wds+exps[pt]+bascore < thresh) {
      ntree=copy_ntree(ntree,pt,pt,lex,exps,pexps,wds,feats);
      ntree->Score+=wds;    /*+p; add pscore every move */
      ntree->WScore+=wds;
      ntree->PScore+=p;
      eval_gold(ntree);
      ptermh=default_theap(ntree,ntree->Score,pt,thish->rank);
    }
  }
  for (i=lex->numNTs-1; i >= 0; i--) {
    if (i < lex->numPTs) i=0;
    if (exps[i] >= lz) continue;
    loox=new_loox(lex->WDs[i],ptree,lex,pt,sword->windex);
    if (loox+exps[i]+bascore >= thresh) continue;
    if (i==0 && stree!=NULL) 
      xtree=copy_ntree(stree,pt,i,lex,exps,pexps,loox,feats);
    else xtree=ntree=copy_ntree(ntree,pt,i,lex,exps,pexps,loox,feats);
    eval_gold(xtree);
    treeh=default_theap(xtree,xtree->Score,pt,thish->rank);
    treeh->score+=xtree->LookScore[pt];
    ctrees=insert_theap(treeh,ctrees);
  }
  return make_ctrees(ctrees,ptermh,thish);
}

TreeHPtr expush(TreeHPtr ctrees,SLPtr string, float thresh, TDParPtr PC)
{
  int min=3, pt=ctrees->preterm;
  float bascore=0;
  LexPtr lex=PC->parser->lexicon;
  FtrPtr feats=PC->parser->features;
  FeatPtr ifeat;
  TreeHPtr thish=ctrees;
  TreePtr tree=thish->tree, ntree, ptree, stree=NULL;
  ctrees=perc_heap(ctrees);
  ntree=newcand_tree(lex,feats,tree,string,pt,pt,PC->parser->ntree);
  ifeat=getminc(PC->parser,ntree,pt,min);
  if (ifeat == NULL) return make_ctrees(ctrees,NULL,thish);
  if (ntree->HChild == NULL && ntree->LSib != NULL)
    stree=newcand_tree(lex,feats,tree,string,0,pt,PC->parser->stree);
  fillexps(PC,ifeat,ntree,stree,min);
  return make_moves(ctrees,thish,ntree,stree,string->word,thresh,PC);
}

TreePtr get_holdtree(TreeHPtr ntrees, TreePtr holdtree)
{
  if (ntrees == NULL) return holdtree;
  if (holdtree!=NULL) {
    holdtree->points--;
    free_trees(holdtree,NULL,0);
  }
  holdtree=ntrees->tree;
  holdtree->points++;
  return holdtree;
}

TreeHPtr advance_parse(TreeHPtr ctrees, SLPtr string, TDParPtr PC)
{
  int ins=0;
  float best=-log(ins), thr=best;
  TreeHPtr ntrees=NULL, treeh, thish;
  while (ctrees!=NULL && (ctrees->tree->NT->terminal || ctrees->score<thr)) {
    if (!ctrees->tree->NT->terminal || ctrees->tree->LookWord != string) {
      ctrees=expush(ctrees,string,thr,PC);
      continue;
    }
    thish=ctrees;
    ctrees=perc_heap(ctrees);
    if (thish->score >= thr) {
      free_theap(thish);
      continue;
    }
    ntrees=insert_theapns(thish,ntrees);
    if (thish->score < best) best=thish->score;
    thr=calc_thresh(best,ins++,PC->thresh);
  }
  free_theaps(ctrees);
  return ntrees;
}

TreeHPtr calc_allparse(TreeHPtr ctrees, SLPtr string, TDParPtr PC, double *plp, double *normp, double *slp, double *snormp)
{
  int ins=0, min=3, i, pt;
  float z=0, lz=-log(z), p, wds;
  double pr, lp;
  WLPtr sword=string->word, tword;
  FeatPtr ifeat;
  TreeHPtr ntrees=NULL, treeh, thish;
  TreePtr ptree;
  LexPtr lex=PC->parser->lexicon;
  for (i=lex->numNTs; i < lex->numWDs-1; i++) 
    PC->wscores[i]=PC->sscores[i]=0;
  plp[0]=normp[0]=slp[0]=snormp[0]=0;
  while (ctrees!=NULL) {
    thish=ctrees;
    ctrees=perc_heap(ctrees);
    p=0;
    pt=thish->preterm;
    ptree=get_ptree(thish->tree);
    lp=thish->tree->Score;
    pr=exp(-lp);
    slp[0]+=pr*lp;
    snormp[0]+=pr;
    for (i=lex->numNTs; i < lex->numWDs-1; i++) {
      string->word=tword=lex->WDs[i];
      if (tword->NTs[pt]>=lz) continue;
      fill_thisarray(thish->tree,PC->parser->features);
      ifeat=getminc(PC->parser,thish->tree,thish->preterm,min);
      wds=getschprob(PC->parser,PC->parser->wsch,thish->tree,1,&p);
      if (tword->oov) wds=new_loox(lex->WDs[pt],ptree,lex,pt,tword->windex);
      lp=thish->tree->Score+wds;
      pr=exp(-lp);
      PC->wscores[i]+=pr;
      PC->sscores[i]+=exp(-thish->tree->Score);
      plp[0]+=pr*lp;
      normp[0]+=pr;
    }
    ntrees=insert_theapns(thish,ntrees);
  }
  string->word=sword;
  return ntrees;
}

TreeHPtr filter_allparse(TreeHPtr ctrees, SLPtr string, TDParPtr PC)
{
  int ins=0, min=3, pt;
  double lp;
  float z=0, lz=-log(z), p, wds;
  WLPtr sword=string->word;
  FeatPtr ifeat;
  TreePtr ptree;
  LexPtr lex=PC->parser->lexicon;
  TreeHPtr ntrees=NULL, treeh, thish;
  while (ctrees!=NULL) {
    thish=ctrees;
    ctrees=perc_heap(ctrees);
    pt=thish->preterm;
    ptree=get_ptree(thish->tree);
    p=0;
    fill_thisarray(thish->tree,PC->parser->features);
    ifeat=getminc(PC->parser,thish->tree,thish->preterm,min);
    wds=getschprob(PC->parser,PC->parser->wsch,thish->tree,1,&p);
    if (sword->oov) wds=new_loox(lex->WDs[pt],ptree,lex,pt,sword->windex);
    if (wds == lz) {
      free_theap(thish);
      continue;      
    }
    thish->tree->Score+=wds;
    thish->tree->WScore+=wds;
    ntrees=insert_theapns(thish,ntrees);
  }
  return ntrees;
}

TreeHPtr percsort(TreeHPtr ctrees)
{
  TreeHPtr ntrees=NULL, thish;
  while (ctrees != NULL) {
    thish=ctrees;
    ctrees=perc_heap(ctrees);
    thish->left=thish->right=thish->next=NULL;
    thish->under=1;
    ntrees=insert_theap(thish,ntrees);
  }
  return ntrees;
}

TreeHPtr show_partial(TreeHPtr xtrees, TDParPtr PC, FILE *fp, SLPtr string, int cnt)
{
  int c=0, openbr;
  float wscore;
  TreePtr tree;
  TreeHPtr ntrees=NULL, thish, ctrees=percsort(xtrees);
  fprintf(fp,"Partial parses for words 1-%d\t(%s)\n",cnt-1,string->outlabel);
  while (ctrees != NULL) {
    thish=ctrees;
    ctrees=perc_heap(ctrees);
    tree=thish->tree;
    thish->left=thish->right=thish->next=NULL;
    thish->under=1;
    ntrees=insert_theap(thish,ntrees);
    if (c < PC->nbest) {
      wscore=calc_wscore(tree);
      openbr=calc_openbr(tree);
      fprintf(fp,"\t%d %f %f %f %f %d  ",c+1,tree->Score,tree->Score+PC->parser->norm,wscore,tree->Score+PC->parser->norm-wscore,openbr);
      show_tree(fp,tree,1);
    }
    c++;
  }
  fflush(fp);
  return ntrees;
}

TreeHPtr calc_condrank(TreeHPtr xtrees, double *tcond)
{
  int rank=1;
  TreeHPtr ntrees=NULL, thish, ctrees=percsort(xtrees);
  if (ctrees == NULL) return ctrees;
  if (ctrees->rank==1) tcond[0]=ctrees->score;
  else tcond[0]=-log(0);
  while (ctrees != NULL) {
    thish=ctrees;
    ctrees=perc_heap(ctrees);
    thish->rank=rank++;
    thish->left=thish->right=thish->next=NULL;
    thish->under=1;
    ntrees=insert_theap(thish,ntrees);
  }
  return ntrees;
}

void show_prefix(FILE *fp,int cnt,char *lbl, ParPtr parser, double lnorm, double amb, double openb, double lcond, double tcond, double dsteps, int wtype)
{
  int clupd=0;
  char wc=':';
  if (wtype==1) wc='+';
  if (wtype==2) wc='-';
  if (cnt > 10) 
    fprintf(fp,"pfix%c%d %13s",wc,cnt-1,lbl);
  else fprintf(fp,"pfix%c%d %14s",wc,cnt-1,lbl);
  fprintf(fp," %7.3f %6.3f %5.3f %6.3f %6.3f %5.2f %5.2f %5.2f %4.1f\n",parser->norm,parser->norm-lnorm,parser->snorm-lnorm,parser->norm-parser->snorm,amb,openb,exp(lcond-parser->cond),exp(lcond-tcond),dsteps);
  parser->pstat->surprisal+=parser->norm-lnorm;
  parser->pstat->nonlex+=parser->snorm-lnorm;
  parser->pstat->lex+=parser->norm-parser->snorm;
  parser->pstat->ambig+=amb;
  parser->pstat->open+=openb;
  parser->pstat->rernk+=exp(lcond-parser->cond);
  parser->pstat->toprr+=exp(lcond-tcond);
  parser->pstat->dsteps+=dsteps;
  if (wtype==0) {
    parser->pstat->owds++;
    parser->pstat->lupd=0;
  }
  else {
    if (wtype==1) {
      parser->pstat->cwds++;
      parser->pstat->lupd=1;
    }
    if (parser->pstat->lupd) {
      parser->pstat->cl_surprisal+=parser->norm-lnorm;
      parser->pstat->cl_nonlex+=parser->snorm-lnorm;
      parser->pstat->cl_lex+=parser->norm-parser->snorm;
      parser->pstat->cl_ambig+=amb;
      parser->pstat->cl_open+=openb;
      parser->pstat->cl_rernk+=exp(lcond-parser->cond);
      parser->pstat->cl_toprr+=exp(lcond-tcond);    
      parser->pstat->cl_dsteps+=dsteps;
    }
  }
}

void show_aprefix(FILE *fp,int cnt,char *lbl, double a, double b, int rank, double lexent, double expsyns, double lograt, int wtype)
{
  int clupd=0;
  char wc=':';
  if (wtype==1) wc='+';
  if (wtype==2) wc='-';
  if (cnt > 10) 
    fprintf(fp,"afix%c%d %13s",wc,cnt-1,lbl);
  else fprintf(fp,"afix%c%d %14s",wc,cnt-1,lbl);
  fprintf(fp," %7.3f %6.3f %6.3f %6.3f %6.3f %6d %7.3f\n",a,b,lexent,expsyns,lexent-expsyns,rank,lograt);
  fflush(fp);
}

void show_norms(FILE *fp,float c,ParPtr parser, int type)
{
  if (type==0)
    fprintf(fp,"%6.3f %5.3f %6.3f %6.3f %5.2f %5.2f %5.2f %4.1f\n",
	    parser->pstat->surprisal/c,
	    parser->pstat->nonlex/c,
	    parser->pstat->lex/c,
	    parser->pstat->ambig/c,
	    parser->pstat->open/c,
	    parser->pstat->rernk/c,
	    parser->pstat->toprr/c,
	    parser->pstat->dsteps/c);
  if (type==1)
    fprintf(fp,"%6.3f %5.3f %6.3f %6.3f %5.2f %5.2f %5.2f %4.1f\n",
	    parser->pstat->cl_surprisal/c,
	    parser->pstat->cl_nonlex/c,
	    parser->pstat->cl_lex/c,
	    parser->pstat->cl_ambig/c,
	    parser->pstat->cl_open/c,
	    parser->pstat->cl_rernk/c,
	    parser->pstat->cl_toprr/c,
	    parser->pstat->cl_dsteps/c);
  if (type==2)
    fprintf(fp,"%6.3f %5.3f %6.3f %6.3f %5.2f %5.2f %5.2f %4.1f\n",
	    (parser->pstat->surprisal-parser->pstat->cl_surprisal)/c,
	    (parser->pstat->nonlex-parser->pstat->cl_nonlex)/c,
	    (parser->pstat->lex-parser->pstat->cl_lex)/c,
	    (parser->pstat->ambig-parser->pstat->cl_ambig)/c,
	    (parser->pstat->open-parser->pstat->cl_open)/c,
	    (parser->pstat->rernk-parser->pstat->cl_rernk)/c,
	    (parser->pstat->toprr-parser->pstat->cl_toprr)/c,
	    (parser->pstat->dsteps-parser->pstat->cl_dsteps)/c);
}

void show_sentnorm(FILE *fp,double derv,TDParPtr PC)
{
  float c=PC->parser->pstat->owds+PC->parser->pstat->cwds;
  fprintf(fp,"pfix sent norm (tot words %2d) ",
	  PC->parser->pstat->owds+PC->parser->pstat->cwds);
  show_norms(fp,c,PC->parser,0);
  if (PC->sfile!=NULL) {
    c=PC->parser->pstat->owds;
    fprintf(fp,"pfix sent norm (opn class %2d) ",PC->parser->pstat->owds);
    show_norms(fp,c,PC->parser,2);  
    c=PC->parser->pstat->cwds;
    fprintf(fp,"pfix sent norm (cls class %2d) ",PC->parser->pstat->cwds);
    show_norms(fp,c,PC->parser,1);  
  }
  c=derv;
  fprintf(fp,"pfix sent norm (steps %5.1f)  ",c);
  show_norms(fp,c,PC->parser,0);  
}

void show_prefixhdr(FILE *fp)
{
  fprintf(fp,"pfix header            ");
  fprintf(fp,"prefix ");
  fprintf(fp,"srprsl ");
  fprintf(fp,"SynSp  ");
  fprintf(fp,"LexSp  ");
  fprintf(fp,"ambig  ");
  fprintf(fp,"open ");
  fprintf(fp,"rernk ");
  fprintf(fp,"toprr ");
  fprintf(fp,"stps");
  fprintf(fp,"\n");
}
 
void show_aprefixhdr(FILE *fp)
{
  fprintf(fp,"afix header           ");
  fprintf(fp,"entropy  ");
  fprintf(fp," SynH");
  fprintf(fp,"  ESurp");
  fprintf(fp,"  ESynS");
  fprintf(fp,"  ELexS");
  fprintf(fp,"   rank ");
  fprintf(fp,"log(p*/p)");
  fprintf(fp,"\n");
}

int get_wtype(SLPtr lastring)
{
  if (lastring->word->cclass==0) return 0;
  if (lastring->word->ccount==0) return 2;
  return 1;
}

TreeHPtr dup_heap(TDParPtr PC, TreeHPtr ntrees, TreeHPtr *strees)
{
  int i=0;
  TreePtr xtree;
  TreeHPtr ctrees=NULL, thish, treeh;
  strees[0]=NULL;
  while (ntrees != NULL) {
    thish=ntrees;
    ntrees=perc_heap(ntrees);
    xtree=copy_tree(thish->tree->NT,thish->tree,PC->parser->features->numfeats);
    xtree->Score=thish->tree->Score;
    xtree->WScore=thish->tree->WScore;
    treeh=default_theap(xtree,xtree->Score,thish->preterm,thish->rank);
    strees[0]=insert_theap(treeh,strees[0]);
    ctrees=insert_theap(thish,ctrees);
  }
  return ctrees;
}

double logcomb(double log1, double log2)
{
  if (log1 < log2) return logcomb(log2,log1);
  return log1+log(exp(log2-log1)+1);
}

int get_hrank(int idx, double *lexent, double *expsyns, double *lograt, TDParPtr PC)
{
  double lnorm, thise, first=1, lbest, lthis, thiss;
  int cnt=1, i;
  if (idx > 0 && PC->wscores[0] > PC->wscores[idx]) cnt++;
  if (PC->wscores[0]>0) {
    lbest=lnorm=log(PC->wscores[0]);
    first=0;
  }
  if (idx==0) lthis=log(PC->wscores[0]);
  for (i=PC->parser->lexicon->numNTs; i < PC->parser->lexicon->numWDs-1; i++) {
    if (PC->wscores[i]>0) {
      if (first) {
	lbest=lnorm=log(PC->wscores[i]);
	first=0;
      }
      else {
	if (log(PC->wscores[i]) > lbest) lbest=log(PC->wscores[i]);
	lnorm=logcomb(lnorm,log(PC->wscores[i]));
      }
    }
    if (i==idx) {
      lthis=log(PC->wscores[i]);
      continue;
    }
    if (PC->wscores[i] > PC->wscores[idx]) cnt++;
  }
  lograt[0]=lbest-lthis;
  lexent[0]=0;
  if (PC->wscores[0]>0) {
    lexent[0]=PC->wscores[0];
    lexent[0]/=exp(lnorm);
    expsyns[0]=PC->sscores[0];
    expsyns[0]/=exp(lnorm);
    expsyns[0]=-log(expsyns[0]);
    expsyns[0]*=lexent[0];
    lexent[0]*=-log(lexent[0]);
  }
  for (i=PC->parser->lexicon->numNTs; i < PC->parser->lexicon->numWDs-1; i++) {
    if (PC->wscores[i]==0) continue;
    thise=PC->wscores[i];
    thise/=exp(lnorm);
    thiss = thise;
    thise*=-log(thise);    
    lexent[0]+=thise;
    thise=PC->sscores[i];
    thise/=exp(lnorm);
    thise=-log(thise);    
    thise*=thiss;
    expsyns[0]+=thise;
  }
  return cnt;
}

TreeHPtr parse_string(SLPtr string, TDParPtr PC, int *fails, FILE *fp)
{
  int cnt=1, wtype, i=0, rank=0;
  double lnorm=0, amb=0, openb, lcond=0, tcond=0, dsteps=0, plp=0, normp=0, 
    slp, snormp, lexent, expsyns, lograt;
  SLPtr lastring=NULL, stopstring;
  WLPtr holdword;
  TreePtr holdtree=NULL;
  TreeHPtr ctrees, ntrees=init_parse(PC,NULL), strees;
  if (PC->prefix) {
    reset_pstat(PC->parser->pstat);
    if (PC->allpossible) show_aprefixhdr(fp);
    else show_prefixhdr(fp);
  }
  if (PC->allpossible) {
    stopstring=default_string(PC->parser->lexicon->allword->label,NULL);
    stopstring->word=PC->parser->lexicon->allword;
  }
  while (string != NULL) {
    if (PC->allpossible) {
      holdword=string->word;
      string->word=PC->parser->lexicon->allword;
    }
    holdtree=get_holdtree(ntrees,holdtree);
    if (lastring != NULL) {
      amb=upnorms(ntrees,PC,&openb);
      ntrees=calc_condrank(ntrees,&tcond);
      dsteps=calc_derv(ntrees,1);
    }
    if (lastring!=NULL && PC->partial) 
      ntrees=show_partial(ntrees,PC,fp,lastring,cnt);
    PC->parser->lcond=ntrees->score;
    ctrees=update_look(ntrees,string,PC);
    if (lastring!=NULL && PC->prefix) {
      wtype=get_wtype(lastring);
      if (PC->allpossible)
	show_aprefix(fp,cnt,lastring->outlabel,(plp/normp)+log(normp),(slp/snormp)+log(snormp),rank,lexent,expsyns,lograt,wtype);
      else show_prefix(fp,cnt,lastring->outlabel,PC->parser,lnorm,amb,openb,
		       lcond,tcond,dsteps,wtype);
      lnorm=PC->parser->norm;
      lcond=PC->parser->lcond;
    }
    if (PC->allpossible) {
      PC->thresh*=2;
      if (lastring != NULL)
	ctrees=dup_heap(PC,ctrees,&strees);
    }
    ntrees=advance_parse(ctrees,string,PC);
    if (PC->allpossible) {
      ntrees=calc_allparse(ntrees,string,PC,&plp,&normp,&slp,&snormp);
      string->word=holdword;
      PC->thresh/=2;
      ntrees=filter_allparse(ntrees,string,PC);
      if (lastring != NULL)
	calc_stop_heap(strees,PC,&plp,&normp,&slp,&snormp);
      rank=get_hrank(string->word->windex,&lexent,&expsyns,&lograt,PC);
    }
    lastring=string;
    string=string->next;
    cnt++;
  }
  amb=upnorms(ntrees,PC,&openb);
  ntrees=calc_condrank(ntrees,&tcond);
  dsteps=calc_derv(ntrees,1);
  if (PC->partial) ntrees=show_partial(ntrees,PC,fp,lastring,cnt);
  if (ntrees != NULL)
    PC->parser->lcond=ntrees->score;
  if (PC->prefix) {
    wtype=get_wtype(lastring);
    if (PC->allpossible)
      show_aprefix(fp,cnt,lastring->outlabel,(plp/normp)+log(normp),(slp/snormp)+log(snormp),rank,lexent,expsyns,lograt,wtype);
    else show_prefix(fp,cnt,lastring->outlabel,PC->parser,lnorm,amb,openb,lcond,
		tcond,dsteps,wtype);
    lnorm=PC->parser->norm;
    lcond=PC->parser->lcond;
  }
  cnt++;
  if (PC->allpossible) {
    PC->thresh*=2;
    ntrees=dup_heap(PC,ntrees,&strees);
    strees=calc_allparse(strees,stopstring,PC,&plp,&normp,&slp,&snormp);
    free_theaps(strees);
  }
  ctrees=stop_heap(ntrees,PC,holdtree,fails,&plp,&normp,&slp,&snormp);
  if (PC->allpossible) {
    rank=get_hrank(0,&lexent,&expsyns,&lograt,PC);
    PC->thresh/=2;
  }

  amb=upnorms(ctrees,PC,&openb);
  if (PC->prefix) {
    dsteps=calc_derv(ctrees,1);
    if (PC->allpossible)
      show_aprefix(fp,cnt,"</s>",(plp/normp)+log(normp),(slp/snormp)+log(snormp),rank,lexent,expsyns,lograt,2);
    else {
      show_prefix(fp,cnt,"</s>",PC->parser,lnorm,amb,openb,lcond,tcond,dsteps,2);
      tcond=calc_derv(ctrees,0);
      show_sentnorm(fp,tcond,PC);
    }
    fprintf(fp,"Full parses for string:\n");
  }
  if (PC->allpossible) free_slist(stopstring);
  return ctrees;
}

void show_copytree(FILE *fp, TreePtr otree, int tabtree, int both)
{
  TreePtr tree=build_copytree(otree), punc=NULL;
  if (both) {
    if (tabtree) show_ttree(fp,tree,1);
    else show_tree(fp,tree,1);
    fprintf(fp,"\t\tDe-transformed tree:");
  }
  lc_detransform(tree);
  if (tabtree) show_ttree(fp,tree,1);
  else show_tree(fp,tree,1);
  free_trees(tree,NULL,0);
}

int countcands(TreeHPtr ctrees)
{
  int cnt=0;
  if (ctrees==NULL) return cnt;
  cnt++;
  cnt+=countcands(ctrees->next);
  cnt+=countcands(ctrees->right);
  cnt+=countcands(ctrees->left);
  return cnt;
}

void fill_cc(FILE *fp, LexPtr lex, char *readstr)
{
  int i, hash;
  WLPtr word, cw;
  while (fgets(readstr,MAXLABLEN,fp) != NULL) {
    i=0;
    while (!ws(readstr[i])) i++;
    readstr[i++]=0;
    word=find_word(readstr,lex,&hash);
    if (word != NULL) {
      word->cclass=1;
      if (readstr[i]=='0') word->ccount=0;
      else word->ccount=1;
    }
    else fprintf(stderr,"=== closed class word not found (%s)\n",readstr);
    if (readstr[0] > 96 && readstr[0] < 123) {
      readstr[0]-=32;
      cw=find_word(readstr,lex,&hash);
      if (cw != NULL) {
	cw->cclass=1;
	if (readstr[i]=='0') cw->ccount=0;
	else cw->ccount=1;
      }
    }
  }
}

void run_parser(char *ifile, TDParPtr PC)
{
  int i, k=0, fails=0, cands, hash;
  char readstr[MAXLABLEN];
  FILE *fp, *ofp;
  LexPtr lex;
  SLPtr string;
  TreeHPtr trees, thish;
  if (PC->verbose) fprintf(stderr,"loading... ");
  fp = gopenfile(ifile,"r");
  PC->parser=init_parser(fp,0);
  lex=PC->parser->lexicon;
  pclose(fp);
  if (PC->allpossible) {
    PC->maxpos=lex->numPTs;
    sprintf(readstr,"OOVALL");
    lex->allword=find_word(&readstr[0],lex,&hash);
    lex->allword=insert_word(&readstr[0],lex,hash,0,1);
    lex->allword->NTs=calloc(lex->numPTs,sizeof(float));
    lex->allword->oov=1;
    PC->wscores=calloc(PC->parser->lexicon->numWDs,sizeof(double));
    PC->sscores=calloc(PC->parser->lexicon->numWDs,sizeof(double));
  }
  if (PC->sfile!=NULL) {
    fp=fopen(PC->sfile,"r");
    fill_cc(fp,lex,&readstr[0]);
    fclose(fp);
  }
  if (PC->verbose) fprintf(stderr,"done\nParsing: ");
  if (PC->pfile == NULL) fp=stdin;
  else fp=fopen(PC->pfile,"r");
  if (PC->ofile == NULL) ofp=stdout;
  else ofp=fopen(PC->ofile,"w");
  while (fgets(readstr,MAXLABLEN,fp) != NULL) {
    string=get_string(&readstr[0],PC);
    trees=parse_string(string,PC,&fails,ofp);
    cands=countcands(trees);
    if (cands > PC->nbest) cands=PC->nbest;
    if (trees != NULL) {
      for (i=0; i < PC->nbest; i++) {
	thish=trees;
	trees=perc_heap(trees);
	thish->tree->WScore=calc_wscore(thish->tree);
	fprintf(ofp,"%d %3d %9.4f %9.4f %9.4f %9.4f\t",i+1,cands,thish->tree->Score,thish->tree->Score+PC->parser->norm,thish->tree->WScore,thish->tree->Score+PC->parser->norm-thish->tree->WScore);
	show_copytree(ofp,thish->tree,PC->tabtree,PC->partial);
	free_theap(thish);
	if (trees==NULL) break;
      }
      if (trees != NULL) free_theaps(trees);
    }
    else show_defstring(ofp,string,lex);
    fflush(ofp);
    free_slist(string);
    verbose_count(PC->verbose,++k,1000);
  } 
  if (PC->pfile != NULL) fclose(fp);
  if (PC->ofile != NULL) fclose(ofp);
  if (PC->verbose) fprintf(stderr," %d\nfailed: %d\n",k,fails);
}

void train_parser(TDTrPtr TC)
{
  int sents, ind, max;
  char readtree[MAXLABLEN];
  compile_feats(TC,&readtree[0]);
  sents=get_lexicon(TC,&readtree[0]);
  collect_counts(TC,&readtree[0],sents+1);
  ind=bucket_counts(TC);
  max=norm_feats(TC,ind,sents+1);
  update_pos(TC->parser);
  if (TC->thresh > 0) prune_feats(TC,max);
  write_parser(TC->parser,TC->ofile);
}

int makeunk(ParPtr par, int wd)
{
  int i=0;
  float lz=-log(i);
  LexPtr lex=par->lexicon;
  for (i=1; i < par->lexicon->numPTs; i++) {
    if (lex->WDs[wd]->NTs[i]>=lz) continue;
    return 0;
  }
  return 1;
}

void makeunks(ParPtr par)
{
  int i=0;
  float lz=-log(i);
  LexPtr lex=par->lexicon;
  for (i=lex->numNTs; i < lex->numWDs; i++) {
    if (lex->WDs[i]->label==par->hlabel) 
      lex->WDs[i]->label=lex->WDs[i]->hlabel;
    if (makeunk(par,i)) lex->WDs[i]->label=par->hlabel;
  }
}

int assign_batch(ParPtr par, int b, int k)
{
  if (k%par->features->scores->lrn->bsize != 0) return b;
  if (b >= par->features->scores->lrn->xvals) {
    par->features->scores->scores=par->features->scores->lrn->scores;
    par->features->scores->pscores=par->features->scores->lrn->ascores;
    par->features->scores->ppsc=par->features->scores->lrn->apsc;
  }
  else {
    par->features->scores->scores=par->features->scores->lrn->fscores[b]->rsc;
    par->features->scores->pscores=par->features->scores->lrn->pscores;
    par->features->scores->ppsc=par->features->scores->lrn->ppsc;
  }
  update_wpos(par);
  makeunks(par);
  return b+1;
}

FUDPtr addfud(CFDataPtr cfeat,int lastidx,float bcost,FUDPtr ifud, float upd)
{
  FUDPtr thisfud, lastfud=NULL, fud=ifud;
  while (fud != NULL && fud->idx < cfeat->idx) {
    lastfud=fud;
    fud=fud->next;
  }
  if (fud != NULL && fud->idx == cfeat->idx) {
    fud->count += upd;
    return ifud;
  }
  thisfud=malloc(sizeof(struct FeatUpd));
  thisfud->idx=cfeat->idx;
  thisfud->bcost=bcost;
  thisfud->lastidx=lastidx;
  thisfud->count=upd;
  if (lastfud==NULL) {
    thisfud->next=ifud;
    return thisfud;
  }
  thisfud->next=fud;
  lastfud->next=thisfud;
  return ifud;
}

FUDPtr update_pschema(FtrPtr feats, int sch, int clbl, TreePtr tree, FUDPtr fud, float upd, LexPtr lex)
{
  int j, k, st, lb, hash, lastidx=-1;
  TWSchPtr fschema=feats->schema[sch];
  FeatPtr currfeat=NULL, prefix;
  ScorePtr scr=feats->scores;
  FHashPtr currhash, holdhash;
  CFDataPtr cfeat;
  for (j=1; j < fschema->numfuncts; j++) {
    prefix=currfeat;
    k=fschema->functs[j]->findex;
    if (prefix==NULL) st=0;
    else st=prefix->dstate;
    lb=tree->featarr[k];
    holdhash=currhash=find_feat(st,lb,j,sch,feats,&hash);
    if (currhash == NULL) {
      holdhash=currhash=insert_feat(st,lb,j,sch,feats,hash,1);
      currhash->feat->data->prefix=prefix;
      currhash->feat->data->bscore=0.0;
    }
    currfeat=currhash->feat;
    currhash=find_feat(currfeat->dstate,clbl,j+1,sch+1,feats,&hash);
    if (currhash == NULL) {
      if (lastidx<0) return fud;
      currhash=insert_feat(currfeat->dstate,clbl,j+1,sch+1,feats,hash,0);
      currhash->feat->dstate=currfeat->data->numc;
      insert_cfeat(feats,currfeat->data,clbl,0);
    }
    cfeat=currfeat->data->cfeat[currhash->feat->dstate];
    fud=addfud(cfeat,lastidx,currfeat->data->bscore,fud,upd);
    lastidx=cfeat->idx;
  }
  return fud;
}

FUDPtr get_fuds(TreePtr tree, ParPtr par, FUDPtr ifud, float upd, TreePtr ltree)
{
  int i,clbl;
  FUDPtr fud=ifud;
  TWSchPtr fschema;
  FtrPtr feats=par->features;
  if (tree==NULL || tree->featarr==NULL || tree==ltree) return fud;
  fud=get_fuds(tree->prev,par,ifud,upd,ltree);
  for (i=0; i < feats->numschema; i++) {
    fschema = feats->schema[i];
    clbl=tree->featarr[fschema->functs[0]->findex];
    if (clbl < 0) continue;
    fud=update_pschema(feats,i,clbl,tree,fud,upd,par->lexicon);
  }
  return fud;
}

int re_lrn(ScorePtr scr, int ind)
{
  int i, asc=0, fsc=scr->lrn->xvals;
  if (scr->lrn->ascores==scr->pscores) asc=1;
  scr->lrn->scores=defscores(scr->lrn->scores,ind,scr->numalloc);
  scr->lrn->ascores=defscores(scr->lrn->ascores,ind,scr->numalloc);
  scr->lrn->pscores=defscores(scr->lrn->pscores,ind,scr->numalloc);
  scr->lrn->lastupd=deflupd(scr->lrn->lastupd,ind,scr->numalloc);
  for (i=0; i < scr->lrn->xvals; i++) {
    if (scr->lrn->fscores[i]->rsc==scr->scores) fsc=i;
    scr->lrn->fscores[i]->rsc=
      defscores(scr->lrn->fscores[i]->rsc,ind,scr->numalloc);
  }
  scr->lrn->backs=realloc(scr->lrn->backs,ind*sizeof(BORPtr));
  for (i=scr->numalloc; i < ind; i++) scr->lrn->backs[i]=NULL;
  if (asc) scr->pscores=scr->lrn->ascores;
  else scr->pscores=scr->lrn->pscores;
  if (fsc==scr->lrn->xvals) scr->scores=scr->lrn->scores;
  else scr->scores=scr->lrn->fscores[fsc]->rsc;
  scr->numalloc=ind;
  return fsc;
}

void newstate(ScorePtr scr, FUDPtr fud)
{
  int i;
  BORPtr back=malloc(sizeof(struct BORev));
  back->state=fud->idx;
  if (fud->lastidx < 0) fprintf(stderr,"oops, lastidx neg: %d\n",fud->lastidx);
  back->next=scr->lrn->backs[fud->lastidx];
  scr->lrn->backs[fud->lastidx]=back;
  scr->lrn->scores[fud->idx]=scr->lrn->scores[fud->lastidx]+fud->bcost;
  scr->lrn->ascores[fud->idx]=scr->lrn->ascores[fud->lastidx];
  scr->lrn->lastupd[fud->idx]=scr->lrn->lastupd[fud->lastidx];
  scr->lrn->pscores[fud->idx]=scr->lrn->pscores[fud->lastidx];
  for (i=0; i < scr->lrn->xvals; i++)
    scr->lrn->fscores[i]->rsc[fud->idx]=
      scr->lrn->fscores[i]->rsc[fud->lastidx]+fud->bcost;
}

void update_avg(LearnPtr lrn, int idx, float upd)
{
  double sum=0, diff=0;
  if (lrn->lastupd[idx]>0 && lrn->ascores[idx] != 0.0) {
    sum=lrn->ascores[idx];
    sum*=lrn->lastupd[idx];
  }
  if (lrn->lastupd[idx] < lrn->iter) {
    diff=lrn->pscores[idx];
    diff*=(lrn->iter-lrn->lastupd[idx]);
    sum+=diff;
    lrn->lastupd[idx]=lrn->iter;
  }
  sum+=upd;
  sum/=lrn->lastupd[idx];
  lrn->ascores[idx]=sum;
}

void upercep(ScorePtr scr, int idx, float upd)
{
  if (upd==0) return;
  update_avg(scr->lrn,idx,upd);
  scr->lrn->pscores[idx]+=upd;
}

void upback(ScorePtr scr,BORPtr back,float upd)
{
  BORPtr nback;
  while (back != NULL) {
    upercep(scr,back->state,upd);
    upback(scr,scr->lrn->backs[back->state],upd);
    back=back->next;
  }
}

void upfud(FUDPtr fud, ScorePtr scr)
{
  upercep(scr,fud->idx,fud->count);
  upback(scr,scr->lrn->backs[fud->idx],fud->count);
}

void update_pfeats(TreePtr gtree, TreePtr btree, ParPtr par, int batch)
{
  ScorePtr scr=par->features->scores;
  TreePtr bgtree=btree, gbtree=NULL;
  int ind=scr->numind, b;
  float upd=par->step;
  FUDPtr fud=NULL, tfud;
  scr->lrn->iter++;
  if (gtree==NULL || btree==NULL) return;
  while (bgtree != NULL && bgtree->gold==NULL) bgtree=bgtree->prev;
  while (bgtree != NULL && !bgtree->NT->terminal) bgtree=bgtree->prev;
  if (bgtree != NULL) gbtree=bgtree->prev->gold->gtree;
  fud=get_fuds(btree,par,fud,upd,bgtree);
  upd*=-1;
  tfud=fud=get_fuds(gtree,par,fud,upd,gbtree);
  if (scr->numind >= scr->numalloc) {
    b=re_lrn(scr,scr->numind+1000);
    if (b != batch-1) {
      fprintf(stderr,"oops, wrong batch %d %d\n",b,batch-1);
      exit(1);
    }
  }
  while (tfud != NULL) {
    if (tfud->idx >= ind) newstate(scr,tfud);
    tfud=tfud->next;
  }
  while (fud != NULL) {
    upfud(fud,scr);
    fud=fud->next;
  }
}

void update_pmodels(ParPtr par, GLPtr gold, int batch)
{
  if (gold==NULL) return;
  update_pfeats(gold->gtree,gold->best,par, batch);
  if (gold->best != NULL) {
    gold->best->points--;
    free_trees(gold->best,NULL,0);
  }
}

GLPtr percep_string(SLPtr string, TDParPtr PC, int *fails, GLPtr gold)
{
  SLPtr lastring=NULL;
  GLPtr rgold=gold;
  double lpr, normp, spr, snormp;
  TreeHPtr ctrees, ntrees=init_parse(PC,gold);
  while (string != NULL) {
    ctrees=update_look(ntrees,string,PC);
    ntrees=advance_parse(ctrees,string,PC);
    lastring=string;
    string=string->next;
    rgold=get_rgold(ntrees,rgold,0);
  }
  ctrees=stop_heap(ntrees,PC,NULL,fails,&lpr,&normp,&spr,&snormp);
  if (rgold != NULL && rgold->best==NULL && ctrees != NULL &&
      ctrees->tree->prev->gold==NULL) {
    while (rgold->next != NULL) rgold=rgold->next;
    rgold->best=ctrees->tree;
    rgold->best->points++;
  }
  if (ctrees != NULL) {
    show_copytree(stdout,ctrees->tree,PC->tabtree,PC->partial);
    free_theaps(ctrees);
  }
  else show_defstring(stdout,string,PC->parser->lexicon);
  return rgold;
}

void update_avgs(ScorePtr scr)
{
  int i;
  for (i=0; i < scr->numind; i++) {
    if (scr->lrn->lastupd[i]==scr->lrn->iter) continue;
    update_avg(scr->lrn,i,0.0);
  }
}

void parse_percep(TDParPtr PC)
{
  int k=0, fails=0, batch=0;
  char readstr[MAXLABLEN];
  FILE *fp;
  SLPtr string;
  GLPtr gold, ugold;
  TreeHPtr trees;
  TreePtr tree, btree;
  if (PC->verbose) fprintf(stderr,"Parsing: ");
  if (PC->pfile == NULL) fp=stdin;
  else fp=fopen(PC->pfile,"r");
  batch=assign_batch(PC->parser,batch,k);
  while (fgets(readstr,MAXLABLEN,fp) != NULL) {
    tree=build_tree(readstr,PC->parser,1);
    string=get_tstring(PC->parser,tree);
    gold=get_gold(tree,NULL);
    ugold=percep_string(string,PC,&fails,gold);
    update_pmodels(PC->parser,ugold,batch);
    free_gold(gold);
    free_trees(tree,NULL,1);
    verbose_count(PC->verbose,++k,1000);
    batch=assign_batch(PC->parser,batch,k);
  } 
  if (PC->pfile != NULL) fclose(fp);
  update_avgs(PC->parser->features->scores);
  assign_batch(PC->parser,PC->parser->features->scores->lrn->xvals,0);
  if (PC->verbose) fprintf(stderr," %d\nfailed: %d\n",k,fails);
}

void train_percep(char *ifile, TDParPtr PC, float step)
{
  FILE *fp;
  if (PC->verbose) fprintf(stderr,"loading... ");
  fp = gopenfile(ifile,"r");
  PC->parser=init_parser(fp,1);
  PC->parser->step=step;
  pclose(fp);
  if (PC->verbose) fprintf(stderr,"done\n");
  PC->parser->hlabel=malloc(MAXLABLEN*sizeof(char));
  sprintf(PC->parser->hlabel,"NoMatchZAAARG01@#$5235");
  fill_borev(PC->parser->features);
  parse_percep(PC);
  write_parser(PC->parser,PC->ofile);
}

