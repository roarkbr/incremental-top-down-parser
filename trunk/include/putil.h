#ifndef PUTILHDR
#define PUTILHDR

#include "pstruct.h"

void failproc(char *msg);
int ws(char n);
int nonalpha(char *str);
int noalpha(char *str);
TDTrPtr TDTrainDefault();
TDParPtr TDParseDefault();
LexPtr LexDefault(int numWH, int numWD, int stop);
ParPtr ParserDefault(int numWH, int numFH, int numWD, int stop);
void LexFree(LexPtr lex);
int prime(int i);
WLPtr find_word(char *word, LexPtr lex, int *hash);
WLPtr insert_word(char *word,LexPtr lex,int hash,int children,int terminal);
WLPtr find_nword(char *word, LexPtr lex, int *hash, WLPtr fword, char *flabel, WLPtr lword, int first, int climit, int oov);
WLPtr find_zword(char *word, LexPtr lex, int *hash, WLPtr fword, char *flabel, WLPtr lword);
SLPtr default_string(char *label, SLPtr lastring);
WLPtr find_nsword(char *word, LexPtr lex, int *hash, char *flabel, int first);
void add_slcprods(TDTrPtr TC,char *label);
void add_slcfile(TDTrPtr TC, char *filename);
void make_slcentry(WLPtr pword, WLPtr lword, WLPtr word);
void add_featschemas(TDTrPtr TC,char *label);
void add_featfile(TDTrPtr TC, char *filename);
void wordtoken(char *readtree, char *readterm, int *i, int len, int STOPs);
void show_tree(FILE *fp, TreePtr tree, int first);
int show_ctree(FILE *fp, TreePtr tree, int first);
TreePtr build_tree(char *readtree, ParPtr par, int oov);
void free_slist(SLPtr slist);
void free_tree(TreePtr tree);
void free_trees(TreePtr tree, SLPtr tstring, int dostr);
void show_lcats(TDTrPtr TC);
void lc_detransform(TreePtr tree);
FHashPtr find_feat(int state, int label, int level, int schema, FtrPtr features, int *hash);
FHashPtr insert_feat(int state, int label, int level, int schema, FtrPtr features, int hash, int afeat);
void insert_cfeat(FtrPtr feats, FDataPtr fdata, int clbl, int lam);
void update_pos(ParPtr par);
void update_wpos(ParPtr par);
int norm_feats(TDTrPtr TC, int ind, int sents);
void prune_feats(TDTrPtr TC, int max);
TreeHPtr insert_theap(TreeHPtr tree, TreeHPtr trees);
TreeHPtr insert_theapns(TreeHPtr tree, TreeHPtr trees);
TreeHPtr default_theap(TreePtr tree, float score, int pt, int rank);
TreePtr default_tree(WLPtr nt, TreePtr lastree);
void default_link(WLPtr nt, TreePtr lastree, TreePtr tree);
void free_theaps(TreeHPtr trees);
void free_theap(TreeHPtr trees);
TreeHPtr init_parse(TDParPtr PC,GLPtr gold);
void get_ntprobs(int lb, int siblb, TDParPtr PC);
TreePtr get_ptree(TreePtr tree);
TreePtr copy_tree(WLPtr nt, TreePtr ntree, int nf);
void fill_farrays(TreePtr tree, TreePtr otree, FtrPtr feats, int all);
void fill_lclfarrays(TreePtr tree, TreePtr otree, FtrPtr feats);
TreePtr build_copytree(TreePtr btree);
TreePtr get_hnode(TreePtr tree);
void show_defstring(FILE *fp, SLPtr string, LexPtr lex);
int bucket_counts(TDTrPtr TC);
void verbose_count(int verbose, int k, int mod);
ParPtr init_parser(FILE *fp, int lrn);
TreePtr gettwop(TreePtr ftree, TreePtr tree, char c, char lc);
void fill_tnull(TreeHPtr tree);
TreeHPtr perc_heap(TreeHPtr ctrees);
void zapentry(FtrPtr feats, int i, FHashPtr last, FHashPtr this, CFDataPtr *cfeats, int xvals);
float getboff(FtrPtr feats, FeatPtr feat, int lbl, float bcost, int bottom, float *p);
SLPtr get_tstring(ParPtr par, TreePtr tree);
void fill_tree(TreePtr tree, ParPtr par);
GLPtr get_rgold(TreeHPtr ntrees,GLPtr gold,int last);
GLPtr get_gold(TreePtr tree, GLPtr pgold);
float getcost(ScorePtr scr, FeatPtr feat, FHashPtr fhash, float bcost, float *p);
float *defscores(float *scores, int numind, int numalloc);
int *deflupd(int *lupd, int numind, int numalloc);
TreePtr goldcand(TreeHPtr ntrees,int last);
void putpuncback(TreePtr tree, TreePtr punc);
TreePtr punctnorm(TreePtr tree);
void reset_pstat(PStats pstat);

#endif /* PUTILHDR */
