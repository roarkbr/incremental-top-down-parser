The incremental top-down parser is a statistical syntactic parser that processes input strings from left-to-right, producing partial derivations in a top-down manner, using beam search as detailed in [Roark (2001)](http://www.aclweb.org/anthology/J/J01/J01-2004.pdf) and [Roark (2004)](http://www.lanzaroark.org/brian-roark/publications/docs/NLE04.pdf?attredirects=0). It can output parser state statistics of utility for psycholinguistic studies, as detailed in [Roark et al.  (2009)](http://www.aclweb.org/anthology-new/D/D09/D09-1034.pdf).

Note that this code is research code, provided to the community without promise of additional technical support. Many researchers have used this code in their research, and were able to successfully compile and use the code on linux and Mac OS X systems. If you follow the instructions below, you should be able to use it as well, but if you encounter problems you are pretty much on your own.

In addition to the parser and the model training code, a model trained on sections 2-21 of the Penn WSJ treebank has been included.  For those who want to train their own models, there are some scripts for converting Penn Treebank files from the 'pretty-print' representation provided by the LDC to the format expected by the model training routine: one parse per line, with no function tags or empty nodes, etc. You are free to change your non-terminal labels (e.g., changing VBs to AUXs, etc.) but the code assumes that POS-tags and non-POS-tag non-terminals are disjoint.

## Quick start usage ##

Download the tar file, untar, change to the incremental-top-down-parser directory and refer to the README for installation instructions.  The tar file does not include the pre-trained model, which must be downloaded separately or you must train your own model with a normalized treebank.

### Parsing ###

To get information on usage, type the following on the command line after installation:

`bin/tdparse -?`

The following is standard usage for parsing a file of strings 'mystrings.txt':

`bin/tdparse -v -F parse.output parse.model.wsj.slc.p05 mystrings.txt`

The strings file 'mystrings.txt' should have one string per line.  See README file for further information.

### Training models ###

To train a model, you must first prepare a training treebank.  Scripts have been provided to assist in converting Penn Treebank style parses from their raw format as distributed by the LDC to the format expected by the model training code, see README for specifics.

To get information on usage, type the following on the command line after installation:

`bin/tdptrain -?`

The following is standard usage for training a model from a given treebank 'treebank.norm.txt':

`bin/tdptrain -t0.05 -p -s"(SINV(S" -v -F parse.model.slc.p05 -m data/tree.functs.sfile treebank.norm.txt`

See README file for further information.

## Citing ##

The basic citation for the parsing approach is:

Brian Roark. 2001. [Probabilistic top-down parsing and language modeling](http://www.aclweb.org/anthology/J/J01/J01-2004.pdf). In Computational Linguistics, 27(2), pages 249-276.

If using the parser to derive parser state statistics for psycholinguistic modeling, it's also a good idea to cite:

Brian Roark, Asaf Bachrach, Carlos Cardenas and Christophe Pallier. 2009. [Deriving lexical and syntactic expectation-based measures for psycholinguistic modeling via incremental top-down parsing](http://www.aclweb.org/anthology-new/D/D09/D09-1034.pdf). In Proceedings of the Conference on Empirical Methods in Natural Language Processing (EMNLP), pp. 324-333.

## References ##

Brian Roark. 2011. [Expected surprisal and entropy](http://www.lanzaroark.org/brian-roark/publications/docs/techrpt-CSLU-11-004.pdf?attredirects=0). Technical Report #CSLU-11-004, Center for Spoken Language Processing, Oregon Health & Science University.

Brian Roark, Asaf Bachrach, Carlos Cardenas and Christophe Pallier. 2009. [Deriving lexical and syntactic expectation-based measures for psycholinguistic modeling via incremental top-down parsing](http://www.aclweb.org/anthology-new/D/D09/D09-1034.pdf). In Proceedings of the Conference on Empirical Methods in Natural Language Processing (EMNLP), pp. 324-333.

Michael Collins and Brian Roark. 2004. [Incremental Parsing with the Perceptron Algorithm](http://www.aclweb.org/anthology/P/P04/P04-1015.pdf). In Proceedings of the 42nd Annual Meeting of the Association for Computational Linguistics (ACL), pages 111-118.

Brian Roark. 2004. [Robust garden path parsing](http://www.cslu.ogi.edu/people/roark/NLE04.pdf). Natural Language Engineering, 10(1), pages 1-24.

Brian Roark. 2001. [Robust Probabilistic Predictive Syntactic Processing: Motivations, Models, and Applications](http://arxiv.org/abs/cs/0105019). Ph.D. Thesis, Department of Cognitive and Linguistic Sciences, Brown University.

Brian Roark. 2001. [Probabilistic top-down parsing and language modeling](http://www.aclweb.org/anthology/J/J01/J01-2004.pdf). In Computational Linguistics, 27(2), pages 249-276.

Mark Johnson and Brian Roark. 2000. [Compact non-left-recursive grammars using the selective left-corner transform and factoring](http://www.aclweb.org/anthology/C/C00/C00-1052.pdf). In Proceedings of the 18th International Conference on Computational Linguistics (COLING), pages 355-361.

Brian Roark and Mark Johnson. 1999. [Efficient probabilistic top-down and left-corner parsing](http://www.aclweb.org/anthology/P/P99/P99-1054.pdf). In Proceedings of the 37th Annual Meeting of the Association for Computational Linguistics, pages 421-428.