#! /bin/sh

TREESURGEONPATH=..
export CLASSPATH=$TREESURGEONPATH/stanford-tregex.jar:$CLASSPATH

# removes empty nodes, normalizes root category to TOP and ensures unary TOP 
# production; finds and removes function tags, coreference indices and 
# ambiguous categories; prunes (CODE trees (used in Switchboard for speaker ID)
java -mx100m edu.stanford.nlp.trees.tregex.tsurgeon.Tsurgeon -s \
     -treeFile $1 pruneEMPTY pruneCODE relabelTOP relabelBINTOP | \
  grep "^(" | sed 's/(^/(/g' | \
  sed 's/^(BINTOP.*$/(TOP &)/g' | sed 's/BINTOP/S/g' | \
  sed 's/^([^T].*$/(TOP &)/g' | \
  sed 's/([^-=^| ][^-=^| ]*[-=^|][^ ]*/&~/g' | \
  sed 's/[-=^|][^ ]*~//g'
