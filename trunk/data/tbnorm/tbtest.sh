#! /bin/sh

tmpdata=${TMPDIR:-/tmp}

cat $1 | sed 's/  *(/~(/g' | sed 's/[^ ]*~//g' | sed 's/  *[^( ]*$//g' | \
 sed 's/ *(/\
/g' | grep . | sort -u >$tmpdata/pre-terms.txt
cat $1 | sed 's/ *(/\
/g' | grep -v " " | grep . | sort -u >$tmpdata/non-terms.txt
cat $tmpdata/pre-terms.txt | while read i; do \
  grep -xF $i $tmpdata/non-terms.txt; done >$tmpdata/matched.txt
wc $tmpdata/matched.txt | while read a b; do 
 if [ $a -lt 1 ]; then echo "POS-tags and non-POS-tag non-terminals are disjoint";
 else echo "The following are both POS-tags and non-POS-tag non-terminals:"; 
 cat $tmpdata/matched.txt; fi; done
rm $tmpdata/matched.txt $tmpdata/non-terms.txt $tmpdata/pre-terms.txt
