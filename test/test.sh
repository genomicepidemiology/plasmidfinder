#!/bin/bash 
plasmidfinder.py -i /test/test.fsa -o /test/ -mp blastn -x -q
file=/test/results_tab.tsv
DIFF=$(diff $file /test/test_results.tsv)
if [ "$DIFF" == "" ] && [ -s $file ] ;
   then     
   echo "TEST SUCCEEDED"; 
else
   echo "TEST FAILED";
fi
