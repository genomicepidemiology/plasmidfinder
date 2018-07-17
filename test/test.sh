#!/bin/bash 
plasmidfinder.py -i /test/test.fsa -o /test/ -mp blastn -x
file=/test/results_tab.txt
DIFF=$(diff $file /test/test_results.txt)
if [ "$DIFF" == "" ] && [ -s $file ] ;
   then     
   echo "TEST SUCCEEDED"; 
else
   echo "TEST FAILED";
fi
