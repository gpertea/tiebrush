#!/bin/env bash
/bin/rm -f tst_* t?/tst_*

valgrind --track-origins=yes --leak-check=full --show-reachable=yes ../tiebrush \
 -o t1/tst_t1.bam t1/t1s[0-9].bam |& tee tst_valgrind_t1.log

valgrind --track-origins=yes --leak-check=full --show-reachable=yes ../tiebrush \
 -o tst_t12.bam t1/t1.bam t2/t2.bam |& tee tst_valgrind_t12.log

valgrind --track-origins=yes --leak-check=full --show-reachable=yes ../tiecov -s t1/tst_t1.sample \
  -c t1/tst_t1.coverage -j t1/tst_t1.junctions t1/t1.bam |& tee tst_valgrind_tiecov.log 

