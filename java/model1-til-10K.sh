#!/bin/sh

#java -Xmx1500m -cp bin cs224n.assignments.WordAlignmentTester -path ../../pa2-data/ -model model3 -data test -verbose > ../results/Model3/sentences-0

java -Xmx1500m -cp bin cs224n.assignments.WordAlignmentTester -path ../../pa2-data/ -model model3 -data test -sentences 10 -verbose  > ../results/Model3/10.txt

java -Xmx1500m -cp bin cs224n.assignments.WordAlignmentTester -path ../../pa2-data/ -model model3 -data test -sentences 1000 -verbose  > ../results/Model3/1000.txt

java -Xmx1500m -cp bin cs224n.assignments.WordAlignmentTester -path ../../pa2-data/ -model model3 -data test -sentences 10000 -verbose  > ../results/Model3/10000.txt


