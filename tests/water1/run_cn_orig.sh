#!/bin/sh

PRE=Input-2.2.water1.xyz.water1.xyz

O1=${PRE}.Graph
O2=${PRE}.GraphGeod

rm -f $O1 $O2
rm -f ${PRE}.0.Graph ${PRE}.0.GraphGeod

#../../src/ChemNetworks-2.2.exe Input-2.2 water1.xyz
#diff -sq $O1 orig/$O1
#diff -sq $O2 orig/$O2


../../src/ChemNetworks-2.2.exe Input-2.2 water1.xyz
diff -sq ${PRE}.Graph     orig/${O1}
diff -sq ${PRE}.GraphGeod orig/${O2}
