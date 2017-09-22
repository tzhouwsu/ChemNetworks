#!/bin/sh

PRE=Input-2.2.water1.xyz.water1.xyz

O1=${PRE}.Graph
O2=${PRE}.GraphGeod

rm -f $O1 $O2
rm -f ${PRE}.*.Graph ${PRE}.*.GraphGeod

#../../simple_driver/simple_driver Input-2.2 water1.xyz
#diff -sq $O1 orig/$O1
#diff -sq $O2 orig/$O2


#../../simple_driver/simple_driver Input-2.2 water1.xyz -new
../../simple_driver2/simple_driver Input-2.2 water1.xyz -new

I=0
while [ $I -lt 5 ]
do
  diff -sq ${PRE}.$I.Graph     orig/${O1}
  diff -sq ${PRE}.$I.GraphGeod orig/${O2}
  I=`expr $I + 1`
done
