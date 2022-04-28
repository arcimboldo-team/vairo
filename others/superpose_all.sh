#!/bin/bash


prog=superpose

# Find do seems to be the best construct to check directories, understands spaces
find atzr-* -maxdepth 0 -type d \( ! -name . \) | while read dir; do
N="$dir";

echo 
echo $N'********************************************'
echo 

t1=$dir'/template.pdb'
t2=$dir'/experimental_template.pdb'
t3=$dir'/experimental_ranked_0.pdb'
t4=$dir'/experimental_ranked_1.pdb'

p0=$dir'/ranked_0.pdb'
p1=$dir'/ranked_1.pdb'

echo $p0 vs $t1
$prog $p0 $t1  -o $dir/fit_${N}_ranked0_template.pdb | grep -A1 "r.m.s.d";
echo $p1 vs $t1
$prog $p1 $t1  -o $dir/fit_${N}_ranked1_template.pdb | grep -A1 "r.m.s.d";

echo $p0 vs $t2
$prog $p0 $t2  -o $dir/fit_${N}_ranked0_experimental_template.pdb | grep -A1 "r.m.s.d";
echo $p1 vs $t2
$prog $p1 $t2  -o $dir/fit_${N}_ranked1_experimental_template.pdb | grep -A1 "r.m.s.d";

echo $p0 vs $t3
$prog $p0 $t3  -o $dir/fit_${N}_ranked0_experimental_ranked0.pdb | grep -A1 "r.m.s.d";
echo $p1 vs $t3
$prog $p1 $t3  -o $dir/fit_${N}_ranked1_experimental_ranked0.pdb | grep -A1 "r.m.s.d";

echo $p0 vs $t4
$prog $p0 $t4  -o $dir/fit_${N}_ranked0_experimental_ranked1.pdb | grep -A1 "r.m.s.d";
echo $p1 vs $t4
$prog $p1 $t4  -o $dir/fit_${N}_ranked1_experimental_ranked1.pdb | grep -A1 "r.m.s.d";


done
exit 0

# | grep -A2 "r.m.s.d"
