#!/bin/sh

paste - - < $1 |
awk '{dup[$2]=dup[$2] ? dup[$2] " duplicate of " $1 : $1} END {for (x in dup) print dup[x]\n x}' |
grep 'duplicate' > multireads.txt;

paste - - < $1 |
awk -F '[_\t]' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' |
awk -F '\t' 'function cp(){a[$10]=$0;b[$10]=($3)^2;c[$10]=($7)^2}{f=$10 in a;d[$10]++}
!f{o[++i]=$10; cp();next}f && b[$10]>($3)^2{cp();next}f && b[$10]==($3)^2 && ($7)^2<c[$10]
{cp();next}END{for(i=1; i in o; i++)print a[o[i]],(d[o[i]]>1?"multi":"unique")}' |
awk -F '[ \t]' '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8"_"$9"_"$11"\n"$10}' 
> ${1%%.fa}_isomir.fa;
