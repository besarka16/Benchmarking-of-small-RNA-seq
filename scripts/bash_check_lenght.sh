
#counts lines longer than 28

grep -A1 "^[@]" Qiag__MirXplore_1_SE.fastq.trim | grep -v "^[@]" | sed -n '2~2!p' | grep '.\{29\}' | wc -l

#counts lines shorter than 16

grep -A1 "^[@]" Qiag__MirXplore_1_SE.fastq.trim | grep -v "^[@]" | sed -n '2~2!p' | grep -E '^.{1,15}$' | wc -l

