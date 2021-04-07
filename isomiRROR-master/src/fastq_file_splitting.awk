#!/usr/bin/awk -f
BEGIN{}
{rec=rec sep $0; sep=ORS} 
!(NR%4){print rec > fn; rec=sep=""} 
NR%4==2{fn = substr(FILENAME,1,length(FILENAME)-6)"_"length($0)"nt.txt"}
END{}

