#!/usr/bin/awk -f
BEGIN{}
{rec=rec sep $0; sep=ORS} 
{fn = substr(FILENAME,1,length(FILENAME)-3)"_"length($0)"nt.txt"}
!(NR%2){print rec > fn; rec=sep=""} 
END{}

