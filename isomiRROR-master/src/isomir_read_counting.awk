#!/usr/bin/awk -f
BEGIN{}
{dups[$0]++} 
END{for (num in dups) {print num,dups[num]}}

