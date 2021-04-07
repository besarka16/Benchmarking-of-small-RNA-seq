#!/usr/bin/awk -f
BEGIN{}
{
	if ($14 =="NM:i:0") {
		print $3"_0mm_0_X_Seed_"substr($10,2,6)"_M8_"substr($10,8,1)"_End3_"substr($10,length($10)-6,7)"_length_"length($10);
	}

	else if ($14 != " " && substr($14,6,1) ~ /^[0-9]/){
		poly_num=0;
		iso_str="";
		split($13,a,":");
		tmp_pos=0;
		for(i=1; i<=length(a[3]);i++){
			if (substr(a[3],i,1) ~ /^[A-Z]/){
				pos_array[poly_num]=i;
				if(poly_num==0){
					tmp_pos+=poly_num+substr(a[3],1,i-1)+1;
				}else{
					tmp_pos+=substr(a[3],int(pos_array[poly_num-1]+1),int((pos_array[poly_num]-pos_array[poly_num-1]-1)))+1;
				}
				iso_str=iso_str""tmp_pos":"substr(a[3],i,1)">"substr($10,tmp_pos,1)",";
				poly_num++;
		}
	}
	print $3"_Poly_"poly_num"_"substr(iso_str,1,length(iso_str)-1)"_Seed_"substr($10,2,6)"_M8_"substr($10,8,1)"_End3_"substr($10,length($10)-6,7)"_length_"length($10);
	}


}
END{}
