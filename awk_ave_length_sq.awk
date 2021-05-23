BEGIN{
	ORS = " ";
}
{
	total[0,FNR] = $1;
	total[1,FNR] += ($4+$7)*($4+$7);
	total[2,FNR] += ($10+$11);
}
END{
	for(i=1 ; i<=FNR ; i++) { 
		for(j=0; j<3; j++) printf "%f ", total[j,i];
		printf "%f ", 1.0;
		printf "\n";
		}
}

