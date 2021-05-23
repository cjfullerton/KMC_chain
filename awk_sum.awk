BEGIN{
	ORS = " ";
}
{
	for(j=1 ; j<=NF ; j++) total[j,FNR] += $j;
}
END{
	for(i=1 ; i<=FNR ; i++) { 
		for(j=1; j<=NF; j++) printf "%f ", total[j,i];
		printf "\n";
	}
}

