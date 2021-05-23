BEGIN{
	ORS = " ";
}
{
	if(FNR < 94) {
	for(j=1 ; j<=NF ; j++) total[j,FNR] += $j;
	tot += 1;
}
}
END{
	tot = tot/93;
	for(i=1 ; i<=93 ; i++) { 
		for(j=1; j<=NF; j++) printf "%f ", total[j,i]/tot;
		printf "\n";
		}
}

