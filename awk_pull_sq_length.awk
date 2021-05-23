{
	if($3 > 0.00) {
		printf "%f ", ($2/$3)/1000.0;
	}
	else {
		printf "%f ", 0.00;
	}
}
