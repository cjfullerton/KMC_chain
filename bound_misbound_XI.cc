#include "bound_misbound_XI.h"


vector <double> bound_misbound_XI::return_rate_vector() {
	
	return rate_vector;

}

vector <int> bound_misbound_XI::return_event(int event_no) { //events now much more complicated...
	
	return event_vector[event_no];

}


void bound_misbound_XI::calculate_rate_vector() {

	double coo_ov_cr, cio_ov_cr, coj_ov_cr, cij_ov_cr, coo_ov_cr_pow_ncg, coo_ov_cr_pow_ncb, rate_acc = 0.0, rate;
	vector <int> event_details;

	rates.resize(0);
	event_vector.resize(0);
	rate_vector.resize(0);

	//nucleate - 0, grow - 1, shrink - 2, new mixed - 3 dissolve - 4
	//good - 0, bad - 1, mixed - 2
	
	//First, deal with nucleation events

	coo_ov_cr = double(clusters[0][0])/double(number_monomers)*ct_ov_cr;

	if(clusters[0][0] >= 2) {

		rate = f0cr*coo_ov_cr*clusters[0][0];
		if(rate > 0.0) {	
			rates.push_back(rate); //NUCLEATE GOOD CLUSTER
			event_details.push_back(0); event_details.push_back(0);
			event_vector.push_back(event_details);
			event_details.resize(0);
		}
	}

	if(clusters[0][0] >= 2 ) {

		
		rate = 0; //TURN OFF BAD NUCLEATION FOR NOW!
		if(rate > 0.0) {
			rates.push_back(rate); //NUCLEATE BAD CLUSTER
			event_details.push_back(0); event_details.push_back(1);
			event_vector.push_back(event_details);
			event_details.resize(0);
		}
	}

	//Now deal with growth and shrinking events
	
	for(int i = 2; i< number_monomers+1; i++) {

		if(clusters[i][0] > 0) {
			
			rate = f0cr*coo_ov_cr*clusters[i][0];
			if(rate > 0.0) {
				rates.push_back(rate); //GROW GOOD CLUSTER
				event_details.push_back(1); event_details.push_back(0); event_details.push_back(i);
				event_vector.push_back(event_details);
				event_details.resize(0);
			}
			
			if(i<=ncg) rate = f0cr*exp(-eps/nu)*clusters[i][0];
			else rate = f0cr*exp(-eps)*clusters[i][0];
			if(rate > 0.0) {
				rates.push_back(rate); //SHRINK/DISSOLVE GOOD CLUSTER
				event_details.push_back(2); event_details.push_back(0); event_details.push_back(i);
				if(i == 2) event_details[0] = 4;
				event_vector.push_back(event_details);
				event_details.resize(0);
			}
		
			if(i>=ncg) {
				rate = m*f0cr*coo_ov_cr*clusters[i][0];
				if(rate > 0.0) {
					rates.push_back(rate); //GOOD GROWS TO MIXED CLUSTER
					event_details.push_back(3); event_details.push_back(0); event_details.push_back(i);
					event_vector.push_back(event_details);
					event_details.resize(0);
				}
			}
		}

		if(clusters[0][i] > 0) {

			rate = m*f0cr*coo_ov_cr*clusters[0][i];
			if(rate > 0.0) {
				rates.push_back(rate); //GROW BAD CLUSTER
				event_details.push_back(1); event_details.push_back(1); event_details.push_back(i);
				event_vector.push_back(event_details);
				event_details.resize(0);
			}
			
			if(i<=ncb) rate = f0cr*exp(-eps/mu/nu)*clusters[0][i];
			else rate = f0cr*exp(-eps/mu)*clusters[0][i];
			if(rate > 0.0) {
				rates.push_back(rate); //SHRINK/DISSOLVE BAD CLUSTER
				event_details.push_back(2); event_details.push_back(1); event_details.push_back(i);
				if(i == 2) event_details[0] = 4;
				event_vector.push_back(event_details);
				event_details.resize(0);
			}
		}

		for(int j = 1; j < number_monomers+1; j++) {

			if(clusters[i][j]>0) {

				rate = m*f0cr*coo_ov_cr*clusters[i][j];
				if(rate > 0.0){
					rates.push_back(rate); //GROW MIXED CLUSTER
					event_details.push_back(1); event_details.push_back(2); event_details.push_back(i); event_details.push_back(j);
					event_vector.push_back(event_details);
					event_details.resize(0);
				}

				rate = f0cr*exp(-eps/mu)*clusters[i][j];
				if(rate > 0.0) {
					rates.push_back(rate); //SHRINK/DISSOLVE MIXED CLUSTER
					event_details.push_back(2); event_details.push_back(2); event_details.push_back(i); event_details.push_back(j);
					event_vector.push_back(event_details);
					event_details.resize(0);
				}
				}
			}
		
		}

	rate_vector.push_back(0);

	for(int  i = 0; i < rates.size(); i++) {
		rate_acc = rate_acc + rates[i];
		rate_vector.push_back(rate_acc);
	}

}

void bound_misbound_XI::perform_event(vector <int> event_details_in) {

	if(event_details_in[0] == 0) nucleate(event_details_in[1]);

	if(event_details_in[0] == 1) {
		
		if(event_details_in.size()  == 3) grow(event_details_in[1], event_details_in[2]);
		else if (event_details_in.size() == 4) grow(event_details_in[1], event_details_in[2], event_details_in[3]);

	}

	if(event_details_in[0] == 2) {
		
		if(event_details_in.size()  == 3) shrink(event_details_in[1], event_details_in[2]);
		else if (event_details_in.size() == 4) shrink(event_details_in[1], event_details_in[2], event_details_in[3]);

	}

	if(event_details_in[0] == 3) new_mixed(event_details_in[2]);

	if(event_details_in[0] == 4) dissolve_cluster(event_details_in[1]);
}

void bound_misbound_XI::nucleate(int cluster_type){

	if(cluster_type == 0){
		clusters[2][0] = clusters[2][0] + 1;
		clusters[0][0] = clusters[0][0] - 2;

		events_array[0] = events_array[0] + 1;
	}
	else
	{
		clusters[0][2] = clusters[0][2] + 1;
		clusters[0][0] = clusters[0][0] - 2;
	}

	nucleate_ev[cluster_type] += 1;


}

void bound_misbound_XI::grow(int cluster_type, int cluster_size) {

	int super_or_sub;

	clusters[0][0] = clusters[0][0] - 1;

	if(cluster_type == 0) {

		if(HoG_flag == 1 && (cluster_size + 1) == HoG_limit) {
			clusters[cluster_size][0] = clusters[cluster_size][0] - 1;
			clusters[0][0] = clusters[0][0] + cluster_size + 1;

		}
		else {
			clusters[cluster_size][0] = clusters[cluster_size][0] - 1;
			clusters[cluster_size + 1][0] = clusters[cluster_size + 1][0] + 1;
		}

		if(cluster_size<=ncg) events_array[cluster_size] = events_array[cluster_size] + 1;
	}
	else {	
		clusters[0][cluster_size] = clusters[0][cluster_size] - 1;
		clusters[0][cluster_size + 1] = clusters[0][cluster_size + 1] + 1;
	}
	
	if(cluster_size < ncg) super_or_sub = 0;
	else super_or_sub = 1;

	grow_ev[cluster_type][super_or_sub] += 1;

}

void bound_misbound_XI::grow(int cluster_type, int cluster_size_1, int cluster_size_2) {

	clusters[0][0] = clusters[0][0] - 1;

	clusters[cluster_size_1][cluster_size_2] = clusters[cluster_size_1][cluster_size_2] - 1;
	clusters[cluster_size_1][cluster_size_2 + 1] = clusters[cluster_size_1][cluster_size_2 + 1] + 1;

	grow_ev[cluster_type][1] += 1;

}

void bound_misbound_XI::shrink(int cluster_type, int cluster_size) {

	int super_or_sub;

	clusters[0][0] = clusters[0][0] + 1;

	if(cluster_type == 0) {
		clusters[cluster_size][0] = clusters[cluster_size][0] - 1;
		clusters[cluster_size - 1][0] = clusters[cluster_size - 1][0] + 1;
	}
	else {	
		clusters[0][cluster_size] = clusters[0][cluster_size] - 1;
		clusters[0][cluster_size - 1] = clusters[0][cluster_size - 1] + 1;
	}
	
	if(cluster_size < ncg) super_or_sub = 0;
	else super_or_sub = 1;

	shrink_ev[cluster_type][super_or_sub] += 1;

}

void bound_misbound_XI::shrink(int cluster_type, int cluster_size_1, int cluster_size_2) {

	clusters[0][0] = clusters[0][0] + 1;

	clusters[cluster_size_1][cluster_size_2] = clusters[cluster_size_1][cluster_size_2] - 1;
	clusters[cluster_size_1][cluster_size_2 -  1] = clusters[cluster_size_1][cluster_size_2 - 1] + 1;

	shrink_ev[cluster_type][1] += 1;

}

void bound_misbound_XI::new_mixed(int cluster_size_1) { 

	clusters[0][0] = clusters[0][0] - 1;

	clusters[cluster_size_1][0] = clusters[cluster_size_1][0] - 1;
	clusters[cluster_size_1][1] = clusters[cluster_size_1][1] + 1;

	grow_ev[2][1] += 1;

}

void bound_misbound_XI::dissolve_cluster(int cluster_type) {

	if(cluster_type == 0) {
		clusters[2][0] = clusters[2][0] - 1;
		clusters[0][0] = clusters[0][0] + 2;
	}
	else {
		clusters[0][2] = clusters[0][2] - 1;
		clusters[0][0] = clusters[0][0] + 2;
	}

	dissolve_ev[cluster_type] += 1;
}

void bound_misbound_XI::print_state(double time) {

	int sub_bad_acc = 0, super_bad_acc = 0, sub_good_no = 0, sub_good_acc = 0, super_good_acc = 0, super_good_no = 0, mixed_good_acc = 0, mixed_bad_acc = 0, good_cluster_count = 0, bad_cluster_count = 0, mixed_cluster_count = 0;

	for(int i = 2; i<number_monomers+1; i++) {
		if(clusters[i][0] > 0) {
			if(i<=ncg) {
				sub_good_acc = sub_good_acc + clusters[i][0]*i;
				sub_good_no = sub_good_no + clusters[i][0];
			}
			else{
				super_good_acc = super_good_acc + clusters[i][0]*i;
				super_good_no = super_good_no + clusters[i][0];
			}
			good_cluster_count++;
		}

		if(clusters[0][i] > 0) {
			if(i<=ncb) sub_bad_acc = sub_bad_acc + clusters[0][i]*i;
			else super_bad_acc = super_bad_acc + clusters[0][i]*i;
			bad_cluster_count++;
		}

		for(int j = 1; j<number_monomers+1; j++) {
			if(clusters[i][j] > 0) {
				mixed_good_acc = mixed_good_acc + clusters[i][j]*i;
				mixed_bad_acc = mixed_bad_acc + clusters[i][j]*j;
				mixed_cluster_count++;
			}
		}
	}

	printf("#state %d %E %d %d %d %d %d %d %d %d %d %d\n", rank_internal, time, clusters[0][0], sub_good_acc, super_good_acc, sub_bad_acc, super_bad_acc, mixed_good_acc, mixed_bad_acc, sub_good_no, super_good_no, mixed_cluster_count);
}

void bound_misbound_XI::print_events(double time) {

	printf("#events %d %E %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", rank_internal, time, ncg, ncb, nucleate_ev[0], nucleate_ev[1], dissolve_ev[0], dissolve_ev[1], grow_ev[0][0], grow_ev[0][1], grow_ev[1][1], grow_ev[2][1], shrink_ev[0][0], shrink_ev[0][1], shrink_ev[1][1], shrink_ev[2][1]);
}

void bound_misbound_XI::print_clusters(double time) {

	for(int i = 0; i<number_monomers; i++) {
		for(int j = 0; j<number_monomers; j++) {
			if(clusters[i][j] > 0) {
				printf("#clusters %d %E %d %d %d\n", rank_internal, time, i, j, clusters[i][j]);
			}
		}
	}
}

void bound_misbound_XI::print_ave_cluster_size(double time) {

	int no_clusters = 0, good_cluster_size = 0, bad_cluster_size = 0, mixed_cluster_size = 0;

	for(int i = 0; i<number_monomers; i++) {
		for(int j = 0; j<number_monomers; j++) {
			if(clusters[i][j] > 0){
				no_clusters++;
				if(j == 0) good_cluster_size = good_cluster_size + i*clusters[i][j];
				else if(i == 0) bad_cluster_size = bad_cluster_size + j*clusters[i][j];
				else mixed_cluster_size = mixed_cluster_size + (i+j)*clusters[i][j];
			}
		}
	}

	printf("#ave_clust_size %d %E %d %d %d %d\n", rank_internal, time, good_cluster_size, bad_cluster_size, mixed_cluster_size, no_clusters);

}

void bound_misbound_XI::change_strength(double new_eps) {

	eps = new_eps;

}

void bound_misbound_XI::set_rank(int new_rank) {

	rank_internal = new_rank;

}

double bound_misbound_XI::ave_cluster_size() {

	double no_in_good_clusters = 0.0, no_good_post_clusters = 0.0, ave_cluster;

	for(int i = ncg+1; i < number_monomers; i++) {
			if(clusters[i][0] > 0) {
				no_good_post_clusters = no_good_post_clusters + 1.0;
				no_in_good_clusters = no_in_good_clusters + i*clusters[i][0];
			}
	}
	
	if(no_good_post_clusters > 0.0) ave_cluster = no_in_good_clusters/no_good_post_clusters;
	else ave_cluster = 0.0;

	return ave_cluster;
				
}

void bound_misbound_XI::print_smashed(double time) {

	printf("#no_smashed %d %E ", rank_internal, time);

	for(int i = 0; i<=ncg; i++) {
//		if(clusters[i][0] > 0) printf("%f ", float(events_array[i])/float(clusters[i][0]));
//		else printf("0.0 ");
		printf("%d ", events_array[i]);
	}

	printf("%d\n", clusters[0][0]);

}
