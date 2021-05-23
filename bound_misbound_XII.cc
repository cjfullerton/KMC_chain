#include "bound_misbound_XII.h"


vector <double> bound_misbound_XII::return_rate_vector() {
	
	return rate_vector;

}

vector <int> bound_misbound_XII::return_event(int event_no) { //events now much more complicated...

	return event_vector[event_no];

}


void bound_misbound_XII::calculate_rate_vector() {

	double coo_ov_cr, rate_acc = 0.0, rate;
	int delta;
	vector <int> event_details;

	rates.resize(0);
	event_vector.resize(0);
	rate_vector.resize(0);

	coo_ov_cr = double(clusters[0][0])/double(number_monomers)*ct_ov_cr;

	for(int i = 0; i < occupy_list.size(); i++) {

		if((clusters[0][0] > 0 && occupy_list[i][0]+occupy_list[i][1] < number_monomers) || (occupy_list[i][0] == 0 && occupy_list[i][1] == 0 && clusters[0][0] > 1)){ //check if can grow

			rate = f0cr*coo_ov_cr*clusters[occupy_list[i][0]][occupy_list[i][1]];

			event_details.push_back(-1);

			if(occupy_list[i][1] == 0) {

				if(occupy_list[i][0] == 0) delta = 2;
				else delta = 1;

				event_details.push_back(occupy_list[i][0]+delta); event_details.push_back(occupy_list[i][1]);
				event_details.push_back(occupy_list[i][0]); event_details.push_back(occupy_list[i][1]);

				rates.push_back(rate);
				event_vector.push_back(event_details);

				event_details.pop_back(); event_details.pop_back(); event_details.pop_back(); event_details.pop_back();

				if(occupy_list[i][0] > ncg) {

					delta = 1;

					event_details.push_back(occupy_list[i][0]); event_details.push_back(occupy_list[i][1]+delta);
					event_details.push_back(occupy_list[i][0]); event_details.push_back(occupy_list[i][1]);

					rates.push_back(m*rate);
					event_vector.push_back(event_details);

					event_details.pop_back(); event_details.pop_back(); event_details.pop_back(); event_details.pop_back();
				}
			}
			else if (occupy_list[i][1] > 0) {

				delta = 1;

				event_details.push_back(occupy_list[i][0]); event_details.push_back(occupy_list[i][1]+delta);
				event_details.push_back(occupy_list[i][0]); event_details.push_back(occupy_list[i][1]);

				rates.push_back(m*rate);
				event_vector.push_back(event_details);

				event_details.pop_back(); event_details.pop_back(); event_details.pop_back(); event_details.pop_back();
			}
				
		event_details.resize(0);


		}

		if(occupy_list[i][0] > 1){ //check if can shrink
			
			event_details.push_back(1);

			if(occupy_list[i][1] == 0) {

				if(occupy_list[i][0] == 2) delta = 2;
				else delta = 1;

				event_details.push_back(occupy_list[i][0]-delta); event_details.push_back(occupy_list[i][1]);
				event_details.push_back(occupy_list[i][0]); event_details.push_back(occupy_list[i][1]);

				if(occupy_list[i][0] <= ncg) rate = f0cr*exp(-eps/nu)*clusters[occupy_list[i][0]][occupy_list[i][1]];
				else rate = f0cr*exp(-eps)*clusters[occupy_list[i][0]][occupy_list[i][1]];

				rates.push_back(rate);
				event_vector.push_back(event_details);

				event_details.pop_back(); event_details.pop_back(); event_details.pop_back(); event_details.pop_back();
			}

			if(occupy_list[i][1] > 0) {

				delta = 1;

				event_details.push_back(occupy_list[i][0]); event_details.push_back(occupy_list[i][1]-delta);
				event_details.push_back(occupy_list[i][0]); event_details.push_back(occupy_list[i][1]);

				rate = f0cr*exp(-eps/mu)*clusters[occupy_list[i][0]][occupy_list[i][1]];

				rates.push_back(rate);
				event_vector.push_back(event_details);

				event_details.pop_back(); event_details.pop_back(); event_details.pop_back(); event_details.pop_back();
			}

			event_details.resize(0);

		}
	}

	rate_vector.push_back(0.0);

	for(int  i = 0; i < rates.size(); i++) {
		rate_acc = rate_acc + rates[i];
		rate_vector.push_back(rate_acc);
	}

}

void bound_misbound_XII::perform_event(vector <int> event_details_in) {

	int end_of_list;

	add_cluster.resize(2,0);

	clusters[0][0] = clusters[0][0] + event_details_in[0];

	clusters[event_details_in[1]][event_details_in[2]] = clusters[event_details_in[1]][event_details_in[2]] + 1;
	clusters[event_details_in[3]][event_details_in[4]] = clusters[event_details_in[3]][event_details_in[4]] - 1;

	if(occupy_list_location[event_details_in[1]][event_details_in[2]] == -1) {
		add_cluster[0] = event_details_in[1];
		add_cluster[1] = event_details_in[2];
		occupy_list.push_back(add_cluster);
		occupy_list_location[event_details_in[1]][event_details_in[2]] = occupy_list.size() - 1;
	}

	if(clusters[event_details_in[3]][event_details_in[4]] == 0) {

		end_of_list = occupy_list.size()-1;

		occupy_list_location[occupy_list[end_of_list][0]][occupy_list[end_of_list][1]] = occupy_list_location[event_details_in[3]][event_details_in[4]];
		occupy_list[occupy_list_location[event_details_in[3]][event_details_in[4]]] = occupy_list[end_of_list];
		occupy_list_location[event_details_in[3]][event_details_in[4]] = -1;
		occupy_list.pop_back();

	}
		
		

}

void bound_misbound_XII::print_state(double time) {

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

void bound_misbound_XII::print_events(double time) {

	printf("#events %d %E %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", rank_internal, time, ncg, ncb, nucleate_ev[0], nucleate_ev[1], dissolve_ev[0], dissolve_ev[1], grow_ev[0][0], grow_ev[0][1], grow_ev[1][1], grow_ev[2][1], shrink_ev[0][0], shrink_ev[0][1], shrink_ev[1][1], shrink_ev[2][1]);
}

void bound_misbound_XII::print_clusters(double time) {

	for(int i = 0; i<number_monomers; i++) {
		for(int j = 0; j<number_monomers; j++) {
			if(clusters[i][j] > 0) {
				printf("#clusters %d %E %d %d %d\n", rank_internal, time, i, j, clusters[i][j]);
			}
		}
	}
}

void bound_misbound_XII::print_ave_cluster_size(double time) {

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

void bound_misbound_XII::change_strength(double new_eps) {

	eps = new_eps;

}

void bound_misbound_XII::set_rank(int new_rank) {

	rank_internal = new_rank;

}

double bound_misbound_XII::ave_cluster_size() {

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

void bound_misbound_XII::print_smashed(double time) {

	printf("#no_smashed %d %E ", rank_internal, time);

	for(int i = 0; i<=ncg; i++) {
//		if(clusters[i][0] > 0) printf("%f ", float(events_array[i])/float(clusters[i][0]));
//		else printf("0.0 ");
		printf("%d ", events_array[i]);
	}

	printf("%d\n", clusters[0][0]);

}
