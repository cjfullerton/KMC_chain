#include "kinetic_monte_carlo.h"

#define KEY_NOT_FOUND -1


int kinetic_monte_carlo::binary_search_recur(vector <double> A, double key, int imin, int imax) {
	

	if(imax < imin)
		//set empty, return value showing not found
		return KEY_NOT_FOUND;
	else
	{
		int imid = find_midpoint(imin, imax);

		if(A[imid] > key)
			return binary_search_recur(A, key, imin, imid-1);
		else if(A[imid] < key)
			return binary_search_recur(A, key, imid+1, imax);
		else
			return imid;
	}
}

int kinetic_monte_carlo::binary_search_iter(vector <double> A, double key, int imin, int imax) {

	while(imax >= imin)
	{
		int imid = find_midpoint(imin, imax);

		if(A[imid] < key && key <= A[imid+1])
			return imid;
		else if(A[imid] <= key)
			imin = imid+1;
		else
			imax = imid -1;
	}

	printf("imax: %d imin: %d \n", imax, imin);

	return KEY_NOT_FOUND;

}

int kinetic_monte_carlo::find_midpoint(int imin, int imax) {
	
	return floor((imin + imax)/2);

}

void kinetic_monte_carlo::run_dynamics(double total_time) {

	int event, event_index, num_poss_events, num_states, diag = 1;
	double rand_to_chose_event, rand_to_advance_time;
	vector <double> rate_vec;

	model_i->print_state(current_time);
	model_i->print_events(current_time);
	model_i->calculate_rate_vector();
	//model_i->print_ave_cluster_size(current_time);

	while(current_time < total_time) {

		rate_vec = model_i->return_rate_vector();
		num_poss_events = rate_vec.size();

		rand_to_chose_event = rand_1->prob_gen()*rate_vec[num_poss_events-1];
		event_index = binary_search_iter(rate_vec, rand_to_chose_event, 0, num_poss_events-1);
		model_i->perform_event(model_i->return_event(event_index));

		rand_to_advance_time = rand_1->prob_gen();
		current_time = current_time + 1.0/rate_vec[num_poss_events-1]*log(1.0/rand_to_advance_time);
		

		model_i->calculate_rate_vector();

		model_i->print_state(current_time);
		model_i->print_events(current_time);
	//	model_i->print_ave_cluster_size(current_time);

	}

	//model_i->print_clusters(current_time);

}

void kinetic_monte_carlo::run_dynamics_II(double cluster_size) {

	int event, event_index, num_poss_events, num_states, diag = 1;
	double rand_to_chose_event, rand_to_advance_time;
	vector <double> rate_vec;

	model_i->calculate_rate_vector();

	while(model_i->ave_cluster_size() < cluster_size) {

		rate_vec = model_i->return_rate_vector();
		num_poss_events = rate_vec.size();

		rand_to_chose_event = rand_1->prob_gen()*rate_vec[num_poss_events-1];
		event_index = binary_search_iter(rate_vec, rand_to_chose_event, 0, num_poss_events-1);
		model_i->perform_event(model_i->return_event(event_index));

		rand_to_advance_time = rand_1->prob_gen();
		current_time = current_time + 1.0/rate_vec[num_poss_events-1]*log(1.0/rand_to_advance_time);

		model_i->calculate_rate_vector();

	}

	printf("%f\n", current_time);

}

void kinetic_monte_carlo::run_dynamics_III(double total_time) {

	int event, event_index, num_poss_events, num_states, diag = 1;
	double rand_to_chose_event, rand_to_advance_time;
	vector <double> rate_vec;

	model_i->print_smashed(current_time);
	model_i->print_state(current_time);

	model_i->calculate_rate_vector();

	while(current_time < total_time) {

		rate_vec = model_i->return_rate_vector();
		num_poss_events = rate_vec.size();

		rand_to_chose_event = rand_1->prob_gen()*rate_vec[num_poss_events-1];
		event_index = binary_search_iter(rate_vec, rand_to_chose_event, 0, num_poss_events-1);
		model_i->perform_event(model_i->return_event(event_index));

		rand_to_advance_time = rand_1->prob_gen();
		current_time = current_time + 1.0/rate_vec[num_poss_events-1]*log(1.0/rand_to_advance_time);
		

		model_i->calculate_rate_vector();

		model_i->print_smashed(current_time);
		model_i->print_state(current_time);

	}

}

void kinetic_monte_carlo::run_dynamics_IV(double total_time, int time_base, int time_div) {

	int event, event_index, num_poss_events, num_states, out_counter, diag = 1;
	double rand_to_chose_event, rand_to_advance_time, time_pow;
	vector <double> rate_vec;
	vector <double> print_times;

	model_i->print_state(current_time);
	model_i->print_events(current_time);
	model_i->calculate_rate_vector();
	
	time_pow = 0;
	out_counter = 0;
	print_times.push_back(0.0);

	while(print_times[out_counter] < total_time) {
		//printf("%f\t%f\n", total_time, print_times[out_counter]);
		for(int i = 0; i<time_div; i++) {
			print_times.push_back(0.001*pow(time_base, time_pow));
			time_pow = time_pow + 1.0/time_div;
			out_counter++;
		}
	}

	out_counter = 0;

	while(current_time < total_time) {

		rate_vec = model_i->return_rate_vector();
		num_poss_events = rate_vec.size();

		rand_to_chose_event = rand_1->prob_gen()*rate_vec[num_poss_events-1];
		while(rand_to_chose_event == 0) rand_to_chose_event = rand_1->prob_gen()*rate_vec[num_poss_events-1]; 
		event_index = binary_search_iter(rate_vec, rand_to_chose_event, 0, num_poss_events-1);
		if(event_index < 0 || event_index > num_poss_events) printf("%d\t%f\t%d\t%d\n", event_index, rand_to_chose_event, 0, num_poss_events-1);
		model_i->perform_event(model_i->return_event(event_index));

		rand_to_advance_time = rand_1->prob_gen();
		while(rand_to_advance_time == 0.0) rand_to_advance_time = rand_1->prob_gen();
		current_time = current_time + 1.0/rate_vec[num_poss_events-1]*log(1.0/rand_to_advance_time);
		if(current_time > total_time) current_time = total_time;
		

		model_i->calculate_rate_vector();
		
		while(current_time >= print_times[out_counter] && out_counter < print_times.size()) {
			printf("%d\t%f\n", out_counter, current_time);
			model_i->print_state(print_times[out_counter]);
			model_i->print_events(print_times[out_counter]);
			out_counter = out_counter + 1;
		}
		

	}

}

void kinetic_monte_carlo::run_dynamics_step(double total_time, int time_base, int time_div, double init_strength, double final_strength, double switch_time) {

	int event, event_index, num_poss_events, num_states, out_counter, diag = 1;
	double rand_to_chose_event, rand_to_advance_time, time_pow;
	vector <double> rate_vec;
	vector <double> print_times;

	model_i->print_state(current_time);
	model_i->print_events(current_time);
	model_i->calculate_rate_vector();
	
	time_pow = 0;
	out_counter = 0;
	print_times.push_back(0.0);

	while(print_times[out_counter] < total_time) {
		//printf("%f\n", print_times[out_counter]);
		for(int i = 0; i<time_div; i++) {
			print_times.push_back(0.001*pow(time_base, time_pow));
			time_pow = time_pow + 1.0/time_div;
			out_counter++;
		}
	}

	model_i->change_strength(init_strength);

	out_counter = 0;

	while(current_time < total_time) {

		rate_vec = model_i->return_rate_vector();
		num_poss_events = rate_vec.size();

		rand_to_chose_event = rand_1->prob_gen()*rate_vec[num_poss_events-1];
		event_index = binary_search_iter(rate_vec, rand_to_chose_event, 0, num_poss_events-1);
		model_i->perform_event(model_i->return_event(event_index));

		rand_to_advance_time = rand_1->prob_gen();
		current_time = current_time + 1.0/rate_vec[num_poss_events-1]*log(1.0/rand_to_advance_time);

		if(current_time >= switch_time) model_i->change_strength(final_strength);

		model_i->calculate_rate_vector();
		
		while(current_time > print_times[out_counter] && out_counter<print_times.size()) {
			model_i->print_state(print_times[out_counter]);
			model_i->print_events(print_times[out_counter]);
			out_counter = out_counter + 1;
		}

	}

}
