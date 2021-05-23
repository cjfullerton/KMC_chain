#ifndef MODKMC_H
#define MODKMC_H

#include <vector>

using namespace std;

//Abstract base class for model to run kinetic monte carlo on.
//Needs to provide a vector of rates to pass to monte carlo class,
//ability to advance time (i.e. make a chosen move) and contain
//the state of the system.

class model_kmc {

	public:
		virtual vector <double> return_rate_vector() = 0;
		virtual vector <int> return_event(int event_no) = 0;
		virtual void calculate_rate_vector() = 0;
		virtual void perform_event(vector <int> event_num) = 0;
		virtual void print_state(double time) = 0;
		virtual void print_events(double time) = 0;
		virtual void print_clusters(double time) = 0;
		virtual void print_smashed(double time) = 0;
		virtual void print_ave_cluster_size(double time) = 0;
		virtual void change_strength(double factor) = 0;
		virtual double ave_cluster_size() = 0;
};

#endif
