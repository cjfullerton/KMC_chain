#ifndef KINETIC_MONTE_CARLO_H
#define KINETIC_MONTE_CARLO_H

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
#include "JRand_dyn.h"
#include "model_kmc.h"

using namespace std;

class kinetic_monte_carlo {

	public:
		kinetic_monte_carlo(model_kmc *i_model_i, JRand_dyn *i_rand_1) {
			model_i = i_model_i; rand_1 = i_rand_1;
			
			current_time = 0;
		}
		~kinetic_monte_carlo() {
			;
		}

		int binary_search_recur(vector <double> A, double key, int imin, int imax);
		int binary_search_iter(vector <double> A, double key, int imin, int imax);
		int find_midpoint(int imin, int imax);
		void run_dynamics(double total_time);
		void run_dynamics_II(double cluster_size);
		void run_dynamics_III(double total_time);
		void run_dynamics_IV(double total_time, int div_time, int base_time);
		void run_dynamics_step(double total_time, int div_time, int base_time, double init_strength, double final_time, double switch_time);

	private:
		model_kmc *model_i;
		JRand_dyn *rand_1;
		double current_time;
};

#endif
