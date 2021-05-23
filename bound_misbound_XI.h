#ifndef BOUND_MISBOUND_X
#define BOUND_MISBOUND_X

#include "model_kmc.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#define F0CR 100.0 //Defines `speed' of all reactions
#define EPS 9.0 //Bond strength
#define MU 4.0 //Good strength over bad strength, more than one
#define NU 5.0 //bonds in sub-critical filament are eps/nu or (eps/mu)/nu
#define DEGEN 4.0 //ways of misbinding vs binding
#define CTOVCR 0.01 //total concentration over reference concentration - small
#define HOGFLAG 0 //turns on 'hammer of god'


using namespace std;

class bound_misbound_XI : public model_kmc {

	protected:	

		int number_monomers;
		
		double f0cr;
		double eps;
		double mu;
		double nu;
		double ct_ov_cr;
		double m;

		int ncg;
		int ncb;

		int HoG_flag; //flag to turn on 'hammer of god'
		int HoG_limit; //'hammer of god' smashes all clusters size greater than this
		vector <int> events_array;
		
		int rank_internal;

		vector < vector <int> > clusters;
		vector <double> rate_vector;
		vector <double> rates;
		vector <int> state_vector;
		vector < vector <int> > event_vector;

		vector <int> nucleate_ev;
		vector <int> dissolve_ev;
		vector < vector <int> > grow_ev;
		vector < vector <int> > shrink_ev;

	public:
		bound_misbound_XI(int i_number_monomers) {

			number_monomers = i_number_monomers;

			f0cr = F0CR;
			eps = EPS;
			mu = MU;
			nu = NU;
			m = DEGEN;
			ct_ov_cr = CTOVCR;
			
			ncg = 6;
			ncb = 6;	

			HoG_flag = HOGFLAG;
			HoG_limit = ncg+1;
			events_array.resize(ncg+1,0);

			rank_internal = 0; //set rank to zero, change through set_rank function before outputing any results...

			create_vectors();

			clusters[0][0] = number_monomers;
		}
		~bound_misbound_XI() {;}

		void create_vectors() {
			clusters.resize(number_monomers+1);
			for(int i = 0; i < number_monomers+1; i++) clusters[i].resize(number_monomers+1,0);
			nucleate_ev.resize(2);
			dissolve_ev.resize(2);
			grow_ev.resize(3);
			for(int i = 0; i<3; i++) grow_ev[i].resize(2);
			shrink_ev.resize(3);
			for(int i = 0; i<3; i++) shrink_ev[i].resize(2);
		}

		vector <double> return_rate_vector();
		vector <int> return_state_vector();
		vector <int> return_event(int event_no);
		void calculate_rate_vector();
		void perform_event(vector <int> event_details_in);
		void nucleate(int cluster_type);
		void grow(int cluster_type, int cluster_size_1);
		void grow(int cluster_type, int cluster_size_1, int cluster_size_2);
		void shrink(int cluster_type, int cluster_size_1);
		void shrink(int cluster_type, int cluster_size_1, int cluster_size_2);
		void new_mixed(int cluster_size_1);
		void dissolve_cluster(int cluster_type);
		void print_state(double time);
		void print_events(double time);
		void print_clusters(double time);
		void print_ave_cluster_size(double time);
		void print_smashed(double time);
		void change_strength(double new_eps);
		void set_rank(int new_rank);
		double ave_cluster_size();

};

#endif
