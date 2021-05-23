#include "b_m_XI_run_step_II.h"

int main(int argc, char *argv[]) {

	double final_time, init_strength, final_strength, switch_time;
	int seed, rank, size, master;

	master = 0;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	printf("RANK: %d\n", rank);
	
	stringstream convert_1(argv[1]);
	convert_1 >> final_time;
	stringstream convert_2(argv[2]);
	convert_2 >> init_strength;
	stringstream convert_3(argv[3]);
	convert_3 >> final_strength;
	stringstream convert_4(argv[4]);
	convert_4 >> switch_time;

	JRand_dyn *rand_1 = new JRand_dyn();
	seed = ((rank+1)*(unsigned int)time(NULL))%1000;
	rand_1->srandom(seed);

	bound_misbound_XI *bm_model = new bound_misbound_XI(NUMBER_PARTICLES);
	kinetic_monte_carlo *kmc = new kinetic_monte_carlo(bm_model, rand_1);

	bm_model->set_rank(rank);

	kmc->run_dynamics_step(final_time, 10, 10, init_strength, final_strength, switch_time);

	MPI_Finalize();

	delete rand_1;
	delete kmc;
	delete bm_model;

	return 0;

}
