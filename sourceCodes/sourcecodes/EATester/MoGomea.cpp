/*
These source codes were created with the use of source codes published for the papers listed below.

 * 1. N.H. Luong, H. La Poutré, and P.A.N. Bosman: Multi-objective Gene-pool
 * Optimal Mixing Evolutionary Algorithms with the Interleaved Multi-start Scheme.
 * In Swarm and Evolutionary Computation, vol. 40, June 2018, pages 238-254,
 * Elsevier, 2018.
 *
 * 2. N.H. Luong, H. La Poutré, and P.A.N. Bosman: Multi-objective Gene-pool
 * Optimal Mixing Evolutionary Algorithms. In Dirk V. Arnold, editor,
 * Proceedings of the Genetic and Evolutionary Computation Conference GECCO 2014:
 * pages 357-364, ACM Press New York, New York, 2014.


 2019.03.21 Michal Przewozniczek: The sources were made objective and the test case objects implemented in this framework were introduced into the original source code

*/


#include "MoGomea.h"
#include "PopulationOptimizer.h"
#include "UIntCommandParam.h"
#include "FloatCommandParam.h"




#define FALSE 0
#define TRUE 1

#define NOT_EXTREME_CLUSTER -1

#define MINIMIZATION 1
#define MAXIMIZATION 2

#define ZEROMAX_ONEMAX 0
#define TRAP5 1
#define KNAPSACK 2
#define LOTZ 3
#define MAXCUT 4




using  namespace MoGomea;


uint32_t CMoGomea::iERROR_PARENT_CMoGomeaOptimizer = CError::iADD_ERROR_PARENT("iERROR_PARENT_CMoGomeaOptimizer");
uint32_t CMoGomea::iERROR_CODE_MOGOMEA_GENOTYPE_LEN_BELOW_0 = CError::iADD_ERROR("iERROR_CODE_MOGOMEA_GENOTYPE_LEN_BELOW_0");


/**
 * Computes the two-log of x.
 */
double math_log_two_ = log(2.0);
double log2_(double x)
{
	return(log(x) / math_log_two_);
}

//---------------------------------------------CMoGomea-------------------------------------------------------
CMoGomea::CMoGomea(CProblem<CBinaryCoding, CBinaryCoding>* pcProblem, CLog * pcLog, uint32_t iRandomSeed)
	: CBinaryMultiObjectiveOptimizer(pcProblem, pcLog, iRandomSeed)
	//: CBinaryOptimizer(pcProblem, pcLog, iRandomSeed)
{
	//pc_problem = (CBinaryMultiObjectiveProblem *) pcProblem;

	pc_genotype = NULL;

};//CMoGomea::CMoGomea(CProblem<CBinaryCoding, CBinaryCoding>* pcProblem, CLog * pcLog, uint32_t iRandomSeed)


CMoGomea::CMoGomea(CMoGomea *pcOther) : CBinaryMultiObjectiveOptimizer(pcOther)//CBinaryOptimizer(pcOther)
{
	::MessageBox(NULL, "No implementation: CMoGomea::CMoGomea(CMoGomea *pcOther) : CBinaryMultiObjectiveOptimizer(pcOther)//CBinaryOptimizer(pcOther)", "Implementation missing", MB_OK);
};//CMoGomea::CMoGomea(CMoGomea *pcOther) : CBinaryMultiObjectiveOptimizer(pcOther)//CBinaryOptimizer(pcOther)


CMoGomea::~CMoGomea()
{
	if (pc_genotype != NULL)  delete  pc_genotype;

	logNumberOfEvaluationsAtVTR();
	writeCurrentElitistArchive(TRUE);
	ezilaitiniArrayOfPopulation();
	ezilaitiniArrayOfParetoFronts();

	//::Tools::vShow("deinit ok");
}//CMoGomea::~CMoGomea()


void CMoGomea::logNumberOfEvaluationsAtVTR()
{
	//in this project the optimal front is supported by CBinaryMultiObjectiveProblem
	return;

	/*int      default_front_size;
	double **default_front, metric_elitist_archive;
	FILE *file;
	char string[1000];

	if (use_vtr == FALSE)
		return;

	if (haveDPFSMetric())
	{
		default_front = getDefaultFront(&default_front_size);
		metric_elitist_archive = computeDPFSMetric(default_front, default_front_size, elitist_archive_objective_values, elitist_archive_size);

		sprintf(string, "number_of_evaluations_when_all_points_found_%d.dat", number_of_parameters);
		file = fopen(string, "a");
		if (metric_elitist_archive <= vtr)
		{
			fprintf(file, "%d\n", number_of_evaluations);
		}
		else
		{
			fprintf(file, "Cannot find all points within current budget!\n");
		}
		fclose(file);
	}*/
}//void CMoGomea::logNumberOfEvaluationsAtVTR()



void CMoGomea::writeCurrentElitistArchive(char final)
{
	return;
	int   i, j, k, index;
	char  string[1000];
	FILE *file;

	/* Elitist archive */
	if (final)
		sprintf(string, "elitist_archive_generation_final.dat");
	else
		sprintf(string, "elitist_archive_at_evaluation_%d.dat", number_of_evaluations);
	file = fopen(string, "w");

	for (i = 0; i < elitist_archive_size; i++)
	{
		for (j = 0; j < number_of_objectives; j++)
		{
			sprintf(string, "%13e ", elitist_archive_objective_values[i][j]);
			fputs(string, file);
		}

		sprintf(string, "   %f     ", elitist_archive_constraint_values[i]);
		fputs(string, file);

		for (j = 0; j < number_of_parameters; j++)
		{
			sprintf(string, "%d ", elitist_archive[i][j]);
			fputs(string, file);
		}
		sprintf(string, "\n");
		fputs(string, file);
	}
	fclose(file);
}//void CMoGomea::writeCurrentElitistArchive(char final)


void CMoGomea::ezilaitiniArrayOfPopulation()
{
	int i;
	for (i = 0; i < number_of_populations; i++)
	{
		population_id = i;
		assignPointersToCorrespondingPopulation();
		ezilaitiniMemoryOfCorrespondingPopulation();
	}
	free(array_of_populations);
	free(array_of_objective_values);
	free(array_of_constraint_values);
	free(array_of_population_sizes);
	free(array_of_objective_ranges);
	free(array_of_t_NIS);
	free(array_of_number_of_generations);
	free(array_of_number_of_evaluations_per_population);
	free(array_of_number_of_clusters);
}//void CMoGomea::ezilaitiniArrayOfPopulation()



void CMoGomea::assignPointersToCorrespondingPopulation()
{
	population = array_of_populations[population_id];
	objective_values = array_of_objective_values[population_id];
	constraint_values = array_of_constraint_values[population_id];
	population_size = array_of_population_sizes[population_id];
	objective_ranges = array_of_objective_ranges[population_id];
	t_NIS = array_of_t_NIS[population_id];
	number_of_generations = array_of_number_of_generations[population_id];
	number_of_mixing_components = array_of_number_of_clusters[population_id];
}//void CMoGomea::assignPointersToCorrespondingPopulation()


void CMoGomea::ezilaitiniMemoryOfCorrespondingPopulation()
{
	int i;

	for (i = 0; i < population_size; i++)
	{
		free(population[i]);
		free(objective_values[i]);
	}
	free(population);
	free(objective_values);
	free(constraint_values);
	free(objective_ranges);
}//void CMoGomea::ezilaitiniMemoryOfCorrespondingPopulation()


void CMoGomea::ezilaitiniArrayOfParetoFronts()
{
	return;
	int i, j;

	FILE *file;
	file = fopen("population_status.dat", "w");
	for (i = 0; i < number_of_populations; i++)
	{
		fprintf(file, "Pop %d: %d\n", ((int)(pow(2, i)))*smallest_population_size, array_of_population_statuses[i]);
	}
	fclose(file);
	for (i = 0; i < maximum_number_of_populations; i++)
	{
		if (array_of_Pareto_front_size_of_each_population[i] > 0)
		{
			for (j = 0; j < array_of_Pareto_front_size_of_each_population[i]; j++)
				free(array_of_Pareto_front_of_each_population[i][j]);
			free(array_of_Pareto_front_of_each_population[i]);
		}
	}
	free(array_of_Pareto_front_of_each_population);
	free(array_of_Pareto_front_size_of_each_population);
	free(array_of_population_statuses);
}//void CMoGomea::ezilaitiniArrayOfParetoFronts()

short CMoGomea::haveDPFSMetric(void)
{
	//in this project the optimal front is supported by CBinaryMultiObjectiveProblem
	return(0);

	
	/*int default_front_size;

	getDefaultFront(&default_front_size);
	if (default_front_size > 0)
		return(1);

	return(0);*/
}//short CMoGomea::haveDPFSMetric(void)


CError CMoGomea::eConfigure(istream *psSettings)
{
	CError c_err(iERROR_PARENT_CMoGomeaOptimizer);
	
	c_err = CBinaryOptimizer::eConfigureMoGomea(psSettings);
	
	return(c_err);
};//CError CMoGomea::eConfigure(istream *psSettings)





void CMoGomea::vInitialize(time_t tStartTime)
{
	CBinaryOptimizer::vInitialize(tStartTime);
	t_start = tStartTime;

	CError  c_err(iERROR_PARENT_CMoGomeaOptimizer);
	i_templ_length = pc_problem->pcGetEvaluation()->iGetNumberOfElements();

	if (i_templ_length <= 0)
	{
		c_err.vSetError(CMoGomea::iERROR_CODE_MOGOMEA_GENOTYPE_LEN_BELOW_0);
		return;
	}//if  (i_templ_length  <=  0)


	pc_log->vPrintLine("Initializing...", true);

	pc_genotype = new CBinaryCoding(i_templ_length);

	number_of_parameters = i_templ_length;

	CBinaryMultiObjectiveProblem  *pc_multi_problem;
	pc_multi_problem = (CBinaryMultiObjectiveProblem *)pc_problem->pcGetEvaluation();
	number_of_objectives = pc_multi_problem->iGetActiveMeasures();

	optimization = new char[number_of_objectives];
	for (int ii = 0; ii < number_of_objectives; ii++)
		optimization[ii] = MAXIMIZATION;


	c_time_counter.vSetStartNow();

	smallest_population_size = 8;

	v_initializeCommonVariables();
	pc_log->vPrintLine("v_initializeCommonVariables();", true);
		
	initializeMemoryForArrayOfPopulations();
	pc_log->vPrintLine("initializeMemoryForArrayOfPopulations();", true);

	initializeArrayOfParetoFronts();
	pc_log->vPrintLine("initializeArrayOfParetoFronts();", true);

	pc_log->vPrintLine("DONE...", true);


	/*CNSGA2Individual  *pc_buf;
	for (int ii = 0; ii < i_pop_size; ii++)
	{
		pc_buf = new CNSGA2Individual(i_templ_length, (CBinaryMultiObjectiveProblem *)pc_problem->pcGetEvaluation(), this);
		pc_buf->vInitialize();
		v_population.push_back(pc_buf);
	}//for (int ii = 0; ii < i_pop_size; ii++)*/

	c_time_counter.vSetStartNow();
};//void CMoGomea::vInitialize(time_t tStartTime)


void CMoGomea::v_initializeCommonVariables()
{
	int i;

	//initializeRandomNumberGenerator();
	generation_base = 2;

	number_of_generations = 0;
	number_of_evaluations = 0;
	objective_discretization_in_effect = 0;
	elitist_archive_size = 0;
	elitist_archive_capacity = 10;
	//elitist_archive = (char **)Malloc(elitist_archive_capacity * sizeof(char *));
	elitist_archive = new char*[elitist_archive_capacity];
	//elitist_archive_objective_values = (double **)Malloc(elitist_archive_capacity * sizeof(double *));
	elitist_archive_objective_values = new  double *[elitist_archive_capacity];
	//elitist_archive_constraint_values = (double *)Malloc(elitist_archive_capacity * sizeof(double));
	elitist_archive_constraint_values = new double[elitist_archive_capacity];

	for (i = 0; i < elitist_archive_capacity; i++)
	{
		//elitist_archive[i] = (char *)Malloc(number_of_parameters * sizeof(char));
		elitist_archive[i] = new char[number_of_parameters];
		//elitist_archive_objective_values[i] = (double *)Malloc(number_of_objectives * sizeof(double));
		elitist_archive_objective_values[i] = new double[number_of_objectives];
	}
	elitist_archive_copy = NULL;
	//objective_discretization = (double *)Malloc(number_of_objectives * sizeof(double));
	objective_discretization = new double[number_of_objectives];

	//MI_matrix = (double **)Malloc(number_of_parameters * sizeof(double *));
	MI_matrix = new double *[number_of_parameters];
	for (i = 0; i < number_of_parameters; i++)
		MI_matrix[i] = new double[number_of_parameters];

	population_indices_of_cluster_members = NULL;
	population_cluster_sizes = NULL;

	offspring = NULL;

	number_of_populations = 0;

	lt = NULL;
}//void CMoGomea::v_initializeCommonVariables()




void CMoGomea::initializeMemoryForArrayOfPopulations()
{
	int i;
	maximum_number_of_populations = 20;

	//array_of_populations = (char***)Malloc(maximum_number_of_populations * sizeof(char**));
	array_of_populations = new char**[maximum_number_of_populations];
	//array_of_objective_values = (double***)Malloc(maximum_number_of_populations * sizeof(double**));
	array_of_objective_values = new double **[maximum_number_of_populations];
	//array_of_constraint_values = (double**)Malloc(maximum_number_of_populations * sizeof(double*));
	array_of_constraint_values = new double *[maximum_number_of_populations];
	//array_of_objective_ranges = (double**)Malloc(maximum_number_of_populations * sizeof(double));   //prw: MISTAKE in the original source codes, it should be sizeof(double *), it works, because sizeof(double) is larger
	array_of_objective_ranges = new double*[maximum_number_of_populations];

	//array_of_t_NIS = (int*)Malloc(maximum_number_of_populations * sizeof(int));
	array_of_t_NIS = new int[maximum_number_of_populations];
	//array_of_number_of_generations = (int*)Malloc(maximum_number_of_populations * sizeof(int));
	array_of_number_of_generations = new int[maximum_number_of_populations];
	for (i = 0; i < maximum_number_of_populations; i++)
	{
		array_of_number_of_generations[i] = 0;
		array_of_t_NIS[i] = 0;
	}

	//array_of_number_of_evaluations_per_population = (long*)Malloc(maximum_number_of_populations * sizeof(long));
	array_of_number_of_evaluations_per_population = new long[maximum_number_of_populations];
	for (i = 0; i < maximum_number_of_populations; i++)
		array_of_number_of_evaluations_per_population[i] = 0;

	/* Popupulation-sizing free scheme. */
	//array_of_population_sizes = (int*)Malloc(maximum_number_of_populations * sizeof(int));
	array_of_population_sizes = new int[maximum_number_of_populations];
	array_of_population_sizes[0] = smallest_population_size;
	for (i = 1; i < maximum_number_of_populations; i++)
		array_of_population_sizes[i] = array_of_population_sizes[i - 1] * 2;

	/* Number-of-clusters parameter-free scheme. */
	//array_of_number_of_clusters = (int*)Malloc(maximum_number_of_populations * sizeof(int));
	array_of_number_of_clusters = new int[maximum_number_of_populations];
	array_of_number_of_clusters[0] = number_of_objectives + 1;
	for (i = 1; i < maximum_number_of_populations; i++)
		array_of_number_of_clusters[i] = array_of_number_of_clusters[i - 1] + 1;
}//void CMoGomea::initializeMemoryForArrayOfPopulations()


void CMoGomea::initializeArrayOfParetoFronts()
{
	int i;

	//array_of_population_statuses = (char*)Malloc(maximum_number_of_populations * sizeof(char));
	array_of_population_statuses = new char[maximum_number_of_populations];
	for (i = 0; i < maximum_number_of_populations; i++)
		array_of_population_statuses[i] = TRUE;

	//array_of_Pareto_front_size_of_each_population = (int*)Malloc(maximum_number_of_populations * sizeof(int));
	array_of_Pareto_front_size_of_each_population = new int[maximum_number_of_populations];
	for (i = 0; i < maximum_number_of_populations; i++)
		array_of_Pareto_front_size_of_each_population[i] = 0;

	//array_of_Pareto_front_of_each_population = (double***)Malloc(maximum_number_of_populations * sizeof(double**));
	array_of_Pareto_front_of_each_population = new double**[maximum_number_of_populations];
	for (i = 0; i < maximum_number_of_populations; i++)
		array_of_Pareto_front_of_each_population[i] = NULL;
}//void CMoGomea::initializeArrayOfParetoFronts()


void CMoGomea::schedule_runMultiplePop_clusterPop_learnPop_improvePop()
{
	int i;
	smallest_population_size = 8;

	initializeMemoryForArrayOfPopulations();
	initializeArrayOfParetoFronts();
	/*while (!checkTerminationCondition())
	{
		population_id = 0;
		do
		{
			if (array_of_number_of_generations[population_id] == 0)
			{
				population_size = array_of_population_sizes[population_id];
				number_of_mixing_components = array_of_number_of_clusters[population_id];

				initialize();

				putInitializedPopulationIntoArray();

				if (stop_population_when_front_is_covered)
				{
					updateParetoFrontForCurrentPopulation(objective_values, constraint_values, population_size);
					checkWhichSmallerPopulationsNeedToStop();
				}

				writeGenerationalStatistics();
			}
			else if (array_of_population_statuses[population_id] == TRUE)
			{
				assignPointersToCorrespondingPopulation();

				learnLinkageOnCurrentPopulation();

				improveCurrentPopulation();

				selectFinalSurvivors();

				computeObjectiveRanges();

				adaptObjectiveDiscretization();

				array_of_t_NIS[population_id] = t_NIS;

				if (stop_population_when_front_is_covered)
				{
					updateParetoFrontForCurrentPopulation(objective_values, constraint_values, population_size);
					checkWhichSmallerPopulationsNeedToStop();
				}

				writeGenerationalStatistics();
			}

			array_of_number_of_generations[population_id]++;
			if (use_print_progress_to_screen)
				printf("%d ", array_of_number_of_generations[population_id]);
			population_id++;
			if (checkTerminationCondition() == TRUE)
				break;
		} while (array_of_number_of_generations[population_id - 1] % generation_base == 0);
		if (use_print_progress_to_screen)
			printf(":   %d\n", number_of_evaluations);
	}

	if (use_print_progress_to_screen)
	{
		printf("Population Status:\n");
		for (i = 0; i < number_of_populations; i++)
			printf("Pop %d: %d\n", ((int)(pow(2, i)))*smallest_population_size, array_of_population_statuses[i]);
	}
	logNumberOfEvaluationsAtVTR();
	writeCurrentElitistArchive(TRUE);
	ezilaitiniArrayOfPopulation();
	ezilaitiniArrayOfParetoFronts();*/
}//void CMoGomea::schedule_runMultiplePop_clusterPop_learnPop_improvePop()



void CMoGomea::schedule()
{
	schedule_runMultiplePop_clusterPop_learnPop_improvePop();
}//void CMoGomea::schedule()



bool CMoGomea::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)
{
	//return(bRunIteration_dummy(iIterationNumber, tStartTime));



	population_id = 0;
	do
	{
		if (array_of_number_of_generations[population_id] == 0)
		{
			population_size = array_of_population_sizes[population_id];
			number_of_mixing_components = array_of_number_of_clusters[population_id];

			initialize();

			putInitializedPopulationIntoArray();

			if (stop_population_when_front_is_covered)
			{
				updateParetoFrontForCurrentPopulation(objective_values, constraint_values, population_size);
				checkWhichSmallerPopulationsNeedToStop();
			}

			writeGenerationalStatistics();
		}
		else if (array_of_population_statuses[population_id] == TRUE)
		{
			assignPointersToCorrespondingPopulation();

			learnLinkageOnCurrentPopulation();

			improveCurrentPopulation();

			selectFinalSurvivors();

			computeObjectiveRanges();

			adaptObjectiveDiscretization();

			array_of_t_NIS[population_id] = t_NIS;

			if (stop_population_when_front_is_covered)
			{
				updateParetoFrontForCurrentPopulation(objective_values, constraint_values, population_size);
				checkWhichSmallerPopulationsNeedToStop();
			}

			writeGenerationalStatistics();
		}

		array_of_number_of_generations[population_id]++;
		if (use_print_progress_to_screen)
			printf("%d ", array_of_number_of_generations[population_id]);
		population_id++;

		//PRW: fix stop condition check here
		//if (checkTerminationCondition() == TRUE)
			//break;
	} while (array_of_number_of_generations[population_id - 1] % generation_base == 0);
	if (use_print_progress_to_screen)
		printf(":   %d\n", number_of_evaluations);

	int  i;
	if (use_print_progress_to_screen)
	{
		printf("Population Status:\n");
		for (i = 0; i < number_of_populations; i++)
			printf("Pop %d: %d\n", ((int)(pow(2, i)))*smallest_population_size, array_of_population_statuses[i]);
	}


	CBinaryMultiObjectiveProblem  *pc_multi_problem;
	pc_multi_problem = (CBinaryMultiObjectiveProblem *)pc_problem->pcGetEvaluation();


	//double  d_fit;
	//d_fit = pc_multi_problem->dEvaluate(&c_genotype);

	double  d_time_passed;
	c_time_counter.bGetTimePassed(&d_time_passed);

	CString  s_buf;
	//s_buf.Format("iteration: %d  pop: %d; GenSize:%d PFsize:%d HyperVolume: %.8lf InvGenDist: %.8lf [time:%.2lf]", iIterationNumber, v_population.size(), v_non_dominated_pfs.at(0).at(0)->pc_genotype->iGetNumberOfBits(), (int)dPFQualityPointNum(), dPFQualityHyperVolume(), dPFQualityInverseGenerationalDistance(), d_time_passed);
	s_buf.Format("iteration: %d  PFsize:%d HyperVolume: %.8lf InvGenDist: %.8lf GenDist:%.8lf MaxSpread:%.8lf DominatedPF: %d [time:%.2lf] [ffe: %.0lf]", iIterationNumber, (int)dPFQualityPointNum(), dPFQualityHyperVolume(), dPFQualityInverseGenerationalDistance(), dPFQualityGenerationalDistance(), dPFQualityMaximumSpread(), (int) dPFQualityDominatedOptimalPoints(), d_time_passed, (double) pc_multi_problem->iGetFFE());
	pc_log->vPrintLine(s_buf, true);

	
	if (pc_multi_problem->pvGetParetoFront()->size() > 0)
	{
		double  d_fit;
		d_fit = pc_multi_problem->dEvaluate(pc_multi_problem->pvGetParetoFront()->at(0)->pcGetGenotype());
		
		b_update_best_individual(iIterationNumber, tStartTime, pc_multi_problem->pvGetParetoFront()->at(0)->pcGetGenotype()->piGetBits(), d_fit);
	}//if (pc_multi_problem->pvGetParetoFront()->size() > 0)


	

	

	return(true);
};//bool CMoGomea::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)




/**
 * Performs initializations that are required before starting a run.
 */
void CMoGomea::initialize()
{
	number_of_populations++;

	initializeMemory();

	initializePopulationAndFitnessValues();

	computeObjectiveRanges();
}//void CMoGomea::initialize()



/**
 * Initializes the memory.
 */
void CMoGomea::initializeMemory(void)
{
	int i;

	//::Tools::vShow(population_size);
	//::Tools::vShow(number_of_parameters);
	//::Tools::vShow(number_of_objectives);

	//objective_ranges = (double *)Malloc(population_size * sizeof(double));
	objective_ranges = new double[population_size];
	//population = (char **)Malloc(population_size * sizeof(char *));
	population = new char*[population_size];
	//objective_values = (double **)Malloc(population_size * sizeof(double *));
	objective_values = new double *[population_size];
	//constraint_values = (double *)Malloc(population_size * sizeof(double));
	constraint_values = new double[population_size];

	for (i = 0; i < population_size; i++)
	{
		//population[i] = (char *)Malloc(number_of_parameters * sizeof(char));
		population[i] = new char[number_of_parameters];
		//objective_values[i] = (double *)Malloc(number_of_objectives * sizeof(double));
		objective_values[i] = new double[number_of_objectives];
	}

	t_NIS = 0;
	number_of_generations = 0;
}//void CMoGomea::initializeMemory(void)


/**
 * Initializes the population and the objective values by randomly
 * generation n solutions.
 */
void CMoGomea::initializePopulationAndFitnessValues()
{
	int i, j;
	for (i = 0; i < population_size; i++)
	{
		for (j = 0; j < number_of_parameters; j++)
			population[i][j] = (randomInt(2) == 1) ? TRUE : FALSE;
		evaluateIndividual(population[i], objective_values[i], &(constraint_values[i]), NOT_EXTREME_CLUSTER);
		updateElitistArchive(population[i], objective_values[i], constraint_values[i]);
	}
}//void CMoGomea::initializePopulationAndFitnessValues()


/**
 * Computes the ranges of all fitness values
 * of all solutions currently in the populations.
 */
void CMoGomea::computeObjectiveRanges(void)
{
	int    i, j;
	double low, high;

	for (j = 0; j < number_of_objectives; j++)
	{
		low = objective_values[0][j];
		high = objective_values[0][j];

		for (i = 0; i < population_size; i++)
		{
			if (objective_values[i][j] < low)
				low = objective_values[i][j];
			if (objective_values[i][j] > high)
				high = objective_values[i][j];
		}

		objective_ranges[j] = high - low;
	}
}//void CMoGomea::computeObjectiveRanges(void)




void CMoGomea::learnLinkageOnCurrentPopulation()
{
	int i, j, k, size_of_one_cluster;

	initializeClusters();

	population_indices_of_cluster_members = clustering(objective_values, population_size, number_of_objectives,
		number_of_mixing_components, &size_of_one_cluster);
	//population_cluster_sizes = (int*)Malloc(number_of_mixing_components * sizeof(int));
	population_cluster_sizes = new int[number_of_mixing_components];
	for (k = 0; k < number_of_mixing_components; k++)
		population_cluster_sizes[k] = size_of_one_cluster;

	// find extreme-region clusters
	determineExtremeClusters();

	// learn linkage tree for every cluster
	for (i = 0; i < number_of_mixing_components; i++)
		learnLinkageTree(i);
}//void CMoGomea::learnLinkageOnCurrentPopulation()




int** CMoGomea::clustering(double **objective_values_pool, int pool_size, int number_of_dimensions,
	int number_of_clusters, int *pool_cluster_size)
{
	int i, j, k, j_min, number_to_select,
		*pool_indices_of_leaders, *k_means_cluster_sizes, **pool_indices_of_cluster_members_k_means,
		**pool_indices_of_cluster_members, size_of_one_cluster;
	double distance, distance_smallest, epsilon,
		**objective_values_pool_scaled, **objective_means_scaled_new, *distances_to_cluster;

	if (number_of_clusters > 1)
		*pool_cluster_size = (2 * pool_size) / number_of_clusters;
	else
	{
		*pool_cluster_size = pool_size;
		//pool_indices_of_cluster_members = (int**)Malloc(number_of_clusters * sizeof(int*));
		pool_indices_of_cluster_members = new int*[number_of_clusters];
		//pool_indices_of_cluster_members[0] = (int*)Malloc(pool_size * sizeof(int));
		pool_indices_of_cluster_members[0] = new int[pool_size];
		for (i = 0; i < pool_size; i++)
			pool_indices_of_cluster_members[0][i] = i;
		return (pool_indices_of_cluster_members);
	}

	size_of_one_cluster = *pool_cluster_size;

	/* Determine the leaders */
	//objective_values_pool_scaled = (double **)Malloc(pool_size * sizeof(double *));
	objective_values_pool_scaled = new double*[pool_size];
	for (i = 0; i < pool_size; i++)
		//objective_values_pool_scaled[i] = (double *)Malloc(number_of_dimensions * sizeof(double));
		objective_values_pool_scaled[i] = new double[number_of_dimensions];
	for (i = 0; i < pool_size; i++)
		for (j = 0; j < number_of_dimensions; j++)
			objective_values_pool_scaled[i][j] = objective_values_pool[i][j] / objective_ranges[j];

	/* Heuristically find k far-apart leaders */
	number_to_select = number_of_clusters;
	pool_indices_of_leaders = greedyScatteredSubsetSelection(objective_values_pool_scaled, pool_size, number_of_dimensions, number_to_select);

	for (i = 0; i < number_of_clusters; i++)
		for (j = 0; j < number_of_dimensions; j++)
			objective_means_scaled[i][j] = objective_values_pool[pool_indices_of_leaders[i]][j] / objective_ranges[j];

	/* Perform k-means clustering with leaders as initial mean guesses */
	//objective_means_scaled_new = (double **)Malloc(number_of_clusters * sizeof(double *));
	objective_means_scaled_new = new double*[number_of_clusters];
	for (i = 0; i < number_of_clusters; i++)
		//objective_means_scaled_new[i] = (double *)Malloc(number_of_dimensions * sizeof(double));
		objective_means_scaled_new[i] = new double[number_of_dimensions];

	//pool_indices_of_cluster_members_k_means = (int **)Malloc(number_of_clusters * sizeof(int *));
	pool_indices_of_cluster_members_k_means = new int*[number_of_clusters];
	for (i = 0; i < number_of_clusters; i++)
		//pool_indices_of_cluster_members_k_means[i] = (int *)Malloc(pool_size * sizeof(int));
		pool_indices_of_cluster_members_k_means[i] = new int[pool_size];

	//k_means_cluster_sizes = (int *)Malloc(number_of_clusters * sizeof(int));
	k_means_cluster_sizes = new int[number_of_clusters];

	epsilon = 1e+308;
	while (epsilon > 1e-10)
	{
		for (j = 0; j < number_of_clusters; j++)
		{
			k_means_cluster_sizes[j] = 0;
			for (k = 0; k < number_of_dimensions; k++)
				objective_means_scaled_new[j][k] = 0.0;
		}

		for (i = 0; i < pool_size; i++)
		{
			j_min = -1;
			distance_smallest = -1;
			for (j = 0; j < number_of_clusters; j++)
			{
				distance = distanceEuclidean(objective_values_pool_scaled[i], objective_means_scaled[j], number_of_dimensions);
				if ((distance_smallest < 0) || (distance < distance_smallest))
				{
					j_min = j;
					distance_smallest = distance;
				}
			}
			pool_indices_of_cluster_members_k_means[j_min][k_means_cluster_sizes[j_min]] = i;
			for (k = 0; k < number_of_dimensions; k++)
				objective_means_scaled_new[j_min][k] += objective_values_pool_scaled[i][k];
			k_means_cluster_sizes[j_min]++;
		}

		for (j = 0; j < number_of_clusters; j++)
			for (k = 0; k < number_of_dimensions; k++)
				objective_means_scaled_new[j][k] /= (double)k_means_cluster_sizes[j];

		epsilon = 0;
		for (j = 0; j < number_of_clusters; j++)
		{
			epsilon += distanceEuclidean(objective_means_scaled[j], objective_means_scaled_new[j], number_of_dimensions);
			for (k = 0; k < number_of_dimensions; k++)
				objective_means_scaled[j][k] = objective_means_scaled_new[j][k];
		}
	}

	/* Shrink or grow the result of k-means clustering to get the final equal-sized clusters */
	//pool_indices_of_cluster_members = (int**)Malloc(number_of_clusters * sizeof(int*));
	pool_indices_of_cluster_members = new int*[number_of_clusters];
	//distances_to_cluster = (double *)Malloc(pool_size * sizeof(double));
	distances_to_cluster = new double[pool_size];
	for (i = 0; i < number_of_clusters; i++)
	{
		for (j = 0; j < pool_size; j++)
			distances_to_cluster[j] = distanceEuclidean(objective_values_pool_scaled[j], objective_means_scaled[i], number_of_dimensions);

		for (j = 0; j < k_means_cluster_sizes[i]; j++)
			distances_to_cluster[pool_indices_of_cluster_members_k_means[i][j]] = 0;

		pool_indices_of_cluster_members[i] = mergeSort(distances_to_cluster, pool_size);
	}

	// Re-calculate clusters' means
	for (i = 0; i < number_of_clusters; i++)
	{
		for (j = 0; j < number_of_dimensions; j++)
			objective_means_scaled[i][j] = 0.0;

		for (j = 0; j < size_of_one_cluster; j++)
		{
			for (k = 0; k < number_of_dimensions; k++)
				objective_means_scaled[i][k] +=
				objective_values_pool_scaled[pool_indices_of_cluster_members[i][j]][k];
		}

		for (j = 0; j < number_of_dimensions; j++)
		{
			objective_means_scaled[i][j] /= (double)size_of_one_cluster;
		}
	}

	free(distances_to_cluster);
	free(k_means_cluster_sizes);
	for (i = 0; i < number_of_clusters; i++)
		free(pool_indices_of_cluster_members_k_means[i]);
	free(pool_indices_of_cluster_members_k_means);
	for (i = 0; i < number_of_clusters; i++)
		free(objective_means_scaled_new[i]);
	free(objective_means_scaled_new);
	for (i = 0; i < pool_size; i++)
		free(objective_values_pool_scaled[i]);
	free(objective_values_pool_scaled);
	free(pool_indices_of_leaders);

	return (pool_indices_of_cluster_members);
}//int** CMoGomea::clustering(double **objective_values_pool, int pool_size, int number_of_dimensions,



/**
 * Selects n points from a set of points. A
 * greedy heuristic is used to find a good
 * scattering of the selected points. First,
 * a point is selected with a maximum value
 * in a randomly selected dimension. The
 * remaining points are selected iteratively.
 * In each iteration, the point selected is
 * the one that maximizes the minimal distance
 * to the points selected so far.
 */
int *CMoGomea::greedyScatteredSubsetSelection(double **points, int number_of_points, int number_of_dimensions, int number_to_select)
{
	int     i, index_of_farthest, random_dimension_index, number_selected_so_far,
		*indices_left, *result;
	double *nn_distances, distance_of_farthest, value;

	if (number_to_select > number_of_points)
	{
		printf("\n");
		printf("Error: greedyScatteredSubsetSelection asked to select %d solutions from set of size %d.", number_to_select, number_of_points);
		printf("\n\n");

		exit(0);
	}

	//result = (int *)Malloc(number_to_select * sizeof(int));
	result = new int[number_to_select];

	//indices_left = (int *)Malloc(number_of_points * sizeof(int));
	indices_left = new int[number_of_points];
	for (i = 0; i < number_of_points; i++)
		indices_left[i] = i;

	/* Find the first point: maximum value in a randomly chosen dimension */
	random_dimension_index = randomInt(number_of_dimensions);

	index_of_farthest = 0;
	distance_of_farthest = points[indices_left[index_of_farthest]][random_dimension_index];
	for (i = 1; i < number_of_points; i++)
	{
		if (points[indices_left[i]][random_dimension_index] > distance_of_farthest)
		{
			index_of_farthest = i;
			distance_of_farthest = points[indices_left[i]][random_dimension_index];
		}
	}

	number_selected_so_far = 0;
	result[number_selected_so_far] = indices_left[index_of_farthest];
	indices_left[index_of_farthest] = indices_left[number_of_points - number_selected_so_far - 1];
	number_selected_so_far++;

	/* Then select the rest of the solutions: maximum minimum
	 * (i.e. nearest-neighbour) distance to so-far selected points */
	//nn_distances = (double *)Malloc(number_of_points * sizeof(double));
	nn_distances = new  double[number_of_points];
	for (i = 0; i < number_of_points - number_selected_so_far; i++)
		nn_distances[i] = distanceEuclidean(points[indices_left[i]], points[result[number_selected_so_far - 1]], number_of_dimensions);

	while (number_selected_so_far < number_to_select)
	{
		index_of_farthest = 0;
		distance_of_farthest = nn_distances[0];
		for (i = 1; i < number_of_points - number_selected_so_far; i++)
		{
			if (nn_distances[i] > distance_of_farthest)
			{
				index_of_farthest = i;
				distance_of_farthest = nn_distances[i];
			}
		}

		result[number_selected_so_far] = indices_left[index_of_farthest];
		indices_left[index_of_farthest] = indices_left[number_of_points - number_selected_so_far - 1];
		nn_distances[index_of_farthest] = nn_distances[number_of_points - number_selected_so_far - 1];
		number_selected_so_far++;

		for (i = 0; i < number_of_points - number_selected_so_far; i++)
		{
			value = distanceEuclidean(points[indices_left[i]], points[result[number_selected_so_far - 1]], number_of_dimensions);
			if (value < nn_distances[i])
				nn_distances[i] = value;
		}
	}

	free(nn_distances);
	free(indices_left);
	return(result);
}//int *CMoGomea::greedyScatteredSubsetSelection(double **points, int number_of_points, int number_of_dimensions, int number_to_select)



void CMoGomea::determineExtremeClusters()
{
	int i, j, index_best, a, b, c, *order;
	// find extreme clusters
	order = createRandomOrdering(number_of_objectives);

	for (i = 0; i < number_of_mixing_components; i++)
		which_extreme[i] = -1;  // not extreme cluster

	if (number_of_mixing_components > 1)
	{
		for (i = 0; i < number_of_objectives; i++)
		{
			index_best = -1;

			for (j = 0; j < number_of_mixing_components; j++)
			{
				if (optimization[order[i]] == MINIMIZATION)
				{
					if (((index_best == -1) || (objective_means_scaled[j][order[i]] < objective_means_scaled[index_best][order[i]])) &&
						(which_extreme[j] == -1))
						index_best = j;
				}
				else if (optimization[order[i]] == MAXIMIZATION)
				{
					if (((index_best == -1) || (objective_means_scaled[j][order[i]] > objective_means_scaled[index_best][order[i]])) &&
						(which_extreme[j] == -1))
						index_best = j;
				}
			}
			which_extreme[index_best] = order[i];
		}
	}

	free(order);
}//void CMoGomea::determineExtremeClusters()

void CMoGomea::initializeClusters()
{
	int i;
	//lt = (int ***)Malloc(number_of_mixing_components * sizeof(int **));
	lt = new int**[number_of_mixing_components];
	//lt_length = (int *)Malloc(number_of_mixing_components * sizeof(int));
	lt_length = new int[number_of_mixing_components];
	//lt_number_of_indices = (int **)Malloc(number_of_mixing_components * sizeof(int *));
	lt_number_of_indices = new int*[number_of_mixing_components];
	for (i = 0; i < number_of_mixing_components; i++)
	{
		lt[i] = NULL;
		lt_number_of_indices[i] = NULL;
		lt_length[i] = 0;
	}

	//which_extreme = (int*)Malloc(number_of_mixing_components * sizeof(int));
	which_extreme = new int[number_of_mixing_components];

	//objective_means_scaled = (double **)Malloc(number_of_mixing_components * sizeof(double *));
	objective_means_scaled = new double*[number_of_mixing_components];
	for (i = 0; i < number_of_mixing_components; i++)
		//objective_means_scaled[i] = (double *)Malloc(number_of_objectives * sizeof(double));
		objective_means_scaled[i] = new double[number_of_objectives];
}//void CMoGomea::initializeClusters()




void CMoGomea::ezilaitiniClusters()
{
	int i, j;

	if (lt == NULL)
		return;

	for (i = 0; i < number_of_mixing_components; i++)
	{
		if (lt[i] != NULL)
		{
			for (j = 0; j < lt_length[i]; j++)
				free(lt[i][j]);
			free(lt[i]);
			free(lt_number_of_indices[i]);
		}
	}

	free(lt); lt = NULL;
	free(lt_length);
	free(lt_number_of_indices);

	free(which_extreme);

	for (i = 0; i < number_of_mixing_components; i++)
		free(objective_means_scaled[i]);
	free(objective_means_scaled);
}//void CMoGomea::ezilaitiniClusters()



void CMoGomea::improveCurrentPopulation(void)
{
	int     i, j, k, j_min, cluster_index, objective_index, number_of_cluster,
		*sum_cluster, *clusters;
	double *objective_values_scaled,
		distance, distance_smallest;

	offspring_size = population_size;
	//offspring = (char**)Malloc(offspring_size * sizeof(char*));
	offspring = new char*[offspring_size];
	//objective_values_offspring = (double**)Malloc(offspring_size * sizeof(double*));
	objective_values_offspring = new double*[offspring_size];
	//constraint_values_offspring = (double*)Malloc(offspring_size * sizeof(double));
	constraint_values_offspring = new double[offspring_size];

	for (i = 0; i < offspring_size; i++)
	{
		//offspring[i] = (char*)Malloc(number_of_parameters * sizeof(char));
		offspring[i] = new char[number_of_parameters];
		//objective_values_offspring[i] = (double*)Malloc(number_of_objectives * sizeof(double));
		objective_values_offspring[i] = new double[number_of_objectives];
	}//for (i = 0; i < offspring_size; i++)


	//objective_values_scaled = (double *)Malloc(number_of_objectives * sizeof(double));
	objective_values_scaled = new double[number_of_objectives];
	//sum_cluster = (int*)Malloc(number_of_mixing_components * sizeof(int));
	sum_cluster = new int[number_of_mixing_components];

	for (i = 0; i < number_of_mixing_components; i++)
		sum_cluster[i] = 0;

	elitist_archive_front_changed = FALSE;
	for (i = 0; i < population_size; i++)
	{
		number_of_cluster = 0;
		//clusters = (int*)Malloc(number_of_mixing_components * sizeof(int));
		clusters = new int[number_of_mixing_components];
		for (j = 0; j < number_of_mixing_components; j++)
		{
			for (k = 0; k < population_cluster_sizes[j]; k++)
			{
				if (population_indices_of_cluster_members[j][k] == i)
				{
					clusters[number_of_cluster] = j;
					number_of_cluster++;
					break;
				}
			}
		}
		if (number_of_cluster > 0)
			cluster_index = clusters[randomInt(number_of_cluster)];
		else
		{
			for (j = 0; j < number_of_objectives; j++)
				objective_values_scaled[j] = objective_values[i][j] / objective_ranges[j];

			distance_smallest = -1;
			j_min = -1;
			for (j = 0; j < number_of_mixing_components; j++)
			{
				distance = distanceEuclidean(objective_values_scaled, objective_means_scaled[j], number_of_objectives);
				if ((distance_smallest < 0) || (distance < distance_smallest))
				{
					j_min = j;
					distance_smallest = distance;
				}

			}

			cluster_index = j_min;
		}

		sum_cluster[cluster_index]++;
		if (which_extreme[cluster_index] == -1)
		{
			performMultiObjectiveGenepoolOptimalMixing(cluster_index, population[i], objective_values[i], constraint_values[i],
				offspring[i], objective_values_offspring[i], &(constraint_values_offspring[i]));
		}
		else
		{
			objective_index = which_extreme[cluster_index];
			performSingleObjectiveGenepoolOptimalMixing(cluster_index, objective_index, population[i], objective_values[i], constraint_values[i],
				offspring[i], objective_values_offspring[i], &(constraint_values_offspring[i]));
		}
		free(clusters);
	}

	free(objective_values_scaled); free(sum_cluster);

	if (!elitist_archive_front_changed)
		t_NIS++;
	else
		t_NIS = 0;
}//void CMoGomea::improveCurrentPopulation(void)




void CMoGomea::copyValuesFromDonorToOffspring(char *solution, char *donor, int cluster_index, int linkage_group_index)
{
	int i, parameter_index;
	for (i = 0; i < lt_number_of_indices[cluster_index][linkage_group_index]; i++)
	{
		parameter_index = lt[cluster_index][linkage_group_index][i];
		solution[parameter_index] = donor[parameter_index];
	}//for (i = 0; i < lt_number_of_indices[cluster_index][linkage_group_index]; i++)
}//void CMoGomea::copyValuesFromDonorToOffspring(char *solution, char *donor, int cluster_index, int linkage_group_index)



void CMoGomea::copyFromAToB(char *solution_a, double *obj_a, double con_a, char *solution_b, double *obj_b, double *con_b)
{
	int i;
	for (i = 0; i < number_of_parameters; i++)
		solution_b[i] = solution_a[i];
	for (i = 0; i < number_of_objectives; i++)
		obj_b[i] = obj_a[i];
	*con_b = con_a;
}//void CMoGomea::copyFromAToB(char *solution_a, double *obj_a, double con_a, char *solution_b, double *obj_b, double *con_b)



void CMoGomea::mutateSolution(char *solution, int lt_factor_index, int cluster_index)
{
	double mutation_rate, prob;
	int i, parameter_index;

	if (use_pre_mutation == FALSE && use_pre_adaptive_mutation == FALSE)
		return;

	mutation_rate = 0.0;
	if (use_pre_mutation == TRUE)
		mutation_rate = 1.0 / ((double)number_of_parameters);
	else if (use_pre_adaptive_mutation == TRUE)
		mutation_rate = 1.0 / ((double)lt_number_of_indices[cluster_index][lt_factor_index]);


	for (i = 0; i < lt_number_of_indices[cluster_index][lt_factor_index]; i++)
	{
		prob = randomRealUniform01();
		if (prob < mutation_rate)
		{
			parameter_index = lt[cluster_index][lt_factor_index][i];
			if (solution[parameter_index] == 0)
				solution[parameter_index] = 1;
			else
				solution[parameter_index] = 0;
		}
	}

}//void CMoGomea::mutateSolution(char *solution, int lt_factor_index, int cluster_index)



/**
 * Multi-objective Gene-pool Optimal Mixing
 * Construct an offspring from a parent solution in a middle-region cluster.
 */
void CMoGomea::performMultiObjectiveGenepoolOptimalMixing(int cluster_index, char *parent, double *parent_obj, double parent_con,
	char *result, double *obj, double *con)
{
	char   *backup, *donor, is_unchanged, changed, is_improved, is_new_nondominated_point, is_dominated_by_archive;
	int     i, j, index, donor_index, position_of_existed_member, *order, linkage_group_index, number_of_linkage_sets;
	double  *obj_backup, con_backup;

	/* Clone the parent solution. */
	copyFromAToB(parent, parent_obj, parent_con, result, obj, con);

	/* Create a backup version of the parent solution. */
	//backup = (char *)Malloc(number_of_parameters * sizeof(char));
	backup = new char[number_of_parameters];
	//obj_backup = (double *)Malloc(number_of_objectives * sizeof(double));
	obj_backup = new double[number_of_objectives];
	copyFromAToB(result, obj, *con, backup, obj_backup, &con_backup);

	number_of_linkage_sets = lt_length[cluster_index] - 1; /* Remove root from the linkage tree. */
	order = createRandomOrdering(number_of_linkage_sets);

	/* Traverse the linkage tree for Gene-pool Optimal Mixing */
	changed = FALSE;
	for (i = 0; i < number_of_linkage_sets; i++)
	{
		linkage_group_index = order[i];

		donor_index = randomInt(population_cluster_sizes[cluster_index]);

		donor = population[population_indices_of_cluster_members[cluster_index][donor_index]];
		copyValuesFromDonorToOffspring(result, donor, cluster_index, linkage_group_index);
		mutateSolution(result, linkage_group_index, cluster_index);

		/* Check if the new intermediate solution is different from the previous state. */
		is_unchanged = TRUE;
		for (j = 0; j < lt_number_of_indices[cluster_index][linkage_group_index]; j++)
		{
			if (backup[lt[cluster_index][linkage_group_index][j]] != result[lt[cluster_index][linkage_group_index][j]])
			{
				is_unchanged = FALSE;
				break;
			}
		}

		if (is_unchanged == FALSE)
		{
			is_improved = FALSE;
			evaluateIndividual(result, obj, con, NOT_EXTREME_CLUSTER);
			updateElitistArchiveWithReplacementOfExistedMember(result, obj, *con, &is_new_nondominated_point, &is_dominated_by_archive);

			/* Check for weak Pareto domination. */
			if (constraintWeaklyParetoDominates(obj, *con, obj_backup, con_backup))
				is_improved = TRUE;

			/* Check if the new intermediate solution is dominated by any solution in the archive.
				Note that if the new intermediate solution is a point in the archive, it is NOT considered to be dominated by the archive.*/
			if (!is_dominated_by_archive)
				is_improved = TRUE;

			if (is_improved)
			{
				changed = TRUE;
				copyFromAToB(result, obj, *con, backup, obj_backup, &con_backup);
			}
			else
				copyFromAToB(backup, obj_backup, con_backup, result, obj, con);
		}
	}
	free(order);

	/* Forced Improvement */
	if ((!changed) || (t_NIS > (1 + floor(log10(population_size)))))
	{
		changed = FALSE;
		order = createRandomOrdering(number_of_linkage_sets);
		/* Perform another round of Gene-pool Optimal Mixing with the donors randomly selected from the archive. */
		for (i = 0; i < number_of_linkage_sets; i++)
		{
			donor_index = randomInt(elitist_archive_size);
			linkage_group_index = order[i];
			copyValuesFromDonorToOffspring(result, elitist_archive[donor_index], cluster_index, linkage_group_index);
			mutateSolution(result, linkage_group_index, cluster_index);

			/* Check if the new intermediate solution is different from the previous state. */
			is_unchanged = TRUE;
			for (j = 0; j < lt_number_of_indices[cluster_index][linkage_group_index]; j++)
			{
				if (backup[lt[cluster_index][linkage_group_index][j]] != result[lt[cluster_index][linkage_group_index][j]])
				{
					is_unchanged = FALSE;
					break;
				}
			}

			if (is_unchanged == FALSE)
			{
				is_improved = FALSE;

				evaluateIndividual(result, obj, con, NOT_EXTREME_CLUSTER);
				updateElitistArchiveWithReplacementOfExistedMember(result, obj, *con, &is_new_nondominated_point, &is_dominated_by_archive);

				/* Check for (strict) Pareto domination. */
				if (constraintParetoDominates(obj, *con, obj_backup, con_backup))
					is_improved = TRUE;

				/* Check if a truly new non-dominated solution is created. */
				if (is_new_nondominated_point)
					is_improved = TRUE;

				if (is_improved)
				{
					changed = TRUE;
					copyFromAToB(result, obj, *con, backup, obj_backup, &con_backup);
					break;
				}
				else
					copyFromAToB(backup, obj_backup, con_backup, result, obj, con);
			}
		}
		free(order);

		if (!changed)
		{
			donor_index = randomInt(elitist_archive_size);

			copyFromAToB(elitist_archive[donor_index], elitist_archive_objective_values[donor_index],
				elitist_archive_constraint_values[donor_index],
				result, obj, con);
		}
	}

	free(backup); free(obj_backup);
}//void CMoGomea::performMultiObjectiveGenepoolOptimalMixing(int cluster_index, char *parent, double *parent_obj, double parent_con,



/**
 * Single-objective Gene-pool Optimal Mixing
 * Construct an offspring from a parent solution in an extreme-region cluster.
 */
void CMoGomea::performSingleObjectiveGenepoolOptimalMixing(int cluster_index, int objective_index,
	char *parent, double *parent_obj, double parent_con,
	char *result, double *obj, double *con)
{
	char   *backup, *donor, *elitist_copy, is_unchanged, changed, is_improved, is_new_nondominated_point, is_dominated_by_archive;
	int     i, j, index, donor_index, number_of_linkage_sets, linkage_group_index, *order;
	double  *obj_backup, con_backup;

	/* Clone the parent solution. */
	copyFromAToB(parent, parent_obj, parent_con, result, obj, con);

	/* Create a backup version of the parent solution. */
	//backup = (char *)Malloc(number_of_parameters * sizeof(char));
	backup = new char[number_of_parameters];
	//obj_backup = (double *)Malloc(number_of_objectives * sizeof(double));
	obj_backup = new double[number_of_objectives];
	copyFromAToB(result, obj, *con, backup, obj_backup, &con_backup);

	number_of_linkage_sets = lt_length[cluster_index] - 1; /* Remove root from the linkage tree. */

	order = createRandomOrdering(number_of_linkage_sets);

	/* Traverse the linkage tree for Gene-pool Optimal Mixing */
	changed = FALSE;
	for (i = 0; i < number_of_linkage_sets; i++)
	{
		linkage_group_index = order[i];
		donor_index = randomInt(population_cluster_sizes[cluster_index]);

		donor = population[population_indices_of_cluster_members[cluster_index][donor_index]];
		copyValuesFromDonorToOffspring(result, donor, cluster_index, linkage_group_index);
		mutateSolution(result, linkage_group_index, cluster_index);

		/* Check if the new intermediate solution is different from the previous state. */
		is_unchanged = TRUE;
		for (j = 0; j < lt_number_of_indices[cluster_index][linkage_group_index]; j++)
		{
			if (backup[lt[cluster_index][linkage_group_index][j]] != result[lt[cluster_index][linkage_group_index][j]])
			{
				is_unchanged = FALSE;
				break;
			}
		}

		if (is_unchanged == FALSE)
		{
			is_improved = FALSE;
			evaluateIndividual(result, obj, con, objective_index);
			updateElitistArchiveWithReplacementOfExistedMember(result, obj, *con, &is_new_nondominated_point, &is_dominated_by_archive);

			if (betterFitness(obj, *con, obj_backup, con_backup, objective_index) ||
				equalFitness(obj, *con, obj_backup, con_backup, objective_index))
				is_improved = TRUE;

			if (is_improved)
			{
				changed = TRUE;
				copyFromAToB(result, obj, *con, backup, obj_backup, &con_backup);
			}
			else
				copyFromAToB(backup, obj_backup, con_backup, result, obj, con);
		}
	}
	free(order);

	//elitist_copy = (char*)Malloc(number_of_parameters * sizeof(char));
	elitist_copy = new char[number_of_parameters];
	/* Forced Improvement*/
	if ((!changed) || (t_NIS > (1 + floor(log10(population_size)))))
	{
		changed = FALSE;
		/* Find in the archive the solution having the best value in the corresponding objective. */
		donor_index = 0;
		for (j = 0; j < elitist_archive_size; j++)
		{
			if (optimization[objective_index] == MINIMIZATION)
			{
				if (elitist_archive_objective_values[j][objective_index] < elitist_archive_objective_values[donor_index][objective_index])
					donor_index = j;
			}
			else if (optimization[objective_index] == MAXIMIZATION)
			{
				if (elitist_archive_objective_values[j][objective_index] > elitist_archive_objective_values[donor_index][objective_index])
					donor_index = j;
			}
		}

		for (j = 0; j < number_of_parameters; j++)
			elitist_copy[j] = elitist_archive[donor_index][j];

		/* Perform Gene-pool Optimal Mixing with the single-objective best-found solution as the donor. */
		order = createRandomOrdering(number_of_linkage_sets);
		for (i = 0; i < number_of_linkage_sets; i++)
		{
			linkage_group_index = order[i];
			copyValuesFromDonorToOffspring(result, elitist_copy, cluster_index, linkage_group_index);
			mutateSolution(result, linkage_group_index, cluster_index);

			/* Check if the new intermediate solution is different from the previous state. */
			is_unchanged = TRUE;
			for (j = 0; j < lt_number_of_indices[cluster_index][linkage_group_index]; j++)
			{
				if (backup[lt[cluster_index][linkage_group_index][j]] != result[lt[cluster_index][linkage_group_index][j]])
				{
					is_unchanged = FALSE;
					break;
				}
			}

			if (is_unchanged == FALSE)
			{
				is_improved = FALSE;
				evaluateIndividual(result, obj, con, objective_index);
				updateElitistArchiveWithReplacementOfExistedMember(result, obj, *con, &is_new_nondominated_point, &is_dominated_by_archive);

				/* Check if strict improvement in the corresponding objective. */
				if (betterFitness(obj, *con, obj_backup, con_backup, objective_index))
					is_improved = TRUE;

				if (is_improved == TRUE)
				{
					changed = TRUE;
					copyFromAToB(result, obj, *con, backup, obj_backup, &con_backup);
					break;
				}
				else
					copyFromAToB(backup, obj_backup, con_backup, result, obj, con);
			}
		}
		free(order);

		if (!changed)
		{
			donor_index = 0;
			for (j = 0; j < elitist_archive_size; j++)
			{
				if (optimization[objective_index] == MINIMIZATION)
				{
					if (elitist_archive_objective_values[j][objective_index] < elitist_archive_objective_values[donor_index][objective_index])
						donor_index = j;
				}
				else if (optimization[objective_index] == MAXIMIZATION)
				{
					if (elitist_archive_objective_values[j][objective_index] > elitist_archive_objective_values[donor_index][objective_index])
						donor_index = j;
				}
			}

			copyFromAToB(elitist_archive[donor_index], elitist_archive_objective_values[donor_index],
				elitist_archive_constraint_values[donor_index],
				result, obj, con);
		}
	}

	free(backup); free(obj_backup); free(elitist_copy);
}//void CMoGomea::performSingleObjectiveGenepoolOptimalMixing(int cluster_index, int objective_index,

/**
 * Determines the solutions that finally survive the generation (offspring only).
 */
void CMoGomea::selectFinalSurvivors()
{
	int i, j;

	for (i = 0; i < population_size; i++)
	{
		for (j = 0; j < number_of_parameters; j++)
			population[i][j] = offspring[i][j];
		for (j = 0; j < number_of_objectives; j++)
			objective_values[i][j] = objective_values_offspring[i][j];
		constraint_values[i] = constraint_values_offspring[i];
	}
}//void CMoGomea::selectFinalSurvivors()



void CMoGomea::freeAuxiliaryPopulations()
{
	int i, k;

	if (population_indices_of_cluster_members != NULL)
	{
		for (k = 0; k < number_of_mixing_components; k++)
			free(population_indices_of_cluster_members[k]);
		free(population_indices_of_cluster_members);
		population_indices_of_cluster_members = NULL;
		free(population_cluster_sizes);
	}

	if (offspring != NULL)
	{
		for (i = 0; i < offspring_size; i++)
		{
			free(offspring[i]);
			free(objective_values_offspring[i]);
		}
		free(offspring);
		free(objective_values_offspring);
		free(constraint_values_offspring);
		offspring = NULL;
	}

	ezilaitiniClusters();
};//void CMoGomea::freeAuxiliaryPopulations()


























bool CMoGomea::bRunIteration_dummy(uint32_t iIterationNumber, time_t tStartTime)
{
	CString  s_buf;
	double  d_time_passed;

	c_time_counter.bGetTimePassed(&d_time_passed);


	CBinaryCoding c_genotype(i_templ_length);

	for (int ii = 0; ii < i_templ_length; ii++)
		c_genotype.piGetBits()[ii] = RandUtils::iRandNumber(0, 1);


	CBinaryMultiObjectiveProblem  *pc_multi_problem;
	pc_multi_problem = (CBinaryMultiObjectiveProblem *)pc_problem->pcGetEvaluation();


	double  d_fit;
	d_fit = pc_multi_problem->dEvaluate(&c_genotype);


	//s_buf.Format("iteration: %d  pop: %d; GenSize:%d PFsize:%d HyperVolume: %.8lf InvGenDist: %.8lf [time:%.2lf]", iIterationNumber, v_population.size(), v_non_dominated_pfs.at(0).at(0)->pc_genotype->iGetNumberOfBits(), (int)dPFQualityPointNum(), dPFQualityHyperVolume(), dPFQualityInverseGenerationalDistance(), d_time_passed);
	s_buf.Format("iteration: %d  PFsize:%d HyperVolume: %.8lf InvGenDist: %.8lf [time:%.2lf]", iIterationNumber, (int)dPFQualityPointNum(), dPFQualityHyperVolume(), dPFQualityInverseGenerationalDistance(), d_time_passed);
	pc_log->vPrintLine(s_buf, true);

	b_update_best_individual(iIterationNumber, tStartTime, c_genotype.piGetBits(), d_fit);

	return(true);
};//bool CMoGomea::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)



CString  CMoGomea::sAdditionalSummaryInfo()
{
	return("test");
};//CString  CMoGomea::sAdditionalSummaryInfo()





void CMoGomea::writeGenerationalStatistics(void)
{
	return;
	int     i;
	char    string[1000];
	FILE   *file;

	file = NULL;
	if ((number_of_generations == 0 && population_id == 0) ||
		(number_of_generations == 0 && population_id == -1))
	{
		file = fopen("statistics.dat", "w");

		sprintf(string, "# Generation  Population  Evaluations   [ Cluster_Index ]\n");
		fputs(string, file);
	}
	else
		file = fopen("statistics.dat", "a");

	sprintf(string, "  %10d %10d %11d     [ ", number_of_generations, population_size, number_of_evaluations);
	fputs(string, file);

	for (i = 0; i < number_of_mixing_components; i++)
	{
		sprintf(string, "%4d", i);
		fputs(string, file);
		if (i < number_of_mixing_components - 1)
		{
			sprintf(string, " ");
			fputs(string, file);
		}
	}
	sprintf(string, " ]\n");
	fputs(string, file);

	fclose(file);

	freeAuxiliaryPopulations();
}//void CMoGomea::writeGenerationalStatistics(void)



void CMoGomea::adaptObjectiveDiscretization(void)
{
	int    i, j, k, na, nb, nc, elitist_archive_size_target_lower_bound, elitist_archive_size_target_upper_bound;
	double low, high, *elitist_archive_objective_ranges;

	elitist_archive_size_target_lower_bound = (int)(0.75*elitist_archive_size_target);
	elitist_archive_size_target_upper_bound = (int)(1.25*elitist_archive_size_target);

	if (objective_discretization_in_effect && (elitist_archive_size < elitist_archive_size_target_lower_bound))
		objective_discretization_in_effect = 0;

	if (elitist_archive_size > elitist_archive_size_target_upper_bound)
	{
		objective_discretization_in_effect = 1;

		//elitist_archive_objective_ranges = (double *)Malloc(number_of_objectives * sizeof(double));
		elitist_archive_objective_ranges = new double[number_of_objectives];
		for (j = 0; j < number_of_objectives; j++)
		{
			low = elitist_archive_objective_values[0][j];
			high = elitist_archive_objective_values[0][j];

			for (i = 0; i < elitist_archive_size; i++)
			{
				if (elitist_archive_objective_values[i][j] < low)
					low = elitist_archive_objective_values[i][j];
				if (elitist_archive_objective_values[i][j] > high)
					high = elitist_archive_objective_values[i][j];
			}

			elitist_archive_objective_ranges[j] = high - low;
		}

		na = 1;
		nb = (int)pow(2.0, 25.0);

		for (k = 0; k < 25; k++)
		{
			nc = (na + nb) / 2;
			for (i = 0; i < number_of_objectives; i++)
				objective_discretization[i] = elitist_archive_objective_ranges[i] / ((double)nc);

			/* Restore the original elitist archive after the first cycle in this loop */
			if (k > 0)
			{
				elitist_archive_size = 0;
				for (i = 0; i < elitist_archive_copy_size; i++)
					addToElitistArchive(elitist_archive_copy[i], elitist_archive_copy_objective_values[i], elitist_archive_copy_constraint_values[i]);
			}

			/* Copy the entire elitist archive */
			if (elitist_archive_copy != NULL)
			{
				for (i = 0; i < elitist_archive_copy_size; i++)
				{
					free(elitist_archive_copy[i]);
					free(elitist_archive_copy_objective_values[i]);
				}
				free(elitist_archive_copy);
				free(elitist_archive_copy_objective_values);
				free(elitist_archive_copy_constraint_values);
			}

			elitist_archive_copy_size = elitist_archive_size;
			//elitist_archive_copy = (char **)Malloc(elitist_archive_copy_size * sizeof(char *));
			elitist_archive_copy = new char*[elitist_archive_copy_size];
			//elitist_archive_copy_objective_values = (double **)Malloc(elitist_archive_copy_size * sizeof(double *));
			elitist_archive_copy_objective_values = new double*[elitist_archive_copy_size];
			//elitist_archive_copy_constraint_values = (double *)Malloc(elitist_archive_copy_size * sizeof(double));
			elitist_archive_copy_constraint_values = new double[elitist_archive_copy_size];

			for (i = 0; i < elitist_archive_copy_size; i++)
			{
				//elitist_archive_copy[i] = (char *)Malloc(number_of_parameters * sizeof(char));
				elitist_archive_copy[i] = new char[number_of_parameters];
				//elitist_archive_copy_objective_values[i] = (double *)Malloc(number_of_objectives * sizeof(double));
				elitist_archive_copy_objective_values[i] = new double[number_of_objectives];
			}
			for (i = 0; i < elitist_archive_copy_size; i++)
			{
				for (j = 0; j < number_of_parameters; j++)
					elitist_archive_copy[i][j] = elitist_archive[i][j];
				for (j = 0; j < number_of_objectives; j++)
					elitist_archive_copy_objective_values[i][j] = elitist_archive_objective_values[i][j];
				elitist_archive_copy_constraint_values[i] = elitist_archive_constraint_values[i];
			}

			/* Clear the elitist archive */
			elitist_archive_size = 0;

			/* Rebuild the elitist archive */
			for (i = 0; i < elitist_archive_copy_size; i++)
				updateElitistArchive(elitist_archive_copy[i], elitist_archive_copy_objective_values[i], elitist_archive_copy_constraint_values[i]);

			if (elitist_archive_size <= elitist_archive_size_target_lower_bound)
				na = nc;
			else
				nb = nc;
		}

		free(elitist_archive_objective_ranges);
	}
}//void CMoGomea::adaptObjectiveDiscretization(void)



/**
 * Returns a random integer, distributed uniformly between 0 and maximum.
 */
int CMoGomea::randomInt(int maximum)
{
	return(RandUtils::iRandNumber(0, maximum-1));

	/*int result;
	result = (int)(((double)maximum)*randomRealUniform01());
	return(result);*/
}//int CMoGomea::randomInt(int maximum)


double CMoGomea::randomRealUniform01(void)
{
	return(RandUtils::dRandNumber(0,1));

	/*int64_t n26, n27;
	double  result;

	random_seed_changing = (random_seed_changing * 0x5DEECE66DLLU + 0xBLLU) & ((1LLU << 48) - 1);
	n26 = (int64_t)(random_seed_changing >> (48 - 26));
	random_seed_changing = (random_seed_changing * 0x5DEECE66DLLU + 0xBLLU) & ((1LLU << 48) - 1);
	n27 = (int64_t)(random_seed_changing >> (48 - 27));
	result = (((int64_t)n26 << 27) + n27) / ((double)(1LLU << 53));

	return(result);*/
}//double CMoGomea::randomRealUniform01(void)


double CMoGomea::distanceEuclidean(double *x, double *y, int number_of_dimensions)
{
	int    i;
	double value, result;

	result = 0.0;
	for (i = 0; i < number_of_dimensions; i++)
	{
		value = y[i] - x[i];
		result += value * value;
	}
	result = sqrt(result);

	return(result);
}//double CMoGomea::distanceEuclidean(double *x, double *y, int number_of_dimensions)



int* CMoGomea::createRandomOrdering(int size_of_the_set)
{
	int *order, a, b, c, i;

	//order = (int *)Malloc(size_of_the_set * sizeof(int));
	order = new int[size_of_the_set];
	for (i = 0; i < size_of_the_set; i++)
		order[i] = i;
	for (i = 0; i < size_of_the_set; i++)
	{
		a = randomInt(size_of_the_set);
		b = randomInt(size_of_the_set);
		c = order[a];
		order[a] = order[b];
		order[b] = c;
	}

	return order;
}//int* CMoGomea::createRandomOrdering(int size_of_the_set)



/**
 * Sorts an array of doubles and returns the sort-order (small to large).
 */
int *CMoGomea::mergeSort(double *array, int array_size)
{
	int i, *sorted, *tosort;

	//sorted = (int *)Malloc(array_size * sizeof(int));
	sorted = new int[array_size];
	//tosort = (int *)Malloc(array_size * sizeof(int));
	tosort = new int[array_size];
	for (i = 0; i < array_size; i++)
		tosort[i] = i;

	if (array_size == 1)
		sorted[0] = 0;
	else
		mergeSortWithinBounds(array, sorted, tosort, 0, array_size - 1);

	free(tosort);

	return(sorted);
}//int *CMoGomea::mergeSort(double *array, int array_size)


/**
 * Subroutine of merge sort, sorts the part of the array between p and q.
 */
void CMoGomea::mergeSortWithinBounds(double *array, int *sorted, int *tosort, int p, int q)
{
	int r;

	if (p < q)
	{
		r = (p + q) / 2;
		mergeSortWithinBounds(array, sorted, tosort, p, r);
		mergeSortWithinBounds(array, sorted, tosort, r + 1, q);
		mergeSortMerge(array, sorted, tosort, p, r + 1, q);
	}
}//void CMoGomea::mergeSortWithinBounds(double *array, int *sorted, int *tosort, int p, int q)


/**
 * Subroutine of merge sort, merges the results of two sorted parts.
 */
void CMoGomea::mergeSortMerge(double *array, int *sorted, int *tosort, int p, int r, int q)
{
	int i, j, k, first;

	i = p;
	j = r;
	for (k = p; k <= q; k++)
	{
		first = 0;
		if (j <= q)
		{
			if (i < r)
			{
				if (array[tosort[i]] < array[tosort[j]])
					first = 1;
			}
		}
		else
			first = 1;

		if (first)
		{
			sorted[k] = tosort[i];
			i++;
		}
		else
		{
			sorted[k] = tosort[j];
			j++;
		}
	}

	for (k = p; k <= q; k++)
		tosort[k] = sorted[k];
}//void CMoGomea::mergeSortMerge(double *array, int *sorted, int *tosort, int p, int r, int q)






void CMoGomea::evaluateIndividual(char *solution, double *obj, double *con, int objective_index_of_extreme_cluster)
{
	number_of_evaluations++;
	if (population_id != -1)
		array_of_number_of_evaluations_per_population[population_id] += 1;



	
	for (int  i = 0; i < number_of_parameters; i++)
	{
		if (solution[i] == 0)
			pc_genotype->piGetBits()[i] = 0;
		else if (solution[i] == 1)
			pc_genotype->piGetBits()[i] = 1;
	}//for (int  i = 0; i < number_of_parameters; i++)



	CBinaryMultiObjectiveProblem  *pc_multi_problem;
	pc_multi_problem = (CBinaryMultiObjectiveProblem *)pc_problem->pcGetEvaluation();


	vector  <double>  v_objectives;
	pc_multi_problem->vEvaluateParetoFront(&v_objectives, pc_genotype);


	if (v_objectives.size() != 2)
	{
		CString  s_buf;

		s_buf.Format("void CMoGomea::evaluateIndividual(char *solution, double *obj, double *con, int objective_index_of_extreme_cluster) :\n if (v_objectives.size() != 2)");
		Tools::vShow(s_buf);

		return;
	}//if (v_objectives.size() != 2)


	obj[0] = v_objectives.at(0);
	obj[1] = v_objectives.at(1);


	

	/*switch (problem_index)
	{
	case ZEROMAX_ONEMAX: onemaxProblemEvaluation(solution, obj, con, objective_index_of_extreme_cluster); break;
	case TRAP5: trap5ProblemEvaluation(solution, obj, con, objective_index_of_extreme_cluster); break;
	case LOTZ: lotzProblemEvaluation(solution, obj, con, objective_index_of_extreme_cluster); break;
	case KNAPSACK: knapsackProblemEvaluation(solution, obj, con, objective_index_of_extreme_cluster); break;
	case MAXCUT: maxcutProblemEvaluation(solution, obj, con, objective_index_of_extreme_cluster); break;
	default:
		printf("Cannot evaluate this problem!\n");
		exit(1);
	}*/

	logElitistArchiveAtSpecificPoints();
}//void CMoGomea::evaluateIndividual(char *solution, double *obj, double *con, int objective_index_of_extreme_cluster)


char CMoGomea::betterFitness(double *objective_value_x, double constraint_value_x, double *objective_value_y, double constraint_value_y, int objective_index)
{
	short result;

	result = FALSE;

	if (constraint_value_x > 0) /* x is infeasible */
	{
		if (constraint_value_y > 0) /* Both are infeasible */
		{
			if (constraint_value_x < constraint_value_y)
				result = TRUE;
		}
	}
	else
	{
		if (constraint_value_y > 0)
			result = TRUE;
		else
		{
			if (optimization[objective_index] == MINIMIZATION)
			{
				if (objective_value_x[objective_index] < objective_value_y[objective_index])
					result = TRUE;
			}
			else if (optimization[objective_index] == MAXIMIZATION)
			{
				if (objective_value_x[objective_index] > objective_value_y[objective_index])
					result = TRUE;
			}
		}
	}

	return (result);
}//char CMoGomea::betterFitness(double *objective_value_x, double constraint_value_x, double *objective_value_y, double constraint_value_y, int objective_index)


char CMoGomea::equalFitness(double *objective_value_x, double constraint_value_x, double *objective_value_y, double constraint_value_y, int objective_index)
{
	short result;

	result = FALSE;

	if (constraint_value_x > 0) /* x is infeasible */
	{
		if (constraint_value_y > 0) /* Both are infeasible */
		{
			if (constraint_value_x == constraint_value_y)
				result = TRUE;
		}
	}
	else
	{
		if (constraint_value_y == 0)
		{
			if (objective_value_x[objective_index] == objective_value_y[objective_index])
				result = TRUE;
		}
	}

	return (result);
}//char CMoGomea::equalFitness(double *objective_value_x, double constraint_value_x, double *objective_value_y, double constraint_value_y, int objective_index)



void CMoGomea::logElitistArchiveAtSpecificPoints()
{
	//we do not log here anymore
	return;

	if (number_of_evaluations%log_progress_interval == 0)
		writeCurrentElitistArchive(FALSE);
}//void CMoGomea::logElitistArchiveAtSpecificPoints()




/**
 * Updates the elitist archive by offering a new solution
 * to possibly be added to the archive. If there are no
 * solutions in the archive yet, the solution is added.
 * Solution A is always dominated by solution B that is
 * in the same domination-box if B dominates A or A and
 * B do not dominate each other. If the solution is not
 * dominated, it is added to the archive and all solutions
 * dominated by the new solution, are purged from the archive.
 */
void CMoGomea::updateElitistArchive(char *solution, double *solution_objective_values, double solution_constraint_value)
{
	short is_dominated_itself;
	int   i, *indices_dominated, number_of_solutions_dominated;

	if (elitist_archive_size == 0)
		addToElitistArchive(solution, solution_objective_values, solution_constraint_value);
	else
	{
		//indices_dominated = (int *)Malloc(elitist_archive_size * sizeof(int));
		indices_dominated = new int[elitist_archive_size];
		number_of_solutions_dominated = 0;
		is_dominated_itself = 0;
		for (i = 0; i < elitist_archive_size; i++)
		{
			if (constraintParetoDominates(elitist_archive_objective_values[i], elitist_archive_constraint_values[i], solution_objective_values, solution_constraint_value))
				is_dominated_itself = 1;
			else
			{
				if (!constraintParetoDominates(solution_objective_values, solution_constraint_value, elitist_archive_objective_values[i], elitist_archive_constraint_values[i]))
				{
					if (sameObjectiveBox(elitist_archive_objective_values[i], solution_objective_values))
						is_dominated_itself = 1;
				}
			}

			if (is_dominated_itself)
				break;
		}

		if (!is_dominated_itself)
		{
			for (i = 0; i < elitist_archive_size; i++)
			{
				if (constraintParetoDominates(solution_objective_values, solution_constraint_value, elitist_archive_objective_values[i], elitist_archive_constraint_values[i]))
				{
					indices_dominated[number_of_solutions_dominated] = i;
					number_of_solutions_dominated++;
				}
			}

			if (number_of_solutions_dominated > 0)
				removeFromElitistArchive(indices_dominated, number_of_solutions_dominated);

			addToElitistArchive(solution, solution_objective_values, solution_constraint_value);
		}

		free(indices_dominated);
	}
}//void CMoGomea::updateElitistArchive(char *solution, double *solution_objective_values, double solution_constraint_value)




void CMoGomea::addToElitistArchive(char *solution, double *solution_objective_values, double solution_constraint_value)
{
	int      i, j, elitist_archive_capacity_new;
	char **elitist_archive_new;
	double **elitist_archive_objective_values_new;
	double *elitist_archive_constraint_values_new;

	if (elitist_archive_capacity == elitist_archive_size)
	{
		elitist_archive_capacity_new = elitist_archive_capacity * 2 + 1;
		//elitist_archive_new = (char **)Malloc(elitist_archive_capacity_new * sizeof(char *));
		elitist_archive_new = new char*[elitist_archive_capacity_new];
		//elitist_archive_objective_values_new = (double **)Malloc(elitist_archive_capacity_new * sizeof(double *));
		elitist_archive_objective_values_new = new double*[elitist_archive_capacity_new];
		//elitist_archive_constraint_values_new = (double *)Malloc(elitist_archive_capacity_new * sizeof(double));
		elitist_archive_constraint_values_new = new double[elitist_archive_capacity_new];

		for (i = 0; i < elitist_archive_capacity_new; i++)
		{
			//elitist_archive_new[i] = (char *)Malloc(number_of_parameters * sizeof(char));
			elitist_archive_new[i] = new char[number_of_parameters];
			//elitist_archive_objective_values_new[i] = (double *)Malloc(number_of_objectives * sizeof(double));
			elitist_archive_objective_values_new[i] = new double[number_of_objectives];
		}

		for (i = 0; i < elitist_archive_size; i++)
		{
			for (j = 0; j < number_of_parameters; j++)
				elitist_archive_new[i][j] = elitist_archive[i][j];
			for (j = 0; j < number_of_objectives; j++)
				elitist_archive_objective_values_new[i][j] = elitist_archive_objective_values[i][j];
			elitist_archive_constraint_values_new[i] = elitist_archive_constraint_values[i];
		}

		for (i = 0; i < elitist_archive_capacity; i++)
		{
			free(elitist_archive[i]);
			free(elitist_archive_objective_values[i]);
		}
		free(elitist_archive);
		free(elitist_archive_objective_values);
		free(elitist_archive_constraint_values);

		elitist_archive_capacity = elitist_archive_capacity_new;
		elitist_archive = elitist_archive_new;
		elitist_archive_objective_values = elitist_archive_objective_values_new;
		elitist_archive_constraint_values = elitist_archive_constraint_values_new;
	}

	for (j = 0; j < number_of_parameters; j++)
		elitist_archive[elitist_archive_size][j] = solution[j];
	for (j = 0; j < number_of_objectives; j++)
		elitist_archive_objective_values[elitist_archive_size][j] = solution_objective_values[j];
	elitist_archive_constraint_values[elitist_archive_size] = solution_constraint_value; // Notice here //

	elitist_archive_size++;
}//void CMoGomea::addToElitistArchive(char *solution, double *solution_objective_values, double solution_constraint_value)




void CMoGomea::removeFromElitistArchive(int *indices, int number_of_indices)
{
	int      i, j, elitist_archive_size_new;
	char **elitist_archive_new;
	double **elitist_archive_objective_values_new;
	double *elitist_archive_constraint_values_new;

	//elitist_archive_new = (char**)Malloc(elitist_archive_capacity * sizeof(char *));
	elitist_archive_new = new char*[elitist_archive_capacity];
	//elitist_archive_objective_values_new = (double **)Malloc(elitist_archive_capacity * sizeof(double *));
	elitist_archive_objective_values_new = new double*[elitist_archive_capacity];
	//elitist_archive_constraint_values_new = (double *)Malloc(elitist_archive_capacity * sizeof(double));
	elitist_archive_constraint_values_new = new double[elitist_archive_capacity];

	for (i = 0; i < elitist_archive_capacity; i++)
	{
		//elitist_archive_new[i] = (char *)Malloc(number_of_parameters * sizeof(char));
		elitist_archive_new[i] = new char[number_of_parameters];
		//elitist_archive_objective_values_new[i] = (double *)Malloc(number_of_objectives * sizeof(double));
		elitist_archive_objective_values_new[i] = new double[number_of_objectives];
	}

	elitist_archive_size_new = 0;
	for (i = 0; i < elitist_archive_size; i++)
	{
		if (!isInListOfIndices(i, indices, number_of_indices))
		{
			for (j = 0; j < number_of_parameters; j++)
				elitist_archive_new[elitist_archive_size_new][j] = elitist_archive[i][j];
			for (j = 0; j < number_of_objectives; j++)
				elitist_archive_objective_values_new[elitist_archive_size_new][j] = elitist_archive_objective_values[i][j];
			elitist_archive_constraint_values_new[elitist_archive_size_new] = elitist_archive_constraint_values[i];

			elitist_archive_size_new++;
		}
	}

	for (i = 0; i < elitist_archive_capacity; i++)
	{
		free(elitist_archive[i]);
		free(elitist_archive_objective_values[i]);
	}
	free(elitist_archive);
	free(elitist_archive_objective_values);
	free(elitist_archive_constraint_values);

	elitist_archive_size = elitist_archive_size_new;
	elitist_archive = elitist_archive_new;
	elitist_archive_objective_values = elitist_archive_objective_values_new;
	elitist_archive_constraint_values = elitist_archive_constraint_values_new;
}//void CMoGomea::removeFromElitistArchive(int *indices, int number_of_indices)


short CMoGomea::isInListOfIndices(int index, int *indices, int number_of_indices)
{
	int i;

	for (i = 0; i < number_of_indices; i++)
		if (indices[i] == index)
			return(1);

	return(0);
}//short CMoGomea::isInListOfIndices(int index, int *indices, int number_of_indices)



void CMoGomea::updateElitistArchiveWithReplacementOfExistedMember(char *solution, double *solution_objective_values, double solution_constraint_value, char *is_new_nondominated_point, char *is_dominated_by_archive)
{
	short is_existed, index_of_existed_member;
	int   i, *indices_dominated, number_of_solutions_dominated;
	int distance_old, distance_new;

	*is_new_nondominated_point = TRUE;
	*is_dominated_by_archive = FALSE;

	if (elitist_archive_size == 0)
		addToElitistArchive(solution, solution_objective_values, solution_constraint_value);
	else
	{
		//indices_dominated = (int *)Malloc(elitist_archive_size * sizeof(int));
		indices_dominated = new int[elitist_archive_size];
		number_of_solutions_dominated = 0;
		is_existed = 0;
		for (i = 0; i < elitist_archive_size; i++)
		{
			if (constraintParetoDominates(elitist_archive_objective_values[i], elitist_archive_constraint_values[i], solution_objective_values, solution_constraint_value))
			{
				*is_dominated_by_archive = TRUE;
				*is_new_nondominated_point = FALSE;
			}
			else
			{
				if (!constraintParetoDominates(solution_objective_values, solution_constraint_value, elitist_archive_objective_values[i], elitist_archive_constraint_values[i]))
				{
					if (sameObjectiveBox(elitist_archive_objective_values[i], solution_objective_values))
					{
						is_existed = 1;
						index_of_existed_member = i;
						*is_new_nondominated_point = FALSE;
					}
				}
			}

			if ((*is_new_nondominated_point) == FALSE)
				break;
		}

		if ((*is_new_nondominated_point) == TRUE)
		{
			for (i = 0; i < elitist_archive_size; i++)
			{
				if (constraintParetoDominates(solution_objective_values, solution_constraint_value, elitist_archive_objective_values[i], elitist_archive_constraint_values[i]))
				{
					indices_dominated[number_of_solutions_dominated] = i;
					number_of_solutions_dominated++;
				}
			}

			if (number_of_solutions_dominated > 0)
				removeFromElitistArchive(indices_dominated, number_of_solutions_dominated);

			addToElitistArchive(solution, solution_objective_values, solution_constraint_value);
			elitist_archive_front_changed = TRUE;
		}

		if (is_existed)
		{
			distance_old = hammingDistanceToNearestNeighborInParameterSpace(elitist_archive[index_of_existed_member], index_of_existed_member);
			distance_new = hammingDistanceToNearestNeighborInParameterSpace(solution, index_of_existed_member);

			if (distance_new > distance_old)
			{
				for (i = 0; i < number_of_parameters; i++)
					elitist_archive[index_of_existed_member][i] = solution[i];
				for (i = 0; i < number_of_objectives; i++)
					elitist_archive_objective_values[index_of_existed_member][i] = solution_objective_values[i];
				elitist_archive_constraint_values[index_of_existed_member] = solution_constraint_value;
			}
		}

		free(indices_dominated);
	}
}//void CMoGomea::updateElitistArchiveWithReplacementOfExistedMember(char *solution, double *solution_objective_values, double solution_constraint_value, char *is_new_nondominated_point, char *is_dominated_by_archive)



short CMoGomea::constraintWeaklyParetoDominates(double *objective_values_x, double constraint_value_x, double *objective_values_y, double constraint_value_y)
{
	short result;

	result = FALSE;

	if (constraint_value_x > 0) /* x is infeasible */
	{
		if (constraint_value_y > 0) /* Both are infeasible */
		{
			if (constraint_value_x <= constraint_value_y)
				result = TRUE;
		}
	}
	else /* x is feasible */
	{
		if (constraint_value_y > 0) /* x is feasible and y is not */
			result = TRUE;
		else /* Both are feasible */
			result = weaklyParetoDominates(objective_values_x, objective_values_y);
	}

	return(result);
}//short CMoGomea::constraintWeaklyParetoDominates(double *objective_values_x, double constraint_value_x, double *objective_values_y, double constraint_value_y)



short CMoGomea::weaklyParetoDominates(double *objective_values_x, double *objective_values_y)
{
	int   i, result;
	result = 1;

	for (i = 0; i < number_of_objectives; i++)
	{
		if (fabs(objective_values_x[i] - objective_values_y[i]) >= 0.00001)
		{
			if (optimization[i] == MINIMIZATION)
			{
				if (objective_values_x[i] > objective_values_y[i])
				{
					result = 0;
					break;
				}
			}
			else if (optimization[i] == MAXIMIZATION)
			{
				if (objective_values_x[i] < objective_values_y[i])
				{
					result = 0;
					break;
				}
			}

		}
	}

	return(result);
}//short CMoGomea::weaklyParetoDominates(double *objective_values_x, double *objective_values_y)



/**
 * Returns 1 if x constraint-Pareto-dominates y, 0 otherwise.
 * x is not better than y unless:
 * - x and y are both infeasible and x has a smaller sum of constraint violations, or
 * - x is feasible and y is not, or
 * - x and y are both feasible and x Pareto dominates y
 */
short CMoGomea::constraintParetoDominates(double *objective_values_x, double constraint_value_x, double *objective_values_y, double constraint_value_y)
{
	short result;

	result = FALSE;

	if (constraint_value_x > 0) /* x is infeasible */
	{
		if (constraint_value_y > 0) /* Both are infeasible */
		{
			if (constraint_value_x < constraint_value_y)
				result = TRUE;
		}
	}
	else /* x is feasible */
	{
		if (constraint_value_y > 0) /* x is feasible and y is not */
			result = TRUE;
		else /* Both are feasible */
			result = paretoDominates(objective_values_x, objective_values_y);
	}

	return(result);
}//short CMoGomea::constraintParetoDominates(double *objective_values_x, double constraint_value_x, double *objective_values_y, double constraint_value_y)



/**
 * Returns 1 if x Pareto-dominates y, 0 otherwise.
 */
short CMoGomea::paretoDominates(double *objective_values_x, double *objective_values_y)
{
	short strict;
	int   i, result;

	result = 1;
	strict = 0;

	for (i = 0; i < number_of_objectives; i++)
	{
		if (fabs(objective_values_x[i] - objective_values_y[i]) >= 0.00001)
		{
			if (optimization[i] == MINIMIZATION)
			{
				if (objective_values_x[i] > objective_values_y[i])
				{
					result = 0;
					break;
				}
				if (objective_values_x[i] < objective_values_y[i])
					strict = 1;
			}
			else if (optimization[i] == MAXIMIZATION)
			{
				if (objective_values_x[i] < objective_values_y[i])
				{
					result = 0;
					break;
				}
				if (objective_values_x[i] > objective_values_y[i])
					strict = 1;
			}

		}
	}

	if (strict == 0 && result == 1)
		result = 0;

	return(result);
}//short CMoGomea::paretoDominates(double *objective_values_x, double *objective_values_y)





char CMoGomea::isDominatedByElitistArchive(double *obj, double con, char *is_new_nondominated_point, int *position_of_existed_member)
{
	int j;

	*is_new_nondominated_point = TRUE;
	*position_of_existed_member = -1;
	for (j = 0; j < elitist_archive_size; j++)
	{
		if (constraintParetoDominates(elitist_archive_objective_values[j], elitist_archive_constraint_values[j], obj, con))
		{
			*is_new_nondominated_point = FALSE;
			return(TRUE);
		}
		else
		{
			if (!constraintParetoDominates(obj, con, elitist_archive_objective_values[j], elitist_archive_constraint_values[j]))
			{
				if (sameObjectiveBox(elitist_archive_objective_values[j], obj))
				{
					*is_new_nondominated_point = FALSE;
					*position_of_existed_member = j;
					return(FALSE);
				}
			}
		}
	}
	return(FALSE);
}//char CMoGomea::isDominatedByElitistArchive(double *obj, double con, char *is_new_nondominated_point, int *position_of_existed_member)


/**
 * Returns 1 if two solutions share the same objective box, 0 otherwise.
 */
short CMoGomea::sameObjectiveBox(double *objective_values_a, double *objective_values_b)
{
	int i;

	if (!objective_discretization_in_effect)
	{
		/* If the solutions are identical, they are still in the (infinitely small) same objective box. */
		for (i = 0; i < number_of_objectives; i++)
		{
			if (objective_values_a[i] != objective_values_b[i])
				return(0);
		}

		return(1);
	}


	for (i = 0; i < number_of_objectives; i++)
	{
		if (((int)(objective_values_a[i] / objective_discretization[i])) != ((int)(objective_values_b[i] / objective_discretization[i])))
			return(0);
	}

	return(1);
}//short CMoGomea::sameObjectiveBox(double *objective_values_a, double *objective_values_b)



int CMoGomea::hammingDistanceInParameterSpace(char *solution_1, char *solution_2)
{
	int i, distance;
	distance = 0;
	for (i = 0; i < number_of_parameters; i++)
	{
		if (solution_1[i] != solution_2[i])
			distance++;
	}

	return distance;
}//int CMoGomea::hammingDistanceInParameterSpace(char *solution_1, char *solution_2)



int CMoGomea::hammingDistanceToNearestNeighborInParameterSpace(char *solution, int replacement_position)
{
	int i, distance_to_nearest_neighbor, distance;
	distance_to_nearest_neighbor = -1;
	for (i = 0; i < elitist_archive_size; i++)
	{
		if (i != replacement_position)
		{
			distance = hammingDistanceInParameterSpace(solution, elitist_archive[i]);
			if (distance < distance_to_nearest_neighbor || distance_to_nearest_neighbor < 0)
				distance_to_nearest_neighbor = distance;
		}
	}

	return distance_to_nearest_neighbor;
}//int CMoGomea::hammingDistanceToNearestNeighborInParameterSpace(char *solution, int replacement_position)





void CMoGomea::checkWhichSmallerPopulationsNeedToStop()
{
	int i;
	for (i = population_id - 1; i >= 0; i--)
	{
		if (array_of_population_statuses[i] == FALSE)
			continue;
		if (checkParetoFrontCover(population_id, i) == TRUE)
			array_of_population_statuses[i] = FALSE;
	}
}//void CMoGomea::checkWhichSmallerPopulationsNeedToStop()



void CMoGomea::updateParetoFrontForCurrentPopulation(double **objective_values_pop, double *constraint_values_pop, int pop_size)
{
	int i, j, index, rank0_size;
	char *isDominated;
	//isDominated = (char*)Malloc(pop_size * sizeof(char));
	isDominated = new char[pop_size * sizeof(char)];
	for (i = 0; i < pop_size; i++)
		isDominated[i] = FALSE;
	for (i = 0; i < pop_size; i++)
	{
		if (isDominated[i] == TRUE)
			continue;
		for (j = i + 1; j < pop_size; j++)
		{
			if (isDominated[j] == TRUE)
				continue;
			if (constraintParetoDominates(objective_values_pop[i], constraint_values_pop[i], objective_values_pop[j], constraint_values_pop[j]) == TRUE)
				isDominated[j] = TRUE;
			else if (constraintParetoDominates(objective_values_pop[j], constraint_values_pop[j], objective_values_pop[i], constraint_values_pop[i]) == TRUE)
			{
				isDominated[i] = TRUE;
				break;
			}
		}
	}

	rank0_size = 0;
	for (i = 0; i < pop_size; i++)
		if (isDominated[i] == FALSE)
			rank0_size++;

	if (array_of_Pareto_front_size_of_each_population[population_id] > 0)
	{
		for (i = 0; i < array_of_Pareto_front_size_of_each_population[population_id]; i++)
		{
			free(array_of_Pareto_front_of_each_population[population_id][i]);
		}
		free(array_of_Pareto_front_of_each_population[population_id]);
	}

	//array_of_Pareto_front_of_each_population[population_id] = (double**)Malloc(rank0_size * sizeof(double*));
	array_of_Pareto_front_of_each_population[population_id] = new double*[rank0_size];
	for (i = 0; i < rank0_size; i++)
		//array_of_Pareto_front_of_each_population[population_id][i] = (double*)Malloc(number_of_objectives * sizeof(double));
		array_of_Pareto_front_of_each_population[population_id][i] = new double[number_of_objectives];
	array_of_Pareto_front_size_of_each_population[population_id] = rank0_size;

	index = 0;
	for (i = 0; i < pop_size; i++)
	{
		if (isDominated[i] == TRUE)
			continue;
		for (j = 0; j < number_of_objectives; j++)
			array_of_Pareto_front_of_each_population[population_id][index][j] = objective_values_pop[i][j];
		index++;
	}
	free(isDominated);
}//void CMoGomea::updateParetoFrontForCurrentPopulation(double **objective_values_pop, double *constraint_values_pop, int pop_size)



char CMoGomea::checkParetoFrontCover(int pop_index_1, int pop_index_2)
{
	int i, j, count;
	count = 0;

	for (i = 0; i < array_of_Pareto_front_size_of_each_population[pop_index_2]; i++)
	{
		for (j = 0; j < array_of_Pareto_front_size_of_each_population[pop_index_1]; j++)
			if ((constraintParetoDominates(array_of_Pareto_front_of_each_population[pop_index_1][j], 0,
				array_of_Pareto_front_of_each_population[pop_index_2][i], 0) == TRUE) ||
				sameObjectiveBox(array_of_Pareto_front_of_each_population[pop_index_1][j], array_of_Pareto_front_of_each_population[pop_index_2][i]) == TRUE)
			{
				count++;
				break;
			}
	}
	// Check if all points in front 2 are dominated by or exist in front 1
	if (count == array_of_Pareto_front_size_of_each_population[pop_index_2])
		return TRUE;
	return FALSE;
}//char CMoGomea::checkParetoFrontCover(int pop_index_1, int pop_index_2)




void CMoGomea::putInitializedPopulationIntoArray()
{
	array_of_objective_ranges[population_id] = objective_ranges;
	array_of_populations[population_id] = population;
	array_of_objective_values[population_id] = objective_values;
	array_of_constraint_values[population_id] = constraint_values;
	array_of_t_NIS[population_id] = 0;
};//void CMoGomea::putInitializedPopulationIntoArray()








void CMoGomea::learnLinkageTree(int cluster_index)
{
	char   done;
	int    i, j, k, a, b, c, r0, r1, *indices, *order,
		lt_index, factor_size, **mpm_new, *mpm_new_number_of_indices, mpm_new_length,
		*NN_chain, NN_chain_length;
	double p, *cumulative_probabilities, **S_matrix, mul0, mul1;

	/* Compute joint entropy matrix */
	for (i = 0; i < number_of_parameters; i++)
	{
		for (j = i + 1; j < number_of_parameters; j++)
		{
			//indices = (int *)Malloc(2 * sizeof(int));
			indices = new int[2];
			indices[0] = i;
			indices[1] = j;
			cumulative_probabilities = estimateParametersForSingleBinaryMarginal(cluster_index, indices, 2, &factor_size);

			MI_matrix[i][j] = 0.0;
			for (k = 0; k < factor_size; k++)
			{
				if (k == 0)
					p = cumulative_probabilities[k];
				else
					p = cumulative_probabilities[k] - cumulative_probabilities[k - 1];
				if (p > 0)
					MI_matrix[i][j] += -p * log2(p);
			}

			MI_matrix[j][i] = MI_matrix[i][j];

			free(indices);
			free(cumulative_probabilities);
		}
		//indices = (int *)Malloc(1 * sizeof(int));
		indices = new int;
		indices[0] = i;
		cumulative_probabilities = estimateParametersForSingleBinaryMarginal(cluster_index, indices, 1, &factor_size);

		MI_matrix[i][i] = 0.0;
		for (k = 0; k < factor_size; k++)
		{
			if (k == 0)
				p = cumulative_probabilities[k];
			else
				p = cumulative_probabilities[k] - cumulative_probabilities[k - 1];
			if (p > 0)
				MI_matrix[i][i] += -p * log2(p);
		}

		free(indices);
		free(cumulative_probabilities);
	}

	/* Then transform into mutual information matrix MI(X,Y)=H(X)+H(Y)-H(X,Y) */
	for (i = 0; i < number_of_parameters; i++)
		for (j = i + 1; j < number_of_parameters; j++)
		{
			MI_matrix[i][j] = MI_matrix[i][i] + MI_matrix[j][j] - MI_matrix[i][j];
			MI_matrix[j][i] = MI_matrix[i][j];
		}


	/* Initialize MPM to the univariate factorization */
	order = createRandomOrdering(number_of_parameters);
	//mpm = (int **)Malloc(number_of_parameters * sizeof(int *));
	mpm = new int*[number_of_parameters];
	//mpm_number_of_indices = (int *)Malloc(number_of_parameters * sizeof(int));
	mpm_number_of_indices = new int[number_of_parameters];
	mpm_length = number_of_parameters;
	for (i = 0; i < number_of_parameters; i++)
	{
		//indices = (int *)Malloc(1 * sizeof(int));
		indices = new int;
		indices[0] = order[i];
		mpm[i] = indices;
		mpm_number_of_indices[i] = 1;
	}
	free(order);

	/* Initialize LT to the initial MPM */
	if (lt[cluster_index] != NULL)
	{
		for (i = 0; i < lt_length[cluster_index]; i++)
			free(lt[cluster_index][i]);
		free(lt[cluster_index]);
		free(lt_number_of_indices[cluster_index]);
	}
	//lt[cluster_index] = (int **)Malloc((number_of_parameters + number_of_parameters - 1) * sizeof(int *));
	lt[cluster_index] = new int*[(number_of_parameters + number_of_parameters - 1)];
	//lt_number_of_indices[cluster_index] = (int *)Malloc((number_of_parameters + number_of_parameters - 1) * sizeof(int));
	lt_number_of_indices[cluster_index] = new int[(number_of_parameters + number_of_parameters - 1)];
	lt_length[cluster_index] = number_of_parameters + number_of_parameters - 1;
	lt_index = 0;
	for (i = 0; i < mpm_length; i++)
	{
		lt[cluster_index][lt_index] = mpm[i];
		lt_number_of_indices[cluster_index][lt_index] = mpm_number_of_indices[i];
		lt_index++;
	}

	/* Initialize similarity matrix */
	//S_matrix = (double **)Malloc(number_of_parameters * sizeof(double *));
	S_matrix = new  double *[number_of_parameters];
	for (i = 0; i < number_of_parameters; i++)
		//S_matrix[i] = (double *)Malloc(number_of_parameters * sizeof(double));
		S_matrix[i] = new double[number_of_parameters];
	for (i = 0; i < mpm_length; i++)
		for (j = 0; j < mpm_length; j++)
			S_matrix[i][j] = MI_matrix[mpm[i][0]][mpm[j][0]];
	for (i = 0; i < mpm_length; i++)
		S_matrix[i][i] = 0;

	//NN_chain = (int *)Malloc((number_of_parameters + 2) * sizeof(int));
	NN_chain = new int[(number_of_parameters + 2)];
	NN_chain_length = 0;
	done = FALSE;
	while (done == FALSE)
	{
		if (NN_chain_length == 0)
		{
			NN_chain[NN_chain_length] = randomInt(mpm_length);
			NN_chain_length++;
		}

		while (NN_chain_length < 3)
		{
			NN_chain[NN_chain_length] = determineNearestNeighbour(NN_chain[NN_chain_length - 1], S_matrix, mpm_length);
			NN_chain_length++;
		}

		while (NN_chain[NN_chain_length - 3] != NN_chain[NN_chain_length - 1])
		{
			NN_chain[NN_chain_length] = determineNearestNeighbour(NN_chain[NN_chain_length - 1], S_matrix, mpm_length);
			if (((S_matrix[NN_chain[NN_chain_length - 1]][NN_chain[NN_chain_length]] == S_matrix[NN_chain[NN_chain_length - 1]][NN_chain[NN_chain_length - 2]])) && (NN_chain[NN_chain_length] != NN_chain[NN_chain_length - 2]))
				NN_chain[NN_chain_length] = NN_chain[NN_chain_length - 2];
			NN_chain_length++;
		}
		r0 = NN_chain[NN_chain_length - 2];
		r1 = NN_chain[NN_chain_length - 1];
		if (r0 > r1)
		{
			a = r0;
			r0 = r1;
			r1 = a;
		}
		NN_chain_length -= 3;

		if (r1 < mpm_length) // This test is required for exceptional cases in which the nearest-neighbor ordering has changed within the chain while merging within that chain
		{
			//indices = (int *)Malloc((mpm_number_of_indices[r0] + mpm_number_of_indices[r1]) * sizeof(int));
			indices = new int[(mpm_number_of_indices[r0] + mpm_number_of_indices[r1])];

			i = 0;
			for (j = 0; j < mpm_number_of_indices[r0]; j++)
			{
				indices[i] = mpm[r0][j];
				i++;
			}
			for (j = 0; j < mpm_number_of_indices[r1]; j++)
			{
				indices[i] = mpm[r1][j];
				i++;
			}

			lt[cluster_index][lt_index] = indices;
			lt_number_of_indices[cluster_index][lt_index] = mpm_number_of_indices[r0] + mpm_number_of_indices[r1];
			lt_index++;

			mul0 = ((double)mpm_number_of_indices[r0]) / ((double)mpm_number_of_indices[r0] + mpm_number_of_indices[r1]);
			mul1 = ((double)mpm_number_of_indices[r1]) / ((double)mpm_number_of_indices[r0] + mpm_number_of_indices[r1]);
			for (i = 0; i < mpm_length; i++)
			{
				if ((i != r0) && (i != r1))
				{
					S_matrix[i][r0] = mul0 * S_matrix[i][r0] + mul1 * S_matrix[i][r1];
					S_matrix[r0][i] = S_matrix[i][r0];
				}
			}

			//mpm_new = (int **)Malloc((mpm_length - 1) * sizeof(int *));
			mpm_new = new int*[(mpm_length - 1)];
			//mpm_new_number_of_indices = (int *)Malloc((mpm_length - 1) * sizeof(int));
			mpm_new_number_of_indices = new int[(mpm_length - 1)];
			mpm_new_length = mpm_length - 1;
			for (i = 0; i < mpm_new_length; i++)
			{
				mpm_new[i] = mpm[i];
				mpm_new_number_of_indices[i] = mpm_number_of_indices[i];
			}

			mpm_new[r0] = indices;
			mpm_new_number_of_indices[r0] = mpm_number_of_indices[r0] + mpm_number_of_indices[r1];
			if (r1 < mpm_length - 1)
			{
				mpm_new[r1] = mpm[mpm_length - 1];
				mpm_new_number_of_indices[r1] = mpm_number_of_indices[mpm_length - 1];

				for (i = 0; i < r1; i++)
				{
					S_matrix[i][r1] = S_matrix[i][mpm_length - 1];
					S_matrix[r1][i] = S_matrix[i][r1];
				}

				for (j = r1 + 1; j < mpm_new_length; j++)
				{
					S_matrix[r1][j] = S_matrix[j][mpm_length - 1];
					S_matrix[j][r1] = S_matrix[r1][j];
				}
			}

			for (i = 0; i < NN_chain_length; i++)
			{
				if (NN_chain[i] == mpm_length - 1)
				{
					NN_chain[i] = r1;
					break;
				}
			}

			free(mpm);
			free(mpm_number_of_indices);
			mpm = mpm_new;
			mpm_number_of_indices = mpm_new_number_of_indices;
			mpm_length = mpm_new_length;

			if (mpm_length == 1)
				done = TRUE;
		}
	}

	free(NN_chain);

	free(mpm_new);
	free(mpm_number_of_indices);

	for (i = 0; i < number_of_parameters; i++)
		free(S_matrix[i]);
	free(S_matrix);
}//void CMoGomea::learnLinkageTree(int cluster_index)



/**
 * Estimates the cumulative probability distribution of a
 * single binary marginal for a cluster (subpopulation).
 */
double *CMoGomea::estimateParametersForSingleBinaryMarginal(int cluster_index, int *indices, int number_of_indices, int *factor_size)
{
	int     i, j, index, power_of_two;
	char *solution;
	double *result;

	*factor_size = (int)pow(2, number_of_indices);
	//result = (double *)Malloc((*factor_size) * sizeof(double));
	result = new double[(*factor_size)];

	for (i = 0; i < (*factor_size); i++)
		result[i] = 0.0;

	for (i = 0; i < population_cluster_sizes[cluster_index]; i++)
	{
		index = 0;
		power_of_two = 1;
		for (j = number_of_indices - 1; j >= 0; j--)
		{
			solution = population[population_indices_of_cluster_members[cluster_index][i]];
			index += (solution[indices[j]] == TRUE) ? power_of_two : 0;
			power_of_two *= 2;
		}

		result[index] += 1.0;
	}

	for (i = 0; i < (*factor_size); i++)
		result[i] /= (double)population_cluster_sizes[cluster_index];

	for (i = 1; i < (*factor_size); i++)
		result[i] += result[i - 1];

	result[(*factor_size) - 1] = 1.0;

	return(result);
}//double *CMoGomea::estimateParametersForSingleBinaryMarginal(int cluster_index, int *indices, int number_of_indices, int *factor_size)



/**
 * Determines nearest neighbour according to similarity values.
 */
int CMoGomea::determineNearestNeighbour(int index, double **S_matrix, int mpm_length)
{
	int i, result;

	result = 0;
	if (result == index)
		result++;
	for (i = 1; i < mpm_length; i++)
	{
		//    if( (S_matrix[index][i] > S_matrix[index][result]) && (i != index) )
		if (((S_matrix[index][i] > S_matrix[index][result]) || ((S_matrix[index][i] == S_matrix[index][result]) && (mpm_number_of_indices[i] < mpm_number_of_indices[result]))) && (i != index))
			result = i;
	}

	return(result);
}//int CMoGomea::determineNearestNeighbour(int index, double **S_matrix, int mpm_length)




void CMoGomea::printLTStructure(int cluster_index)
{
	int i, j;

	for (i = 0; i < lt_length[cluster_index]; i++)
	{
		printf("[");
		for (j = 0; j < lt_number_of_indices[cluster_index][i]; j++)
		{
			printf("%d", lt[cluster_index][i][j]);
			if (j < lt_number_of_indices[cluster_index][i] - 1)
				printf(" ");
		}
		printf("]\n");
	}
	printf("\n");
	fflush(stdout);
};//void CMoGomea::printLTStructure(int cluster_index)





