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



#ifndef CMOGOMEA_OPTIMIZER_H
#define CMOGOMEA_OPTIMIZER_H


#include "MultiObjectiveOptimizer.h"
#include "BinaryCoding.h"
#include "BinaryOptimizer.h"
#include "BinaryEvaluationMultiObjective.h"
#include "Error.h"
#include "Log.h"
#include "Optimizer.h"
#include "Problem.h"
#include  "util\timer.h"
#include  "util\tools.h"
#include "RandUtils.h"


#include <istream>
#include <algorithm>



namespace MoGomea
{
	//class CNSGA2Individual;

	class CMoGomea : public CBinaryMultiObjectiveOptimizer  //CBinaryOptimizer
	{
		//friend class CNSGA2Individual;
	public:
		static uint32_t iERROR_PARENT_CMoGomeaOptimizer;
		static uint32_t iERROR_CODE_MOGOMEA_GENOTYPE_LEN_BELOW_0;

		CMoGomea(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed);
		CMoGomea(CMoGomea *pcOther);
		~CMoGomea();

		virtual COptimizer<CBinaryCoding, CBinaryCoding> *pcCopy() { return new CMoGomea(this); };

		virtual CError eConfigure(istream *psSettings);

		virtual void vInitialize(time_t tStartTime);
		virtual bool bRunIteration(uint32_t iIterationNumber, time_t tStartTime);
		bool bRunIteration_dummy(uint32_t iIterationNumber, time_t tStartTime);

		CString  sAdditionalSummaryInfo();



		//mo gomea functions
		void v_initializeCommonVariables();

		void schedule();
		void schedule_runMultiplePop_clusterPop_learnPop_improvePop();
		
	private:
		void initializeMemoryForArrayOfPopulations();
		void initializeArrayOfParetoFronts();

		//at the end:
		void logNumberOfEvaluationsAtVTR();
		short haveDPFSMetric(void);
		void writeCurrentElitistArchive(char final);
		void ezilaitiniArrayOfPopulation();

		void assignPointersToCorrespondingPopulation();
		void ezilaitiniMemoryOfCorrespondingPopulation();
		void ezilaitiniArrayOfParetoFronts();

		void writeGenerationalStatistics(void);
		void adaptObjectiveDiscretization(void);

		void putInitializedPopulationIntoArray();
		void updateParetoFrontForCurrentPopulation(double **objective_values_pop, double *constraint_values_pop, int pop_size);
		void checkWhichSmallerPopulationsNeedToStop();
		char checkParetoFrontCover(int pop_index_1, int pop_index_2);


		//mogomea
		void initialize();
		void initializeMemory(void);
		void initializePopulationAndFitnessValues();
		void computeObjectiveRanges(void);
		void learnLinkageOnCurrentPopulation();
		int** clustering(double **objective_values_pool, int pool_size, int number_of_dimensions, int number_of_clusters, int *pool_cluster_size);
		int *greedyScatteredSubsetSelection(double **points, int number_of_points, int number_of_dimensions, int number_to_select);
		void determineExtremeClusters();
		void initializeClusters();
		void ezilaitiniClusters();
		void improveCurrentPopulation(void);
		void copyValuesFromDonorToOffspring(char *solution, char *donor, int cluster_index, int linkage_group_index);
		void copyFromAToB(char *solution_a, double *obj_a, double con_a, char *solution_b, double *obj_b, double *con_b);
		void mutateSolution(char *solution, int lt_factor_index, int cluster_index);
		void performMultiObjectiveGenepoolOptimalMixing(int cluster_index, char *parent, double *parent_obj, double parent_con, char *result, double *obj, double *con);
		void performSingleObjectiveGenepoolOptimalMixing(int cluster_index, int objective_index,
			char *parent, double *parent_obj, double parent_con,
			char *result, double *obj, double *con);
		void selectFinalSurvivors();
		void freeAuxiliaryPopulations();


		//utils
		int randomInt(int maximum);
		double randomRealUniform01(void);
		int* createRandomOrdering(int size_of_the_set);
		double distanceEuclidean(double *x, double *y, int number_of_dimensions);
		int *mergeSort(double *array, int array_size);
		void mergeSortWithinBounds(double *array, int *sorted, int *tosort, int p, int q);
		void mergeSortMerge(double *array, int *sorted, int *tosort, int p, int r, int q);


		void evaluateIndividual(char *solution, double *obj, double *con, int objective_index_of_extreme_cluster);
		char betterFitness(double *objective_value_x, double constraint_value_x, double *objective_value_y, double constraint_value_y, int objective_index);
		char equalFitness(double *objective_value_x, double constraint_value_x, double *objective_value_y, double constraint_value_y, int objective_index);
		void logElitistArchiveAtSpecificPoints();
		void updateElitistArchive(char *solution, double *solution_objective_values, double solution_constraint_value);
		void addToElitistArchive(char *solution, double *solution_objective_values, double solution_constraint_value);
		void removeFromElitistArchive(int *indices, int number_of_indices);
		void updateElitistArchiveWithReplacementOfExistedMember(char *solution, double *solution_objective_values, double solution_constraint_value, char *is_new_nondominated_point, char *is_dominated_by_archive);
		short constraintWeaklyParetoDominates(double *objective_values_x, double constraint_value_x, double *objective_values_y, double constraint_value_y);
		short weaklyParetoDominates(double *objective_values_x, double *objective_values_y);
		short constraintParetoDominates(double *objective_values_x, double constraint_value_x, double *objective_values_y, double constraint_value_y);
		short paretoDominates(double *objective_values_x, double *objective_values_y);

		char isDominatedByElitistArchive(double *obj, double con, char *is_new_nondominated_point, int *position_of_existed_member);
		short sameObjectiveBox(double *objective_values_a, double *objective_values_b);
		short isInListOfIndices(int index, int *indices, int number_of_indices);
		int hammingDistanceInParameterSpace(char *solution_1, char *solution_2);
		int hammingDistanceToNearestNeighborInParameterSpace(char *solution, int replacement_position);
		


		//linkage tree learning
		void learnLinkageTree(int cluster_index);
		double *estimateParametersForSingleBinaryMarginal(int cluster_index, int *indices, int number_of_indices, int *factor_size);
		int determineNearestNeighbour(int index, double **S_matrix, int mpm_length);
		void printLTStructure(int cluster_index);




		TimeCounters::CTimeCounter  c_time_counter;
		time_t t_start;

		int  i_templ_length;
		CBinaryCoding *pc_genotype;


		//mo gomea variables

		/*-=-=-=-=-=-=-=-=-=-=-=- Section Global Variables -=-=-=-=-=-=-=-=-=-=-=-=-*/
		char    **population,               /* The population containing the solutions. */
			***array_of_populations,    /* The array containing all populations in the parameter-less scheme. */
			**offspring,                /* Offspring solutions. */
			**elitist_archive,          /* Archive of elitist solutions. */
			**elitist_archive_copy;     /* Copy of the elitist archive. */

		int     problem_index,                          /* The index of the optimization problem. */
			number_of_parameters,                   /* The number of parameters to be optimized. */
			number_of_generations,                  /* The current generation count. */
			*array_of_number_of_generations,        /* The array containing generation counts of all populations in the parameter-less scheme.*/
			generation_base,                        /* The number of iterations that population of size N_i is run before running 1 iteration of population of size N_(i+1). */
			population_size,                        /* The size of each population. */
			*array_of_population_sizes,             /* The array containing population sizes of all populations in the parameter-less scheme. */
			smallest_population_size,               /* The size of the first population. */
			population_id,                          /* The index of the population that is currently operating. */
			offspring_size,                         /* The size of the offspring population. */

			number_of_objectives,                   /* The number of objective functions. */
			elitist_archive_size,                   /* Number of solutions in the elitist archive. */
			elitist_archive_size_target,            /* The lower bound of the targeted size of the elitist archive. */
			elitist_archive_copy_size,              /* Number of solutions in the elitist archive copy. */
			elitist_archive_capacity,               /* Current memory allocation to elitist archive. */
			number_of_mixing_components,            /* The number of components in the mixture distribution. */
			*array_of_number_of_clusters,           /* The array containing the number-of-cluster of each population in the parameter-less scheme. */
			*population_cluster_sizes,              /* The size of each cluster. */
			**population_indices_of_cluster_members,/* The corresponding index in the population of each solution in a cluster. */
			*which_extreme,                         /* The corresponding objective of an extreme-region cluster. */

			t_NIS,                          /* The number of subsequent generations without an improvement (no-improvement-stretch). */
			*array_of_t_NIS,                /* The array containing the no-improvement-stretch of each population in the parameter-less scheme. */
			maximum_number_of_populations,  /* The maximum number of populations that can be run (depending on memory budget). */
			number_of_populations,          /* The number of populations that have been initialized. */

			**mpm,                          /* The marginal product model. */
			*mpm_number_of_indices,         /* The number of variables in each factor in the mpm. */
			mpm_length,                     /* The number of factors in the mpm. */

			***lt,                          /* The linkage tree, one for each cluster. */
			**lt_number_of_indices,         /* The number of variables in each factor in the linkage tree of each cluster. */
			*lt_length;                     /* The number of factors in the linkage tree of each cluster. */

		long    number_of_evaluations,            /* The current number of times a function evaluation was performed. */
			log_progress_interval,            /* The interval (in terms of number of evaluations) at which the elitist archive is logged. */
			maximum_number_of_evaluations,    /* The maximum number of evaluations. */
			*array_of_number_of_evaluations_per_population; /* The array containing the number of evaluations used by each population in the parameter-less scheme. */

		double  **objective_values,                 /* Objective values for population members. */
			***array_of_objective_values,       /* The array containing objective values of all populations in the parameter-less scheme. */
			*constraint_values,                 /* Constraint values of population members. */
			**array_of_constraint_values,       /* The array containing constraint values of all populations in the parameter-less scheme. */

			**objective_values_offspring,       /* Objective values of offspring solutions. */
			*constraint_values_offspring,       /* Constraint values of offspring solutions. */

			**elitist_archive_objective_values,         /* Objective values of solutions stored in elitist archive. */
			**elitist_archive_copy_objective_values,    /* Copy of objective values of solutions stored in elitist archive. */
			*elitist_archive_constraint_values,         /* Constraint values of solutions stored in elitist archive. */
			*elitist_archive_copy_constraint_values,    /* Copy of constraint values of solutions stored in elitist archive. */

			*objective_ranges,                          /* Ranges of objectives observed in the current population. */
			**array_of_objective_ranges,                /* The array containing ranges of objectives observed in each population in the parameter-less scheme. */
			**objective_means_scaled,                   /* The means of the clusters in the objective space, linearly scaled according to the observed ranges. */
			*objective_discretization,                  /* The length of the objective discretization in each dimension (for the elitist archive). */
			vtr,                              /* The value-to-reach (in terms of Inverse Generational Distance). */
			**MI_matrix;                      /* Mutual information between any two variables */

		int64_t random_seed,                      /* The seed used for the random-number generator. */
			random_seed_changing;             /* Internally used variable for randomly setting a random seed. */

		char    use_pre_mutation,                   /* Whether to use weak mutation. */
			use_pre_adaptive_mutation,          /* Whether to use strong mutation. */
			use_print_progress_to_screen,       /* Whether to print the progress of the optimization to screen. */
			use_repair_mechanism,               /* Whether to use a repair mechanism (provided by users) if the problem is constrained. */
			*optimization,                      /* Maximization or Minimization for each objective. */
			print_verbose_overview,             /* Whether to print a overview of settings (0 = no). */
			use_vtr,                            /* Whether to terminate at the value-to-reach (VTR) (0 = no). */
			objective_discretization_in_effect, /* Whether the objective space is currently being discretized for the elitist archive. */
			elitist_archive_front_changed;      /* Whether the Pareto front formed by the elitist archive is changed in this generation. */


		char    *array_of_population_statuses;
		double  ***array_of_Pareto_front_of_each_population;
		int     *array_of_Pareto_front_size_of_each_population;
		char    stop_population_when_front_is_covered;
		

	};//class CNSGA2 : public CBinaryOptimizer




}//namespace Nsga2



#endif//CMOGOMEA_OPTIMIZER_H