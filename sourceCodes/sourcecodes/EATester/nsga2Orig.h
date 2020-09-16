/*
These source codes were partially downloaded from:
http://www.iitk.ac.in/kangal/codes.shtml
They were integrated with this project by Michal Przewozniczek 2019 mail:



The commentary in the main header file, of the original source codes is:
/* This is a Multi-Objective GA program.
**********************************************************************
*  This program is the implementation of the NSGA-2 proposed by      *
*                                                                    *
*  Prof. Kalyanmoy Deb and his students .                            *
*                                                                    *
*  copyright Kalyanmoy Deb
**********************************************************************

18.08.2003: The keepaliven.h file is modified to have normalized
			crowding distance calculation. The previous version of
			the code did not have this feature. This way, maintaining
			a good distribution of solutions in problems having quite
			a different range of objective functions were difficult.
			Hopefully, with this modification, such difficulties will
			not appear. --  K. Deb
18.08.2003: Also the dfit.h file is deleted. It was not needed any way.

The user have to give the input manualy or through a data file.

The user needs to enter objective functions in func-con.h
The code can also take care of the constraints. Enter the constraints
in the space provided in the func-con.h file.
Constraints must be of the following type:
g(x) >= 0.0
Also normalize all constraints (see the example problem in func-con.h)

If your program asks you to increase the values of some parameters in the
program come to main program and accordingly changed the values which are
defined against #define ...

The program generates few output files. These are described as
1.output.out
*           This file has the detailed record for all the variables,
*           the fitness values, constraint values, overall constraint
			violation (penalty)  and their ranks for all the members
*           of old population in the left hand side of the |**|
*           and of new population in the right hand side.

2.all_fitness.out
*         This file prints the record of all the fitness values for
*         different individual of new popultion created at all
*         generations.

3.g_rank_record.out
*        This file maintains the record of individuals in global pop-
*        -ulation at different ranks for all the generations.

4.ranks.out
*         This file prints the number of individual at different ranks
*          in old and new population and finds rank ratios

5.final_fitness.out
*                 This file has the fitness value of all feasible and
				  non-dominated individuals at the final generation

6.final_var.out
*                 This file has the all the variables of the feasible
				  and non-dominated individuals at the final generation.
				  The i-th solutions here corresponds to the i-th solution
				  in the final_fitness.out file.

7.plot.out        This file contains gnuplot-based file for plotting
				  the non-dominated feasible solutions obtained by the code.
*************************************************************************
*         This is recommended to delete or rename all the *.out files
*         obtained from the previous runs as some files are opened in
*         append mode so they give false resemblence of data if the
*         user is not careful

Compilation procedure:  gcc nsga2.c -lm
Run ./a.out with or without an input file

Input data files: Three files are included, but at one time one is needed
depending on the type of variables used:
inp-r (template file input-real)  : All variables are real-coded
inp-b (template file input-binary): All variables are binary-coded
inp-rb(template file input-rl+bin): Some variables are real and some are binary
*/


#ifndef CNSGA2_ORIG_OPTIMIZER_H
#define CNSGA2_ORIG_OPTIMIZER_H


#define NSGA2_ORIG_ARGUMENT_CROSS_TYPE "crossover_type (1-singlePoint;2-uniform)"
#define NSGA2_ORIG_ARGUMENT_PROB_CROSS "ProbCross"
#define NSGA2_ORIG_ARGUMENT_PROB_MUT "ProbMut"


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

#define square(x) ((x)*(x))
#define maxvar 1//#define maxvar    20  /*Max no. of variables*/


namespace Nsga2Orig
{
	class CNSGA2origIndividual : public CMultiIndividual
	{
	public:
		CNSGA2origIndividual();
		virtual ~CNSGA2origIndividual();

		void  vSetGenotype(int  iGenSize) { pc_genotype = new  CBinaryCoding(iGenSize); }
		CBinaryCoding* pcGetGenotype() { return(pc_genotype); }

		void  vSetProblem(CBinaryMultiObjectiveProblem  *pcMultiObjProblem);
		CMultiIndividual  *pcClone();

		double  dEvaluate();
		

		void  vRate();

		float  *fitness;

		int rank,              /*Rank of the individual*/
			flag;              /*Flag for ranking*/

		float	cub_len,             /*crowding distance of the individual*/
			error;              /* overall constraint violation for the individual*/
	};//class CNSGA2origIndividual



	class CNSGA2origPopulation
	{
	public:
		CNSGA2origPopulation(int  iPopSize);
		~CNSGA2origPopulation();

		int maxrank;            /*Maximum rank present in the population*/
		float *rankrat;  /*Rank Ratio*/
		int *rankno;     /*Individual at different ranks*/
		CNSGA2origIndividual *ind, /*Different Individuals*/
			*ind_ptr;

		int  i_pop_size;

	};//class CNSGA2origPopulation



	class  globpop
	{
	public:
		globpop(int  iPopSize,  int  iFuncNum,  int  iProblemLength);
		~globpop();

		int maxrank,   /*Max rank of the global population*/
			//rankar[2 * maxpop][2 * maxpop], /*record of array of individual numbers at a particular rank */
			**rankar,

			//rankno[2 * maxpop];           /*record of no. of individuals at a particular rank*/
			*rankno;

		//int genes[2 * maxpop][maxchrom],
		//	rank[2 * maxpop],            /*rank of different individuals*/
		//	flag[2 * maxpop];            /*Setting the flag */
		int  **genes;
		int  *rank;
		int  *flag;

		//float fitness[2 * maxpop][maxfun], /*Fitness function values for the different
		//				   individuals*/
		//	cub_len[2 * maxpop],              /*Dummyfitness*/
		//	xreal[2 * maxpop][maxvar],       /*value of the decoded variables for different individuals */
		//	xbin[2 * maxpop][maxvar],   /* binray-coded variables */
		//	error[2 * maxpop],               /*Error Values of the individuals*/
		//	constr[2 * maxpop][maxcons];
		float **fitness; 
		float *cub_len;

		int  i_pop_size;
		int  i_func_num;
		int  i_problem_length;
	};



	class CNSGA2orig : public CBinaryMultiObjectiveOptimizer  //CBinaryOptimizer
	{
		friend class CNSGA2Individual;
	public:
		static uint32_t iERROR_PARENT_CNSGA2origOptimizer;
		static uint32_t iERROR_CODE_NSGA2_ORIG_GENOTYPE_LEN_BELOW_0;

		CNSGA2orig(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed);
		CNSGA2orig(CNSGA2orig *pcOther);
		~CNSGA2orig();

		virtual COptimizer<CBinaryCoding, CBinaryCoding> *pcCopy() { return new CNSGA2orig(this); };

		virtual CError eConfigure(istream *psSettings);

		virtual void vInitialize(time_t tStartTime);
		bool bRunIteration_dummy(uint32_t iIterationNumber, time_t tStartTime);
		virtual bool bRunIteration(uint32_t iIterationNumber, time_t tStartTime);

		CString  sAdditionalSummaryInfo();

	private:
		void init(CNSGA2origPopulation *pop_ptr);
		void nselect(CNSGA2origPopulation *old_pop_ptr, CNSGA2origPopulation *pop2_ptr);
		float randomperc();
		void func(CNSGA2origPopulation *pop_ptr);
		void ranking(CNSGA2origPopulation *pop_ptr);
		void rankcon(CNSGA2origPopulation *pop_ptr);
		int indcmp(float *ptr1, float *ptr2);
		int indcmp3(float *ptr1, float *ptr2);
		int indcmp1(float *ptr1, float *ptr2);
		void unicross(CNSGA2origPopulation *new_pop_ptr, CNSGA2origPopulation *mate_pop_ptr);
		void crossover(CNSGA2origPopulation *new_pop_ptr, CNSGA2origPopulation *mate_pop_ptr);
		void mutate(CNSGA2origPopulation *new_pop_ptr);
		void mutate_bitwise(CNSGA2origPopulation *new_pop_ptr);
		void keepalive(CNSGA2origPopulation *pop1_ptr, CNSGA2origPopulation *pop2_ptr, CNSGA2origPopulation *pop3_ptr, int gen);
		void grank(int gen);
		void gshare(int rnk);
		void sort(int m1);
		void gsort(int rnk, int sel);

		CNSGA2origPopulation *old_pop_ptr;
		CNSGA2origPopulation  *new_pop_ptr;
		CNSGA2origPopulation  *mate_pop_ptr;

		globpop  *globalpop;


		TimeCounters::CTimeCounter  c_time_counter;
		time_t t_start;

		int  i_templ_length;//#define maxchrom 200  /*Max chromosome length*/


		int  i_pop_size;//#define maxpop   500  /*Max population */

		
//#define maxfun    10  /*Max no. of functions */
//#define maxcons   20  /*Max no. of Constraints*/

		int gener,       /*No of generations*/
			//nvar, nchrom,          /*No of variables*/
			ncons,         /*No of Constraints*/
			vlen[maxvar],  /*Array to store no of bits for each variable*/
			nmut,          /* No of Mutations */
			ncross,        /*No of crossovers*/
			ans;
		float seed,      /*Random Seed*/
			pcross,        /*Cross-over Probability*/
			pmut_b, pmut_r,          /*Mutation Probability*/
			lim_b[maxvar][2], lim_r[maxvar][2];/*Limits of variable in array*/
		float di,        /*Distribution Index for the Cross-over*/
			dim,           /*Distribution Index for the Mutation*/
			delta_fit,     /* variables required forfitness for fitness sharing */
			min_fit,
			front_ratio;
		int optype,      /*Cross-over type*/
			nfunc,         /*No of functions*/
			sharespace;    /*Sharing space (either parameter or fitness)*/

		double coef[maxvar]; /*Variable used for decoding*/

		//float fpara1[2 * maxpop][2];
		//keepalive
		float  **fpara1;
		int left, Lastrank;




		//float  array[2 * maxpop][2]
		float  **array;//void CNSGA2orig::gsort(int rnk, int sel)

		//float length[2 * maxpop][2] 
		float  **length;//void CNSGA2orig::gshare(int rnk)

		//int   gflg[2 * maxpop] 
		int   *gflg;//void CNSGA2orig::grank(int gen)
		int  *front_pop;//keepalive

		
		//static int popsize,  /*Population Size*/
			//chrom;             /*Chromosome size*/

	};//class CNSGA2 : public CBinaryOptimizer
}//namespace Nsga2Orig

#endif//CNSGA2_ORIG_OPTIMIZER_H