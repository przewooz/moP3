/////////////////////////////////////////////////////////////////////
//  Implementation of Multiobjective Evolutionary Algorithm Based on 
//  Decomposition (MOEA/D) For Multiobjective Knapsack Problem (2006)
//
//  See the details of MOEA/D in the following paper
//  Q. Zhang and H. Li, MOEA/D: A Multi-objective Evolutionary Algorithm Based on Decomposition, 
//  IEEE Trans. on Evolutionary Computation, in press, 2007
//
//  The source code of MOEA/D was implemented by Hui Li and Qingfu Zhang 
//
//  If you have any questions about the codes, please contact 
//  Qingfu Zhang at qzhang@essex.ac.uk  or Hui Li at hzl@cs.nott.ac.uk
/////////////////////////////////////////////////////////////////////


//*********************************************************************
//
//  (1) txt files to be loaded
//      The code reads the test instances from the file "TestInstances.txt". The data of these instances
//      are included in the folder "Instances". The predefined weight  vectors can be loaded from the 
//      files in the folder "Weights". 

//  (2) parameter settings
//      (a) The population size equals to the number of weight vectors. It is specified in the files containing
//      weight vectors.
//      (b) The neighborhood size can be set in "MOEAD.cpp".
//
//  (3) operators for knapsack
//      The implementation of crossover operator and greedy repair heuristics can be found in "Individual.cpp"
//
//*********************************************************************

//the original source codes were downloaded from https://github.com/ZhenkunWang/MOEADCODES
//2020.03.05 Michal Przewozniczek : The sources were modified to make them fully objective and the test case objects implemented in this framework were introduced into the original source code




#ifndef CMOEAD_OPTIMIZER_H
#define CMOEAD_OPTIMIZER_H


#define MOEAD_ARGUMENT_IPR  "IPr"



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




#define MOEAD_IM1 2147483563
#define MOEAD_IM2 2147483399
#define MOEAD_AM (1.0/IM1)
#define MOEAD_IMM1 (IM1-1)
#define MOEAD_IA1 40014
#define MOEAD_IA2 40692
#define MOEAD_IQ1 53668
#define MOEAD_IQ2 52774
#define MOEAD_IR1 12211
#define MOEAD_IR2 3791
#define MOEAD_NTAB 32
#define MOEAD_NDIV (1+IMM1/NTAB)
#define MOEAD_EPS 1.2e-7
#define MOEAD_RNMX (1.0-EPS)
#define MOEAD_SIGMA 1


using namespace std;

namespace MoeaD
{
	class CMoeadIndividual;
	class  CMoeadSubProblem;
	class CMoeadWeightVector;



	class  CMoead : public CBinaryMultiObjectiveOptimizer  //CBinaryOptimizer
	{
		static uint32_t iERROR_PARENT_CMoeaDOptimizer;
		static uint32_t iERROR_CODE_MOEAD_GENOTYPE_LEN_BELOW_0;

		friend class CMoeadIndividual;
		friend class CMoeadSubProblem;


	public:

		CMoead(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed);
		CMoead(CMoead *pcOther);
		~CMoead();

		CString  sAdditionalSummaryInfo();

		virtual COptimizer<CBinaryCoding, CBinaryCoding> *pcCopy() { return new CMoead(this); };

		virtual CError eConfigure(istream *psSettings);

		virtual void vInitialize(time_t tStartTime);
		virtual bool bRunIteration(uint32_t iIterationNumber, time_t tStartTime);
		bool bRunIteration_dummy(uint32_t iIterationNumber, time_t tStartTime);

	private:
		TimeCounters::CTimeCounter  c_time_counter;
		time_t t_start;

		int  i_templ_length;
		CBinaryCoding *pc_genotype;

		CBinaryMultiObjectiveProblem  *pc_multi_problem;


	public:
		vector<CMoeadSubProblem *> CurrentPopulation;  //  population of subproblems
		vector<CMoeadIndividual *> SecondPopulation;   //  external population containing all nondominated solutions found 

		double  d_mut_prob;
		int i_population_size;                     //  population size
		int NeighborhoodSize_Tr;                //  neighborhood size for replacement
		int NeighborhoodSize_Tm;                //  neighborhood size for mating

		int GlobalBestVectorIndex;              //  Global best vector for the new offspring
		double Pr;                              //  parameter control the replacement size

		vector<double> v_reference_point;

		void SaveSecondPopulation();
		void Show();
		void MinFastSort(vector<double> &dist, vector<int> &indx, int n, int m);

		void InitializePopulation();
		void InitializeWeightVector();
		void InitializeNeighborhood();
		void InitializeReferencePoint();

		void UpdateReferencePoint(CMoeadIndividual &ind);
		void UpdateSecondPopulation(CMoeadIndividual &ind);
		void UpdateNeighboringSolution(CMoeadIndividual &ind, int iPop);

		void Run(int fevals);

		void FindGlobalBestVectorIndex(CMoeadIndividual &ind);

	};//class  CMoead : public CBinaryMultiObjectiveOptimizer  //CBinaryOptimizer





	class CMoeadIndividual
	{
	public:
		CMoeadIndividual(CMoead  *pcParent);
		virtual ~CMoeadIndividual();


	private:
		CMoead  *pc_parent;

	public:
		vector<int>    v_genotype;     // items in knapsack - binary

		vector<double>  v_objectives;//measure values


	public:
		bool  dominated;

	public:
		void   SinglePointXover(CMoeadIndividual &parent1, CMoeadIndividual &parent2);
		double ComputingFitnessValue(vector<double> &lambda, int  iWeightingType);
		void   Randomize();

		void   Show();
		void   ComputeObjectives();          // function computing
		void   ComputingCapacityValue__();        // capacity computing

		//void   GreedyRepairHeuristic(vector<double> &lambda, char* strFuncType);
		void   GreedyRepairHeuristic(vector<double> &lambda, int  iWeightingType);

		bool   IsFeasible();           //check the constraints


		bool operator==(const CMoeadIndividual& ind);
		bool operator>(const  CMoeadIndividual& ind);
		bool operator<(const  CMoeadIndividual& ind);

	};//class CIndividual




	class CMoeadSubProblem
	{
	public:
		CMoeadSubProblem(CMoead  *pcParent);
		virtual ~CMoeadSubProblem();

	public:
		CMoeadIndividual    *pc_current_solution;
		CMoeadWeightVector  *pc_weight_vector;
		vector<int>    IndexOfNeighbor;

	private:
		CMoead  *pc_parent;
	};//class CMoeadSubProblem



	class CMoeadWeightVector
	{
	public:
		CMoeadWeightVector(CBinaryMultiObjectiveProblem  *pcMultiProblem);
		virtual ~CMoeadWeightVector();

		void Show();

	public:
		vector<double> lambda;

		CBinaryMultiObjectiveProblem  *pc_multi_problem;

	public:
		double DistanceTo(CMoeadWeightVector &weight);
	};

}//namespace MoeaD


#endif//CMOEAD_OPTIMIZER_H