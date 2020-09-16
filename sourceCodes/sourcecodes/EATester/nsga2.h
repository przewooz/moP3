#ifndef CNSGA2_OPTIMIZER_H
#define CNSGA2_OPTIMIZER_H

#define NSGA2_ARGUMENT_TOUR_SIZE "tournament_size"
#define NSGA2_ARGUMENT_PROB_CROSS "ProbCross"
#define NSGA2_ARGUMENT_PROB_MUT "ProbMut"
#define NSGA2_FULL_GENOTYPE_MUT "FullGenotypeMut"

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



namespace Nsga2
{
	class CNSGA2Individual;

	class CNSGA2 : public CBinaryMultiObjectiveOptimizer  //CBinaryOptimizer
	{
		friend class CNSGA2Individual;
	public:
		static uint32_t iERROR_PARENT_CNSGA2Optimizer;
		static uint32_t iERROR_CODE_NSGA2_GENOTYPE_LEN_BELOW_0;

		CNSGA2(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed);
		CNSGA2(CNSGA2 *pcOther);
		~CNSGA2();

		virtual COptimizer<CBinaryCoding, CBinaryCoding> *pcCopy() { return new CNSGA2(this); };

		virtual CError eConfigure(istream *psSettings);

		virtual void vInitialize(time_t tStartTime);
		virtual bool bRunIteration(uint32_t iIterationNumber, time_t tStartTime);

		CString  sAdditionalSummaryInfo();

	private:
		void  v_evolution();
		void  v_non_dominated_sorting();
		void  v_non_dominated_sorting_old();
		void  v_compute_crowding_distance();
		void  v_add_to_non_dominated_pf(vector<CMultiIndividual *>  *pvFrontToFill, CNSGA2Individual  *pcInd);
		void  v_report_pareto_front(int  iId);
		void  v_get_half_ofChildren_parent_pop();
		//void  v_update_global_pareto_front();
		//int   i_join_the_global_pareto(CNSGA2Individual *pcInd);//returns -1 if ind was not added, or the number of thrown-out individuals

		CNSGA2Individual *pc_parent_tournament(int iTournamentSize);



		TimeCounters::CTimeCounter  c_time_counter;
		time_t t_start;

		int  i_templ_length;

		
		int  i_pop_size;// = 1000;
		int  i_tournament_size;// = 2;
		double  d_prob_cross;// = 0.8;
		double  d_prob_mut;// = 0.05;

		bool  b_2_point_cross = false;
		bool  b_whole_genotype_mut;// = true;


		//vector<CNSGA2Individual  *>  v_global_pareto_front;
		vector<vector<CNSGA2Individual  *>>  v_non_dominated_pfs;

		vector<CNSGA2Individual  *>  v_population;

	};//class CNSGA2 : public CBinaryOptimizer



	class  CNSGA2Individual : public CMultiIndividual
	{
		friend class CNSGA2;
	public:
		CNSGA2Individual(int  iTemplLength, CBinaryMultiObjectiveProblem *pcProblem, CNSGA2 *pcParent);
		CNSGA2Individual(CNSGA2Individual &pcOther);
		~CNSGA2Individual();

		void  vInitialize();
		//virtual void  vRate();
		double  dEvaluate();

		virtual CMultiIndividual  *pcClone() { return(new CNSGA2Individual(*this)); };

		//vector<double>  *pvGetFitness();
		//int  iCheckDomination(CNSGA2Individual  *pcOther);
		void  vCompCrowdingDistance(CNSGA2Individual  *pcPrev, CNSGA2Individual  *pcNext);
		double  dGetCrowdingDistance() { return(d_crowding); };

		void  vCross(CNSGA2Individual  *pcOtherParent, vector<CNSGA2Individual *> *pvCrossResult);
		void  vCross2point(CNSGA2Individual  *pcOtherParent, vector<CNSGA2Individual *> *pvCrossResult);
		void  vMutate();
		void  vMutateWholeGenotype();
		
	private:
		CNSGA2  *pc_parent;
		//CBinaryMultiObjectiveProblem *pc_problem;

		int  i_templ_length;
		//CBinaryCoding *pc_genotype;

		//vector<double>  v_fitness;
		double  d_crowding;
		int  i_front;

		//bool  b_rated;
	};//class  CNSGA2Individual

}//namespace Nsga2



#endif//CNSGA2_OPTIMIZER_H