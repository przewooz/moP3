#ifndef MULTI_OBJECTIVE_OPTIMIZER
#define MULTI_OBJECTIVE_OPTIMIZER

#include "BinaryOptimizer.h"
#include "BinaryEvaluationMultiObjective.h"
#include "BinaryCoding.h"
#include "Individual.h"
#include "Log.h"
#include "Optimizer.h"
#include "Problem.h"

#include <ctime>
#include <cstdint>

using namespace std;


//class  CMultiIndividual;

class CBinaryMultiObjectiveOptimizer : public CBinaryOptimizer
{
public:
	CBinaryMultiObjectiveOptimizer(CProblem<CBinaryCoding, CBinaryCoding> *pcProblem, CLog *pcLog, uint32_t iRandomSeed);
	CBinaryMultiObjectiveOptimizer(CBinaryMultiObjectiveOptimizer *pcOther);


	double  dPFQualityHyperVolume();
	double  dPFQualityPointNum();
	double  dPFQualityInverseGenerationalDistance();
	double  dPFQualityGenerationalDistance();
	double  dPFQualityDominatedOptimalPoints();
	double  dPFQualityMaximumSpread();

	virtual ~CBinaryMultiObjectiveOptimizer();

protected:
	CBinaryMultiObjectiveProblem  *pc_get_multi_problem();

};//class class CBinaryMultiObjectiveOptimizer : public CBinaryOptimizer


/*class  CMultiIndividual
{
public:

	virtual void  vRate() {};
	vector<double>  *pvGetFitness();
	bool  bFrontsDiffer(CMultiIndividual  *pcOther);

protected:
	vector<double>  v_fitness;
	bool  b_rated;

};//class  CMultiIndividual*/


#endif//MULTI_OBJECTIVE_OPTIMIZER