#include "MultiObjectiveOptimizer.h"



//---------------------------------------------CBinaryMultiObjectiveOptimizer-------------------------------------------------------
CBinaryMultiObjectiveOptimizer::CBinaryMultiObjectiveOptimizer(CProblem<CBinaryCoding, CBinaryCoding>* pcProblem, CLog * pcLog, uint32_t iRandomSeed)
	: CBinaryOptimizer(pcProblem, pcLog, iRandomSeed)
{
	//pc_problem = (CBinaryMultiObjectiveProblem *) pcProblem;

}//CBinaryMultiObjectiveOptimizer::CBinaryMultiObjectiveOptimizer(CProblem<CBinaryCoding, CBinaryCoding>* pcProblem, CLog * pcLog, uint32_t iRandomSeed)


CBinaryMultiObjectiveOptimizer::CBinaryMultiObjectiveOptimizer(CBinaryMultiObjectiveOptimizer *pcOther) : CBinaryOptimizer(pcOther)
{
	::MessageBox(NULL, "Brak implementacji: CNSGA2::CNSGA2(CNSGA2 *pcOther) : CBinaryOptimizer(pcOther)", "BRAK", MB_OK);
}//CBinaryMultiObjectiveOptimizer::CBinaryMultiObjectiveOptimizer(CBinaryMultiObjectiveOptimizer *pcOther) : CBinaryOptimizer(pcOther)


CBinaryMultiObjectiveOptimizer::~CBinaryMultiObjectiveOptimizer()
{
	
}//CBinaryMultiObjectiveOptimizer::~CBinaryMultiObjectiveOptimizer()


CBinaryMultiObjectiveProblem  *CBinaryMultiObjectiveOptimizer::pc_get_multi_problem()
{
	if (pc_problem == NULL)  return(NULL);

	CBinaryMultiObjectiveProblem *pc_problem_multi;
	pc_problem_multi = (CBinaryMultiObjectiveProblem *)pc_problem->pcGetEvaluation();

	return(pc_problem_multi);
}//CBinaryMultiObjectiveProblem  *CBinaryMultiObjectiveOptimizer::pc_get_multi_problem()



double  CBinaryMultiObjectiveOptimizer::dPFQualityHyperVolume()
{
	CBinaryMultiObjectiveProblem *pc_problem_multi;
	pc_problem_multi = pc_get_multi_problem();

	if (pc_problem_multi == NULL)  return(-1);

	return(pc_problem_multi->dPFQualityHyperVolume());
}//double  CBinaryMultiObjectiveOptimizer::dPFQualityHyperVolume()



double  CBinaryMultiObjectiveOptimizer::dPFQualityInverseGenerationalDistance()
{
	CBinaryMultiObjectiveProblem *pc_problem_multi;
	pc_problem_multi = pc_get_multi_problem();

	if (pc_problem_multi == NULL)  return(-1);

	return(pc_problem_multi->dPFQualityInverseGenerationalDistance());
}//double  CBinaryMultiObjectiveOptimizer::dPFQualityInverseGenerationalDistance()



double  CBinaryMultiObjectiveOptimizer::dPFQualityGenerationalDistance()
{
	CBinaryMultiObjectiveProblem *pc_problem_multi;
	pc_problem_multi = pc_get_multi_problem();

	if (pc_problem_multi == NULL)  return(-1);

	return(pc_problem_multi->dPFQualityGenerationalDistance());
}//double  CBinaryMultiObjectiveOptimizer::dPFQualityGenerationalDistance()


double  CBinaryMultiObjectiveOptimizer::dPFQualityMaximumSpread()
{
	CBinaryMultiObjectiveProblem *pc_problem_multi;
	pc_problem_multi = pc_get_multi_problem();

	if (pc_problem_multi == NULL)  return(-1);

	return(pc_problem_multi->dPFQualityMaximumSpread());
}//double  CBinaryMultiObjectiveOptimizer::dPFQualityMaximumSpread()




double  CBinaryMultiObjectiveOptimizer::dPFQualityPointNum()
{
	CBinaryMultiObjectiveProblem *pc_problem_multi;
	pc_problem_multi = pc_get_multi_problem();

	if (pc_problem_multi == NULL)  return(-1);

	return(pc_problem_multi->dPFQualityPointNum());
}//double  CBinaryMultiObjectiveOptimizer::dPFQualityPointNum()



double  CBinaryMultiObjectiveOptimizer::dPFQualityDominatedOptimalPoints()
{
	CBinaryMultiObjectiveProblem *pc_problem_multi;
	pc_problem_multi = pc_get_multi_problem();

	if (pc_problem_multi == NULL)  return(-1);

	return(pc_problem_multi->dPFQualityDominatedOptimalPoints());
}//double  CBinaryMultiObjectiveOptimizer::dPFQualityDominatedOptimalPoints()



//---------------------------------------------CMultiIndividual-------------------------------------------------------

vector<double>  *CMultiIndividual::pvGetFitness()
{
	if (b_rated == true)  return(&v_fitness);
	vRate();
	return(&v_fitness);
};//vector<double>  *CMultiIndividual::pvGetFitness()



bool  CMultiIndividual::bFrontsDiffer(CMultiIndividual  *pcOther)
{
	vRate();

	for (int ii = 0; ii < v_fitness.size(); ii++)
	{
		if (v_fitness.at(ii) != pcOther->pvGetFitness()->at(ii))  return(false);
	}//for (int ii = 0; ii < v_fitness.size(); ii++)

	return(true);
};//bool  CMultiIndividual::bFrontsDiffer(CMultiIndividual  *pcOther)

