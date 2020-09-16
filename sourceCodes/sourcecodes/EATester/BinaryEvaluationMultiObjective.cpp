#include "BinaryEvaluationMultiObjective.h"

#include "FloatCommandParam.h"
#include "StringCommandParam.h"
#include "UIntCommandParam.h"

#include "EvaluationUtils.h"

#include <atlstr.h>
#include <cfloat>
#include <functional>



CBinaryMultiObjectiveProblem::CBinaryMultiObjectiveProblem(const CBinaryMultiObjectiveProblem &pcOther)
{
	for (int i_measure = 0; i_measure < pcOther.v_measures.size(); i_measure++)
	{
		v_measures.push_back(pcOther.v_measures.at(i_measure)->pcClone());
	}//for (int i_measure = 0; i_measure < v_mesaures.size(); i_measure++)

	i_weghting_type = pcOther.i_weghting_type;
}//CBinaryMultiObjectiveProblem::CBinaryMultiObjectiveProblem(CBinaryMultiObjectiveProblem &pcOther)


CBinaryMultiObjectiveProblem::~CBinaryMultiObjectiveProblem()
{
	for (int i_measure = 0; i_measure < v_measures.size(); i_measure++)
		delete v_measures.at(i_measure);

	for (int ii = 0; ii < v_global_pareto_front.size(); ii++)
		delete  v_global_pareto_front.at(ii);

	for (int ii = 0; ii < v_theoretical_optimal_pareto_front.size(); ii++)
		delete  v_theoretical_optimal_pareto_front.at(ii);

}//CBinaryMultiObjectiveProblem::~CBinaryMultiObjectiveProblem()



void  CBinaryMultiObjectiveProblem::vReportPF(vector<CString> *pvPFReport)
{
	CString  s_line;

	s_line.Format("%d", v_global_pareto_front.size());
	pvPFReport->push_back(s_line);


	for (int ii = 0; ii < v_global_pareto_front.size(); ii++)
	{
		s_line = v_global_pareto_front.at(ii)->sReportObj();
		pvPFReport->push_back(s_line);
	}//for (int ii = 0; ii < v_global_pareto_front.size(); ii++)

}//void  CBinaryMultiObjectiveProblem::vReportPF(vector<CString> *pvPFReport)



CString  CBinaryMultiObjectiveProblem::sMultiObjectiveReportIter() 
{
	CString  s_result;
	
	s_result.Format
		(
			"OptPFfound: %.2lf (%d/%d) PFsize:%d HyperVolume: %.8lf InvGenDist: %.8lf GenDist:%.8lf MaxSpread:%.8lf DominatedPF: %d [ffe: %.0lf]", 
			dPFQualityOptimalPointsPerc(), iPFQualityOptimalPointsFound(), iPFQualityOptimalPointsNum(),
			(int)dPFQualityPointNum(), dPFQualityHyperVolume(), dPFQualityInverseGenerationalDistance(), dPFQualityGenerationalDistance(), dPFQualityMaximumSpread(), (int)dPFQualityDominatedOptimalPoints(), (double)iGetFFE()
		);

	return(s_result);
};//CString  CBinaryMultiObjectiveProblem::sMultiObjectiveReportIter() 



void  CBinaryMultiObjectiveProblem::vFlushGlobalPareto()
{
	for (int ii = 0; ii < v_global_pareto_front.size(); ii++)
		delete  v_global_pareto_front.at(ii);
	v_global_pareto_front.clear();

	c_time_counter.vSetStartNow();
}//void  CBinaryMultiObjectiveProblem::vFlushGlobalPareto()



void  CBinaryMultiObjectiveProblem::vEvaluateParetoFront(vector<double> *pvPF, CBinaryCoding *pcFenotype)
{
	v_evaluate_pareto_front(pvPF, pcFenotype);
}//void  CBinaryMultiObjectiveProblem::vEvaluateParetoFront(vector<double> *pvPF, CBinaryCoding *pcFenotype)



void  CBinaryMultiObjectiveProblem::v_evaluate_pareto_front(vector<double> *pvPF, CBinaryCoding *pcFenotype)
{
	double  d_time;
	if (c_time_counter.bGetTimePassed(&d_time) == false)  c_time_counter.vSetStartNow();

	i_ffe++;

	if (pvPF->size() != v_measures.size())
	{
		pvPF->clear();
		for (int ii = 0; ii < v_measures.size(); ii++)
			pvPF->push_back(0);
	}//if (pvPF->size() != v_mesaures.size())



	for (int i_measure = 0; i_measure < v_measures.size(); i_measure++)
	{
		pvPF->at(i_measure) = v_measures.at(i_measure)->dEvaluate(pcFenotype);
	}//for (int i_measure = 0; i_measure < v_mesaures.size(); i_measure++)


	CMultiIndividualPFpoints  c_dummy_pf_point;//the class and object are used only to get to CMultiIndividual methods
	c_dummy_pf_point.vConfigure(pcFenotype, false, pvPF, true, this);
	i_join_the_global_pareto(&c_dummy_pf_point);
}//void  CBinaryMultiObjectiveProblem::vEvaluateParetoFront(vector<double> *pvPF);



double CBinaryMultiObjectiveProblem::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)
{
	CString  s_buf;
	vEvaluateParetoFront(&v_pf_buffer, pcFenotype);

	if (i_weghting_type == MULTI_OBJECTIVE_PROBLEM_WEIGHTING_WEIGHT_VECTOR)
	{
		double  d_weighted_result = 0;

		for (int i_measure = 0; i_measure < v_measures.size(); i_measure++)
		{
			d_weighted_result += v_pf_buffer.at(i_measure) * v_measures.at(i_measure)->dWeight;

			//s_buf.Format("%.2lf * %.2lf = %.2lf", v_pf_buffer.at(i_measure), v_mesaures.at(i_measure)->dWeight, v_pf_buffer.at(i_measure) * v_mesaures.at(i_measure)->dWeight);
			//::Tools::vReportInFile("fit.txt", s_buf);
		}//for (int i_measure = 0; i_measure < v_mesaures.size(); i_measure++)

		//::Tools::vReportInFile("fit.txt", "");
		//::Tools::vShow("abbbb");

		return(d_weighted_result);
	}//if (i_weghting_type == MULTI_OBJECTIVE_PROBLEM_WEIGHTING_WEIGHT_VECTOR)


	if (i_weghting_type == MULTI_OBJECTIVE_PROBLEM_WEIGHTING_TCHEBYCHEV)
	{
		double  d_result_max, d_result_current;

		/*//MOCK
		v_measures.at(0)->dWeight = 0.8888888888888890;
		v_measures.at(1)->dWeight = 0.1111111111111110;

		v_pf_buffer.at(0) = 7;
		v_pf_buffer.at(1) = 18;
		//MOCK END*/


		d_result_max = 0;
		for (int i_measure = 0; i_measure < v_measures.size(); i_measure++)
		{
			d_result_current = (v_measures.at(i_measure)->dWeight * MULTI_OBJECTIVE_PROBLEM_WEIGHTING_TCHEBYCHEV_MULTIPLIER /*2.0 * v_measures.at(i_measure)->dGetMax()*/ - v_pf_buffer.at(i_measure)) * v_measures.at(i_measure)->dWeight;

			if (i_measure == 0)  d_result_max = d_result_current;
			if (d_result_max < d_result_current)  d_result_max = d_result_current;

			//CString  s_buf;
			//s_buf.Format("max: %.16lf current: %.16lf \t weight: %.4lf \t value: %.4lf", d_result_max, d_result_current, v_measures.at(i_measure)->dWeight, v_pf_buffer.at(i_measure));
			//::Tools::vReportInFile("____etst.txt", s_buf);
		}//for (int i_measure = 0; i_measure < v_mesaures.size(); i_measure++)

		d_result_max *= -1;
		//::Tools::vShow(d_result_max);
		return(d_result_max);
	}//if (i_weghting_type == MULTI_OBJECTIVE_PROBLEM_WEIGHTING_TCHEBYCHEV)


	return(-1);
}//double CBinaryMultiObjectiveProblem::d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift)



CString  CBinaryMultiObjectiveProblem::sAdditionalSummaryInfo()
{
	CString  s_result;
	s_result.Format
		(
			"\t GenLen: \t %d \t PFSize: \t %d \t OptPfFound: \t %d \t %d \t OptPfPerc: \t %.2lf \t HypVol: \t %.8lf \t IGD: \t %.8lf \t GD: \t %.8lf \t MaxSpread: \t %.8lf \t DomOptimalPoints: \t %d \t LastPfUpdate: \t %u \t LastPfUpdateTime: \t %.2lf", 
			iGetNumberOfElements(), (int)dPFQualityPointNum(), iPFQualityOptimalPointsFound(), iPFQualityOptimalPointsNum(), dPFQualityOptimalPointsPerc(),
			dPFQualityHyperVolume(), dPFQualityInverseGenerationalDistance(), dPFQualityGenerationalDistance(), dPFQualityMaximumSpread(), (int)dPFQualityDominatedOptimalPoints(), i_last_pareto_front_update_ffe, d_last_pareto_front_update_time
		);


	return(s_result);

	for (int ii = 0; ii < v_global_pareto_front.size(); ii++)
	{
		s_result += "(" + v_global_pareto_front.at(ii)->sReportObj() + ")";
	}//for (int ii = 0; ii < v_global_pareto_front.size(); ii++)

	return(s_result);
}//CString  CBinaryMultiObjectiveProblem::sAdditionalSummaryInfo()


double  CBinaryMultiObjectiveProblem::dPFQualityHyperVolume()
{
	//CBinaryMultiObjectiveProblem *pc_problem_multi;
	//pc_problem_multi = (CBinaryMultiObjectiveProblem *)pc_problem->pcGetEvaluation();


	vector<double>  v_perfect;
	for (int i_measure = 0; i_measure < pvGetMeasures()->size(); i_measure++)
	{
		v_perfect.push_back(pvGetMeasures()->at(i_measure)->dGetMax());
	}//for (int i_measure = 0; i_measure < pc_problem_multi->vGetMeasures()->size(); i_measure++)



	double  d_hyper_volume, d_piece_volume;

	d_hyper_volume = 0;
	for (int ii = 0; ii < v_global_pareto_front.size(); ii++)
	{

		if (ii == 0)
		{
			d_piece_volume = (v_global_pareto_front.at(0)->pvGetFitness()->at(0) - 0) * (v_perfect.at(1) - v_global_pareto_front.at(0)->pvGetFitness()->at(1));
			d_hyper_volume += d_piece_volume;
		}//if (ii == 0)

		if (ii > 0)
		{
			d_piece_volume = (v_global_pareto_front.at(ii)->pvGetFitness()->at(0) - v_global_pareto_front.at(ii - 1)->pvGetFitness()->at(0)) * (v_perfect.at(1) - v_global_pareto_front.at(ii)->pvGetFitness()->at(1));
			d_hyper_volume += d_piece_volume;
		}//if (ii > 0)

		if (ii == v_global_pareto_front.size() - 1)
		{
			d_piece_volume = (v_perfect.at(0) - v_global_pareto_front.at(ii)->pvGetFitness()->at(0)) * (v_perfect.at(1) - 0);
			d_hyper_volume += d_piece_volume;
		}//if (ii == v_global_pareto_front.size() - 1)

	}//for (int ii = 0; ii < v_global_pareto_front.size(); ii++)


	return(d_hyper_volume);
}//double  CBinaryMultiObjectiveProblem::dPFQualityHyperVolume()



double  CBinaryMultiObjectiveProblem::dPFQualityPointNum()
{
	return(v_global_pareto_front.size());
}//double  CBinaryMultiObjectiveProblem::dPFQualityPointNum()




double  CBinaryMultiObjectiveProblem::dPFQualityInverseGenerationalDistance()
{
	if (v_theoretical_optimal_pareto_front.size() == 0)
	{
		return(0);
	}//if (v_theoretical_optimal_pareto_front.size() == 0)



	double  d_inversed_generational_distance, d_cur_dist;
	d_inversed_generational_distance = 0;

	for (int ii = 0; ii < v_theoretical_optimal_pareto_front.size(); ii++)
	{
		d_cur_dist = ((CMultiIndividualPFpoints *)v_theoretical_optimal_pareto_front.at(ii))->dGetInversedGenerationalDistance(&v_global_pareto_front);

		d_inversed_generational_distance += d_cur_dist;
	}//for (int ii = 0; ii < v_theoretical_optimal_pareto_front.size(); ii++)

	d_inversed_generational_distance = d_inversed_generational_distance / v_theoretical_optimal_pareto_front.size();

	return(d_inversed_generational_distance);
}//double  CBinaryMultiObjectiveProblem::dPFQualityInverseGenerationalDistance()




int  CBinaryMultiObjectiveProblem::iPFQualityOptimalPointsFound()
{
	int  i_result;

	if (v_theoretical_optimal_pareto_front.size() == 0)
	{
		return(0);
	}//if (v_theoretical_optimal_pareto_front.size() == 0)



	double  d_cur_dist;
	i_result = 0;

	for (int ii = 0; ii < v_theoretical_optimal_pareto_front.size(); ii++)
	{
		d_cur_dist = ((CMultiIndividualPFpoints *)v_theoretical_optimal_pareto_front.at(ii))->dGetInversedGenerationalDistance(&v_global_pareto_front);
		
		if (d_cur_dist == 0) i_result++;
	}//for (int ii = 0; ii < v_theoretical_optimal_pareto_front.size(); ii++)
	
	return(i_result);
}//double  CBinaryMultiObjectiveProblem::iPFQualityOptimalPointsFound()




double  CBinaryMultiObjectiveProblem::dPFQualityOptimalPointsPerc()
{
	double  d_result;

	d_result = iPFQualityOptimalPointsFound();
	d_result = d_result / iPFQualityOptimalPointsNum();

	return(d_result);
}//double  CBinaryMultiObjectiveProblem::dPFQualityOptimalPointsPerc()




double  CBinaryMultiObjectiveProblem::dPFQualityGenerationalDistance()
{
	if (v_theoretical_optimal_pareto_front.size() == 0)
	{
		return(0);
	}//if (v_theoretical_optimal_pareto_front.size() == 0)



	double  d_generational_distance, d_cur_dist;
	d_generational_distance = 0;

	for (int ii = 0; ii < v_global_pareto_front.size(); ii++)
	{
		d_cur_dist = ((CMultiIndividualPFpoints *)v_global_pareto_front.at(ii))->dGetInversedGenerationalDistance(&v_theoretical_optimal_pareto_front, false);

		d_generational_distance += d_cur_dist;
	}//for (int ii = 0; ii < v_theoretical_optimal_pareto_front.size(); ii++)

	d_generational_distance = d_generational_distance / v_global_pareto_front.size();

	return(d_generational_distance);
}//double  CBinaryMultiObjectiveProblem::dPFQualityGenerationalDistance()




double  CBinaryMultiObjectiveProblem::dPFQualityMaximumSpread()
{
	if (v_global_pareto_front.size() == 0)  return(VERY_LARGE_VALUE);
	if  (v_global_pareto_front.at(0)->pvGetFitness()->size() != 2)  return(VERY_LARGE_VALUE);

	double  d_max_spread;
	d_max_spread = v_global_pareto_front.at(0)->dGetObjectiveSpaceDistance(v_global_pareto_front.at(v_global_pareto_front.size() -1));

	return(d_max_spread);
}//double  CBinaryMultiObjectiveProblem::dPFQualityMaximumSpread()



double  CBinaryMultiObjectiveProblem::dPFQualityDominatedOptimalPoints()
{
	if (v_theoretical_optimal_pareto_front.size() == 0)
	{
		return(0);
	}//if (v_theoretical_optimal_pareto_front.size() == 0)



	int i_dominated_pf_optimal_points;
	i_dominated_pf_optimal_points = 0;

	for (int i_theoretical = 0; i_theoretical < v_theoretical_optimal_pareto_front.size(); i_theoretical++)
	{
		for (int i_current = 0; i_current < v_global_pareto_front.size(); i_current++)
		{
			if (
				v_theoretical_optimal_pareto_front.at(i_theoretical)->iCheckDomination(v_global_pareto_front.at(i_current)) < 0
				)
				i_dominated_pf_optimal_points++;
		}//for (int i_current = 0; i_current < v_global_pareto_front.size(); i_current++)
		
	}//for (int ii = 0; ii < v_theoretical_optimal_pareto_front.size(); ii++)


	return(i_dominated_pf_optimal_points);
}//double  CBinaryMultiObjectiveProblem::dPFQualityDominatedOptimalPoints()




//returns -1 if ind was not added, or the number of thrown-out individuals
int   CBinaryMultiObjectiveProblem::i_join_the_global_pareto(CMultiIndividual *pcInd)
{
	int  i_dominated_individuals;
	int  i_domination;

	i_dominated_individuals = 0;
	for (int ii = 0; ii < v_global_pareto_front.size(); ii++)
	{
		i_domination = pcInd->iCheckDomination(v_global_pareto_front.at(ii));

		if (i_domination < 0)
		{
			if (i_dominated_individuals > 0)  ::Tools::vShow("Dominated candidated, that was dominateing part of the current global-front");
			return(-1);
		}//if (i_domination < 0)


		if (i_domination == 0)
		{
			//if (pcInd->pc_genotype->bGenotypesDiffer(((CNSGA2Individual *) v_global_pareto_front.at(ii))->pc_genotype) == true)
			if (pcInd->bFrontsDiffer(v_global_pareto_front.at(ii)) == true)
			{
				if (i_dominated_individuals > 0)  ::Tools::vShow("The same candidate, that was dominateing part of the current global-front");
				return(-1);
			}//if (pcInd->pc_genotype->bGenotypesDiffer(v_global_pareto_front.at(ii)->pc_genotype) == true)
		}//if (i_domination == 0)

		if (i_domination > 0)
		{
			i_dominated_individuals++;
			delete  v_global_pareto_front.at(ii);
			v_global_pareto_front.erase(v_global_pareto_front.begin() + ii);
			ii--;
		}//if (i_domination > 0)
	}//for (int ii = 0; ii < v_global_pareto_front.size(); ii++)


	CMultiIndividual  *pc_copy_for_global_front;
	pc_copy_for_global_front = pcInd->pcClone();
	v_add_to_non_dominated_pf(&v_global_pareto_front, pc_copy_for_global_front);
	i_last_pareto_front_update_ffe = iGetFFE();
	d_last_pareto_front_update_time = c_time_counter.dGetTimePassed();
	//::Tools::vShow(d_last_pareto_front_update_time);

	return(i_dominated_individuals);
}//int   CBinaryMultiObjectiveProblem::i_join_the_global_pareto(CMultiIndividual *pcInd)



void  CBinaryMultiObjectiveProblem::v_add_to_non_dominated_pf(vector<CMultiIndividual *>  *pvFrontToFill, CMultiIndividual  *pcInd)
{
	for (int ii = 0; ii < pvFrontToFill->size(); ii++)
	{
		if (pcInd->pvGetFitness()->at(0) < pvFrontToFill->at(ii)->pvGetFitness()->at(0))
		{
			pvFrontToFill->insert(pvFrontToFill->begin() + ii, pcInd);
			return;
		}//if (pcInd->pvGetFitness()->at(0) < pvFrontToFill->at(ii)->pvGetFitness()->at(0))
		else
		{
			if (pcInd->pvGetFitness()->at(0) == pvFrontToFill->at(ii)->pvGetFitness()->at(0))
			{
				if (pcInd->pvGetFitness()->size() > 1)
				{
					if (pcInd->pvGetFitness()->at(1) > pvFrontToFill->at(ii)->pvGetFitness()->at(1))
					{
						pvFrontToFill->insert(pvFrontToFill->begin() + ii, pcInd);
						return;
					}//if (pcInd->pvGetFitness()->at(1) > pvFrontToFill->at(ii)->pvGetFitness()->at(1))
				}//if  (pcInd->pvGetFitness()->size() > 1)
			}//if (pcInd->pvGetFitness()->at(0) == pvFrontToFill->at(ii)->pvGetFitness()->at(0))
		}//else  if (pcInd->pvGetFitness()->at(0) > pvFrontToFill->at(ii)->pvGetFitness()->at(0))

	}//for (int ii = 0; ii < pvFrontToFill->size(); ii++)

	pvFrontToFill->push_back(pcInd);
};//void  void  CBinaryMultiObjectiveProblem::v_add_to_non_dominated_pf(vector<CMultiIndividual *>  *pvFrontToFill, CBinaryMultiObjectiveProblem  *pcInd)




CError CBinaryMultiObjectiveProblem::eConfigure(istream *psSettings)
{
	CError  c_err;
	CString  s_buf;


	//not every problem needs that
	//c_err = CBinaryFileConfigEvaluation::eConfigure(psSettings);
	//if (c_err)  return(c_err);


	CString  s_weghting_type;
	CStringCommandParam p_optimal_pf_weight_type(EVALUATION_ARGUMENT_BINARY_MULTI_WEIGHTING, true);
	s_weghting_type = p_optimal_pf_weight_type.sGetValue(psSettings, &c_err);
	
	bool  b_weight_found = false;
	if (s_weghting_type == MULTI_OBJECTIVE_PROBLEM_WEIGHTING_WEIGHT_VECTOR_TEXT)
	{
		vSetWeghtingType(MULTI_OBJECTIVE_PROBLEM_WEIGHTING_WEIGHT_VECTOR);
		b_weight_found = true;
	}//if (s_weghting_type == MULTI_OBJECTIVE_PROBLEM_WEIGHTING_WEIGHT_VECTOR_TEXT)

	if (s_weghting_type == MULTI_OBJECTIVE_PROBLEM_WEIGHTING_TCHEBYCHEV_TEXT)
	{
		vSetWeghtingType(MULTI_OBJECTIVE_PROBLEM_WEIGHTING_TCHEBYCHEV);
		b_weight_found = true;
	}//if (s_weghting_type == MULTI_OBJECTIVE_PROBLEM_WEIGHTING_WEIGHT_VECTOR_TEXT)

	if (b_weight_found == false)
	{
		s_buf.Format("Unknown weighting type '%s'", s_weghting_type);
		c_err.vSetError(s_buf);
	}//if (b_weight_found == false)
	   

	if (c_err)  return(c_err);


	CStringCommandParam p_optimal_pf_file_not_really_optimal(EVALUATION_ARGUMENT_BINARY_MULTI_OPTIMAL_PF_NOT_OPTIMAL, false);
	p_optimal_pf_file_not_really_optimal.sGetValue(psSettings, &c_err);//we only check if the line exists
	if (c_err)  return(c_err);

	if (p_optimal_pf_file_not_really_optimal.bHasValue() == true)  b_theoretical_pareto_front_really_optimal = false;	 


	return(c_err);
}//CError CBinaryMultiObjectiveProblem::eConfigure(istream *psSettings)


//-----------------------------------------CMultiIndividualPFpoints------------------------------------------------------------

CMultiIndividualPFpoints::~CMultiIndividualPFpoints()
{
	if (b_own_solution == false) pc_genotype = NULL; //if genotype is not owned we can not delete it (in ~CMultiIndividual)
}//CMultiIndividualPFpoints::~CMultiIndividualPFpoints()


CMultiIndividual  *CMultiIndividualPFpoints::pcClone()
{
	CMultiIndividualPFpoints  *pc_result;

	pc_result = new CMultiIndividualPFpoints();
	pc_result->vConfigure(pc_genotype, true, &v_fitness, b_rated, pc_problem);
	pc_result->b_matched_by_pf = b_matched_by_pf;

	return(pc_result);
}//CMultiIndividual  *CMultiIndividualPFpoints::pcClone()



void  CMultiIndividualPFpoints::vConfigure(CBinaryCoding *pcRatedSolution, bool  bOwnSolution, vector<double>  *pvPfBuffer, bool  bRated, CBinaryMultiObjectiveProblem *pcProblem)
{
	b_own_solution = bOwnSolution;

	if ( (b_own_solution == true)&&(pcRatedSolution != NULL) )
		pc_genotype = new CBinaryCoding(pcRatedSolution);
	else
		pc_genotype = pcRatedSolution;

	b_rated = bRated;
	if (b_rated == true)  v_fitness = *pvPfBuffer;

	pc_problem = pcProblem;

}//void  CMultiIndividualPFpoints::vConfigure(CBinaryCoding *pcRatedSolution, bool  bOwnSolution, vector<double>  *pvPfBuffer, bool  bRated)



bool CMultiIndividualPFpoints::bLoadFromFilePRW(FILE  *pfSource)
{
	CString  s_line;

	s_line = ::Tools::sReadLine(pfSource);

	if (s_line == "")  return(false);


	double  d_obj_value;
	int  i_index;
	i_index = 0;

	d_obj_value = ::Tools::dExtractFromString(s_line, i_index, &i_index);
	v_fitness.push_back(d_obj_value);
	d_obj_value = ::Tools::dExtractFromString(s_line, i_index, &i_index);
	v_fitness.push_back(d_obj_value);

	//::Tools::vShow(d_obj_value);

	if (i_index < 2)  return(false);

	if (v_fitness.at(0) == 0)
	{
		::Tools::vShow("break");
	}

	//s_line.Format("TEST: %.4lf  %.4lf", v_fitness.at(0), v_fitness.at(1));
	//::Tools::vShow(s_line);

	pc_genotype = NULL;
	b_rated = true;
	b_own_solution = true;
}//bool CMultiIndividualPFpoints::bLoadFromFilePRW(FILE  *pfSource)




bool CMultiIndividualPFpoints::bLoadFromFile_MaxCutMoGomea(FILE  *pfSource)
{
	CString  s_line;

	s_line = ::Tools::sReadLine(pfSource);

	if (s_line == "")  return(false);


	double  d_obj_value;
	int  i_index;
	i_index = 0;

	d_obj_value = ::Tools::dExtractFromString(s_line, i_index, &i_index);
	v_fitness.push_back(d_obj_value);
	d_obj_value = ::Tools::dExtractFromString(s_line, i_index, &i_index);
	v_fitness.push_back(d_obj_value);

	if (i_index < 2)  return(false);

	if (v_fitness.at(0) == 0)
	{
		::Tools::vShow("break");
	}

	//s_line.Format("TEST: %.2lf  %.2lf", v_fitness.at(0), v_fitness.at(1));
	//::Tools::vShow(s_line);

	pc_genotype = NULL;
	b_rated = true;
	b_own_solution = true;
	


	return(true);
}//bool CMultiIndividualPFpoints::bLoadFromFile_MaxCutMoGomea(FILE  *pfSource)



void  CMultiIndividualPFpoints::vRepDominance(CMultiIndividual *pcDominator)
{
	CString  s_report, s_buf;

	s_report.Format("PF [%s]  NEW DOMINATING: [%s]", sReportObj(), pcDominator->sReportObj());

	::Tools::vReportInFile("zzz_report_new_pf.txt", s_report);
}//void  CMultiIndividualPFpoints::vRepDominance(CMultiIndividual *pcDominating)


double  CMultiIndividualPFpoints::dGetInversedGenerationalDistance(vector  <CMultiIndividual *>  *pvGlobalParetoFront, bool  bCheckDominance)
{
	//if (b_matched_by_pf == true)  return(0);
	if (pvGlobalParetoFront->size() < 1)  return(VERY_LARGE_VALUE);

	double d_igd_cur, d_igd_min;

	d_igd_min = dGetObjectiveSpaceDistance(pvGlobalParetoFront->at(0));

	if (bCheckDominance == true)
	{
		if (iCheckDomination(pvGlobalParetoFront->at(0)) < 0)
		{
			d_igd_min *= -1;
			//vRepDominance(pvGlobalParetoFront->at(0));
		}//if (iCheckDomination(pvGlobalParetoFront->at(0)) < 0)
	}//if (bCheckDominance == true)

	for (int ii = 1; ii < pvGlobalParetoFront->size(); ii++)
	{
		d_igd_cur = dGetObjectiveSpaceDistance(pvGlobalParetoFront->at(ii));

		if (bCheckDominance == true)
		{
			if (iCheckDomination(pvGlobalParetoFront->at(ii)) < 0)
			{
				d_igd_cur *= -1;
				//vRepDominance(pvGlobalParetoFront->at(ii));
			}//if (iCheckDomination(pvGlobalParetoFront->at(ii)) < 0)
		}//if (bCheckDominance == true)

		if (d_igd_min > d_igd_cur)  d_igd_min = d_igd_cur;
	}//for (int ii = 0; ii < pvGlobalParetoFront->size(); ii++)


	if (d_igd_min == 0)  b_matched_by_pf = true;


	return(d_igd_min);
};//double  CMultiIndividualPFpoints::dGetInversedGenerationalDistance(vector  <CMultiIndividual *>  *pvGlobalParetoFront)





//-----------------------------------------CMultiIndividual------------------------------------------------------------


CMultiIndividual::~CMultiIndividual()
{
	if  (pc_genotype != NULL)  delete  pc_genotype;
};//CMultiIndividual::~CMultiIndividual()


int  CMultiIndividual::iCheckDomination(CMultiIndividual  *pcOther)
{
	bool  b_better = false;
	bool  b_worse = false;

	vRate();
	pcOther->vRate();

	for (int i_obj = 0; i_obj < v_fitness.size(); i_obj++)
	{
		if (v_fitness.at(i_obj) < pcOther->pvGetFitness()->at(i_obj))
			b_worse = true;
		else
		{
			if (v_fitness.at(i_obj) > pcOther->pvGetFitness()->at(i_obj))  b_better = true;
		}//else  if (v_fitness.at(i_obj) < pcOther->pvGetFitness()->at(i_obj))

		if ((b_better == true) && (b_worse == true))  return(0);
	}//for (int i_obj = 0; i_obj < v_fitness.size(); i_obj++)


	if ((b_better == true) && (b_worse == true))  return(0);
	if ((b_better == true) && (b_worse == false))  return(1);
	if ((b_better == false) && (b_worse == true))  return(-1);

	return(0);
};//CMultiIndividual::iCheckDomination(CNSGA2Individual  *pcOther)



void  CMultiIndividual::vRate()
{
	if (b_rated == true)  return;
	pc_problem->vEvaluateParetoFront(&v_fitness, pc_genotype);
	b_rated = true;
};//void  CMultiIndividual::vRate()



CString  CMultiIndividual::sReport()
{
	CString  s_result, s_buf;

	for (int ii = 0; ii < pc_genotype->iGetNumberOfBits(); ii++)
	{
		s_buf.Format("%d", pc_genotype->piGetBits()[ii]);
		s_result += s_buf;
	}//for (int ii = 0; ii < pc_genotype->iGetNumberOfBits(); ii++)

	vRate();

	for (int ii = 0; ii < v_fitness.size(); ii++)
	{
		s_buf.Format(" M[%d]: %.4lf", ii, v_fitness.at(ii));
		s_result += s_buf;
	}//for (int ii = 0; ii < v_fitness.size(); ii++)

	return(s_result);
}//CString  CMultiIndividual::sReport()



CString  CMultiIndividual::sReportObj()
{
	CString  s_result, s_buf;

	vRate();


	for (int ii = 0; ii < v_fitness.size(); ii++)
	{
		s_buf.Format("%.4lf \t", v_fitness.at(ii));

		s_result += s_buf;
	}//for (int ii = 0; ii < v_fitness.size(); ii++)

	if  (pc_genotype  !=  NULL)  s_result += " \t " + pc_genotype->sToString();


	return(s_result);
}//CString  CMultiIndividual::sReportObj()



double  CMultiIndividual::dGetObjectiveSpaceDistance(CMultiIndividual  *pcOther)
{
	vRate();
	pcOther->vRate();

	if (v_fitness.size() != pcOther->v_fitness.size())  return(-1);

	double  d_os_distance, d_buf;
	d_os_distance = 0;
	for (int  ii = 0; ii < v_fitness.size(); ii++)
	{
		d_buf = v_fitness.at(ii) - pcOther->v_fitness.at(ii);
		d_buf = d_buf * d_buf;
		//d_buf = ::sqrt(d_buf);
	
		d_os_distance += d_buf;
	}//for (int  ii = 0; ii < v_fitness.size(); ii++)

	d_buf = ::sqrt(d_os_distance);


	return(d_os_distance);
}//double  CMultiIndividual::dGetIGD_Distance(CMultiIndividual  *pcOther)



//-----------------------------------------CBinaryMultiPaintsMeasureMakespan------------------------------------------------------------


double CBinaryMultiPaintsMeasureMakespan::dEvaluate(CBinaryCoding *pcFenotype)
{
	int  i_mach_number;

	i_mach_number = ((CBinaryMultiPaints*)pc_parent_problem)->i_machine_number;

	while (v_machines_makespan_buffer.size() < i_mach_number)
		v_machines_makespan_buffer.push_back(0);


	for (int ii = 0; ii < v_machines_makespan_buffer.size(); ii++)
		v_machines_makespan_buffer.at(ii) = 0;
	

	int  i_genes_per_recipe;
	int  i_current_number_of_elements = 0;
	for (int i_recipe = 0; i_recipe < ((CBinaryMultiPaints*)pc_parent_problem)->v_recipes.size(); i_recipe++)
	{
		i_genes_per_recipe = (int)(((CBinaryMultiPaints*)pc_parent_problem)->v_paint_demands.at(((CBinaryMultiPaints*)pc_parent_problem)->v_recipes.at(i_recipe).i_paint_produced_offset) / ((CBinaryMultiPaints*)pc_parent_problem)->v_recipes.at(i_recipe).d_produced_amount);
		//upper bound...
		if (((CBinaryMultiPaints*)pc_parent_problem)->v_paint_demands.at(((CBinaryMultiPaints*)pc_parent_problem)->v_recipes.at(i_recipe).i_paint_produced_offset) > ((CBinaryMultiPaints*)pc_parent_problem)->v_recipes.at(i_recipe).d_produced_amount * i_genes_per_recipe)  i_genes_per_recipe++;

		

		for (int i_mach_off = 0; i_mach_off < ((CBinaryMultiPaints*)pc_parent_problem)->v_recipes.at(i_recipe).v_machines.size(); i_mach_off++)
		{
			for (int i_gene = i_current_number_of_elements; i_gene < i_current_number_of_elements + i_genes_per_recipe; i_gene++)
			{
				if (pcFenotype->piGetBits()[i_gene] == 1)  
					v_machines_makespan_buffer.at(((CBinaryMultiPaints*)pc_parent_problem)->v_recipes.at(i_recipe).v_machines.at(i_mach_off))
					+= ((CBinaryMultiPaints*)pc_parent_problem)->v_recipes.at(i_recipe).d_duration;
			}//for (int i_gene = i_current_number_of_elements; i_gene < i_current_number_of_elements + i_genes_per_recipe; i_genes_per_recipe++)

			i_current_number_of_elements += i_genes_per_recipe;
		}//for  (int  i_mach_off = 0; i_mach_off < ((CBinaryMultiPaints*)pc_parent_problem)->v_recipes.at(i_recipe).v_machines.size(); i_mach_off++)
				
	}//for (int i_recipe = 0; i_recipe < v_recipes.size(); i_recipe++)

	
	double  d_max;
	d_max = v_machines_makespan_buffer.at(0);

	for (int ii = 1; ii < v_machines_makespan_buffer.size(); ii++)
	{
		if (d_max < v_machines_makespan_buffer.at(ii))  d_max = v_machines_makespan_buffer.at(ii);
	}//for (int ii = 1; v_machines_makespan_buffer.size(); ii++)


	return(1.0/(d_max + 1));
}//double CBinaryMultiPaintsMeasureMakespan::dEvaluate(CBinaryCoding *pcFenotype)



double CBinaryMultiPaintsMeasureOverhead::dEvaluate(CBinaryCoding *pcFenotype)
{
	int  i_paints;

	i_paints = ((CBinaryMultiPaints*)pc_parent_problem)->v_paint_demands.size();

	while (v_paint_amount_buffer.size() < i_paints)
		v_paint_amount_buffer.push_back(0);


	for (int ii = 0; ii < v_paint_amount_buffer.size(); ii++)
		v_paint_amount_buffer.at(ii) = 0;

	//constraint handling first - we must satisfy all demands
	int  i_current_number_of_elements;
	int  i_genes_per_recipe;

	i_current_number_of_elements = 0;

	for (int i_recipe = 0; i_recipe < ((CBinaryMultiPaints*)pc_parent_problem)->v_recipes.size(); i_recipe++)
	{
		i_genes_per_recipe = (int)(((CBinaryMultiPaints*)pc_parent_problem)->v_paint_demands.at(((CBinaryMultiPaints*)pc_parent_problem)->v_recipes.at(i_recipe).i_paint_produced_offset) / ((CBinaryMultiPaints*)pc_parent_problem)->v_recipes.at(i_recipe).d_produced_amount);
		//upper bound...
		if (((CBinaryMultiPaints*)pc_parent_problem)->v_paint_demands.at(((CBinaryMultiPaints*)pc_parent_problem)->v_recipes.at(i_recipe).i_paint_produced_offset) > ((CBinaryMultiPaints*)pc_parent_problem)->v_recipes.at(i_recipe).d_produced_amount * i_genes_per_recipe)  i_genes_per_recipe++;

		i_genes_per_recipe = i_genes_per_recipe * ((CBinaryMultiPaints*)pc_parent_problem)->v_recipes.at(i_recipe).v_machines.size();


		//if (v_recipes.at(i_recipe).i_paint_produced_offset == i_paint)
		for (int i_gene = i_current_number_of_elements; i_gene < i_current_number_of_elements + i_genes_per_recipe; i_gene++)
		{
			if (pcFenotype->piGetBits()[i_gene] == 1)  v_paint_amount_buffer.at(((CBinaryMultiPaints*)pc_parent_problem)->v_recipes.at(i_recipe).i_paint_produced_offset) += ((CBinaryMultiPaints*)pc_parent_problem)->v_recipes.at(i_recipe).d_produced_amount;
		}//for (int i_gene = i_current_number_of_elements; i_gene < i_current_number_of_elements + i_genes_per_recipe; i_genes_per_recipe++)

		i_current_number_of_elements += i_genes_per_recipe;
	}//for (int i_recipe = 0; i_recipe < v_recipes.size(); i_recipe++)



	double  d_max_overhead;

	d_max_overhead = v_paint_amount_buffer.at(0) - ((CBinaryMultiPaints*)pc_parent_problem)->v_paint_demands.at(0);

	/*CString  s_buf;
	//::Tools::vReportInFile("over.txt", pcFenotype->sToString());
	s_buf.Format("%.2lf - %.2lf = %.2lf", v_paint_amount_buffer.at(0), ((CBinaryMultiPaints*)pc_parent_problem)->v_paint_demands.at(0), v_paint_amount_buffer.at(0) - ((CBinaryMultiPaints*)pc_parent_problem)->v_paint_demands.at(0));
	::Tools::vReportInFile("over.txt", s_buf);
	s_buf.Format("%.2lf - %.2lf = %.2lf", v_paint_amount_buffer.at(1), ((CBinaryMultiPaints*)pc_parent_problem)->v_paint_demands.at(1), v_paint_amount_buffer.at(1) - ((CBinaryMultiPaints*)pc_parent_problem)->v_paint_demands.at(1));
	::Tools::vReportInFile("over.txt", s_buf);
	::Tools::vShow("b");*/

	for (int ii = 1; ii < v_paint_amount_buffer.size(); ii++)
	{
		if (d_max_overhead < v_paint_amount_buffer.at(ii) - ((CBinaryMultiPaints*)pc_parent_problem)->v_paint_demands.at(ii))
			d_max_overhead = v_paint_amount_buffer.at(ii) - ((CBinaryMultiPaints*)pc_parent_problem)->v_paint_demands.at(ii);
	}//for (int ii = 1; ii < v_paint_amount_buffer.size(); ii++)

	return(1.0 / (d_max_overhead + 1));
}//double CBinaryMultiPaintsMeasureOverhead::dEvaluate(CBinaryCoding *pcFenotype)







uint32_t CBinaryMultiFuncTool::iERROR_PARENT_CBinaryMultiFuncTool = CError::iADD_ERROR_PARENT("iERROR_PARENT_CBinaryMultiFuncLinear");
uint32_t CBinaryMultiFuncTool::iERROR_CODE_MNULTI_OBJECTIVE_PAINTS_FUNC_UNKNOWN_TYPE = CError::iADD_ERROR("iERROR_CODE_MNULTI_OBJECTIVE_PAINTS_FUNC_UNKNOWN_TYPE");



CError  CBinaryMultiFuncToolFuncLinear::eConfig(CString  sLine)
{
	CError  c_err(CBinaryMultiFuncTool::iERROR_PARENT_CBinaryMultiFuncTool);
	CString  s_marker;

	s_marker = s_MULTI_PAINTS_CONFIG_FUNC_TYPE_LINEAR;
	d_multiplier = ::Tools::dExtractFromString(sLine, s_marker.GetLength());

	return(c_err);
}//CError  CBinaryMultiFuncToolFuncLinear::eConfig(CString  sLine)


CError  CBinaryMultiFuncToolFuncExp::eConfig(CString  sLine)
{
	CError  c_err(CBinaryMultiFuncTool::iERROR_PARENT_CBinaryMultiFuncTool);
	CString  s_marker;

	s_marker = s_MULTI_PAINTS_CONFIG_FUNC_TYPE_LINEAR;
	d_base = ::Tools::dExtractFromString(sLine, s_marker.GetLength());

	return(c_err);
}//CError  CBinaryMultiFuncToolFuncExp::eConfig(CString  sLine)


CError  CBinaryMultiFuncToolFuncPow::eConfig(CString  sLine)
{
	CError  c_err(CBinaryMultiFuncTool::iERROR_PARENT_CBinaryMultiFuncTool);
	CString  s_marker;

	s_marker = s_MULTI_PAINTS_CONFIG_FUNC_TYPE_LINEAR;
	d_pow = ::Tools::dExtractFromString(sLine, s_marker.GetLength());

	return(c_err);
}//CError  CBinaryMultiFuncToolFuncPow::eConfig(CString  sLine)



double  CBinaryMultiFuncToolFuncExp::dGetValue(double  dArg)
{
	double  d_res;
	int  i_pow;


	i_pow = (int) dArg;

	d_res = 1;
	for (int ii = 0; ii < i_pow; ii++)
		d_res *= d_base;

	/*CString  s_buf;
	s_buf.Format("res: %.2lf  base:%.2lf  pow: %d", d_res, d_base, i_pow);
	::Tools::vShow(s_buf);*/

	return(d_res);
}//double  CBinaryMultiFuncToolFuncExp::dGetValue(double  dArg)




double  CBinaryMultiFuncToolFuncPow::dGetValue(double  dArg)
{
	double  d_res;
	int  i_pow;


	i_pow = (int) d_pow;

	d_res = 1;
	for (int ii = 0; ii < i_pow; ii++)
		d_res *= dArg;

	return(d_res);
}//double  CBinaryMultiFuncToolFuncPow::dGetValue(double  dArg)



CError  CBinaryMultiFuncTool::eConfig(CString  sLine)
{
	CError  c_err(iERROR_PARENT_CBinaryMultiFuncTool);
	double  d_buf;

	if (pc_func != NULL)  delete  pc_func;

	if (sLine.Find(s_MULTI_PAINTS_CONFIG_FUNC_TYPE_LINEAR) >= 0)
	{
		pc_func = new CBinaryMultiFuncToolFuncLinear();
		c_err = pc_func->eConfig(sLine);
		return(c_err);
	}//if (i_offset = sLine.Find(s_MULTI_PAINTS_CONFIG_FUNC_TYPE_LINEAR) >= 0)


	if (sLine.Find(s_MULTI_PAINTS_CONFIG_FUNC_TYPE_EXP) >= 0)
	{
		pc_func = new CBinaryMultiFuncToolFuncExp();
		c_err = pc_func->eConfig(sLine);
		return(c_err);
	}//if (sLine.Find(s_MULTI_PAINTS_CONFIG_FUNC_TYPE_EXP) >= 0)


	if (sLine.Find(s_MULTI_PAINTS_CONFIG_FUNC_TYPE_POW) >= 0)
	{
		pc_func = new CBinaryMultiFuncToolFuncPow();
		c_err = pc_func->eConfig(sLine);
		return(c_err);
	}//if (sLine.Find(s_MULTI_PAINTS_CONFIG_FUNC_TYPE_POW) >= 0)


	return(c_err);
}//CError  CBinaryMultiFuncLinear::eConfig(CString  sLine)


CString  CBinaryMultiFuncTool::sToString()
{
	if (pc_func == NULL)  return("NO FNCTION");
	return(pc_func->sToString());
}//CString  CBinaryMultiFuncLinear::sToString()



double  CBinaryMultiFuncTool::dGetValue(double  dArg)
{
	if (pc_func == NULL)  return(0);
	return(pc_func->dGetValue(dArg));
}//double  CBinaryMultiFuncTool::dGetValue(double  dArg)







uint32_t CBinaryMultiFuncAmountTool::iERROR_PARENT_CBinaryMultiFuncAmountTool = CError::iADD_ERROR_PARENT("iERROR_PARENT_CBinaryMultiFuncAmountTool");
uint32_t CBinaryMultiFuncAmountTool::iERROR_CODE_MNULTI_OBJECTIVE_PAINTS_AMOUNT_FUNC_UNKNOWN_TYPE = CError::iADD_ERROR("iERROR_CODE_MNULTI_OBJECTIVE_PAINTS_AMOUNT_FUNC_UNKNOWN_TYPE");



CError  CBinaryMultiFuncToolAmountFuncSum::eConfig(CString  sLine)
{
	CError  c_err(CBinaryMultiFuncTool::iERROR_PARENT_CBinaryMultiFuncTool);
	CString  s_marker;

	s_marker = s_MULTI_PAINTS_CONFIG_AMOUNT_FUNC_TYPE_SUM;
	d_multiplier = ::Tools::dExtractFromString(sLine, s_marker.GetLength());

	return(c_err);
}//CError  CBinaryMultiFuncToolAmountFuncSum::eConfig(CString  sLine)


CError  CBinaryMultiFuncToolAmountFuncLargest::eConfig(CString  sLine)
{
	CError  c_err(CBinaryMultiFuncTool::iERROR_PARENT_CBinaryMultiFuncTool);
	CString  s_marker;

	s_marker = s_MULTI_PAINTS_CONFIG_AMOUNT_FUNC_TYPE_LARGEST;
	d_multiplier = ::Tools::dExtractFromString(sLine, s_marker.GetLength());

	return(c_err);
}//CError  CBinaryMultiFuncToolAmountFuncLargest::eConfig(CString  sLine)


double  CBinaryMultiFuncToolAmountFuncSum::dGetValue(vector<double>  *pvProducedAmounts)
{
	double  d_result;

	d_result = 0;
	for (int ii = 0; ii < pvProducedAmounts->size(); ii++)
		d_result += pvProducedAmounts->at(ii);

	d_result *= d_multiplier;

	return(d_result);
}//double  CBinaryMultiFuncToolAmountFuncSum::dGetValue(vector<double>  *pvProducedAmounts)



double  CBinaryMultiFuncToolAmountFuncLargest::dGetValue(vector<double>  *pvProducedAmounts)
{
	double  d_result;

	d_result = 0;
	for (int ii = 0; ii < pvProducedAmounts->size(); ii++)
	{
		if  (d_result < pvProducedAmounts->at(ii))  d_result = pvProducedAmounts->at(ii);
	}//for (int ii = 0; ii < pvProducedAmounts->size(); ii++)

	d_result *= d_multiplier;

	return(d_result);
}//double  CBinaryMultiFuncToolAmountFuncLargest::dGetValue(vector<double>  *pvProducedAmounts)




CError  CBinaryMultiFuncAmountTool::eConfig(CString  sLine)
{
	CError  c_err(iERROR_PARENT_CBinaryMultiFuncAmountTool);
	double  d_buf;

	if (pc_func != NULL)  delete  pc_func;

	if (sLine.Find(s_MULTI_PAINTS_CONFIG_AMOUNT_FUNC_TYPE_SUM) >= 0)
	{
		pc_func = new CBinaryMultiFuncToolAmountFuncSum();
		c_err = pc_func->eConfig(sLine);
		return(c_err);
	}//if (i_offset = sLine.Find(s_MULTI_PAINTS_CONFIG_FUNC_TYPE_LINEAR) >= 0)


	if (sLine.Find(s_MULTI_PAINTS_CONFIG_AMOUNT_FUNC_TYPE_LARGEST) >= 0)
	{
		pc_func = new CBinaryMultiFuncToolAmountFuncLargest();
		c_err = pc_func->eConfig(sLine);
		return(c_err);
	}//if (i_offset = sLine.Find(s_MULTI_PAINTS_CONFIG_FUNC_TYPE_LINEAR) >= 0)


	return(c_err);
}//CError  CBinaryMultiFuncLinear::eConfig(CString  sLine)


CString  CBinaryMultiFuncAmountTool::sToString()
{
	if (pc_func == NULL)  return("NO FNCTION");
	return(pc_func->sToString());
}//CString  CBinaryMultiFuncLinear::sToString()



double  CBinaryMultiFuncAmountTool::dGetValue(vector<double>  *pvProducedAmounts)
{
	if (pc_func == NULL)  return(0);
	return(pc_func->dGetValue(pvProducedAmounts));
}//double  CBinaryMultiFuncAmountTool::dGetValue(vector<double>  *pvProducedAmounts);







uint32_t CBinaryMultiPaints::iERROR_PARENT_CBinaryMultiPaints = CError::iADD_ERROR_PARENT("iERROR_PARENT_CBinaryMultiPaints");
uint32_t CBinaryMultiPaints::iERROR_CODE_MNULTI_OBJECTIVE_PAINTS_UNKNOWN_CONFIG_FILE = CError::iADD_ERROR("iERROR_CODE_MNULTI_OBJECTIVE_PAINTS_UNKNOWN_CONFIG_FILE");

CBinaryMultiPaints::CBinaryMultiPaints(): CBinaryMultiObjectiveProblem::CBinaryMultiObjectiveProblem()
{
	pc_solution_buf = NULL;
}//CBinaryMultiPaints::CBinaryMultiPaints()


CBinaryMultiPaints::CBinaryMultiPaints(const CBinaryMultiPaints &pcOther) : CBinaryMultiObjectiveProblem::CBinaryMultiObjectiveProblem(pcOther)
{
	if (pcOther.pc_solution_buf != NULL)
		pc_solution_buf = new CBinaryCoding(pcOther.pc_solution_buf->iGetNumberOfBits());
	else
		pc_solution_buf = NULL;
}//CBinaryMultiPaints::CBinaryMultiPaints(const CBinaryMultiPaints &pcOther)




CBinaryMultiPaints::CBinaryMultiPaints(FILE *pfConfig, CError *pcError)
{
	*pcError = e_init(pfConfig);
}//CBinaryMultiPaints::CBinaryMultiPaints(FILE *pfConfig, CError *pcError)

CBinaryMultiPaints::~CBinaryMultiPaints()
{
	if (pc_solution_buf != NULL)  delete  pc_solution_buf;
	//delete pc_deceptive_concatenation_function;
}//CBinaryMultiPaints::~CBinaryMultiPaints()



void  CBinaryMultiPaints::vEvaluateParetoFront(vector<double> *pvPF, CBinaryCoding *pcFenotype)
{
	v_prepare_solution(pcFenotype);
	v_evaluate_pareto_front(pvPF, pc_solution_buf);
}//void  CBinaryMultiPaints::vEvaluateParetoFront(vector<double> *pvPF, CBinaryCoding *pcFenotype)






void  CBinaryMultiPaints::v_prepare_solution(CBinaryCoding *pcFenotype)
{
	if (pc_solution_buf == NULL)  pc_solution_buf = new CBinaryCoding(pcFenotype->iGetNumberOfBits());

	for (int ii = 0; ii < pcFenotype->iGetNumberOfBits(); ii++)
		pc_solution_buf->piGetBits()[ii] = pcFenotype->piGetBits()[ii];
		//pc_solution_buf->piGetBits()[ii] = 0;

	/*pc_solution_buf->piGetBits()[0] = 1;
	pc_solution_buf->piGetBits()[1] = 1;
	pc_solution_buf->piGetBits()[2] = 1;
	pc_solution_buf->piGetBits()[8] = 1;
	pc_solution_buf->piGetBits()[14] = 1;
	pc_solution_buf->piGetBits()[15] = 1;
	pc_solution_buf->piGetBits()[16] = 1;*/




	//constraint handling first - we must satisfy all demands
	int  i_current_number_of_elements;
	int  i_genes_per_recipe;
	double  d_paint_sum;

	for (int i_paint = 0; i_paint < v_paint_demands.size(); i_paint++)
	{
		i_current_number_of_elements = 0;
		d_paint_sum = 0;
		for (int i_recipe = 0; i_recipe < v_recipes.size(); i_recipe++)
		{
			i_genes_per_recipe = (int)(v_paint_demands.at(v_recipes.at(i_recipe).i_paint_produced_offset) / v_recipes.at(i_recipe).d_produced_amount);
			//upper bound...
			if (v_paint_demands.at(v_recipes.at(i_recipe).i_paint_produced_offset) > v_recipes.at(i_recipe).d_produced_amount * i_genes_per_recipe)  i_genes_per_recipe++;

			i_genes_per_recipe = i_genes_per_recipe * v_recipes.at(i_recipe).v_machines.size();


			if (v_recipes.at(i_recipe).i_paint_produced_offset == i_paint)
			{
				for (int i_gene = i_current_number_of_elements; i_gene < i_current_number_of_elements + i_genes_per_recipe; i_gene++)
				{
					if (pc_solution_buf->piGetBits()[i_gene] == 1)  d_paint_sum += v_recipes.at(i_recipe).d_produced_amount;
				}//for (int i_gene = i_current_number_of_elements; i_gene < i_current_number_of_elements + i_genes_per_recipe; i_genes_per_recipe++)
			}//if (v_recipes.at(i_recipe).i_paint_produced_offset == i_paint)

			i_current_number_of_elements += i_genes_per_recipe;
		}//for (int i_recipe = 0; i_recipe < v_recipes.size(); i_recipe++)


		//check if you can decrease paint amount
		if (d_paint_sum > v_paint_demands.at(i_paint))
		{
			i_current_number_of_elements = 0;
			for (int i_recipe = 0; i_recipe < v_recipes.size(); i_recipe++)
			{
				i_genes_per_recipe = (int)(v_paint_demands.at(v_recipes.at(i_recipe).i_paint_produced_offset) / v_recipes.at(i_recipe).d_produced_amount);
				//upper bound...
				if (v_paint_demands.at(v_recipes.at(i_recipe).i_paint_produced_offset) > v_recipes.at(i_recipe).d_produced_amount * i_genes_per_recipe)  i_genes_per_recipe++;

				i_genes_per_recipe = i_genes_per_recipe * v_recipes.at(i_recipe).v_machines.size();


				if (v_recipes.at(i_recipe).i_paint_produced_offset == i_paint)
				{
					for (int i_gene = i_current_number_of_elements; i_gene < i_current_number_of_elements + i_genes_per_recipe; i_gene++)
					{
						if (pc_solution_buf->piGetBits()[i_gene] == 1)
						{
							if (d_paint_sum - v_recipes.at(i_recipe).d_produced_amount >= v_paint_demands.at(i_paint))
							{
								pc_solution_buf->piGetBits()[i_gene] = 0;
								d_paint_sum -= v_recipes.at(i_recipe).d_produced_amount;
							}//if (d_paint_sum - v_recipes.at(i_recipe).d_produced_amount > v_paint_demands.at(i_paint))							
						}//if (c_solution_buf.piGetBits()[i_gene] == 1)
					}//for (int i_gene = i_current_number_of_elements; i_gene < i_current_number_of_elements + i_genes_per_recipe; i_genes_per_recipe++)
				}//if (v_recipes.at(i_recipe).i_paint_produced_offset == i_paint)

				i_current_number_of_elements += i_genes_per_recipe;
			}//for (int i_recipe = 0; i_recipe < v_recipes.size(); i_recipe++)
		}//if (d_paint_sum > v_paint_demands.at(i_paint))


		//increase paint amount
		if (d_paint_sum < v_paint_demands.at(i_paint))
		{
			i_current_number_of_elements = 0;
			for (int i_recipe = 0; i_recipe < v_recipes.size(); i_recipe++)
			{
				i_genes_per_recipe = (int)(v_paint_demands.at(v_recipes.at(i_recipe).i_paint_produced_offset) / v_recipes.at(i_recipe).d_produced_amount);
				//upper bound...
				if (v_paint_demands.at(v_recipes.at(i_recipe).i_paint_produced_offset) > v_recipes.at(i_recipe).d_produced_amount * i_genes_per_recipe)  i_genes_per_recipe++;

				i_genes_per_recipe = i_genes_per_recipe * v_recipes.at(i_recipe).v_machines.size();


				if (v_recipes.at(i_recipe).i_paint_produced_offset == i_paint)
				{
					for (int i_gene = i_current_number_of_elements; i_gene < i_current_number_of_elements + i_genes_per_recipe; i_gene++)
					{
						if (pc_solution_buf->piGetBits()[i_gene] == 0)
						{
							if (d_paint_sum < v_paint_demands.at(i_paint))
							{
								pc_solution_buf->piGetBits()[i_gene] = 1;
								d_paint_sum += v_recipes.at(i_recipe).d_produced_amount;
							}//if (d_paint_sum < v_recipes.at(i_recipe).d_produced_amount)
						}//if (c_solution_buf.piGetBits()[i_gene] == 1)
					}//for (int i_gene = i_current_number_of_elements; i_gene < i_current_number_of_elements + i_genes_per_recipe; i_genes_per_recipe++)
				}//if (v_recipes.at(i_recipe).i_paint_produced_offset == i_paint)

				i_current_number_of_elements += i_genes_per_recipe;
			}//for (int i_recipe = 0; i_recipe < v_recipes.size(); i_recipe++)
		}//if (d_paint_sum < v_paint_demands.at(i_paint))

	}//for  (int  i_paint = 0; i_paint < v_paint_demands.size(); i_paint++)

	//::Tools::vReportInFile("zz_results.txt", pc_solution_buf->sToString());
	//::Tools::vShow("aaa");
	
}//void  CBinaryMultiPaints::v_prepare_solution()



CError CBinaryMultiPaints::e_init(FILE *pfConfig)
{
	CError  c_err;


	i_number_of_elements = 1;
	d_max_value = 1;

	c_err = e_load_settings_from_file(pfConfig);

	//eSave("zzz_paints_gen_test_base.txt");
	//aaaaaa
	//eSave("zzz_paints_gen_test_multiplied.txt");

	//::Tools::vShow("aaaaa");


	if (!c_err)
	{
		int  i_genes_per_recipe;

		i_number_of_elements = 0;
		for (int i_recipe = 0; i_recipe < v_recipes.size(); i_recipe++)
		{
			i_genes_per_recipe = (int) (v_paint_demands.at(v_recipes.at(i_recipe).i_paint_produced_offset) / v_recipes.at(i_recipe).d_produced_amount);

			//upper bound...
			if (v_paint_demands.at(v_recipes.at(i_recipe).i_paint_produced_offset) > v_recipes.at(i_recipe).d_produced_amount * i_genes_per_recipe)  i_genes_per_recipe++;
			i_genes_per_recipe = i_genes_per_recipe * v_recipes.at(i_recipe).v_machines.size();

			i_number_of_elements += i_genes_per_recipe;
		}//for (int i_recipe = 0; i_recipe < v_recipes.size(); i_recipe++)
				
		d_max_value = 1;
	}//if (!c_error)*/


	//::Tools::vShow(i_number_of_elements);


	//::Tools::vShow("Mwahahaha");
	fclose(pfConfig);

	return c_err;
}//CError CBinaryMultiPaints::e_init(FILE *pfConfig)



void  CBinaryMultiPaints::v_get_files_names_from_summary(CString  sDir, CString  sSummary, vector<CString>  *pvFileNames)
{
	int  i_index;
	FILE  *pf_summary;
	CString  s_buf, s_line;
	
	pf_summary = fopen(sDir + sSummary, "r");

	while (!feof(pf_summary))
	{
		s_line = ::Tools::sReadLine(pf_summary);
		s_line.Replace(" ", "");

		if (s_line != "")
		{
			i_index = 0;
			s_buf = Tools::sExtractFromString(s_line, i_index, &i_index);
			while ((s_buf != "LogFile:") && (s_buf != ""))
			{
				s_buf = Tools::sExtractFromString(s_line, i_index, &i_index);
				//::Tools::vShow("o" + s_buf + "o");
			}//while (s_buf != "LogFile:")


			if (s_buf != "")
			{
				s_buf = Tools::sExtractFromString(s_line, i_index, &i_index);
				s_buf = sDir + s_buf + "pf.txt";

				pvFileNames->push_back(s_buf);
				//::Tools::vShow(s_buf);
			}//if (s_buf != "")
		}//if (s_buf != "")
	}//while (!feof(pf_problem_list))

	fclose(pf_summary);

	//for (int ii = 0; ii < pvFileNames->size(); ii++)
		//Tools::vShow(pvFileNames->at(ii));

}//void  CBinaryMultiPaints::vGetFilesNamesFromSummary(CString  sDir, CString  sSummary, vector<CString>  *pvFileNames)



CError  CBinaryMultiPaints::eLoadPFFromManyFileFiles()
{
	CError  c_err(iERROR_PARENT_CBinaryMultiPaints);

	CBinaryMultiPaintsMeasureMakespan  *pc_makespan;
	pc_makespan = new CBinaryMultiPaintsMeasureMakespan(this);
	pc_makespan->dWeight = 0.5;
	v_measures.push_back(pc_makespan);

	CBinaryMultiPaintsMeasureOverhead  *pc_overhead;
	pc_overhead = new CBinaryMultiPaintsMeasureOverhead(this);
	pc_overhead->dWeight = 0.5;
	v_measures.push_back(pc_overhead);



	CString  s_buf, s_line;
	FILE  *pf_problem_list;
	vector<CString>  v_problems;
	vector<CString>  v_pf_files_list;
	int i_index;


	//pf_problem_list = fopen("zzz_serie_paints_list.txt","r");
	//pf_problem_list = fopen("C:\\Projekty\\projekty_pwr\\island_ga\\x64\\Release\\zzz_serie_paints_list_old.txt", "r");
	
	//pf_problem_list = fopen("C:\\Projekty\\projekty_pwr\\island_ga\\x64\\Release\\zzz_serie_paints_list_all3.txt", "r");
	pf_problem_list = fopen("C:\\Projekty\\projekty_pwr\\island_ga\\x64\\Release\\zzz_serie_paints_list_all3_multiple.txt", "r");


	while (!feof(pf_problem_list))
	{
		s_buf = ::Tools::sReadLine(pf_problem_list);
		s_buf.Replace(" ", "");

		if  ( (s_buf != "")&&(s_buf != "a") )  v_problems.push_back(s_buf);
	}//while (!feof(pf_problem_list))

	fclose(pf_problem_list);

	//for (int ii = 0; ii < v_problems.size(); ii++)
		//Tools::vShow(v_problems.at(ii));

	/*::Tools::vReportInFile("zzzz_problems.txt", "", true);
	for (int ii = 0; ii < v_problems.size(); ii++)
		::Tools::vReportInFile("zzzz_problems.txt", v_problems.at(ii));//*/


	/*v_get_files_names_from_summary
	(
		"D:\\studia\\doktoranckie\\publikacje\\X 2019.02 Multi-objective paint problem (York)\\205_badania_paints_mogo\\paints_nsga2\\",
		"summary_paints_nsga400pop.txt", &v_pf_files_list
	);*/


	
	
	
	
	
	v_get_files_names_from_summary
		(
		//"D:\\studia\\doktoranckie\\publikacje\\R 2019.02 Multi-objective paint problem (York)\\205_badania_paints_mogo\\paints_mogo\\",
		//"D:\\studia\\doktoranckie\\publikacje\\X 2019.02 Multi-objective paint problem (York)\\research\\experiments\\paints_ok\\",


		//"D:\\studia\\doktoranckie\\publikacje\\R 2019.02 Multi-objective paint problem (York)\\research\\experiments\\paints_mogo\\",
		//"summary_paints_mogomea.txt", 

		
		"D:\\studia\\doktoranckie\\publikacje\\R 2019.02 Multi-objective paint problem (York)\\205_badania_paints_mogo\\multiplied_paints_many_halls\\multiplied_paints_many_halls_mogo\\",
		"summary_paints_mogomea.txt",
	    &v_pf_files_list
		//"D:\\studia\\doktoranckie\\publikacje\\X 2019.02 Multi-objective paint problem (York)\\205_badania_paints_mogo\\paints_mop3_smart\\",
		//"summary_paints_mop3new.txt", &v_pf_files_list


		);
	

	v_get_files_names_from_summary
		(
			"D:\\studia\\doktoranckie\\publikacje\\R 2019.02 Multi-objective paint problem (York)\\205_badania_paints_mogo\\multiplied_paints_many_halls\\multiplied_paints_many_halls_mop3\\",
			"summary_paints_mop3.txt", &v_pf_files_list
		);

	/*v_get_files_names_from_summary
	(
		"D:\\studia\\doktoranckie\\publikacje\\R 2019.02 Multi-objective paint problem (York)\\205_badania_paints_mogo\\multiplied_paints_many_halls\\multiplied_paints_many_halls_mop3_smart\\",
		"summary_paints_mop3_smart.txt", &v_pf_files_list
	);*/




	v_get_files_names_from_summary
	(
		"D:\\studia\\doktoranckie\\publikacje\\R 2019.02 Multi-objective paint problem (York)\\205_badania_paints_mogo\\multiplied_paints_many_halls\\multiplied_paints_many_halls_moead\\",
		"summary_paints_moead.txt", &v_pf_files_list

		//"D:\\studia\\doktoranckie\\publikacje\\R 2019.02 Multi-objective paint problem (York)\\205_badania_paints_mogo\\paints_moead\\",
		//"summary_paints_moead400.txt", &v_pf_files_list
	);


	v_get_files_names_from_summary
	(
		"D:\\studia\\doktoranckie\\publikacje\\R 2019.02 Multi-objective paint problem (York)\\205_badania_paints_mogo\\multiplied_paints_many_halls\\multiplied_paints_many_halls_nsga2\\",
		"summary_paints_nsga2.txt", &v_pf_files_list

		//"D:\\studia\\doktoranckie\\publikacje\\R 2019.02 Multi-objective paint problem (York)\\205_badania_paints_mogo\\paints_moead\\",
		//"summary_paints_moead400.txt", &v_pf_files_list
	);


	


	/*::Tools::vReportInFile("zzzz_test.txt", "", true);
	for (int ii = 0; ii < v_pf_files_list.size(); ii++)
		::Tools::vReportInFile("zzzz_test.txt", v_pf_files_list.at(ii));//*/

	vector<CString>  v_result;
	vector<CString>  v_pseudo_theoretical_report_file;
	FILE  *pf_loc_pareto;
	for (int i_problem = 0; i_problem < v_problems.size(); i_problem++)
	{
		//first - create pareto front
		v_theoretical_optimal_pareto_front.clear();
		v_global_pareto_front.clear();


		for (int i_pf_file = 0; i_pf_file < v_pf_files_list.size(); i_pf_file++)
		{
			if (v_pf_files_list.at(i_pf_file).Find(v_problems.at(i_problem)) > 0)
			{
				::Tools::vReportInFile("zzzz_test223.txt", v_pf_files_list.at(i_pf_file));
				pf_loc_pareto = fopen(v_pf_files_list.at(i_pf_file), "r");
				eLoadPFFromFile(pf_loc_pareto);
				fclose(pf_loc_pareto);
			}//if (v_pf_files_list.at(i_pf_file).Find(v_problems.at(i_problem)) > 0)

		}//for (int i_pf_file = 0; i_pf_file < v_pf_files_list.size(); i_pf_file++)
		

		/*v_pseudo_theoretical_report_file.clear();
		vReportPF(&v_pseudo_theoretical_report_file);

		::Tools::vReportInFile("zzzz_pseudo_theoretical.txt", v_problems.at(i_problem), true);
		for (int ii = 0; ii < v_pseudo_theoretical_report_file.size(); ii++)
			::Tools::vReportInFile("zzzz_pseudo_theoretical.txt", v_pseudo_theoretical_report_file.at(ii));//*/

		for (int ii = 0; ii < v_global_pareto_front.size(); ii++)
			v_theoretical_optimal_pareto_front.push_back(v_global_pareto_front.at(ii));

		v_global_pareto_front.clear();


		
		//second get IGD for all test cases
		for (int i_pf_file = 0; i_pf_file < v_pf_files_list.size(); i_pf_file++)
		{
			if (v_pf_files_list.at(i_pf_file).Find(v_problems.at(i_problem)) > 0)
			{
				v_global_pareto_front.clear();

				pf_loc_pareto = fopen(v_pf_files_list.at(i_pf_file), "r");
				eLoadPFFromFile(pf_loc_pareto);

				//::Tools::vShow(2.5);
				//::Tools::vReportInFile("zzzz_result.txt", v_pf_files_list.at(i_pf_file) + "\t" + this->sAdditionalSummaryInfo());
				v_result.push_back(v_pf_files_list.at(i_pf_file) + "\t" + this->sAdditionalSummaryInfo());
				//::Tools::vShow(4.55625);
				fclose(pf_loc_pareto);
			}//if (v_pf_files_list.at(i_pf_file).Find(v_problems.at(i_problem)) > 0)

		}//for (int i_pf_file = 0; i_pf_file < v_pf_files_list.size(); i_pf_file++)


		//::Tools::vShow("a");

	}//for (int i_problem = 0; i_problem < v_problems.size(); i_problem++)


	for (int i_pf_file = 0; i_pf_file < v_pf_files_list.size(); i_pf_file++)
	{
		for (int ii = 0; ii < v_result.size(); ii++)
		{
			if  (v_result.at(ii).Find(v_pf_files_list.at(i_pf_file)) >= 0)
				::Tools::vReportInFile("zzzz_result22.txt", v_result.at(ii));
		}//for (int ii = 0; ii < v_result.size(); ii++)
	}//for (int i_pf_file = 0; i_pf_file < v_pf_files_list.size(); i_pf_file++)

	return(c_err);
}//CError  CBinaryMultiPaints::eLoadPFFromManyFileFiles()


CError  CBinaryMultiPaints::eLoadPFFromFile(FILE  *pfSource)
{
	CError  c_err(iERROR_PARENT_CBinaryMultiPaints);

	CString  s_buf;

	s_buf = ::Tools::sReadLine(pfSource);

	bool  b_loaded, b_joined_pf;
	CMultiIndividualPFpoints *pc_pf_point_buf;
	b_loaded = true;

	while ((!feof(pfSource)) && (b_loaded == true))
	{
		pc_pf_point_buf = new CMultiIndividualPFpoints();
		b_loaded = pc_pf_point_buf->bLoadFromFilePRW(pfSource);
		b_joined_pf = false;

		if (b_loaded == true)
		{
			//v_theoretical_optimal_pareto_front.push_back(pc_pf_point_buf);
			i_join_the_global_pareto(pc_pf_point_buf);
		}//if (b_loaded == true)

		delete  pc_pf_point_buf;
	}//while ((!feof(pfSource)) && (i_pareto_points_num < v_theoretical_optimal_pareto_front.size()))


	return(c_err);
};//CError  CBinaryMultiPaints::eLoadPFFromFile(FILE  *pfSource)


CError CBinaryMultiPaints::e_load_settings_from_file(FILE *pfConfig)
{	
	CError  c_err(iERROR_PARENT_CBinaryMultiPaints);
	CString  s_buf, s_line;
	CString  s_setting_name;
	double  d_measure_weight;

	CString  s_settings_type;


	//c_err = eLoadPFFromManyFileFiles();
	//::Tools::vShow("DONE");

	int  i_task_multiplier;

	i_task_multiplier = Tools::iReadLine(pfConfig, &s_setting_name);
	if (s_setting_name != s_MULTI_PAINTS_TASK_MULTIPLE)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", s_MULTI_PAINTS_TASK_MULTIPLE, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (s_setting_name != s_MULTI_PAINTS_TASK_MULTIPLE)


	s_settings_type = Tools::sReadLine(pfConfig, &s_setting_name);
	if (s_setting_name != s_MULTI_PAINTS_CONFIG_TYPE)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", s_MULTI_PAINTS_CONFIG_TYPE, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (s_setting_name != s_MULTI_PAINTS_MEASURE_MAKESPAN)


	bool  b_recognized = false;
	if (s_settings_type == s_MULTI_PAINTS_CONFIG_TYPE_HAND_MADE)
	{
		b_recognized = true;
		c_err = e_load_settings_from_file_hand_made(pfConfig);
	}//if (s_settings_type == s_MULTI_PAINTS_CONFIG_TYPE_HAND_MADE)

	if (s_settings_type == s_MULTI_PAINTS_CONFIG_TYPE_GENERATOR)
	{
		b_recognized = true;
		c_err = e_load_settings_from_file_generator(pfConfig);
	}//if (s_settings_type == s_MULTI_PAINTS_CONFIG_TYPE_GENERATOR)


	if (b_recognized == false)
	{
		s_buf.Format("Unknown multi objective paints config file type (%s)", s_settings_type);
		c_err.vSetError(iERROR_CODE_MNULTI_OBJECTIVE_PAINTS_UNKNOWN_CONFIG_FILE, s_buf);

		return(c_err);
	}//if (b_recognized == false)


	if  (i_task_multiplier > 1)  c_err = e_multiply_task(i_task_multiplier);


	//eSave("zzz_paints_gen_test.txt");
	//::Tools::vShow("aaaaa");


	//eGenerateTestCaseSerie("zzz_serie_paints_list.txt");

	//eSave("zzz_report.txt");


	return(c_err);
}//CError CBinaryMultiPaints::eLoadSettingsFromFile(FILE *pfConfig)




CError CBinaryMultiPaints::e_load_settings_from_file_hand_made(FILE *pfConfig)
{
	CError  c_err(iERROR_PARENT_CBinaryMultiPaints);
	CString  s_buf, s_line;
	CString  s_setting_name;
	double  d_measure_weight;


	d_measure_weight = Tools::dReadLine(pfConfig, &s_setting_name);
	if (s_setting_name != s_MULTI_PAINTS_MEASURE_MAKESPAN)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", s_MULTI_PAINTS_MEASURE_MAKESPAN, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (s_setting_name != s_MULTI_PAINTS_MEASURE_MAKESPAN)

	if (d_measure_weight >= 0)
	{
		CBinaryMultiPaintsMeasureMakespan  *pc_makespan;

		pc_makespan = new CBinaryMultiPaintsMeasureMakespan(this);
		pc_makespan->dWeight = d_measure_weight;
		v_measures.push_back(pc_makespan);
	}//if (d_measure_weight >= 0)



	d_measure_weight = Tools::dReadLine(pfConfig, &s_setting_name);
	if (s_setting_name != s_MULTI_PAINTS_MEASURE_OVERHEAD)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", s_MULTI_PAINTS_MEASURE_OVERHEAD, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (s_setting_name != s_MULTI_PAINTS_MEASURE_OVERHEAD)

	if (d_measure_weight >= 0)
	{
		CBinaryMultiPaintsMeasureOverhead  *pc_makespan;

		pc_makespan = new CBinaryMultiPaintsMeasureOverhead(this);
		pc_makespan->dWeight = d_measure_weight;
		v_measures.push_back(pc_makespan);
	}//if (d_measure_weight >= 0)






	int  i_paints_num;
	i_paints_num = Tools::iReadLine(pfConfig, &s_setting_name);

	if (s_setting_name != s_MULTI_PAINTS_PAINTS)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", s_MULTI_PAINTS_PAINTS, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (s_setting_name != s_MULTI_PAINTS_PAINTS)


	if (i_paints_num <= 0)
	{
		s_buf.Format("Invalid paints number (%d)", i_paints_num);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (i_paints_num <= 0)


	double  d_paint_demand;
	for (int i_paint_demand = 0; i_paint_demand < i_paints_num; i_paint_demand++)
	{
		d_paint_demand = Tools::dReadLine(pfConfig, &s_setting_name);

		v_paint_demands.push_back(d_paint_demand);
	}//for (int i_paint_demand = 0; i_paint_demand < i_paints_num; i_paint_demand++)


	i_machine_number = Tools::iReadLine(pfConfig, &s_setting_name);
	if (s_setting_name != s_MULTI_PAINTS_MACHINES)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", s_MULTI_PAINTS_MACHINES, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (s_setting_name != s_MULTI_PAINTS_MACHINES)


	if (i_machine_number <= 0)
	{
		s_buf.Format("Invalid machines number (%d)", i_machine_number);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (i_paints_num <= 0)


	int  i_recipe_number;
	i_recipe_number = Tools::iReadLine(pfConfig, &s_setting_name);
	if (s_setting_name != s_MULTI_PAINTS_RECIPES)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", s_MULTI_PAINTS_RECIPES, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (s_setting_name != s_MULTI_PAINTS_MACHINES)



	if (i_recipe_number <= 0)
	{
		s_buf.Format("Invalid recipes number (%d)", i_recipe_number);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (i_paints_num <= 0)


	CBinaryMultiPaintsRecipe  c_recipe_buf;
	for (int i_recipe = 0; i_recipe < i_recipe_number; i_recipe++)
	{
		s_line = Tools::sReadLine(pfConfig, &s_setting_name);
		c_err = c_recipe_buf.eReadFromLine(s_line, v_recipes.size());
		if (c_err)  return(c_err);
		v_recipes.push_back(c_recipe_buf);
	}//for (int i_recipe = 0; i_recipe < i_recipe_number; i_recipe++)


	/*i_bit_length = Tools::iReadLine(pfSettings, &s_setting_name);
	if (s_setting_name != COMP_PROBLEM_PARAM_BIT_LENGTH)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_BIT_LENGTH, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)



	d_max_val = Tools::dReadLine(pfSettings, &s_setting_name);
	if (s_setting_name != COMP_PROBLEM_PARAM_MAX_VAL)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", COMP_PROBLEM_PARAM_MAX_VAL, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if  (s_setting_name  !=  VGA_PARAM_GENERATIONS)*/


	eSave("paints_test_control.txt");

	return(c_err);
}//CError CBinaryMultiPaints::e_load_settings_from_file_hand_made(FILE *pfConfig)







CError CBinaryMultiPaints::e_load_settings_from_file_generator(FILE *pfConfig)
{
	CError  c_err(iERROR_PARENT_CBinaryMultiPaints);
	CString  s_buf, s_line;
	CString  s_setting_name;
	double  d_measure_weight;


	d_measure_weight = Tools::dReadLine(pfConfig, &s_setting_name);
	if (s_setting_name != s_MULTI_PAINTS_MEASURE_MAKESPAN)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", s_MULTI_PAINTS_MEASURE_MAKESPAN, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (s_setting_name != s_MULTI_PAINTS_MEASURE_MAKESPAN)

	if (d_measure_weight >= 0)
	{
		CBinaryMultiPaintsMeasureMakespan  *pc_makespan;

		pc_makespan = new CBinaryMultiPaintsMeasureMakespan(this);
		pc_makespan->dWeight = d_measure_weight;
		v_measures.push_back(pc_makespan);
	}//if (d_measure_weight >= 0)



	d_measure_weight = Tools::dReadLine(pfConfig, &s_setting_name);
	if (s_setting_name != s_MULTI_PAINTS_MEASURE_OVERHEAD)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", s_MULTI_PAINTS_MEASURE_OVERHEAD, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (s_setting_name != s_MULTI_PAINTS_MEASURE_OVERHEAD)

	if (d_measure_weight >= 0)
	{
		CBinaryMultiPaintsMeasureOverhead  *pc_makespan;

		pc_makespan = new CBinaryMultiPaintsMeasureOverhead(this);
		pc_makespan->dWeight = d_measure_weight;
		v_measures.push_back(pc_makespan);
	}//if (d_measure_weight >= 0)






	int  i_paints_num;
	i_paints_num = Tools::iReadLine(pfConfig, &s_setting_name);

	if (s_setting_name != s_MULTI_PAINTS_PAINTS)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", s_MULTI_PAINTS_PAINTS, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (s_setting_name != s_MULTI_PAINTS_PAINTS)


	if (i_paints_num <= 0)
	{
		s_buf.Format("Invalid paints number (%d)", i_paints_num);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (i_paints_num <= 0)



	i_machine_number = Tools::iReadLine(pfConfig, &s_setting_name);
	if (s_setting_name != s_MULTI_PAINTS_MACHINES)
	{
		s_buf.Format("Error at reading (%s) parameter read:(%s)", s_MULTI_PAINTS_MACHINES, s_setting_name);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (s_setting_name != s_MULTI_PAINTS_MACHINES)


	if (i_machine_number <= 0)
	{
		s_buf.Format("Invalid machines number (%d)", i_machine_number);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (i_paints_num <= 0)




	
	bool  b_line_for_all;
	b_line_for_all = false;
	int  i_line_index;
	CString  s_to_produce, s_base_amount, s_eff;
	CBinaryMultiFuncTool  c_func_buf;
	CBinaryMultiFuncAmountTool  c_func_amount_buf;
	CBinaryMultiPaintsRecipe  c_recipe_buf;
	int  i_recipe_offset;
	vector<double>  v_amount_per_machine;

	i_recipe_offset = 0;
	for (int i_paint = 0; i_paint < i_paints_num; i_paint++)
	{
		if (b_line_for_all == false)
		{
			s_line = ::Tools::sReadLine(pfConfig, &s_setting_name);
			if (s_setting_name.Find(s_MULTI_PAINTS_CONFIG_ALL_PAINTS_CONFIG) == 0)  b_line_for_all = true;

			i_line_index = 0;
			s_to_produce = ::Tools::sExtractFromString(s_line, i_line_index, &i_line_index);
			s_base_amount = ::Tools::sExtractFromString(s_line, i_line_index, &i_line_index);
			s_eff = ::Tools::sExtractFromString(s_line, i_line_index, &i_line_index);

		}//if (b_line_for_all == false)


		//c_err = c_func_buf.eConfig(s_base_amount);
		//if (c_err)  return(c_err);

		//c_err = c_func_amount_buf.eConfig(s_to_produce);
		//if (c_err)  return(c_err);

		//::Tools::vShow(c_func_buf.sToString());
		//::Tools::vShow(c_func_amount_buf.sToString());

		v_amount_per_machine.clear();

		for (int i_mach = 0; i_mach < i_machine_number; i_mach++)
		{
			c_recipe_buf.v_clear();

			c_recipe_buf.i_offset = i_recipe_offset;
			i_recipe_offset++;

			c_recipe_buf.i_paint_produced_offset = i_paint;
			c_recipe_buf.v_machines.push_back(i_mach);
			c_recipe_buf.d_duration = 1;

			c_err = c_func_buf.eConfig(s_base_amount);
			if (c_err)  return(c_err);

			c_recipe_buf.d_produced_amount = c_func_buf.dGetValue(i_mach + 1);
			c_recipe_buf.d_duration = c_recipe_buf.d_duration * c_recipe_buf.d_produced_amount;
			v_amount_per_machine.push_back(c_recipe_buf.d_produced_amount);

			if (i_mach > 0)
			{
				c_err = c_func_buf.eConfig(s_eff);
				if (c_err)  return(c_err);

				c_recipe_buf.d_duration = c_recipe_buf.d_duration / c_func_buf.dGetValue(i_mach);
			}//if (i_mach > 0)

			v_recipes.push_back(c_recipe_buf);
		}//for (int i_mach = 0; i_mach < i_machine_number; i_mach++)


		c_err = c_func_amount_buf.eConfig(s_to_produce);
		if (c_err)  return(c_err);

		v_paint_demands.push_back(c_func_amount_buf.dGetValue(&v_amount_per_machine));
		

	
	}//for (int i_paint = 0; i_paint < i_paints_num; i_paint++)
	





	return(c_err);
}//CError CBinaryMultiPaints::e_load_settings_from_file_generator(FILE *pfConfig)




CError CBinaryMultiPaints::e_multiply_task(int  iMultiplier)
{
	CError  c_err(iERROR_PARENT_CBinaryMultiPaints);

	int  i_orig_paint_dem_size;
	i_orig_paint_dem_size = v_paint_demands.size();

	for (int i_replication = 1; i_replication < iMultiplier; i_replication++)
	{
		for (int ii = 0; ii < i_orig_paint_dem_size; ii++)
			v_paint_demands.push_back(v_paint_demands.at(ii));
	}//for (int i_replication = 0; i_replication < iMultiplier; i_replication++)


	int  i_orig_machine_numer;
	i_orig_machine_numer = i_machine_number;
	i_machine_number *= iMultiplier;


	int  i_orig_recipe_size;
	i_orig_recipe_size = v_recipes.size();

	CBinaryMultiPaintsRecipe  c_recipe_buf;


	for (int i_replication = 1; i_replication < iMultiplier; i_replication++)
	{
		for (int i_recipe = 0; i_recipe < i_orig_recipe_size; i_recipe++)
		{
			c_recipe_buf = v_recipes.at(i_recipe);

			for (int i_mach_offset = 0; i_mach_offset < c_recipe_buf.v_machines.size(); i_mach_offset++)
				c_recipe_buf.v_machines.at(i_mach_offset) += i_replication * i_orig_machine_numer;

			c_recipe_buf.i_paint_produced_offset += i_replication * i_orig_paint_dem_size;

			v_recipes.push_back(c_recipe_buf);
		}//for (int i_recipe = 0; i_recipe < i_orig_recipe_size; i_recipe++)
			
	}//for (int i_replication = 0; i_replication < iMultiplier; i_replication++)

	


	return(c_err);
}//CError CBinaryMultiPaints::e_multiply_task(int  iMultiplier)



CError  CBinaryMultiPaints::eGenerateTestCaseSerie(CString  sFileNameList)
{
	CError  c_err(iERROR_PARENT_CBinaryMultiPaints);

	vector<int>  v_machines_num;
	vector<int>  v_paints_num;

	vector<CString>  v_amount_type;
	vector<double>  v_amount_mud;

	vector<CString>  v_base_type;
	vector<double>  v_base_multi;

	vector<CString>  v_eff_type;
	vector<double>  v_eff_mud;

	FILE  *pf_filename_list;

	pf_filename_list = fopen(sFileNameList, "w+");
	if (pf_filename_list == NULL)
	{
		c_err.vSetError("CError  CBinaryMultiPaints::eGenerateTestCaseSerie(CString  sFileNameList)  - Can not open FileNameList File");
		return(c_err);
	}//if (pf_filename_list == NULL)


	v_machines_num.push_back(4);
	v_machines_num.push_back(8);
	v_machines_num.push_back(12);

	v_paints_num.push_back(4);
	v_paints_num.push_back(4);
	v_paints_num.push_back(8);

	v_amount_type.push_back(s_MULTI_PAINTS_CONFIG_AMOUNT_FUNC_TYPE_SUM);
	v_amount_mud.push_back(1.5);
	v_amount_type.push_back(s_MULTI_PAINTS_CONFIG_AMOUNT_FUNC_TYPE_LARGEST);
	v_amount_mud.push_back(2.98);


	v_base_type.push_back(s_MULTI_PAINTS_CONFIG_FUNC_TYPE_EXP);
	v_base_multi.push_back(3);
	v_base_type.push_back(s_MULTI_PAINTS_CONFIG_FUNC_TYPE_EXP);
	v_base_multi.push_back(2);
	v_base_type.push_back(s_MULTI_PAINTS_CONFIG_FUNC_TYPE_POW);
	v_base_multi.push_back(2);
	v_base_type.push_back(s_MULTI_PAINTS_CONFIG_FUNC_TYPE_LINEAR);
	v_base_multi.push_back(2);

	v_eff_type.push_back(s_MULTI_PAINTS_CONFIG_FUNC_TYPE_EXP);
	v_eff_mud.push_back(2);
	v_eff_type.push_back(s_MULTI_PAINTS_CONFIG_FUNC_TYPE_POW);
	v_eff_mud.push_back(2);
	v_eff_type.push_back(s_MULTI_PAINTS_CONFIG_FUNC_TYPE_LINEAR);
	v_eff_mud.push_back(2);
	v_eff_type.push_back(s_MULTI_PAINTS_CONFIG_FUNC_TYPE_LINEAR);
	v_eff_mud.push_back(2);

	



	for (int i_mach = 0; i_mach < v_machines_num.size(); i_mach++)
	{
		int i_paint = i_mach;
		//for (int i_paint = 0; i_paint < v_paints_num.size(); i_paint++)
		{
			for (int i_amount = 0; i_amount < v_amount_type.size(); i_amount++)
			{
				for (int i_base = 0; i_base < v_base_type.size(); i_base++)
				{
					int i_eff = i_base;
					//for (int i_eff = 0; i_eff < v_eff_type.size(); i_eff++)
					{
						c_err = eGenerateTestCase
							(
								pf_filename_list,
								v_paints_num.at(i_paint), v_machines_num.at(i_mach),
								v_amount_type.at(i_amount), v_amount_mud.at(i_amount),
								v_base_type.at(i_base), v_base_multi.at(i_base),
								v_eff_type.at(i_eff), v_eff_mud.at(i_eff)
							);
					}//for (int i_eff = 0; i_eff < v_eff_type.size(); i_eff++)
				}//for (int i_amount = 0; i_amount < v_amount_type.size(); i_amount++)
			}//for (int i_sum = 0; i_sum < v_sum_type.size(); i_sum++)
		}//for (int i_paint = 0; i_paint < v_paints_num.size(); i_paint++)
	}//for (int i_mach = 0; i_mach < v_machines_num.size(); i_mach++)

	fclose(pf_filename_list);

	return(c_err);
}//CError  CBinaryMultiPaints::eGenerateTestCaseSerie(CString  sFileNameList)




CError  CBinaryMultiPaints::eGenerateTestCase
	(
		FILE  *pfFileNameList,
		int  iPaints, int iMachines,
		CString sAmountType, double dAmountMulti,
		CString sBaseType, double  dBaseVal,
		CString sEffType, double  dEffVal
	)
{
	CError  c_err;
	CString  s_buf;
	CString  s_file_name;

	CString s_amount_file;
	CString s_base_file;
	CString s_eff_type;
	
	s_amount_file = sAmountType;
	s_base_file = sBaseType;
	s_eff_type = sEffType;

	s_amount_file.Replace("*", "");
	s_base_file.Replace("*", "");
	s_base_file.Replace("^", "");
	s_eff_type.Replace("*", "");
	s_eff_type.Replace("^", "");

	s_buf.Format("%.0lf", dAmountMulti*10);
	s_amount_file += s_buf;

	s_buf.Format("%.0lf", dBaseVal * 10);
	s_base_file += s_buf;

	s_buf.Format("%.0lf", dEffVal * 10);
	s_eff_type += s_buf;

	s_file_name.Format("paint_%.2d_%.2d_%s_%s_%s", iPaints, iMachines, s_amount_file, s_base_file, s_eff_type);


	FILE  *pf_dest;
	pf_dest = fopen(s_file_name + ".txt", "w+");
	if (pf_dest == NULL)
	{
		s_buf.Format("Can not create file (%s)", s_file_name + ".txt");
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (pf_dest == NULL)
	
	if (pfFileNameList != NULL)  fprintf(pfFileNameList, s_file_name + "\n");



	fprintf(pf_dest, "PaintsGenerator\\\\PaintsConfigType\n");

	fprintf(pf_dest, "0.5\\\\%s\n", s_MULTI_PAINTS_MEASURE_MAKESPAN);
	fprintf(pf_dest, "0.5\\\\%s\n", s_MULTI_PAINTS_MEASURE_OVERHEAD);

	fprintf(pf_dest, "%d\\\\%s\n", iPaints, s_MULTI_PAINTS_PAINTS);
	fprintf(pf_dest, "%d\\\\%s\n", iMachines, s_MULTI_PAINTS_MACHINES);

	fprintf
		(
			pf_dest, "%s%lf %s%lf %s%lf\\\\%s      amount_to_produce; mach_base_amount; mach_eff\n", 
			sAmountType, dAmountMulti,
			sBaseType, dBaseVal,
			sEffType, dEffVal,
			s_MULTI_PAINTS_CONFIG_ALL_PAINTS_CONFIG
		);

	
	

	fclose(pf_dest);


	return(c_err);
}//CError  CBinaryMultiPaints::eGenerateTestCase



CError  CBinaryMultiPaints::eSave(CString  sDest)
{
	CError  c_err;
	CString  s_buf;
	FILE  *pf_dest;

	pf_dest = fopen(sDest, "w+");

	if (pf_dest == NULL)
	{
		s_buf.Format("Can not open destination file (CError  CBinaryMultiPaints::eSave(CString  sDest)) (%s)", sDest);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (pf_dest == NULL)

	c_err = eSave(pf_dest);

	fclose(pf_dest);

	return(c_err);
}//void  CBinaryMultiPaints::vSave(CString  sDest)




CError  CBinaryMultiPaints::eSave(FILE  *pfDest)
{
	CError  c_err;

	CString  s_line, s_buf;

	s_buf.Format("%d\\\\%s\n", v_paint_demands.size(), s_MULTI_PAINTS_PAINTS);
	fprintf(pfDest, s_buf);


	for (int i_paint_demand = 0; i_paint_demand < v_paint_demands.size(); i_paint_demand++)
	{
		s_buf.Format("%.2lf\\\\%d\n", v_paint_demands.at(i_paint_demand), i_paint_demand);
		fprintf(pfDest, s_buf);
	}//for (int i_paint_demand = 0; i_paint_demand < i_paints_num; i_paint_demand++)

	s_buf.Format("%d\\\\%s\n", i_machine_number, s_MULTI_PAINTS_MACHINES);
	fprintf(pfDest, s_buf);


	s_buf.Format("%d\\\\%s\n", v_recipes.size(), s_MULTI_PAINTS_RECIPES);
	fprintf(pfDest, s_buf);


	for (int i_recipe = 0; i_recipe < v_recipes.size(); i_recipe++)
	{
		s_buf.Format("%s\\\\%d\n", v_recipes.at(i_recipe).sGetRecipe(), i_recipe);
		fprintf(pfDest, s_buf);
	}//for (int i_paint_demand = 0; i_paint_demand < i_paints_num; i_paint_demand++)





	return(c_err);
}//void  CBinaryMultiPaints::vSave(FILE  *pfDest)


void  CBinaryMultiPaintsRecipe::v_clear()
{
	i_offset = -1;
	v_machines.clear();
	d_duration = 1;
	i_paint_produced_offset = -1;
	d_produced_amount = 0;
	//d_product = 1;
}//void  CBinaryMultiPaintsRecipe::v_clear()

CError  CBinaryMultiPaintsRecipe::eReadFromLine(CString  sLine, int iOffset)
{
	CError  c_err;

	CString  s_buf;
	int  i_index;
	CString  s_driving_sequence;

	v_clear();

	i_offset = iOffset;
	i_index = 0;

	s_driving_sequence = ::Tools::sExtractFromString(sLine, i_index, &i_index);
	if (s_driving_sequence != s_MULTI_PAINTS_RECIPE_MACHINES_START)
	{
		s_buf.Format("Searching for machines start in recipe offset (%d) found(%s) expected(%s)", i_offset, s_driving_sequence, s_MULTI_PAINTS_RECIPE_MACHINES_START);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (s_buf != s_MULTI_PAINTS_RECIPE_MACHINES_START)


	int  i_machines_num;
	i_machines_num = ::Tools::iExtractFromString(sLine, i_index, &i_index);
	if (i_machines_num < 0)
	{
		s_buf.Format("Invalid machine offset (%d)", i_machines_num);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (i_machines_num <= 0)


	int  i_machine_offset;
	for (int i_machine = 0; i_machine < i_machines_num; i_machine++)
	{
		i_machine_offset = ::Tools::iExtractFromString(sLine, i_index, &i_index);
		v_machines.push_back(i_machine_offset);
	}//for (int i_machine = 0; i_machine < i_machines_num; i_machine++)


	s_driving_sequence = ::Tools::sExtractFromString(sLine, i_index, &i_index);
	if (s_driving_sequence != s_MULTI_PAINTS_RECIPE_MACHINES_END)
	{
		s_buf.Format("Searching for machines end in recipe offset (%d) found(%s) expected(%s)", i_offset, s_driving_sequence, s_MULTI_PAINTS_RECIPE_MACHINES_END);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (s_buf != s_MULTI_PAINTS_RECIPE_MACHINES_START)


	d_duration = ::Tools::dExtractFromString(sLine, i_index, &i_index);
	i_paint_produced_offset = ::Tools::iExtractFromString(sLine, i_index, &i_index);
	d_produced_amount = ::Tools::dExtractFromString(sLine, i_index, &i_index);


	//int  iReadLine(FILE  *pfSource, CString  *psComment = NULL);
	//double  dReadLine(FILE  *pfSource, CString  *psComment = NULL);
	//CString  sExtractFromString(CString  sLineToSearch, int  iIndex /*= 0*/, int  *piFinishIndex = NULL);
	//int  iExtractFromString(CString  sLineToSearch, int  iIndex /*= 0*/, int  *piFinishIndex = NULL);
	//double  dExtractFromString(CString  sLineToSearch, int  iIndex, int  *piFinishIndex = NULL);




	return(c_err);
}//CEror  CBinaryMultiPaintsRecipe::eReadFromLine(CString  sLine)


CString  CBinaryMultiPaintsRecipe::sGetRecipe()
{
	CString  s_buf;
	CString  s_recipe;

	s_recipe = s_MULTI_PAINTS_RECIPE_MACHINES_START;

	s_buf.Format(" %d", v_machines.size());
	s_recipe += s_buf;

	for (int ii = 0; ii < v_machines.size(); ii++)
	{
		s_buf.Format(" %d", v_machines.at(ii));
		s_recipe += s_buf;
	}//for (int ii = 0; ii < v_machines.size(); ii++)

	s_recipe += " ";
	s_recipe += s_MULTI_PAINTS_RECIPE_MACHINES_END;

	s_buf.Format(" %.2lf %d %.2lf", d_duration, i_paint_produced_offset, d_produced_amount);
	s_recipe += s_buf;


	return(s_recipe);
}//CString  CBinaryMultiPaintsRecipe::sGetRecipe()



CError CBinaryMultiPaints::eConfigure(istream *psSettings)
{
	CError  c_err;

	c_err = CBinaryFileConfigEvaluation::eConfigure(psSettings);
	if (c_err)  return(c_err);

	c_err = CBinaryMultiObjectiveProblem::eConfigure(psSettings);
	if (c_err)  return(c_err);


	return(c_err);
}//CError CBinaryMultiPaints::eConfigure(istream *psSettings)


CError CBinaryMultiPaints::eReport(FILE *pfReport)
{
	CError  c_err;

	fprintf(pfReport, "dupa dupa");
	//return pc_deceptive_concatenation_function->eCreateReport(pfReport);

	return(c_err);
}//CError CBinaryDeceptiveConcatenationEvaluation::eReport(FILE *pfReport)




//---------------------------------------------------CBinaryMultiReversing---------------------------------------------------
double CBinaryMultiReversingMeasure1s::dEvaluate(CBinaryCoding *pcFenotype)
{
	return(((CBinaryMultiReversing*)pc_parent_problem)->pc_evaluation_inner->dEvaluate(pcFenotype));
}//double CBinaryMultiReversingMeasure1s::dEvaluate(CBinaryCoding *pcFenotype)




CBinaryMultiReversingMeasure0s::CBinaryMultiReversingMeasure0s(CBinaryMultiObjectiveProblem  *pcParentProblem, double  dMax) : CMultiObjectiveMeasure(pcParentProblem)
{
	pc_phenotype_inner = NULL;
	s_name = s_MULTI_REVERSING_MEASURE_1s;
	d_max = dMax;
}//CBinaryMultiReversingMeasure0s::CBinaryMultiReversingMeasure0s(CBinaryMultiObjectiveProblem  *pcParentProblem) : CMultiObjectiveMeasure(pcParentProblem)


CBinaryMultiReversingMeasure0s::CBinaryMultiReversingMeasure0s(const CBinaryMultiReversingMeasure0s &pcOther) : CMultiObjectiveMeasure(pcOther)
{
	b_reverse_bit_order_for_0_measure = pcOther.b_reverse_bit_order_for_0_measure;

	if (pcOther.pc_phenotype_inner == NULL)
		pc_phenotype_inner = NULL;
	else
		pc_phenotype_inner = new CBinaryCoding(pcOther.pc_phenotype_inner);
}//CBinaryMultiReversingMeasure0s::CBinaryMultiReversingMeasure0s(const CBinaryMultiReversingMeasure0s &pcOther)

CBinaryMultiReversingMeasure0s::~CBinaryMultiReversingMeasure0s()
{
	if (pc_phenotype_inner != NULL)  delete  pc_phenotype_inner;
}//CBinaryMultiReversingMeasure0s::~CBinaryMultiReversingMeasure0s()



double CBinaryMultiReversingMeasure0s::dEvaluate(CBinaryCoding *pcFenotype)
{
	if (pc_phenotype_inner == NULL)  pc_phenotype_inner = new CBinaryCoding(pcFenotype);

	if (b_reverse_bit_order_for_0_measure == false)
	{
		for (int ii = 0; ii < pcFenotype->iGetNumberOfBits(); ii++)
			pc_phenotype_inner->piGetBits()[ii] = pcFenotype->piGetBits()[ii] * (-1) + 1;
	}//if (b_reverse_bit_order_for_0_measure == false)
	else
	{
		for (int ii = 0; ii < pcFenotype->iGetNumberOfBits(); ii++)
			pc_phenotype_inner->piGetBits()[ii] = pcFenotype->piGetBits()[pcFenotype->iGetNumberOfBits() - ii - 1] * (-1) + 1;
	}//else  if (b_reverse_bit_order_for_0_measure == false)


	return(((CBinaryMultiReversing*)pc_parent_problem)->pc_evaluation_inner->dEvaluate(pc_phenotype_inner));
}//double CBinaryMultiReversingMeasure0s::dEvaluate(CBinaryCoding *pcFenotype)



CBinaryMultiReversing::CBinaryMultiReversing() : CBinaryMultiObjectiveProblem::CBinaryMultiObjectiveProblem()
{
	
}//CBinaryMultiReversing::CBinaryMultiReversing()

CBinaryMultiReversing::CBinaryMultiReversing(FILE *pfConfig, CError *pcError)
{
	*pcError = e_init(pfConfig);
}//CBinaryMultiPaints::CBinaryMultiPaints(FILE *pfConfig, CError *pcError)


CBinaryMultiReversing::CBinaryMultiReversing(const CBinaryMultiReversing &pcOther) : CBinaryMultiObjectiveProblem::CBinaryMultiObjectiveProblem(pcOther)
{
	::Tools::vShow("CBinaryMultiReversing::CBinaryMultiReversing(const CBinaryMultiReversing &pcOther)   CA NOT COPY INNER PROBLEM");
	//if (pcOther.pc_evaluation_inner != NULL)
		//pc_evaluation_inner = new pcOther.pc_evaluation_inner.
	//else
		//pc_evaluation_inner = NULL;
}//CBinaryMultiPaints::CBinaryMultiPaints(const CBinaryMultiPaints &pcOther)

CBinaryMultiReversing::~CBinaryMultiReversing()
{
	//delete pc_deceptive_concatenation_function;
}//CBinaryMultiPaints::~CBinaryMultiPaints()


CError CBinaryMultiReversing::eConfigure(istream *psSettings)
{
	CError  c_err;
	   

	c_err  =  CBinaryMultiObjectiveProblem::eConfigure(psSettings);
	if (c_err)  return(c_err);


	CUIntCommandParam p_reversed_it_order_for_0_measure(EVALUATION_ARGUMENT_REVERSE_BIT_ORDER);
	int i_reverse_bit_order_for_0_measure;
	i_reverse_bit_order_for_0_measure = p_reversed_it_order_for_0_measure.iGetValue(psSettings, &c_err);
	if (c_err)  return(c_err);
	

	CUIntCommandParam p_optimal_pf_bit_step(EVALUATION_ARGUMENT_BINARY_MULTI_OPTIMAL_PF_BIT_STEP);
	int  i_create_optimal_pf_bit_step;
	i_create_optimal_pf_bit_step = p_optimal_pf_bit_step.iGetValue(psSettings, &c_err);
	if (c_err)  return(c_err);

	

	pc_evaluation_inner = EvaluationUtils::pcGetEvaluation<CBinaryCoding>(psSettings, &c_err);
	if (c_err)  return(c_err);
	

	i_number_of_elements = pc_evaluation_inner->iGetNumberOfElements();
	d_max_value = pc_evaluation_inner->dGetMaxValue();


	//this MUST be after pc_evaluation_inner creation
	CBinaryMultiReversingMeasure1s  *pc_measure_1s;
	pc_measure_1s = new CBinaryMultiReversingMeasure1s(this, pc_evaluation_inner->dGetMaxValue());
	pc_measure_1s->dWeight = 0.5;
	v_measures.push_back(pc_measure_1s);

	CBinaryMultiReversingMeasure0s  *pc_measure_0s;
	pc_measure_0s = new CBinaryMultiReversingMeasure0s(this, pc_evaluation_inner->dGetMaxValue());
	pc_measure_0s->dWeight = 0.5;
	v_measures.push_back(pc_measure_0s);
	pc_measure_0s->vSetReversingBitOrder(i_reverse_bit_order_for_0_measure);


	//this MUST be after measures
	vFillOptimalPf(i_create_optimal_pf_bit_step);
		

	return(c_err);
}//virtual CError CBinaryMultiReversing::eConfigure(istream *psSettings)



void  CBinaryMultiReversing::vFillOptimalPf(int  iCreateOptimalPfBitStep)
{
	if (iCreateOptimalPfBitStep <= 0)  return;


	CBinaryCoding c_genotype(i_number_of_elements);

	for (int ii = 0; ii < i_number_of_elements; ii++)
		c_genotype.piGetBits()[ii] = 0;
	   	

	CMultiIndividualPFpoints  *pc_pf_point_buf;

	pc_pf_point_buf = new CMultiIndividualPFpoints();
	pc_pf_point_buf->vConfigure(&c_genotype, true, &v_pf_buffer, false, this);
	pc_pf_point_buf->vRate();
	v_theoretical_optimal_pareto_front.push_back(pc_pf_point_buf);

	int  i_cur_start_bit;
	i_cur_start_bit = 0;

	while (i_cur_start_bit < i_number_of_elements)
	{
		for (int ii = i_cur_start_bit; (ii < i_cur_start_bit + iCreateOptimalPfBitStep)&&(ii < i_number_of_elements); ii++)
			c_genotype.piGetBits()[ii] = 1;

		pc_pf_point_buf = new CMultiIndividualPFpoints();
		pc_pf_point_buf->vConfigure(&c_genotype, true, &v_pf_buffer, false, this);
		pc_pf_point_buf->vRate();
		v_theoretical_optimal_pareto_front.push_back(pc_pf_point_buf);

		i_cur_start_bit += iCreateOptimalPfBitStep;
	}//while (i_cur_start_bit < i_number_of_elements)


	/*::Tools::vReportInFile("zzz_opt_pf.txt", "", true);
	for (int ii = 0; ii < v_theoretical_optimal_pareto_front.size(); ii++)
	{
		::Tools::vReportInFile("zzz_opt_pf.txt", v_theoretical_optimal_pareto_front.at(ii)->sReport());		
	}//for (int ii = 0; ii < v_theoretical_optimal_pareto_front.size(); ii++)*/


	vFlushGlobalPareto();//creating and rating v_theoretical_optimal_pareto_front, makes object storing the globaly best solutions
	vSetZeroFFE();//no ffe used by a method yet

	//::Tools::vShow(0);
}//void  CBinaryMultiReversing::vFillOptimalPf(int  iCreateOptimalPfBitStep)




CError CBinaryMultiReversing::e_init(FILE *pfConfig)
{
	CError  c_err;

	::Tools::vShow("I should not be here!    (CError CBinaryMultiReversing::e_init(FILE *pfConfig))");
	
	return(c_err);
};//CError CBinaryMultiReversing::e_init(FILE *pfConfig)




//---------------------------------------------------CBinaryMultiMaxcutMoGomeaMeasure---------------------------------------------------
CBinaryMultiMaxcutMoGomeaMeasure::CBinaryMultiMaxcutMoGomeaMeasure(CBinaryMultiObjectiveProblem  *pcParentProblem, double  dMax) : CMultiObjectiveMeasure(pcParentProblem)
{
	s_name = s_MULTI_MAXCUT_MEASURE;
	d_max = dMax;
}//CBinaryMultiMaxcutMoGomeaMeasure::CBinaryMultiMaxcutMoGomeaMeasure(CBinaryMultiObjectiveProblem  *pcParentProblem, double  dMax) : CMultiObjectiveMeasure(pcParentProblem)


CBinaryMultiMaxcutMoGomeaMeasure::~CBinaryMultiMaxcutMoGomeaMeasure()
{

}//CBinaryMultiMaxcutMoGomeaMeasure::~CBinaryMultiMaxcutMoGomeaMeasure()


double  CBinaryMultiMaxcutMoGomeaMeasure::dEvaluate(CBinaryCoding *pcFenotype)
{
	double  d_result;

	d_result = 0;
	for (int i_edge = 0; i_edge < v_edges.size(); i_edge++)
	{
		if (pcFenotype->piGetBits()[v_edges.at(i_edge).i_start_node] != pcFenotype->piGetBits()[v_edges.at(i_edge).i_end_node])
			d_result += v_edges.at(i_edge).d_weight;
	}//for (int i_edge = 0; i_edge < v_edges.size(); i_edge++)


	return(d_result);
}//double  CBinaryMultiMaxcutMoGomeaMeasure::dEvaluate(CBinaryCoding *pcFenotype)


CError  CBinaryMultiMaxcutMoGomeaMeasure::eLoadFromFile(FILE  *pfSource)
{
	CError  c_err;

	CString  s_buf;
	int  i_index;
	int  i_nodes, i_links;
	CString  s_line;
	s_line = Tools::sReadLine(pfSource);

	i_index = 0;
	i_nodes = ::Tools::iExtractFromString(s_line, i_index, &i_index);

	if (pc_parent_problem->iGetNumberOfElements() != i_nodes)
	{
		s_buf.Format("CError  CBinaryMultiMaxcutMoGomeaMeasure::eLoadFromFile(FILE  *pfSource)      if (pc_parent_problem->iGetNumberOfElements() != i_nodes)");
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (pc_parent_problem->iGetNumberOfElements() != i_nodes)
	
	i_links = ::Tools::iExtractFromString(s_line, i_index, &i_index);


	//s_buf.Format("%d  %d", i_nodes, i_links);
	//::Tools::vShow(s_buf);

	CBinaryMultiMaxcutMoGomeaMeasureEdge  c_buf;

	while ((!feof(pfSource)) && (v_edges.size() < i_links))
	{
		s_line = Tools::sReadLine(pfSource);
		c_buf.eLoad(s_line);
		v_edges.push_back(c_buf);
	}//while ((!feof(pfSource)) && (v_edges.size() < i_links))



	if (v_edges.size() != i_links)
	{
		s_buf.Format("CError  CBinaryMultiMaxcutMoGomeaMeasure::eLoadFromFile(FILE  *pfSource)      if (v_edges.size() != i_links)");
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (pc_parent_problem->iGetNumberOfElements() != i_nodes)



	d_max = 0;
	for (int ii = 0; ii < v_edges.size(); ii++)
		d_max += v_edges.at(ii).d_weight;


	return(c_err);
};//CError  CBinaryMultiMaxcutMoGomeaMeasure::eLoadFromFile(FILE  *pfSource)



CError  CBinaryMultiMaxcutMoGomeaMeasure::eSaveToFile(FILE  *pfSource)
{
	CError  c_err;

	CString  s_buf;
	int  i_index;
	CString  s_line;
	s_line.Format("%d %d\n", pc_parent_problem->iGetNumberOfElements(), v_edges.size());
	fprintf(pfSource, s_line);
	
	for (int ii = 0; ii < v_edges.size(); ii++)
	{
		fprintf(pfSource, v_edges.at(ii).sToString() + "\n");
	}//for (int ii = 0; ii < v_edges.size(); ii++)


	return(c_err);
}//CError  CBinaryMultiMaxcutMoGomeaMeasure::eSaveToFile(FILE  *pfSource)



CBinaryMultiMaxcutMoGomeaMeasureEdge::CBinaryMultiMaxcutMoGomeaMeasureEdge()
{
	i_start_node = -1;
	i_end_node = -1;
	d_weight = 0;
}//CBinaryMultiMaxcutMoGomeaMeasureEdge::CBinaryMultiMaxcutMoGomeaMeasureEdge()



CError  CBinaryMultiMaxcutMoGomeaMeasureEdge::eLoad(CString  sSource)
{
	CError  c_err;

	int  i_index;
	i_index = 0;

	i_start_node = ::Tools::iExtractFromString(sSource, i_index, &i_index);
	i_start_node--;//nodes are enumerated from 1
	i_end_node = ::Tools::iExtractFromString(sSource, i_index, &i_index);
	i_end_node--;//nodes are enumerated from 1
	d_weight = ::Tools::dExtractFromString(sSource, i_index, &i_index);

	return(c_err);
}//CError  CBinaryMultiMaxcutMoGomeaMeasureEdge::eLoad(CString  sSource)


CString  CBinaryMultiMaxcutMoGomeaMeasureEdge::sToString()
{
	CString  s_result;

	s_result.Format("%d %d %lf", i_start_node, i_end_node, d_weight);

	return(s_result);
}//CString  CBinaryMultiMaxcutMoGomeaMeasureEdge::sToString()



//---------------------------------------------------CBinaryMultiMaxcutMoGomea---------------------------------------------------

CBinaryMultiMaxcutMoGomea::CBinaryMultiMaxcutMoGomea() : CBinaryMultiObjectiveProblem::CBinaryMultiObjectiveProblem()
{

}//CBinaryMultiMaxcutMoGomea::CBinaryMultiMaxcutMoGomea()

CBinaryMultiMaxcutMoGomea::CBinaryMultiMaxcutMoGomea(FILE *pfConfig, CError *pcError)
{
	*pcError = e_init(pfConfig);
}//CBinaryMultiMaxcutMoGomea::CBinaryMultiMaxcutMoGomea(FILE *pfConfig, CError *pcError)


CBinaryMultiMaxcutMoGomea::CBinaryMultiMaxcutMoGomea(const CBinaryMultiMaxcutMoGomea &pcOther) : CBinaryMultiObjectiveProblem::CBinaryMultiObjectiveProblem(pcOther)
{
	::Tools::vShow("CBinaryMultiMaxcutMoGomea::CBinaryMultiMaxcutMoGomea(const CBinaryMultiMaxcutMoGomea &pcOther)   CAN NOT COPY INNER PROBLEM");
	//if (pcOther.pc_evaluation_inner != NULL)
		//pc_evaluation_inner = new pcOther.pc_evaluation_inner.
	//else
		//pc_evaluation_inner = NULL;
}//CBinaryMultiMaxcutMoGomea::CBinaryMultiMaxcutMoGomea(const CBinaryMultiMaxcutMoGomea &pcOther) : CBinaryMultiObjectiveProblem::CBinaryMultiObjectiveProblem(pcOther)

CBinaryMultiMaxcutMoGomea::~CBinaryMultiMaxcutMoGomea()
{
	//delete pc_deceptive_concatenation_function;
}//CBinaryMultiMaxcutMoGomea::~CBinaryMultiMaxcutMoGomea()



CError  CBinaryMultiMaxcutMoGomea::eLoadPFFromFile(FILE  *pfSource)
{
	CError  c_err;

	CMultiIndividualPFpoints  *pc_pf_point_buf;

	CString  s_buf;
	int  i_pareto_points_num;

	i_pareto_points_num = ::Tools::iReadLine(pfSource);

	//::Tools::vShow(i_pareto_points_num);

	while ((!feof(pfSource)) && (v_theoretical_optimal_pareto_front.size() < i_pareto_points_num))
	{
		pc_pf_point_buf = new CMultiIndividualPFpoints();
		pc_pf_point_buf->bLoadFromFile_MaxCutMoGomea(pfSource);

		v_theoretical_optimal_pareto_front.push_back(pc_pf_point_buf);
	}//while ((!feof(pfSource)) && (i_pareto_points_num < v_theoretical_optimal_pareto_front.size()))


	return(c_err);
}//CError  CBinaryMultiMaxcutMoGomea::eLoadPFFromFile(FILE  *pfSource)



CError CBinaryMultiMaxcutMoGomea::eConfigure(istream *psSettings)
{
	CError  c_err;


	c_err = CBinaryMultiObjectiveProblem::eConfigure(psSettings);
	if (c_err)  return(c_err);

	
	CUIntCommandParam p_config_arg_num(EVALUATION_ARGUMENT_BINARY_MULTI_MAXCUT_MOGOMEA_ARG_NUM);
	
	i_number_of_elements = p_config_arg_num.iGetValue(psSettings, &c_err);
	//d_max_value = pc_evaluation_inner->dGetMaxValue();
	d_max_value = 1;
	

	CBinaryMultiMaxcutMoGomeaMeasure  *pc_measure;

	for (int i_measure = 0; i_measure < 2; i_measure++)
	{
		pc_measure = new CBinaryMultiMaxcutMoGomeaMeasure(this, 1);
		pc_measure->dWeight = 0.5;
		v_measures.push_back(pc_measure);

		CString  s_config_file_path;
		s_config_file_path.Format(EVALUATION_ARGUMENT_BINARY_MULTI_MAXCUT_MOGOMEA_FILE_TEMPLATE, i_number_of_elements, i_measure);

		if (!c_err)
		{
			FILE *pf_config = fopen(s_config_file_path, "r");

			if (pf_config != NULL)
			{
				pc_measure->eLoadFromFile(pf_config);
				fclose(pf_config);

				/*FILE  *pf_test;
				CString  s_buf;
				s_buf.Format("zzz_max_cut_test_%d.txt", i_measure);
				pf_test = fopen(s_buf, "w+");
				pc_measure->eSaveToFile(pf_test);
				fclose(pf_test);//*/
			}//if (pf_config != NULL)
			else
			{
				c_err.vSetError(CError::iERROR_CODE_SYSTEM_FILE_NOT_FOUND, s_config_file_path);
			}//else if (pf_config != NULL)
		}//if (!c_error)
	}//for  (int  i_measure = 0; i_measure < 2; i_measure++)


	
	
	CStringCommandParam p_optimal_pf_file(EVALUATION_ARGUMENT_BINARY_MULTI_MAXCUT_OPTIMAL_PF_FILE, false);
	CString  s_file_opt_pf;
	s_file_opt_pf = p_optimal_pf_file.sGetValue(psSettings, &c_err);
	if (c_err)  return(c_err);

	if (p_optimal_pf_file.bHasValue() == true)
	{
		FILE  *pf_opt_pf;
		pf_opt_pf = fopen(s_file_opt_pf, "r");
		if (pf_opt_pf == NULL)
		{
			c_err.vSetError(CError::iERROR_CODE_SYSTEM_FILE_NOT_FOUND, s_file_opt_pf);
			return(c_err);
		}//if (pf_opt_pf == NULL)

		eLoadPFFromFile(pf_opt_pf);

		//::Tools::vShow(v_theoretical_optimal_pareto_front.size());

		fclose(pf_opt_pf);
	}//if  (p_optimal_pf_file.bHasValue() == true)







	/*CUIntCommandParam p_reversed_it_order_for_0_measure(EVALUATION_ARGUMENT_REVERSE_BIT_ORDER);
	int i_reverse_bit_order_for_0_measure;
	i_reverse_bit_order_for_0_measure = p_reversed_it_order_for_0_measure.iGetValue(psSettings, &c_err);
	if (c_err)  return(c_err);


	CUIntCommandParam p_optimal_pf_bit_step(EVALUATION_ARGUMENT_BINARY_MULTI_OPTIMAL_PF_BIT_STEP);
	int  i_create_optimal_pf_bit_step;
	i_create_optimal_pf_bit_step = p_optimal_pf_bit_step.iGetValue(psSettings, &c_err);
	if (c_err)  return(c_err);


	pc_evaluation_inner = EvaluationUtils::pcGetEvaluation<CBinaryCoding>(psSettings, &c_err);
	if (c_err)  return(c_err);


	i_number_of_elements = pc_evaluation_inner->iGetNumberOfElements();
	d_max_value = pc_evaluation_inner->dGetMaxValue();


	//this MUST be after pc_evaluation_inner creation
	CBinaryMultiReversingMeasure1s  *pc_measure_1s;
	pc_measure_1s = new CBinaryMultiReversingMeasure1s(this, pc_evaluation_inner->dGetMaxValue());
	pc_measure_1s->dWeight = 0.5;
	v_mesaures.push_back(pc_measure_1s);

	CBinaryMultiReversingMeasure0s  *pc_measure_0s;
	pc_measure_0s = new CBinaryMultiReversingMeasure0s(this, pc_evaluation_inner->dGetMaxValue());
	pc_measure_0s->dWeight = 0.5;
	v_mesaures.push_back(pc_measure_0s);
	pc_measure_0s->vSetReversingBitOrder(i_reverse_bit_order_for_0_measure);


	//this MUST be after measures
	vFillOptimalPf(i_create_optimal_pf_bit_step);*/


	return(c_err);
};//virtual CError CBinaryMultiMaxcutMoGomea::eConfigure(istream *psSettings)





//---------------------------------------------------CBinaryMultiKnapsackMoGomeaMeasure---------------------------------------------------
CBinaryMultiKnapsackMoGomeaMeasure::CBinaryMultiKnapsackMoGomeaMeasure(CBinaryMultiObjectiveProblem  *pcParentProblem, double  dMax) : CMultiObjectiveMeasure(pcParentProblem)
{
	s_name = s_MULTI_KNAPSACK_MEASURE;
	d_max = dMax;
}//CBinaryMultiKnapsackMoGomeaMeasure::CBinaryMultiKnapsackMoGomeaMeasure(CBinaryMultiObjectiveProblem  *pcParentProblem, double  dMax) : CMultiObjectiveMeasure(pcParentProblem)


CBinaryMultiKnapsackMoGomeaMeasure::~CBinaryMultiKnapsackMoGomeaMeasure()
{

}//CBinaryMultiKnapsackMoGomeaMeasure::~CBinaryMultiKnapsackMoGomeaMeasure()


double  CBinaryMultiKnapsackMoGomeaMeasure::dEvaluate(CBinaryCoding *pcFenotype)
{
	double  d_result;

	d_result = 0;
	for (int ii = 0; ii < pcFenotype->iGetNumberOfBits(); ii++)
	{
		if (pcFenotype->piGetBits()[ii] == 1)
			d_result += v_items.at(ii).d_profit;
	}//for (int ii = 0; ii < pcFenotype->iGetNumberOfBits(); ii++)

	
	return(d_result);
}//double  CBinaryMultiKnapsackMoGomeaMeasure::dEvaluate(CBinaryCoding *pcFenotype)


void  CBinaryMultiKnapsackMoGomeaMeasure::vConstraintSatisfied(CBinaryCoding *pcFenotype, vector<CBinaryMultiKnapsackMoGomeaMeasureItemDataRatio>  *pvItemsRatioOrder)
{
	double  d_total_weight;
	d_total_weight = 0;

	for (int ii = 0; ii < pcFenotype->iGetNumberOfBits(); ii++)
	{
		if (pcFenotype->piGetBits()[ii] == 1)
			d_total_weight += v_items.at(ii).d_weight;
	}//for (int ii = 0; ii < pcFenotype->iGetNumberOfBits(); ii++)


	for (
		int  i_item_to_remove = 0;
		(i_item_to_remove < pvItemsRatioOrder->size()) && (d_total_weight > d_capacity);
		i_item_to_remove++
		)
	{
		if (pcFenotype->piGetBits()[pvItemsRatioOrder->at(i_item_to_remove).i_item_offset] == 1)
		{
			pcFenotype->piGetBits()[pvItemsRatioOrder->at(i_item_to_remove).i_item_offset] = 0;
			d_total_weight -= v_items.at(pvItemsRatioOrder->at(i_item_to_remove).i_item_offset).d_weight;
		}//if (pcFenotype->piGetBits()[pvItemsRatioOrder->at(i_item_to_remove).i_item_offset] == 1)
	}//for (
	

}//void  CBinaryMultiKnapsackMoGomeaMeasure::vConstraintSatisfied(CBinaryCoding *pcFenotype, vector<CBinaryMultiKnapsackMoGomeaMeasureItemDataRatio>  *pvItemsRatioOrder)



CError  CBinaryMultiKnapsackMoGomeaMeasureItemData::eLoad(FILE  *pfSource, int iItemOffset)
{
	CError  c_err;

	CString  s_line, s_buf;
	int  i_buf;


	i_item_offset = ::Tools::iReadLine(pfSource);
	i_item_offset--;
	if (i_item_offset != iItemOffset)
	{
		s_buf.Format("WRONG item offset: (%d != %d)   CError  CBinaryMultiKnapsackMoGomeaMeasureItemData::eLoad(FILE  *pfSource, int iItemOffset)", i_item_offset, iItemOffset);
		c_err.vSetError(s_buf);
		return(c_err);
	}//if (i_item_offset != iItemOffset)

	s_line = ::Tools::sReadLine(pfSource);
	d_weight = ::Tools::dExtractFromString(s_line, 1);

	s_line = ::Tools::sReadLine(pfSource);
	d_profit = ::Tools::dExtractFromString(s_line, 1);


	return(c_err);
}//CError  CBinaryMultiKnapsackMoGomeaMeasureItemData::eLoad(FILE  *pfSource)



CError  CBinaryMultiKnapsackMoGomeaMeasureItemData::eSave(FILE  *pfDest)
{
	CError  c_err;

	fprintf(pfDest, "%d\n", i_item_offset);
	fprintf(pfDest, "+%.2lf\n", d_weight);
	fprintf(pfDest, "+%.2lf\n", d_profit);

	return(c_err);
}//CError  CBinaryMultiKnapsackMoGomeaMeasureItemData::eSave(FILE  *pfDest)



CError  CBinaryMultiKnapsackMoGomeaMeasure::eSave(FILE  *pfSource)
{
	CError  c_err;


	fprintf(pfSource, "+%.2lf\n", d_capacity);

	for (int ii = 0; ii < v_items.size(); ii++)
		v_items.at(ii).eSave(pfSource);

	fprintf(pfSource, "\n", d_capacity);

	return(c_err);
}//CError  CBinaryMultiKnapsackMoGomeaMeasure::eSave(FILE  *pfSource)



CError  CBinaryMultiKnapsackMoGomeaMeasure::eLoadFromFile(FILE  *pfSource, int  iArgNumber)
{
	CError  c_err;

	CString  s_line, s_buf;

	s_line = ::Tools::sReadLine(pfSource);

	while  ( (feof(pfSource) == false)&&(s_line.GetLength() == 0) )
		s_line = ::Tools::sReadLine(pfSource);

	d_capacity = ::Tools::dExtractFromString(s_line, 1);

	//::Tools::vShow(d_capacity);

	int  i_item;
	CBinaryMultiKnapsackMoGomeaMeasureItemData  c_item_buf;
	i_item = 0;
	while ((feof(pfSource) == false) && (i_item < iArgNumber))
	{
		c_err = c_item_buf.eLoad(pfSource, i_item);
		if (c_err)  return(c_err);

		v_items.push_back(c_item_buf);
		i_item++;
	}//while ((feof(pfSource) == false) && (i_item < iArgNumber))


	d_max = 0;
	for (int ii = 0; ii < v_items.size(); ii++)
		d_max += v_items.at(ii).d_profit;

	return(c_err);
}//CError  CBinaryMultiKnapsackMoGomeaMeasure::eLoadFromFile(FILE  *pfSource, int  iArgNumber)



//---------------------------------------------------CBinaryMultiKnapsackMoGomea---------------------------------------------------

CBinaryMultiKnapsackMoGomea::CBinaryMultiKnapsackMoGomea() : CBinaryMultiObjectiveProblem::CBinaryMultiObjectiveProblem()
{
	pc_solution_buf = NULL;
}//CBinaryMultiMaxcutMoGomea::CBinaryMultiMaxcutMoGomea()

CBinaryMultiKnapsackMoGomea::CBinaryMultiKnapsackMoGomea(FILE *pfConfig, CError *pcError)
{
	pc_solution_buf = NULL;
	*pcError = e_init(pfConfig);
}//CBinaryMultiKnapsackMoGomea::CBinaryMultiKnapsackMoGomea(FILE *pfConfig, CError *pcError)


CBinaryMultiKnapsackMoGomea::CBinaryMultiKnapsackMoGomea(const CBinaryMultiKnapsackMoGomea &pcOther) : CBinaryMultiObjectiveProblem::CBinaryMultiObjectiveProblem(pcOther)
{
	pc_solution_buf = NULL;
	::Tools::vShow("CBinaryMultiKnapsackMoGomea::CBinaryMultiKnapsackMoGomea(const CBinaryMultiKnapsackMoGomea &pcOther) : CBinaryMultiObjectiveProblem::CBinaryMultiObjectiveProblem(pcOther)   CAN NOT COPY INNER PROBLEM");
	//if (pcOther.pc_evaluation_inner != NULL)
		//pc_evaluation_inner = new pcOther.pc_evaluation_inner.
	//else
		//pc_evaluation_inner = NULL;
}//CBinaryMultiMaxcutMoGomea::CBinaryMultiMaxcutMoGomea(const CBinaryMultiMaxcutMoGomea &pcOther) : CBinaryMultiObjectiveProblem::CBinaryMultiObjectiveProblem(pcOther)


CBinaryMultiKnapsackMoGomea::~CBinaryMultiKnapsackMoGomea()
{
	if (pc_solution_buf != NULL)  delete  pc_solution_buf;
}//CBinaryMultiKnapsackMoGomea::~CBinaryMultiKnapsackMoGomea()


CError CBinaryMultiKnapsackMoGomea::eConfigure(istream *psSettings)
{
	CError  c_err;


	c_err = CBinaryMultiObjectiveProblem::eConfigure(psSettings);
	if (c_err)  return(c_err);



	CUIntCommandParam p_config_arg_num(EVALUATION_ARGUMENT_BINARY_MULTI_KNAPSACK_MOGOMEA_ARG_NUM);

	i_number_of_elements = p_config_arg_num.iGetValue(psSettings, &c_err);
	//d_max_value = pc_evaluation_inner->dGetMaxValue();
	d_max_value = 1;

	if (c_err)  return(c_err);


	//CBinaryMultiKnapsackMoGomeaMeasure  *pc_measure;

	//pc_measure = new CBinaryMultiKnapsackMoGomeaMeasure(this, 1);
	//pc_measure->dWeight = 0.5;
	//v_mesaures.push_back(pc_measure);

	CString  s_config_file_path;
	s_config_file_path.Format(EVALUATION_ARGUMENT_BINARY_MULTI_KNAPSACK_MOGOMEA_FILE_TEMPLATE, i_number_of_elements);



	if (!c_err)
	{
		FILE *pf_config = fopen(s_config_file_path, "r");

		if (pf_config != NULL)
		{
			eLoadFromFile(pf_config);
			fclose(pf_config);
		}//if (pf_config != NULL)
		else
		{
			c_err.vSetError(CError::iERROR_CODE_SYSTEM_FILE_NOT_FOUND, s_config_file_path);
		}//else if (pf_config != NULL)
	}//if (!c_error)




	CStringCommandParam p_optimal_pf_file(EVALUATION_ARGUMENT_BINARY_MULTI_KNAPSACK_OPTIMAL_PF_FILE, false);
	
		
	CString  s_file_opt_pf;
	s_file_opt_pf = p_optimal_pf_file.sGetValue(psSettings, &c_err);
	if (c_err)  return(c_err);

	if (p_optimal_pf_file.bHasValue() == true)
	{
		FILE  *pf_opt_pf;
		pf_opt_pf = fopen(s_file_opt_pf, "r");
		if (pf_opt_pf == NULL)
		{
			c_err.vSetError(CError::iERROR_CODE_SYSTEM_FILE_NOT_FOUND, s_file_opt_pf);
			return(c_err);
		}//if (pf_opt_pf == NULL)

		eLoadPFFromFile(pf_opt_pf);

		//::Tools::vShow(v_theoretical_optimal_pareto_front.size());

		fclose(pf_opt_pf);
	}//if  (p_optimal_pf_file.bHasValue() == true)


	return(c_err);
};//CError CBinaryMultiKnapsackMoGomea::eConfigure(istream *psSettings)



CError  CBinaryMultiKnapsackMoGomea::eLoadFromFile(FILE  *pfSource)
{
	CError  c_err;

	CString  s_line, s_buf;
	int  i_objectives_num, i_arg_num;
	int  i_index;

	s_line = ::Tools::sReadLine(pfSource);

	i_index = 0;
	i_objectives_num = ::Tools::iExtractFromString(s_line, i_index, &i_index);
	i_arg_num = ::Tools::iExtractFromString(s_line, i_index, &i_index);

	//s_buf.Format("arg: %d,  obj: %d", i_arg_num, i_objectives_num);
	//::Tools::vShow(s_buf);

	if (i_objectives_num != 2)
	{
		s_buf.Format("CError  CBinaryMultiKnapsackMoGomea::eLoadFromFile(FILE  *pfSource)     if (i_objectives_num != 2)");
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (i_objectives_num != 2)


	if (i_arg_num != i_number_of_elements)
	{
		s_buf.Format("CError  CBinaryMultiKnapsackMoGomea::eLoadFromFile(FILE  *pfSource)     if (i_arg_num != i_number_of_elements)");
		c_err.vSetError(s_buf);

		return(c_err);
	}//if (i_objectives_num != 2)



	CBinaryMultiKnapsackMoGomeaMeasure  *pc_measure;
	for (int i_mesaure = 0; i_mesaure < 2; i_mesaure++)
	{
		pc_measure = new CBinaryMultiKnapsackMoGomeaMeasure(this, 1);
		pc_measure->dWeight = 0.5;
		v_measures.push_back(pc_measure);

		pc_measure->eLoadFromFile(pfSource, i_number_of_elements);
	}//for (int i_mesaure = 0; i_mesaure < 2; i_mesaure++)


	CBinaryMultiKnapsackMoGomeaMeasureItemDataRatio  c_ratio_buf;
	for (int i_item = 0; i_item < i_number_of_elements; i_item++)
	{
		c_ratio_buf.i_item_offset = i_item;
		c_ratio_buf.d_ratio = 0;
		
		for (int i_mesaure = 0; i_mesaure < v_measures.size(); i_mesaure++)
			c_ratio_buf.d_ratio += ((CBinaryMultiKnapsackMoGomeaMeasure *)v_measures.at(i_mesaure))->v_items.at(i_item).dGetRatio();

		c_ratio_buf.d_ratio = c_ratio_buf.d_ratio / v_measures.size();
		v_items_ratio_order.push_back(c_ratio_buf);
	}//for  (int  i_item = 0; i_item < i_number_of_elements; i_item++)

	std::sort(v_items_ratio_order.begin(), v_items_ratio_order.end());


	/*//checkup
	FILE  *pf_test;
	pf_test = fopen("zzz_knapsack_test.txt", "w+");
	eSave(pf_test);
	fclose(pf_test);*/

	

	return(c_err);
};//CError  CBinaryMultiKnapsackMoGomea::eLoadFromFile(FILE  *pfSource)



CError  CBinaryMultiKnapsackMoGomea::eSave(FILE  *pfSource)
{
	CError  c_err;

	for (int i_mesaure = 0; i_mesaure < 2; i_mesaure++)
	{
		((CBinaryMultiKnapsackMoGomeaMeasure *)v_measures.at(i_mesaure))->eSave(pfSource);
	}//for (int i_mesaure = 0; i_mesaure < 2; i_mesaure++)

	return(c_err);
};//CError  CBinaryMultiKnapsackMoGomea::eSave(FILE  *pfSource)



void  CBinaryMultiKnapsackMoGomea::vEvaluateParetoFront(vector<double> *pvPF, CBinaryCoding *pcFenotype)
{
	v_prepare_solution(pcFenotype);
	v_evaluate_pareto_front(pvPF, pc_solution_buf);
}//void  CBinaryMultiKnapsackMoGomea::vEvaluateParetoFront(vector<double> *pvPF, CBinaryCoding *pcFenotype)



void  CBinaryMultiKnapsackMoGomea::v_prepare_solution(CBinaryCoding *pcFenotype)
{
	if (pc_solution_buf == NULL)  pc_solution_buf = new CBinaryCoding(pcFenotype->iGetNumberOfBits());

	for (int ii = 0; ii < pcFenotype->iGetNumberOfBits(); ii++)
		pc_solution_buf->piGetBits()[ii] = pcFenotype->piGetBits()[ii];

	for (int i_measure = 0; i_measure < v_measures.size(); i_measure++)
	{
		((CBinaryMultiKnapsackMoGomeaMeasure *)v_measures.at(i_measure))->vConstraintSatisfied(pc_solution_buf, &v_items_ratio_order);
	}//for (int i_measure = 0; i_measure < v_measures.size(); i_measure++)

};//void  CBinaryMultiKnapsackMoGomea::v_prepare_solution(CBinaryCoding *pcFenotype)



CError  CBinaryMultiKnapsackMoGomea::eLoadPFFromFile(FILE  *pfSource)
{
	CError  c_err;

	CString  s_buf;

	s_buf = ::Tools::sReadLine(pfSource);

	bool  b_loaded;
	CMultiIndividualPFpoints *pc_pf_point_buf;
	b_loaded = true;
	
	while ((!feof(pfSource)) && (b_loaded == true))
	{
		pc_pf_point_buf = new CMultiIndividualPFpoints();
		b_loaded = pc_pf_point_buf->bLoadFromFile_MaxCutMoGomea(pfSource);

		if (b_loaded == true)
			v_theoretical_optimal_pareto_front.push_back(pc_pf_point_buf);
		else
			delete  pc_pf_point_buf;			
	}//while ((!feof(pfSource)) && (i_pareto_points_num < v_theoretical_optimal_pareto_front.size()))
	

	return(c_err);
};//CError  CBinaryMultiKnapsackMoGomea::eLoadPFFromFile(FILE  *pfSource)
