#include "nsga2.h"
#include "PopulationOptimizer.h"
#include "UIntCommandParam.h"
#include "FloatCommandParam.h"



using  namespace Nsga2;


uint32_t CNSGA2::iERROR_PARENT_CNSGA2Optimizer = CError::iADD_ERROR_PARENT("iERROR_PARENT_CNSGA2Optimizer");
uint32_t CNSGA2::iERROR_CODE_NSGA2_GENOTYPE_LEN_BELOW_0 = CError::iADD_ERROR("iERROR_CODE_NSGA2_GENOTYPE_LEN_BELOW_0");


//---------------------------------------------CNSGA2-------------------------------------------------------
CNSGA2::CNSGA2(CProblem<CBinaryCoding, CBinaryCoding>* pcProblem, CLog * pcLog, uint32_t iRandomSeed)
	: CBinaryMultiObjectiveOptimizer(pcProblem, pcLog, iRandomSeed)
	//: CBinaryOptimizer(pcProblem, pcLog, iRandomSeed)
{
	//pc_problem = (CBinaryMultiObjectiveProblem *) pcProblem;

}//CNSGA2::CNSGA2(CProblem<CBinaryCoding, CBinaryCoding>* pcProblem, CLog * pcLog, uint32_t iRandomSeed)


CNSGA2::CNSGA2(CNSGA2 *pcOther) : CBinaryMultiObjectiveOptimizer(pcOther)//CBinaryOptimizer(pcOther)
{
	::MessageBox(NULL, "No implementation: CNSGA2::CNSGA2(CNSGA2 *pcOther) : CBinaryMultiObjectiveOptimizer(pcOther)//CBinaryOptimizer(pcOther)", "Implementation missing", MB_OK);
}//CNSGA2::CNSGA2(CNSGA2 *pcOther)


CError CNSGA2::eConfigure(istream *psSettings)
{
	CError c_err(iERROR_PARENT_CNSGA2Optimizer);

	c_err = CBinaryOptimizer::eConfigureNSGA2(psSettings);



	if (!c_err)
	{
		CUIntCommandParam p_population_size(POPULATION_OPTIMIZER_ARGUMENT_POPULATION_SIZE);
		i_pop_size = p_population_size.iGetValue(psSettings, &c_err);
	}//if (!c_err)


	if (!c_err)
	{
		CUIntCommandParam p_tournament_size(NSGA2_ARGUMENT_TOUR_SIZE);
		i_tournament_size = p_tournament_size.iGetValue(psSettings, &c_err);
	}//if (!c_err)

	if (!c_err)
	{
		CFloatCommandParam p_prob_cross(NSGA2_ARGUMENT_PROB_CROSS);
		d_prob_cross = p_prob_cross.fGetValue(psSettings, &c_err);
	}//if (!c_err)

	if (!c_err)
	{
		CFloatCommandParam p_prob_mut(NSGA2_ARGUMENT_PROB_MUT);
		d_prob_mut = p_prob_mut.fGetValue(psSettings, &c_err);
	}//if (!c_err)

	if (!c_err)
	{
		CUIntCommandParam p_full_genotype_mut(NSGA2_FULL_GENOTYPE_MUT);
		if (p_full_genotype_mut.iGetValue(psSettings, &c_err) == 1)
			b_whole_genotype_mut = true;
		else
			b_whole_genotype_mut = false;
	}//if (!c_err)


	return c_err;
}//CError CNSGA2::eConfigure(istream *psSettings)



CNSGA2::~CNSGA2()
{
	//v_non_dominated_pfs; is not an owner - no deletion

	for (int ii = 0; ii < v_population.size(); ii++)
		delete  v_population.at(ii);
}//CNSGA2::~CNSGA2()




CString  CNSGA2::sAdditionalSummaryInfo()
{
	CString  s_result;
		
	s_result.Format("CString  CNSGA2::sAdditionalSummaryInfo() - JUST FILL ME UP!");

	return(s_result);
}//CString  CNSGA2::sAdditionalSummaryInfo()




void CNSGA2::vInitialize(time_t tStartTime)
{
	CBinaryOptimizer::vInitialize(tStartTime);
	t_start = tStartTime;

	CError  c_err(iERROR_PARENT_CNSGA2Optimizer);
	i_templ_length = pc_problem->pcGetEvaluation()->iGetNumberOfElements();

	if (i_templ_length <= 0)
	{
		c_err.vSetError(CNSGA2::iERROR_CODE_NSGA2_GENOTYPE_LEN_BELOW_0);
		return;
	}//if  (i_templ_length  <=  0)


	pc_log->vPrintLine("Initializing...", true);


	c_time_counter.vSetStartNow();
	pc_log->vPrintLine("DONE...", true);


	CNSGA2Individual  *pc_buf;
	for (int ii = 0; ii < i_pop_size; ii++)
	{
		pc_buf = new CNSGA2Individual(i_templ_length, (CBinaryMultiObjectiveProblem *) pc_problem->pcGetEvaluation(), this);
		pc_buf->vInitialize();
		v_population.push_back(pc_buf);
	}//for (int ii = 0; ii < i_pop_size; ii++)

	c_time_counter.vSetStartNow();	


	
}//void CNSGA2::vInitialize(time_t tStartTime)




bool CNSGA2::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)
{
	CString  s_buf;
	double  d_time_passed;

	c_time_counter.bGetTimePassed(&d_time_passed);

	for (int ii = 0; ii < v_population.size(); ii++)
		v_population.at(ii)->vRate();

	v_non_dominated_sorting();
	v_compute_crowding_distance();
	//v_report_pareto_front(0);
	//::Tools::vShow("asd");

	v_evolution();

	v_non_dominated_sorting();
	v_compute_crowding_distance();
	v_report_pareto_front(1);
	
	v_get_half_ofChildren_parent_pop();

	//v_update_global_pareto_front();

	CBinaryMultiObjectiveProblem  *pc_multi_problem;
	pc_multi_problem = (CBinaryMultiObjectiveProblem *)pc_problem->pcGetEvaluation();

	c_time_counter.bGetTimePassed(&d_time_passed);
	//s_buf.Format("iteration: %d  pop: %d; GenSize:%d PFsize:%d HyperVolume: %.8lf InvGenDist: %.8lf DominatedPF: %d [time:%.2lf]", iIterationNumber, v_population.size(), v_non_dominated_pfs.at(0).at(0)->pc_genotype->iGetNumberOfBits(), (int) dPFQualityPointNum(), dPFQualityHyperVolume(), dPFQualityInverseGenerationalDistance(), (int)dPFQualityDominatedOptimalPoints(), d_time_passed);
	s_buf.Format("iteration: %d  pop: %d; GenSize:%d PFsize:%d HyperVolume: %.8lf InvGenDist: %.8lf GenDist:%.8lf MaxSpread:%.8lf DominatedPF: %d [time:%.2lf] [ffe: %.0lf]", iIterationNumber, v_population.size(), v_non_dominated_pfs.at(0).at(0)->pc_genotype->iGetNumberOfBits(), (int)dPFQualityPointNum(), dPFQualityHyperVolume(), dPFQualityInverseGenerationalDistance(), dPFQualityGenerationalDistance(), dPFQualityMaximumSpread(), (int)dPFQualityDominatedOptimalPoints(), d_time_passed, (double)pc_multi_problem->iGetFFE());
	pc_log->vPrintLine(s_buf, true);

	b_update_best_individual(iIterationNumber, tStartTime, v_non_dominated_pfs.at(0).at(0)->pc_genotype->piGetBits(), v_non_dominated_pfs.at(0).at(0)->dEvaluate());

	return(true);
}//bool CNSGA2::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)





void  CNSGA2::v_get_half_ofChildren_parent_pop()
{
	v_population.clear();

	for (int i_front = 0; (i_front < v_non_dominated_pfs.size()) && (v_population.size() < i_pop_size); i_front++)
	{
		if (v_population.size() + v_non_dominated_pfs.at(i_front).size() <= i_pop_size)
		{
			for (int ii = 0; ii < v_non_dominated_pfs.at(i_front).size(); ii++)
			{
				v_population.push_back(v_non_dominated_pfs.at(i_front).at(ii));
			}//for (int ii = 0; ii < v_non_dominated_pfs.at(i_front).size(); ii++)
		}//if (v_population.size() + v_non_dominated_pfs.size() <= i_pop_size)
		else
		{
			if (v_population.size() >= i_pop_size)
			{
				for (int ii = 0; ii < v_non_dominated_pfs.at(i_front).size(); ii++)
					delete  v_non_dominated_pfs.at(i_front).at(ii);
			}//if (v_population.size() >= i_pop_size)
			else
			{
				std::sort
				(
					v_non_dominated_pfs.at(i_front).begin(), v_non_dominated_pfs.at(i_front).end(),
					[](CNSGA2Individual *pc0, CNSGA2Individual *pc1) {return pc0->d_crowding > pc1->d_crowding; }
				);

				while ((0 < v_non_dominated_pfs.at(i_front).size()) && (v_population.size() < i_pop_size))
				{
					v_population.push_back(v_non_dominated_pfs.at(i_front).at(0));
					v_non_dominated_pfs.at(i_front).erase(v_non_dominated_pfs.at(i_front).begin());
				}//while ( (0 < v_non_dominated_pfs.at(i_front).size())&&(v_population.size() < i_pop_size) )
			}//else  if (v_population.size() >= i_pop_size)
		}//else if (v_population.size() + v_non_dominated_pfs.size() <= i_pop_size)

	}//for (int i_front = 0; i_front < v_non_dominated_pfs.size(); i_front++)
}//void  CNSGA2::v_get_half_ofChildren_parent_pop()




CNSGA2Individual *CNSGA2::pc_parent_tournament(int iTournamentSize)
{
	int  i_parent_offset;
	CNSGA2Individual *pc_parent;
	CNSGA2Individual *pc_parent_other;

	i_parent_offset = RandUtils::iRandNumber(0, v_population.size() - 1);
	pc_parent = v_population.at(i_parent_offset);

	for (int ii = 0; ii < iTournamentSize - 1; ii++)
	{
		i_parent_offset = RandUtils::iRandNumber(0, v_population.size() - 1);
		pc_parent_other = v_population.at(i_parent_offset);

		if (pc_parent_other->i_front < pc_parent->i_front)
			pc_parent = pc_parent_other;
		else
		{
			if (pc_parent_other->i_front == pc_parent->i_front)
			{
				if (pc_parent_other->d_crowding > pc_parent->d_crowding)  pc_parent = pc_parent_other;
			}//if (pc_parent_other->i_front == pc_parent->i_front)
		}//else  if (pc_parent_other->i_front < pc_parent->i_front)

	}//for (int ii = 0; ii < iTournamentSize - 1; ii++)


	return(pc_parent);
}//CNSGA2Individual *CNSGA2::pc_parent_tournament(int iTournamentSize)



void  CNSGA2::v_evolution()
{
	vector <CNSGA2Individual *> v_offsprings;

	CNSGA2Individual  *pc_parent_0, *pc_parent_1;
	vector<CNSGA2Individual *> v_cross_result;

	while (v_offsprings.size() < i_pop_size)
	{
		pc_parent_0 = pc_parent_tournament(i_tournament_size);
		pc_parent_1 = pc_parent_tournament(i_tournament_size);

		if  (b_2_point_cross == true)
			pc_parent_0->vCross2point(pc_parent_1, &v_cross_result);
		else
			pc_parent_0->vCross(pc_parent_1, &v_cross_result);


		for (int i_child = 0; i_child < v_cross_result.size(); i_child++)
		{
			if (b_whole_genotype_mut == true)
				v_cross_result.at(i_child)->vMutateWholeGenotype();
			else
				v_cross_result.at(i_child)->vMutate();


			v_cross_result.at(i_child)->vRate();
			v_offsprings.push_back(v_cross_result.at(i_child));
		}//for (int i_child = 0; i_child < v_cross_result.size(); i_child++)

		v_cross_result.clear();
	}//for (int ii = 0; ii < i_pop_size; ii++)

	for (int ii = 0; ii < v_offsprings.size(); ii++)
		v_population.push_back(v_offsprings.at(ii));
}//void  CNSGA2::v_evolution()



void  CNSGA2::v_non_dominated_sorting()
{
	v_non_dominated_pfs.clear();

	vector<CNSGA2Individual  *>  v_new_front;
	vector<CNSGA2Individual  *>  v_ind_to_add;
	v_ind_to_add = v_population;


	bool  b_dominated;
	int  i_domination;

	while (v_ind_to_add.size() > 0)
	{
		for (int i_cand = 0; i_cand < v_ind_to_add.size(); i_cand++)
		{
			b_dominated = false;

			for (int i_opp = 0; (i_opp < v_new_front.size()) && (b_dominated == false); i_opp++)
			{
				i_domination = v_ind_to_add.at(i_cand)->iCheckDomination(v_new_front.at(i_opp));
				if (i_domination < 0)  b_dominated = true;
			}//for (int i_opp = 0; i_opp < v_ind_to_add.size(); i_opp++)

			for (int i_opp = 0; (i_opp < v_ind_to_add.size()) && (b_dominated == false); i_opp++)
			{
				i_domination = v_ind_to_add.at(i_cand)->iCheckDomination(v_ind_to_add.at(i_opp));
				if (i_domination < 0)  b_dominated = true;
			}//for (int i_opp = 0; i_opp < v_ind_to_add.size(); i_opp++)


			if (b_dominated == false)
			{
				v_new_front.push_back(v_ind_to_add.at(i_cand));
				v_ind_to_add.erase(v_ind_to_add.begin() + i_cand);
				i_cand--;
			}//if (b_dominated == false)
		}//for (int i_cand = 0; i_cand < v_ind_to_add.size(); i_cand++)

		v_non_dominated_pfs.push_back(vector<CNSGA2Individual *>());

		for (int ii = 0; ii < v_new_front.size(); ii++)
			v_add_to_non_dominated_pf((vector<CMultiIndividual *> *) &(v_non_dominated_pfs.at(v_non_dominated_pfs.size() - 1)), v_new_front.at(ii));

		v_new_front.clear();
	}//while  (v_ind_to_add.size() > 0)


	//set front info in individuals
	for (int i_front = 0; i_front < v_non_dominated_pfs.size(); i_front++)
	{
		for (int i_ind = 0; i_ind < v_non_dominated_pfs.at(i_front).size(); i_ind++)
		{
			v_non_dominated_pfs.at(i_front).at(i_ind)->i_front = i_front;
		}//for (int i_ind = 1; i_ind < v_non_dominated_pfs.at(i_front).size() - 1; i_ind++)
	}//for (int i_front = 0; i_front < v_non_dominated_pfs.size(); i_front++)


}//void  CNSGA2::v_non_dominated_sorting()



void  CNSGA2::v_add_to_non_dominated_pf(vector<CMultiIndividual *>  *pvFrontToFill, CNSGA2Individual  *pcInd)
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
};//void  void  CNSGA2::v_add_to_non_dominated_pf(vector<CMultiIndividual *>  *pvFrontToFill, CBinaryMultiObjectiveProblem  *pcInd)





void  CNSGA2::v_compute_crowding_distance()
{
	double  d_max_crowding_distance;
	
	for (int i_front = 0; i_front < v_non_dominated_pfs.size(); i_front++)
	{
		d_max_crowding_distance = 0;

		for (int i_ind = 1; i_ind < v_non_dominated_pfs.at(i_front).size() - 1; i_ind++)
		{
			v_non_dominated_pfs.at(i_front).at(i_ind)->vCompCrowdingDistance(v_non_dominated_pfs.at(i_front).at(i_ind - 1), v_non_dominated_pfs.at(i_front).at(i_ind + 1));

			if (d_max_crowding_distance < v_non_dominated_pfs.at(i_front).at(i_ind)->dGetCrowdingDistance())  d_max_crowding_distance = v_non_dominated_pfs.at(i_front).at(i_ind)->dGetCrowdingDistance();
		}//for (int i_ind = 1; i_ind < v_non_dominated_pfs.at(i_front).size() - 1; i_ind++)

		v_non_dominated_pfs.at(i_front).at(0)->d_crowding = d_max_crowding_distance + 1;
		v_non_dominated_pfs.at(i_front).at(v_non_dominated_pfs.at(i_front).size() - 1)->d_crowding = d_max_crowding_distance + 1;
	}//for (int i_front = 0; i_front < v_non_dominated_pfs.size(); i_front++)

}//void  CNSGA2::v_compute_crowding_distance()




void  CNSGA2::v_non_dominated_sorting_old()
{
	v_non_dominated_pfs.clear();


	CNSGA2Individual  *pc_cur_ind;

	int  i_current_front;


	for (int i_ind = 0; i_ind < v_population.size(); i_ind++)
	{
		pc_cur_ind = v_population.at(i_ind);

		if (v_non_dominated_pfs.size() == 0)
		{
			v_non_dominated_pfs.push_back(vector<CNSGA2Individual *>());
			v_non_dominated_pfs.at(0).push_back(pc_cur_ind);
		}//if (i_current_front < 0)
		else
		{
			vector<CNSGA2Individual *>  v_dropped_individuals;
			int  i_domination;
			bool  b_add_to_this_level = true;


			i_current_front = 0;

			while (b_add_to_this_level == true)
			{
				for (int i_pf_ind = 0; (i_pf_ind < v_non_dominated_pfs.at(i_current_front).size())&&(b_add_to_this_level == true); i_pf_ind++)
				{
					i_domination = pc_cur_ind->iCheckDomination(v_non_dominated_pfs.at(i_current_front).at(i_pf_ind));

					//form new level
					if (i_domination > 0)
					{
						v_dropped_individuals.push_back(v_non_dominated_pfs.at(i_current_front).at(i_pf_ind));
						v_non_dominated_pfs.at(i_current_front).erase(v_non_dominated_pfs.at(i_current_front).begin() + i_pf_ind);
						i_pf_ind--;
					}//if (i_domination > 0)

					if (i_domination < 0)
					{
						b_add_to_this_level = false;
						if (v_dropped_individuals.size() > 0)  ::Tools::vShow("void  CNSGA2::v_non_dominated_sorting() - PROBLEM");
					}//if (i_domination < 0)
				}//for (int i_pf_ind = 0; i_pf_ind < v_non_dominated_pfs.size(); i_pf_ind++)

				if (b_add_to_this_level == true)
				{
					v_add_to_non_dominated_pf((vector<CMultiIndividual *> *) &(v_non_dominated_pfs.at(i_current_front)), pc_cur_ind);
					b_add_to_this_level = false;//to break the while loop

					if (v_non_dominated_pfs.size() > i_current_front + 1)
						v_non_dominated_pfs.insert(v_non_dominated_pfs.begin() + i_current_front + 1, v_dropped_individuals);
					else
						v_non_dominated_pfs.push_back(v_dropped_individuals);

					
				}//if (b_add_to_this_level == true)
				else
				{
					i_current_front++;
					if (i_current_front >= v_non_dominated_pfs.size())
					{
						v_non_dominated_pfs.push_back(vector<CNSGA2Individual *>());
						v_non_dominated_pfs.at(i_current_front).push_back(pc_cur_ind);
					}//if (i_current_front >= v_non_dominated_pfs.size())

					b_add_to_this_level = true;//to continue the while loop
				}//else  if (b_add_to_this_level == true)

			}//while  (b_add_to_this_level  ==  true)
		}//else  if (i_current_front < 0)

	}//for (int i_ind = 0; i_ind < v_population.size(); i_ind++)
	
}//void  CNSGA2::v_non_dominated_sorting()






/*void  CNSGA2::v_add_to_non_dominated_pf(vector<CNSGA2Individual *>  *pvFrontToFill, CNSGA2Individual  *pcInd)
{
	for (int ii = 0; ii < pvFrontToFill->size(); ii++)
	{
		if (pcInd->pvGetFitness()->at(0) > pvFrontToFill->at(ii)->pvGetFitness()->at(0))
		{
			pvFrontToFill->insert(pvFrontToFill->begin() + ii, pcInd);
			return;
		}//if (pcInd->pvGetFitness()->at(0) > pvFrontToFill->at(ii)->pvGetFitness()->at(0))
		else
		{
			if (pcInd->pvGetFitness()->at(0) == pvFrontToFill->at(ii)->pvGetFitness()->at(0))
			{
				if (pcInd->pvGetFitness()->size() > 1)
				{
					if (pcInd->pvGetFitness()->at(1) < pvFrontToFill->at(ii)->pvGetFitness()->at(1))
					{
						pvFrontToFill->insert(pvFrontToFill->begin() + ii, pcInd);
						return;
					}//if (pcInd->pvGetFitness()->at(1) < pvFrontToFill->at(ii)->pvGetFitness()->at(1))
				}//if  (pcInd->pvGetFitness()->size() > 1)
			}//if (pcInd->pvGetFitness()->at(0) == pvFrontToFill->at(ii)->pvGetFitness()->at(0))
		}//else  if (pcInd->pvGetFitness()->at(0) > pvFrontToFill->at(ii)->pvGetFitness()->at(0))

	}//for (int ii = 0; ii < pvFrontToFill->size(); ii++)

	pvFrontToFill->push_back(pcInd);
};//void  CNSGA2::v_add_to_non_dominated_pf(vector<CNSGA2Individual *>  *pvFrontToFill, CNSGA2Individual  *pcInd)*/



void  CNSGA2::v_report_pareto_front(int  iId)
{
	CString  s_filename;
	CString  s_line, s_buf;

	s_filename.Format("zz_pareto_%.2d.txt", iId);
	::Tools::vReportInFile(s_filename, "", true);
	
	for (int i_front = 0; i_front < v_non_dominated_pfs.size(); i_front++)
	{
		s_line.Format("FRONT: %d", i_front);
		::Tools::vReportInFile(s_filename, s_line);

		for (int i_obj = 0; i_obj < ((CBinaryMultiObjectiveProblem *)pc_problem->pcGetEvaluation())->iGetActiveMeasures(); i_obj++)
		{
			s_line.Format("OBJ (%.2d): \t", i_obj);

			for (int i_ind = 0; i_ind < v_non_dominated_pfs.at(i_front).size(); i_ind++)
			{
				s_buf.Format("%.8lf\t", v_non_dominated_pfs.at(i_front).at(i_ind)->pvGetFitness()->at(i_obj));
				s_line += s_buf;				
			}//for (int i_ind = 0; i_ind < v_pareto_fronts.at(i_front).size(); i_ind++)

			s_line += "\n";
			::Tools::vReportInFile(s_filename, s_line);
		}//for (int i_obj = 0; i_obj < ((CBinaryMultiObjectiveProblem *) pc_problem)->iGetActiveMeasures(); i_obj++)


		s_line.Format("Crowding: \t");

		for (int i_ind = 0; i_ind < v_non_dominated_pfs.at(i_front).size(); i_ind++)
		{
			s_buf.Format("%.8lf\t", v_non_dominated_pfs.at(i_front).at(i_ind)->dGetCrowdingDistance());
			s_line += s_buf;
		}//for (int i_ind = 0; i_ind < v_pareto_fronts.at(i_front).size(); i_ind++)

		//s_line += "\n";
		::Tools::vReportInFile(s_filename, s_line);

		::Tools::vReportInFile(s_filename, "");
	}//for (int i_front = 0; i_front < v_pareto_fronts.size(); i_front++)
}//void  CNSGA2::v_report_pareto_front()







//---------------------------------------------CNSGA2Individual-------------------------------------------------------

CNSGA2Individual::CNSGA2Individual(int  iTemplLength, CBinaryMultiObjectiveProblem *pcProblem, CNSGA2  *pcParent)
{
	i_templ_length = iTemplLength;
	pc_parent = pcParent;


	pc_genotype = new CBinaryCoding(i_templ_length);

	pc_problem = pcProblem;

	b_rated = false;
}//CNSGA2Individual::CNSGA2Individual(int  iTemplLength, CBinaryMultiObjectiveProblem *pcProblem, CNSGA2  *pcParent)



CNSGA2Individual::CNSGA2Individual(CNSGA2Individual &pcOther)
{
	i_templ_length = pcOther.i_templ_length;
	pc_parent = pcOther.pc_parent;
	pc_problem = pcOther.pc_problem;
	b_rated = pcOther.b_rated;
	v_fitness = pcOther.v_fitness;

	pc_genotype = new CBinaryCoding(pcOther.pc_genotype);

	

	//pi_genotype = new int[i_templ_length];
	//for (int ii = 0; ii < i_templ_length; ii++)
		//pi_genotype[ii] = pcOther.pi_genotype[ii];

}//CNSGA2Individual::CNSGA2Individual(CNSGA2Individual &pcOther)


CNSGA2Individual::~CNSGA2Individual()
{
	//delete  pc_genotype;  //it is done in CMultiIndividual
};//CNSGA2Individual::~CNSGA2Individual()



void  CNSGA2Individual::vInitialize()
{
	for (int ii = 0; ii < i_templ_length; ii++)
	{
		pc_genotype->piGetBits()[ii] = RandUtils::iRandNumber(0, 1);
	}//for (int ii = 0; ii < i_templ_length; ii++)

	d_crowding = 0;
	i_front = 0;
	b_rated = false;
};//void  CNSGA2Individual::vInitialize()


double  CNSGA2Individual::dEvaluate()
{
	double  d_res;

	vRate();

	d_res = 1;
	for (int ii = 0; ii < v_fitness.size(); ii++)
		d_res *= v_fitness.at(ii);

	return(d_res);
}//double  CNSGA2Individual::dEvaluate()




/*void  CNSGA2Individual::vRate()
{
	if (b_rated == true)  return;
	pc_problem->vEvaluateParetoFront(&v_fitness, pc_genotype);
	b_rated = true;
};//void  CNSGA2Individual::vRate()*/


/*vector<double>  *CNSGA2Individual::pvGetFitness()
{
	if (b_rated == true)  return(&v_fitness);
	vRate();
	return(&v_fitness);
};//vector<double>  *CNSGA2Individual::pvGetFitness()*/





void  CNSGA2Individual::vCompCrowdingDistance(CNSGA2Individual  *pcPrev, CNSGA2Individual  *pcNext)
{
	double  d_buf;
	d_crowding = 1;

	for (int i_obj = 0; i_obj < pcPrev->pvGetFitness()->size(); i_obj++)
	{
		d_buf = ::abs(pcPrev->pvGetFitness()->at(i_obj) - pcNext->pvGetFitness()->at(i_obj));
		d_crowding *= d_buf;
	}//for (int i_obj = 0; i_obj < pcPrev->pvGetFitness()->size(); i_obj++)
	
};//void  CNSGA2Individual::vCompCrowdingDistance(CNSGA2Individual  *pcPrev, CNSGA2Individual  *pcNext)



void  CNSGA2Individual::vCross2point(CNSGA2Individual  *pcOtherParent, vector<CNSGA2Individual *> *pvCrossResult)
{
	CNSGA2Individual  *pc_child_0, *pc_child_1;

	pc_child_0 = new CNSGA2Individual(*this);
	pc_child_1 = new CNSGA2Individual(*pcOtherParent);

	pvCrossResult->push_back(pc_child_0);
	pvCrossResult->push_back(pc_child_1);

	if (RandUtils::dRandNumber(0, 1) < pc_parent->d_prob_cross)
	{
		int  i_buf;
		int  i_point_0, i_point_1;

		pc_child_0->b_rated = false;
		pc_child_1->b_rated = false;

		i_point_0 = RandUtils::dRandNumber(0, i_templ_length - 1);
		i_point_1 = RandUtils::dRandNumber(0, i_templ_length - 1);

		if (i_point_0 > i_point_1)
		{
			i_buf = i_point_0;
			i_point_0 = i_point_1;
			i_point_1 = i_buf;
		}//if (i_point_0 > i_point_1)

		
		for (int ii = 0; ii < i_point_0; ii++)
		{
			pc_child_0->pc_genotype->piGetBits()[ii] = pc_genotype->piGetBits()[ii];
			pc_child_1->pc_genotype->piGetBits()[ii] = pcOtherParent->pc_genotype->piGetBits()[ii];
		}//for (int ii = 0; ii < i_point_0; ii++)


		for (int ii = i_point_0; ii < i_point_1; ii++)
		{
			pc_child_0->pc_genotype->piGetBits()[ii] = pcOtherParent->pc_genotype->piGetBits()[ii]; 
			pc_child_1->pc_genotype->piGetBits()[ii] = pc_genotype->piGetBits()[ii];
		}//for (int ii = i_point_0; ii < i_point_1; ii++)


		for (int ii = i_point_1; ii < i_templ_length; ii++)
		{
			pc_child_0->pc_genotype->piGetBits()[ii] = pc_genotype->piGetBits()[ii];
			pc_child_1->pc_genotype->piGetBits()[ii] = pcOtherParent->pc_genotype->piGetBits()[ii];
		}//for (int ii = i_point_1; ii < i_templ_length; ii++)

	}//if (RandUtils::dRandNumber(0, 1) < pc_parent->d_prob_cross)
}//void  CNSGA2Individual::vCross2point(CNSGA2Individual  *pcOtherParent, vector<CNSGA2Individual *> *pvCrossResult)




void  CNSGA2Individual::vCross(CNSGA2Individual  *pcOtherParent, vector<CNSGA2Individual *> *pvCrossResult)
{
	CNSGA2Individual  *pc_child_0, *pc_child_1;

	pc_child_0 = new CNSGA2Individual(*this);
	pc_child_1 = new CNSGA2Individual(*pcOtherParent);

	pvCrossResult->push_back(pc_child_0);
	pvCrossResult->push_back(pc_child_1);

	if (RandUtils::dRandNumber(0, 1) < pc_parent->d_prob_cross)
	{
		pc_child_0->b_rated = false;
		pc_child_1->b_rated = false;

		for (int ii = 0; ii < i_templ_length; ii++)
		{
			if (RandUtils::dRandNumber(0, 1) < 0.5)
			{
				pc_child_0->pc_genotype->piGetBits()[ii] = pc_genotype->piGetBits()[ii];
				pc_child_1->pc_genotype->piGetBits()[ii] = pcOtherParent->pc_genotype->piGetBits()[ii];
			}//if (RandUtils::dRandNumber(0, 1) < 0.5)
			else
			{
				pc_child_0->pc_genotype->piGetBits()[ii] = pcOtherParent->pc_genotype->piGetBits()[ii]; 
				pc_child_1->pc_genotype->piGetBits()[ii] = pc_genotype->piGetBits()[ii];
			}//else  if (RandUtils::dRandNumber(0, 1) < 0.5)

		}//for (int ii = 0; ii < i_templ_length; ii++)

	}//if (RandUtils::dRandNumber(0, 1) < pc_parent->d_prob_cross)

};//void  CNSGA2Individual::vCross(CNSGA2Individual  *pcOtherParent, vector<CNSGA2Individual *> vCrossResult)




void  CNSGA2Individual::vMutateWholeGenotype()
{
	double  d_eff_prob_mut;

	d_eff_prob_mut = pc_parent->d_prob_mut;
	d_eff_prob_mut = d_eff_prob_mut / pc_problem->iGetNumberOfElements();

	//::Tools::vShow(d_eff_prob_mut);

	for (int i_gene = 0; i_gene < i_templ_length; i_gene++)
	{
		if (RandUtils::dRandNumber(0, 1) < d_eff_prob_mut)
		{
			pc_genotype->piGetBits()[i_gene] = pc_genotype->piGetBits()[i_gene] * (-1) + 1;
		}//if (RandUtils::dRandNumber(0, 1) < pc_parent->d_prob_mut)
	}//for (int i_gene = 0; i_gene < i_templ_length; i_gene++)
};//void  CNSGA2Individual::vMutateWholeGenotype()


void  CNSGA2Individual::vMutate()
{
	if (RandUtils::dRandNumber(0, 1) < pc_parent->d_prob_mut)
	{
		int  i_gene;

		i_gene = RandUtils::iRandNumber(0, i_templ_length - 1);
		pc_genotype->piGetBits()[i_gene] = pc_genotype->piGetBits()[i_gene] * (-1) + 1;
	}//if (RandUtils::dRandNumber(0, 1) < pc_parent->d_prob_mut)

};//void  CNSGA2Individual::vMutate()

