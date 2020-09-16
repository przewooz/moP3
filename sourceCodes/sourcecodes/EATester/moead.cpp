



#include "moead.h"
#include "PopulationOptimizer.h"
#include "UIntCommandParam.h"
#include "FloatCommandParam.h"


using namespace MoeaD;


uint32_t CMoead::iERROR_PARENT_CMoeaDOptimizer = CError::iADD_ERROR_PARENT("iERROR_PARENT_CMoeaDOptimizer");
uint32_t CMoead::iERROR_CODE_MOEAD_GENOTYPE_LEN_BELOW_0 = CError::iADD_ERROR("iERROR_CODE_MOEAD_GENOTYPE_LEN_BELOW_0");



//---------------------------------------------CMoead-------------------------------------------------------
CMoead::CMoead(CProblem<CBinaryCoding, CBinaryCoding>* pcProblem, CLog * pcLog, uint32_t iRandomSeed)
	: CBinaryMultiObjectiveOptimizer(pcProblem, pcLog, iRandomSeed)
	//: CBinaryOptimizer(pcProblem, pcLog, iRandomSeed)
{
	//pc_problem = (CBinaryMultiObjectiveProblem *) pcProblem;

	pc_genotype = NULL;


};//CMoead::CMoead(CProblem<CBinaryCoding, CBinaryCoding>* pcProblem, CLog * pcLog, uint32_t iRandomSeed)


CMoead::CMoead(CMoead *pcOther) : CBinaryMultiObjectiveOptimizer(pcOther)//CBinaryOptimizer(pcOther)
{
	::MessageBox(NULL, "No implementation: CMoead::CMoead(CMoead *pcOther) : CBinaryMultiObjectiveOptimizer(pcOther)//CBinaryOptimizer(pcOther)", "Implementation missing", MB_OK);
};//CMoead::CMoead(CMoead *pcOther) : CBinaryMultiObjectiveOptimizer(pcOther)//CBinaryOptimizer(pcOther)


CMoead::~CMoead()
{
	if (pc_genotype != NULL)  delete  pc_genotype;

	for (int ii = 0; ii < CurrentPopulation.size(); ii++)
		delete  CurrentPopulation.at(ii);

	for (int ii = 0; ii < SecondPopulation.size(); ii++)
		delete  SecondPopulation.at(ii);

}//CMoead::~CMoead()


CError CMoead::eConfigure(istream *psSettings)
{
	CError c_err(iERROR_PARENT_CMoeaDOptimizer);

	c_err = CBinaryOptimizer::eConfigureMoead(psSettings);

	if (!c_err)
	{
		CUIntCommandParam p_population_size(POPULATION_OPTIMIZER_ARGUMENT_POPULATION_SIZE);
		i_population_size = p_population_size.iGetValue(psSettings, &c_err);
	}//if (!c_err)


	if (!c_err)
	{
		CFloatCommandParam p_prob_cross(MOEAD_ARGUMENT_IPR);
		Pr = p_prob_cross.fGetValue(psSettings, &c_err);
	}//if (!c_err)
	   

	
	return(c_err);
};//CError CMoead::eConfigure(istream *psSettings)



void CMoead::vInitialize(time_t tStartTime)
{
	CBinaryOptimizer::vInitialize(tStartTime);
	t_start = tStartTime;

	CError  c_err(iERROR_PARENT_CMoeaDOptimizer);
	i_templ_length = pc_problem->pcGetEvaluation()->iGetNumberOfElements();

	if (i_templ_length <= 0)
	{
		c_err.vSetError(CMoead::iERROR_CODE_MOEAD_GENOTYPE_LEN_BELOW_0);
		return;
	}//if  (i_templ_length  <=  0)

	pc_multi_problem = (CBinaryMultiObjectiveProblem *)pc_problem->pcGetEvaluation();


	vector<double> ref(vector<double>(pc_multi_problem->iGetActiveMeasures(), 0));
	v_reference_point = ref;

	//CString  s_buf;
	//s_buf.Format("measures active = %d", pc_multi_problem->iGetActiveMeasures());
	//::Tools::vShow(s_buf);

	d_mut_prob = 1;
	d_mut_prob = d_mut_prob / i_templ_length;




	pc_log->vPrintLine("Initializing...", true);

	pc_genotype = new CBinaryCoding(i_templ_length);

	
	pc_log->vPrintLine("InitializeWeightVector...", true);
	InitializeWeightVector();

	pc_log->vPrintLine("InitializePopulation...", true);
	InitializePopulation();

	pc_log->vPrintLine("InitializeNeighborhood...", true);
	InitializeNeighborhood();


	CBinaryCoding c_genotype(i_templ_length);

	for (int ii = 0; ii < i_templ_length; ii++)
		c_genotype.piGetBits()[ii] = ::RandUtils::iRandNumber(0, 1);

	double  d_fit;
	d_fit = pc_multi_problem->dEvaluate(&c_genotype);

	b_update_best_individual(1, tStartTime, c_genotype.piGetBits(), d_fit);


	pc_log->vPrintLine("DONE...", true);



	c_time_counter.vSetStartNow();
};//void CMoead::vInitialize(time_t tStartTime)



bool CMoead::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)
{
	//return(bRunIteration_dummy(iIterationNumber, tStartTime));

	for (int iPop = 0; iPop < i_population_size; iPop++)
	{
		// select two neighboring solutions randomly
		int p1, p2;
		if (::RandUtils::dRandNumber(0, 1) < MOEAD_SIGMA)
		{
			int id1 = int(::RandUtils::dRandNumber(0, 1)*NeighborhoodSize_Tm);// NeighborhoodSize_Tm
			int id2 = int(::RandUtils::dRandNumber(0, 1)*NeighborhoodSize_Tm);
			p1 = CurrentPopulation[iPop]->IndexOfNeighbor[id1];
			p2 = CurrentPopulation[iPop]->IndexOfNeighbor[id2];
		}
		else
		{
			int id1 = int(::RandUtils::dRandNumber(0, 1)*i_population_size);// NeighborhoodSize_Tm
			int id2 = int(::RandUtils::dRandNumber(0, 1)*i_population_size);
			p1 = id1;
			p2 = id2;
		}

		// produce an offspring solution
		CMoeadIndividual offspring(this);
		offspring.SinglePointXover(*(CurrentPopulation[p1]->pc_current_solution), *(CurrentPopulation[p2]->pc_current_solution));
		offspring.GreedyRepairHeuristic(CurrentPopulation[iPop]->pc_weight_vector->lambda, pc_multi_problem->iGetActiveMeasures());

		/*
		for(int Obj=0;Obj<NumberOfKnapsacks;Obj++)
			{
				fout<<offspring.profit[Obj]<<" ";
				foutP_1<<CurrentPopulation[p1].CurrentSolution.profit[Obj]<<" ";
				foutP_2<<CurrentPopulation[p2].CurrentSolution.profit[Obj]<<" ";
			}
		fout<<"\n";foutP_1<<"\n";foutP_2<<"\n";
		*/
		// update neighboring solutions, reference point, and second population 
		//this->FindGlobalBestVectorIndex(offspring);
		this->UpdateNeighboringSolution(offspring, iPop);//GlobalBestVectorIndex
		this->UpdateReferencePoint(offspring);
		this->UpdateSecondPopulation(offspring);
	}//for (int iPop = 0; iPop < i_population_size; iPop++)


	double  d_time_passed;
	c_time_counter.bGetTimePassed(&d_time_passed);

	CString  s_buf;
	//s_buf.Format("iteration: %d  pop: %d; GenSize:%d PFsize:%d HyperVolume: %.8lf InvGenDist: %.8lf [time:%.2lf]", iIterationNumber, v_population.size(), v_non_dominated_pfs.at(0).at(0)->pc_genotype->iGetNumberOfBits(), (int)dPFQualityPointNum(), dPFQualityHyperVolume(), dPFQualityInverseGenerationalDistance(), d_time_passed);
	s_buf.Format("iteration: %d  PFsize:%d HyperVolume: %.8lf InvGenDist: %.8lf GenDist:%.8lf MaxSpread:%.8lf DominatedPF: %d [time:%.2lf] [ffe: %.0lf]", iIterationNumber, (int)dPFQualityPointNum(), dPFQualityHyperVolume(), dPFQualityInverseGenerationalDistance(), dPFQualityGenerationalDistance(), dPFQualityMaximumSpread(), (int)dPFQualityDominatedOptimalPoints(), d_time_passed, (double)pc_multi_problem->iGetFFE());
	pc_log->vPrintLine(s_buf, true);


	if (pc_multi_problem->pvGetParetoFront()->size() > 0)
	{
		double  d_fit;
		d_fit = pc_multi_problem->dEvaluate(pc_multi_problem->pvGetParetoFront()->at(0)->pcGetGenotype());

		b_update_best_individual(iIterationNumber, tStartTime, pc_multi_problem->pvGetParetoFront()->at(0)->pcGetGenotype()->piGetBits(), d_fit);
	}//if (pc_multi_problem->pvGetParetoFront()->size() > 0)

	return(true);
}//bool CMoead::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)



bool CMoead::bRunIteration_dummy(uint32_t iIterationNumber, time_t tStartTime)
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
};//bool CMoead::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)




// read predefined weight vectors from txt files
void CMoead::InitializeWeightVector()
{
	/*char filename[1024];
	sprintf(filename,"weights/wks%d%d.dat", NumberOfItems, NumberOfKnapsacks);

	ifstream indata;
	indata.open(filename); // opens the file
	if(!indata) { // file couldn't be opened
	  cerr << "Error: file could not be opened" << endl;
	  exit(1);
	}

	indata>>PopulationSize;*/

	double  d_step;
	d_step = 1;
	d_step = d_step / (i_population_size - 1);

	//CSubProblem  subp;
	CMoeadSubProblem  *pc_new_subp;
	for (int i = 1; i <= i_population_size; i++)
	{
		/*for(int iKnap=0; iKnap < pc_problem->iGetObjectiveNumber(); iKnap++)
		{
			indata>>subp.WeightVector.lambda[iKnap];
		}*/

		/*subp.WeightVector.lambda[0] = d_step * (i-1);
		subp.WeightVector.lambda[1] = 1.0 - subp.WeightVector.lambda[0];

		CurrentPopulation.push_back(subp);*/

		pc_new_subp = new CMoeadSubProblem(this);

		pc_new_subp->pc_weight_vector->lambda[0] = d_step * (i - 1);
		pc_new_subp->pc_weight_vector->lambda[1] = 1.0 - pc_new_subp->pc_weight_vector->lambda[0];

		CurrentPopulation.push_back(pc_new_subp);
	}//for(int i=1; i<= i_population_size; i++)

	//indata.close();
	NeighborhoodSize_Tm = 0.1*i_population_size;
	NeighborhoodSize_Tr = Pr * i_population_size;
}//void CMoead::InitializeWeightVector()



// initialize population with N solutions
void CMoead::InitializePopulation()
{
	for (int iPop = 0; iPop < CurrentPopulation.size(); iPop++)
	{
		CurrentPopulation[iPop]->pc_current_solution->Randomize();
		CurrentPopulation[iPop]->pc_current_solution->GreedyRepairHeuristic(CurrentPopulation[iPop]->pc_weight_vector->lambda, pc_multi_problem->iGetWeghtingType());
		UpdateReferencePoint(*(CurrentPopulation[iPop]->pc_current_solution));
	}//for(int iPop=0; iPop<CurrentPopulation.size(); iPop++)
}//void CMoead::InitializePopulation()




// compare the offspring solution with its neighhoring solutions
void CMoead::UpdateNeighboringSolution(CMoeadIndividual &offspring, int iPop)
{
	for (int n = 0; n < NeighborhoodSize_Tr; n++)
	{
		double f1, f2;
		double id = CurrentPopulation[iPop]->IndexOfNeighbor[n];    // the index of neighboring subproblem
		//f1 = offspring.ComputingFitnessValue(CurrentPopulation[id]->WeightVector.lambda, strFitnessType);  // fitness of the offspring
		f1 = offspring.ComputingFitnessValue(CurrentPopulation[id]->pc_weight_vector->lambda, pc_multi_problem->iGetWeghtingType());  // fitness of the offspring
		//f2 = CurrentPopulation[id].CurrentSolution.ComputingFitnessValue(CurrentPopulation[id].WeightVector.lambda, strFitnessType);  // fitness of neighbors
		f2 = CurrentPopulation[id]->pc_current_solution->ComputingFitnessValue(CurrentPopulation[id]->pc_weight_vector->lambda, pc_multi_problem->iGetWeghtingType());  // fitness of neighbors
		// if offspring is better, then update the neighbor
		if (f1 < f2)
		{
			*(CurrentPopulation[id]->pc_current_solution) = offspring;
		}//if(f1<f2)
	}
}//void CMoead::UpdateNeighboringSolution(CMoeadIndividual &offspring, int iPop)



// update the reference point with the best value for each objective
void CMoead::UpdateReferencePoint(CMoeadIndividual &ind)
{
	for (int iObj = 0; iObj < pc_multi_problem->iGetActiveMeasures(); iObj++)
	{
		if (ind.v_genotype[iObj] > v_reference_point[iObj])
		{
			v_reference_point[iObj] = ind.v_genotype[iObj];
		}
	}
}//void CMoead::UpdateReferencePoint(CMoeadIndividual &ind)


// determine the neighboring relationship between subproblems according to 
// the distances between weight vectors
void CMoead::InitializeNeighborhood()
{

	double MaxSize = NeighborhoodSize_Tr > NeighborhoodSize_Tm ? NeighborhoodSize_Tr : NeighborhoodSize_Tm;
	for (int iPop = 0; iPop < i_population_size; iPop++)
	{
		vector<int>    indx;
		vector<double> dist;

		for (int iPop2 = 0; iPop2 < i_population_size; iPop2++)
		{
			indx.push_back(iPop2);
			double tp = CurrentPopulation[iPop]->pc_weight_vector->DistanceTo(*(CurrentPopulation[iPop2]->pc_weight_vector));
			dist.push_back(tp);
		}

		this->MinFastSort(dist, indx, i_population_size, MaxSize);

		for (int i = 0; i < MaxSize; i++)
		{
			CurrentPopulation[iPop]->IndexOfNeighbor.push_back(indx[i]);
		}

		indx.clear();
		dist.clear();
	}
}//void CMoead::InitializeNeighborhood()


void CMoead::InitializeReferencePoint()
{

	vector<double> lambda(vector<double>(pc_multi_problem->iGetActiveMeasures(), 0));

	//for(int iObj=0; iObj<NumberOfKnapsacks; iObj++)
	for (int iObj = 0; iObj < pc_multi_problem->iGetActiveMeasures(); iObj++)
	{
		lambda[iObj] = 1.0;
		CMoeadIndividual ind(this);
		ind.Randomize();
		//ind.GreedyRepairHeuristic(lambda,"_WEIGHTEDSUM");
		ind.GreedyRepairHeuristic(lambda, MULTI_OBJECTIVE_PROBLEM_WEIGHTING_WEIGHT_VECTOR);
		this->UpdateReferencePoint(ind);
		lambda[iObj] = 0.0;
	}//for (int iObj = 0; iObj < pc_problem->iGetObjectiveNumber(); iObj++)
}//void CMOEAD::InitializeReferencePoint()



void CMoead::MinFastSort(vector<double> &dist, vector<int> &indx, int n, int m)
{
	for (int i = 0; i < m; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			if (dist[i] > dist[j])
			{
				double temp = dist[i];
				dist[i] = dist[j];
				dist[j] = temp;
				int id = indx[i];
				indx[i] = indx[j];
				indx[j] = id;
			}
		}
	}
}//void CMoead::MinFastSort(vector<double> &dist, vector<int> &indx, int n, int m)



void CMoead::UpdateSecondPopulation(CMoeadIndividual &ind)
{
	//*
	int iCount = 0;
	for (int n = 0; n < SecondPopulation.size(); n++)
	{

		if (ind == *(SecondPopulation[n]))
			return;

		if (*(SecondPopulation[n]) > ind)
			return;

		if (ind > *(SecondPopulation[n]))
		{
			SecondPopulation[n]->dominated = true;
			iCount++;
		}
		else
		{
			SecondPopulation[n]->dominated = false;
		}
	}

	for (int n = 0; n < SecondPopulation.size(); n++)
	{
		if (SecondPopulation[n]->dominated)
		{
			delete  SecondPopulation[n];
			SecondPopulation.erase(SecondPopulation.begin() + n);
			n--;
		}
	}

	CMoeadIndividual  *pc_new_ind;
	pc_new_ind = new CMoeadIndividual(this);

	*pc_new_ind = ind;
	SecondPopulation.push_back(pc_new_ind);
}//void CMOEAD::UpdateSecondPopulation(CIndividual &ind)



void CMoead::Show()
{
	for (int iPop = 0; iPop < SecondPopulation.size(); iPop++)
	{
		printf("\n");
		SecondPopulation[iPop]->Show();
	}

	printf("\n");
}//void CMoead::Show()



void CMoead::SaveSecondPopulation()
{
	char saveFilename[1024];
	int Temp = 100 * Pr;
	//sprintf(saveFilename,"POF/POF_MOEAD_KS%d%d_R%d_%d.dat", NumberOfItems, NumberOfKnapsacks, run_id,Temp);
	int  i_run_id = 1;
	sprintf(saveFilename, "POF/POF_MOEAD_KS%d%d_R%d_%d.dat", pc_multi_problem->iGetActiveMeasures(), pc_multi_problem->iGetActiveMeasures(), i_run_id, Temp);
	std::fstream fout;
	fout.open(saveFilename, std::ios::out);
	for (int n = 0; n < SecondPopulation.size(); n++)
	{
		//for(int k=0; k<NumberOfKnapsacks; k++)
		for (int k = 0; k < pc_multi_problem->iGetActiveMeasures(); k++)
		{
			fout << SecondPopulation[n]->v_objectives[k] << " ";
		}
		fout << "\n";
	}
	fout.close();
}//void CMoead::SaveSecondPopulation()




void CMoead::FindGlobalBestVectorIndex(CMoeadIndividual &ind)
{
	double MinFitness = MOEAD_IM1;                      //Initialized to very large number
	for (int BestIndex = 0; BestIndex < i_population_size; BestIndex++)
	{
		double f1 = ind.ComputingFitnessValue(CurrentPopulation[BestIndex]->pc_weight_vector->lambda, pc_multi_problem->iGetActiveMeasures());//fitness of the offspring
		if (f1 < MinFitness)
		{
			MinFitness = f1;
			GlobalBestVectorIndex = BestIndex;
		}
	}
}//void CMoead::FindGlobalBestVectorIndex(CMoeadIndividual &ind)









//---------------------------------------------------CMoeadIndividual-------------------------------


CMoeadIndividual::CMoeadIndividual(CMoead  *pcParent)
{
	pc_parent = pcParent;
	vector<int> item_vec(vector<int>(pc_parent->i_templ_length, 0));
	vector<double> obj_vec(vector<double>(pc_parent->pc_multi_problem->iGetActiveMeasures(), 0));

	v_genotype = item_vec;
	v_objectives = obj_vec;
	//weight = knap_vec;
	//profit = knap_vec;
}//CMoeadIndividual::CMoeadIndividual(CMoead  *pcParent)



CMoeadIndividual::~CMoeadIndividual()
{

}//CMoeadIndividual::~CMoeadIndividual()



// randomly initialize solution
void CMoeadIndividual::Randomize()
{

	for (int ii = 0; ii < v_genotype.size(); ii++)
	{
		if (::RandUtils::dRandNumber(0, 1) < 0.5)
			v_genotype.at(ii) = 0;
		else
			v_genotype.at(ii) = 1;
	}//for (int ii = 0; ii < v_genotype.size(); ii++)


	for (int ii = 0; ii < v_genotype.size(); ii++)
		pc_parent->pc_genotype->piGetBits()[ii] = v_genotype.at(ii);

	pc_parent->pc_multi_problem->vEvaluateParetoFront(&v_objectives, pc_parent->pc_genotype);
}//void CMoeadIndividual::Randomize()


// Greedy repair method in MOMHLib++ (Andrzej Jaszkiewicz)
//void CMoeadIndividual::GreedyRepairHeuristic(vector<double>& lambda, char* strFuncType)
void CMoeadIndividual::GreedyRepairHeuristic(vector<double>& lambda, int  iWeightingType)
{
	//we are always feasible here

}//void CMoeadIndividual::GreedyRepairHeuristic(vector<double>& lambda, int  iWeightingType)




bool CMoeadIndividual::IsFeasible()
{
	/*for (int i = 0; i < NumberOfKnapsacks; i++)
		if (weight[i] > knapsack_capacity[i])
			return false;*/
	return true;
}//bool CMoeadIndividual::IsFeasible()


void CMoeadIndividual::ComputingCapacityValue__()
{
	//not necessary
	/*int sum;
	for (int i = 0; i < NumberOfKnapsacks; i++)
	{
		sum = 0;
		for (int j = 0; j < NumberOfItems; j++)
			sum = sum + knapsack_weight[i][j] * items[j];
		weight[i] = sum;
	}*/
}//void CMoeadIndividual::ComputingCapacityValue__()



void CMoeadIndividual::ComputeObjectives()
{
	for (int ii = 0; ii < v_genotype.size(); ii++)
		pc_parent->pc_genotype->piGetBits()[ii] = v_genotype.at(ii);

	pc_parent->pc_multi_problem->vEvaluateParetoFront(&v_objectives, pc_parent->pc_genotype);
}//void CMoeadIndividual::ComputeObjectives()


bool CMoeadIndividual::operator>(const CMoeadIndividual& indiv)
{
	for (int i = 0; i < v_objectives.size(); i++)
		if (indiv.v_objectives.at(i) > v_objectives.at(i))
			return false;
	return true;
}//bool CMoeadIndividual::operator>(const CMoeadIndividual& indiv)




bool CMoeadIndividual::operator<(const CMoeadIndividual& indiv)
{
	for (int i = 0; i < v_objectives.size(); i++)
		if (indiv.v_objectives[i] < v_objectives[i])
			return false;
	return true;
}//bool CMoeadIndividual::operator<(const CMoeadIndividual& indiv)



bool CMoeadIndividual::operator==(const CMoeadIndividual& ind)
{
	for (int i = 0; i < v_objectives.size(); i++)
		if (ind.v_objectives[i] != v_objectives[i])
			return false;
	return true;
}//bool CMoeadIndividual::operator==(const CMoeadIndividual& ind)


void CMoeadIndividual::Show()
{
	std::cout << "profit: ";
	for (int i = 0; i < v_objectives.size(); i++)
		std::cout << v_objectives[i] << " ";

	std::cout << "\n";
}//void CMoeadIndividual::Show()




double CMoeadIndividual::ComputingFitnessValue(vector<double> &lambda, int  iWeightingType)
{
	double fitness = 0;

	//if (!strcmp(strFuncType, "_WEIGHTEDSUM"))
	if (iWeightingType == MULTI_OBJECTIVE_PROBLEM_WEIGHTING_WEIGHT_VECTOR)
	{
		fitness = 0;
		/*for (int i = 0; i < NumberOfKnapsacks; i++)
			fitness += lambda[i] * (ReferencePoint[i] - this->profit[i]);*/
		for (int ii = 0; ii < v_objectives.size(); ii++)
			fitness += lambda[ii] * (pc_parent->v_reference_point.at(ii) - this->v_objectives.at(ii));
	}

	//if (!strcmp(strFuncType, "_TCHEBYCHEFF"))
	if (iWeightingType == MULTI_OBJECTIVE_PROBLEM_WEIGHTING_TCHEBYCHEV)
	{
		fitness = -1e30;
		/*for (int i = 0; i < NumberOfKnapsacks; i++)
		{
			double dif = 1.1*ReferencePoint[i] - this->profit[i];
			double s = lambda[i] * (dif > 0 ? dif : -dif);
			if (s > fitness) 	fitness = s;
		}*/
		for (int ii = 0; ii < v_objectives.size(); ii++)
		{
			double dif = 1.1*pc_parent->v_reference_point.at(ii) - this->v_objectives.at(ii);
			double s = lambda[ii] * (dif > 0 ? dif : -dif);
			if (s > fitness) 	fitness = s;
		}//for (int ii = 0; ii < v_objectives.size(); ii++)
	}//if (!strcmp(strFuncType, "_TCHEBYCHEFF"))

	return fitness;
}//double CMoeadIndividual::ComputingFitnessValue(vector<double> &lambda, int  iWeightingType)


void CMoeadIndividual::SinglePointXover(CMoeadIndividual &parent1, CMoeadIndividual &parent2)
{
	//int pos = int(Rnd.GetNumber()* v_genotype.size());
	int pos = RandUtils::iRandNumber(0, v_genotype.size() - 1);


	for (int i_gene = 0; i_gene < v_genotype.size(); i_gene++)
	{
		if (i_gene < pos)
			this->v_genotype[i_gene] = parent1.v_genotype[i_gene];
		else
			this->v_genotype[i_gene] = parent2.v_genotype[i_gene];
		//if (Rnd.GetNumber() < pc_parent->d_mut_prob)
		if  (::RandUtils::dRandNumber(0, 1)  < pc_parent->d_mut_prob)
			this->v_genotype[i_gene] = 1 - this->v_genotype[i_gene];
	}

	ComputeObjectives();

}//void CMoeadIndividual::SinglePointXover(CIndividual &parent1, CIndividual &parent2)


//---------------------------------------------------CMoeadSubProblem-------------------------------

CMoeadSubProblem::CMoeadSubProblem(CMoead  *pcParent)
{
	pc_parent = pcParent;
	pc_current_solution = new CMoeadIndividual(pc_parent);
	pc_weight_vector = new CMoeadWeightVector(pc_parent->pc_multi_problem);
}//CMoeadSubProblem::CMoeadSubProblem(CMoead  *pcParent)

CMoeadSubProblem::~CMoeadSubProblem()
{
	delete  pc_current_solution;
	delete  pc_weight_vector;
}//CMoeadSubProblem::~CMoeadSubProblem()




//---------------------------------------------------CMoeadWeightVector-------------------------------

CMoeadWeightVector::CMoeadWeightVector(CBinaryMultiObjectiveProblem  *pcMultiProblem)
{
	pc_multi_problem = pcMultiProblem;

	//vector<double> w_vec(vector<double>(NumberOfKnapsacks, 0));
	vector<double> w_vec(vector<double>(pc_multi_problem->iGetActiveMeasures(), 0));
	lambda = w_vec;
}//CMoeadWeightVector::CMoeadWeightVector(CBinaryMultiObjectiveProblem  *pcMultiProblem)


CMoeadWeightVector::~CMoeadWeightVector()
{

}//CMoeadWeightVector::~CMoeadWeightVector()


double CMoeadWeightVector::DistanceTo(CMoeadWeightVector &weight)
{
	double dist = 0;
	//for(int iObj=0; iObj<NumberOfKnapsacks; iObj++)
	for (int iObj = 0; iObj < pc_multi_problem->iGetActiveMeasures(); iObj++)
	{
		double diff = this->lambda[iObj] - weight.lambda[iObj];
		dist += diff * diff;
	}
	return sqrt(dist);
}//double CMoeadWeightVector::DistanceTo(CMoeadWeightVector &weight)


void CMoeadWeightVector::Show()
{
	//for(int iObj=0; iObj<NumberOfKnapsacks; iObj++)
	for (int iObj = 0; iObj < pc_multi_problem->iGetActiveMeasures(); iObj++)
	{
		printf(" %f ", lambda[iObj]);
	}
	printf("\n");
}//void CMoeadWeightVector::Show()
