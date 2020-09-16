#include "Optimizer.h"

#include "BinaryCoding.h"
//#include "GenePatternReplacementUtils.h"
#include "GenerationUtils.h"
//#include "LinkageUtils.h"
//#include "PermutationCoding.h"
//#include "RealCoding.h"
#include "StopConditionUtils.h"
#include "UIntCommandParam.h"


template <class TGenotype, class TFenotype>
COptimizer<TGenotype, TFenotype>::COptimizer(CProblem<TGenotype, TFenotype> *pcProblem, CLog *pcLog, uint32_t iRandomSeed)
{
	pc_problem = pcProblem;

	pc_stop_condition = nullptr;

//	pl_gene_patterns = nullptr;
//	pc_linkage = nullptr;

	pc_log = pcLog;

	pc_empty_generation = GenerationUtils::pcGetEmptyGeneration(pcProblem);
	pc_best_individual = nullptr;

	i_random_seed = iRandomSeed;

	b_own_params = true;
	b_own_gene_patterns = true;

//	pc_linkage_analyzer = NULL;
}//COptimizer<TGenotype, TFenotype>::COptimizer(CProblem<TGenotype, TFenotype> *pcProblem, CLog *pcLog, uint32_t iRandomSeed)

template <class TGenotype, class TFenotype>
COptimizer<TGenotype, TFenotype>::COptimizer(COptimizer<TGenotype, TFenotype> *pcOther)
{
	pc_problem = pcOther->pc_problem;

	pc_stop_condition = pcOther->pc_stop_condition;

	i_gene_pattern_min_length = pcOther->i_gene_pattern_min_length;
	i_linkage_frequency = pcOther->i_linkage_frequency;

	//pl_gene_patterns = pcOther->pl_gene_patterns;
	//pc_linkage = pcOther->pc_linkage;

	pc_log = pcOther->pc_log;

	pc_empty_generation = pcOther->pc_empty_generation;
	pc_best_individual = nullptr;

	i_random_seed = pcOther->i_random_seed;

//	pc_linkage_analyzer = pcOther->pc_linkage_analyzer;

	b_own_params = false;
	b_own_gene_patterns = false;
}//COptimizer<TGenotype, TFenotype>::COptimizer(COptimizer<TGenotype, TFenotype> *pcOther)

template <class TGenotype, class TFenotype>
COptimizer<TGenotype, TFenotype>::~COptimizer()
{
	v_clear_params();

	if (b_own_params)
	{
		delete pc_empty_generation;
	}//if (b_own_params)

	vResetBestIndividual();


/*	if (pc_linkage_analyzer != NULL)
	{
		//vReportLinkage();

		vector<vector <int>>  v_dependent_genes;
		pc_problem->vGetTrueGeneDependencies(&v_dependent_genes);
		pc_linkage_analyzer->vSetTrueLinkage(&v_dependent_genes);

		if  (b_linkage_improver == true)
			pc_linkage_analyzer->vSaveLinkageReport(pc_log->sGetLogFile(), "Improver");
		else
			pc_linkage_analyzer->vSaveLinkageReport(pc_log->sGetLogFile(), "");
	}//if (pc_linkage_analyzer != NULL)*/
}//COptimizer::~COptimizer()



template <class TGenotype, class TFenotype>
CError COptimizer<TGenotype, TFenotype>::eConfigure(istream *psSettings)
{
	CError c_error;

	v_clear_params();

	pc_stop_condition = StopConditionUtils::pcGetStopCondition<TGenotype, TFenotype>(pc_problem->pcGetEvaluation(), psSettings, &c_error);


	if (!c_error)
	{
		CUIntCommandParam p_max_number_of_gene_patterns(OPTIMIZER_ARGUMENT_MAX_NUMBER_OF_GENE_PATTERNS, false);
		uint32_t i_max_number_of_gene_patterns = p_max_number_of_gene_patterns.iGetValue(psSettings, &c_error);

		/*if (!c_error && p_max_number_of_gene_patterns.bHasValue())
		{
			CGenePatternReplacement *pc_gene_pattern_replacement = GenePatternReplacementUtils::pcGetGenePatternReplacement(psSettings, &c_error);
			pl_gene_patterns = new CGenePatternList(i_max_number_of_gene_patterns, pc_gene_pattern_replacement);

			if (!c_error && pl_gene_patterns)
			{
				CUIntCommandParam p_gene_pattern_min_length(OPTIMIZER_ARGUMENT_GENE_PATTERN_MIN_LENGTH, 1, UINT32_MAX);
				i_gene_pattern_min_length = p_gene_pattern_min_length.iGetValue(psSettings, &c_error);
			}//if (!c_error && pl_gene_patterns)
		}//if (!c_error && p_max_number_of_gene_patterns.bHasValue())*/
	}//if (!c_error)

	/*if (!c_error)
	{
		pc_linkage = LinkageUtils::pcGetLinkage<TGenotype>(psSettings, &c_error, false);

		if (!c_error && pc_linkage)
		{
			CUIntCommandParam p_linkage_frequency(OPTIMIZER_ARGUMENT_LINKAGE_FREQUENCY, 1, UINT32_MAX);
			i_linkage_frequency = p_linkage_frequency.iGetValue(psSettings, &c_error);
		}//if (!c_error && pc_linkage)
	}//if (!c_error)


	int  i_linkage_quality = 0;
	if (!c_error)
	{
		CUIntCommandParam p_linkage_quality(OPTIMIZER_ARGUMENT_LINKAGE_ANALYZER, false);
		i_linkage_quality = p_linkage_quality.iGetValue(psSettings, &c_error);

		if (p_linkage_quality.bHasValue() == false)  i_linkage_quality = 0;

		if (i_linkage_quality > 0)
		{
			pc_linkage_analyzer = new CLinkageAnalyzer();//new CLinkageAnalyzer(pc_log);
			pc_linkage_analyzer->eConfigure(psSettings);

			CUIntCommandParam p_linkage_quality(OPTIMIZER_ARGUMENT_LINKAGE_ANALYZER, false);
			i_linkage_quality = p_linkage_quality.iGetValue(psSettings, &c_error);

		}//if (i_linkage_quality > 0)
	}//if (!c_error)*/

	return c_error;
}//CError COptimizer<TGenotype, TFenotype>::eConfigure(istream *psSettings)



template <class TGenotype, class TFenotype>
CError COptimizer<TGenotype, TFenotype>::eConfigure3LO(istream *psSettings)
{
	CError c_err;

	v_clear_params();

	pc_stop_condition = StopConditionUtils::pcGetStopCondition<TGenotype, TFenotype>(pc_problem->pcGetEvaluation(), psSettings, &c_err);

	return(c_err);
}//template <class TGenotype, class TFenotype>


template <class TGenotype, class TFenotype>
CError COptimizer<TGenotype, TFenotype>::eConfigureLinkageImprover(istream *psSettings)
{
	CError c_err;

	v_clear_params();

	pc_stop_condition = StopConditionUtils::pcGetStopCondition<TGenotype, TFenotype>(pc_problem->pcGetEvaluation(), psSettings, &c_err);

	return(c_err);
}//template <class TGenotype, class TFenotype>


template <class TGenotype, class TFenotype>
CError COptimizer<TGenotype, TFenotype>::eConfigureNSGA2(istream *psSettings)
{
	CError c_err;

	v_clear_params();

	pc_stop_condition = StopConditionUtils::pcGetStopCondition<TGenotype, TFenotype>(pc_problem->pcGetEvaluation(), psSettings, &c_err);

	return(c_err);
}//template <class TGenotype, class TFenotype>



template <class TGenotype, class TFenotype>
CError COptimizer<TGenotype, TFenotype>::eConfigureMoGomea(istream *psSettings)
{
	CError c_err;

	v_clear_params();

	pc_stop_condition = StopConditionUtils::pcGetStopCondition<TGenotype, TFenotype>(pc_problem->pcGetEvaluation(), psSettings, &c_err);

	return(c_err);
}//template <class TGenotype, class TFenotype>



template <class TGenotype, class TFenotype>
CError COptimizer<TGenotype, TFenotype>::eConfigureMoead(istream *psSettings)
{
	CError c_err;

	v_clear_params();

	pc_stop_condition = StopConditionUtils::pcGetStopCondition<TGenotype, TFenotype>(pc_problem->pcGetEvaluation(), psSettings, &c_err);

	return(c_err);
}//template <class TGenotype, class TFenotype>





template <class TGenotype, class TFenotype>
void COptimizer<TGenotype, TFenotype>::vInitialize(time_t tStartTime)
{
	vResetBestIndividual();

	/*if (pc_linkage_analyzer != NULL)
	{
		vector<vector <int>>  v_dependent_genes;
		pc_problem->vGetTrueGeneDependencies(&v_dependent_genes);
		pc_linkage_analyzer->vSetTrueLinkage(&v_dependent_genes);
		//pc_linkage_analyzer->vSaveLinkageReport(pc_log->sGetLogFile());
	}//if (pc_linkage_analyzer != NULL)*/
}//void COptimizer<TGenotype, TFenotype>::vInitialize(time_t tStartTime)


template <class TGenotype, class TFenotype>
void COptimizer<TGenotype, TFenotype>::vRun()
{
	CTimeCounter  c_time_counter;
	time_t t_start_time = time(nullptr);

	CString  s_buf;
	
	double  d_best_fitness;
	uint32_t i_iteration_number;
	double  d_time;
	uint64_t i_ffe;

	vInitialize(t_start_time);
	c_time_counter.vSetStartNow();

	i_iteration_number = 0;

	while (!pc_stop_condition->bStop(t_start_time, i_iteration_number, pc_problem->pcGetEvaluation()->iGetFFE(), pc_best_individual))
	{
		bRunIteration(i_iteration_number, t_start_time);
		i_iteration_number++;

		d_best_fitness = pc_best_individual->dGetFitnessValue();
		c_time_counter.bGetTimePassed(&d_time);
		i_ffe = pc_problem->pcGetEvaluation()->iGetFFE();

		CString log_message;
		log_message.AppendFormat("[PRW LOG] best fitness: \t %.8lf \t ffe: \t %u \t time: \t %.4lf", d_best_fitness, i_ffe, d_time);
		//pc_log->vPrintLine(log_message, true, PRW_LOG_SYSTEM);

		if (pc_problem->pcGetEvaluation()->bMultiObjective() == true)
		{
			pc_log->vPrintLine(pc_problem->pcGetEvaluation()->sMultiObjectiveReportIter(), true);
			//CBinaryMultiObjectiveProblem *pc_problem_multi;
			//pc_problem_multi = (CBinaryMultiObjectiveProblem *) pc_problem->pcGetEvaluation();

			//return(pc_problem_multi->dPFQualityInverseGenerationalDistance());

			//s_buf.Format("iteration: %d  PFsize:%d HyperVolume: %.8lf InvGenDist: %.8lf GenDist:%.8lf MaxSpread:%.8lf DominatedPF: %d [time:%.2lf] [ffe: %.0lf]", i_iteration_number, (int)dPFQualityPointNum(), dPFQualityHyperVolume(), dPFQualityInverseGenerationalDistance(), dPFQualityGenerationalDistance(), dPFQualityMaximumSpread(), (int)dPFQualityDominatedOptimalPoints(), d_time_passed, (double)pc_multi_problem->iGetFFE());
			//  ///pc_log->vPrintLine(s_buf, true);
		}//if (pc_problem->pcGetEvaluation()->bMultiObjective() == true)

	}//while (!pc_stop_condition->bStop(t_start_time, i_iteration_number, pc_problem->pcGetEvaluation()->iGetFFE(), pc_best_individual))

	vExecuteBeforeEnd();

	if (pc_problem->pcGetEvaluation()->bMultiObjective() == true)
	{
		//pc_log->vPrintLine("PARETO FRONT:", true, LOG_SYSTEM_PARETO_FRONT);

		vector<CString>  v_pf_report;
		pc_problem->pcGetEvaluation()->vReportPF(&v_pf_report);

		for (int ii = 0; ii < v_pf_report.size(); ii++)
		{
			pc_log->vPrintLine(v_pf_report.at(ii), false, LOG_SYSTEM_PARETO_FRONT);
		}//for (int ii = 0; ii < v_pf_report.size(); ii++)				
	}//if (pc_problem->pcGetEvaluation()->bMultiObjective() == true)
	
}//void COptimizer<TGenotype, TFenotype>::vRun()


template <class TGenotype, class TFenotype>
CString  COptimizer<TGenotype, TFenotype>::sGetLogName()
{
	if (pc_log == NULL)  return("<no log>");

	return(pc_log->sGetLogFile());
}//CString  COptimizer<TGenotype, TFenotype>::sGetLogName()





template <class TGenotype, class TFenotype>
void COptimizer<TGenotype, TFenotype>::vResetBestIndividual()
{
	delete pc_best_individual;
	pc_best_individual = nullptr;
}//void COptimizer<TGenotype, TFenotype>::vResetBestIndividual()



template <class TGenotype, class TFenotype>
void COptimizer<TGenotype, TFenotype>::vSetBestIndividual(CIndividual<TGenotype, TFenotype> *pcBestIndividual, bool bCopy)
{
	vResetBestIndividual();
	pc_best_individual = bCopy ? new CIndividual<TGenotype, TFenotype>(pcBestIndividual) : pcBestIndividual;
}//void COptimizer<TGenotype, TFenotype>::vSetBestIndividual(CIndividual<TGenotype, TFenotype> *pcBestIndividual, bool bCopy)

template <class TGenotype, class TFenotype>
bool COptimizer<TGenotype, TFenotype>::b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime, CIndividual<TGenotype, TFenotype> *pcIndividual)
{
	bool b_updated = false;

	if (!pc_best_individual || pc_problem->bIsBetterIndividual(pcIndividual, pc_best_individual))
	{
		delete pc_best_individual;
		pc_best_individual = new CIndividual<TGenotype, TFenotype>(pcIndividual);

		v_update_statistics_of_best(iIterationNumber, tStartTime);

		b_updated = true;
	}//if (!pc_best_individual || pc_problem->bIsBetterIndividual(pcIndividual, pc_best_individual))

	return b_updated;
}//bool COptimizer<TGenotype, TFenotype>::b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime, CIndividual<TGenotype, TFenotype> *pcIndividual)

template <class TGenotype, class TFenotype>
bool COptimizer<TGenotype, TFenotype>::b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime, double dCurrentBestFitnessValue, function<void(TGenotype*)> &&fUpdateBestGenotype)
{
	bool b_updated = false;

	if (!pc_best_individual || pc_problem->bIsBetterFitnessValue(dCurrentBestFitnessValue, pc_best_individual->dGetFitnessValue()))
	{
		if (!pc_best_individual)
		{
			pc_best_individual = pc_create_individual(pc_empty_generation->pcGenerateEmpty());
		}//if (!pc_best_individual)

		fUpdateBestGenotype(pc_best_individual->pcGetGenotype());
		pc_best_individual->vSetFitnessValue(dCurrentBestFitnessValue);

		v_update_statistics_of_best(iIterationNumber, tStartTime);

		b_updated = true;
	}//if (!pc_best_individual || pc_problem->bIsBetterFitnessValue(dCurrentBestFitnessValue, pc_best_individual->dGetFitnessValue()))

	return b_updated;
}//bool COptimizer<TGenotype, TFenotype>::b_update_best_individual(uint32_t iIterationNumber, time_t tStartTime, double dCurrentBestFitnessValue, function<void(TGenotype*)> &&fUpdateBestGenotype)

template <class TGenotype, class TFenotype>
void COptimizer<TGenotype, TFenotype>::v_clear_params()
{
	if (b_own_params)
	{
		delete pc_stop_condition;
		pc_stop_condition = nullptr;

		//delete pc_linkage;
		//pc_linkage = nullptr;
	}//if (b_own_params)

	if (b_own_gene_patterns)
	{
		//delete pl_gene_patterns;
		//pl_gene_patterns = nullptr;
	}//if (b_own_gene_patterns)
}//void COptimizer<TGenotype, TFenotype>::v_clear_params()

template <class TGenotype, class TFenotype>
void COptimizer<TGenotype, TFenotype>::v_update_statistics_of_best(uint32_t iIterationNumber, time_t tStartTime)
{
	t_best_time = time(nullptr) - tStartTime;
	i_best_ffe = pc_problem->pcGetEvaluation()->iGetFFE();
}//void COptimizer<TGenotype, TFenotype>::v_update_statistics_of_best(uint32_t iIterationNumber, time_t tStartTime)

template class COptimizer<CBinaryCoding, CBinaryCoding>;
//template class COptimizer<CRealCoding, CPermutationCoding>;
//template class COptimizer<CRealCoding, CRealCoding>;