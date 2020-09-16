#pragma once

#ifndef BINARY_EVALUATION_MULTI_OBJECTIVE_H
#define BINARY_EVALUATION_MULTI_OBJECTIVE_H


#include "BinaryCoding.h"
#include "CompProblem.h"
#include "Error.h"
#include "Evaluation.h"
#include "BinaryEvaluation.h"
#include "FiberNet.h"

#include "../DSMGA2/chromosome.h"
#include "../P3/Evaluation.h"

#include <cstdint>
#include <iostream>
#include <vector>

using namespace FiberNets;
using namespace ProblemTools;
using namespace std;




#define EVALUATION_ARGUMENT_BINARY_MULTI_WEIGHTING "multi_objective_weight_type"
#define MULTI_OBJECTIVE_PROBLEM_WEIGHTING_WEIGHT_VECTOR_TEXT "weight_vector"
#define MULTI_OBJECTIVE_PROBLEM_WEIGHTING_TCHEBYCHEV_TEXT "tchebychev_vector"


#define EVALUATION_ARGUMENT_BINARY_MULTI_OPTIMAL_PF_BIT_STEP "optimal_pf_bit_step"
#define EVALUATION_ARGUMENT_BINARY_MULTI_OPTIMAL_PF_NOT_OPTIMAL "optimal_pf_not_really_optimal"
#define VERY_LARGE_VALUE    999999998888

#define MULTI_OBJECTIVE_PROBLEM_WEIGHTING_WEIGHT_VECTOR    0
#define MULTI_OBJECTIVE_PROBLEM_WEIGHTING_TCHEBYCHEV    1
#define MULTI_OBJECTIVE_PROBLEM_WEIGHTING_TCHEBYCHEV_MULTIPLIER    10


class  CMultiIndividual;
class  CMultiObjectiveMeasure;

class CBinaryMultiObjectiveProblem : public CBinaryFileConfigEvaluation
{
public:
	CBinaryMultiObjectiveProblem() { b_theoretical_pareto_front_really_optimal = true; i_last_pareto_front_update_ffe = 0; d_last_pareto_front_update_time = 0; i_weghting_type = MULTI_OBJECTIVE_PROBLEM_WEIGHTING_WEIGHT_VECTOR;};
	CBinaryMultiObjectiveProblem(const CBinaryMultiObjectiveProblem &pcOther);
	virtual ~CBinaryMultiObjectiveProblem();

	int  iGetActiveMeasures() { return(v_measures.size()); };
	vector<CMultiObjectiveMeasure*>  *pvGetMeasures() { return(&v_measures); };
	int  iGetWeghtingType() { return(i_weghting_type); };
	void  vSetWeghtingType(int  iNewWeghtingType) { i_weghting_type = iNewWeghtingType; };
	virtual void  vEvaluateParetoFront(vector<double> *pvPF, CBinaryCoding *pcFenotype);

	//PF quality measures
	double  dPFQualityHyperVolume();
	double  dPFQualityPointNum();
	double	dPFQualityOptimalPointsPerc();
	int		iPFQualityOptimalPointsFound();
	int		iPFQualityOptimalPointsNum() { return(v_theoretical_optimal_pareto_front.size()); };
	double  dPFQualityInverseGenerationalDistance();
	double  dPFQualityGenerationalDistance();
	double  dPFQualityDominatedOptimalPoints();
	double  dPFQualityMaximumSpread();

	void  vFlushGlobalPareto();

	vector<CMultiIndividual  *>  *pvGetParetoFront() { return(&v_global_pareto_front); }

	virtual CString  sAdditionalSummaryInfo();
	virtual bool  bMultiObjective() { return(true); }
	virtual void  vReportPF(vector<CString> *pvPFReport);
	virtual bool bOptimalFound() { if ((dPFQualityInverseGenerationalDistance() == 0) && (b_theoretical_pareto_front_really_optimal == true))  return(true); return(false); }
	//virtual bool bOptimalFound() { if ((dPFQualityInverseGenerationalDistance() == 0) && (b_theoretical_pareto_front_really_optimal == true) && (v_theoretical_optimal_pareto_front.size() > 0))  return(true); return(false); }
	virtual CString  sMultiObjectiveReportIter();

	CError eConfigure(istream *psSettings);

	vector<CMultiIndividual  *>  *pvGetGlobalParetoFront() { return(&v_global_pareto_front); };
protected:

	int  i_weghting_type;
	vector<CMultiObjectiveMeasure *>  v_measures;
	virtual void  v_prepare_solution(CBinaryCoding *pcFenotype) = 0;//prepares constructs solution for rating by measures - dependent on the problem
	virtual double d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift);
	void  v_evaluate_pareto_front(vector<double> *pvPF, CBinaryCoding *pcFenotype);

	int   i_join_the_global_pareto(CMultiIndividual *pcInd);
	void  v_add_to_non_dominated_pf(vector<CMultiIndividual *>  *pvFrontToFill, CMultiIndividual  *pcInd);



	vector<CMultiIndividual  *>  v_global_pareto_front;
	vector<CMultiIndividual  *>  v_theoretical_optimal_pareto_front;//by default it is empty. problems can fill it during configuration

	
	TimeCounters::CTimeCounter  c_time_counter;
	uint64_t  i_last_pareto_front_update_ffe;
	double  d_last_pareto_front_update_time;
	vector<double>  v_pf_buffer;
	bool  b_theoretical_pareto_front_really_optimal;

};//class CBinaryMultiObjectiveProblem



class  CMultiIndividual
{
public:
	CMultiIndividual() { pc_genotype = NULL; };
	virtual ~CMultiIndividual();

	virtual void  vRate();
	vector<double>  *pvGetFitness();
	bool  bFrontsDiffer(CMultiIndividual  *pcOther);

	int  iCheckDomination(CMultiIndividual  *pcOther);

	virtual CMultiIndividual  *pcClone() = 0;

	CString  sReport();
	CString  sReportObj();

	CBinaryCoding *pcGetGenotype() { return(pc_genotype); }


	double  dGetObjectiveSpaceDistance(CMultiIndividual  *pcOther);

protected:
	vector<double>  v_fitness;
	bool  b_rated;

	CBinaryMultiObjectiveProblem *pc_problem;
	CBinaryCoding *pc_genotype;

};//class  CMultiIndividual



class  CMultiIndividualPFpoints : public CMultiIndividual
{
public:
	CMultiIndividualPFpoints() { pc_genotype = NULL; b_matched_by_pf = false; pc_problem = NULL; };
	virtual ~CMultiIndividualPFpoints();

	virtual CMultiIndividual  *pcClone();
	void  vConfigure(CBinaryMultiObjectiveProblem *pcProblem) { pc_problem = pcProblem; }
	void  vConfigure(CBinaryCoding *pcRatedSolution,  bool  bOwnSolution, vector<double>  *pvPfBuffer, bool  bRated, CBinaryMultiObjectiveProblem *pcProblem);
	double  dGetInversedGenerationalDistance(vector  <CMultiIndividual *>  *pvGlobalParetoFront, bool  bCheckDominance = true);

	void  vRepDominance(CMultiIndividual *pcDominator);

	bool  bLoadFromFile_MaxCutMoGomea(FILE  *pfSource);
	bool bLoadFromFilePRW(FILE  *pfSource);

	
protected:
	bool  b_own_solution;
	bool  b_matched_by_pf;
};//class  CMultiIndividualPFpoints : public CMultiIndividual






class  CMultiObjectiveMeasure
{
public:
	CMultiObjectiveMeasure(CBinaryMultiObjectiveProblem  *pcParentProblem) { pc_parent_problem = pcParentProblem; };
	
	virtual double  dGetMax() = 0;
	virtual double  dEvaluate(CBinaryCoding *pcFenotype) = 0;
	virtual CMultiObjectiveMeasure  *pcClone() = 0;
	CString   sGetName() { return(s_name); }

	double  dWeight;

protected:

	CString  s_name;
	CBinaryMultiObjectiveProblem  *pc_parent_problem;
};//class  CMultiObjectiveMeasure




#define EVALUATION_ARGUMENT_REVERSE_BIT_ORDER "reverse_bit_for_0_measure"

#define  s_MULTI_REVERSING_MEASURE_1s    "MEASURE_1s"
class  CBinaryMultiReversingMeasure1s : public CMultiObjectiveMeasure
{
public:
	CBinaryMultiReversingMeasure1s(CBinaryMultiObjectiveProblem  *pcParentProblem, double  dMax) : CMultiObjectiveMeasure(pcParentProblem)
	{
		s_name = s_MULTI_REVERSING_MEASURE_1s;
		d_max = dMax;
	}

	virtual double  dGetMax() { return(d_max); };
	double  dEvaluate(CBinaryCoding *pcFenotype);
	CMultiObjectiveMeasure  *pcClone() { return(new CBinaryMultiReversingMeasure1s(*this)); };

private:
	double  d_max;

};//class  CBinaryMultiReversingMeasure1s : public CMultiObjectiveMeasure


#define  s_MULTI_REVERSING_MEASURE_0s    "MEASURE_0s"
class  CBinaryMultiReversingMeasure0s : public CMultiObjectiveMeasure
{
public:
	CBinaryMultiReversingMeasure0s(CBinaryMultiObjectiveProblem  *pcParentProblem, double  dMax);
	CBinaryMultiReversingMeasure0s(const CBinaryMultiReversingMeasure0s &pcOther);

	~CBinaryMultiReversingMeasure0s();

	virtual double  dGetMax() { return(d_max); };
	double  dEvaluate(CBinaryCoding *pcFenotype);
	CMultiObjectiveMeasure  *pcClone() { return(new CBinaryMultiReversingMeasure0s(*this)); };


	void  vSetReversingBitOrder(bool  bReverseBitOrder) {b_reverse_bit_order_for_0_measure = bReverseBitOrder;}
private:
	bool  b_reverse_bit_order_for_0_measure;
	CBinaryCoding *pc_phenotype_inner;
	double  d_max;
};//class  CBinaryMultiReversingMeasure1s : public CMultiObjectiveMeasure



class CBinaryMultiReversing : public CBinaryMultiObjectiveProblem
{
friend class CBinaryMultiReversingMeasure1s;
friend class CBinaryMultiReversingMeasure0s;

public:
	CBinaryMultiReversing();
	CBinaryMultiReversing(FILE *pfConfig, CError *pcError);
	CBinaryMultiReversing(const CBinaryMultiReversing &pcOther);

	virtual ~CBinaryMultiReversing();

	CError eReport(FILE *pfReport);

	CError  eSave(CString  sDest);
	CError  eSave(FILE  *pfDest);

	//mutliobjective methods
	//void  vEvaluateParetoFront(vector<double> *pvPF, CBinaryCoding *pcFenotype);//inherited is enough

	virtual CError eConfigure(istream *psSettings);
	void  vFillOptimalPf(int  iCreateOptimalPfBitStep);


protected:
	virtual CError e_init(FILE *pfConfig);
	CError e_load_settings_from_file(FILE *pfConfig) {};


	CEvaluation<CBinaryCoding> *pc_evaluation_inner;


	double d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift) { return(CBinaryMultiObjectiveProblem::d_evaluate(pcFenotype, iShift)); };
	void  v_prepare_solution(CBinaryCoding *pcFenotype) {};

};//class CBinaryMultiReversing : public CBinaryMultiObjectiveProblem



#define EVALUATION_ARGUMENT_BINARY_MULTI_MAXCUT_MOGOMEA_ARG_NUM "maxcut_mogomea_arg_num"
#define EVALUATION_ARGUMENT_BINARY_MULTI_MAXCUT_OPTIMAL_PF_FILE  "maxcut_mogomea_opt_pf"

#define EVALUATION_ARGUMENT_BINARY_MULTI_MAXCUT_MOGOMEA_FILE_TEMPLATE "maxcut\\maxcut_instance_%d_%d.txt"


#define  s_MULTI_MAXCUT_MEASURE    "MAXCUT_Measure"
class CBinaryMultiMaxcutMoGomeaMeasure;

class  CBinaryMultiMaxcutMoGomeaMeasureEdge
{
	friend class CBinaryMultiMaxcutMoGomeaMeasure;
public:
	CBinaryMultiMaxcutMoGomeaMeasureEdge();

	CError  eLoad(CString  sSource);
	CString  sToString();


private:
	int  i_start_node;
	int  i_end_node;
	double  d_weight;
};//class  CBinaryMultiMaxcutMoGomeaMeasureEdge


class  CBinaryMultiMaxcutMoGomeaMeasure : public CMultiObjectiveMeasure
{
public:
	CBinaryMultiMaxcutMoGomeaMeasure(CBinaryMultiObjectiveProblem  *pcParentProblem, double  dMax);
	

	virtual ~CBinaryMultiMaxcutMoGomeaMeasure();

	virtual double  dGetMax() { return(d_max); };
	double  dEvaluate(CBinaryCoding *pcFenotype);
	CMultiObjectiveMeasure  *pcClone() { return(new CBinaryMultiMaxcutMoGomeaMeasure(*this)); };

	CError  eLoadFromFile(FILE  *pfSource);
	CError  eSaveToFile(FILE  *pfSource);
	

private:
	double  d_max;
	vector<CBinaryMultiMaxcutMoGomeaMeasureEdge>  v_edges;

};//class  CBinaryMultiMaxcutMoGomeaMeasure : public CMultiObjectiveMeasure


class CBinaryMultiMaxcutMoGomea : public CBinaryMultiObjectiveProblem
{
public:
	CBinaryMultiMaxcutMoGomea();
	CBinaryMultiMaxcutMoGomea(FILE *pfConfig, CError *pcError);
	CBinaryMultiMaxcutMoGomea(const CBinaryMultiMaxcutMoGomea &pcOther);

	virtual ~CBinaryMultiMaxcutMoGomea();

	
	CError eReport(FILE *pfReport);

	CError  eSave(CString  sDest);
	CError  eSave(FILE  *pfDest);

	//mutliobjective methods
	//void  vEvaluateParetoFront(vector<double> *pvPF, CBinaryCoding *pcFenotype);//inherited is enough

	virtual CError eConfigure(istream *psSettings);
	CError  eLoadPFFromFile(FILE  *pfSource);

protected:
	virtual CError e_init(FILE *pfConfig) { CError  c_err; return(c_err); };

	double d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift) { return(CBinaryMultiObjectiveProblem::d_evaluate(pcFenotype, iShift)); };
	void  v_prepare_solution(CBinaryCoding *pcFenotype) {};

};//class CBinaryMultiMaxcutMoGomea : public CBinaryMultiObjectiveProblem





#define EVALUATION_ARGUMENT_BINARY_MULTI_KNAPSACK_MOGOMEA_ARG_NUM "knapsack_mogomea_arg_num"
#define EVALUATION_ARGUMENT_BINARY_MULTI_KNAPSACK_OPTIMAL_PF_FILE  "knapsack_mogomea_opt_pf"

#define EVALUATION_ARGUMENT_BINARY_MULTI_KNAPSACK_MOGOMEA_FILE_TEMPLATE "knapsack\\knapsack.%d.2.txt"


#define  s_MULTI_KNAPSACK_MEASURE    "KNAPSACK_Measure"

class CBinaryMultiKnapsackMoGomeaMeasure;

class  CBinaryMultiKnapsackMoGomeaMeasureItemData
{
	friend class CBinaryMultiKnapsackMoGomeaMeasure;
public:
	CBinaryMultiKnapsackMoGomeaMeasureItemData() { d_weight = 0; d_profit = 0; };

	CError  eLoad(FILE  *pfSource, int iItemOffset);
	CError  eSave(FILE  *pfDest);

	double  dGetRatio() { return(d_profit / d_weight); }

private:
	int  i_item_offset;
	double  d_weight;
	double  d_profit;
};//class  CBinaryMultiKnapsackMoGomeaMeasureItemData


class  CBinaryMultiKnapsackMoGomeaMeasureItemDataRatio
{
public:
	double  d_ratio;
	int i_item_offset;


	bool  operator<(CBinaryMultiKnapsackMoGomeaMeasureItemDataRatio &pcOther) const { return(d_ratio < pcOther.d_ratio); }
};//class  CBinaryMultiKnapsackMoGomeaMeasureItemDataRatio



class CBinaryMultiKnapsackMoGomea;
class  CBinaryMultiKnapsackMoGomeaMeasure : public CMultiObjectiveMeasure
{
	friend class CBinaryMultiKnapsackMoGomea;
public:
	CBinaryMultiKnapsackMoGomeaMeasure(CBinaryMultiObjectiveProblem  *pcParentProblem, double  dMax);


	virtual ~CBinaryMultiKnapsackMoGomeaMeasure();

	virtual double  dGetMax() { return(d_max); };
	double  dEvaluate(CBinaryCoding *pcFenotype);
	void  vConstraintSatisfied(CBinaryCoding *pcFenotype, vector<CBinaryMultiKnapsackMoGomeaMeasureItemDataRatio>  *pvItemsRatioOrder);

	CMultiObjectiveMeasure  *pcClone() { return(new CBinaryMultiKnapsackMoGomeaMeasure(*this)); };

	CError  eLoadFromFile(FILE  *pfSource, int  iArgNumber);
	CError  eSave(FILE  *pfSource);


private:
	double  d_max;

	double  d_capacity;
	vector<CBinaryMultiKnapsackMoGomeaMeasureItemData>  v_items;

};//class  CBinaryMultiKnapsackMoGomeaMeasure : public CMultiObjectiveMeasure





class CBinaryMultiKnapsackMoGomea : public CBinaryMultiObjectiveProblem
{
public:
	CBinaryMultiKnapsackMoGomea();
	CBinaryMultiKnapsackMoGomea(FILE *pfConfig, CError *pcError);
	CBinaryMultiKnapsackMoGomea(const CBinaryMultiKnapsackMoGomea &pcOther);

	virtual ~CBinaryMultiKnapsackMoGomea();


	
	//mutliobjective methods
	void  vEvaluateParetoFront(vector<double> *pvPF, CBinaryCoding *pcFenotype);//inherited is enough

	virtual CError eConfigure(istream *psSettings);
	CError  eLoadFromFile(FILE  *pfSource);
	CError  eSave(FILE  *pfSource);

	CError  eLoadPFFromFile(FILE  *pfSource);

protected:
	virtual CError e_init(FILE *pfConfig) { CError  c_err; return(c_err); };

	double d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift) { return(CBinaryMultiObjectiveProblem::d_evaluate(pcFenotype, iShift)); };
	void  v_prepare_solution(CBinaryCoding *pcFenotype);

	vector<CBinaryMultiKnapsackMoGomeaMeasureItemDataRatio>  v_items_ratio_order;

	CBinaryCoding *pc_solution_buf;
};//class CBinaryMultiMaxcutMoGomea : public CBinaryMultiObjectiveProblem









//--------------------------------------PAINTS PROBLEM--------------------------------------------------------------

class  CBinaryMultiPaintsMeasureMakespan : public CMultiObjectiveMeasure
{
#define  s_MULTI_PAINTS_MEASURE_MAKESPAN    "MEASURE_Makespan"
public:
	CBinaryMultiPaintsMeasureMakespan(CBinaryMultiObjectiveProblem  *pcParentProblem) : CMultiObjectiveMeasure(pcParentProblem)
	{s_name = s_MULTI_PAINTS_MEASURE_MAKESPAN;}

	virtual double  dGetMax() { return(1); };
	double  dEvaluate(CBinaryCoding *pcFenotype);
	CMultiObjectiveMeasure  *pcClone() { return(new CBinaryMultiPaintsMeasureMakespan(*this)); };

private:
	vector<double>  v_machines_makespan_buffer;

};//class  CBinaryMultiPaintsMeasureMakespan



class  CBinaryMultiPaintsMeasureOverhead : public CMultiObjectiveMeasure
{
#define  s_MULTI_PAINTS_MEASURE_OVERHEAD    "MEASURE_Overhead"
public:
	CBinaryMultiPaintsMeasureOverhead(CBinaryMultiObjectiveProblem  *pcParentProblem) : CMultiObjectiveMeasure(pcParentProblem)
	{s_name = s_MULTI_PAINTS_MEASURE_OVERHEAD;}

	virtual double  dGetMax() { return(1); };
	double  dEvaluate(CBinaryCoding *pcFenotype);
	CMultiObjectiveMeasure  *pcClone() { return(new CBinaryMultiPaintsMeasureOverhead(*this)); };

private:
	vector<double>  v_paint_amount_buffer;
};//class  CBinaryMultiPaintsMeasureOverhead


class  CBinaryMultiPaints;

class  CBinaryMultiPaintsRecipe
{
friend class CBinaryMultiPaints;
friend class CBinaryMultiPaintsMeasureMakespan;
friend class CBinaryMultiPaintsMeasureOverhead;
#define  s_MULTI_PAINTS_RECIPE_MACHINES_START    "Machines["
#define  s_MULTI_PAINTS_RECIPE_MACHINES_END      "]"
public:
	CBinaryMultiPaintsRecipe() { v_clear(); };
	CError  eReadFromLine(CString  sLine, int iOffset);
	CString  sGetRecipe();

private:
	int  i_offset;
	vector<int>  v_machines;
	double  d_duration;
	int  i_paint_produced_offset;
	double  d_produced_amount;
	//double  d_product;


	void  v_clear();
};//class  CBinaryMultiPaintsRecipe




#define  s_MULTI_PAINTS_CONFIG_FUNC_TYPE_LINEAR    "lin*"
#define  s_MULTI_PAINTS_CONFIG_FUNC_TYPE_EXP    "exp^"
#define  s_MULTI_PAINTS_CONFIG_FUNC_TYPE_POW    "pow^"


class CBinaryMultiFuncToolFunc
{
public:
	virtual  CString  sToString() = 0;
	virtual double  dGetValue(double  dArg) = 0;
	virtual CError  eConfig(CString  sLine) = 0;
};//class CBinaryMultiFuncToolFunc


class CBinaryMultiFuncToolFuncLinear : public CBinaryMultiFuncToolFunc
{
public:
	CBinaryMultiFuncToolFuncLinear() { d_multiplier = 0; };

	CString  sToString() { CString  s_res; s_res.Format("%s%.2lf", s_MULTI_PAINTS_CONFIG_FUNC_TYPE_LINEAR, d_multiplier);  return(s_res); };
	double  dGetValue(double  dArg) { return(d_multiplier * dArg); };
	void  vConfigure(double  dMultiplier) {d_multiplier = dMultiplier;};
	CError  eConfig(CString  sLine);
private:
	double  d_multiplier;
};//class CBinaryMultiFuncToolFuncLinear : public CBinaryMultiFuncToolFunc


class CBinaryMultiFuncToolFuncExp : public CBinaryMultiFuncToolFunc
{
public:
	CBinaryMultiFuncToolFuncExp() { d_base = 0; };

	CString  sToString() { CString  s_res; s_res.Format("%s%.2lf", s_MULTI_PAINTS_CONFIG_FUNC_TYPE_EXP, d_base);  return(s_res); };
	double  dGetValue(double  dArg);
	void  vConfigure(double  dBase) { d_base = dBase; };
	CError  eConfig(CString  sLine);
private:
	double  d_base;
};//class CBinaryMultiFuncToolFuncExp : public CBinaryMultiFuncToolFunc


class CBinaryMultiFuncToolFuncPow : public CBinaryMultiFuncToolFunc
{
public:
	CBinaryMultiFuncToolFuncPow() { d_pow = 0; };

	CString  sToString() { CString  s_res; s_res.Format("%s%.2lf", s_MULTI_PAINTS_CONFIG_FUNC_TYPE_POW, d_pow);  return(s_res); };
	double  dGetValue(double  dArg);
	void  vConfigure(double  dPow) { d_pow = dPow; };
	CError  eConfig(CString  sLine);
private:
	double  d_pow;
};//class CBinaryMultiFuncToolFuncPow : public CBinaryMultiFuncToolFunc


class CBinaryMultiFuncTool
{
public:
	static uint32_t iERROR_PARENT_CBinaryMultiFuncTool;
	static uint32_t iERROR_CODE_MNULTI_OBJECTIVE_PAINTS_FUNC_UNKNOWN_TYPE;

	CBinaryMultiFuncTool() { pc_func = NULL; };
	~CBinaryMultiFuncTool() { if (pc_func != NULL)  delete  pc_func; };
	CError  eConfig(CString  sLine);
	CString  sToString();

	double  dGetValue(double  dArg);

private:
	CBinaryMultiFuncToolFunc  *pc_func;
};



#define  s_MULTI_PAINTS_CONFIG_AMOUNT_FUNC_TYPE_SUM    "sum*"
#define  s_MULTI_PAINTS_CONFIG_AMOUNT_FUNC_TYPE_LARGEST    "largest*"

class CBinaryMultiFuncToolAmountFunc
{
public:
	virtual  CString  sToString() = 0;
	virtual double  dGetValue(vector<double>  *pvProducedAmounts) = 0;
	virtual CError  eConfig(CString  sLine) = 0;
};//class CBinaryMultiFuncToolAmountFunc



class CBinaryMultiFuncToolAmountFuncSum : public CBinaryMultiFuncToolAmountFunc
{
public:
	CBinaryMultiFuncToolAmountFuncSum() { d_multiplier = 0; };

	CString  sToString() { CString  s_res; s_res.Format("%s%.2lf", s_MULTI_PAINTS_CONFIG_AMOUNT_FUNC_TYPE_SUM, d_multiplier);  return(s_res); };
	double  dGetValue(vector<double>  *pvProducedAmounts);
	void  vConfigure(double  dMultiplier) { d_multiplier = dMultiplier; };
	CError  eConfig(CString  sLine);
private:
	double  d_multiplier;
};//class CBinaryMultiFuncToolAmountFuncSum



class CBinaryMultiFuncToolAmountFuncLargest : public CBinaryMultiFuncToolAmountFunc
{
public:
	CBinaryMultiFuncToolAmountFuncLargest() { d_multiplier = 0; };

	CString  sToString() { CString  s_res; s_res.Format("%s%.2lf", s_MULTI_PAINTS_CONFIG_AMOUNT_FUNC_TYPE_LARGEST, d_multiplier);  return(s_res); };
	double  dGetValue(vector<double>  *pvProducedAmounts);
	void  vConfigure(double  dMultiplier) { d_multiplier = dMultiplier; };
	CError  eConfig(CString  sLine);
private:
	double  d_multiplier;
};//class CBinaryMultiFuncToolAmountFuncLargest



class CBinaryMultiFuncAmountTool
{
public:
	static uint32_t iERROR_PARENT_CBinaryMultiFuncAmountTool;
	static uint32_t iERROR_CODE_MNULTI_OBJECTIVE_PAINTS_AMOUNT_FUNC_UNKNOWN_TYPE;

	CBinaryMultiFuncAmountTool() { pc_func = NULL; };
	~CBinaryMultiFuncAmountTool() { if (pc_func != NULL)  delete  pc_func; };
	CError  eConfig(CString  sLine);
	CString  sToString();

	double  dGetValue(vector<double>  *pvProducedAmounts);

private:
	CBinaryMultiFuncToolAmountFunc  *pc_func;
};



class CBinaryMultiPaints : public CBinaryMultiObjectiveProblem
{
#define  s_MULTI_PAINTS_CONFIG_TYPE    "PaintsConfigType"
#define  s_MULTI_PAINTS_CONFIG_TYPE_HAND_MADE    "PaintsHandMade"
#define  s_MULTI_PAINTS_CONFIG_TYPE_GENERATOR    "PaintsGenerator"

#define  s_MULTI_PAINTS_CONFIG_ALL_PAINTS_CONFIG    "FOR_ALL_PAINTS"


#define  s_MULTI_PAINTS_TASK_MULTIPLE    "MultiplyTask"

#define  s_MULTI_PAINTS_PAINTS    "PaintsNum"
#define  s_MULTI_PAINTS_MACHINES  "MachinesNum"
#define  s_MULTI_PAINTS_RECIPES   "RecipesNum"

friend class CBinaryMultiPaintsMeasureMakespan;
friend class CBinaryMultiPaintsMeasureOverhead;

public:
	static uint32_t iERROR_PARENT_CBinaryMultiPaints;
	static uint32_t iERROR_CODE_MNULTI_OBJECTIVE_PAINTS_UNKNOWN_CONFIG_FILE;

	CBinaryMultiPaints();
	CBinaryMultiPaints(FILE *pfConfig, CError *pcError);
	CBinaryMultiPaints(const CBinaryMultiPaints &pcOther);

	virtual ~CBinaryMultiPaints();

	CError eConfigure(istream *psSettings);
	CError eReport(FILE *pfReport);

	CError  eSave(CString  sDest);
	CError  eSave(FILE  *pfDest);

	CError  eGenerateTestCaseSerie(CString  sFileNameList);
	CError  eGenerateTestCase
		(
			FILE  *pfFileNameList,
			int  iPaints, int iMachines,
			CString sAmountType, double dAmountMulti,
			CString sBaseType, double  dBaseVal,
			CString sEffType, double  dEffVal

		);

	//mutliobjective methods
	void  vEvaluateParetoFront(vector<double> *pvPF, CBinaryCoding *pcFenotype);

protected:
	virtual CError e_init(FILE *pfConfig);
	CError e_load_settings_from_file(FILE *pfConfig);
	CError e_load_settings_from_file_hand_made(FILE *pfConfig);
	CError e_load_settings_from_file_generator(FILE *pfConfig);

	CError e_multiply_task(int  iMultiplier);


	CError  eLoadPFFromManyFileFiles();
	CError  eLoadPFFromFile(FILE  *pfSource);
	void  v_get_files_names_from_summary(CString  sDir, CString  sSummary, vector<CString>  *pvFileNames);


	double d_evaluate(CBinaryCoding *pcFenotype, uint16_t iShift) { return(CBinaryMultiObjectiveProblem::d_evaluate(pcFenotype, iShift)); };
	void  v_prepare_solution(CBinaryCoding *pcFenotype);

private:

	vector<double>  v_paint_demands;
	int  i_machine_number;
	vector<CBinaryMultiPaintsRecipe>  v_recipes;

	//int  i_instance_repeater;

	CBinaryCoding *pc_solution_buf;
	//CBinaryCoding *pc_feno_part;

};//class CBinaryMultiPaints


//--------------------------------------PAINTS PROBLEM--------------------------------------------------------------


#endif//BINARY_EVALUATION_H