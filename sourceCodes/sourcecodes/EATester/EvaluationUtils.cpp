#include "EvaluationUtils.h"

#include "BinaryCoding.h"
#include "BinaryEvaluation.h"
#include "BinaryEvaluationMultiObjective.h"
#include "EnumCommandParam.h"
//#include "PermutationCoding.h"
//#include "PermutationEvaluation.h"
//#include "RealCoding.h"
//#include "RealEvaluation.h"
#include "StringUtils.h"

#include <atlstr.h>
#include <unordered_map>
#include <utility>

template <class TFenotype>
CEvaluation<TFenotype> * EvaluationUtils::pcGetEvaluation(istream *psSettings, CError *pcError)
{
	CEvaluation<TFenotype> *pc_evaluation = nullptr;

	size_t i_fenotype_type_hash_code = typeid(TFenotype).hash_code();

	unordered_map<CString, EEvaluationType> m_evaluation_types;

	m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_CONCATENATION, EVALUATION_CONCATENATION));

	if (i_fenotype_type_hash_code == typeid(CBinaryCoding).hash_code())
	{
		m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_BINARY_CYCLIC_FENOTYPE, EVALUATION_BINARY_CYCLIC_FENOTYPE));
		m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_BINARY_DECEPTIVE, EVALUATION_BINARY_DECEPTIVE));
		m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_BINARY_DECEPTIVE_CONCATENATION, EVALUATION_BINARY_DECEPTIVE_CONCATENATION));
		m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_BINARY_KNAPSACK, EVALUATION_BINARY_KNAPSACK));
		m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_BINARY_MAX_SAT, EVALUATION_BINARY_MAX_SAT));
		m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_BINARY_DECEPTIVE_STEP_TRAP, EVALUATION_BINARY_DECEPTIVE_STEP_TRAP));
		m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_BINARY_NEAREST_NEIGHBOR_NK, EVALUATION_BINARY_NEAREST_NEIGHBOR_NK));
		m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_BINARY_ISING_SPIN_GLASS, EVALUATION_BINARY_ISING_SPIN_GLASS));
		m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_BINARY_MAX_3_SAT, EVALUATION_BINARY_MAX_3_SAT));
		m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_BINARY_DISCRETIZED_RASTRIGIN, EVALUATION_BINARY_DISCRETIZED_RASTRIGIN));
		m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_BINARY_DISCRETIZED_ROSENBROCK, EVALUATION_BINARY_DISCRETIZED_ROSENBROCK));
		m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_BINARY_HIFF, EVALUATION_BINARY_HIFF));
		m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_BINARY_CYCLIC_TRAP, EVALUATION_BINARY_CYCLIC_TRAP));
		m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_BINARY_LEADING_1, EVALUATION_BINARY_LEADING_1));
		m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_BINARY_JUMP_FUNC, EVALUATION_BINARY_JUMP_FUNC));

		m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_BINARY_MULTI_REVERSING, EVALUATION_BINARY_MULTI_REVERSING));
		m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_BINARY_MULTI_MAX_CUT_MOGOMEA, EVALUATION_BINARY_MULTI_MAX_CUT_MOGOMEA));
		m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_BINARY_MULTI_KNAPSACK_MOGOMEA, EVALUATION_BINARY_MULTI_KNAPSACK_MOGOMEA));
		m_evaluation_types.insert(pair<const CString, EEvaluationType>(EVALUATION_ARGUMENT_TYPE_BINARY_MULTI_PAINTS, EVALUATION_BINARY_MULTI_PAINTS));
	}//if (i_fenotype_type_hash_code == typeid(CBinaryCoding).hash_code())
	

	CEnumCommandParam<EEvaluationType> p_type(EVALUATION_ARGUMENT_TYPE, &m_evaluation_types);
	EEvaluationType e_type = p_type.eGetValue(psSettings, pcError);

	if (!*pcError)
	{
		switch (e_type)
		{
			case EVALUATION_CONCATENATION:
			{
				pc_evaluation = new CConcatenationEvaluation<TFenotype>();
				break;
			}//case EVALUATION_CONCATENATION
			case EVALUATION_BINARY_CYCLIC_FENOTYPE:
			{
				pc_evaluation = (CEvaluation<TFenotype>*)new CBinaryCyclicFenotypeEvaluation();
				break;
			}//case EVALUATION_BINARY_CYCLIC_FENOTYPE
			case EVALUATION_BINARY_DECEPTIVE:
			{
				pc_evaluation = (CEvaluation<TFenotype>*)new CBinaryDeceptiveEvaluation();
				break;
			}//case EVALUATION_BINARY_DECEPTIVE
			case EVALUATION_BINARY_DECEPTIVE_CONCATENATION:
			{
				pc_evaluation = (CEvaluation<TFenotype>*)new CBinaryDeceptiveConcatenationEvaluation();
				break;
			}//case EVALUATION_BINARY_DECEPTIVE_CONCATENATION
			case EVALUATION_BINARY_KNAPSACK:
			{
				pc_evaluation = (CEvaluation<TFenotype>*)new CBinaryKnapsackEvaluation();
				break;
			}//case EVALUATION_BINARY_KNAPSACK
			case EVALUATION_BINARY_MAX_SAT:
			{
				pc_evaluation = (CEvaluation<TFenotype>*)new CBinaryMaxSatEvaluation();
				break;
			}//case EVALUATION_BINARY_MAX_SAT
			case EVALUATION_BINARY_DECEPTIVE_STEP_TRAP:
			{
				pc_evaluation = (CEvaluation<TFenotype>*)new CBinaryDeceptiveStepTrapEvaluation();
				break;
			}//case EVALUATION_BINARY_DECEPTIVE_STEP_TRAP
			case EVALUATION_BINARY_NEAREST_NEIGHBOR_NK:
			{
				pc_evaluation = (CEvaluation<TFenotype>*)new CBinaryNearestNeighborNKEvaluation();
				break;
			}//case EVALUATION_BINARY_NEAREST_NEIGHBOR_NK
			case EVALUATION_BINARY_ISING_SPIN_GLASS:
			{
				pc_evaluation = (CEvaluation<TFenotype>*)new CBinaryIsingSpinGlassEvaluation();
				break;
			}//case EVALUATION_BINARY_ISING_SPIN_GLASS
			case EVALUATION_BINARY_MAX_3_SAT:
			{
				pc_evaluation = (CEvaluation<TFenotype>*)new CBinaryMax3SatEvaluation();
				break;
			}//case EVALUATION_BINARY_MAX_3_SAT
			case EVALUATION_BINARY_DISCRETIZED_RASTRIGIN:
			{
				pc_evaluation = (CEvaluation<TFenotype>*)new CBinaryDiscretizedRastriginEvaluation();
				break;
			}//case EVALUATION_BINARY_DISCRETIZED_RASTRIGIN
			case EVALUATION_BINARY_DISCRETIZED_ROSENBROCK:
			{
				pc_evaluation = (CEvaluation<TFenotype>*)new CBinaryDiscretizedRosenbrockEvaluation();
				break;
			}//case EVALUATION_BINARY_DISCRETIZED_ROSENBROCK
			case EVALUATION_BINARY_HIFF:
			{
				pc_evaluation = (CEvaluation<TFenotype>*)new CBinaryHIFFEvaluation();
				break;
			}//case EVALUATION_BINARY_HIFF
			case EVALUATION_BINARY_CYCLIC_TRAP:
			{
				pc_evaluation = (CEvaluation<TFenotype>*)new CBinaryCyclicTrapEvaluation();
				break;
			}//case EVALUATION_BINARY_CYCLIC_TRAP
			case EVALUATION_BINARY_LEADING_1:
			{
				pc_evaluation = (CEvaluation<TFenotype>*)new CBinaryLeading1();
				break;
			}//case EVALUATION_BINARY_MULTI_PAINTS
			case EVALUATION_BINARY_JUMP_FUNC:
			{
				pc_evaluation = (CEvaluation<TFenotype>*)new CBinaryJumpFunc();
				break;
			}//case EVALUATION_BINARY_JUMP_FUNC

			case EVALUATION_BINARY_MULTI_REVERSING:
			{
				pc_evaluation = (CEvaluation<TFenotype>*)new CBinaryMultiReversing();
				break;
			}//case EVALUATION_BINARY_MULTI_REVERSING
			case EVALUATION_BINARY_MULTI_MAX_CUT_MOGOMEA:
			{
				pc_evaluation = (CEvaluation<TFenotype>*)new CBinaryMultiMaxcutMoGomea();
				break;
			}//case EVALUATION_BINARY_MULTI_MAX_CUT_MOGOMEA:

			case EVALUATION_BINARY_MULTI_KNAPSACK_MOGOMEA:
			{
				pc_evaluation = (CEvaluation<TFenotype>*)new CBinaryMultiKnapsackMoGomea();
				break;
			}//case EVALUATION_BINARY_MULTI_KNAPSACK_MOGOMEA:

			case EVALUATION_BINARY_MULTI_PAINTS:
			{
				pc_evaluation = (CEvaluation<TFenotype>*)new CBinaryMultiPaints();
				break;
			}//case EVALUATION_BINARY_MULTI_PAINTS
			default:
			{
				pcError->vSetError(CError::iERROR_CODE_OPERATOR_NOT_FOUND, "evaluation");
				break;
			}//default
		}//switch (pcParams->eGetType())

		if (!*pcError)
		{
			*pcError = pc_evaluation->eConfigure(psSettings);
		}//if (!*pcError)
	}//if (!*pcError)

	return pc_evaluation;
}//CEvaluation<TFenotype> * EvaluationUtils::pcGetEvaluation(istream *psSettings, CError *pcError)




template CEvaluation<CBinaryCoding> * EvaluationUtils::pcGetEvaluation(istream*, CError*);
