#include "nsga2Orig.h"
#include "PopulationOptimizer.h"
#include "UIntCommandParam.h"
#include "FloatCommandParam.h"



using  namespace Nsga2Orig;


uint32_t CNSGA2orig::iERROR_PARENT_CNSGA2origOptimizer = CError::iADD_ERROR_PARENT("iERROR_PARENT_CNSGA2origOptimizer");
uint32_t CNSGA2orig::iERROR_CODE_NSGA2_ORIG_GENOTYPE_LEN_BELOW_0 = CError::iADD_ERROR("iERROR_CODE_NSGA2_ORIG_GENOTYPE_LEN_BELOW_0");


//---------------------------------------------CNSGA2orig-------------------------------------------------------
CNSGA2orig::CNSGA2orig(CProblem<CBinaryCoding, CBinaryCoding>* pcProblem, CLog * pcLog, uint32_t iRandomSeed)
	: CBinaryMultiObjectiveOptimizer(pcProblem, pcLog, iRandomSeed)
	//: CBinaryOptimizer(pcProblem, pcLog, iRandomSeed)
{
	//pc_problem = (CBinaryMultiObjectiveProblem *) pcProblem;

	old_pop_ptr = NULL;
	new_pop_ptr = NULL;
	mate_pop_ptr = NULL;
	globalpop = NULL;

	fpara1 = NULL;
	array = NULL;
	length = NULL;
	gflg = NULL;
	front_pop = NULL;

}//CNSGA2::CNSGA2(CProblem<CBinaryCoding, CBinaryCoding>* pcProblem, CLog * pcLog, uint32_t iRandomSeed)


CNSGA2orig::CNSGA2orig(CNSGA2orig *pcOther) : CBinaryMultiObjectiveOptimizer(pcOther)//CBinaryOptimizer(pcOther)
{
	::MessageBox(NULL, "No implementation: CNSGA2orig::CNSGA2orig(CNSGA2orig *pcOther) : CBinaryMultiObjectiveOptimizer(pcOther)//CBinaryOptimizer(pcOther)", "Implementation missing", MB_OK);
}//CNSGA2::CNSGA2(CNSGA2 *pcOther)


CNSGA2orig::~CNSGA2orig()
{
	if  (old_pop_ptr != NULL) delete old_pop_ptr;
	if (new_pop_ptr != NULL) delete new_pop_ptr;
	if (mate_pop_ptr != NULL) delete mate_pop_ptr;
	if (globalpop != NULL) delete globalpop;

	if (fpara1 != NULL)
	{
		for (int ii = 0; ii < 2 * i_pop_size; ii++)
			delete  fpara1[ii];

		delete  fpara1;
	}//if (fpara1 != NULL)


	if (array != NULL)
	{
		for (int ii = 0; ii < 2 * i_pop_size; ii++)
			delete  array[ii];

		delete  array;
	}//if (array != NULL)


	if (length != NULL)
	{
		for (int ii = 0; ii < 2 * i_pop_size; ii++)
			delete  length[ii];

		delete  length;
	}//if (length != NULL)


	if (gflg != NULL)  delete  gflg;
	if (front_pop != NULL)  delete  front_pop;

	
}//CNSGA2orig::~CNSGA2orig()



CError CNSGA2orig::eConfigure(istream *psSettings)
{
	CError c_err(iERROR_PARENT_CNSGA2origOptimizer);

	c_err = CBinaryOptimizer::eConfigureNSGA2(psSettings);


	if (!c_err)
	{
		CUIntCommandParam p_population_size(POPULATION_OPTIMIZER_ARGUMENT_POPULATION_SIZE);
		i_pop_size = p_population_size.iGetValue(psSettings, &c_err);
	}//if (!c_err)


	if (!c_err)
	{
		CUIntCommandParam p_tournament_size(NSGA2_ORIG_ARGUMENT_CROSS_TYPE);
		optype = p_tournament_size.iGetValue(psSettings, &c_err);
	}//if (!c_err)




	if (!c_err)
	{
		CFloatCommandParam p_prob_cross(NSGA2_ORIG_ARGUMENT_PROB_CROSS);
		pcross = p_prob_cross.fGetValue(psSettings, &c_err);
	}//if (!c_err)

	if (!c_err)
	{
		CFloatCommandParam p_prob_mut(NSGA2_ORIG_ARGUMENT_PROB_MUT);
		pmut_b = p_prob_mut.fGetValue(psSettings, &c_err);
	}//if (!c_err)

	/*if (!c_err)
	{
		CUIntCommandParam p_full_genotype_mut(NSGA2_FULL_GENOTYPE_MUT);
		if (p_full_genotype_mut.iGetValue(psSettings, &c_err) == 1)
			b_whole_genotype_mut = true;
		else
			b_whole_genotype_mut = false;
	}//if (!c_err)*/


	return c_err;
}//CError CNSGA2orig::eConfigure(istream *psSettings)



void CNSGA2orig::vInitialize(time_t tStartTime)
{
	CBinaryOptimizer::vInitialize(tStartTime);
	t_start = tStartTime;

	CError  c_err(iERROR_PARENT_CNSGA2origOptimizer);
	i_templ_length = pc_problem->pcGetEvaluation()->iGetNumberOfElements();

	if (i_templ_length <= 0)
	{
		c_err.vSetError(CNSGA2orig::iERROR_CODE_NSGA2_ORIG_GENOTYPE_LEN_BELOW_0);
		return;
	}//if  (i_templ_length  <=  0)


	nmut = 0;
	ncross = 0;
	ncons = 0;
		
	nfunc = pc_get_multi_problem()->iGetActiveMeasures();

	pc_log->vPrintLine("Initializing...", true);


	c_time_counter.vSetStartNow();
	pc_log->vPrintLine("DONE...", true);


	old_pop_ptr = new CNSGA2origPopulation(i_pop_size);
	new_pop_ptr = new CNSGA2origPopulation(i_pop_size);
	mate_pop_ptr = new CNSGA2origPopulation(i_pop_size);

	pc_log->vPrintLine("pops done...", true);

	globalpop = new globpop(i_pop_size, pc_get_multi_problem()->iGetActiveMeasures(), i_templ_length);

	pc_log->vPrintLine("globpop done...", true);


	init(old_pop_ptr);
	init(new_pop_ptr);
	init(mate_pop_ptr);

	pc_log->vPrintLine("pops init done...", true);


	fpara1 = new float*[2*i_pop_size];
	for  (int  ii = 0; ii < 2 * i_pop_size; ii++)
		fpara1[ii] = new float[2];


	//float  array[2 * maxpop][2] //void CNSGA2orig::gsort(int rnk, int sel)
	array = new float*[2 * i_pop_size];
	for (int ii = 0; ii < 2 * i_pop_size; ii++)
		array[ii] = new float[2];

	//float length[2 * maxpop][2]  //void CNSGA2orig::gshare(int rnk)	
	length = new float*[2 * i_pop_size];
	for (int ii = 0; ii < 2 * i_pop_size; ii++)
		length[ii] = new float[2];

	//int   gflg[2 * maxpop]  //void CNSGA2orig::grank(int gen)
	gflg = new int[2 * i_pop_size];

	front_pop = new int[i_pop_size];


	pc_log->vPrintLine("buffers init done...", true);

	for (int j = 0; j < i_pop_size; j++)
	{
		/*Initializing the Rank array having different individuals
	  at a particular  rank to zero*/
		old_pop_ptr->rankno[j] = 0;
		new_pop_ptr->rankno[j] = 0;
	}

	c_time_counter.vSetStartNow();

	func(old_pop_ptr);
	func(new_pop_ptr);
	func(mate_pop_ptr);


	pc_log->vPrintLine("func done...", true);


}//void CNSGA2orig::vInitialize(time_t tStartTime)



float CNSGA2orig::randomperc()
{
	return((float) RandUtils::dRandNumber(0, 1));
}//float CNSGA2orig::randomperc()


void CNSGA2orig::init(CNSGA2origPopulation *pop_ptr)
{
	int i, j, r;
	float d;
	pop_ptr->ind_ptr = &(pop_ptr->ind[0]);


	CBinaryMultiObjectiveProblem  *pc_multi_obj_problem;
	pc_multi_obj_problem = pc_get_multi_problem();

	/*Loop Over the population size*/
	for (i = 0; i < i_pop_size; i++)
	{

		pop_ptr->ind_ptr->vSetGenotype(i_templ_length);
		pop_ptr->ind_ptr->vSetProblem(pc_multi_obj_problem);

		/*Loop over the chromosome length*/
		for (j = 0; j < i_templ_length; j++)
		{
			/*Generate a Random No. if it is less than 0.5 it
			  generates a 0 in the string otherwise 1*/
			d = randomperc();
			if (d >= 0.5)
			{
				//pop_ptr->ind_ptr->genes[j] = 1;
				pop_ptr->ind_ptr->pcGetGenotype()->piGetBits()[j] = 1;
			}
			else
			{
				//pop_ptr->ind_ptr->genes[j] = 0; 
				pop_ptr->ind_ptr->pcGetGenotype()->piGetBits()[j] = 0;
			}
		}
		pop_ptr->ind_ptr = &(pop_ptr->ind[i + 1]);
	}
	pop_ptr->ind_ptr = &(pop_ptr->ind[0]);
	return;
}//void CNSGA2orig::init(CNSGA2origPopulation *pop_ptr)



void CNSGA2orig::func(CNSGA2origPopulation *pop_ptr)
{

	/*File ptr to the file to store the value of the g for last iteration
	  g is the parameter required for a particular problem
	  Every problem is not required*/

	//float *realx_ptr, /*Pointer to the array of x values*/
	//	*binx_ptr,      /* Pointer to the binary variables */
	//	*fitn_ptr,      /*Pointer to the array of fitness function*/
	//	x[2 * maxvar],     /* problem variables */
	//	f[maxfun],     /*array of fitness values*/
	//	*err_ptr,      /*Pointer to the error */
	//	cstr[maxcons];

	int i, j, k;
	float error, cc;

	pop_ptr->ind_ptr = &(pop_ptr->ind[0]);

	/*Initializing the max rank to zero*/
	pop_ptr->maxrank = 0;
	for (i = 0; i < i_pop_size; i++)
	{
		pop_ptr->ind_ptr = &(pop_ptr->ind[i]);
		//realx_ptr = &(pop_ptr->ind_ptr->xreal[0]);
		//binx_ptr = &(pop_ptr->ind_ptr->xbin[0]);

		/*for (j = 0; j < nvar; j++)
		{ // Real-coded variables 
			x[j] = *realx_ptr++;
		}*/

		/*for (j = 0; j < nchrom; j++)
		{ // Binary-codced variables
			x[nvar + j] = *binx_ptr++;
		}*/

		//fitn_ptr = &(pop_ptr->ind_ptr->fitness[0]);
		//err_ptr = &(pop_ptr->ind_ptr->error);



		/*   DO NOT CHANGE ANYTHING ABOVE   */
		/*----------------------CODE YOUR OBJECTIVE FUNCTIONS HERE------------*/
		/*All functions must be of minimization type, negate maximization
			  functions */
			  /*============Start Coding Your Function From This Point=============*/
			  // First fitness function
		//f[0] = x[0];
		// Second Fitness Function
		//f[1] = x[1];
		pop_ptr->ind_ptr->vRate();

		/*=========End Your Coding Upto This Point===============*/

		/******************************************************************/
		/*              Put The Constraints Here                          */
		/******************************************************************/
		// g(x) >= 0 type (normalize g(x) as in the cstr[1] below)
		/*===========Start Coding Here=============*/

		//cstr[0] = x[0] * x[0] + x[1] * x[1] - 1.0 - 0.1*cos(16.0*atan(x[0] / x[1]));
		//cstr[1] = (-square(x[0] - 0.5) - square(x[1] - 0.5) + 0.5) / 0.5;

		/*===========Constraints Are Coded Upto Here=============*/
		/*   DO NOT CHANGE ANYTHING BELOW  */



		/*for (k = 0; k < nfunc; k++)
		{
			*fitn_ptr++ = f[k];
		}

		for (k = 0; k < ncons; k++)
		{
			pop_ptr->ind_ptr->constr[k] = cstr[k];
		}
		error = 0.0;
		for (k = 0; k < ncons; k++)
		{
			cc = cstr[k];
			if (cc < 0.0)
				error = error - cc;
		}
		*err_ptr = error;*/
	}

	/*---------------------------* RANKING *------------------------------*/

	if (ncons == 0)
		ranking(pop_ptr);
	else
		rankcon(pop_ptr);

	return;
}//void CNSGA2orig::func(CNSGA2origPopulation *pop_ptr)






void CNSGA2orig::ranking(CNSGA2origPopulation *pop_ptr)
{
	int i, j, k,       /*counters*/
		rnk,           /*rank*/
		val,           /*value obtained after comparing two individuals*/
		nondom,        /*no of non dominated members*/
		maxrank1,      /*Max rank of the population*/
		//rankarr[maxpop], /*Array storing the individual number at a rank*/
		q;

	int  *rankarr; /*Array storing the individual number at a rank*/
	rankarr = new int[i_pop_size];



	float *ptr1, *ptr2;

	/*------------------------------* RANKING *------------------------------*/

	/*Initializing the ranks to zero*/
	rnk = 0;

	nondom = 0;
	maxrank1 = 0;

	/*min_fit is initialize to start distributing the dummy fitness =
	  popsize to the rank one individuals and keeping the record such
	  that the minimum fitness of the better rank individual is always
	  greater than max fitness of the relatively worse rank*/


	  /*Difference in the fitness of minimum dummy fitness of better rank
		and max fitness of the next ranked individuals*/

		/*Initializing all the flags to 2*/

	for (j = 0; j < i_pop_size; j++)
	{
		pop_ptr->ind[j].flag = 2;
	}

	q = 0;

	for (k = 0; k < i_pop_size; k++, q = 0)
	{
		for (j = 0; j < i_pop_size; j++)
		{
			if (pop_ptr->ind[j].flag != 1)break;
			/*Break if all the individuals are assigned a rank*/
		}
		if (j == i_pop_size)break;

		rnk = rnk + 1;

		for (j = 0; j < i_pop_size; j++)
		{
			if (pop_ptr->ind[j].flag == 0) pop_ptr->ind[j].flag = 2;
			/*Set the flag of dominated individuals to 2*/
		}

		for (i = 0; i < i_pop_size; i++)
		{
			/*Select an individual which rank to be assigned*/

			pop_ptr->ind_ptr = &(pop_ptr->ind[i]);

			if (pop_ptr->ind_ptr->flag != 1 && pop_ptr->ind_ptr->flag != 0)
			{
				ptr1 = &(pop_ptr->ind_ptr->fitness[0]);

				for (j = 0; j < i_pop_size; j++)
				{

					/*Select the other individual which has not got a rank*/
					if (i != j)
					{
						if (pop_ptr->ind[j].flag != 1)
						{
							pop_ptr->ind_ptr = &(pop_ptr->ind[j]);
							ptr2 = &(pop_ptr->ind_ptr->fitness[0]);

							/*Compare the two individuals for fitness*/
							val = indcmp(ptr1, ptr2);

							/*VAL = 2 for dominated individual which rank to be given*/
							/*VAL = 1 for dominating individual which rank to be given*/

							/*VAL = 3 for non comparable individuals*/

							if (val == 2)
							{
								pop_ptr->ind[i].flag = 0;/* individual 1 is dominated */
								break;
							}

							if (val == 1)
							{
								pop_ptr->ind[j].flag = 0;/* individual 2 is dominated */
							}

							if (val == 3)
							{
								nondom++;/* individual 1 & 2 are non dominated */
								if (pop_ptr->ind[j].flag != 0)
									pop_ptr->ind[j].flag = 3;
							}

						}   /*if loop ends*/
					}       /* i != j loop ends*/
				}           /*loop over j ends*/
				if (j == i_pop_size)
				{

					/*Assign the rank and set the flag*/
					pop_ptr->ind[i].rank = rnk;
					pop_ptr->ind[i].flag = 1;
					rankarr[q] = i;
					q++;
				}
			}       /*Loop over flag check ends*/
		}           /*Loop over i ends */
		pop_ptr->rankno[rnk - 1] = q;
	}
	maxrank1 = rnk;


	/* Find Max Rank of the population    */
	for (i = 0; i < i_pop_size; i++)
	{
		rnk = pop_ptr->ind[i].rank;

		if (rnk > maxrank1)maxrank1 = rnk;

	}

	pop_ptr->maxrank = maxrank1;

	delete  rankarr;
	return;
};//void CNSGA2orig::ranking(CNSGA2origPopulation *pop_ptr)


int CNSGA2orig::indcmp(float *ptr1, float *ptr2)
{
	/*float fit1[maxfun], fit2[maxfun];
	int i, value, m, n;
	for (i = 0; i < nfunc; i++)
	{
		fit1[i] = *ptr1++;
		fit2[i] = *ptr2++;
	}*/

	float *fit1, *fit2;
	int i, value, m, n;
	fit1 = ptr1;
	fit2 = ptr2;

	m = 0;
	n = 0;
	while (m < nfunc && fit1[m] <= fit2[m])
	{
		if ((fit2[m] - fit1[m]) < 1e-7) n++;
		m++;
	}
	if (m == nfunc)
	{
		if (n == nfunc) value = 3;
		else value = 1;             /*value = 1 for dominationg*/
	}
	else
	{
		m = 0;
		n = 0;
		while (m < nfunc && fit1[m] >= fit2[m])
		{
			if ((fit1[m] - fit2[m]) < 1e-7) n++;
			m++;
		}
		if (m == nfunc)
		{
			if (n != nfunc)
				value = 2;                       /*value =  2 for dominated */
			else value = 3;
		}
		else value = 3;                   /*value = 3 for incomparable*/
	}

	return value;
}







void CNSGA2orig::rankcon(CNSGA2origPopulation *pop_ptr)
{
	int i, j, k,       /*counters*/
		rnk,           /*rank*/
		val,           /*value obtained after comparing two individuals*/
		nondom,        /*no of non dominated members*/
		maxrank1,      /*Max rank of the population*/
		//rankarr[maxpop], /*Array storing the individual number at a rank*/
		q;

	int  *rankarr; /*Array storing the individual number at a rank*/
	rankarr = new int[i_pop_size];

	float *ptr1, *ptr2, *err_ptr1, *err_ptr2;

	/*------------------------------* RANKING *------------------------------*/

	/*Initializing the ranks to zero*/
	rnk = 0;

	nondom = 0;
	maxrank1 = 0;
	/*min_fit is initialize to start distributing the dummy fitness =
	  popsize to the rank one individuals and keeping the record such
	  that the minimum fitness of the better rank individual is always
	  greater than max fitness of the relatively worse rank*/

	min_fit = i_pop_size;


	/*Difference in the fitness of minimum dummy fitness of better rank
	  and max fitness of the next ranked individuals*/

	delta_fit = 0.1 * i_pop_size;

	/*Initializing all the flags to 2*/

	for (j = 0; j < i_pop_size; j++)
	{
		pop_ptr->ind[j].flag = 2;
	}

	q = 0;

	for (k = 0; k < i_pop_size; k++, q = 0)
	{
		for (j = 0; j < i_pop_size; j++)
		{
			if (pop_ptr->ind[j].flag != 1)break;
			/*Break if all the individuals are assigned a rank*/
		}
		if (j == i_pop_size)break;

		rnk = rnk + 1;

		for (j = 0; j < i_pop_size; j++)
		{
			if (pop_ptr->ind[j].flag == 0) pop_ptr->ind[j].flag = 2;
			/*Set the flag of dominated individuals to 2*/
		}

		for (i = 0; i < i_pop_size; i++)
		{
			/*Select an individual which rank to be assigned*/

			pop_ptr->ind_ptr = &(pop_ptr->ind[i]);

			if (pop_ptr->ind_ptr->flag != 1 && pop_ptr->ind_ptr->flag != 0)
			{
				ptr1 = &(pop_ptr->ind_ptr->fitness[0]);
				err_ptr1 = &(pop_ptr->ind_ptr->error);

				for (j = 0; j < i_pop_size; j++)
				{

					/*Select the other individual which has not got a rank*/
					if (i != j)
					{
						if (pop_ptr->ind[j].flag != 1)
						{
							pop_ptr->ind_ptr = &(pop_ptr->ind[j]);
							ptr2 = &(pop_ptr->ind_ptr->fitness[0]);
							err_ptr2 = &(pop_ptr->ind_ptr->error);

							if (*err_ptr1 < 1.0e-6 && *err_ptr2 > 1.0e-6)
							{
								/*first ind is feasible second individaul
							  is infeasible*/
								pop_ptr->ind[j].flag = 0;
							}
							else
							{
								if (*err_ptr1 > 1.0e-6 && *err_ptr2 < 1.0e-6)
								{
									/*first individual is infeasible and
									  second is feasible*/
									pop_ptr->ind[i].flag = 0;
									break;
								}
								else
								{
									/*both are feasible or both are infeasible*/
									if (*err_ptr1 > *err_ptr2)
									{
										pop_ptr->ind[i].flag = 0;
										/*first individual is more infeasible*/
										break;
									}
									else
									{
										if (*err_ptr1 < *err_ptr2)
										{
											pop_ptr->ind[j].flag = 0;
											/*second individual is more
											  infeasible*/
										}
										else
										{
											/*Compare the two individuals for
											  fitness*/
											val = indcmp3(ptr1, ptr2);

											/*VAL = 2 for dominated individual
											  which rank to be given*/

											  /*VAL = 1 for dominating individual
												which rank to be given*/

												/*VAL = 3 for non comparable
												  individuals*/

											if (val == 2)
											{
												pop_ptr->ind[i].flag = 0;
												/* individual 1 is dominated */
												break;
											}

											if (val == 1)
											{
												pop_ptr->ind[j].flag = 0;
												/* individual 2 is dominated */
											}

											if (val == 3)
											{
												nondom++;
												/* individual 1 & 2 are
											   non dominated */
												if (pop_ptr->ind[j].flag != 0)
													pop_ptr->ind[j].flag = 3;
											}

										}   /*if loop ends*/
									}       /* i != j loop ends*/
								}
							}
						}
					}
				}        /*loop over j ends*/
				if (j == i_pop_size)
				{
					/*Assign the rank and set the flag*/
					pop_ptr->ind[i].rank = rnk;
					pop_ptr->ind[i].flag = 1;
					rankarr[q] = i;
					q++;
				}
			}       /*Loop over flag check ends*/
		}           /*Loop over i ends */
		pop_ptr->rankno[rnk - 1] = q;
	}
	maxrank1 = rnk;


	/*     Find Max Rank of the population    */
	for (i = 0; i < i_pop_size; i++)
	{
		rnk = pop_ptr->ind[i].rank;

		if (rnk > maxrank1)maxrank1 = rnk;

	}

	pop_ptr->maxrank = maxrank1;

	delete  rankarr;
	return;
}//void CNSGA2orig::rankcon(CNSGA2origPopulation *pop_ptr)





int CNSGA2orig::indcmp3(float *ptr1, float *ptr2)
{
	/*float fit1[maxfun], fit2[maxfun];
	int i, value, m, n;
	for (i = 0; i < nfunc; i++)
	{
		fit1[i] = *ptr1++;
		fit2[i] = *ptr2++;
	}*/

	float *fit1, *fit2;
	int i, value, m, n;
	fit1 = ptr1;
	fit2 = ptr2;

	m = 0;
	n = 0;
	while (m < nfunc && fit1[m] <= fit2[m])
	{
		if (fit1[m] == fit2[m]) n++;
		m++;
	}
	if (m == nfunc)
	{
		if (n == nfunc) value = 3;
		else value = 1;             /*value = 1 for dominationg*/
	}
	else
	{
		m = 0;
		n = 0;
		while (m < nfunc && fit1[m] >= fit2[m])
		{
			if (fit1[m] == fit2[m]) n++;
			m++;
		}
		if (m == nfunc)
		{
			if (n != nfunc)
				value = 2;                       /*value =  2 for dominated */
			else value = 3;
		}
		else value = 3;                   /*value = 3 for incomparable*/
	}

	return value;
}


int CNSGA2orig::indcmp1(float *ptr1, float *ptr2)
{
	/*float fit1[maxfun], fit2[maxfun];
	int i, value, m, n;
	for (i = 0; i < nfunc; i++)
	{
		fit1[i] = *ptr1++;
		fit2[i] = *ptr2++;
	}*/

	float *fit1, *fit2;
	int i, value, m, n;
	fit1 = ptr1;
	fit2 = ptr2;


	m = 0; n = 0;
	while (m < nfunc && fit1[m] <= fit2[m])
	{
		if ((fit2[m] - fit1[m]) < 1e-7) n++;
		m++;
	}
	if (m == nfunc)
	{
		if (n == nfunc) value = 3;
		else value = 1;                    /*value = 1 for dominating*/
	}
	else
	{
		m = 0; n = 0;
		while (m < nfunc && fit1[m] >= fit2[m])
		{
			if ((fit1[m] - fit2[m]) < 1e-7) n++;
			m++;
		}
		if (m == nfunc)
		{
			if (n != nfunc)
				value = 2;                       /*value =  2 for dominated */
			else value = 3;
		}
		else value = 3;                   /*value = 3 for incomparable*/
	}
	return value;
}









CString  CNSGA2orig::sAdditionalSummaryInfo()
{
	CString  s_result;

	s_result.Format("CString  CNSGA2::sAdditionalSummaryInfo() - JUST FILL ME UP!");

	return(s_result);
}//CString  CNSGA2orig::sAdditionalSummaryInfo()



bool CNSGA2orig::bRunIteration_dummy(uint32_t iIterationNumber, time_t tStartTime)
{
	CString  s_buf;
	double  d_time_passed;

	c_time_counter.bGetTimePassed(&d_time_passed);

	CBinaryMultiObjectiveProblem  *pc_multi_problem;
	pc_multi_problem = (CBinaryMultiObjectiveProblem *)pc_problem->pcGetEvaluation();

	

	CBinaryCoding  *pc_genotype;
	pc_genotype = new CBinaryCoding(i_templ_length);


	for (int ii = 0; ii < i_templ_length; ii++)
	{
		pc_genotype->piGetBits()[ii] = RandUtils::iRandNumber(0, 1);
	}//for (int ii = 0; ii < i_templ_length; ii++)


	vector<double>  v_fitness;
	pc_multi_problem->vEvaluateParetoFront(&v_fitness, pc_genotype);

	double  d_res = 1;
	for (int ii = 0; ii < v_fitness.size(); ii++)
		d_res *= v_fitness.at(ii);

	b_update_best_individual(iIterationNumber, tStartTime, pc_genotype->piGetBits(), d_res);


	/*for (int ii = 0; ii < v_population.size(); ii++)
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



	c_time_counter.bGetTimePassed(&d_time_passed);
	//s_buf.Format("iteration: %d  pop: %d; GenSize:%d PFsize:%d HyperVolume: %.8lf InvGenDist: %.8lf DominatedPF: %d [time:%.2lf]", iIterationNumber, v_population.size(), v_non_dominated_pfs.at(0).at(0)->pc_genotype->iGetNumberOfBits(), (int) dPFQualityPointNum(), dPFQualityHyperVolume(), dPFQualityInverseGenerationalDistance(), (int)dPFQualityDominatedOptimalPoints(), d_time_passed);
	s_buf.Format("iteration: %d  pop: %d; GenSize:%d PFsize:%d HyperVolume: %.8lf InvGenDist: %.8lf GenDist:%.8lf MaxSpread:%.8lf DominatedPF: %d [time:%.2lf] [ffe: %.0lf]", iIterationNumber, v_population.size(), v_non_dominated_pfs.at(0).at(0)->pc_genotype->iGetNumberOfBits(), (int)dPFQualityPointNum(), dPFQualityHyperVolume(), dPFQualityInverseGenerationalDistance(), dPFQualityGenerationalDistance(), dPFQualityMaximumSpread(), (int)dPFQualityDominatedOptimalPoints(), d_time_passed, (double)pc_multi_problem->iGetFFE());
	pc_log->vPrintLine(s_buf, true);

	b_update_best_individual(iIterationNumber, tStartTime, v_non_dominated_pfs.at(0).at(0)->pc_genotype->piGetBits(), v_non_dominated_pfs.at(0).at(0)->dEvaluate());*/

	return(true);
};//bool CNSGA2orig::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)



void CNSGA2orig::nselect(CNSGA2origPopulation *old_pop_ptr, CNSGA2origPopulation *pop2_ptr)
{
	int *fit_ptr1, *fit_ptr2;

	float rnd2, *f1_ptr, *f2_ptr;

	int *s1_ptr, *s2_ptr, *select_ptr, *ptry;
	float /**select_ptr_r,*/ *s1_ptr_r, *s2_ptr_r;

	void *j, *j1;

	int c, i, rnd, rnd1, k, n, j2, r, /*s,*/ r1;  //PRW 2019.06.23 - unsed

	old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[0]);

	pop2_ptr->ind_ptr = &(pop2_ptr->ind[0]);

	j = &(old_pop_ptr->ind[i_pop_size - 1]);

	old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[0]);
	j2 = 0;
	r = i_pop_size;
	//s = chrom;  //PRW 2019.06.23 - unsed

	for (n = 0, k = 0; n < i_pop_size; n++, k++)
	{
		pop2_ptr->ind_ptr = &(pop2_ptr->ind[k]);
		select_ptr = &(pop2_ptr->ind_ptr->pcGetGenotype()->piGetBits()[0]);
		//select_ptr_r = &(pop2_ptr->ind_ptr->xreal[0]);

		rnd2 = randomperc();

		rnd2 = i_pop_size * rnd2;

		rnd = floor(rnd2);

		if (rnd == 0)
			rnd = i_pop_size - k;

		if (rnd == i_pop_size)
			rnd = (i_pop_size - 2) / 2;

		/*Select first parent randomly*/
		j = &(old_pop_ptr->ind[rnd - 1]);

		rnd2 = randomperc();

		rnd2 = i_pop_size * rnd2;

		rnd1 = floor(rnd2);

		if (rnd1 == 0)
			rnd1 = i_pop_size - n;

		if (rnd1 == i_pop_size)
			rnd1 = (i_pop_size - 4) / 2;


		/*Select second parent randomly*/
		j1 = &(old_pop_ptr->ind[rnd1 - 1]);

		old_pop_ptr->ind_ptr = (CNSGA2origIndividual*) j;

		s1_ptr = &(old_pop_ptr->ind_ptr->pcGetGenotype()->piGetBits()[0]);
		//s1_ptr_r = &(old_pop_ptr->ind_ptr->xreal[0]);
		fit_ptr1 = &(old_pop_ptr->ind_ptr->rank);
		f1_ptr = &(old_pop_ptr->ind_ptr->cub_len);

		old_pop_ptr->ind_ptr = (CNSGA2origIndividual*) j1;
		s2_ptr = &(old_pop_ptr->ind_ptr->pcGetGenotype()->piGetBits()[0]);
		//s2_ptr_r = &(old_pop_ptr->ind_ptr->xreal[0]);
		fit_ptr2 = &(old_pop_ptr->ind_ptr->rank);
		f2_ptr = &(old_pop_ptr->ind_ptr->cub_len);
		/*--------------------------------------------------------------------------*/

		/*------------------SELECTION PROCEDURE------------------------------------*/

			  /*Comparing the fitnesses*/

		if (*fit_ptr1 > *fit_ptr2)
		{
			//for (i = 0; i < chrom; i++)
			for (i = 0; i < i_templ_length; i++)
				*select_ptr++ = *s2_ptr++;
			//for (i = 0; i < nvar; i++)
				//*select_ptr_r++ = *s2_ptr_r++;
		}
		else
		{
			if (*fit_ptr1 < *fit_ptr2)
			{
				//for (i = 0; i < chrom; i++)
				for (i = 0; i < i_templ_length; i++)
					*select_ptr++ = *s1_ptr++;
				//for (i = 0; i < nvar; i++)
					//*select_ptr_r++ = *s1_ptr_r++;
			}
			else
			{
				if (*f1_ptr < *f2_ptr)
				{
					//for (i = 0; i < chrom; i++)
					for (i = 0; i < i_templ_length; i++)
						*select_ptr++ = *s2_ptr++;
					//for (i = 0; i < nvar; i++)
						//*select_ptr_r++ = *s2_ptr_r++;
				}
				else
				{
					//for (i = 0; i < chrom; i++)
					for (i = 0; i < i_templ_length; i++)
						*select_ptr++ = *s1_ptr++;
					//for (i = 0; i < nvar; i++)
						//*select_ptr_r++ = *s1_ptr_r++;
				}
			}
		}
	}
	return;
}





void CNSGA2orig::mutate(CNSGA2origPopulation *new_pop_ptr)
{
	int i, *ptr, j, r;
	float rand1, *rand_float_ptr;

	rand1 = randomperc();
	new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[0]);

	for (j = 0; j < i_pop_size; j++)
	{
		ptr = &(new_pop_ptr->ind_ptr->pcGetGenotype()->piGetBits()[0]);
		new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[j + 1]);

		/*Select bit */
		for (i = 0; i < i_templ_length; i++)
		{
			rand1 = randomperc();

			/*Check whether to do mutation or not*/
			if (rand1 <= pmut_b)
			{
				if (*ptr == 0)
					*ptr = 1;
				else
					*ptr = 0;
				nmut++;
			}
			ptr++;
		}
	}
	return;
}



void CNSGA2orig::mutate_bitwise(CNSGA2origPopulation *new_pop_ptr)
{
	int i, *ptr, j, r;
	float rand1, *rand_float_ptr;

	rand1 = randomperc();
	new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[0]);


	float  f_effective_pmut_b;
	f_effective_pmut_b = pmut_b;
	f_effective_pmut_b = f_effective_pmut_b / i_templ_length;

	for (j = 0; j < i_pop_size; j++)
	{
		ptr = &(new_pop_ptr->ind_ptr->pcGetGenotype()->piGetBits()[0]);
		new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[j + 1]);

		/*Select bit */
		for (i = 0; i < i_templ_length; i++)
		{
			rand1 = randomperc();

			/*Check whether to do mutation or not*/
			if (rand1 <= f_effective_pmut_b)
			{
				if (*ptr == 0)
					*ptr = 1;
				else
					*ptr = 0;
				nmut++;
			}
			ptr++;
		}
	}
	return;
}





bool CNSGA2orig::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)
{
	//return(bRunIteration_dummy(iIterationNumber, tStartTime));

	CString  s_buf;
	double  d_time_passed;

	c_time_counter.bGetTimePassed(&d_time_passed);

	CBinaryMultiObjectiveProblem  *pc_multi_problem;
	pc_multi_problem = (CBinaryMultiObjectiveProblem *)pc_problem->pcGetEvaluation();

	nfunc = pc_multi_problem->iGetActiveMeasures();

	int  maxrank1, j, l;
	float tot;


	/*s_buf.Format("Generation = %d\n", iIterationNumber);
	pc_log->vPrintLine(s_buf, true);

	s_buf.Format("Population at generation no. -->%d\n", iIterationNumber);
	pc_log->vPrintLine(s_buf, true);
	s_buf.Format("#Generation No. -->%d\n", iIterationNumber);
	pc_log->vPrintLine(s_buf, true);
	s_buf.Format("#Variable_vector  Fitness_vector Constraint_violation Overall_penalty\n");
	pc_log->vPrintLine(s_buf, true);*/


	/*--------SELECT----------------*/
	nselect(old_pop_ptr, mate_pop_ptr);

	//new_pop_ptr = &(newpop);
	//mate_pop_ptr = &(matepop);


	b_update_best_individual(iIterationNumber, tStartTime, new_pop_ptr->ind[0].pcGetGenotype()->piGetBits(), new_pop_ptr->ind[0].dEvaluate()); 
	

	/*CROSSOVER----------------------------*/
	//if (nchrom > 0)
	if (i_templ_length > 0)
	{

		if (optype == 1)
		{
			crossover(new_pop_ptr, mate_pop_ptr);
			/*Binary Cross-over*/
		}

		if (optype == 2)
		{
			unicross(new_pop_ptr, mate_pop_ptr);
			/*Binary Uniform Cross-over*/
		}
	}
	//if (nvar > 0)
		//realcross(new_pop_ptr, mate_pop_ptr);
	///*Real Cross-over*/

	

	/*------MUTATION-------------------*/
	//new_pop_ptr = &(newpop);

	//if (nchrom > 0)
	if (i_templ_length > 0)
		//mutate(new_pop_ptr);
		mutate_bitwise(new_pop_ptr);
	/*Binary Mutation */

	
	//if (nvar > 0)
		//real_mutate(new_pop_ptr);
	///*Real Mutation*/

	//new_pop_ptr = &(newpop);

	///*-------DECODING----------*/
	//if (nchrom > 0)
	//	decode(new_pop_ptr);
	///*Decoding for binary strings*/

	/*----------FUNCTION EVALUATION-----------*/
	//new_pop_ptr = &(newpop);
	func(new_pop_ptr);

	
	/*-------------------SELECTION KEEPING FRONTS ALIVE--------------*/
	//old_pop_ptr = &(oldpop);
	//new_pop_ptr = &(newpop);
	//mate_pop_ptr = &(matepop);

	/*Elitism And Sharing Implemented*/
	keepalive(old_pop_ptr, new_pop_ptr, mate_pop_ptr, iIterationNumber);


	//mate_pop_ptr = &(matepop);
	//if (nchrom > 0)
	//	decode(mate_pop_ptr);

	//mate_pop_ptr = &(matepop);
	///*------------------REPORT PRINTING--------------------------------*/
	//report(i, old_pop_ptr, mate_pop_ptr, rep_ptr, gen_ptr, lastit);

	///*==================================================================*/

	/*----------------Rank Ratio Calculation------------------------*/
	//new_pop_ptr = &(matepop);
	//old_pop_ptr = &(oldpop);

	/*Finding the greater maxrank among the two populations*/

	if (old_pop_ptr->maxrank > new_pop_ptr->maxrank)
		maxrank1 = old_pop_ptr->maxrank;
	else
		maxrank1 = new_pop_ptr->maxrank;

	//fprintf(rep2_ptr, "--------RANK AT GENERATION %d--------------\n", i + 1);
	//fprintf(rep2_ptr, "Rank old ranks   new ranks     rankratio\n");

	for (j = 0; j < maxrank1; j++)
	{
		/*Sum of the no of individuals at any rank in old population
		  and the new populaion*/

		tot = (old_pop_ptr->rankno[j]) + (new_pop_ptr->rankno[j]);

		/*Finding the rank ratio for new population at this rank*/

		new_pop_ptr->rankrat[j] = (new_pop_ptr->rankno[j]) / tot;

		/*Printing this rank ratio to a file called ranks.dat*/

	//	fprintf(rep2_ptr, " %d\t  %d\t\t %d\t %f\n", j + 1, old_pop_ptr->rankno[j], new_pop_ptr->rankno[j], new_pop_ptr->rankrat[j]);

	}

	//fprintf(rep2_ptr, "-----------------Rank Ratio-------------------\n");
	///*==================================================================*/

	///*=======Copying the new population to old population======*/

	//old_pop_ptr = &(oldpop);
	//new_pop_ptr = &(matepop);

	for (j = 0; j < i_pop_size; j++)
	{
		old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[j]);
		new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[j]);
		//if (nchrom > 0)
		{
			/*For Binary GA copying of the chromosome*/

			for (l = 0; l < i_templ_length; l++)
				old_pop_ptr->ind_ptr->pcGetGenotype()->piGetBits()[l] = new_pop_ptr->ind_ptr->pcGetGenotype()->piGetBits()[l];

			//for (l = 0; l < nchrom; l++)
				//old_pop_ptr->ind_ptr->xbin[l] = new_pop_ptr->ind_ptr->xbin[l];
		}
		//if (nvar > 0)
		{
			/*For Real Coded GA copying of the chromosomes*/
			//for (l = 0; l < nvar; l++)
				//old_pop_ptr->ind_ptr->xreal[l] = new_pop_ptr->ind_ptr->xreal[l];
		}

		/*Copying the fitness vector */
		for (l = 0; l < nfunc; l++)
			old_pop_ptr->ind_ptr->fitness[l] = new_pop_ptr->ind_ptr->fitness[l];

		/*Copying the dummy fitness*/
		old_pop_ptr->ind_ptr->cub_len = new_pop_ptr->ind_ptr->cub_len;

		/*Copying the rank of the individuals*/
		old_pop_ptr->ind_ptr->rank = new_pop_ptr->ind_ptr->rank;

		/*Copying the error and constraints of the individual*/

		old_pop_ptr->ind_ptr->error = new_pop_ptr->ind_ptr->error;
		/*for (l = 0; l < ncons; l++)
		{
			old_pop_ptr->ind_ptr->constr[l] = new_pop_ptr->ind_ptr->constr[l];
		}*/

		/*Copying the flag of the individuals*/
		old_pop_ptr->ind_ptr->flag = new_pop_ptr->ind_ptr->flag;
	}   // end of j

	maxrank1 = new_pop_ptr->maxrank;

	/*Copying the array having the record of the individual
 // at different ranks */
	for (l = 0; l < i_pop_size; l++)
	{
		old_pop_ptr->rankno[l] = new_pop_ptr->rankno[l];
	}

	///*Copying the maxrank */
	old_pop_ptr->maxrank = new_pop_ptr->maxrank;

	///*Printing the fitness record for last generation in a file last*/
	//if (i == gener - 1)
	//{  // for the last generation 
	//	old_pop_ptr = &(matepop);
	//	for (f = 0; f < popsize; f++) // for printing
	//	{
	//		old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[f]);

	//		if ((old_pop_ptr->ind_ptr->error <= 0.0) && (old_pop_ptr->ind_ptr->rank == 1))  // for all feasible solutions and non-dominated solutions
	//		{
	//			for (l = 0; l < nfunc; l++)
	//				fprintf(end_ptr, "%f\t", old_pop_ptr->ind_ptr->fitness[l]);
	//			for (l = 0; l < ncons; l++)
	//			{
	//				fprintf(end_ptr, "%f\t", old_pop_ptr->ind_ptr->constr[l]);
	//			}
	//			if (ncons > 0)
	//				fprintf(end_ptr, "%f\t", old_pop_ptr->ind_ptr->error);
	//			fprintf(end_ptr, "\n");

	//			if (nvar > 0)
	//			{
	//				for (l = 0; l < nvar; l++)
	//				{
	//					fprintf(g_var, "%f\t", old_pop_ptr->ind_ptr->xreal[l]);
	//				}
	//				fprintf(g_var, "  ");
	//			}

	//			if (nchrom > 0)
	//			{
	//				for (l = 0; l < nchrom; l++)
	//				{
	//					fprintf(g_var, "%f\t", old_pop_ptr->ind_ptr->xbin[l]);
	//				}
	//			}
	//			fprintf(g_var, "\n");
	//		}  // feasibility check
	//	} // end of f (printing)

	//} // for the last generation


	//b_update_best_individual(iIterationNumber, tStartTime, pc_genotype->piGetBits(), d_res);


	///*for (int ii = 0; ii < v_population.size(); ii++)
	//	v_population.at(ii)->vRate();

	//v_non_dominated_sorting();
	//v_compute_crowding_distance();
	////v_report_pareto_front(0);
	////::Tools::vShow("asd");

	//v_evolution();

	//v_non_dominated_sorting();
	//v_compute_crowding_distance();
	//v_report_pareto_front(1);

	//v_get_half_ofChildren_parent_pop();

	////v_update_global_pareto_front();

	//

	//c_time_counter.bGetTimePassed(&d_time_passed);
	////s_buf.Format("iteration: %d  pop: %d; GenSize:%d PFsize:%d HyperVolume: %.8lf InvGenDist: %.8lf DominatedPF: %d [time:%.2lf]", iIterationNumber, v_population.size(), v_non_dominated_pfs.at(0).at(0)->pc_genotype->iGetNumberOfBits(), (int) dPFQualityPointNum(), dPFQualityHyperVolume(), dPFQualityInverseGenerationalDistance(), (int)dPFQualityDominatedOptimalPoints(), d_time_passed);
	//s_buf.Format("iteration: %d  pop: %d; GenSize:%d PFsize:%d HyperVolume: %.8lf InvGenDist: %.8lf GenDist:%.8lf MaxSpread:%.8lf DominatedPF: %d [time:%.2lf] [ffe: %.0lf]", iIterationNumber, v_population.size(), v_non_dominated_pfs.at(0).at(0)->pc_genotype->iGetNumberOfBits(), (int)dPFQualityPointNum(), dPFQualityHyperVolume(), dPFQualityInverseGenerationalDistance(), dPFQualityGenerationalDistance(), dPFQualityMaximumSpread(), (int)dPFQualityDominatedOptimalPoints(), d_time_passed, (double)pc_multi_problem->iGetFFE());
	//pc_log->vPrintLine(s_buf, true);

	//b_update_best_individual(iIterationNumber, tStartTime, v_non_dominated_pfs.at(0).at(0)->pc_genotype->piGetBits(), v_non_dominated_pfs.at(0).at(0)->dEvaluate());*/

	return(true);
};//bool CNSGA2orig::bRunIteration(uint32_t iIterationNumber, time_t tStartTime)




void CNSGA2orig::unicross(CNSGA2origPopulation *new_pop_ptr, CNSGA2origPopulation *mate_pop_ptr)
{
	CString  s_buf;

	int i, j, r, *gene, y, n, *par1, *par2, *chld1, *chld2;
	float rnd;
	for (i = 0, y = 0, n = 0; i < i_pop_size/2; i++)
	{
		for (j = 0; j < i_templ_length; j++)
		{

			/*Select a bit for doing cross-over*/
			new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[y]);
			chld1 = &(new_pop_ptr->ind_ptr->pcGetGenotype()->piGetBits()[j]);

			new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[y + 1]);
			chld2 = &(new_pop_ptr->ind_ptr->pcGetGenotype()->piGetBits()[j]);

			mate_pop_ptr->ind_ptr = &(mate_pop_ptr->ind[n]);
			par1 = &(mate_pop_ptr->ind_ptr->pcGetGenotype()->piGetBits()[j]);

			mate_pop_ptr->ind_ptr = &(mate_pop_ptr->ind[n + 1]);
			par2 = &(mate_pop_ptr->ind_ptr->pcGetGenotype()->piGetBits()[j]);

			rnd = randomperc();

			/*Checking whether to do cross-over or not*/
			if (rnd <= pcross)
			{
				ncross++;
				*chld1 = *par2;
				*chld2 = *par2;
			}

			else
			{
				*chld1 = *par1;
				*chld2 = *par2;
			}
		}

		y = y + 2;
		n = n + 2;

		//s_buf.Format("i=%d    y=%d      n=%d", i, y, n);
		//Tools::vShow(s_buf);

	}

	for (i = 0; i < i_pop_size; i++)
	{
		new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[i]);
		gene = &(new_pop_ptr->ind_ptr->pcGetGenotype()->piGetBits()[0]);
		for (j = 0; j < i_templ_length; j++)
		{
			gene = &(new_pop_ptr->ind_ptr->pcGetGenotype()->piGetBits()[j]);
		}
	}
	return;
}





void CNSGA2orig::crossover(CNSGA2origPopulation *new_pop_ptr, CNSGA2origPopulation *mate_pop_ptr)
{
	int i, j, k, l, m, n, y, mating_site, *par1, *par2, *chld1, *chld2, c;
	float rnd;
	int r;
	rnd = randomperc();

	new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[0]);

	mate_pop_ptr->ind_ptr = &(mate_pop_ptr->ind[0]);

	for (i = 0, y = 0, n = 0; i < i_pop_size / 2; i++)
	{
		new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[n]);
		chld1 = &(new_pop_ptr->ind_ptr->pcGetGenotype()->piGetBits()[0]);
		n = n + 1;

		new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[n]);
		chld2 = &(new_pop_ptr->ind_ptr->pcGetGenotype()->piGetBits()[0]);
		n = n + 1;

		mate_pop_ptr->ind_ptr = &(mate_pop_ptr->ind[y]);
		par1 = &(mate_pop_ptr->ind_ptr->pcGetGenotype()->piGetBits()[0]);
		y = y + 1;

		mate_pop_ptr->ind_ptr = &(mate_pop_ptr->ind[y]);
		par2 = &(mate_pop_ptr->ind_ptr->pcGetGenotype()->piGetBits()[0]);
		y = y + 1;

		rnd = randomperc();
		if (rnd < pcross)
		{
			ncross++;
			rnd = randomperc();
			c = floor(rnd*(i_templ_length + 10));
			mating_site = c;
			if (mating_site >= i_templ_length)
			{
				mating_site = mating_site / 2;
			}

			for (k = 0; k < i_templ_length; k++)
			{
				if (k > mating_site - 1)
				{
					*chld1++ = *par2++;
					*chld2++ = *par1++;
				}
				else
				{
					*chld1++ = *par1++;
					*chld2++ = *par2++;
				}
			}
		}
		else
		{
			for (k = 0; k < i_templ_length; k++)
			{
				*chld1++ = *par1++;
				*chld2++ = *par2++;
			}
		}
	}
	return;
}



globpop::globpop(int  iPopSize, int  iFuncNum, int  iProblemLength)
{
	i_pop_size = iPopSize;
	i_func_num = iFuncNum;
	i_problem_length = iProblemLength;


	//rankar[2 * maxpop][2 * maxpop], /*record of array of individual numbers at a particular rank */
	rankar = new int*[2 * i_pop_size];
	for (int ii = 0; ii < 2 * i_pop_size; ii++)
		rankar[ii] = new int[2 * i_pop_size];

	//rankno[2 * maxpop];           /*record of no. of individuals at a particular rank*/
	rankno = new int[2 * i_pop_size];


	//int genes[2 * maxpop][maxchrom],
	//	rank[2 * maxpop],            /*rank of different individuals*/
	//	flag[2 * maxpop];            /*Setting the flag */
	genes = new int*[2 * i_pop_size];
	for (int ii = 0; ii < 2 * i_pop_size; ii++)
		genes[ii] = new int[i_problem_length];

	rank = new int[2 * i_pop_size];
	flag = new int[2 * i_pop_size];
	cub_len = new float[2 * i_pop_size];



	//float fitness[2 * maxpop][maxfun], /*Fitness function values for the different
	fitness = new float*[2 * i_pop_size];
	for (int ii = 0; ii < 2 * i_pop_size; ii++)
		fitness[ii] = new float[i_func_num];

}//globpop::globpop(int  iPopSize, int  iFuncNum, int  iProblemLength)


globpop::~globpop()
{
	for (int ii = 0; ii < 2 * i_pop_size; ii++)
	{
		delete  rankar[ii];
		delete  genes[ii];
		delete  fitness[ii];
	}//for (int ii = 0; ii < 2 * i_pop_size; ii++)
	delete  rankar;
	delete  genes;
	delete  fitness;

	delete  rankno;
	delete  rank;
	delete  flag;
}//globpop::~globpop()



void CNSGA2orig::grank(int gen)
{
	int i, j, k, rnk, val, nondom, popsize1, /*gflg[2 * maxpop],*/ q;
	//int  *gflg;
	float *ptr1, *ptr2;
	//FILE *gr;
	//gr = fopen("g_rank_record.out", "a");
	//fprintf(gr, "Genration no. = %d\n", gen);

	//gflg = new int[2 * i_pop_size];

	/*----------------------------* RANKING *---------------------------------*/
	rnk = 0;
	nondom = 0;
	popsize1 = 2 * i_pop_size;

	for (i = 0; i < popsize1; i++)
	{
		gflg[i] = 2;
	}

	for (k = 0; k < popsize1; k++)
	{
		q = 0;
		for (j = 0; j < popsize1; j++)
		{
			if (gflg[j] != 1) break;
		}
		if (j == popsize1) break;
		rnk = rnk + 1;
		for (j = 0; j < popsize1; j++)
		{
			if (gflg[j] == 0) gflg[j] = 2;
		}
		for (i = 0; i < popsize1; i++)
		{
			if (gflg[i] != 1 && gflg[i] != 0)
			{
				ptr1 = &(globalpop->fitness[i][0]);
				for (j = 0; j < popsize1; j++)
				{
					if (i != j)
					{
						if (gflg[j] != 1)
						{
							ptr2 = &(globalpop->fitness[j][0]);
							val = indcmp1(ptr1, ptr2);
							if (val == 2)
							{
								gflg[i] = 0;/* individual 1 is dominated */
								break;
							}
							if (val == 1)
							{
								gflg[j] = 0;/* individual 2 is dominated */
							}
							if (val == 3)
							{
								nondom++;/* individual 1 & 2 are non dominated */
								if (gflg[j] != 0)gflg[j] = 3;
							}
						}
					}
				}
				if (j == popsize1)
				{
					globalpop->rank[i] = rnk;
					gflg[i] = 1;
					globalpop->rankar[rnk - 1][q] = i;
					q++;
				}
			}
		}
		globalpop->rankno[rnk - 1] = q;
	}
	globalpop->maxrank = rnk;

	/*fprintf(gr, "   RANK     No Of Individuals\n");
	for (i = 0; i < rnk; i++)
		fprintf(gr, "\t%d\t%d\n", i + 1, globalpop.rankno[i]);*/

	//fclose(gr);
	//delete  gflg;
	return;

}//void CNSGA2orig::grank(int gen)



void CNSGA2orig::gshare(int rnk)
{
	float /*length[2 * maxpop][2],*/ max;
	int i, j, m1, a;
	float min, Diff;  // Added 18.08.2003

	//float  **length;
	//length = new float*[2*i_pop_size];
	//for (int ii = 0; ii < 2 * i_pop_size; ii++)
//		length[ii] = new float[2];

	m1 = globalpop->rankno[rnk - 1];

	for (j = 0; j < nfunc; j++)
	{
		for (i = 0; i < m1; i++)
		{
			fpara1[i][0] = 0;
			fpara1[i][1] = 0;
		}

		for (i = 0; i < m1; i++)
		{
			a = globalpop->rankar[rnk - 1][i];
			fpara1[i][0] = (float)a;
			fpara1[i][1] = globalpop->fitness[a][j];
		}

		sort(m1); /*Sort the arrays in ascending order of the fitness*/

		max = fpara1[m1 - 1][1];
		min = fpara1[0][1];  // Added 18.08.2003
		Diff = max - min;      // Added 18.08.2003 and 5 subsequent lines
		if (Diff < 0.0)
		{
			pc_log->vPrintLine("Something wrong in keepaliven.h\n", true);
			exit(1);
		}
		for (i = 0; i < m1; i++)
		{
			if (i == 0 || i == (m1 - 1))
			{
				length[i][0] = fpara1[i][0];
				length[i][1] = 100 * max;
			}
			else
			{
				length[i][0] = fpara1[i][0];
				length[i][1] = fabs(fpara1[i + 1][1] - fpara1[i - 1][1]) / Diff; // crowding distances are normalized 18.08.2003
			}
		}
		for (i = 0; i < m1; i++)
		{
			a = length[i][0];
			globalpop->cub_len[a] += length[i][1];
		}
	}


	//for (int ii = 0; ii < 2 * i_pop_size; ii++)
		//delete  length[ii];
	//delete  length;

	return;
}


void CNSGA2orig::sort(int m1)
{
	float temp, temp1;
	int i1, j1, k1;
	for (k1 = 0; k1 < m1 - 1; k1++)
	{
		for (i1 = k1 + 1; i1 < m1; i1++)
		{
			if (fpara1[k1][1] > fpara1[i1][1])
			{
				temp = fpara1[k1][1];
				temp1 = fpara1[k1][0];
				fpara1[k1][1] = fpara1[i1][1];
				fpara1[k1][0] = fpara1[i1][0];
				fpara1[i1][1] = temp;
				fpara1[i1][0] = temp1;
			}
		}
	}
	return;
}



void CNSGA2orig::keepalive(CNSGA2origPopulation *pop1_ptr, CNSGA2origPopulation *pop2_ptr, CNSGA2origPopulation *pop3_ptr, int gen)
{
	int i, j, jj, k, m, a1, l, /*front_pop,*/ rec;

	int sum, st, str, pool, poolf, sel, r1;

	int *gene1_ptr, *gene2_ptr, leftsum, x;

	float rnd, a, *gene3_ptr, x3, *gene4_ptr, *xbin1_ptr, *xbin2_ptr;

	/*Forming the global mating pool*/
		
	nfunc = pc_get_multi_problem()->iGetActiveMeasures();
	
	for (i = 0; i < i_pop_size; i++)
	{
		//if (nchrom > 0)
		{
			/*Binary Coded GA genes are copied*/
			for (k = 0; k < i_templ_length; k++)
			{
				globalpop->genes[i][k] = pop1_ptr->ind[i].pcGetGenotype()->piGetBits()[k];
				globalpop->genes[i + i_pop_size][k] = pop2_ptr->ind[i].pcGetGenotype()->piGetBits()[k];
			}
			/*for (k = 0; k < nchrom; k++)
			{
				globalpop.xbin[i][k] = pop1_ptr->ind[i].xbin[k];
				globalpop.xbin[i + popsize][k] = pop2_ptr->ind[i].xbin[k];
			}*/
		}
		//if (nvar > 0)
		//{
		//	/*For Real Coded GA x values are copied */
		//	for (k = 0; k < nvar; k++)
		//	{
		//		globalpop.xreal[i][k] = pop1_ptr->ind[i].xreal[k];
		//		globalpop.xreal[i + popsize][k] = pop2_ptr->ind[i].xreal[k];
		//	}
		//}

		/*Fitness is copied to the global pool */
		for (l = 0; l < nfunc; l++)
		{
			globalpop->fitness[i][l] = pop1_ptr->ind[i].fitness[l];
			globalpop->fitness[i + i_pop_size][l] = pop2_ptr->ind[i].fitness[l];
		}

		/*Initial;ising the dummyfitness to zero */
		/*globalpop.cub_len[i] = 0;
		globalpop.cub_len[i + popsize] = 0;
		globalpop.error[i] = pop1_ptr->ind[i].error;
		globalpop.error[i + popsize] = pop2_ptr->ind[i].error;
		for (jj = 0; jj < ncons; jj++)
		{
			globalpop.constr[i][jj] = pop1_ptr->ind[i].constr[jj];
			globalpop.constr[i + popsize][jj] = pop2_ptr->ind[i].constr[jj];
		}*/
	}
	
	//global_pop_ptr = &(globalpop);

	/*Finding the global ranks */
	//if (ncons == 0)
		grank(gen);
	//else
		//grankc(gen);

	
	m = globalpop->maxrank;

	/* Sharing the fitness to get the dummy fitness */
	for (i = 0; i < m; i++)
	{
		gshare(i + 1);
	}

	poolf = i_pop_size;
	pool = 0;


	/*Initializing the flags of population to zero */
	for (i = 0; i < 2 * i_pop_size; i++)
	{
		globalpop->flag[i] = 0;
	}

	
	// decide which all solutions belong to the pop3 
	rec = 0;
	st = 0;
	for (i = 0; i < m; i++)
	{
		/*    Elitism Applied Here     */
		st = pool;
		pool += globalpop->rankno[i];

		if (pool <= i_pop_size)
		{
			for (k = 0; k < 2 * i_pop_size; k++)
			{
				if (globalpop->rank[k] == i + 1)
					globalpop->flag[k] = 1;
			}
			pop3_ptr->rankno[i] = globalpop->rankno[i];
		}
		else
		{
			sel = i_pop_size - st;
			Lastrank = i + 1;
			pop3_ptr->rankno[i] = sel;
			gsort(i + 1, sel);
			break;
		}
	}

	
	k = 0;
	for (i = 0, k = 0; i < 2 * i_pop_size && k < i_pop_size; i++)
	{
		//if (nchrom > 0)
		{
			if (globalpop->flag[i] == 1)
			{
				gene1_ptr = &(globalpop->genes[i][0]);
				//xbin1_ptr = &(globalpop->xbin[i][0]);
				pop3_ptr->ind_ptr = &(pop3_ptr->ind[k]);
				//gene2_ptr = &(pop3_ptr->ind_ptr->genes[0]);
				gene2_ptr = &(pop3_ptr->ind_ptr->pcGetGenotype()->piGetBits()[0]);
				//xbin2_ptr = &(pop3_ptr->ind_ptr->xbin[0]);

				//for (j = 0; j < chrom; j++)
				for (j = 0; j < i_templ_length; j++)
				{
					*gene2_ptr++ = *gene1_ptr++;
				}
				//for (j = 0; j < nchrom; j++)
					//*xbin2_ptr++ = *xbin1_ptr++;
			}
		}
		/*if (nvar > 0)
		{
			if (globalpop.flag[i] == 1)
			{
				gene3_ptr = &(globalpop.xreal[i][0]);
				pop3_ptr->ind_ptr = &(pop3_ptr->ind[k]);
				gene4_ptr = &(pop3_ptr->ind_ptr->xreal[0]);

				for (j = 0; j < nvar; j++)
				{
					*gene4_ptr++ = *gene3_ptr++;
				}
			}
		}*/
		if (globalpop->flag[i] == 1)
		{
			for (j = 0; j < nfunc; j++)
				pop3_ptr->ind[k].fitness[j] = globalpop->fitness[i][j];
			pop3_ptr->ind[k].cub_len = globalpop->cub_len[i];

			//if (ncons != 0)
				//pop3_ptr->ind[k].error = globalpop->error[i];
			//for (jj = 0; jj < ncons; jj++)
				//pop3_ptr->ind[k].constr[jj] = globalpop->constr[i][jj];
			pop3_ptr->ind[k].rank = globalpop->rank[i];

			k++;  // increment the pop3 counter
		}
	}

	pop3_ptr->maxrank = Lastrank;

	return;
}



void CNSGA2orig::gsort(int rnk, int sel)
{
	int i, j, a, q;
	float /*array[2 * maxpop][2],*/ temp, temp1;

	q = globalpop->rankno[rnk - 1];

	for (i = 0; i < q; i++)
	{
		array[i][0] = globalpop->rankar[rnk - 1][i];
		a = globalpop->rankar[rnk - 1][i];
		array[i][1] = globalpop->cub_len[a];
	}
	for (i = 0; i < q; i++)
	{
		for (j = i + 1; j < q; j++)
		{
			if (array[i][1] < array[j][1])
			{
				temp = array[i][1];
				temp1 = array[i][0];
				array[i][1] = array[j][1];
				array[i][0] = array[j][0];

				array[j][1] = temp;
				array[j][0] = temp1;
			}
		}
	}

	for (i = 0; i < sel; i++)
	{
		a = array[i][0];
		globalpop->flag[a] = 1;
	}
	return;
}





//---------------------------------------------CNSGA2origPopulation-------------------------------------------------------

CNSGA2origPopulation::CNSGA2origPopulation(int  iPopSize) 
{ 
	i_pop_size = iPopSize;
	rankrat = new float[i_pop_size]; 

	rankno = new int[i_pop_size];
	ind = new CNSGA2origIndividual[i_pop_size];
}//CNSGA2origPopulation::CNSGA2origPopulation(int  iPopSize) 


CNSGA2origPopulation::~CNSGA2origPopulation()
{
	delete rankrat;

	delete  rankno;
	delete  [] ind;
}//CNSGA2origPopulation::CNSGA2origPopulation(int  iPopSize) 



//---------------------------------------------CNSGA2origIndividual-------------------------------------------------------
CNSGA2origIndividual::CNSGA2origIndividual()
{
	pc_genotype = NULL;
	fitness = NULL;
};//CNSGA2origIndividual::CNSGA2origIndividual()


CNSGA2origIndividual::~CNSGA2origIndividual()
{
	if (fitness != NULL)  delete  fitness;
}//CNSGA2origIndividual::~CNSGA2origIndividual()


void  CNSGA2origIndividual::vRate()
{
	b_rated = false;

	CMultiIndividual::vRate();

	if (fitness == NULL)
	{
		fitness = new float[v_fitness.size()];
	}//if (fitness == NULL)

	for (int ii = 0; ii < v_fitness.size(); ii++)
		fitness[ii] = (float) v_fitness.at(ii);

}//void  CNSGA2origIndividual::vRate()


void  CNSGA2origIndividual::vSetProblem(CBinaryMultiObjectiveProblem  *pcMultiObjProblem)
{
	pc_problem = pcMultiObjProblem;
}//void  CNSGA2origIndividual::vSetProblem(CBinaryMultiObjectiveProblem  *pcMultiObjProblem)


double  CNSGA2origIndividual::dEvaluate()
{
	double  d_res;

	vRate();

	d_res = 1;
	for (int ii = 0; ii < v_fitness.size(); ii++)
		d_res *= v_fitness.at(ii);

	return(d_res);
}//double  CNSGA2origIndividual::dEvaluate()


CMultiIndividual  *CNSGA2origIndividual::pcClone()
{
	::MessageBox(NULL, "No implementation: CMultiIndividual  *CNSGA2origPopulation::pcClone()", "Implementation missing", MB_OK);
	return(NULL);
};//CMultiIndividual  *CNSGA2origPopulation::pcClone()