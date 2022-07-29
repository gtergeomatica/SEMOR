#include "rtklib.h"

double ComputeDynamicRisk(double* _break, double* risk, double a, int nRisks)
{
	/* Risk function(Koulouri's method)
		: param _break : list of sorted values giving the break points of the variable range(including max)
		: param risk : list of sorted values to assign to each interval between break points
		: param a : float value to evaluate
		: return r risk associated to a
		*/

	int i;

	for (i = 0; i < nRisks; i++)
	{
		/*if (_break[i] != _break[i]) /* _break[i] is NAN, occhio ai settings di compilazione perch� potrebbe non funzionare */
		/*							/* https://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c) */
		/*	continue; */

		if (isnan(_break[i]))
			continue;

		if (a < _break[i])
			return risk[i];
	}

	return risk[i];
}

double* ComputeWeights(double* scintParms, const prcopt_t opt)
{
	/*
		Calculates Koulouri's weights on the satellites based on associated scint params values
		:param sat_info : dictionary with association id satellites - scintillation params
		: param params : list of scintillation params
		: param _break : dictionary of lists of sorted values giving the break points for each param(including max)
		: param risk : dictionary of lists of sorted values to assign to each interval between break points for each param
		: return : wk dictionary with association id satellites - list of weights as sorted in params list
		*/

	int i;
	double* rval = (double*)malloc(6 * sizeof(double));
	double* wght = (double*)malloc(6 * sizeof(double));
	double scpa[6];
	scpa[0] = scintParms[0]; // S4 L1 codice
	scpa[1] = scintParms[4]; // S4 L2 codice
	scpa[2] = scintParms[8]; // S4 L5 codice
	scpa[3] = scintParms[1]; // SP L1 fase
	scpa[4] = scintParms[5]; // SP L2 fase
	scpa[5] = scintParms[9]; // SP L5 fase

	for (i = 0; i < 6; i++)
	{
		rval[i] = ComputeDynamicRisk(opt.gter_koulouri_breaks[i], opt.gter_koulouri_risks[i], scpa[i], opt.gter_n_valid_koulouri_risks[i]);
		wght[i] = (1.0 - rval[i]) * (1.0 - rval[i]); /* parametrizzare esponente */
	}

	free(rval);
	return wght;
}

void WeightObservableWithKoulouri(const rtk_t* rtk, obsd_t* obs, double* recPos, double* satAzel)
{
	int mapIdx;
	int frqIdx;
	double* satPp;
	double* scintParms;
	double* tmp;
	satPp = CalculatePPGeo(recPos[0], recPos[1], rtk->opt.gter_im_h, satAzel);
	mapIdx = SearchOnGeoGrid(rtk->gterIonoGrid, rtk->gterNGridRow, rtk->gterNGridStep, satPp);

	if (mapIdx < 0 || rtk->gterCurrentIonoMap == NULL)
	{
		for (frqIdx = 0; frqIdx < NFREQ; frqIdx++) {

			obs->ionoWeight[frqIdx] = obs->ionoWeight[frqIdx + NFREQ] = 0.0;
			obs->s4[frqIdx] = 999.99; /*LB 999.99 is the null value for scint parameters*/
			obs->sigmaphi[frqIdx] = 999.99;
			obs->T[frqIdx] = 999.99;
			obs->p[frqIdx] = 999.99;
		}
			
	}
	else
	{
		scintParms = rtk->gterCurrentIonoMap[mapIdx];
		tmp = ComputeWeights(scintParms, rtk->opt);
		for (frqIdx = 0; frqIdx < NFREQ; frqIdx++)
		{
			obs->s4[frqIdx] = scintParms[frqIdx * 4];
			obs->sigmaphi[frqIdx] = scintParms[frqIdx * 4+1];
			obs->T[frqIdx] = scintParms[frqIdx * 4+2];
			obs->p[frqIdx] = scintParms[frqIdx * 4 + 3];
			
			if (obs->L[frqIdx] > 0.0)
				obs->ionoWeight[frqIdx] = fmin(tmp[0], tmp[1]);//tmp[frqIdx + NFREQ]; /* store carrier phases measures before, then code-pseudoranges*/

			if (obs->P[frqIdx] > 0.0)
				obs->ionoWeight[frqIdx + NFREQ] = fmin(tmp[0], tmp[1]); //tmp[frqIdx];
		}

		free(tmp);
		//free(scintParms);
	}

	free(satPp);
}