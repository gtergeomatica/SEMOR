#include "rtklib.h"


double weightMatrixDllL1(double S4L1, double cN0L1, prcopt_t opt, int roundDec)
{
	/* Calculation of varianceand weight DLL L1
	: param S4_L1 : S4(L1)
	: param c_n0_L1 : SNR(L1)
	: param config : constant params for Conker's equations
	: return : sigma ^ 2 variance DLL L1(Conker's equation)
	w = 1 / sigma ^ 2 weight DLL L1
	*/

	double rval = NAN;

	if (fabs(cN0L1) > 0.0 && S4L1 < opt.gter_park_s4)
	{
		double S4L1Sq = S4L1 * S4L1;
		double mult = pow(10.0, roundDec);

		rval = opt.gter_park_bn_dll_l1 * opt.gter_park_d * (1.0 + 1.0 / (opt.gter_park_eta_l1 * cN0L1 * (1.0 - 2.0 * S4L1Sq))) / (2.0 * cN0L1 * (1.0 - S4L1Sq));
		rval *= opt.gter_park_w_l1 * opt.gter_park_w_l1; /* meters conversion */
		rval = round(mult / rval) / mult;
	}

	return rval;
}

double weightMatrixDllL2(double S4L1, double S4L2, double cN0L1, double cN0L2, prcopt_t opt, int roundDec)
{
	/* Calculation of varianceand weight DLL L2
	:param S4_L1 : S4(L1)
	: param S4_L2 : S4(L2)
	: param c_n0_L1 : SNR(L1)
	: param c_n0_L2 : SNR(L2)
	: param config : constant params for Conker's equations
	: return : sigma ^ 2 variance DLL L2(Conker's equation)
	w = 1 / sigma ^ 2 weight DLL L2
	*/

	double rval = NAN;

	if (fabs(cN0L1) > 0.0 && S4L1 < opt.gter_park_s4 && fabs(cN0L2) > 0.0 && S4L2 < opt.gter_park_s4)
	{
		double mult = pow(10.0, roundDec);
		rval = opt.gter_park_bn_dll_l2 * (1.0 + 1.0 / (2.0 * opt.gter_park_eta_l2 * cN0L1 * (1.0 - S4L1 * S4L1))) / (2.0 * cN0L2 * (1.0 - S4L2 * S4L2));
		rval *= opt.gter_park_w_l2 * opt.gter_park_w_l2; /* meters conversion */
		rval = round(mult / rval) / mult;
	}

	return rval;
}

double weightMatrixPllL1(double S4L1, double cN0L1, double tL1, double pL1, prcopt_t opt, int roundDec)
{
	/* Calculation of varianceand weight PLL L1
	: param S4_L1 : S4(L1)
	: param c_n0_L1 : SNR(L1)
	: param T_L1 : T(L1)
	: param p_L1 : p(L1)
	: param config : constant params for Conker's equations
	: return : sigma ^ 2 variance PLL L1(Conker's equation)
	w = 1 / sigma ^ 2 weight PLL L1
	*/

	double rval = NAN;

	if (fabs(cN0L1) > 0.0 && S4L1 < opt.gter_park_s4)
	{
		double S4L1Sq = S4L1 * S4L1;
		double mult = pow(10.0, roundDec);
		double sigmaT = (opt.gter_park_bn_pll_l1 * (1.0 + 1.0 / (2.0 * opt.gter_park_eta_l1 * cN0L1 * (1.0 - 2.0 * S4L1Sq))) / (cN0L1 * (1.0 - S4L1Sq)));
		double sigmaS = PI * tL1 / ((opt.gter_park_k_pll_l1 * pow(opt.gter_park_fn_pll_l1, pL1 - 1.0) * sin(fabs(2.0 * opt.gter_park_k_pll_l1 + 1.0 - pL1) * PI / (2.0 * opt.gter_park_k_pll_l1))));
		rval = sigmaT + sigmaS + opt.gter_park_sigma_osc;
		rval = sqrt(rval);
		rval *= 180.0 / PI;
		rval *= opt.gter_park_lambda1 / 360.0;
		rval *= rval;
		rval = round(mult / rval) / mult;
	}

	return rval;
}

double weightMatrixPllL2(double S4L1, double S4L2, double cN0L1, double cN0L2, double tL2, double pL2, prcopt_t opt, int roundDec)
{
	/* Calculation of varianceand weight PLL L2
		: param S4_L1 : S4(L1)
		: param c_n0_L1 : SNR(L1)
		: param S4_L2 : S4(L2)
		: param c_n0_L2 : SNR(L2)
		: param T_L2 : T(L2)
		: param p_L2 : p(L2)
		: param config : constant params for Conker's equations
		: return : sigma ^ 2 variance PLL L2(Conker's equation)
		w = 1 / sigma ^ 2 weight PLL L2
		*/

	double rval = NAN;

	if (fabs(cN0L2) > 0.0 && S4L1 < opt.gter_park_s4 && fabs(cN0L2) > 0.0 && S4L2 < opt.gter_park_s4)
	{
		double mult = pow(10.0, roundDec);
		double aSq = opt.gter_park_a * opt.gter_park_a;
		double sigmaT = (opt.gter_park_bn_pll_l2 * (1.0 + 1.0 / (2.0 * opt.gter_park_eta_l2 * cN0L2 * (1.0 - 2.0 * S4L2 * S4L2))) / (cN0L2 * (1.0 - S4L2 * S4L2)));
		double sigmaS = (aSq + 1.0 / aSq) * (PI * tL2 / (opt.gter_park_k_pll_l2 * (pow(opt.gter_park_fn_pll_l2, pL2 - 1.0)) * sin(fabs(2.0 * opt.gter_park_k_pll_l2 + 1.0 - pL2) * PI / (2.0 * opt.gter_park_k_pll_l2))));
		rval = sigmaT + sigmaS + opt.gter_park_sigma_osc;
		rval = sqrt(rval);
		rval *= 180.0 / PI;
		rval *= opt.gter_park_lambda2 / 360.0;
		rval *= rval;
		rval = round(mult / rval) / mult;
	}

	return rval;
}


double* ComputeWeightsWithConker(obsd_t* obs, double* scintParms, prcopt_t opt, int roundDec)
{
	/* sat_obs observables must be in the following order
			Pseudorange_f1, CarrierPhase_f1, SNR_f1, Pseudorange_f2, CarrierPhase_f2, SNR_f2
	   scint_parms observables must be in the following order
			S4_L1, S4_L2, S4_L5, SIG_L1, SIG_L2, SIG_L5, T_L1, T_L2, T_L5, p_L1, p_L2, p_L5
		rval (weights) return value must be in the following order
			w_DLL_L1, w_DLL_L2, w_DLL_L5, w_PLL_L1, w_PLL_L2, w_PLL_L5
	   */

	double* rval = (double*)malloc(MAX_PARK_PARAMS * sizeof(double));

	double cN0L1 = pow(10.0, 0.1 * obs->SNR[0]/1000.0);
	double cN0L2 = pow(10.0, 0.1 * obs->SNR[1]/1000.0);

	double S4L1 = scintParms[0];
	double S4L2 = scintParms[1];
	double tL1 = scintParms[6];
	double tL2 = scintParms[7];
	double pL1 = scintParms[9];
	double pL2 = scintParms[10];

	rval[0] = weightMatrixPllL1(S4L1, cN0L1, tL1, pL1, opt, roundDec);
	rval[1] = weightMatrixPllL2(S4L1, S4L2, cN0L1, cN0L2, tL2, pL2, opt, roundDec);
	rval[2] = NAN; /* PLL L5 still missing */
	rval[3] = weightMatrixDllL1(S4L1, cN0L1, opt, roundDec);
	rval[4] = weightMatrixDllL2(S4L1, S4L2, cN0L1, cN0L2, opt, roundDec);
	rval[5] = NAN; /* DLL L5 still missing */

	return rval;
}


void WeightObservableWithPark(const rtk_t* rtk, obsd_t* obs, double* recPos, double* satAzel)
{
	int mapIdx;
	int frqIdx;
	double* satPp;
	double* scintParms;
	double* tmp;

	/* warning! test only
	satAzel[0] = 36 * PI / 180.0;
	satAzel[1] = 15 * PI / 180.0;
	obs->SNR[0] = 34250;
	obs->SNR[1] = 42500;
	 end of warning! test only */

	satPp = CalculatePPGeo(recPos[0], recPos[1], rtk->opt.gter_im_h, satAzel);
	mapIdx = SearchOnGeoGrid(rtk->gterIonoGrid, rtk->gterNGridRow, rtk->gterNGridStep, satPp);

	if (mapIdx < 0 || rtk->gterCurrentIonoMap == NULL)
	{
		for (frqIdx = 0; frqIdx < NFREQ; frqIdx++){
			obs->ionoWeight[frqIdx] = 0.0;
			obs->s4[frqIdx] = 999.99; /*LB 999.99 is the null value for scint parameters*/
			obs->sigmaphi[frqIdx] = 999.99;
			obs->T[frqIdx] = 999.99;
			obs->p[frqIdx] = 999.99;
		}
	}
	else
	{
		scintParms = rtk->gterCurrentIonoMap[mapIdx];
		tmp = ComputeWeightsWithConker(obs, scintParms, rtk->opt, 6);
		for (frqIdx = 0; frqIdx < NFREQ; frqIdx++)
		{
			obs->s4[frqIdx] = scintParms[frqIdx * 4];
			obs->sigmaphi[frqIdx] = scintParms[frqIdx * 4 + 1];
			obs->T[frqIdx] = scintParms[frqIdx * 4 + 2];
			obs->p[frqIdx] = scintParms[frqIdx * 4 + 3];

			if (obs->L[frqIdx] > 0.0)
				obs->ionoWeight[frqIdx] = tmp[frqIdx]; /* store carrier phases measures before, then code-pseudoranges*/

			if (obs->P[frqIdx] > 0.0)
				obs->ionoWeight[frqIdx + NFREQ] = tmp[frqIdx + NFREQ];
		}

		free(tmp);
		/*free(scintParms);*/
	}
	free(satPp);
}