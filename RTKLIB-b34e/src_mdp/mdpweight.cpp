/*------------------------------------------------------------------------------
* ionomitigation.c : functions required for multipath-based observables 
*                    weighting (GTER addin 2021)
*-----------------------------------------------------------------------------*/

#include "rtklib.h"
#include <cfloat>

void FillMdpStructWithCurrentEpochData(mdpStruct* mdp, const obsd_t* obs, int nObs)
{
	int obsIdx, freqIdx, satIdx;
	int sysIdx;
	double deltaCodePhase;
	double lambda;
	double lambdaGPS[NFREQ];
	double lambdaGAL[NFREQ];
	double lambdaGLO[NFREQ];

	ShiftMdpStructEpochContent(mdp);
	lambdaGPS[0] = CLIGHT / FREQ1;
	lambdaGPS[1] = CLIGHT / FREQ2;
	lambdaGPS[2] = CLIGHT / FREQ5;

	lambdaGAL[0] = CLIGHT / FREQ1;
	lambdaGAL[1] = CLIGHT / FREQ7;
	lambdaGAL[2] = CLIGHT / FREQ5;

	lambdaGLO[0] = CLIGHT / FREQ1_GLO; /* warning: approximation! */
	lambdaGLO[1] = CLIGHT / FREQ2_GLO; /* warning: approximation! */
	lambdaGLO[2] = CLIGHT / FREQ3_GLO;

	for (obsIdx = 0; obsIdx < nObs; obsIdx++)
	{
		if (obsIdx == 0)
			mdp->times[0] = obs[obsIdx].time;

		for (freqIdx = 0; freqIdx < NFREQ; freqIdx++)
		{
			sysIdx = satsys(obs[obsIdx].sat, NULL);
			lambda = (sysIdx == SYS_GPS) ? lambdaGPS[freqIdx] : (sysIdx == SYS_GLO) ? lambdaGLO[freqIdx] : (sysIdx == SYS_GAL) ? lambdaGAL[freqIdx] : 0.0;

			if (obs[obsIdx].L[freqIdx] == 0.0 || obs[obsIdx].P[freqIdx] == 0.0)
				deltaCodePhase = NAN;
			else
				deltaCodePhase = obs[obsIdx].P[freqIdx] - obs[obsIdx].L[freqIdx] * lambda;

			mdp->deltaCp[obs[obsIdx].rcv - 1][freqIdx][0][obs[obsIdx].sat] = deltaCodePhase;
		}

	}
}

mdpStruct* InitializeMdpStructure(int mdpEpochs, double thSnr, double thMdp)
{
	int statIdx, epIdx, satIdx, freqIdx;
	mdpStruct* mdp = new mdpStruct();
	mdp->nEpochs = min(mdpEpochs, 10);
	mdp->thMdp = thMdp;
	mdp->thSnr = thSnr;

	for (epIdx = 0; epIdx < mdp->nEpochs; epIdx++) /* cycle on epochs */
	{
		mdp->times[epIdx].time = 0;
		mdp->times[epIdx].sec = 0.0;

		for (statIdx = 0; statIdx < 2; statIdx++) /* cycle on rover (0) and reference (1)*/
		{
			for (freqIdx = 0; freqIdx < NFREQ; freqIdx++) /* cycle on frequencies */
			{
				for (satIdx = 0; satIdx < MAXSAT; satIdx++) /* cycle on satellites*/
				{
					mdp->deltaCp[statIdx][freqIdx][epIdx][satIdx] = (double)NAN;
				}
			}
		}
	}

	return mdp;
}

void ShiftMdpStructEpochContent(mdpStruct* mdp)
{
	int statIdx, epIdx, freqIdx, satIdx;
	for (epIdx = mdp->nEpochs - 1; epIdx > 0; epIdx--) /* cycle on epochs */
	{
		mdp->times[epIdx] = mdp->times[epIdx - 1];

		for (statIdx = 0; statIdx < 2; statIdx++) /* cycle on rover (0) and reference (1)*/
		{
			for (freqIdx =0; freqIdx<NFREQ; freqIdx++) /* cycle on frequencies */
			{
				for (satIdx = 0; satIdx < MAXSAT; satIdx++)
				{
					mdp->deltaCp[statIdx][freqIdx][epIdx][satIdx] = mdp->deltaCp[statIdx][freqIdx][epIdx - 1][satIdx];
				}
			}
		}
	}
}

void StartDebugFile(double snrTh, double mdpTh)
{
	char path[MAXSTRPATH];
	sprintf(path, "C:\\Workarea\\gter\\data\\output\\SNR_%5.1lf_%__MDP_%3.1lf.log", snrTh, mdpTh);
	FILE* f = fopen(path, "w");
	fclose(f);
}

void WriteDebugFile(obsd_t* obs, int n, double snrTh, double mdpTh)
{
	int i, prn, sysId;
	char sys;
	char path[MAXSTRPATH];
	sprintf(path, "C:\\Workarea\\gter\\data\\output\\SNR_%5.1lf_%__MDP_%3.1lf.log", snrTh, mdpTh);
	FILE* f = fopen(path, "a");
	fprintf(f, ">\n");
	for (i = 0; i < n; i++)
	{
		if (obs[i].weight[0] < 100.0)
		{
			sysId = satsys(obs[i].sat, &prn);
			sys = sysId == SYS_GPS ? 'G' : sysId == SYS_GLO ? 'R' : sysId == SYS_GAL ? 'E' : 'U';
			fprintf(f, "%c%02d\n", sys, prn);
		}
			
	}
	fclose(f);
}

int Weight(obsd_t* obs, int satIdx, int frqIdx, int criterion, double coeff, double mdp)
{
	int cond;

	switch (criterion)
	{
		case 1:
			cond = obs[satIdx].snrFlag[frqIdx] == 0;
			break;
		case 2:
			cond = obs[satIdx].mdpFlag[frqIdx] == 0;
			break;
		case 3:
			cond = obs[satIdx].snrFlag[frqIdx] == 0 && obs[satIdx].mdpFlag[frqIdx] == 0;
			break;
		case 4:
			cond = obs[satIdx].snrFlag[frqIdx] == 0 || obs[satIdx].mdpFlag[frqIdx] == 0;
			break;
		default:
			return -1;
	}

	if (cond)
		obs[satIdx].weight[frqIdx] = coeff * fabs(mdp) * pow(10.0, obs[satIdx].SNR[frqIdx] / 10.0);
	else
		obs[satIdx].weight[frqIdx] = 1.0;

	return 0;
}

void WeightObservablesWithMdp(mdpStruct* mdp, obsd_t* obs, int n, int criterion, double coeff)
{
	int obsIdx, freqIdx, epIdx, k;
	double incrMdp;
	double lastMdp;
	int validMdp;
	for (obsIdx = 0; obsIdx < n; obsIdx++) /* cycle on observables */
	{
		for (freqIdx = 0; freqIdx < NFREQ; freqIdx++) /* cycle on frequencies */
		{
			incrMdp = 0.0;
			lastMdp = DBL_MAX;
			validMdp = 1;
			for (epIdx = 0; epIdx < mdp->nEpochs - 1; epIdx++) /* cycle on epochs */
			{
				if (isnan(mdp->deltaCp[obs[obsIdx].rcv - 1][freqIdx][epIdx][obs[obsIdx].sat]) || isnan(mdp->deltaCp[obs[obsIdx].rcv - 1][freqIdx][epIdx + 1][obs[obsIdx].sat]))
				{
					validMdp = 0;
					break;
				}
				else
				{
					incrMdp += (mdp->deltaCp[obs[obsIdx].rcv - 1][freqIdx][epIdx][obs[obsIdx].sat] - mdp->deltaCp[obs[obsIdx].rcv - 1][freqIdx][epIdx + 1][obs[obsIdx].sat]);
					if (epIdx == 0)
						lastMdp = incrMdp;
				}
			}

			if (validMdp)
			{
				double mdpthreshold = mdp->thMdp;
				if (mdp->nEpochs > 2)     /* MDP dynamic --> mdpthreshold represents not a threshold per se but the increment to the average mdp value to set the threshold*/
					mdpthreshold += mdp->nEpochs / (mdp->nEpochs - 1);

				obs[obsIdx].snrFlag[freqIdx] = obs[obsIdx].SNR[freqIdx]/1000.0 > mdp->thSnr;
				obs[obsIdx].mdpFlag[freqIdx] = fabs(lastMdp) < mdp->thMdp;

				if (Weight(obs, obsIdx, freqIdx, criterion, coeff, lastMdp) != 0)
					validMdp = 0;
			}

			if (!validMdp)
				obs[obsIdx].weight[freqIdx] = WEIGHT_INC_OBS; /* check if correct! */
		}
	}
}
