/*------------------------------------------------------------------------------
* ionomitigation.c : functions required for ionosphere scintillation mitigation
*                    (GTER addin 2021)
*-----------------------------------------------------------------------------*/

#include "rtklib.h"

void decodeKeyValue(char key[], float value, grid_map* mtx) {
	if (strcmp("ncols", key) == 0) {
		mtx->ncols = (int)value;
	}
	else if (strcmp("nrows", key) == 0) {
		mtx->nrows = (int)value;
	}
	else if (strcmp("xllcorner", key) == 0) {
		mtx->xllcorner = (int)value;
	}
	else if (strcmp("yllcorner", key) == 0) {
		mtx->yllcorner = (int)value;
	}
	else if (strcmp("cellsize", key) == 0) {
		mtx->cellsize = value;
	}
	else if (strcmp("NODATA_value", key) == 0) {
		mtx->NODATA_value = value;
	}
	return;
}

int fillMap(char* token, grid_map* mat, int indLine) {
	int count = 0;
	int flag = 0;
	while (token != NULL) {
		mat->map[(indLine * mat->nrows) + count] = atof(token);
		token = strtok(NULL, " ");
		count++;
	}
	if (count < (mat->ncols-1)) { flag = 1; }
	return flag;
}

/* void ReadGridMapFile(char* filename, grid_map* map)
* Description:
* This function reads the file with grid parameters (S4 or SigmaPhi) according to naming convention and formats provided by GTER 
* 
* Return value: 
* void
* 
* Parameters
* I		char* filename		file (including path) to be read
* O		grid_map* map			pointer to the grid_map structure filled with data read from the file
* 
* Notes: 
* if the file does not exists, it cannot be read or has a wrong format, a structure with default (zero) values is returned */
void ReadGridMapFile(char* filename, grid_map* map)
{
	char buf[1024];
	char key[50];
	float value;
	int count = 0;
	char* token;
	int indLine = 0;

	FILE *fl = fopen(filename, "r");
	if (fl == NULL) {
		map->isValid = 0;
		return;
	}

	while (fgets(buf, sizeof(buf), fl) != NULL) {
		if (count < 6) {
			sscanf(buf, "%s %f", key, &value);
			decodeKeyValue(key, value, map);
			count++;
		}
		else {
			if (indLine == 0) { map->map = mat((map->nrows - 1), (map->ncols - 1)); }
			token = strtok(buf, " ");
			if (fillMap(token, map, indLine)) {
				map->isValid = 0;
				break;
			}
			else {
				indLine++;
			}
		}
	}
	map->isValid = 1;

	/* Debug: print map */
	/*int k;
	int j;
	for (j = 0; j < indLine; j++) {
		for (k = 0; k < map->nrows; k++) {
			printf("%lf ", map->map[(j * map->nrows) + k]);
		}
		printf("\n");
	}*/

	fclose(fl);
	return;
}

/* int CountValidSats(int len, int* indexes)
* Description:
* This function scans the 'indexes' vector, which has length 'len', and returns the number of elements with 
* value >=0 at the beginning of the vector. It breaks the search when the first <0 value is found
*
* Return value:
* int:						the number of >=0 elements at the beginning of the vector
*
* Parameters
* I		int len				length of the vector to scan
* I		int* indexes		pointer to the vector to scan
*
* Notes:
* --- */
int CountValidSats(int len, int* indexes)
{
	int retVal = 0;
	int i;
	for (i = 0; i < len; i++)
	{
		if (indexes[i] >= 0)
			retVal++;
		else
			break;
	}

	return retVal;
}

/* double* GetCoordinatesOfValidSatellites(double* rs, int nValid, int* indexes)
* Description:
* This function builds a vector of doubles with the ECEF coordinates of the valid satellites at emission time
*  (whose indexes are taken from the 'indexes' vector), as read from the rs vector (array produced by RTKLIB,
* which contains coordinates and velocity of all tracked satellites)
*
* Return value:
* double*:					length 3 x nValid; the vector with the ECEF coordinates of the valid satellites, 
*							stored according to the following schema:
*                           [X_sV0, Y_sV0, Z_sV0, X_sV1, ..., X_sV(N-1), Y_sV(N-1), Z_sV(N-1)]
*
* Parameters
* I		double* rs			RTKLIB vector containing coordinates and velocity of all tracked satellites
*							at emission time (length 6 x n_tracked_sats)
*							[X_s0, Y_s0, Z_s0, VelX_s0, VelY_s0, VelZ_s0, ..., X_s(N-1), Y_s(N-1), Z_s(N-1), Vel_Xs(N-1), Vel_Ys(N-1), Vel_Zs(N-1)]
* I		int nValid			length of valid elements of the 'indexes' vector
* I		int* indexes		vector containing the indexes of the valid satellites
*
* Notes:
* The function must be called separately for rover and reference satellites */
double* GetCoordinatesOfValidSatellites(double* rs, int nValid, int* indexes)
{
	int i, j;
	double* retCoords = mat(nValid, 3);
	for (i = 0; i < nValid; i++)
		for (j = 0; j < 3; j++)
			retCoords[i*3 + j] = rs[indexes[i] * 6 + j];

	return retCoords;
}

/* double* GetAzElOfValidSatellites(double* azel, int nValid, int* indexes)
* Description:
* This function builds a vector of doubles with the azimuth and elevation of the valid satellites at emission time
*  (whose indexes are taken from the 'indexes' vector), as read from the azel vector (array produced by RTKLIB,
* which contains azimuth and elevation all tracked satellites)
*
* Return value:
* double*:					length 2 x nValid; the vector with the azimuth and elevation of the valid satellites,
*							stored according to the following schema:
*                           [az_sV0, el_sV0, az_sV1, ..., az_sV(N-1), el_sV(N-1)]
*
* Parameters
* I		double* azel		RTKLIB vector containing azimuth and elevation of all tracked satellites
*							at emission time (length 2 x n_tracked_sats)
*							[az_s0, el_s0, ..., az_s(N-1), el_s(N-1)]
* I		int nValid			length of valid elements of the 'indexes' vector
* I		int* indexes		vector containing the indexes of the valid satellites
*
* Notes:
* The function must be called separately for rover and reference satellites */
double* GetAzElOfValidSatellites(double* azel, int nValid, int* indexes)
{
	int i, j;
	double* retAzEl = mat(nValid, 2);
	for (i = 0; i < nValid; i++)
		for (j = 0; j < 2; j++)
			retAzEl[i * 2 + j] = azel[indexes[i] * 2 + j];

	return retAzEl;
}


/* int* GetSnrOfValidSatellites(const obsd_t* obs, int nValid, int* indexes)
* Description:
* This function builds a vector of integers with the SNR values of all tracked frequencies of all 
* valid satellites (whose indexes are taken from the 'indexes' vector), read from the obs structure 
* (structure produced by RTKLIB, which contains all data read from the obs RINEX files or equivalent
* objects)
*
* Return value:
* int*:						vector with the SNR values of the valid satellites, stored according to the
*							following schema:
*							[SNR_l1_sV0, SNR_l2_sV0, SNR_l5_sV0, SNR_l1_sV1, ..., SNR_l1_sV(N-1), SNR_l2_sV(N-1), SNR_l5_sV(N-1)]
*
* Parameters
* I		const obsd_t* obs	RTKLIB sctructure containing observables of all tracked satellites
* I		int nValid			length of valid elements of the 'indexes' vector
* I		int* indexes		vector containing the indexes of the valid satellites
*
* Notes:
* The function must be called separately for rover and reference observables 
* SNR are stored in input and output structures according to RTKLIB convention (1000 x dBHz_value)*/
int* GetSnrOfValidSatellites(const obsd_t* obs, int nValid, int* indexes)
{
	int i, j;
	int* retSnrs = mat(nValid, 3);
	for (i = 0; i < nValid; i++)
		for (j = 0; j < 3; j++)
			retSnrs[i * 3 + j] = obs[indexes[i]].SNR[j];

	return retSnrs;
}

/* double* GetApproximateRoverCoords(rtk_t* rtk)
*  Description:
* This function reads the approximate rover coordinates (output of the single point positioning)
* from the 'rtk' structure, expresssed in meters in ECEF Cartesian reference system
*
* Return value:
* double*:					vector with X, Y and Z coordinates of the rover receiver
*
* Parameters
* I		rtk_t* rtk	RTKLIB sctructure containing the results of RTK processing
*
* Notes:
* --- */
double* GetApproximateRoverCoords(rtk_t* rtk)
{
	int i;
	double retCoords[3];
	for (i = 0; i < 3; i++)
		retCoords[i] = rtk->sol.rr[i];

	return retCoords;
}

int sys2progIdx(int sys)
{
	switch (sys)
	{
		case SYS_GPS:
			return 0;
		case SYS_SBS:
			return 1;
		case SYS_GLO:
			return 2;
		case SYS_GAL:
			return 3;
		case SYS_QZS:
			return 4;
		case SYS_CMP:
			return 5;
		case SYS_IRN:
			return 6;
		default:
			return 0;
	}
}

/* int* GetSysGivenObservationIndex(int nSats, int* obsIdxs, const obsd_t* obs)
*  Description:
* For each valid satellite (whose progressive index in the 'obs' structure is stored in obsIdxs), the function
* retrieved the navigation system index (0:GPS, 1:SBS, 2:GLO, 3:GAL, 4:QZS, 5:CMP, 6:IRN)

* Parameters
* I		int nSats				total number of valid satellites (excluding hubs for each satellite system)
* I		int* indexes		    vector containing the progressive indexes of the valid satellites inside the OBS structure
* I		int* satSyss			integer vector with the id of the navigation system of each not-hub satellite
*								(0: GPS, 1: SBS, 2: GLO, 3: GAL, 4: QZS, 5: CMP, 6: IRN - length nSats)
* I		const obsd_t* obs		RTKLIB sctructure containing observables of all tracked satellites
*
* Notes:
* It is currently empty and just a "placeholder". Please modify name, set of input parameters, description, etc. */
int* GetSysGivenObservationIndex(int nSats, int* obsIdxs, const obsd_t* obs)
{
	int i;
	int* retSyss = mat(nSats, 1);
	for (i = 0; i < nSats; i++)
		retSyss[i] = sys2progIdx(satsys(obs[obsIdxs[i]].sat, NULL));
	return retSyss;
}

/* void FillRimMatrixWithModelOne(int nSats, int nf, double* rovSatsCoord, double* refSatsCoord, double* rovHubCoord,
								double* refHubCoord, double* rovRecCoords, double* rovSatsSnrs, double* refSatsSnrs,
								double* rovHubSnrs, double* refHubSnrs, grid_map s4, grid_map sigmaPhi, double* Rim)
*  Description:
* This function implements the first model for ionosphere scintillation mitigation, based on the input parameters

* Parameters
* I		int nSats				total number of valid satellites (excluding hubs for each satellite system)
* I     int nf					number of frequencies used in processing
* I		int* satSyss			integer vector with the id of the navigation system of each not-hub satellite
*								(0:GPS, 1:SBS, 2:GLO, 3:GAL, 4:QZS, 5:CMP, 6:IRN - length nSats)
* I		double* rovSatsAzel		double vector with the ECEF coordinates of all non-hub satellites at transmission time
*                               of the signal tracked by the rover (length: 3 x nSats)
* I		double* refSatsAzel		double vector with the ECEF coordinates of all non-hub satellites at transmission time
*								of the signal tracked by the reference (length: 3 x nSats)
* I		double* rovHubAzel		double vector with the ECEF coordinates of all hub satellites at transmission time
*								of the signal tracked by the rover (length: 3 x 6)
* I		double* refHubAzel		double vector with the ECEF coordinates of all hub satellites at transmission time
*								of the signal tracked by the reference (length: 3 x 6)
* I		double* rovRecCoords	double vector with the approximate ECEF coordinates of the rover receiver
* I		double* rovSatsSnrs		double vector with the SNR values at all frequencies of all non-hub satellites tracked by the
*								rover (length: 3 x nSats)
* I		double* refSatsSnrs		double vector with the SNR values at all frequencies of all non-hub satellites tracked by the
*								rover (length: 3 x nSats)
* I		double* rovHubSnrs		double vector with the SNR values at all frequencies of all hub satellites tracked by the
*								rover (length: 3 x 6)
* I		double* refHubSnrs		double vector with the SNR values at all frequencies of all hub satellites tracked by the
*								rover (length: 3 x 6)
* I		grid_map s4				structure with the s4 parameters read from the S4 map (.asc) file
* I		grid_map sigmaPhi		structure with the sigmaPhi parameters read from the sigmaPhi map (.asc) file
* O		double* Rim				Covariance matrix output of the function (length: (2 x nf x nSats) x (2 x nf x nSats))
*
* Notes:
* It is currently empty and just a "placeholder". Please modify name, set of input parameters, description, etc. */
void FillRimMatrixWithModelOne(int nSats, int nf, int* satSyss, double* rovSatsAzel, double* refSatsAzel, double* rovHubAzel,
								double* refHubAzel, double* rovRecCoords, double* rovSatsSnrs, double* refSatsSnrs,
								double* rovHubSnrs, double* refHubSnrs, grid_map* s4, grid_map sigmaPhi, double* Rim)
{
	int nv = 2 * nf * nSats; /* number of double-differenced observables, 2 for each non-hub satellite 
							    (one pseudorange and one carrier phase for each frequency of each satellite)*/

	if(s4[0].isValid && sigmaPhi.isValid) /* modify this check to consider the case where all s4 maps are valid */
	{ }
	/* NB structure of Rim to be discussed, to met the same format of the current 'R' matrix produced by RTKLIB */
}

void FillRimMatrixWithModelTwo(int nSats, int nf, int* satSyss, double* rovSatsAzel, double* refSatsAzel, double* rovHubAzel,
								double* refHubAzel, double* rovRecCoords, double* rovSatsSnrs, double* refSatsSnrs,
								double* rovHubSnrs, double* refHubSnrs, grid_map* s4, grid_map sigmaPhi, double* Rim)
{
}

void FillRimMatrixWithModelThree(int nSats, int nf, int* satSyss, double* rovSatsAzel, double* refSatsAzel, double* rovHubAzel,
								double* refHubAzel, double* rovRecCoords, double* rovSatsSnrs, double* refSatsSnrs,
								double* rovHubSnrs, double* refHubSnrs, grid_map* s4, grid_map sigmaPhi, double* Rim)
{
}

void FillRimMatrixWithModelFour(int nSats, int nf, int* satSyss, double* rovSatsAzel, double* refSatsAzel, double* rovHubAzel,
								double* refHubAzel, double* rovRecCoords, double* rovSatsSnrs, double* refSatsSnrs,
								double* rovHubSnrs, double* refHubSnrs, grid_map* s4, grid_map sigmaPhi, double* Rim)
{
}

/* void ComputeVarCovarMatrix(int* iuValid, int* irValid, int* hubs, rtk_t* rtk, const obsd_t* obs, double* rs, double* Rim)

* Description:
* "Main" function for computing the Rim variance-covariance matrix

* Parameters
* I		int* iuValid			vector with the indexes of the valid satellites (non-hub satellites used to compute
*								the double differences) tracked by the rover receiver (length: MAXSAT)
* I		int* irValid			vector with the indexes of the valid satellites (non-hub satellites used to compute
*								the double differences) tracked by the rover receiver (length: MAXSAT)
* I		int* hubs			    vector with the indexes of the hub satellites for each satellite system (length 12: 
*								hubs[0-5] contains the hubs sats indexes tracked by the rover for the six
*								allowed satellite systems (0:GPS/SBS,1:GLO,2:GAL,3:BDS,4:QZS,5:IRN ), hubs[6-11] the 
*								hubs sats indexes tracked by the reference
* I		rtk_t* rtk				RTKLIB sctructure containing the results of RTK processing
* I		const obsd_t* obs		RTKLIB sctructure containing observables of all tracked satellites
* I		double* azel			RTKLIB vector containing azimuth and elevation of all tracked satellites
*								at emission time (length 2 x n_tracked_sats)
*								[az_s0, el_s0, ..., az_s(N-1), el_s(N-1)]
* O		double* Rim				covariance matrix output of the function (length: (2 x rtk->opt.nf x nSats) x (2 x rtk->opt.nf x nSats))
*
*Notes:
*---*/
void ComputeVarCovarMatrix(int* iuValid, int* irValid, int* hubs, rtk_t* rtk, const obsd_t* obs, double* azel, double* Rim)
{
	int nSats;
	int i, j;
	nSats = CountValidSats(MAXSAT, iuValid);
	for (i = 0; i < 2*rtk->opt.nf*nSats; i++)
		for (j = 0; j < 2 * rtk->opt.nf * nSats; j++)
			Rim[i * 2 * rtk->opt.nf * nSats + j] = i == j ? 1.0 : 0.0;

	double* rovSatsAzel = GetAzElOfValidSatellites(azel, nSats, iuValid);
	double* refSatsAzel = GetAzElOfValidSatellites(azel, nSats, irValid);
	double* rovHubAzel = GetAzElOfValidSatellites(azel, 2, hubs);
	double* refHubAzel = GetAzElOfValidSatellites(azel, 2, &hubs[6]);

	/*double* rovSatsCoord = GetCoordinatesOfValidSatellites(rs, nSats, iuValid);
	double* refSatsCoord = GetCoordinatesOfValidSatellites(rs, nSats, irValid);
	double* rovHubCoord = GetCoordinatesOfValidSatellites(rs, 6, hubs);
	double* refHubCoord = GetCoordinatesOfValidSatellites(rs, 6, &hubs[6]);*/

	int* satSyss = GetSysGivenObservationIndex(nSats, iuValid, obs);

	int* rovSatsSnrs = GetSnrOfValidSatellites(obs, nSats, iuValid);
	int* refSatsSnrs = GetSnrOfValidSatellites(obs, nSats, irValid);
	int* rovHubSnrs = GetSnrOfValidSatellites(obs, 6, hubs);
	int* refHubSnrs = GetSnrOfValidSatellites(obs, 6, &hubs[6]);

	double* rovRecCoords = GetApproximateRoverCoords(rtk);

	switch (rtk->opt.im_model)
	{
		case 1:
			FillRimMatrixWithModelOne(nSats, rtk->opt.nf, satSyss, rovSatsAzel, refSatsAzel, rovHubAzel, refHubAzel, rovRecCoords,
									  rovSatsSnrs, refSatsSnrs, rovHubSnrs, refHubSnrs, rtk->s4, rtk->sigmaPhi, Rim);
			break;
		case 2:
			FillRimMatrixWithModelTwo(nSats, rtk->opt.nf, satSyss, rovSatsAzel, refSatsAzel, rovHubAzel, refHubAzel, rovRecCoords,
									  rovSatsSnrs, refSatsSnrs, rovHubSnrs, refHubSnrs, rtk->s4, rtk->sigmaPhi, Rim);
			break;
		case 3:
			FillRimMatrixWithModelThree(nSats, rtk->opt.nf, satSyss, rovSatsAzel, refSatsAzel, rovHubAzel, refHubAzel, rovRecCoords,
				                      rovSatsSnrs, refSatsSnrs, rovHubSnrs, refHubSnrs, rtk->s4, rtk->sigmaPhi, Rim);
			break;
		case 4:
		default:
			FillRimMatrixWithModelFour(nSats, rtk->opt.nf, satSyss, rovSatsAzel, refSatsAzel, rovHubAzel, refHubAzel, rovRecCoords,
				                      rovSatsSnrs, refSatsSnrs, rovHubSnrs, refHubSnrs, rtk->s4, rtk->sigmaPhi, Rim);
			break;
	}

	free(rovSatsAzel); free(refSatsAzel); free(rovHubAzel); free(refHubAzel);
	free(rovSatsSnrs); free(refSatsSnrs); free(rovHubSnrs); free(refHubSnrs);
	free(satSyss);
}