/*------------------------------------------------------------------------------
* ionomitigation.c : functions required for ionosphere scintillation mitigation
*                    (GTER addin 2021)
*-----------------------------------------------------------------------------*/

#include "rtklib.h"

/* void decodeKeyValue(char key[], float value, grid_map* mtx) {
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
} */

/* int fillMap(char* token, grid_map* mat, int indLine) {
	int count = 0;
	int flag = 0;
	while (token != NULL) {
		mat->map[(indLine * mat->nrows) + count] = atof(token);
		token = strtok(NULL, " ");
		count++;
	}
	if (count < (mat->ncols-1)) { flag = 1; }
	return flag;
} */

/* void ReadGridMapFile(char* filename, grid_map* map)
* Description:
* This function reads the file with grid parameters (
or SigmaPhi) according to naming convention and formats provided by GTER 
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
/* void ReadGridMapFile(char* filename, grid_map* map)
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

	// Debug: print map
	//int k;
	//int j;
	//for (j = 0; j < indLine; j++) {
	//	for (k = 0; k < map->nrows; k++) {
	//		printf("%lf ", map->map[(j * map->nrows) + k]);
	//	}
	//	printf("\n");
	//}

	fclose(fl);
	return;
} */



double Average(double* coord, int nCoord, int isLat)
{
	double inc = 0.0;

	isLat = isLat > 0;

	for (int i = isLat; i < nCoord; i += 2)
		inc += coord[i];

	return inc / 4.0;
}

void SetGterCorrectionStatus(rtk_t* rtk)
{
	rtk->sol.gter_corr_status = rtk->opt.gter_mdp_enable;
	if (rtk->opt.gter_im_model > 0 && rtk->gterIonoGrid != NULL) rtk->sol.gter_corr_status += 2;
	if (rtk->opt.gter_im_model > 0 && rtk->gterCurrentIonoMap != NULL) rtk->sol.gter_corr_status += 4;
}

double* CalculatePPGeo(double obsLat, double obsLon, double h, double* satHorizCoord)
{
	/*
		Given observer lat(deg), lon(deg) and satellite azimuth(deg), elevation(deg)
		it calculates the Pierce Point of link observer - satellite at height h
		in geodetic coordinates
		: param obs_lat : observer latitude in deg
		: param obs_lon : observer longitude in deg
		: param h : height at which to compute the PP
		: param az : satellite azimuth
		: param el : satellite elevation
		: return : latitude, longitude(deg) of the Pierce Point at level h
		*/

	double Re = 6373.0; /* Earth's radius in km: warning! replace with RE_WGS84/1000.0*/
	double zen = PI / 2.0 - satHorizCoord[1]; /* zenith angle */

	double sbet = Re * sin(zen) / (Re + h);
	double bet = asin(sbet);
	double alf = zen - bet;
	double phi = asin(sin(obsLat) * cos(alf) + cos(obsLat) * sin(alf) * cos(satHorizCoord[0]));
	double slam = sin(alf) * sin(satHorizCoord[0]) / cos(phi);
	double clam = (cos(alf) - sin(obsLat) * sin(phi)) / cos(obsLat) / cos(obsLat);
	double lam = obsLon + atan2(slam, clam);

	/* back to degree */
	double* geod_coord = (double*)malloc(2 * sizeof(double));

	geod_coord[0] = lam * (180.0 / PI); /* pierce point longitude */
	geod_coord[1] = phi * (180.0 / PI); /* pierce point latitude */

	return geod_coord;
}

int SearchOnGeoGrid(double** grid, int nRow, double step, double* ppCoord)
{
	/*
	It searches on a grid defined with lon, lat of each cell,
	where a given PP(lat, lon) is located, i.e.to which cell it belongs, and returns the id of the cell
	This is needed to associate a PP(lat, lon) to a grid cell.
	:param grid : dataframe with coordinates(lat, lon) for each extreme of each cell in a grid
	: param PP_lon : longitude of Pierce Point
	: param PP_lat : latitude of Pierce Point
	: return : id of the cell containing the PP(NaN if not found in any cell)
	*/

	int rval = -1;
	double lonAvg, latAvg;
	double dist = 999999.9, tmpDist;
	double lonDist, latDist;

	for (int i = 0; i < nRow; i++)
	{
		lonAvg = Average(grid[i], MAX_GRID_COL, 0);
		latAvg = Average(grid[i], MAX_GRID_COL, 1);

		lonDist = ppCoord[0] - lonAvg;
		latDist = ppCoord[1] - latAvg;

		if (fabs(lonDist) < step && fabs(latDist) < step)
		{
			tmpDist = sqrt(lonDist * lonDist + latDist * latDist);
			if (tmpDist < dist)
			{
				dist = tmpDist;
				rval = i;
			}
		}
	}

	return rval;
}

double* FillCSValues(char* line, int nCol)
{
	const char* tok;
	int i;
	double* rval;

	i = 0;
	rval = (double*)malloc(nCol * sizeof(double));

	/* Skip first element (index) */
	tok = strtok(line, ";");
	tok = strtok(NULL, ";\n");

	for (; tok && *tok; tok = strtok(NULL, ";\n"))
		sscanf(tok, "%lf", &rval[i++]);

	return rval;
}

double** ReadCsv(char* fPpath, int nCol, int* nRow, double* step, int skipFirstRow)
{
	FILE* fMap;
	char line[1024];
	double** rval;
	int j;

	fMap = fopen(fPpath, "r+");
	if (fMap == NULL)
		return NULL;

	if (skipFirstRow) fgets(line, 1024, fMap); /* skip header */

	*nRow = 0;

	while (fgets(line, 1024, fMap))
		(*nRow)++;

	if (*nRow == 0)
	{
		fclose(fMap);
		return NULL;
	}

	rewind(fMap);
	rval = (double**)malloc(*nRow * sizeof(double*));

	*nRow = 0;
	if (skipFirstRow) fgets(line, 1024, fMap); /* skip header */

	while (fgets(line, 1024, fMap))
	{
		rval[*nRow] = (double*)malloc(nCol * sizeof(double));
		for (j = 0; j < nCol; j++)
			rval[*nRow][j] = NAN;
		rval[(*nRow)++] = FillCSValues(line, nCol);
	}

	fclose(fMap);

	if (step != NULL)
		*step = fabs(rval[1][0] - rval[0][0]);

	return rval;
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
* ... warning! complete ths description!
* O		double* Rim				Covariance matrix output of the function (length: (2 x nf x nSats) x (2 x nf x nSats))
*
* Notes:
* It is currently empty and just a "placeholder". Please modify name, set of input parameters, description, etc. */
void FillRimMatrixWithMdpSnrWeight(int nSats, int nValid, int* iuValid, int nf, double* w, double* Rim)
{
	int i,j,k=0;
	int matrixSize = nValid * 2;
	if (w != NULL)
	{
		for (i = 0; i < nSats; i++)
		{
			for (j = 0; j < nf; j++)
			{
				if (w[i * nf + j] > 0.0)
				{
					Rim[k * (matrixSize + 1)] *=  w[i * nf + j];
					Rim[(k + nValid) * (matrixSize + 1)] *= w[i * nf + j];
					k++;
				}
			}
		}
	}
}

void FillRimMatrixWithIonoWeight(int nSats, int nValid, int* iuValid, int nf, double* w, double* Rim)
{
	int i, j, k = 0;
	int matrixSize = nValid * 2;
	if (w != NULL)
	{
		for (i = 0; i < nSats; i++)
		{
			for (j = 0; j < nf; j++)
			{
				if (w[2 * i * nf + j] >= 0.0)
				{
					Rim[k * (matrixSize + 1)] *=  w[2* i * nf + j];
					Rim[(k + nValid) * (matrixSize + 1)] *= w[2 * i * nf + j + nf];
					k++;
				}
			}
		}
	}
}

/*void FillRimMatrixWithModelTwo(int nSats, int nf, int* satSyss, double* rovSatsAzel, double* refSatsAzel, double* rovHubAzel,
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
}*/

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
void ComputeVarCovarMatrix(int totSat, int* iuValid, int* irValid, int* hubs, rtk_t* rtk, const obsd_t* obs, double* azel, double* wm, double* wi, double* Rim)
{
	int nValid;
	int i, j;
	nValid = CountValidSats(MAXSAT, iuValid);
	for (i = 0; i < 2*rtk->opt.nf* nValid; i++)
		for (j = 0; j < 2 * rtk->opt.nf * nValid; j++)
			Rim[i * 2 * rtk->opt.nf * nValid + j] = i == j ? 1.0 : 0.0;


	double* rovSatsAzel = GetAzElOfValidSatellites(azel, nValid, iuValid);
	double* refSatsAzel = GetAzElOfValidSatellites(azel, nValid, irValid);
	double* rovHubAzel = GetAzElOfValidSatellites(azel, 2, hubs);
	double* refHubAzel = GetAzElOfValidSatellites(azel, 2, &hubs[6]);

	/*double* rovSatsCoord = GetCoordinatesOfValidSatellites(rs, nSats, iuValid);
	double* refSatsCoord = GetCoordinatesOfValidSatellites(rs, nSats, irValid);
	double* rovHubCoord = GetCoordinatesOfValidSatellites(rs, 6, hubs);
	double* refHubCoord = GetCoordinatesOfValidSatellites(rs, 6, &hubs[6]);*/

	int* satSyss = GetSysGivenObservationIndex(nValid, iuValid, obs);

	int* rovSatsSnrs = GetSnrOfValidSatellites(obs, nValid, iuValid);
	int* refSatsSnrs = GetSnrOfValidSatellites(obs, nValid, irValid);
	int* rovHubSnrs = GetSnrOfValidSatellites(obs, 6, hubs);
	int* refHubSnrs = GetSnrOfValidSatellites(obs, 6, &hubs[6]);

	double* rovRecCoords = GetApproximateRoverCoords(rtk);

	if (wm != NULL)
		FillRimMatrixWithMdpSnrWeight(totSat, nValid, iuValid, rtk->opt.nf, wm, Rim);

	if (wi != NULL)
		FillRimMatrixWithIonoWeight(totSat, nValid, iuValid, rtk->opt.nf, wi, Rim);
		
		/*	
		switch (rtk->opt.gter_im_model)
		{
			case 1:
				FillRimMatrixWith(nValid, rtk->opt.nf, satSyss, rovSatsAzel, refSatsAzel, rovHubAzel, refHubAzel, rovRecCoords,
					rovSatsSnrs, refSatsSnrs, rovHubSnrs, refHubSnrs, rtk->s4, rtk->sigmaPhi, Rim);
				break;
					case 2:

							break;
						case 3:
							FillRimMatrixWithModelThree(nValid, rtk->opt.nf, satSyss, rovSatsAzel, refSatsAzel, rovHubAzel, refHubAzel, rovRecCoords,
													  rovSatsSnrs, refSatsSnrs, rovHubSnrs, refHubSnrs, rtk->s4, rtk->sigmaPhi, Rim);
							break;
						case 4:
						default:
							FillRimMatrixWithModelFour(nValid, rtk->opt.nf, satSyss, rovSatsAzel, refSatsAzel, rovHubAzel, refHubAzel, rovRecCoords,
													  rovSatsSnrs, refSatsSnrs, rovHubSnrs, refHubSnrs, rtk->s4, rtk->sigmaPhi, Rim);
							break;
		}*/


	free(rovSatsAzel); free(refSatsAzel); free(rovHubAzel); free(refHubAzel);
	free(rovSatsSnrs); free(refSatsSnrs); free(rovHubSnrs); free(refHubSnrs);
	free(satSyss);
}

void MultiplyElementsOnDiagonalOfMatrix(int nv, double* Rim, double* R, double* Rout)
{
	int i,j;
	for (i = 0; i < nv; i++)
		for (j = 0; j < nv; j++)
			Rout[i * nv + j] = R[i * nv + j];

	for (i = 0; i < nv; i++)
		Rout[i * (nv + 1)] *= Rim[i * (nv + 1)];
}


void StartDebugFile(int crit, double snrTh, double mdpTh)
{
	char path[MAXSTRPATH];
	/*sprintf(path, "C:\\Workarea\\gter\\data\\output\\crit_%d__SNR_%5.1lf_%__MDP_%4.2lf.log", crit, snrTh, mdpTh);*/
#ifdef WIN32
	sprintf(path, "C:\\Workarea\\gter\\data\\output\\crit_%d__SNR_%5.1lf_%__MDP_%4.2lf.log", crit, snrTh, mdpTh);
#else
	sprintf(path, "/mnt/c/Workarea/gter/data/output/dbn_crit_%d__SNR_%5.1lf_%__MDP_%4.2lf.log", crit, snrTh, mdpTh);
#endif
	FILE* f = fopen(path, "w");
	fclose(f);
}

void WriteDebugFile(obsd_t* obs, int n, int crit, double snrTh, double mdpTh)
{
	int i, prn, sysId;
	char sys;
	double ep[6];
	char path[MAXSTRPATH];
	/*sprintf(path, "C:\\Workarea\\gter\\data\\output\\crit_%d__SNR_%5.1lf_%__MDP_%4.2lf.log", crit, snrTh, mdpTh);*/
#ifdef WIN32
	sprintf(path, "C:\\Workarea\\gter\\data\\output\\crit_%d__SNR_%5.1lf_%__MDP_%4.2lf.log", crit, snrTh, mdpTh);
#else
	sprintf(path, "/mnt/c/Workarea/gter/data/output/dbn_crit_%d__SNR_%5.1lf_%__MDP_%4.2lf.log", crit, snrTh, mdpTh);
#endif
	FILE* f = fopen(path, "a");
	time2epoch(obs[0].time, ep);
	fprintf(f, "> %4d %02d %02d %02d %02d %010.7lf\n", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5]);
	for (i = 0; i < n; i++)
	{
		if (obs[i].mdpWeight[0] == 1.0)
		{
			sysId = satsys(obs[i].sat, &prn);
			sys = sysId == SYS_GPS ? 'G' : sysId == SYS_GLO ? 'R' : sysId == SYS_GAL ? 'E' : 'U';
			fprintf(f, "%c%02d\n", sys, prn);
		}

	}
	fclose(f);
}

void WriteDebugMatrix(const char* fileName, int nRows, double* mtrx)
{
	FILE* r; 
	int i, j;
	char row[2560];
	r = fopen(fileName, "w");
	for (i = 0; i < nRows; i++)
	{
		row[0] = 0;

		for (j = 0; j < nRows; j++)
		{
			sprintf(row, "%s %20.15lf", row, mtrx[i * nRows + j]);
		}
		fprintf(r, "%s\n", row);
	}
		
	fclose(r);
}