/*
* Loosely.h
* Loosely Coupled Integration of GPS & IMU
*  Created on: Sept 10, 2017
*      Author: Aaron Boda
*/

#include "pch.h"
#include "FileIO.h"
#include "LooseKF.h"
#include "ReaderIMU.h"
#include "ReaderGNSS.h"
#include "MathUtils.h"
#include "GeoUtils.h"
#include "InitializeIMU.h"
#include "semor.h" // ../../../src/semor.h

#ifndef LOOSELY_H_
#define LOOSELY_H_

class Loosely
{
public:
	// CONSTRUCTOR
	Loosely();
	Loosely(gnss_sol_t fst_pos);
	Loosely(std::string filePathGNSS, std::string filePathIMU, std::string filePathOUT);

	//##CHANGED##
	void init_imu(gnss_sol_t fst_pos);
	void get_imu_sol(gnss_sol_t* int_sol); //In client.c fare una copia del dato gnss, quella copia è la struttura in cui sarà immagazzinato int_sol
	void close_out_file();

	//##CHANGED##
	struct GNSSpos {
		double epoch;
		std::vector<double> posXYZ;
		std::vector<double> velXYZ;
		std::vector<double> accXYZ;
	};

	struct IMUpos {
		double epoch;
		std::vector<double> posXYZ;
		std::vector<double> velXYZ;
		std::vector<double> attXYZ;
	};

	struct INTpos {
		double epoch;
		std::vector<double> posXYZ;
		std::vector<double> velXYZ;
		std::vector<double> attXYZ;
		std::vector<double> df;
		std::vector<double> dw;
	};

	// Attributes
	IMUpos IMUsol;
	INTpos INTsol;
	GNSSpos GNSSsol;

	// For Timing Information
	double _dT;
	double _epochTime;
	double _dTgps;
	double _dTimu;
	double _epochIMU;
	double _epochGNSS;

	// For position
	std::vector<double> _LLH_o; // Reference Position in Geodetic
	std::vector<double> _ECEF_o; // Reference Position Origin in ECEF
	std::vector<double> _ECEF_i; // Rover Position in ECEF
	std::vector<double> _ECEF_imu; // Rover Position in ECEF
	double _Heading_gps;
	double _Heading_imu;

	// Functions
	void read_imu();
	void epochOutput(std::ofstream& fout);
	void epochOutputIMU(std::ofstream& fout);
	void epochOutputGPS(std::ofstream& fout);
	void LooseCoupling(ReaderGNSS &GNSS, IMUmechECEF& Mech);
	void SolutionGNSS(ReaderGNSS GNSS);
	void SolutionIMU(ReaderIMU IMU, IMUmechECEF& Mech);

private:
	// Functions
	bool isValid(std::string observation_filepath);

};

#endif /* LOOSELY_H_ */
