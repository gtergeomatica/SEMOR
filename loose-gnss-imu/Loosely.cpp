/*
* Loosely.cpp
* Loosely Coupled Integration of GPS & IMU
*  Created on: Sept 10, 2018
*      Author: Aaron Boda
*/

#include "pch.h"
#include "Loosely.h"
#include "semor.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/time.h>

using namespace std;
using namespace Eigen;

//raw imu data log file
ofstream out;

// PI
const double PI = 3.1415926535898;

//imu timestamp of the last read sample
double epochAfterRead;

//structure used to make the imu mechanization and so the part that calculates a position given a start position and a sample of imu data
IMUmechECEF MechECEF;

//Used in debug mode, raw imu data input file
ifstream fimu;

//structure used to initialize imu biases and errors
InitializeIMU iniIMU;

//timestamp that tells when the initialization of the imu ends
double IMU_INI_TIME_END;

//structure that contains last read gnss data
ReaderGNSS OBSgnss;

//structure that contains last read imu data
ReaderIMU OBSimu;


pthread_mutex_t lock;
pthread_t imu_thread_id;

//buffer containing raw imu data
//the main thread reads data from this
//the loop_imu_thread() writes data in this
char imu_data[IMUBUF_CAPACITY][IMU_LENGTH];

//cursors for reading and writing in imu_data
int imu_count_write = 0;
int imu_count_read = 0;

Loosely::Loosely(){


}

//close raw imu data log file
void Loosely::close_out_file(){
	out.close();
}

//Used for debugging the collection of imu data
void print_time(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	int week = ((tv.tv_sec+LEAP_SECONDS-GPS_EPOCH))/(7*24*3600);
	double sec = (double)(((tv.tv_sec+LEAP_SECONDS-GPS_EPOCH))%(7*24*3600))+(tv.tv_usec / 1000000.0);
	printf("%lf\n", sec);
}

//read one sample at a time
void Loosely::read_imu(){
	string line;
	if(debug){
		getline(fimu, line);
		if (line.find("GPS") != std::string::npos) {
			getline(fimu, line);
		}
		OBSimu.clearObs();
		OBSimu.obsEpoch(line);
	}
	else{
		char buf[IMU_LENGTH];
		do{
			pthread_mutex_lock(&lock);
			strcpy(buf, imu_data[imu_count_read]);
			pthread_mutex_unlock(&lock);
			if(strlen(buf) != 0){
				line = string(buf);
				OBSimu.clearObs();
				OBSimu.obsEpoch(line);
			}
		}while(OBSimu._IMUdata.imuTime <= _epochIMU);

		imu_count_read++;
		imu_count_read = imu_count_read % IMUBUF_CAPACITY;
	}

	if(logs){
		out << line << endl;
		out.flush();
	}

}


// A routine to facilitate IMU mechanization in ECEF
void Loosely::SolutionIMU(ReaderIMU IMU, IMUmechECEF& MechECEF) {
	// Compute time interval
	_dTimu = IMU._IMUdata.imuTime - _epochIMU; //Difference between current epoch and previous
	// IMU Mechanization
	MechECEF.MechanizerECEF(_dTimu, IMU._IMUdata.Acc, IMU._IMUdata.Gyr, _LLH_o); //Update ECEF position adding to it accelerometer and gyroscope data (previous data is accumulated)
	// Update solution
	_epochIMU = IMU._IMUdata.imuTime; //Update epocIMU with current epoch
	IMUsol.posXYZ = MechECEF._pos;		//Update IMU solution
	IMUsol.velXYZ = MechECEF._vel;		//Update IMU solution
	IMUsol.attXYZ = MechECEF._att;		//Update IMU solution
	_Heading_imu = normalise(IMUsol.attXYZ.at(2), 0, 2 * PI);
}

VectorXd double2eigVector(double a, double b, double c){
	VectorXd v = VectorXd::Zero(3);
	v(0) = a;
	v(1) = b;
	v(2) = c;
	return v;
}

//dedicated thread that reads raw imu data and place it in imu_data in a thread-safe manner
void* loop_imu_thread(void* arg){
	char buf[IMU_LENGTH];
	string line;

	//Initialize lock
	if (pthread_mutex_init(&lock, NULL) != 0)
    {
        printf("\n mutex init failed\n");
		perror("SEMOR: pthread_mutex_init()");
        //TODO: close_semor(1);
		return NULL;
    }

	while(1){
		get_imu_data(buf);

		//modifica matrice "imu_data" in modo thread-safe
		pthread_mutex_lock(&lock);
		strcpy(imu_data[imu_count_write++], buf);
		pthread_mutex_unlock(&lock);
		imu_count_write = imu_count_write % IMUBUF_CAPACITY;
	}
}

double radianToDegree(double r) {
  return r * (180 / PI);
}

gnss_sol_t ecef2geo(gnss_sol_t gnss){

    //1
    // Output vector - Lat, Long, Height
	// Variables
	double x, y, z;
	x = gnss.a; y = gnss.b; z = gnss.c;
	// Semi Major Axis and Eccentricity
	const double a = 6378137; const double e = 0.08181979;
	// Compute Longitude
	double lambda = atan2(y, x);
	// Physical radius of the point
	double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	// Radius in the x-y plane
	double p = sqrt(pow(x, 2) + pow(y, 2));
	// GEOcentric latitude (Initial Approx)
	double phi_o = atan2(p, z); double phi_i = phi_o;
	// Radius of curvature in the prime vertical
	double Rn;
	// Height
	double h;
	// Loop
	for (unsigned i = 0; i < 3; i++) {
		// Recalculate Radius of curvature in the prime vertical
		Rn = a / sqrt(1 - (e * e * sin(phi_i) * sin(phi_i)));
		// Recalculate Height
		h = (p / cos(phi_i)) - (Rn);
		// Recalculate Latitude
		phi_i = atan((z / p) * (pow((1 - ((pow(e, 2))*(Rn / (Rn + h)))), (-1))));
	}
	// Recalculate Height
	h = (p / cos(phi_i)) - Rn;
	// Populate output vector
	/* gnss.a = radianToDegree(phi_i);
	gnss.b = radianToDegree(lambda); */

	//The library takes radian values, not degrees so radianToDegree is not what we want to use here
	gnss.a = phi_i;
	gnss.b = lambda;
	gnss.c = h;
	return gnss;
}

//input: previous epoch best solution (this position will be propagated using imu data until next epoch), output: gnss+imu solution of the next(this) epoch
void Loosely::get_imu_sol(gnss_sol_t* int_sol){
	//Check if imu is initializing
	if(imu_ready == 0){
		int n = 0;
		while(OBSimu._IMUdata.imuTime < (*int_sol).time.sec){ //read whole imu date of a particular gnss epoch
			if(imu_ready == 0 && iniIMU.stepInitializeIMU(OBSimu, IMU_INI_TIME_END, _LLH_o) == 1){ //it calculates _ACCbias and _GYRbias, _RPY
				// Initialize IMU Mechanization
				MechECEF.InitializeMechECEF(_ECEF_imu, _LLH_o, GNSSsol.velXYZ, iniIMU._RPY, iniIMU._ACCbias, iniIMU._GYRbias); //it initializes MechECEF with _ACCbias, _GYRbias and _RPY

				imu_ready = 1; //Tells client.c that it can read imu data
				printf("SEMOR: End initialization\n");
			}
			_epochIMU = OBSimu._IMUdata.imuTime;
			read_imu();
			if(imu_ready == 1){
				imu_ready = 2;
				(*int_sol).time.week = OBSimu._IMUdata.week;
			}
			epochAfterRead = OBSimu._IMUdata.imuTime;
			n++;
		}
		return;
	}
	//Here IMU is ready:

	//initialize MechECEF fields with the input position
	MechECEF._pos.at(0) = (*int_sol).a;
	MechECEF._pos.at(1) = (*int_sol).b;
	MechECEF._pos.at(2) = (*int_sol).c;

	MechECEF._vel.at(0) = (*int_sol).va;
	MechECEF._vel.at(1) = (*int_sol).vb;
	MechECEF._vel.at(2) = (*int_sol).vc;

	double group_time = 0.1;
	double next_time = epochAfterRead + group_time;

	double avg_ax = 0, avg_ay = 0, avg_az = 0;
	double avg_gx = 0, avg_gy = 0, avg_gz = 0;
	int n = 0, nlocal = 0;
	//here we already have the first sample of IMU data that we need (in the last iteration of the next while we read a sample that we use in the next epoch)
	do {
		gnss_sol_t llh = ecef2geo(*int_sol);
		_LLH_o = eigVector2std(double2eigVector(llh.a, llh.b, llh.c));
		SolutionIMU(OBSimu, MechECEF);	//this is the gnss+first imu data //Get position from current MechECEF state and IMU data just read - This has effects on: MechECEF e IMUsol

		(*int_sol).a = IMUsol.posXYZ.at(0);
		(*int_sol).b = IMUsol.posXYZ.at(1);
		(*int_sol).c = IMUsol.posXYZ.at(2);
		(*int_sol).va = IMUsol.velXYZ.at(0);
		(*int_sol).vb = IMUsol.velXYZ.at(1);
		(*int_sol).vc = IMUsol.velXYZ.at(2);

		_epochIMU = OBSimu._IMUdata.imuTime;
		read_imu();
		epochAfterRead = OBSimu._IMUdata.imuTime;
		n++;

	} while (epochAfterRead <= (*int_sol).time.sec);

	//mark the imu solution as available (before (*int_sol).time.week was equal to 0)
	(*int_sol).time.week = OBSimu._IMUdata.week; 
}

//this is called only once
void Loosely::init_imu(gnss_sol_t fst_pos){
	FileIO FIO;

	time_t rawtime;
    struct tm info;
    time( &rawtime );
    info = *localtime( &rawtime );
    char str_time[13];
	char buf[IMU_LENGTH];
	string line;

    sprintf(str_time, "%d_%02d_%02d_%02d_%02d", (info.tm_year+1900), (info.tm_mon+1), info.tm_mday, info.tm_hour, info.tm_min);

	stringstream ss, ss1;
	ss << root_path << "test/imu.csv";
	ss1 << log_dir << "/imu_raw" << ".log";
	FIO.fileSafeIn(ss.str(), fimu);

	if(logs){
		out.open(ss1.str());
	}

	OBSgnss.readEpoch(fst_pos);

	int err = pthread_create(&imu_thread_id, NULL, &loop_imu_thread, NULL);

	if(err != 0){
		printf("ERRORE THREAD\n");
		//TODO: close_semor(1)
	}

	read_imu(); //fill OBSimu

	_epochIMU = OBSimu._IMUdata.imuTime;
	_epochGNSS = fst_pos.time.sec;

	// Initial ECEF position for vehicle from GNSS     It puts the first position of the gnss file (first line) in _ECEF_o (and in _ECEF_imu e in GNSSsol.posXYZ)
	_ECEF_o = eigVector2std(double2eigVector(fst_pos.a, fst_pos.b, fst_pos.c));
	_ECEF_imu = _ECEF_o; GNSSsol.posXYZ = _ECEF_o;
	GNSSsol.velXYZ = eigVector2std(double2eigVector(fst_pos.va, fst_pos.vb, fst_pos.vc));

	_epochTime = _epochGNSS;

	// Initial Position in Geodetic and ENU
	Vector3d h;
	h << _ECEF_o.at(0), _ECEF_o.at(1), _ECEF_o.at(2);
	_LLH_o = eigVector2std(ecef2geo(h));

	IMU_INI_TIME_END = _epochIMU+imu_init_epochs; // Time taken to initialize the imu  (first imu epoch + 300) (for example)

	memset(imu_data, 0, IMUBUF_CAPACITY*IMU_LENGTH);

}
