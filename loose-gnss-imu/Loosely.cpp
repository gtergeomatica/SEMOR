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

ofstream out;

int first = 1;

// PI
const double PI = 3.1415926535898;

double epochAfterRead;

double lastGroupReadEpoch;


IMUmechECEF MechECEF;

ifstream fimu;

InitializeIMU iniIMU;
double IMU_INI_TIME_END;

ReaderGNSS OBSgnss;
ReaderIMU OBSimu;


int last_imu_epoch = 0;

pthread_mutex_t lock;
pthread_t imu_thread_id;
char imu_data[IMUBUF_CAPACITY][IMU_LENGTH];
int imu_count_write = 0;
int imu_count_read = 0;

Loosely::Loosely(){


}

void Loosely::close_out_file(){
	out.close();
}

void print_time(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	int week = ((tv.tv_sec+LEAP_SECONDS-GPS_EPOCH))/(7*24*3600);
	double sec = (double)(((tv.tv_sec+LEAP_SECONDS-GPS_EPOCH))%(7*24*3600))+(tv.tv_usec / 1000000.0);
	printf("%lf\n", sec);
}

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
			strcpy(buf, imu_data[imu_count_read]); //forse qui non serve il lock
			pthread_mutex_unlock(&lock);
			if(strlen(buf) != 0){
				//printf("%s\n", buf);
				line = string(buf);
				OBSimu.clearObs();
				OBSimu.obsEpoch(line);
				//printf("%lf | %lf\n", OBSimu._IMUdata.imuTime, _epochIMU);
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

void closeF(){
	out.close();
}


// A routine to facilitate IMU mechanization in ECEF
void Loosely::SolutionIMU(ReaderIMU IMU, IMUmechECEF& MechECEF) {
	// Compute time interval
	_dTimu = IMU._IMUdata.imuTime - _epochIMU; //Difference between current epoch and previous
	//printf("%lf\n", _dTimu);
	// IMU Mechanization
  //cout << "2-->" << iniIMU._RPY.at(0)*180.0/PI << " | " << iniIMU._RPY.at(0)*180.0/PI << " | " << iniIMU._RPY.at(0)*180.0/PI << endl;
	MechECEF.MechanizerECEF(_dTimu, IMU._IMUdata.Acc, IMU._IMUdata.Gyr, _LLH_o); //Update ECEF position adding to it accelerometer and gyroscope data (previous data is accumulated)
	// Update solution
	_epochIMU = IMU._IMUdata.imuTime; //Update epocIMU with current epoch
	IMUsol.posXYZ = MechECEF._pos;		//Update IMU solution
	IMUsol.velXYZ = MechECEF._vel;		//Update IMU solution
	IMUsol.attXYZ = MechECEF._att;		//Update IMU solution
	_Heading_imu = normalise(IMUsol.attXYZ.at(2), 0, 2 * PI);
	//printf("%lf\n", IMU._IMUdata.imuTime);
}

VectorXd double2eigVector(double a, double b, double c){
	VectorXd v = VectorXd::Zero(3);
	v(0) = a;
	v(1) = b;
	v(2) = c;
	return v;
}

int deb_enabled = 1;
void deb(char message[300]){
    if(deb_enabled){
        printf("%s\n", message);
    }
}

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
		//while(get_imu_data(buf) != 0){} //old version (for get_imu_data_old(buf))

		//modificare matrice "imu_data"
		pthread_mutex_lock(&lock);
		strcpy(imu_data[imu_count_write++], buf);
		pthread_mutex_unlock(&lock);
		imu_count_write = imu_count_write % IMUBUF_CAPACITY;
		//printf("imu_count_write: %d\n", imu_count_write);
		//pthread_mutex_unlock(&lock);
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
	gnss.a = phi_i;
	gnss.b = lambda;
	gnss.c = h;
	return gnss;
}

void Loosely::get_imu_sol(gnss_sol_t* int_sol){
	 /* struct timeval tv;
            gettimeofday(&tv, NULL);
            int week = ((tv.tv_sec+LEAP_SECONDS-GPS_EPOCH))/(7*24*3600);
            double sec = (double)(((tv.tv_sec+LEAP_SECONDS-GPS_EPOCH))%(7*24*3600))+(tv.tv_usec / 1000000.0);
	printf("to process: %d | now: %lf\n", (*int_sol).time.sec, sec); */
	//Check if imu is initializing
	if(imu_ready == 0){
		int n = 0;
		while(/* n < 104 */ OBSimu._IMUdata.imuTime < (*int_sol).time.sec){ //read whole imu date of a particular gnss epoch
			if(imu_ready == 0 && iniIMU.stepInitializeIMU(OBSimu, IMU_INI_TIME_END, _LLH_o) == 1){ //it calculates _ACCbias and _GYRbias, _RPY
				//_dT = 1.0;
				// Initialize IMU Mechanization
        //cout << "1-->" << iniIMU._RPY.at(0)*180.0/PI << " | " << iniIMU._RPY.at(0)*180.0/PI << " | " << iniIMU._RPY.at(0)*180.0/PI << endl;
				MechECEF.InitializeMechECEF(_ECEF_imu, _LLH_o, GNSSsol.velXYZ, iniIMU._RPY, iniIMU._ACCbias, iniIMU._GYRbias); //it initializes MechECEF with _ACCbias, _GYRbias and _RPY

				imu_ready = 1; //Tells client.c that it can read imu data
				printf("SEMOR: End initialization\n");
			}
			_epochIMU = OBSimu._IMUdata.imuTime;
			lastGroupReadEpoch = OBSimu._IMUdata.imuTime;
			read_imu();
			if(imu_ready == 1){
				imu_ready = 2;
				(*int_sol).time.week = OBSimu._IMUdata.week;
				/*avg_acc_x = tot_acc_x/((double)nsamples);
				avg_acc_y = tot_acc_y/((double)nsamples);
				avg_acc_z = tot_acc_z/((double)nsamples);

				printf("%lf - %lf - %lf\n", avg_acc_x, avg_acc_y, avg_acc_z);*/
			}
			epochAfterRead = OBSimu._IMUdata.imuTime;
			n++;
		}
		//cout << n << endl;
		//cout.flush();
		return;
	}
	//IMU is ready:

	//Initialize IMU mechanization with initial position, initial velocity and biases
	/*_ECEF_o = eigVector2std(double2eigVector((*int_sol).a, (*int_sol).b, (*int_sol).c));
	_LLH_o = eigVector2std(ecef2geo(double2eigVector((*int_sol).a, (*int_sol).b, (*int_sol).c)));
	_ECEF_imu = _ECEF_o;
	GNSSsol.velXYZ = eigVector2std(double2eigVector((*int_sol).va, (*int_sol).vb, (*int_sol).vc));*/

	MechECEF._pos.at(0) = (*int_sol).a;
	MechECEF._pos.at(1) = (*int_sol).b;
	MechECEF._pos.at(2) = (*int_sol).c;

	MechECEF._vel.at(0) = (*int_sol).va;
	MechECEF._vel.at(1) = (*int_sol).vb;
	MechECEF._vel.at(2) = (*int_sol).vc;

	/*if(1){
		MechECEF.InitializeMechECEF(_ECEF_imu, _LLH_o, GNSSsol.velXYZ, iniIMU._RPY, iniIMU._ACCbias, iniIMU._GYRbias);
		first = 0;
	}*/
	//here we have previous gnss position



	//print_time();

	double group_time = 0.1;
	double next_time = epochAfterRead + group_time;
	//epochAfterRead = _epochIMU;

	double avg_ax = 0, avg_ay = 0, avg_az = 0;
	double avg_gx = 0, avg_gy = 0, avg_gz = 0;
	int n = 0, nlocal = 0;
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

		/*avg_ax +=OBSimu._IMUdata.Acc.at(0);
		avg_ay +=OBSimu._IMUdata.Acc.at(1);
		avg_az +=OBSimu._IMUdata.Acc.at(2);

		avg_gx +=OBSimu._IMUdata.Gyr.at(0) ;
		avg_gy +=OBSimu._IMUdata.Gyr.at(1) ;
		avg_gz +=OBSimu._IMUdata.Gyr.at(2) ;

		n++;
		nlocal++;
		//printf("a\n");

		if(epochAfterRead >= next_time || epochAfterRead > (*int_sol).time.sec){

			avg_ax /= nlocal;
			avg_ay /= nlocal;
			avg_az /= nlocal;

			avg_gx /= nlocal;
			avg_gy /= nlocal;
			avg_gz /= nlocal;

			OBSimu._IMUdata.Acc.at(0) = avg_ax;
			OBSimu._IMUdata.Acc.at(1) = avg_ay;
			OBSimu._IMUdata.Acc.at(2) = avg_az;

			OBSimu._IMUdata.Gyr.at(0) = avg_gx;
			OBSimu._IMUdata.Gyr.at(1) = avg_gy;
			OBSimu._IMUdata.Gyr.at(2) = avg_gz;

			gnss_sol_t llh = ecef2geo(*int_sol);
			_LLH_o = eigVector2std(double2eigVector(llh.a, llh.b, llh.c));

      //cout << avg_gz << endl ;

			SolutionIMU(OBSimu, MechECEF);
			(*int_sol).a = IMUsol.posXYZ.at(0);
			(*int_sol).b = IMUsol.posXYZ.at(1);
			(*int_sol).c = IMUsol.posXYZ.at(2);
			(*int_sol).va = IMUsol.velXYZ.at(0);
			(*int_sol).vb = IMUsol.velXYZ.at(1);
			(*int_sol).vc = IMUsol.velXYZ.at(2);

			lastGroupReadEpoch = OBSimu._IMUdata.imuTime;

      //cout << "gx : " << avg_gx << " gy = " << avg_gy << " gz = " << avg_gz << endl;
      //cout << "ax : " << avg_ax << " ay = " << avg_ay << " az = " << avg_az << endl;

			//Reset variables for next iterations
			next_time += group_time;
			avg_ax = avg_ay = avg_az = 0;
			avg_gx = avg_gy = avg_gz = 0;
			nlocal = 0;
		}
		//Update Time
		_epochIMU = OBSimu._IMUdata.imuTime;
		read_imu();
		epochAfterRead = OBSimu._IMUdata.imuTime;*/
	} while (epochAfterRead <= (*int_sol).time.sec);

		 /*cout << n << endl;
		cout.flush();
cout << "---------------------------------------------------" << endl;
		cout.flush();*/


	/* if(logs){
		cout << n << endl;
		cout.flush();
		out << "---------------------------------------------------" << endl;
		out.flush();
	} */

	//printf("next_time: %lf\n", next_time);

	/* (*int_sol).a = IMUsol.posXYZ.at(0);
	(*int_sol).b = IMUsol.posXYZ.at(1);
	(*int_sol).c = IMUsol.posXYZ.at(2);

	//printf("%lf, %lf, %lf\n", (*int_sol).a, (*int_sol).b, (*int_sol).c);

  //cout << "-->" << IMUsol.velXYZ.at(0) << endl;

	(*int_sol).va = IMUsol.velXYZ.at(0);
	(*int_sol).vb = IMUsol.velXYZ.at(1);
	(*int_sol).vc = IMUsol.velXYZ.at(2); */

	//cout << n << endl;
	//cout.flush();

	//

	//here we have next gnss+imu position

	//We have position and STD of this epoch


	//Mark imu solution as usable for the comparison (in order to get the best solution for SEMOR)
	(*int_sol).time.week = OBSimu._IMUdata.week;
}

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
	if(relative){
		ss << "test/imu.csv";
		ss1 << log_dir << "/imu_raw" << ".log";
	}
	else{
		ss << root_path << "test/imu.csv";
		ss1 << log_dir << "/imu_raw" << ".log";
	}
	FIO.fileSafeIn(ss.str(), fimu);

	if(logs){
		out.open(ss1.str());
	}

	OBSgnss.readEpoch(fst_pos);

	int err = pthread_create(&imu_thread_id, NULL, &loop_imu_thread, NULL);

	if(err != 0){
		printf("ERRORE THREAD\n");
	}

	read_imu(); //fill OBSimu
	/* get_imu_data(buf);
	line = string(buf);
	OBSimu.clearObs();
	OBSimu.obsEpoch(line);
 */

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
