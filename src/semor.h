#ifndef SEMOR_H
#define SEMOR_H

#include <unistd.h>
#include <stdint.h>
#include <limits.h>
#include <pthread.h>

/*
#ifdef __cplusplus
#include "Loosely.h"
#endif
*/

#define MAXSTR 512
#define LEAP_SECONDS 18
#define IMU_LENGTH 200 //length of imu data string
#define GPS_EPOCH 315964800 /* 6th January 1980*/

extern char root_path[PATH_MAX-200];

extern int logs;
extern int debug;
extern double init_bg_unc;
extern double init_ba_unc;
extern double psd_gyro;
extern double psd_acce;     // acce noise PSD (m^2/s^3)  
extern double psd_bg;     // gyro bias random walk PSD (rad^2/s^3)
extern double psd_ba;
extern int sample_rate;
extern int imu_init_epochs;//in seconds
extern int imu_drift; //in seconds
extern char rtk_port_rtkrcv[6];
extern char ppp_port_rtkrcv[6];

extern double init_x;
extern double init_y;
extern double init_z;

extern double avg_acc_x;
extern double avg_acc_y;
extern double avg_acc_z;

extern double tot_acc_x;
extern double tot_acc_y;
extern double tot_acc_z;

extern int nsamples;


extern pthread_t imu_thread_id;
extern char log_dir[PATH_MAX/2];

#define IMUBUF_CAPACITY 312 //3*imu frequency

extern char imu_data[IMUBUF_CAPACITY][IMU_LENGTH];
extern int imu_count_write; //to store samples in imu matrix
extern int imu_count_read; //to read samples from imu matrix
extern pthread_mutex_t lock;


extern pid_t str2str_pid, rtkrcv1_pid, rtkrcv2_pid;

//extern char imu_data[IMU_DATA_LAST_N_SEC * IMU_HZ][75]; /* Last 2 seconds of imu data (imu send data at 100Hz and it's a string of about 75 chars)*/

typedef struct {
    double dLon;
    double dLat;
    double dAcc;
    double dGyro;
    double dHeigth;
    double dSdn;
    double dSde;
    double dSdu;
    /*double dRobInd;*/
    uint8_t ui8FixQual;
    uint8_t ui8TS[50+1]; //+1 \0
}LocData_t;

extern LocData_t get_data();

typedef struct tow{ /* Time Of Week */
    int week;
    int sec;
}tow_t;

typedef struct gnss_solution{ /* a = x | n ,    b = y | e ,     c = z | u*/
    tow_t time;
    double a; //x
    double b; //y
    double c; //z
    short Q;
    short ns;
    double sda; //sdx
    double sdb; //sdy
    double sdc; //sdz
    double sdab; //sdxy
    double sdbc; //sdyz
    double sdca; //sdxz
    float age;
    float ratio;
    double va; //velocity x
    double vb;//velocity y
    double vc;//velocity z

    //int compared; //If it has no matches for its epoch, it will be printed alone with no comparisons
} gnss_sol_t;

//IMU
extern gnss_sol_t first_pos;
extern int imu_ready;
extern gnss_sol_t best;


//Functions
void start_processing(void);
void close_semor(int status);
#ifdef __cplusplus
//extern Loosely *imu; //##CHANGED##
extern "C"{
#endif
int get_imu_data(char line[IMU_LENGTH]);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
//extern Loosely *imu; //##CHANGED##
extern "C"{
#endif
void init_imu(gnss_sol_t fst_pos);
void imu_sol(gnss_sol_t* cur_gnss); /* Called for each GNSS epoch*/
void close_ctocpp(void); //Delete loosely
#ifdef __cplusplus
}
#endif



#endif /* SEMOR_H */
