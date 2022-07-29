#include "semor.h"
#include "Loosely.h"

/*
This file put c and c++ sources in communication
*/
Loosely imu;

void imu_sol(gnss_sol_t* cur_gnss){ //called every second at start for initializing IMU (calculating biases and errors) and then to take the solution
    imu.get_imu_sol(cur_gnss);
}

void init_imu(gnss_sol_t fst_pos){
    imu.init_imu(fst_pos);
}

void close_ctocpp(void){
    imu.close_out_file();
}
