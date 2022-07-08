#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <errno.h>
#include <unistd.h>
#include <fcntl.h>
#include <linux/input.h>
#include <linux/i2c-dev.h>
#include <sys/time.h>
#include "semor.h"
#include <time.h>

// https://www.sparkfun.com/products/18020
// https://github.com/sparkfun/SparkFun_Qwiic_6DoF_LSM6DSO_Arduino_Library/blob/3e126baaa3b2a4874af9bc6ab9c57712e0395ee2/src/SparkFunLSM6DSO.h#L252
enum LSM6DSO_REGISTERS
{
	FUNC_CFG_ACCESS			= 0x01,
	LSM6DO_PIN_CTRL			= 0x02,

	FIFO_CTRL1				= 0x07,
	FIFO_CTRL2				= 0x08,
	FIFO_CTRL3				= 0x09,
	FIFO_CTRL4				= 0x0A,

	COUNTER_BDR_REG1		= 0x0B,
	COUNTER_BDR_REG2		= 0x0C,

	INT1_CTRL				= 0x0D,
	INT2_CTRL				= 0x0E,
	WHO_AM_I_REG			= 0x0F,
	CTRL1_XL				= 0x10,
	CTRL2_G					= 0x11,
	CTRL3_C					= 0x12,
	CTRL4_C					= 0x13,
	CTRL5_C					= 0x14,
	CTRL6_C					= 0x15,
	CTRL7_G					= 0x16,
	CTRL8_XL				= 0x17,
	CTRL9_XL				= 0x18,
	CTRL10_C				= 0x19,
	ALL_INT_SRC				= 0x1A,
	WAKE_UP_SRC				= 0x1B,
	TAP_SRC					= 0x1C,
	D6D_SRC					= 0x1D,
	STATUS_REG				= 0x1E,
	OUT_TEMP_L				= 0x20,
	OUT_TEMP_H				= 0x21,
	OUTX_L_G				= 0x22,
	OUTX_H_G				= 0x23,
	OUTY_L_G				= 0x24,
	OUTY_H_G				= 0x25,
	OUTZ_L_G				= 0x26,
	OUTZ_H_G				= 0x27,

	OUTX_L_A				= 0x28,
	OUTX_H_A				= 0x29,
	OUTY_L_A				= 0x2A,
	OUTY_H_A				= 0x2B,
	OUTZ_L_A				= 0x2C,
	OUTZ_H_A				= 0x2D,

	EMB_FUNC_STATUS_MP		= 0x35,
	FSM_FUNC_STATUS_A_MP	= 0x36,
	FSM_FUNC_STATUS_B_MP	= 0x37,
	STATUS_MASTER_MAINPAGE	= 0x39,

	FIFO_STATUS1			= 0x3A,
	FIFO_STATUS2			= 0x3B,

	TIMESTAMP0_REG			= 0x40,
	TIMESTAMP1_REG			= 0x41,
	TIMESTAMP2_REG			= 0x42,
	TIMESTAMP3_REG			= 0x43,

	TAP_CFG0				= 0x56,
	TAP_CFG1				= 0x57,
	TAP_CFG2				= 0x58,
	TAP_THS_6D				= 0x59,
	INT_DUR2				= 0x5A,
	WAKE_UP_THS				= 0x5B,
	WAKE_UP_DUR				= 0x5C,
	FREE_FALL				= 0x5D,
	MD1_CFG					= 0x5E,
	MD2_CFG					= 0x5F,

	I3C_BUS_AVB				= 0x62,
	INTERNAL_FREQ_FINE		= 0x63,


	INT_OIS					= 0x6F,
	CTRL1_OIS				= 0x70,
	CTRL2_OIS				= 0x71,
	CTRL3_OIS				= 0x72,
	X_OFS_USR				= 0x73,
	Y_OFS_USR				= 0x74,
	Z_OFS_USR				= 0x75,

	FIFO_DATA_OUT_TAG		= 0x78,
	FIFO_DATA_OUT_X_L		= 0x79,
	FIFO_DATA_OUT_X_H		= 0x7A,
	FIFO_DATA_OUT_Y_L		= 0x7B,
	FIFO_DATA_OUT_Y_H		= 0x7C,
	FIFO_DATA_OUT_Z_L		= 0x7D,
	FIFO_DATA_OUT_Z_H		= 0x7E,
};

enum LSM6DSO_FS_XL
{
	FS_XL_2g			= 0x00,
	FS_XL_16g			= 0x04,
	FS_XL_4g			= 0x08,
	FS_XL_8g			= 0x0C,
	FS_XL_MASK			= 0xF3
};
enum LSM6DSO_ODR_XL
{
	ODR_XL_DISABLE		= 0x00,
	ODR_XL_1_6Hz		= 0xB0, // Low Power only
	ODR_XL_12_5Hz		= 0x10, // Low Power only
	ODR_XL_26Hz			= 0x20, // Low Power only
	ODR_XL_52Hz			= 0x30, // Low Power only
	ODR_XL_104Hz		= 0x40, // Normal Mode
	ODR_XL_208Hz		= 0x50, // Normal Mode
	ODR_XL_416Hz		= 0x60, // High performance
	ODR_XL_833Hz		= 0x70, // High Performance
	ODR_XL_1660Hz		= 0x80, // High Performance
	ODR_XL_3330Hz		= 0x90, // High Performance
	ODR_XL_6660Hz		= 0xA0, // High Performance
	ODR_XL_MASK			= 0x0F
};
enum LSM6DSO_FS_G
{
	FS_G_125dps			= 0x02,
	FS_G_250dps			= 0x00,
	FS_G_500dps			= 0x04,
	FS_G_1000dps		= 0x08,
	FS_G_2000dps		= 0x0C,
	FS_G_MASK			= 0xF0
};
enum LSM6DSO_ODR_GYRO_G
{
	ODR_GYRO_DISABLE	= 0x00,
	ODR_GYRO_12_5Hz		= 0x10, // Low Power only
	ODR_GYRO_26Hz		= 0x20, // Low Power only
	ODR_GYRO_52Hz		= 0x30, // Low Power only
	ODR_GYRO_104Hz		= 0x40, // Normal Mode
	ODR_GYRO_208Hz		= 0x50, // Normal Mode
	ODR_GYRO_416Hz		= 0x60, // High performance
	ODR_GYRO_833Hz		= 0x70, // High Performance
	ODR_GYRO_1660Hz		= 0x80, // High Performance
	ODR_GYRO_3330Hz		= 0x90, // High Performance
	ODR_GYRO_6660Hz		= 0xA0, // High Performance
	ODR_GYRO_MASK		= 0x0F
};
enum LSM6DSO_SW_RESET
{
	SW_RESET_NORMAL_MODE= 0x00,
	SW_RESET_DEVICE		= 0x01,
};
enum LSM6DSO_IF_INC
{
	IF_INC_DISABLED		= 0x00,
	IF_INC_ENABLED		= 0x04,
};
enum LSM6DSO_BDU
{
	BDU_CONTINUOS		= 0x00,
	BDU_BLOCK_UPDATE	= 0x40,
	BDU_MASK			= 0xBF
};
enum LSM6DSO_BOOT
{
	BOOT_NORMAL_MODE	= 0x00,
	BOOT_REBOOT_MODE	= 0x80,
};
enum LSM6DSO_TIMESTAMP
{
	TIMESTAMP_ENABLED	= 0x20,
};
enum LSM6DSO_BDR_GY_FIFO
{
	FIFO_BDR_GYRO_NOT_BATCHED	= 0x00,
	FIFO_BDR_GYRO_12_5Hz		= 0x10,
	FIFO_BDR_GYRO_26Hz			= 0x20,
	FIFO_BDR_GYRO_52Hz			= 0x30,
	FIFO_BDR_GYRO_104Hz			= 0x40,
	FIFO_BDR_GYRO_208Hz			= 0x50,
	FIFO_BDR_GYRO_417Hz			= 0x60,
	FIFO_BDR_GYRO_833Hz			= 0x70,
	FIFO_BDR_GYRO_1667Hz		= 0x80,
	FIFO_BDR_GYRO_3333Hz		= 0x90,
	FIFO_BDR_GYRO_6667Hz		= 0xA0,
	FIFO_BDR_GYRO_6_5Hz			= 0xB0,
	FIFO_BDR_GYRO_MASK			= 0x0F
};
enum LSM6DSO_BDR_XL_FIFO
{
	FIFO_BDR_ACC_NOT_BATCHED   = 0x00,
	FIFO_BDR_ACC_12_5Hz		   = 0x01,
	FIFO_BDR_ACC_26Hz		   = 0x02,
	FIFO_BDR_ACC_52Hz		   = 0x03,
	FIFO_BDR_ACC_104Hz		   = 0x04,
	FIFO_BDR_ACC_208Hz		   = 0x05,
	FIFO_BDR_ACC_417Hz		   = 0x06,
	FIFO_BDR_ACC_833Hz		   = 0x07,
	FIFO_BDR_ACC_1667Hz		   = 0x08,
	FIFO_BDR_ACC_3333Hz		   = 0x09,
	FIFO_BDR_ACC_6667Hz		   = 0x0A,
	FIFO_BDR_ACC_1_6Hz		   = 0x0B,
	FIFO_BDR_ACC_MASK		   = 0xF0
};
enum LSM6DSO_FIFO_TS_DEC
{
	FIFO_TS_DEC_DISABLED = 0x00,
	FIFO_TS_DEC_BY_1	 = 0x40,
	FIFO_TS_DEC_BY_8	 = 0x80,
	FIFO_TS_DEC_BY_32	 = 0xC0,
};
enum LSM6DSO_FIFO_MODE
{
	FIFO_MODE_DISABLED		 = 0x00,
	FIFO_MODE_STOP_WHEN_FULL = 0x01,
	FIFO_MODE_CONT_TO_FIFO	 = 0x03,
	FIFO_MODE_BYPASS_TO_CONT = 0x04,
	FIFO_MODE_CONTINUOUS	 = 0x06,
	FIFO_MODE_BYPASS_TO_FIFO = 0x07,
	FIFO_MODE_MASK			 = 0xF0
};

enum LSM6DSO_FIFO_TAGS
{
	GYROSCOPE_DATA = 0x01,
	ACCELEROMETER_DATA,
	TEMPERATURE_DATA,
	TIMESTAMP_DATA,
	CFG_CHANGE_DATA,
	ACCELERTOMETER_DATA_T_2,
	ACCELERTOMETER_DATA_T_1,
	ACCELERTOMETER_DATA_2xC,
	ACCELERTOMETER_DATA_3xC,
	GYRO_DATA_T_2,
	GYRO_DATA_T_1,
	GYRO_DATA_2xC,
	GYRO_DATA_3xC,
};

// https://raspberry-projects.com/pi/programming-in-c/i2c/using-the-i2c-interface
int i2c_fd;
uint8_t i2c_buf[16];
int i2c_open()
{
	i2c_fd = open("/dev/i2c-1", O_RDWR);
	if (i2c_fd < 0)
		perror("Could not open /dev/i2c-1");
	return i2c_fd;
}
int i2c_addr(uint8_t addr7)
{
	int e = ioctl(i2c_fd, I2C_SLAVE, addr7);
	if (e < 0)
		perror("Could not set slave address");
	return e;
}
int i2c_read(int length)
{
	int n = read(i2c_fd, i2c_buf, length);
	if (n != length)
	{
		if (n < 0)
			perror("Could not read from i2c device");
		else fprintf(stderr, "Could not read %d bytes, got %d\n", length, n);
	}
	return n;
}
int i2c_write(int length)
{
	int n = write(i2c_fd, i2c_buf, length);
	if (n != length)
	{
		if (n < 0)
			perror("Could not write to i2c device");
		else fprintf(stderr, "Could not write %d bytes, did %d\n", length, n);
	}
	return n;
}

#define USE_FIFO 1
#define ADD_STAMP 1

uint32_t last_stamp = 0;

int week;
double sec = 0;
int first_read = 1;
int need_setup = 1;

double avg_acc_x = 0;
double avg_acc_y = 0;
double avg_acc_z = 0;
int nsamples = 0;

double tot_acc_x = 0;
double tot_acc_y = 0;
double tot_acc_z = 0;

double gravity = 9.80665;

int setup(){
	//	open /dev/i2c-1 and set slave address
	if (i2c_open() < 0)
		return 1;
	if (i2c_addr(0x6B) < 0)
		return 2;

//	setup registers
	int setup_guard = 0;
	for (;;)
	{
		if (++setup_guard > 3)
		{
			fputs("Setup failed\n", stderr);
			return 3;
		}

		// test
		i2c_buf[0] = WHO_AM_I_REG;
		if (i2c_write(1) != 1)
			continue;
		if (i2c_read(1) != 1)
			continue;
		if (i2c_buf[0] != 0x6C)
			fprintf(stderr, "WHO AM I: %02X!\n", i2c_buf[0]);

		// reset
		i2c_buf[0] = CTRL3_C;
		i2c_buf[1] = IF_INC_ENABLED | BOOT_REBOOT_MODE | SW_RESET_DEVICE;
		if (i2c_write(2) != 2)
			continue;

		// setup or disable FIFO
		i2c_buf[0] = FIFO_CTRL1;
	#if USE_FIFO
		i2c_buf[1] = 0;
		i2c_buf[2] = 0;
		i2c_buf[3] = FIFO_BDR_ACC_104Hz | FIFO_BDR_GYRO_104Hz;
	#if ADD_STAMP
		i2c_buf[4] = FIFO_MODE_CONTINUOUS | FIFO_TS_DEC_BY_1;
	#else
		i2c_buf[4] = FIFO_MODE_CONTINUOUS;
	#endif
	#else
		i2c_buf[1] = 0;
		i2c_buf[2] = 0;
		i2c_buf[3] = 0;
		i2c_buf[4] = FIFO_MODE_DISABLED;
	#endif
		if (i2c_write(5) != 5)
			continue;

		// start 8g, 500dps, 104Hz
		i2c_buf[0] = CTRL1_XL;
		i2c_buf[1] = FS_XL_8g|ODR_XL_104Hz;
		i2c_buf[2] = FS_G_500dps|ODR_GYRO_104Hz;
	#if USE_FIFO
		i2c_buf[3] = IF_INC_ENABLED;
	#else
		i2c_buf[3] = IF_INC_ENABLED | BDU_BLOCK_UPDATE;
	#endif
		if (i2c_write(4) != 4)
			continue;

		i2c_buf[0] = CTRL10_C;
	#if ADD_STAMP
		i2c_buf[1] = TIMESTAMP_ENABLED;
	#else
		i2c_buf[1] = 0;
	#endif
		if (i2c_write(2) != 2)
			continue;

		break;
	}
}

//	unit conversion constants
	double acclUnits = 1.0/4096;// bit weight for 8g range (8.0/(1<<(16-1)))
	double gyroUnits = 0.0175;	// bit weight for 500dps range

//	read and print data
	#if USE_FIFO
	int16_t ax,ay,az;
	int16_t gx,gy,gz;
	uint32_t stamp;
	uint8_t flags;
	#endif

int get_imu_data_old(char buf[IMU_LENGTH]){
	if(need_setup){
		setup();
		need_setup = 0;
	}

	i2c_buf[0] = FIFO_STATUS1;
	if (i2c_write(1) != 1)
		return 1;
	if (i2c_read(2) != 2)
		return 1;
//	number of unread data
	int n = i2c_buf[0] | (i2c_buf[1] & 3) << 8;
	if(n == 0)
		return 1; //try again, no data right now
	i2c_buf[0] = FIFO_DATA_OUT_TAG;
	if (i2c_write(1) != 1)
		return 1;
	if (i2c_read(7) != 7)
		return 1;
	//printf("%02X:%02X%02X:%02X%02X:%02X%02X\n", i2c_buf[0], i2c_buf[1], i2c_buf[2], i2c_buf[3], i2c_buf[4], i2c_buf[5], i2c_buf[6]);
	switch (i2c_buf[0] >> 3)
	{
	case ACCELEROMETER_DATA:
	case ACCELERTOMETER_DATA_T_1:
	case ACCELERTOMETER_DATA_T_2:
	case ACCELERTOMETER_DATA_2xC:
	case ACCELERTOMETER_DATA_3xC:
		ax = i2c_buf[1] | i2c_buf[2] << 8;
		ay = i2c_buf[3] | i2c_buf[4] << 8;
		az = i2c_buf[5] | i2c_buf[6] << 8;
		flags |= 1;
		break;
	case GYROSCOPE_DATA:
	case GYRO_DATA_T_1:
	case GYRO_DATA_T_2:
	case GYRO_DATA_2xC:
	case GYRO_DATA_3xC:
		gx = i2c_buf[1] | i2c_buf[2] << 8;
		gy = i2c_buf[3] | i2c_buf[4] << 8;
		gz = i2c_buf[5] | i2c_buf[6] << 8;
		flags |= 2;
		break;
	case TIMESTAMP_DATA:
		stamp = i2c_buf[1] | i2c_buf[2] << 8 | i2c_buf[3] << 16 | i2c_buf[4] << 24;
		flags |= 4;
		break;
	default:
		printf("TAG %02X\n", i2c_buf[0]);
		break;
	}
	if (flags == 7)
	{
		int week1;
		double sec1;
		flags = 0;
		if(first_read){
			struct timeval tv;
			gettimeofday(&tv, NULL);
			week = ((tv.tv_sec+LEAP_SECONDS-GPS_EPOCH))/(7*24*3600);
			sec = (double)(((tv.tv_sec+LEAP_SECONDS-GPS_EPOCH))%(7*24*3600))+(tv.tv_usec / 1000000.0);
			first_read = 0;
		}
		else{
			sec += (stamp-last_stamp) * 2.5e-5;
			struct timeval tvv;
			gettimeofday(&tvv, NULL);
			week1 = ((tvv.tv_sec+LEAP_SECONDS-GPS_EPOCH))/(7*24*3600);
			sec1 = (double)(((tvv.tv_sec+LEAP_SECONDS-GPS_EPOCH))%(7*24*3600))+(tvv.tv_usec / 1000000.0);
		}
		/*nsamples += 1;
		tot_acc_x += ax*acclUnits*gravity;
		tot_acc_y += ay*acclUnits*gravity;
		tot_acc_z += az*acclUnits*gravity;*/
		sprintf(buf, "%d %lf %7.4f %7.4f %7.4f %8.3f %8.3f %8.3f | diff: %10u | n: %d | preso: %lf", week,
			sec,
			/*ax,ay,az, gx,gy,gz,*/
			ax*acclUnits*gravity/*-avg_acc_x*/, ay*acclUnits*gravity/*-avg_acc_y*/, az*acclUnits*gravity/*-avg_acc_z*/,
			gx*gyroUnits, gy*gyroUnits, gz*gyroUnits, stamp-last_stamp, n, sec1);

		last_stamp = stamp;

		return 0;

		//printf("%10u: %6d %6d %6d | %6d %6d %6d || %7.4f %7.4f %7.4f | %8.3f %8.3f %8.3f\n",

	}
	return 1;
}

int get_imu_data(char line[IMU_LENGTH]){
	int n;
	if(need_setup){
		setup();
		need_setup = 0;
	}
	while(1){
		i2c_buf[0] = FIFO_STATUS1;
		if (i2c_write(1) != 1)
			continue;
		if (i2c_read(2) != 2)
			continue;
	//	number of unread data
		n = i2c_buf[0] | (i2c_buf[1] & 3) << 8;
		if(n < 3)
			continue; //try again, no data right now
		break;
	}
	while(1){ //exit when return
		i2c_buf[0] = FIFO_DATA_OUT_TAG;
		if (i2c_write(1) != 1)
			continue;
		if (i2c_read(7) != 7)
			continue;
		//printf("%02X:%02X%02X:%02X%02X:%02X%02X\n", i2c_buf[0], i2c_buf[1], i2c_buf[2], i2c_buf[3], i2c_buf[4], i2c_buf[5], i2c_buf[6]);
		switch (i2c_buf[0] >> 3)
		{
		case ACCELEROMETER_DATA:
		case ACCELERTOMETER_DATA_T_1:
		case ACCELERTOMETER_DATA_T_2:
		case ACCELERTOMETER_DATA_2xC:
		case ACCELERTOMETER_DATA_3xC:
			ax = i2c_buf[1] | i2c_buf[2] << 8;
			ay = i2c_buf[3] | i2c_buf[4] << 8;
			az = i2c_buf[5] | i2c_buf[6] << 8;
			flags |= 1;
			break;
		case GYROSCOPE_DATA:
		case GYRO_DATA_T_1:
		case GYRO_DATA_T_2:
		case GYRO_DATA_2xC:
		case GYRO_DATA_3xC:
			gx = i2c_buf[1] | i2c_buf[2] << 8;
			gy = i2c_buf[3] | i2c_buf[4] << 8;
			gz = i2c_buf[5] | i2c_buf[6] << 8;
			flags |= 2;
			break;
		case TIMESTAMP_DATA:
			stamp = i2c_buf[1] | i2c_buf[2] << 8 | i2c_buf[3] << 16 | i2c_buf[4] << 24;
			flags |= 4;
			break;
		default:
			printf("TAG %02X\n", i2c_buf[0]);
			break;
		}
		if (flags == 7)
		{
			flags = 0;

			struct timeval tv;
			gettimeofday(&tv, NULL);
			week = ((tv.tv_sec+LEAP_SECONDS-GPS_EPOCH))/(7*24*3600);
			sec = (double)(((tv.tv_sec+LEAP_SECONDS-GPS_EPOCH))%(7*24*3600))+(tv.tv_usec / 1000000.0);

			/*nsamples += 1;
			tot_acc_x += ax*acclUnits*gravity;
			tot_acc_y += ay*acclUnits*gravity;
			tot_acc_z += az*acclUnits*gravity;*/
			sprintf(line, "%d %lf %7.4f %7.4f %7.4f %8.3f %8.3f %8.3f | %10u %d", week,
				sec,
				/*ax*acclUnits*gravity, ay*acclUnits*gravity, az*acclUnits*gravity,
				gx*gyroUnits, gy*gyroUnits, gz*gyroUnits, stamp, n);*/
        -ax*acclUnits*gravity, -ay*acclUnits*gravity, -az*acclUnits*gravity,
				gx*gyroUnits*3.1415926535898/180.0, gy*gyroUnits*3.1415926535898/180.0, gz*gyroUnits*3.1415926535898/180.0, stamp, n);

			last_stamp = stamp;

			return 0;

			//printf("%10u: %6d %6d %6d | %6d %6d %6d || %7.4f %7.4f %7.4f | %8.3f %8.3f %8.3f\n",

		}
	}
	return 1;
}

void imu_log(void){ //place the call to this function in the flow after the initialization of str_time in client.c (start_processing())

    char path[PATH_MAX];

    char buf[300];

    FILE* f;

    time_t rawtime;
    struct tm info;
    time( &rawtime );
    info = *localtime( &rawtime );
    char str_time[13];

    sprintf(str_time, "%d_%02d_%02d_%02d_%02d", (info.tm_year+1900), (info.tm_mon+1), info.tm_mday, info.tm_hour, info.tm_min);

    struct stat st = {0};
    log_dir[PATH_MAX/2];

    if(logs){
        if(relative){
            sprintf(log_dir, "logs/logs_%s", str_time);
            if (stat(log_dir, &st) == -1) {
                mkdir(log_dir, 0777);
            }
        }
        else{
            sprintf(log_dir, "%slogs/logs_%s", root_path, str_time);
            if (stat(log_dir, &st) == -1) {
                mkdir(log_dir, 0777);
            }
        }
        sprintf(path, "%s/imu_%s.log", log_dir, str_time);
    }
    f = fopen(path, "w");

    while(1){
        get_imu_data(buf);

        fprintf(f, "%s\n", buf);
        
        fflush(f);
    }
}