#include <stdio.h>
#include <signal.h>
#include <string.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "semor.h"

#define MAX_LINE 256 //max length of a SEMOR config's line

//Set default configuration
int logs=1;
int debug=0; //1: read from files (test/galileo.pos, test/gps.pos, test/imu.csv)

//START Used in IMU standard deviation calculation
double init_bg_unc = 2.42406840554768e-05;
double init_ba_unc = 0.048901633857000000;
double psd_gyro  = 3.38802348178723e-09;
double psd_acce      =     2.60420170553977e-06;     // acce noise PSD (m^2/s^3)  
double psd_bg        =     2.61160339323310e-14;     // gyro bias random walk PSD (rad^2/s^3)
double psd_ba        =     1.66067346797506e-09;
int sample_rate = 104; //in hz, sample rate of the IMU
//END Used in IMU standard deviation calculation

int imu_init_epochs =  20;//in seconds, in this interval SEMOR collects IMU data to initialize it by calculating precision errors and biases
int imu_drift = 60; //in seconds, if during this interval only IMU solutions are available, SEMOR will terminate and it needs to be re-initialized (probably
 //something is wrong in RTK or PPP solutions)

char rtk_port_rtkrcv[6] = "8090";
char ppp_port_rtkrcv[6] = "8091";

//Next 3 variables represent the initial position of the device, the are needed for the IMU initialization process (leave the PID in that position
//until the IMU is initialized)
double init_x;
double init_y;
double init_z;
int init_pos = 0; //if init_x, init_y and init_z are assigned with a value this value will be 7, if the value is not 7 SEMOR stops 

//Paths of str2str and rtkrcv executables
char str2str_path[PATH_MAX];
char rtkrcv_path[PATH_MAX];

//Path of RTK and PPP conf files
char rtkconf[PATH_MAX];
char pppconf[PATH_MAX];

//Path to the root SEMOR project directory
char root_path[PATH_MAX-200];

//-in argument of str2str
char str2str_in[60] = "serial://ttyACM0#ubx $tcpcli://192.168.2.91:8081";

//-out arguments of str2str
char str2str_out1_port[16] = "tcpsvr://:8085";
char str2str_out2_port[16] = "tcpsvr://:8086";


//Close standard input and output of the RTKLIB processes (created with fork)
void close_io(void){
    close(STDIN_FILENO);
    close(STDOUT_FILENO);
    close(STDERR_FILENO);
}

//Read one line of the config file
void read_conf_line(char line[MAX_LINE]){
    int i;
    char *token;
    char *eptr;
    char delim[8] = "$= \t"; // '$' indica l'inizio di un commento nel file di configurazione
    token = strtok(line, delim);

    if(strstr(token, "rtkrcv-port-rtk")){
        strcpy(rtk_port_rtkrcv, strtok(NULL, delim));
        if(rtk_port_rtkrcv[strlen(rtk_port_rtkrcv)-1] == '\n')
            rtk_port_rtkrcv[strlen(rtk_port_rtkrcv)-1] = '\0'; //remove new line
        return;
    }

    if(strstr(token, "rtkrcv-port-ppp")){
        strcpy(ppp_port_rtkrcv, strtok(NULL, delim));
        if(ppp_port_rtkrcv[strlen(ppp_port_rtkrcv)-1] == '\n')
            ppp_port_rtkrcv[strlen(ppp_port_rtkrcv)-1] = '\0';
        return;
    }

    if(strstr(token, "init-x")){
        init_x = strtod(strtok(NULL, delim), &eptr);
        if(init_x != 0){
            init_pos |= 1;
        }
        return;
    }
    if(strstr(token, "init-y")){
        init_y = strtod(strtok(NULL, delim), &eptr);
        if(init_y != 0){
            init_pos |= 2;
        }
        return;
    }
    if(strstr(token, "init-z")){
        init_z = strtod(strtok(NULL, delim), &eptr);
        if(init_z != 0){
            init_pos |= 4;
        }
        return;
    }

    if(strstr(token, "imu-drift")){
        imu_drift = atoi(strtok(NULL, delim));
        return;
    }
    
    if(strstr(token, "debug")){
        debug = atoi(strtok(NULL, delim));
        return;
    }
    if(strstr(token, "logs")){
        logs = atoi(strtok(NULL,delim));
        return;
    }
    if(strstr(token, "str2str-in")){
        strcpy(str2str_in, strtok(NULL, delim));
        if(str2str_in[strlen(str2str_in)-1] == '\n')
            str2str_in[strlen(str2str_in)-1] = '\0';
        return;
    }
    if(strstr(token, "str2str-out1-port")){
        strcpy(str2str_out1_port, strtok(NULL, delim));
        if(str2str_out1_port[strlen(str2str_out1_port)-1] == '\n')
            str2str_out1_port[strlen(str2str_out1_port)-1] = '\0';
        return;
    }
    if(strstr(token, "str2str-out2-port")){
        strcpy(str2str_out2_port, strtok(NULL, delim));
        if(str2str_out2_port[strlen(str2str_out2_port)-1] == '\n')
            str2str_out2_port[strlen(str2str_out2_port)-1] = '\0';
        return;
    }
    if(strstr(token, "str2str-path")){
        strcpy(str2str_path, strtok(NULL, delim));
        if(str2str_path[strlen(str2str_path)-1] == '\n')
            str2str_path[strlen(str2str_path)-1] = '\0';
        return;
    }
    if(strstr(token, "rtkrcv-path")){
        strcpy(rtkrcv_path, strtok(NULL, delim));
        if(rtkrcv_path[strlen(rtkrcv_path)-1] == '\n')
            rtkrcv_path[strlen(rtkrcv_path)-1] = '\0';
        return;
    }
    if(strstr(token, "rtkconf")){
        strcpy(rtkconf, strtok(NULL, delim));
        if(rtkconf[strlen(rtkconf)-1] == '\n')
            rtkconf[strlen(rtkconf)-1] = '\0';
        return;
    }
    if(strstr(token, "pppconf")){
        strcpy(pppconf, strtok(NULL, delim));
        if(pppconf[strlen(pppconf)-1] == '\n')
            pppconf[strlen(pppconf)-1] = '\0';
        return;
    }
    if(strstr(token, "sample-rate")){
        sample_rate = atoi(strtok(NULL, delim));
        return;
    }
    if(strstr(token, "imu-init-epochs")){
        imu_init_epochs = atoi(strtok(NULL, delim));
        return;
    }
    if(strstr(token, "init-bg-unc")){
        init_bg_unc = strtod(strtok(NULL, delim), &eptr);
        return;
    }
    if(strstr(token, "init-ba-unc")){
        init_ba_unc = strtod(strtok(NULL, delim), &eptr);
        return;
    }
    if(strstr(token, "psd-gyro")){
        psd_gyro = strtod(strtok(NULL, delim), &eptr);
        return;
    }
    if(strstr(token, "psd-acce")){
        psd_acce = strtod(strtok(NULL, delim), &eptr);
        return;
    }
    if(strstr(token, "psd-bg")){
        psd_bg = strtod(strtok(NULL, delim), &eptr);
        return;
    }
    if(strstr(token, "psd-ba")){
        psd_ba = strtod(strtok(NULL, delim), &eptr);
        return;
    }
}

int main(int argc, char *argv[]){

    char cmd;
    int imu_ready;

    //Setup from configuration
    int i;
    char line[MAX_LINE];
	char cwd[PATH_MAX-400];
    char semor_conf_path[PATH_MAX];
    
    //pids.txt file to manually kill previous executed processes if needed (in case someone forgots to soft-terminate SEMOR entering 'q')
    char pids_file[PATH_MAX];

    char path[60];
    {
        sprintf(path, "/proc/%d/exe", getpid());
        if(readlink(path, root_path, PATH_MAX) == -1){
            perror("SEMOR: readlink()");
            printf("\n");
            return 0;
        }
        for(i=strlen(root_path)-1; i >= 0; i--){
            if(root_path[i] != '/'){
                root_path[i] = '\0';
            }
            else
                break;
        }
        root_path[strlen(root_path)-1] = '\0';
        for(i=strlen(root_path)-1; i >= 0; i--){
            if(root_path[i] != '/'){
                root_path[i] = '\0';
            }
            else
                break;
        }
        //Set default paths
        sprintf(pids_file, "%spids.txt", root_path);
        sprintf(str2str_path, "%sRTKLIB-b34e/app/consapp/str2str/gcc/str2str", root_path);
        sprintf(rtkrcv_path, "%sRTKLIB-b34e/app/consapp/rtkrcv/gcc/rtkrcv", root_path);
        sprintf(rtkconf, "%sconf/rtk4pid.conf", root_path);
        sprintf(pppconf, "%sconf/ppp4pid_navcast.conf", root_path);
        sprintf(semor_conf_path, "%ssemor.conf", root_path);
    }
    FILE* f;
    if( access( semor_conf_path, F_OK ) == 0 ) {
        // read file
        f = fopen(semor_conf_path, "r");
        while(fgets(line, sizeof(line), f)){
            read_conf_line(line);
        }
    } else {
        // create file and use default parameters
        f = fopen(semor_conf_path, "w");

        fprintf(f, "$General\n");
        fprintf(f, "debug=%d   $0:disabled (get realtime data from rtkrcv and imu), 1:enabled (get data from files (in test folder))\n", debug);
        fprintf(f, "logs=%d   $0:disabled, 1:enabled\n", logs);
        fprintf(f, "str2str-path=%s\n", str2str_path);
        fprintf(f, "str2str-in=%s\n", str2str_in);
        fprintf(f, "str2str-out1-port=%s\n", str2str_out1_port);
        fprintf(f, "str2str-out2-port=%s\n", str2str_out2_port);
        fprintf(f, "rtkrcv-path=%s\n", rtkrcv_path);
        fprintf(f, "rtkconf=%s\n", rtkconf);
        fprintf(f, "pppconf=%s\n", pppconf);
        fprintf(f, "rtkrcv-port-rtk=%s\n", rtk_port_rtkrcv);
        fprintf(f, "rtkrcv-port-ppp=%s\n", ppp_port_rtkrcv);
        fprintf(f, "\n$IMU parameters\n");
        fprintf(f, "imu-init-epochs(sec)=%d\n", imu_init_epochs);
        fprintf(f, "imu-drift=%d\n", imu_drift);
        fprintf(f, "sample-rate(hz)=%d\n", sample_rate);
        fprintf(f, "init-bg-unc=%.15g\n", init_bg_unc);
        fprintf(f, "init-ba-unc=%.15g\n", init_ba_unc);
        fprintf(f, "psd-gyro=%.15g\n", psd_gyro);
        fprintf(f, "psd-acce=%.15g\n", psd_acce);
        fprintf(f, "psd-bg=%.15g\n", psd_bg);
        fprintf(f, "psd-ba=%.15g\n", psd_ba);
        fprintf(f, "$Initialization coordinates\n", psd_ba);
        fprintf(f, "init-x=0\n");
        fprintf(f, "init-y=0\n");
        fprintf(f, "init-z=0\n");
    }
    fclose(f);

    if(init_pos != 7){
        printf("Please, set the initialization coordinates in the configuration file.\n");
        close_semor(1);
    }

    //Initialize shared variables
    str2str_pid = rtkrcv1_pid = rtkrcv2_pid = -1;
    imu_ready = 0;

    if(!debug){
        //Set arguments for execv function
        char *const str2str_args[] = {str2str_path, "-in", str2str_in, "-out", str2str_out1_port, "-out", str2str_out2_port, NULL};
        char *const rtkrcv1_args[] = {rtkrcv_path, "-s", "-o", rtkconf, "-d", "/dev/null", NULL}; //"/dev/null" needed to remove rtkrcv terminal issues (it takes control
        //of standard input and output)
        char *const rtkrcv2_args[] = {rtkrcv_path, "-s", "-o", pppconf, "-d", "/dev/null", NULL};

        //Execute str2str
        if ((str2str_pid = fork()) == -1){
            perror("SEMOR: fork error: str2str");
            close_semor(1);
        }
        else if (str2str_pid == 0) {
            close_io();
            execv(str2str_args[0], str2str_args);
            printf("\nSEMOR: execv error (str2str)"); //if this line is executed, something is wrong
            close_semor(1);
        }

        //Execute first rtkrcv instance
        if ((rtkrcv1_pid = fork()) == -1){
            perror("SEMOR: fork error: rtkrcv(1)");
            close_semor(1);
        }
        else if (rtkrcv1_pid == 0) {
            close_io();
            execv(rtkrcv1_args[0], rtkrcv1_args);
            printf("\nSEMOR: execv error (rtkrcv(1))");
            close_semor(1);
        }
        //Execute second rtkrcv instance

        if ((rtkrcv2_pid = fork()) == -1){
            
            perror("SEMOR: fork error: rtkrcv(2)");
            close_semor(1);
        }
        else if (rtkrcv2_pid == 0) {
            close_io();
            execv(rtkrcv2_args[0], rtkrcv2_args); //da cambiare con bin/rtkrcv
            printf("\nSEMOR: execv error (rtkrcv(2))");
            close_semor(1);
        }

        //Create pids.txt file to manually kill previous executed processes if needed (in case someone forgots to soft-terminate SEMOR entering 'q')
        FILE *pids = fopen(pids_file, "w");
        fprintf(pids, "str2str: %d\nrtkrcv(1): %d\nrtkrcv(2): %d", str2str_pid, rtkrcv1_pid, rtkrcv2_pid);
        fclose(pids);

    }

    struct stat st = {0};


    sprintf(semor_conf_path, "%slogs", root_path);
    if (stat(semor_conf_path, &st) == -1) {
        mkdir(semor_conf_path, 0777);
    }

   

    printf("\n");
    printf("SEMOR: press q to stop\n");

    /*
    SEMOR core starts here
    */
    start_processing();



    return 0;
}