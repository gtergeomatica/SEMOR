#include <stdio.h>
#include <signal.h>
#include <string.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "semor.h"

#define MAX_LINE 256

//Set default configuration
int relative = 0;
int logs=1;
int debug=1;
double init_bg_unc = 2.42406840554768e-05;
double init_ba_unc = 0.048901633857000000;
double psd_gyro  = 3.38802348178723e-09;
double psd_acce      =     2.60420170553977e-06;     // accelerator noise PSD (m^2/s^3)  
double psd_bg        =     2.61160339323310e-14;     // gyro bias random walk PSD (rad^2/s^3)
double psd_ba        =     1.66067346797506e-09;
int sample_rate = 104; //in hz
int imu_init_epochs =  20;//in seconds
int imu_drift = 60;

char rtk_port_rtkrcv[6] = "8090";
char ppp_port_rtkrcv[6] = "8091";

double init_x;
double init_y;
double init_z;
int init_pos = 0; //7 if init_x, init_y and init_z are assigned with a value, semor stops otherwise

char str2str_path[PATH_MAX];
char rtkrcv_path[PATH_MAX];
char rtkconf[PATH_MAX];
char pppconf[PATH_MAX];

char root_path[PATH_MAX-200];

char str2str_in[30] = "tcpcli://192.168.2.91:8081";
char str2str_out1_port[16] = "tcpsvr://:8085";
char str2str_out2_port[16] = "tcpsvr://:8086";


void close_io(void){
    close(STDIN_FILENO);
    close(STDOUT_FILENO);
    close(STDERR_FILENO);
}

void read_conf_line(char line[MAX_LINE]){
    int i;
    char *token;
    char *eptr;
    char delim[8] = "$= \t";
    token = strtok(line, delim);

    if(strstr(token, "rtkrcv-port-rtk")){
        strcpy(rtk_port_rtkrcv, strtok(NULL, delim));
        if(rtk_port_rtkrcv[strlen(rtk_port_rtkrcv)-1] == '\n')
            rtk_port_rtkrcv[strlen(rtk_port_rtkrcv)-1] = '\0';
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
    //pid_t str2str_pid, rtkrcv1_pid, rtkrcv2_pid;
    int isset_first_pos;
    int imu_ready;

    //Setup from configuration
    int i;
    char line[MAX_LINE];
	char cwd[PATH_MAX-400];
    char semor_conf_path[PATH_MAX];
    char pids_file[PATH_MAX];

    char path[60];
    /*if(argc == 2 && argv[1][0] == '-' && argv[1][1] == 'r'){
        relative = 1;
        sprintf(str2str_path, "RTKLIB-b34e/app/consapp/str2str/gcc/str2str");
        sprintf(rtkrcv_path, "RTKLIB-b34e/app/consapp/rtkrcv/gcc/rtkrcv");
        sprintf(rtkconf, "conf/rtk4pid.conf");
        sprintf(pppconf, "conf/ppp4pid_navcast.conf");
        sprintf(semor_conf_path, "bin/semor.conf");
        sprintf(pids_file, "pids.txt");
    }*/
    {
        sprintf(path, "/proc/%d/exe", getpid());
        if(readlink(path, root_path, PATH_MAX) == -1){
            perror("SEMOR: readlink()");
            printf("\n");
            printf("Possible workaround:\n");
            printf("Go exactly in the SEMOR root folder and run semor with the argument '-r' (use relative paths)\n");
            printf("Example:\n>cd SEMOR\n>bin/semor -r");
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
    isset_first_pos = 0;
    imu_ready = 0;
    if(!debug){
        char *const str2str_args[] = {str2str_path, "-in", str2str_in, "-out", str2str_out1_port, "-out", str2str_out2_port, NULL};
        char *const rtkrcv1_args[] = {rtkrcv_path, "-s", "-o", rtkconf, "-d", "/dev/null", NULL}; //rtk4pid.conf
        char *const rtkrcv2_args[] = {rtkrcv_path, "-s", "-o", pppconf, "-d", "/dev/null", NULL}; //ppp4pid_navcast.conf

        //Execute str2str
        if ((str2str_pid = fork()) == -1){
            perror("SEMOR: fork error: str2str");
            close_semor(1);
        }
        else if (str2str_pid == 0) {
            close_io();
            execv(str2str_args[0], str2str_args);
            printf("\nSEMOR: execv error (str2str)");
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

        //Create pids.txt file to manually kill previous executed processes if needed
        FILE *pids = fopen(pids_file, "w");
        fprintf(pids, "str2str: %d\nrtkrcv(1): %d\nrtkrcv(2): %d", str2str_pid, rtkrcv1_pid, rtkrcv2_pid);
        fclose(pids);

    }

    struct stat st = {0};


    sprintf(semor_conf_path, "%slogs", root_path);
    if (stat(semor_conf_path, &st) == -1) {
        mkdir(semor_conf_path, 0777);
    }

    /*
    Semor starts here
    */
   

    printf("\n");
    printf("SEMOR: press q to stop\n");
    start_processing();



    return 0;
}