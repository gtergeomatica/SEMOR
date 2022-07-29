#include "semor.h"
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <fcntl.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <signal.h>
#include <string.h>
#include <math.h>
#include <sys/poll.h>
#include <time.h>
#include <sys/stat.h>


//File containing the output positions from SEMOR
#define FILE_PATH "output.txt"

//Input file of RTK and PPP for debugging (in future developments their definition can be added to the conf file)
#define GPS_FILE "test/gps.pos"
#define GALILEO_FILE "test/galileo.pos"

//Constants to improve readability
#define GPS 0
#define GALILEO 1
#define IMU 2

//GPS=0, GALILEO=1, IMU=2 --- array in which new (input) data from GNSS and IMU is written and ready to be processed
gnss_sol_t sol[3];

//Structure in which the "best" solution, so the output of SEMOR for each epoch, is stored (taken by SiConsulting algorithm using get_data() function)
gnss_sol_t best;

//Initial position for imu initialization
gnss_sol_t initial_pos;

//Log files for RTK, PPP and IMU
FILE* sol_file[3];

//output file (with path: FILE_PATH)
FILE *file;

//Declaration of pids variables of str2str and rtkrcv (2 instances)
pid_t str2str_pid, rtkrcv1_pid, rtkrcv2_pid;

int imu_ready; //0: the imu biases are not already calculated so SEMOR can't process any solutions right now, 1: SEMOR is processing solutions 
gnss_sol_t first_pos;

//Used in debug mode
//Set to 1 when one or both RTK and PPP solutions are ahead of the seconds variable (defined later in the code) -> no new data is read from the file/files until
//the seconds variable arrive to the same epochs as both RTK and PPP. This is needed for the synchronization of the solutions and so for the processing
int wait_read[2];

//Used in debug mode
//SEMOR set "sol[<index>].time.week = 0" to indicate that a solution will be ignored in the output production. We set sol[<index>].time.week = 0, for example, when
//we want to ignore a solution A that is waiting for one other solution B since B is behind A. When B reaches again A, A.time.week will be reset to last_week[<index of A>]
//and it won't be ignored for the next rounds
int last_week[2];

//Used in debug mode to simulate the time passing (since input files of RTK and PPP may contain discontinuities among epochs)
int seconds;

//how many times the imu solution has been used consecutively
//when n_imu reaches imu_drift, SEMOR terminates and needs to be re-initialized
int n_imu;

//used to execute init_imu(sol[IMU]) only one/the first time
int first_time = 1;

//path to de log directory
char log_dir[PATH_MAX/2];

const double PI = 3.1415926535898;
double radianToDegree(double r) {
  return r * (180 / PI);
}

//Conversion from ecef to geodetic coordinates
gnss_sol_t ecef2geo_(gnss_sol_t gnss){
  //gives std deviations in m with format E N U
	// Variables
	double x, y, z;
  double lat_P1, lng_P1, h_P1;
  double slat, slng, sh;


  double x_P1 = gnss.a;
  double y_P1 = gnss.b;
  double z_P1 = gnss.c;
  double sx = gnss.sda;
  double sy = gnss.sdb;
  double sz = gnss.sdc;

  // printf("%lf | ", gnss.sda);
  // printf("%lf | ", gnss.sdb);
  // printf("%lf ||", gnss.sdc);
/* Point n°1 */
  x = x_P1; y = y_P1; z = z_P1;
	// Semi Major Axis and Eccentricity
	const double a = 6378137; const double e = 0.08181919;
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
	lat_P1 = radianToDegree(phi_i);
	lng_P1 = radianToDegree(lambda);
	h_P1 = h;
  gnss.a = lat_P1;
  gnss.b = lng_P1;
  gnss.c = h_P1;

  double phi = phi_i;
  double lam = lambda;

  gnss.sda = fabs(-sin(lam)*sx + cos(lam)*sy);
  gnss.sdb = fabs(-cos(lam)*sin(phi)*sx - sin(lam)*sin(phi)*sy +cos(phi)*sz);
  gnss.sdc = fabs(cos(lam)*cos(phi)*sx + sin(lam)*cos(phi)*sy + sin(phi)*sz);
  //
  // printf("%lf | ", gnss.sda);
  // printf("%lf | ", gnss.sdb);
  // printf("%lf \n", gnss.sdc);

  return gnss;

}


LocData_t get_data(){ //SiConsulting
    LocData_t data;
    gnss_sol_t geo_best;

    int sec_today, sec_today_utc, h_utc, h, m, s;

    time_t rawtime, ti;
    struct tm info;

    geo_best = ecef2geo_(best);

    data.dLat = geo_best.a;
    data.dLon = geo_best.b;
    data.dHeigth = geo_best.c;

    data.dSdn = geo_best.sda;
    data.dSde = geo_best.sdb;
    data.dSdu = geo_best.sdc;
    data.ui8FixQual = geo_best.Q;



    sec_today = (best.time.sec-LEAP_SECONDS)%(24*3600); //7200 because GPS is shifted

    /*time( &rawtime );
    info = *localtime( &rawtime );
    h_utc = info.tm_hour; */

    h = sec_today/3600;
	
	m = (sec_today -(3600*h))/60;
	
	s = (sec_today -(3600*h)-(m*60));

    sprintf(data.ui8TS, "%02d:%02d:%02d", h, m, s);
    
    return data;
}



//Close SEMOR and all processes started by it (rtkrcv and str2str)
void close_semor(int status){
    printf("\nSEMOR terminated.\n");

    //Kill RTKLIB executables
    if(str2str_pid != -1 && kill(str2str_pid, SIGKILL) == -1){
        perror("SEMOR: Error killing str2str process");
    }
    if(rtkrcv1_pid != -1 && kill(rtkrcv1_pid, SIGKILL) == -1){
        perror("SEMOR: Error killing rtkrcv(1) process");
    }
    if(rtkrcv2_pid != -1 && kill(rtkrcv2_pid, SIGKILL) == -1){
        perror("SEMOR: Error killing rtkrcv(2) process");
    }

    //lock used in loose-gnss-imu/Loosely.cpp to read and write the imu_data matrix
    pthread_mutex_destroy(&lock);
    
    //close raw imu data log file
    close_ctocpp();

    //close log files
    if(logs){
        fclose(sol_file[GPS]);
        fclose(sol_file[GALILEO]);
        fclose(sol_file[IMU]);
    }

    //close output file
    fclose(file);
    exit(status);
}


gnss_sol_t str2gnss(char str[MAXSTR]){
    gnss_sol_t gnss;
    char *eptr;
    char copy[MAXSTR];

    strcpy(copy, str);

    gnss.time.week = atoi(strtok(copy, " "));
    gnss.time.sec = atoi(strtok(NULL, " "));
    gnss.a = strtod(strtok(NULL, " "), &eptr);
    gnss.b = strtod(strtok(NULL, " "), &eptr);
    gnss.c = strtod(strtok(NULL, " "), &eptr);
    gnss.Q = atoi(strtok(NULL, " "));
    gnss.ns = atoi(strtok(NULL, " "));
    gnss.sda = strtod(strtok(NULL, " "), &eptr);
    gnss.sdb = strtod(strtok(NULL, " "), &eptr);
    gnss.sdc = strtod(strtok(NULL, " "), &eptr);
    gnss.sdab = strtod(strtok(NULL, " "), &eptr);
    gnss.sdbc = strtod(strtok(NULL, " "), &eptr);
    gnss.sdca = strtod(strtok(NULL, " "), &eptr);
    gnss.age = strtof(strtok(NULL, " "), &eptr);
    gnss.ratio = strtof(strtok(NULL, " "), &eptr);
    gnss.va = strtod(strtok(NULL, " "), &eptr);
    gnss.vb = strtod(strtok(NULL, " "), &eptr);
    gnss.vc = strtod(strtok(NULL, " "), &eptr);

    return gnss;
}

void gnss2str(char* str, gnss_sol_t gnss){
    sprintf(str, "sec: %d, pos: %lf %lf %lf, std: %lf %lf %lf", gnss.time.sec, gnss.a, gnss.b, gnss.c, gnss.sda, gnss.sdb, gnss.sdc);
}

void gnsscopy(gnss_sol_t *dest, gnss_sol_t src){
    (*dest).time.week = src.time.week;
    (*dest).time.sec = src.time.sec;
    (*dest).a = src.a;
    (*dest).b = src.b;
    (*dest).c = src.c;
    (*dest).Q = src.Q;
    (*dest).ns = src.ns;
    (*dest).sda = src.sda;
    (*dest).sdb = src.sdb;
    (*dest).sdc = src.sdc;
    (*dest).sdab = src.sdab;
    (*dest).sdbc = src.sdbc;
    (*dest).sdca = src.sdca;
    (*dest).age = src.age;
    (*dest).ratio = src.ratio;
    (*dest).va = src.va;
    (*dest).vb = src.vb;
    (*dest).vc = src.vc;
}

//print solution to output file
void output(gnss_sol_t sol){
    char sol1[MAXSTR];
    gnss2str(sol1, sol);
    fprintf(file, "%s\n", (sol.time.week == 0) ? "no data" : sol1);
    fflush(file);
}

//print solution to log file
void print_solution(int sol_index){
    char sol1[MAXSTR];
    gnss2str(sol1, sol[sol_index]);
    fprintf(sol_file[sol_index], "%s\n", (sol[sol_index].time.week == 0) ? "no data" : sol1);
    fflush(sol_file[sol_index]);
}

//check if 2 positions are similar (if the range of one position plus or minus its standard deviation overlaps with the range of the other)
int similar_pos(gnss_sol_t p1, gnss_sol_t p2){

    //Check a
    if(p1.a < p2.a){
        if(p1.a+3*p1.sda < p2.a-3*p2.sda){
            return 0;
        }
    }
    else{
        if(p2.a+3*p2.sda < p1.a-3*p1.sda){
            return 0;
        }
    }

    //Check b
    if(p1.b < p2.b){
        if(p1.b+3*p1.sdb < p2.b-3*p2.sdb){
            return 0;
        }
    }
    else{
        if(p2.b+3*p2.sdb < p1.b-3*p1.sdb){
            return 0;
        }
    }

    //Check c
    if(p1.c < p2.c){
        if(p1.c+3*p1.sdc < p2.c-3*p2.sdc){
            return 0;
        }
    }
    else{
        if(p2.c+3*p2.sdc < p1.c-3*p1.sdc){
            return 0;
        }
    }

    return 1;
}

//Compute the mean of 2 positions (used when only 2 positions, out of the 3, are similar)
void gnss_avg_2(gnss_sol_t sol1, gnss_sol_t sol2){
    //Position
    best.a = (sol1.a + sol2.a)/2;
    best.b = (sol1.b + sol2.b)/2;
    best.c = (sol1.c + sol2.c)/2;

    //Velocity
    best.va = (sol1.va + sol2.va)/2;
    best.vb = (sol1.vb + sol2.vb)/2;
    best.vc = (sol1.vc + sol2.vc)/2;

    //Standard deviation
    best.sda = sqrt(pow(sol1.sda/2, 2) + pow(sol2.sda/2, 2));
    best.sdb = sqrt(pow(sol1.sdb/2, 2) + pow(sol2.sdb/2, 2));
    best.sdc = sqrt(pow(sol1.sdc/2, 2) + pow(sol2.sdc/2, 2));

    best.sdab = sqrt(pow(sol1.sdab/2, 2) + pow(sol2.sdab/2, 2));
    best.sdbc = sqrt(pow(sol1.sdbc/2, 2) + pow(sol2.sdbc/2, 2));
    best.sdca = sqrt(pow(sol1.sdca/2, 2) + pow(sol2.sdca/2, 2));

    best.time.week = sol1.time.week;
    best.time.sec = sol1.time.sec;
}

//Compute the mean of 3 positions (when all the positions are similar)
void gnss_avg_3(gnss_sol_t sol1, gnss_sol_t sol2, gnss_sol_t sol3){
    //Position
    best.a = (sol1.a + sol2.a + sol3.a)/3;
    best.b = (sol1.b + sol2.b + sol3.b)/3;
    best.c = (sol1.c + sol2.c + sol3.c)/3;

    //Velocity
    best.va = (sol1.va + sol2.va + sol3.va)/3;
    best.vb = (sol1.vb + sol2.vb + sol3.vb)/3;
    best.vc = (sol1.vc + sol2.vc + sol3.vc)/3;

    //Standard deviation
    best.sda = sqrt(pow(sol1.sda/3, 2) + pow(sol2.sda/3, 2) + pow(sol3.sda/3, 2));
    best.sdb = sqrt(pow(sol1.sdb/3, 2) + pow(sol2.sdb/3, 2) + pow(sol3.sdb/3, 2));
    best.sdc = sqrt(pow(sol1.sdc/3, 2) + pow(sol2.sdc/3, 2) + pow(sol3.sdc/3, 2));

    best.sdab = sqrt(pow(sol1.sdab/3, 2) + pow(sol2.sdab/3, 2) + pow(sol3.sdab/3, 2));
    best.sdbc = sqrt(pow(sol1.sdbc/3, 2) + pow(sol2.sdbc/3, 2) + pow(sol3.sdbc/3, 2));
    best.sdca = sqrt(pow(sol1.sdca/3, 2) + pow(sol2.sdca/3, 2) + pow(sol3.sdca/3, 2));

    best.time.week = sol1.time.week;
    best.time.sec = sol1.time.sec;
}

int get_best_sol2(int sol1_idx, int sol2_idx){ //if 2 gnss solutions available - 0: no best found, 1: best found
    if(similar_pos(sol[sol1_idx], sol[sol2_idx])){
        gnss_avg_2(sol[sol1_idx], sol[sol2_idx]);
        return 1;
    }
    return 0;
}

int  get_best_sol_3(){ //if 3 solutions available - 0: no best found, 1: best found (if no best found -> use IMU)

    if(similar_pos(sol[GPS], sol[IMU]) && similar_pos(sol[IMU], sol[GALILEO]) && similar_pos(sol[GPS], sol[GALILEO])){
        gnss_avg_3(sol[GPS], sol[GALILEO], sol[IMU]);
        return 1;
    }

    if(similar_pos(sol[GPS], sol[GALILEO])){
        gnss_avg_2(sol[GPS], sol[GALILEO]);
        return 1;
    }

    if(similar_pos(sol[GPS], sol[IMU])){
        gnss_avg_2(sol[GPS], sol[IMU]);
        return 1;
    }

    if(similar_pos(sol[GALILEO], sol[IMU])){
        gnss_avg_2(sol[GALILEO], sol[IMU]);
        return 1;
    }
    return 0;
}

//core processing of semor here (comparison of the available solutions for each epoch and the calcultation of the output)
void process_solutions(int chk_sols){
    int i; //counter variable
    int is_best_found; //it tells if a best is found: the only time SEMOR can't find the best solution is when we have 2 or 3 solutions available but none of them are similar

    if(debug){
        for(i = 0; i < 2; i++){ //for each gnss solution check if its epoch is higher than the "seconds" variable, if so wait next iterations 
        //before displaying the solution until epochs are equal
            if(sol[i].time.sec > seconds){
                if(!wait_read[i])
                    last_week[i] = sol[i].time.week;
                wait_read[i] = 1;
                sol[i].time.week = 0;

                chk_sols -= pow(2, i); //this position belongs to a next epoch, so it won't be utilized
            }
            else{
                wait_read[i] = 0;
                if(last_week[i] != 0){
                    sol[i].time.week = last_week[i];
                    last_week[i] = 0;
                }
            }
        }
    }

    //while imu is not ready, keep calling imu_sol and return (in this condition imu_sol won't calculate the next imu position but will initialize the imu)
    if(!imu_ready){
        sol[IMU].time.sec++;
        imu_sol(&sol[IMU]); //this takes 1 second
        return;
    }

    //sol[IMU].time.week == 0, no solution will be available (actually I think it's impossibile for this condition to be true ever)
    if(sol[IMU].time.week != 0) 
        chk_sols |= 4;

    //Get best solution between GPS, GALILEO, IMU
    switch(chk_sols){
        case 0: //no solutions, return
            return;
        case 1: //only GPS
            gnsscopy(&best, sol[GPS]);
            is_best_found = 1;
            break;
        case 2: // only GALILEO
            gnsscopy(&best, sol[GALILEO]);
            is_best_found = 1;
            break;
        case 3: //GPS and GALILEO
            is_best_found = get_best_sol2(GPS, GALILEO);
            best.time.week = sol[GPS].time.week;
            best.time.sec = sol[GPS].time.sec;
            break;
        case 4: //only IMU
            gnsscopy(&best, sol[IMU]);
            is_best_found = 1;
            break;
        case 5: //GPS and IMU
            is_best_found = get_best_sol2(GPS, IMU);
            break;
        case 6: //GALILEO and IMU
            is_best_found = get_best_sol2(GALILEO, IMU);
            break;
        case 7: //all solutions available
            is_best_found = get_best_sol_3();
            break;
    }

    
    //print logs
    if(logs){
        print_solution(GPS);
        print_solution(GALILEO);
        print_solution(IMU);
    }

    //if is_best_found == 0 (IMU solution is used)
    if(!is_best_found){
        gnsscopy(&best, sol[IMU]);

        if(n_imu == imu_drift){
            //REINITIALIZE SEMOR
            printf("Re-initialize SEMOR\n");
            close_semor(1);
        }
        else
            n_imu++;
    }
    else{
        //reset n_imu since at least two solutions are similar and the output solution is not drifting
        n_imu = 0;
    }

    //Here we have the output solutions stored in best

    //Print the solution in the output file
    output(best);

    //Flag solutions as already used
    sol[GPS].time.week = 0;
    sol[GALILEO].time.week = 0;

    //Post comparison and output
    //So let's generate the next imu position propagating "best" with the raw acceleration and rotation data until the next epoch
    gnsscopy(&sol[IMU], best);
    sol[IMU].time.week = 0;
    sol[IMU].time.sec += 1; //Get imu position of the next second
    imu_sol(&sol[IMU]); //this runs at least until raw IMU timestamp < next second

}

//Check if the user entered 'q' to terminate SEMOR
void check_termination(){
    char cmd;
    cmd = getchar();
    if(cmd == 'q' || cmd == 'Q'){
        close_semor(0);
    }
}

//Establish connection to the a rtkrcv sockets (executed twice in SEMOR, one for RTK and one for PPP)
void setup_tcp_socket(int* fd, char port[6]){
    int status;
    struct addrinfo hints, *res;

    memset(&hints, 0, sizeof hints);
    hints.ai_socktype = SOCK_STREAM;
    hints.ai_flags = AI_PASSIVE;
    hints.ai_family = AF_INET; 

    
    if ((status = getaddrinfo(NULL, port, &hints, &res)) != 0) {
        fprintf(stderr, "SEMOR: getaddrinfo: %s\n", gai_strerror(status));
        close_semor(1);
    }
    *fd = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
    
    if(*fd == -1){
        perror("SEMOR socket()");
        close_semor(1);
    }

    int flags = fcntl(*fd, F_GETFL, 0);

    if(fcntl(*fd, F_SETFL, flags | O_NONBLOCK) == -1){
        perror("SEMOR: can't set non blocking socket()");
        close_semor(1);
    }

    do{
        connect(*fd, res->ai_addr, res->ai_addrlen);
        if(errno != EISCONN){
            continue; //Wait for rtkrcv to start and to send data
        }
        break; //Stop while
    }while(1);

    free(res);
}

//read a gnss solution from a socket
int read_gnss(int fd, int* offset, char* buf, int buf_length){
    int nbytes = 0;

    if(debug){
        do{
            nbytes = read(fd, buf+*offset, 1);
            if(nbytes == -1){
                perror("SEMOR read file");
                close_semor(1);
            }else if(nbytes == 0){ // Nothing to read
                continue;
            }
            (*offset)++;
        }while(buf[(*offset)-1] != '\r' && buf[(*offset)-1] != '\n');
    }else{
        if ((nbytes = recv(fd, buf + *offset, buf_length - *offset, 0)) < 0) { //(sizeof buf) /*strlen(buf)*/-*offset
            if(errno != EAGAIN && errno != EWOULDBLOCK && errno != ECONNREFUSED){
                perror("SEMOR: read");
                close_semor(1);
            }
        }
    *offset = *offset+nbytes;
    }

    return nbytes;
}

//This function contains the main loop in which: 
//1)SEMOR reads GNSS data from the two instances of rtkrcv
//2)it reads IMU data
//3)it processes the best solution (that can be one of the 3 solutions in input or the mean of at least 2 of them)
//4)it calculates the imu position for the next epoch
void handle_connection(){
    struct addrinfo hints, *res;
    int status;
    char buf[2][MAXSTR];
    int offset[2] = {0, 0};
    int nbytes;
    int socketfd[2];
    char port[6];
    int ret;
    int i;
    struct pollfd fds[3];
    struct pollfd fds1[2][1];
    int timeout_msecs = 1000;
    int check_sols = 0; //3 if both rtk and ppp are read, 1 if only rtk read, 2 if only ppp and 0 if none
    int galileo_ready = 0;
    int new_gnss_data = 0;
    int first_input = 0;

    printf("SEMOR: Start initialization. Please wait...\n");
    
    //Initialize input file descriptors
    if(debug){
        socketfd[GPS] = open(GPS_FILE, O_RDONLY);
        socketfd[GALILEO] = open(GALILEO_FILE, O_RDONLY);
    }
    else{
        setup_tcp_socket(&socketfd[GPS], rtk_port_rtkrcv);
        setup_tcp_socket(&socketfd[GALILEO], ppp_port_rtkrcv);
    }

    //Initialize the structures needed by poll()
    fds[0].fd = socketfd[GPS];
    fds[1].fd = socketfd[GALILEO];
    fds[2].fd = STDIN_FILENO;
    fds[0].events = fds[1].events = fds[2].events = POLLIN;

    fds1[GPS][0].fd = socketfd[GPS];
    fds1[GALILEO][0].fd = socketfd[GALILEO];
    fds1[GPS][0].events = fds1[GALILEO][0].events = POLLIN;

    //Initialize to 0 some variables
    wait_read[0] = wait_read[1] = 0;
    n_imu = 0;

    sol[IMU].time.week = 0;


    while(1){
        check_sols = 0;
        ret = poll(fds, 3, timeout_msecs); //wait for events on the 3 fds
        if (ret == -1){
            perror("SEMOR: poll");
            close_semor(1);
        }
        if(fds[2].revents & POLLIN){
            check_termination(); //Check if user requested the termination of the process
        }
        for(i = 0; i < 2; i++){ //Get GPS and GALILEO solutions
            offset[i] = 0;

            usleep(150000); //gives time to the solutions to be read (if one of them is late)

            if(debug && wait_read[i]){ //Don't read solution i (0:GPS, 1:GALILEO) if the solution in the previous iterations has a higher epoch than "seconds" variable
                continue;
            }
            nbytes = read_gnss(socketfd[i], &offset[i], buf[i], sizeof buf[i]);

            if(strstr(buf[i], "lat") || strstr(buf[i], "latitude") || strstr(buf[i], "ecef") || strlen(buf[i]) < 5){ //Check if the input string contains gnss data or if it is an empty line or header
                offset[i] = 0;

                //Clean up the buffer
                for(int j=0; j<MAXSTR;j++){
                    buf[i][j] = 0;
                }
                continue;
            }

            if(offset[i] != 0 && (buf[i][offset[i]-1] == '\r' || buf[i][offset[i]-1] == '\n')){ //If incoming data is a full gnss measurement string
                buf[i][offset[i]-1] = '\0'; //Replace new line with end of string
                check_sols |= i+1;
                first_input |= i+1;
                sol[i] = str2gnss(buf[i]);//Parse string to gnss_sol_t structure

                //Clean up the buffer
                for(int j=0; j<MAXSTR;j++){
                    buf[i][j] = 0;
                }
            }
        }

        if(first_input < 1){ //at least one GNSS solution available
            continue;
        }

        if(debug){
            if(seconds == 0){
                seconds = (sol[GPS].time.sec <= sol[GALILEO].time.sec) ? sol[GPS].time.sec : sol[GALILEO].time.sec; //Set current second to the minimum of the epochs of the 2 gnss solution
            }
        }

        if(first_time){
            printf("SEMOR: connection established with rtkrcv\nwait %d seconds for the IMU to initialize...\n", imu_init_epochs);
            //Initialize imu epoch

            sol[IMU].time.week = sol[GPS].time.week != 0 ? sol[GPS].time.week : sol[GALILEO].time.week;
            sol[IMU].time.sec = sol[GPS].time.sec != 0 ? sol[GPS].time.sec : sol[GALILEO].time.sec;

            if(debug){
                sol[IMU].time.sec = seconds;
            }
            init_imu(sol[IMU]);
            first_time = 0;
        }

        process_solutions(check_sols); //Get best solution, output it and use it to calculate next imu position

        if(debug)
            seconds++; //Update current second
    }


}

//Initialize structure for the initial position (initial input of the IMU), open log and output files and call handle_connection()
void start_processing(void){
    char path[3][PATH_MAX];

    //Initialize structure of initial_pos
    sol[IMU].a = init_x;
    sol[IMU].b = init_y;
    sol[IMU].c = init_z;

    sol[IMU].sda = 0;
    sol[IMU].sdb = 0;
    sol[IMU].sdc = 0;

    sol[IMU].va = 0;
    sol[IMU].vb = 0;
    sol[IMU].vc = 0;

    //Create log folder and files for this very execution of SEMOR
    time_t rawtime;
    struct tm info;
    time( &rawtime );
    info = *localtime( &rawtime );
    char str_time[13];

    sprintf(str_time, "%d_%02d_%02d_%02d_%02d", (info.tm_year+1900), (info.tm_mon+1), info.tm_mday, info.tm_hour, info.tm_min);

    struct stat st = {0};
    log_dir[PATH_MAX/2];

    if(logs){
        sprintf(log_dir, "%slogs/logs_%s", root_path, str_time);
        if (stat(log_dir, &st) == -1) {
            mkdir(log_dir, 0777);
        }
        sprintf(path[0], "%s/gps_%s.log", log_dir, str_time);
        sprintf(path[1], "%s/galileo_%s.log", log_dir, str_time);
        sprintf(path[2], "%s/imu_%s.log", log_dir, str_time);

        sol_file[GPS] = fopen(path[0], "w");

        sol_file[GALILEO] = fopen(path[1], "w");

        sol_file[IMU] = fopen(path[2], "w");
    }

    char imu_raw_log[PATH_MAX/2];

    sprintf(imu_raw_log, "%s/imu_raw.log", log_dir, str_time);

    FILE* fimu_raw = fopen(imu_raw_log, "w");
    fclose(fimu_raw);

    //open output file
    file = fopen(FILE_PATH, "w");
    if(file == NULL){
        perror("SEMOR fopen()");
    }
    fflush(file);

    handle_connection();
}

