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


#define DEBUG_FILE 1 //If semor takes data from file

#define FILE_PATH "output.txt"

#define GPS_FILE "test/gps.pos"
#define GALILEO_FILE "test/galileo.pos"

#define GPS 0
#define GALILEO 1
#define IMU 2

#define THRESHOLD 20000000

//Solutions

//char rtk_port_rtkrcv[6];
//char ppp_port_rtkrcv[6];

pthread_t imu_thread_id;

gnss_sol_t sol[2]; /* GPS=0, GALILEO=1 */
gnss_sol_t best;
//gnss_sol_t best_old;

gnss_sol_t initial_pos; //initial position for imu initialization

FILE* sol_file[2];

FILE *file;

pid_t str2str_pid, rtkrcv1_pid, rtkrcv2_pid;

int isset_first_pos;
gnss_sol_t first_pos;

int wait_read[2];
int last_week[2];
int seconds;

int first_time = 1;

char log_dir[PATH_MAX/2];

const double PI = 3.1415926535; //3.14159;


double radianToDegree(double r) {
  return r * (180 / PI);
}
/* HUGO */
/*gnss_sol_t ecef2geo_(gnss_sol_t gnss){
  //gives std deviations in ° for lat and lng, and in m for h
	// Variables
	double x, y, z;
  double lat_P1, lng_P1, h_P1;
  double lat_P2, lng_P2, h_P2;
  double slat, slng, sh;


  double x_P1 = gnss.a;
  double y_P1 = gnss.b;
  double z_P1 = gnss.c;
  double sx = gnss.sda;
  double sy = gnss.sdb;
  double sz = gnss.sdc;

  // printf("%lf | ", gnss.a);
  // printf("%lf | ", gnss.b);
  // printf("%lf | ", gnss.c);
  printf("%lf | ", gnss.sda);
  printf("%lf | ", gnss.sdb);
  printf("%lf || ", gnss.sdc);

  double x_P2 = x_P1 + sx;
  double y_P2 = y_P1 + sy;
  double z_P2 = z_P1 + sz;
// Point n°1
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

// Point n°2
  x = x_P2; y = y_P2; z = z_P2;
	// Compute Longitude
    lambda = atan2(y, x);
	// Physical radius of the point
    r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	// Radius in the x-y plane
	p = sqrt(pow(x, 2) + pow(y, 2));
	// GEOcentric latitude (Initial Approx)
	phi_o = atan2(p, z); phi_i = phi_o;
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
	lat_P2 = radianToDegree(phi_i);
	lng_P2 = radianToDegree(lambda);
	h_P2 = h;

  slat = fabs(lat_P2 - lat_P1);
  slng = fabs(lng_P2 - lng_P1);
  sh = fabs(h_P2 - h_P1);

  // printf("%lf | ", lat_P1);
  // printf("%lf | ", lng_P1);
  // printf("%lf | ", h_P1);
  printf("%lf | ", slat);
  printf("%lf | ", slng);
  printf("%lf \n", sh);
  gnss.sda = slat;
  gnss.sdb = slng;
  gnss.sdc = sh;

  return gnss;

}*/

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

/* HUGO END */



/*gnss_sol_t ecef2geo_(gnss_sol_t gnss){

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
	gnss.a = radianToDegree(phi_i);
	gnss.b = radianToDegree(lambda);
	gnss.c = h;
	return gnss;
}*/

gnss_sol_t ecef2geo(gnss_sol_t gnss){
    //2
    // WGS84 constants
    double a = 6378137.0;
    double f = 1.0 / 298.257223563;
    // derived constants
    double b = a - f*a;
    double e = sqrt(pow(a,2)-pow(b,2))/a;
    double clambda = atan2(gnss.b,gnss.a);
    double p = sqrt(pow(gnss.a,2)+pow(gnss.b,2));
    double h_old = 0.0;
    // first guess with h=0 meters
    double theta = atan2(gnss.c,p*(1.0-pow(e,2.0)));
    double cs = cos(theta);
    double sn = sin(theta);
    double N = pow(a,2)/sqrt(pow(a*cs,2)+pow(b*sn,2.0));
    double h = p/cs - N;
    while (fabs(h-h_old) > 1.0e-6){
        h_old = h;
        theta = atan2(gnss.c,p*(1.0-pow(e,2.0)*N/(N+h)));
        cs = cos(theta);
        sn = sin(theta);
        N = pow(a,2.0)/sqrt(pow(a*cs,2.0)+pow(b*sn,2.0));
        h = p/cs - N;
    }
    //llh = {'lon':clambda, 'lat':theta, 'height': h}
    gnss.a = radianToDegree(theta);
    gnss.b = radianToDegree(clambda);
    gnss.c = h;
    return gnss;
}


LocData_t get_data(){ //SiConsulting
    LocData_t data;
    gnss_sol_t geo_best;

    int sec_today, sec_today_utc, h_utc, h, m, s;

    time_t rawtime, ti;
    struct tm info;

    geo_best = best;

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
    if(str2str_pid != -1 && kill(str2str_pid, SIGKILL) == -1){
        perror("SEMOR: Error killing str2str process");
    }
    if(rtkrcv1_pid != -1 && kill(rtkrcv1_pid, SIGKILL) == -1){
        perror("SEMOR: Error killing rtkrcv(1) process");
    }
    if(rtkrcv2_pid != -1 && kill(rtkrcv2_pid, SIGKILL) == -1){
        perror("SEMOR: Error killing rtkrcv(2) process");
    }
    if(logs){
        fclose(sol_file[GPS]);
        fclose(sol_file[GALILEO]);
    }
    fclose(file);
    exit(status);
}

/*LocData_t get_data(){

}*/

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
    /*sprintf(str, "%d %d %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %f %f %lf %lf %lf", gnss.time.week,
    gnss.time.sec, gnss.a, gnss.b, gnss.c, gnss.Q, gnss.ns,
    gnss.sda, gnss.sdb, gnss.sdc, gnss.sdab, gnss.sdbc, gnss.sdca,
    gnss.age, gnss.ratio, gnss.va, gnss.vb, gnss.vc);*/
    sprintf(str, "sec: %d, pos: %.9f %.9f %.9f, std: %.9f %.9f %.9f", gnss.time.sec, gnss.a, gnss.b, gnss.c, gnss.sda, gnss.sdb, gnss.sdc);
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

void output(gnss_sol_t sol, int thrust){
    char sol1[MAXSTR];
    gnss2str(sol1, sol);
    if(sol.time.week == 0){
        fprintf(file, "no_data\n", (sol.time.week == 0) ? "no data" : sol1);
    }else{
        fprintf(file, "%s, thrust:%d\n", sol1, thrust); //if thrust > 0  OK, if thrust = 0 NOT OK
    }
    fflush(file);
}

void print_solution(int sol_index){
    char sol1[MAXSTR];
    gnss2str(sol1, sol[sol_index]);
    fprintf(sol_file[sol_index], "%s\n", (sol[sol_index].time.week == 0) ? "no data" : sol1);
    fflush(sol_file[sol_index]);
}

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
}


int get_best_sol2(int sol1_idx, int sol2_idx){ //if 2 gnss solutions available - 0: no best found, 1: best found
    if(similar_pos(sol[sol1_idx], sol[sol2_idx])){
        gnss_avg_2(sol[sol1_idx], sol[sol2_idx]);
        return 1;
    }
    return 0;
}


void process_solutions(int chk_sols){
    int i;
    int is_best_found = 0;
    //gnss_sol_t tmp;

    if(debug){
        for(i = 0; i < 2; i++){ //for each gnss solution check if its epoch is higher than the "seconds" variable, if so wait next iterations
        //before displaying the solution till epochs are equal
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

    switch(chk_sols){ //Get best solution between GPS, GALILEO, IMU
        case 0: //no solutions for this epoch
            /*if(best_old.time.week != 0){ //if i have 2 previous solutions
                gnsscopy(&tmp, best);
                best.a = best.a + best.a-best_old.a;
                best.b = best.b + best.b-best_old.b;
                best.c = best.c + best.c-best_old.c;
                best.time.sec++;
                if(best.time.sec == 604801){
                    best.time.week++;
                    best.time.sec = 1;
                }
                //Best_old is previous best
                best_old = tmp;
                gnsscopy(&best_old, tmp);
            }else{
                return; //no previous solutions, return
            }*/

            //If no solutions I use the previous best solution but with a mark saying we can't thrust it
            if(best.time.week != 0)
                is_best_found = 1;
            //Do not convert best to geodetic since it is already converted from previous epoch
            printf("No solutions for this epoch, previous best used\n");
            break;
        case 1: //only GPS
            printf("Only RTK for epoch: %d\n", sol[GPS].time.sec);
            gnsscopy(&best, sol[GPS]);
            is_best_found = 1;
            if(coord_type){
                gnsscopy(&sol[GPS], ecef2geo_(sol[GPS]));
                gnsscopy(&best, sol[GPS]);
            }
            break;
        case 2: // only GALILEO
        printf("Only RTK for epoch: %d\n", sol[GALILEO].time.sec);
            gnsscopy(&best, sol[GALILEO]);
            is_best_found = 1;
            if(coord_type){
                gnsscopy(&sol[GALILEO], ecef2geo_(sol[GALILEO]));
                gnsscopy(&best, sol[GALILEO]);
            }
            break;
        case 3: //GPS and GALILEO
            is_best_found = get_best_sol2(GPS, GALILEO);
            if(is_best_found){
                printf("Both RTK and PPP for epoch: %d, similar\n", sol[GPS].time.sec);
            }
            else{
                printf("Both RTK and PPP for epoch: %d, NOT similar\n", sol[GPS].time.sec);
            }
            best.time.week = sol[GPS].time.week;
            best.time.sec = sol[GPS].time.sec;

            if(coord_type){
                if(sol[GPS].time.week != 0) gnsscopy(&sol[GPS], ecef2geo_(sol[GPS]));
                if(sol[GALILEO].time.week != 0) gnsscopy(&sol[GALILEO], ecef2geo_(sol[GALILEO]));
                if(is_best_found){ //avoid to convert multiple times
                    gnsscopy(&best, ecef2geo_(best));
                }
            }
            break;
    }

    if(logs){
        print_solution(GPS);
        print_solution(GALILEO);
    }

    //Here we have the solution
    if(is_best_found) //if best is found, that is used
        output(best, chk_sols);
    else if(sol[GPS].time.week != 0){ //else RTK is used if present since it's more accurate
        gnsscopy(&best, sol[GPS]);
        output(best, chk_sols);

    }
    else if(sol[GALILEO].time.week != 0){ //else PPP is used if present
        gnsscopy(&best, sol[GALILEO]);
        output(best, chk_sols);

    }

    /* if(is_best_found && chk_sols != 1){
        output(best, chk_sols);
    }else if(sol[GPS].time.week != 0){
        //DA CAPIRE MEGLIO COSA FARE QUI
        gnsscopy(&best, sol[GPS]);
        output(best, chk_sols);

    }else{
        best.a = -1;
        fprintf(file, "no-data\n");
        fflush(file);
    } */

    //Flag solutions as already used
    sol[GPS].time.week = 0;
    sol[GALILEO].time.week = 0;
    best.time.week = 0;
    best.time.sec++;
    if(best.time.sec == 604801){
        best.time.week++;
        best.time.sec = 1;
    }


}
void check_termination(){
    char cmd;
    cmd = getchar();
    if(cmd == 'q' || cmd == 'Q'){
        close_semor(0);
    }
}

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
        if(errno != EISCONN){ //(errno == ECONNREFUSED || errno == EAGAIN || errno == EINPROGRESS || errno == EALREADY)
            continue; //Wait for rtkrcv to start and to send data
        }
        break; //Stop while
    }while(1);

    free(res);
}

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

    fds[0].fd = socketfd[GPS];
    fds[1].fd = socketfd[GALILEO];
    fds[2].fd = STDIN_FILENO;
    fds[0].events = fds[1].events = fds[2].events = POLLIN;

    //fcntl(0, F_SETFL, fcntl(0, F_GETFL) | O_NONBLOCK);

    wait_read[0] = wait_read[1] = 0;



    while(1){
        first_input = 0;
        check_sols = 0;
        //usleep(300000);
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
                //strcpy(buf[i], "");
                //memset(buf[i], 0, sizeof buf);
                //strncpy(dest_string,"",strlen(dest_string));
            usleep(300000); //gives time to the solutions to be read (if one of them is late)
            //ret = poll(fds, 3, timeout_msecs);
            if(1/*fds[i].revents & POLLIN*/){
                if(debug && wait_read[i]){ //Don't read solution i (0:GPS, 1:GALILEO) if the solution in the previous iterations has a higher epoch than "seconds" variable
                    continue;
                }
                nbytes = read_gnss(socketfd[i], &offset[i], buf[i], sizeof buf[i]);

                if(strstr(buf[i], "lat") || strstr(buf[i], "latitude") || strstr(buf[i], "ecef") || strlen(buf[i]) < 5){ //Check if the input string contains gnss data or if it is an empty line or header
                    offset[i] = 0;
                    //strcpy(buf[i], "");
                    //memset(buf[tmp], 0, sizeof buf);
                    for(int j=0; j<MAXSTR;j++){
                        buf[i][j] = 0;
                    }
                    continue;
                }

                if(offset[i] != 0 && (buf[i][offset[i]-1] == '\r' || buf[i][offset[i]-1] == '\n')){ //If incoming data is a full gnss measurement string
                    buf[i][offset[i]-1] = '\0'; //Replace new line with end of string
                    //offset[i] = 0; //Reset offset for next reads
                    check_sols |= i+1;
                    //printf("%s\n", buf[i]);
                    first_input |= i+1;
                    sol[i] = str2gnss(buf[i]);//Parse string to gnss_sol_t structure
                    //strcpy(buf[i], "");
                    //memset(buf[tmp], 0, sizeof buf);
                    for(int j=0; j<MAXSTR;j++){
                        buf[i][j] = 0;
                    }
                }
            }
        }

        if(first_time && sol[GPS].time.week == 0){
            continue;
        }

        if(first_time){
            first_time = 0;
            printf("SEMOR: End initialization\n");
            fflush(stdout);
        }

        if(debug){
            if(seconds == 0){
                seconds = (sol[GPS].time.sec <= sol[GALILEO].time.sec) ? sol[GPS].time.sec : sol[GALILEO].time.sec; //Set current second to the minimum of the epochs of the 2 gnss solution
            }
        }

        process_solutions(check_sols); //Get best solution, output it and use it to calculate next imu position

        if(debug)
            seconds++; //Update current second
    }


}

void start_processing(void){
    int ret;
    int i;
    char path[3][PATH_MAX];
    int err;


    //DEBUG
    //initial_pos.time.week = 2032;
    //initial_pos.time.sec = 273375;

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
        sprintf(path[0], "%s/gps_%s.log", log_dir, str_time);
        sprintf(path[1], "%s/galileo_%s.log", log_dir, str_time);

        sol_file[GPS] = fopen(path[0], "w");

        sol_file[GALILEO] = fopen(path[1], "w");
    }

    sprintf(log_dir, "%s%s", root_path, FILE_PATH);
    file = fopen(log_dir, "w");
    if(file == NULL){
        perror("SEMOR fopen()");
    }
    fflush(file);

    err = pthread_create(&imu_thread_id, NULL, &imu_log, NULL);


    handle_connection();
}
