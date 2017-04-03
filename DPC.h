#include <pthread.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/socket.h>
#include <sys/wait.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <netdb.h>
#include <signal.h>
#include <unistd.h>
#include <errno.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define BUFSIZE 256
#define MAXSLEEP 10
#define INTERVAL 1000


FILE *reza_debug; // info regarding the DPC algo.
FILE *log_file_handle; // info regarding the DPC agent logged here. 
double * constraints; // constraint array
double conv_thr; // conergence threshhold
double old_D; // old budget
double mu;
double epsilon;
double power_target;
double C1;
double C2;
double vl;
double vh;
double node_priority;
double v_old;
double v_new;
double ** et;
double *** eij;
char ** hostname; 
int port; // base port: server listens to port+INDEX
int ** Nh; // neighborhood matrix
int ** A_matrix;
int converge; 
int iter_debug = 1;
int * my_sockets;
int comm_error;
int server_sock;
int serversockfd;
int clientsockfd;
int done; // done == 1 stops the main loop 
int NN; // total number of nodes.
int number_constraints;
int my_index; // servers must be indexed from 0 to NN-1

struct packet {
pthread_t thread_id;
int recv_index; 
};

double cal_ext();
void cal_init();
void camm_init();
void transpose (int rows, int cols, double mat[rows][cols], double mat_trans[cols][rows]);
void zeros(int num, int cols, double ** mat);
void value_vec(int num, double * Vec, double value);
double sum(double * array, int length);
int sum_int(int * array, int length);
double sum_col(double ** mat, int column, int length);
double dot_mult(double * a, int * c, int length); // result = a*c'
void element_mult(double * a, int * c, int length); // a = a.*c
void * rec_send(void * arg);
void  * send_rec(void * arg) ;
void comm_init();
void communicate();
void print_error(const char *msg);
double absolute (double value);
int myNanoSleep(time_t sec, long nanosec);
double my_MAX(double value1, double value2);
double my_MIN(double value1, double value2);
void * DPC_opt(void * arg);
void * WL_monitor(void * arg);
void * power_controller(void * arg);

