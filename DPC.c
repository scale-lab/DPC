/*
 * DPC.c
 * 
 * Copyright 2017 scaleLab <reza_azimi@brown.edu>, <mbadieik@seas.harvard.edu>
 * <xin_zhan@brown.edu>, <nali@seas.harvard.edu>, and <sherief_reda@brown.edu> 
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */
#include "DPC.h"


/*
 * initialized the optimization
 * 
 */
void cal_init() {
	int i,j;
	double initial_error;
	v_old = vl;
	et = (double **) malloc(number_constraints * sizeof(double *));
	for (i = 0; i < number_constraints; i++) {
		et[i] = malloc(NN * sizeof(double));
		//initial_error = (dot_mult(vl, A_matrix[i], NN) - constraints[i])/sum_int(A_matrix[i],NN);// for vector vl
		initial_error = (vl*sum_int(A_matrix[i],NN) - constraints[i])/sum_int(A_matrix[i],NN);
		value_vec(NN, et[i], initial_error);
		element_mult(et[i], A_matrix[i], NN);
	}
	eij = (double ***) malloc (number_constraints * sizeof(double **));
	for (i = 0; i < number_constraints; i++) {
		eij[i] = (double **) malloc (NN * sizeof(double *));
		for (j =0; j < NN; j++) {
			eij[i][j] = malloc (NN * sizeof(double)); 
		}
	}
	converge = 0;
	old_D = constraints[0];
}

/*
 * initialized the communication between agents
 * 
 */
void comm_init() {
	int num_servers =0; 
	int num_clients = 0; 
	for (int i = 0; i < my_index; i ++) {
		if (Nh[my_index][i] == 1) {
			num_servers++;
		}
	}
	for (int i = my_index+1; i < NN; i++) {
		if (Nh[my_index][i] == 1) {
			num_clients++;
		}
	}
	fprintf(log_file_handle, "%d want to connect to me and I have to connect to %d nodes!\n", num_servers, num_clients);
	if (num_servers != 0) {
		struct sockaddr_in serv_addr;
		server_sock = socket(AF_INET, SOCK_STREAM, 0);
		if (server_sock < 0) 
			print_error("ERROR opening socket");
		bzero((char *) &serv_addr, sizeof(serv_addr));
		serv_addr.sin_family = AF_INET;
		serv_addr.sin_addr.s_addr = INADDR_ANY;
		serv_addr.sin_port = htons(port+my_index);
		if (bind(server_sock, (struct sockaddr *) &serv_addr,sizeof(serv_addr)) < 0) 
			print_error("ERROR on binding");
		listen(server_sock,NN);
		struct sockaddr_in client_address;
		char clntName[INET_ADDRSTRLEN];
		socklen_t client_len = sizeof(client_address);
		for (int i = 0; i < num_servers; i++) {
			int temp = accept(server_sock, (struct sockaddr *)&client_address, &client_len);
			if (temp < 0) print_error("ERROR on accept");
			if(inet_ntop(AF_INET,&client_address.sin_addr.s_addr,clntName,sizeof(clntName))!=NULL){
				fprintf(log_file_handle,"client is %s%c%d\n",clntName,'/',ntohs(client_address.sin_port));
				for (int j = 0; j < NN; j++) {
					if (strcmp (clntName,hostname[j]) == 0) {
						my_sockets[j] = temp;
					}
				}
			} else {
				fprintf(log_file_handle, "Unable to get address\n"); 
			}
		}	
	}
	if (num_clients != 0) {
		struct sockaddr_in serveraddr;
		struct hostent *server;
		for (int i = my_index+1; i < NN; i++) {
			if (Nh[my_index][i] == 1) {
				fprintf(log_file_handle, "trying to connect to %s\n",hostname[i]);
				my_sockets[i] = socket(AF_INET, SOCK_STREAM, 0);
				if (my_sockets[i] < 0) 
					print_error("ERROR opening socket");
				server = gethostbyname(hostname[i]);
				if (server == NULL) {
				       fprintf(log_file_handle,"ERROR, no such host as %s\n", hostname[i]);
					exit(0);
				}
				bzero((char *) &serveraddr, sizeof(serveraddr));
				serveraddr.sin_family = AF_INET;
				bcopy((char *)server->h_addr, (char *)&serveraddr.sin_addr.s_addr, server->h_length);
				serveraddr.sin_port = htons(port+i);
				// retry
				for (int nsec = 1; nsec<= MAXSLEEP; nsec = nsec+1) {
					if (connect(my_sockets[i], (struct sockaddr *) &serveraddr, sizeof(serveraddr)) == 0) {
						break; 
					} else if (nsec <= MAXSLEEP) {
						fprintf(log_file_handle, "connection refused, wait for %d and connect again\n",nsec);
						sleep(1);			
					} else {
						print_error("ERROR connecting comm init");
					}
				}
			}
		} 
	}
	fprintf(log_file_handle, "DPC Conections Complete.\n");
}

/*
 * The core iteration of DPC
 * 
 */
double cal_ext() {
	struct timespec stop, start, comm_start, comm_stop;
	clock_gettime(CLOCK_REALTIME, &start);
	double vhatt;
	double error_sum= 0;
	double et_new[number_constraints];
	int tot_neigh =0; 
	int conv_neigh = 0;
	int constraint_counter = 0;
	int i, j, ii, comm_time_result;
	for (i = 0; i < number_constraints; i++ ) {
		error_sum += my_MAX(0,et[i][my_index])*A_matrix[i][my_index];
	}
	 
	if ( old_D !=  constraints[0]) {
		fprintf(log_file_handle,"iter %d and new budget is %f and old budget was %f\n",iter_debug,constraints[0], old_D);
		fprintf(log_file_handle, "old et_old = %f \n" , et[0][my_index]);		
		et[0][my_index] = et[0][my_index] - (constraints[0] - old_D);//inject the difference
		old_D = constraints[0]; 
		fprintf(log_file_handle, "new et_old = %f \n" , et[0][my_index]);
	}
	for (i = 0; i < number_constraints; i++ ) {
		error_sum += my_MAX(0,et[i][my_index])*A_matrix[i][my_index];
	}
	vhatt = -epsilon * ((2*C1*v_old+C2)*pow(2, node_priority) + 2*mu*error_sum);		
	clock_gettime(CLOCK_REALTIME, &comm_start);
	communicate();
	if (comm_error != 0) {
		return v_old; //repeat
	} 
	clock_gettime(CLOCK_REALTIME, &comm_stop);
	comm_time_result = comm_stop.tv_sec * 1000000000 + comm_stop.tv_nsec -(comm_start.tv_sec * 1000000000 + comm_start.tv_nsec);
	for (i = 0; i < number_constraints; i++ ) {
			zeros(NN, NN, eij[i]);
	}
	for (i = 0; i < NN; i++) {
		for (j = 0; j < NN; j++) {
			if (Nh[i][j] == 1) {
				for (ii = 0; ii < number_constraints; ii++) {
				 eij[ii][i][j] = epsilon*2*mu*(et[ii][i]-et[ii][j]);
				}
			}	
		}
	}
	v_new = v_old+vhatt;
	v_new = my_MAX(my_MIN(v_new, vh),vl);
	vhatt = v_new-v_old;
	for (i = 0; i < number_constraints; i++) {
			et_new[i] = et[i][my_index] + sum_col(eij[i],my_index,NN) - sum (eij[i][my_index],NN) + vhatt * A_matrix[i][my_index];
	}
	for (i = 0; i < NN; i++) {
		if (Nh[my_index][i] == 1) {
			tot_neigh ++;
			constraint_counter = 0;
			for (ii = 0; ii < number_constraints; ii++ ) {
				if (absolute(et[ii][my_index] - et[ii][i]) < conv_thr){
					constraint_counter++;
				}
			}
			if (constraint_counter == number_constraints) {
				conv_neigh ++;
			}
		}
	}
	if ((conv_neigh == tot_neigh) && vhatt < .05) {
		if (converge == 0) {
			fprintf(log_file_handle,"converged at %d\n",iter_debug);
		}
		converge = 1;  
	} else {
		if (converge == 1) {
			fprintf(log_file_handle,"uconverged at %d\n",iter_debug);
		}
		converge = 0; 
	}
	for (i = 0; i < number_constraints; i++) {
		et[i][my_index] = et_new[i];
	}
	power_target = v_new;
	v_old = v_new;
	clock_gettime(CLOCK_REALTIME, &stop);
	int time_result = stop.tv_sec * 1000000000 + stop.tv_nsec -(start.tv_sec * 1000000000 + start.tv_nsec);
	fprintf(reza_debug,"%f\t%f\t%f\t%d\t%d\t%lld\t%d\t%d\t%d\n",v_new,et[0][my_index],constraints[0], time_result, comm_time_result, (long long)time(NULL),converge, tot_neigh, conv_neigh);
	iter_debug ++;
	struct timespec delay;
	if (converge == 1) {
		delay.tv_sec = 0;
		int temp = stop.tv_nsec-start.tv_nsec;
		if ((stop.tv_sec > start.tv_sec && temp > 0) || (stop.tv_sec > start.tv_sec+1)) {
			//fprintf(reza_debug, "took more than it should!!\n");
			//continue;
			return v_new;
		} else if (stop.tv_sec > start.tv_sec && temp <= 0) {
			stop.tv_nsec = stop.tv_nsec + 1000000000;
		}
		delay.tv_nsec = 100*1000000-(stop.tv_nsec-start.tv_nsec);
		if((myNanoSleep(delay.tv_sec, delay.tv_nsec)) == -1) {
			fprintf(log_file_handle,"sleep error:%s\n",strerror(errno));
		}
	}
	return v_new;
}


/*
 * communicate with the neighbors
 * 
 */
void communicate() {	
	struct packet * tinfo;
	tinfo = calloc(NN, sizeof(struct packet));
	if (tinfo == NULL)
               print_error("calloc");
	int tnum = 0;
	comm_error = 0;
	for (int i = 0; i < my_index; i++) {
		if (Nh[my_index][i] == 1) {		
			tinfo[tnum].recv_index =  i;
			pthread_create(&tinfo[tnum].thread_id,0, rec_send, &tinfo[tnum]);
			tnum++;
		}	
	}
	for (int i = my_index+1; i < NN; i++) {
		if (Nh[my_index][i] == 1) {	
			tinfo[tnum].recv_index =  i;	 
			pthread_create(&tinfo[tnum].thread_id,0, send_rec, &tinfo[tnum]);
			tnum++;
		}	
	}
	for (int i = 0; i < tnum; i++) {
		pthread_join(tinfo[i].thread_id, 0);
	}
	free(tinfo);
}

/*
 * rec and sern a message
 * 
 */
void * rec_send(void * arg)
{
	struct packet * input = arg;
	int n1,i;
	float temp_float;
	char send_buf[BUFSIZE]="";
	char * recv_ptr;
	char buffer[BUFSIZE]="";
	char * temp_buf;
	const char sp[2] = "a";
	for (i = 0; i < number_constraints; i++) {
        if (sprintf(buffer, "%f%s", (float) et[i][my_index], sp) <= 0) {
		    fprintf(log_file_handle, "%d,%d: server send1 sprintf error!\n",iter_debug,i);
        }
        strcat(send_buf, buffer);
	}
	//buffer="";
	if ((n1=recv(my_sockets[input->recv_index],buffer,BUFSIZE,0)) <= 0) {
		fprintf(log_file_handle, "%d: server recv:%s %lld\n",iter_debug,strerror(errno), (long long)time(NULL));
		Nh[my_index][input->recv_index] = 0; 
		Nh[input->recv_index][my_index] = 0;
		comm_error = 1; 
		return (0);
	}
	if ((send(my_sockets[input->recv_index],send_buf,strlen(send_buf),0)) <= 0) {
		fprintf(log_file_handle,"%d: server send:%s %lld\n",iter_debug,strerror(errno), (long long)time(NULL));
		Nh[my_index][input->recv_index] = 0; 
		Nh[input->recv_index][my_index] = 0; 
		comm_error = 1;
		return (0);
	}
	buffer[n1] = 0;
	recv_ptr = buffer;
	temp_buf = strtok_r(recv_ptr,sp, &recv_ptr); 
	for (i = 0; i < number_constraints; i++) {
	    if (sscanf(temp_buf,"%f",&temp_float) <= 0) {
		    fprintf(log_file_handle,"%d,%d: server send1 scanf error:%s\n",iter_debug,i, temp_buf);
	    }   
        et[i][input->recv_index]=(double) temp_float;
        temp_buf = strtok_r(recv_ptr,sp, &recv_ptr);
	}
	return (0);
}


/*
 * send and recieve a message
 * 
 */
void  * send_rec(void * arg) {
	struct packet * input = arg;
	int n1,i;
	float temp_float;
	char send_buf[BUFSIZE]="";
	char buffer[BUFSIZE]="";
	char * recv_ptr;
	char * temp_buf;
	const char sp[2] = "a";
	//char sp1, sp2; 
	for (i = 0; i < number_constraints; i++) {
        if (sprintf(buffer, "%f%s", (float) et[i][my_index], sp) <= 0) {
		//if (sprintf(buffer, "%f%s", (float) 0.2, sp) <= 0) {	
		    fprintf(log_file_handle,"%d,%d: server send1 sprintf error!\n",iter_debug,i);
        }
        strcat(send_buf, buffer);
	}
	if ((send(my_sockets[input->recv_index],send_buf,strlen(send_buf),0)) <= 0) {
		fprintf(log_file_handle, "%d: client send1: %s %lld \n",iter_debug,strerror(errno),(long long)time(NULL));	
		Nh[my_index][input->recv_index] = 0; 
		Nh[input->recv_index][my_index] = 0;
		comm_error = 1; 
	}
	//buffer="";
	if((n1=recv(my_sockets[input->recv_index],buffer,BUFSIZE,0))<= 0) {
		fprintf(log_file_handle, "%d:  rec1:%s %lld\n",iter_debug,strerror(errno),(long long)time(NULL));
		Nh[my_index][input->recv_index] = 0; 
		Nh[input->recv_index][my_index] = 0;
		comm_error = 1; 
	}
	buffer[n1]=0;
	recv_ptr = buffer; 
	temp_buf = strtok_r(recv_ptr,sp, &recv_ptr); 
	for (i = 0; i < number_constraints; i++) {
        if (sscanf(temp_buf,"%f",&temp_float) <= 0) {
		    fprintf(log_file_handle, "%d,%d: server send1 scanf error:%s\n",iter_debug,i, temp_buf);
	    }   
        et[i][input->recv_index]=(double) temp_float;
		temp_buf = strtok_r(recv_ptr,sp, &recv_ptr);
	}
	return (0);
}



double absolute(double value) {
  if (value < 0) {
    return -value;
  }
  else {
    return value;  
  }
}

double my_MAX(double value1, double value2) {
  if (value1 < value2) {
    return value2;
  }
  else {
    return value1;  
  }
}

double my_MIN(double value1, double value2) {
  if (value1 < value2) {
    return value1;
  }
  else {
    return value2;  
  }
}

void transpose (int rows, int cols, double mat[rows][cols], double mat_trans[cols][rows]) {
    int i, j;
    for (i = 0; i< cols; i++){
        for (j = 0; j< rows; j++){
            mat_trans[i][j] = mat[j][i];
        }   
    }
}

void zeros(int num, int cols, double ** mat) {
	int i, j;
	for (i = 0; i < num; i++) {
		for (j = 0; j < cols; j++) {
			mat[i][j] = 0.0;	
		}
	}
}

void value_vec(int num, double * Vec, double value) {
	int i;
	for (i = 0; i < num; i++) {
		Vec[i] = value;
	}
}

/*
 * get the sum of an array
 * 
 */
double sum(double * array, int length) {
	double result = 0.0;
	int i;
	for (i = 0; i < length; i++) {
			result += array[i];
	}
	return result;
}

/*
 * get the sum of a column of a matrix 
 * 
 */
double sum_col(double ** mat, int column, int length) {
	double result = 0.0;
	int i;
	for (i = 0; i < length; i++) {
			result += mat[i][column];
	}
	return result;
}

int sum_int(int * array, int length) {
	int result = 0;
	int i;
	for (i = 0; i < length; i++) {
			result += array[i];
	}
	return result;
}

double dot_mult(double * a, int * c, int length) {
	int i;
	double result = 0.0; 
	for (i = 0; i < length; i++) {
			result += (double) (a[i] * c[i]); 
	}
	return result; 
}

/*
 * element wise mult of two arrays
 * 
 */
void element_mult(double * a, int * c, int length) {
	int i;
	for (i = 0; i < length; i++) {
			a[i] = (double) a[i] * c[i]; 
	}
}

int myNanoSleep(time_t sec, long nanosec) {
	/* Setup timespec */
	struct timespec req;
	req.tv_sec = sec;
	req.tv_nsec = nanosec;
	/* Loop until we've slept long enough */
	while(req.tv_sec > 0 || req.tv_nsec > 0) {
		/* Store remainder back on top of the original required time */
		if(nanosleep(&req, &req) != 0) {
		  /* If any error other than a signal interrupt occurs, return an error */
		  if(errno != EINTR)
		     return -1; 
		} else {
		  /* nanosleep succeeded, so exit the loop */
		  break;
		}
	} 
	return 0; /* Return success */
}

/*
 * print error and exit
 * 
 */
void print_error(const char *msg)
{
    fprintf(log_file_handle, "Exit on ERROR:%s\n", msg);
    exit(1);
}

void * DPC_opt(void * arg) {
	mu = .05;
	epsilon = 2;
	conv_thr = .01;
	comm_init();
	cal_init();
	int i,j;
	while (!done) { // main loop
		power_target = cal_ext(0,0,0);
		
	}
	fprintf(log_file_handle,"closing sockets at %lld\n",(long long)time(NULL));	
	for (i = 0; i < NN; i++) {
		if (Nh[my_index][i]) {
			close(my_sockets[i]);
		}
	}
	for (i = 0; i < number_constraints; i++) {
		for (j =0; j < NN; j++) {
			free(eij[i][j]);
		}
		free(eij[i]);
		free(et[i]);
	}
	free(eij);
	free(et);
	close(server_sock);
	close(serversockfd);
	close(clientsockfd);
	return(0);
}

/*
 * Monitor the workload and give this info to DPC 
 * 
 */
void * WL_monitor(void * arg) {
	// IN DPC originial implemnetation Wl_monitor got these information
	// from slrumd but for simplicity we read it from file in case SLURM 
	// is not installed on your cluster. 
	// for the budget we insert the error using server/client model
	// here again we are reading it form a file to keep it simple.  
	FILE * wl_file;
	FILE * budget_file;
	struct timespec tps, tpe;
	struct timespec delay = {INTERVAL/1000,(INTERVAL%1000)*10e6};
	double C1_lookup[] = {0.0006, 0.0001, 0.0006, 0.0013};
	double C2_lookup[] = {-0.2389, -0.0783, -0.2846, -0.4834};
	double max_power = 200;
	double min_power = 130;
	float float_temp;
	int wl_index, priority;
	while (!done) {
		clock_gettime(CLOCK_REALTIME, &tps);
		wl_file=fopen("workload.txt", "r");
		budget_file = fopen("budget.txt","r");
		if (wl_file == NULL || budget_file == NULL) {
        	continue;// try again 
		}
		if (!fscanf(wl_file, "%d\n%d", &wl_index, &priority)) {
			continue;//try again
		}
		if (!fscanf(budget_file, "%f",  &float_temp)) {
			continue;//try again
		}
		C1 = C1_lookup[wl_index];
		C2 = C2_lookup[wl_index];
		node_priority = priority;
		if (constraints[0] != float_temp){
			fprintf(log_file_handle, "budget updated, old=%f, new=%f\n",constraints[0],float_temp);
			constraints[0] = (double) float_temp;
		}
		vl= min_power;
		vh= max_power;
		clock_gettime(CLOCK_REALTIME, &tpe);
		delay.tv_sec = 0;
		long temp = tpe.tv_nsec-tps.tv_nsec;
		if ((tpe.tv_sec > tps.tv_sec && temp > 0) || (tpe.tv_sec > tps.tv_sec+1)) {
			fprintf(log_file_handle, "took more than it should!!\n");
			continue;
		} else if (tpe.tv_sec > tps.tv_sec && temp <= 0) {
			tpe.tv_nsec = tpe.tv_nsec + 1000000000;
		}
		delay.tv_nsec = INTERVAL*1000000-(tpe.tv_nsec-tps.tv_nsec);
		if((myNanoSleep(delay.tv_sec, delay.tv_nsec)) == -1) {
			fprintf(log_file_handle,"sleep error:%s\n",strerror(errno));
		}
	}
	return 0; 
	
}

void * power_controller(void * arg) {
	// this is place holder for the power controller, based on your infrastructer
	// read the power consumption and enfore the power_target. Here we just write to a file;
	FILE * pc_file=fopen("powerCap.txt", "a"); 
	while (pc_file == NULL) {
			fprintf(log_file_handle, "Can't open the powerCap.txt, try again ...\n");
			sleep(1);
        	pc_file=fopen("powerCap.txt", "a") ;// try again 
	}
	while (!done) {
		fprintf(pc_file, "%f\n",(float) power_target);
		sleep(1);
	}
	fclose(pc_file);
	return(0);
}

int hostname_to_ip(char * hostname , char* ip)
{
    struct hostent *he;
    struct in_addr **addr_list;
    int i;
         
    if ( (he = gethostbyname( hostname ) ) == NULL) 
    {
        // get the host info
        print_error("ERROR gethostbyname");
        return 1;
    }
 
    addr_list = (struct in_addr **) he->h_addr_list; 
    for(i = 0; addr_list[i] != NULL; i++) 
    {
        strcpy(ip , inet_ntoa(*addr_list[i]) );
        return 0;
    }
     
    return 1;
}


void read_config(char * file_name){
	FILE * topo_file;
	int value, n1, i, j;
	float D_float;
    char ip_buffer[30];
	topo_file=fopen(file_name, "r");
	if (topo_file == NULL) {
        	print_error("Error on opening configuration file!");
	}
	n1 = fscanf(topo_file, "%d", &my_index);
	if (!n1) {
		print_error("ERROR reading my_index!");
	}
	fprintf(log_file_handle,"my index is: %d\n",my_index);
	n1 = fscanf(topo_file, "%d", &NN);
	if (!n1) {
		print_error("ERROR reading number of Nodes!");
	}
	fprintf(log_file_handle,"Number of Node is: %d\n",NN);
	my_sockets = (int *) malloc (NN * sizeof(int));
	Nh = (int **)malloc(NN * sizeof(int *));
	for (i=0; i<NN; i++) {
         Nh[i] = (int *)malloc(NN * sizeof(int));
    }
	for (i = 0; i < NN; i++) {
		for(j = 0; j < NN; j++) {
			n1 = fscanf(topo_file, "%d", &value);
			if (!n1) {
				print_error("ERROR reading topology");
				break; 
			} else {
				Nh[i][j] = value;
			} 
		}
	} 
	for (i = 0; i < NN; i++) {
		for(j = 0; j < NN; j++) {
			fprintf(log_file_handle,"%d\t", Nh[i][j]); 
		}
		fprintf(log_file_handle,"\n");
	}	
	hostname = (char **)malloc(NN * sizeof(char *));
	for (i = 0; i < NN; i++) {
		char * temp = malloc (30 * sizeof(char));
		n1 = fscanf(topo_file, "%s", ip_buffer);
		if (!n1) {
			print_error("ERROR reading the ip of the neighbors"); 
		}
		if (hostname_to_ip(ip_buffer , temp)) {
		 	print_error("ERROR host to ip");
		}
		hostname[i] = temp;
		fprintf(log_file_handle, "ip %i, %s %s \n",i, ip_buffer, hostname[i]);
	}
	n1 = fscanf(topo_file, "%d",&port);
        if (!n1) {
        	print_error("ERROR reading the port base");
        }
        fprintf(log_file_handle, "base port is %d \n",port );



	n1 = fscanf(topo_file, "%d", &number_constraints);
	if (!n1) {
		print_error("ERROR reading number of constraints!");
	}
	fprintf(log_file_handle,"Number of constraints is: %d\n",number_constraints);
	constraints = malloc(number_constraints * sizeof(double));
	for (i= 0; i < number_constraints; i++) {
		n1 = fscanf(topo_file, "%f", &D_float);
		if (!n1) {
			print_error("ERROR reading initial constraints!");
		}
		constraints[i] = (double) D_float; 
		fprintf(log_file_handle, "CB%d has budget %f\n",i,constraints[i]);
	}
	A_matrix = (int **)malloc(number_constraints * sizeof(int *));
	for (i=0; i<number_constraints; i++) {
         A_matrix[i] = (int *)malloc(NN * sizeof(int));
    }
	for (i = 0; i < number_constraints; i++) {
		for(j = 0; j < NN; j++) {
			n1 = fscanf(topo_file, "%d", &value);
			if (!n1) {
				print_error("ERROR reading constraint matrix");
				break; 
			} else {
				A_matrix[i][j] = value;
			} 
		}
	} 
	for (i = 0; i < number_constraints; i++) {
		for(j = 0; j < NN; j++) {
			fprintf(log_file_handle,"%d\t", Nh[i][j]); 
		}
		fprintf(log_file_handle,"\n");
	}	
	
	fclose(topo_file);
}

void clean_up() {
	int i; 
	for (i = 0; i < NN; i++) {
		free(Nh[i]);
		free(hostname[i]);
	}
	for (i = 0; i < number_constraints; i++) {
		free(A_matrix[i]);
	}
	free(my_sockets);
	free(hostname);
	free(Nh);
	free(A_matrix);
	free(constraints);
}

int main(int argc, char **argv) {
	log_file_handle = fopen("log.txt", "w");
	reza_debug = fopen("DPC_debug.txt", "w");
	read_config(argv[1]);
	done = 0; 
	pthread_t opt_thread, wlmonitor_thread, pc_thread; 
	pthread_create(&opt_thread, NULL, DPC_opt, NULL);
	pthread_create(&wlmonitor_thread, NULL, WL_monitor, NULL);
	pthread_create(&pc_thread, NULL, power_controller, NULL);
	sleep (30);//DPC agents run forever but here for the sake of example run it for 30 sec.
	done = 1;
	pthread_cancel(opt_thread);
	pthread_cancel(wlmonitor_thread);
	pthread_cancel(pc_thread);
	clean_up();
	fprintf(log_file_handle,"finishing DPC %lld\n",(long long)time(NULL));
	fclose(log_file_handle);
	fclose(reza_debug);
	
}
