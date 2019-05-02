#include <iostream>
#include <random>
#include <sys/time.h>
#include <ctime>
#include <mpi.h>
#include <stdexcept>
#include <math.h>
#include <typeinfo>


#define NP 3

using namespace std;


double mysecond(){

	struct timeval tp;
	struct timezone tzp;
	int i;
	i = gettimeofday(&tp, &tzp);
	return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}



void initialize_arr (float* arr, long dim_size = 98306) {
//	cout << "Inside initialize_arr function" << endl;

	uint64_t n = dim_size;
//	cout << "Size of N: " << n << endl;
	for (uint64_t i =0; i < n*n; i += n){
		for (uint64_t j = i; j < i+n; j++){
			arr[j] = 1. * rand()/RAND_MAX;
//			cout << arr[j] << " ";
		}
//		cout << endl;
	}

//	cout << "Exiting intialize_arr function" << endl;
}



void init_send_arr (float* x1, uint64_t & xdim_size, int & nranks, MPI_Comm & comm_old) {
//	cout << "Inside initialize_arr function" << endl;
	float *  x0;
	MPI_Request req_init;
	MPI_Status stat_init;


	x0 = (float *) malloc(xdim_size * xdim_size * sizeof(float));
//	Iniitializing first array before entering for loop:;
	initialize_arr(x0, xdim_size);
	for (int rank_num = 1; rank_num < nranks; rank_num++){
		if (rank_num % 2 == 1){
			MPI_Isend(&x0[0],xdim_size*xdim_size, MPI_FLOAT, rank_num, 1, comm_old, &req_init);
			initialize_arr(x1,xdim_size);
 			MPI_Wait(&req_init, &stat_init);
		}			
		else if (rank_num % 2 == 0){
			MPI_Isend(&x1[0],xdim_size*xdim_size, MPI_FLOAT, rank_num, 1, comm_old, &req_init);
			initialize_arr(x0,xdim_size);
 			MPI_Wait(&req_init, &stat_init);
		}			
		
		if (rank_num == nranks){
			if (rank_num % 2 == 0){ x1 = x0;}
		
		}	
	}



//	cout << "Exiting intialize_arr function" << endl;
}




void smooth (float* x, float* y, uint64_t x_dim_size=98306, float a=0.05, float b=0.1, float c=0.4) {
//	cout << "Entered smooth function" << endl;
	uint64_t n = x_dim_size;
	int ub_i, lb_i, ub_j, lb_j;
	ub_i = 0; ub_j = 0;
	lb_i = n - 1; lb_j = n - 1;
	float new_val;
	for (uint64_t i = 0, k = 0; i < n*n; i += n){
		for (uint64_t j = i; j < i + n; j++){
			if (j < n || j >= n*(n-1) || j == i || j == i + n - 1){ 
				//y[j] = x[j];
				//cout << "Entered if statement! " << "";
				continue;
			}
			
			y[k] = b * (x[j-n+0] + x[j+n+0] + x[j-1] + x[j+1]) + c * (x[j]);
			//cout << "index: " << j << " value " << y[j] << endl;
			k++;
		}
//		cout << endl;	
	}
/*	
	for (uint64_t i = 0; i < n*n; i += n){
		for (uint64_t j = i; j < i + n; j++){
			if (j < n || j >= n*(n-1) || j == i || j == i + n -1){ 
				continue;
			}
			
			y[j] = b * (x[j-n+0] + x[j+n+0] + x[j-1] + x[j+1]) +
			       c * (x[j]);
//			cout << y[j] << " ";
		}
//		cout << endl;	
	}
*/
//	cout << endl;

}



void count (float* arr, uint64_t dim_size = 98306, uint64_t * num_below=0, float threshold=0.5) {
	uint64_t tot_count=0; uint64_t bel_thres_count=0;
	uint64_t n = dim_size;
	
	for (uint64_t i = 0; i < n*n; i += n){
		for (uint64_t j = i; j < i + n; j++){
			if (arr[j] < threshold){
				bel_thres_count++;
			}
			tot_count++;
		}

	}


	*num_below = bel_thres_count;
//	cout << "Total count after parallel region: " << tot_count << endl;
//	cout << "items below threshold: " << bel_thres_count << endl;
//	cout << "Total items in array:  " << tot_count << endl;

}


void print_x_array(float * x, uint64_t & x_n, int & my_rank, int rank=0){
	if (my_rank == rank){
		cout << "Print x for rank: " << my_rank << endl;
		for (uint64_t i = 0; i <= x_n*x_n; i++){
			cout << x[i] << " " ;
			if ((i+1) % x_n == 0)
				cout << "\n\n";
		}
		cout << "\n\n" << endl;
	}
}	




void print_y_array(float * y, uint64_t & n, int & my_rank, int rank = 0){
	if (my_rank == rank){
		cout << "Print y for rank: " << my_rank << endl;
		for (uint64_t i = 0; i <= n*n; i++){
			cout << y[i] << " " ;
			if ((i+1) % n == 0)
				cout << "\n\n";
		}
	}
}



void main (int argc, char* argv[]) {

	//Initialize MPI
	MPI_Comm comm_old = MPI_COMM_WORLD;
	MPI_Comm comm_cart;
	MPI_Status status;
	MPI_Request request;


	int nranks, rank = -1, ierr, recv_buff= -1;
        ierr = MPI_Init(&argc, &argv);
        ierr = MPI_Comm_size(comm_old, &nranks);
        ierr = MPI_Comm_rank(comm_old, &rank);
	


	//Make sure number of task can make a perfect square
	
	float ps_sq_root = sqrt(nranks);
	if (rank == 0){
		if (ps_sq_root == floor(ps_sq_root)){
			#undef NP
			#define NP ps_sq_root
		}
		//cout << "Type of sqrt val: " << typeid(sqrt(nranks)).name() << endl;
		// This alwyas returns type Decimal even if it is an integer, does not help.
		else {
			cout << "Number of tasks: " << nranks << endl;
			cout << "Sqrt = " << ps_sq_root << endl;
			throw invalid_argument( "Number of tasks does not make a perfect square.");
		}
	}
	



/*	int* dim[2]; 
	ierr = MPI_Dims_create( nranks, 2, *dim);
	
	for (int i = 0; i < 2; i++){
		cout << "dimensions: " << dim[i] << endl;	
	}
*/



	srand(time(NULL));
	float * x; float * y;
//	uint64_t n = 1000;
//	uint64_t n = 98306;	
//	uint64_t n = 49152; // for 2 x 2 matrix
//	uint64_t n = 32768; // for 3 x 3 matrix
//	uint64_t n = 12288; // for 8 x 8 matrix
	uint64_t n = 16384;


	uint64_t x_n = n + 2;	
	double array_size = (double)sizeof(float)*n*n/1073741824;

	float b = 0.1, c = 0.4, t = 0.1;
	uint64_t elm_bel_thres_x_ct, elm_bel_thres_y_ct;
	uint64_t total_elm_bel_thres_x_ct, total_elm_bel_thres_y_ct;
	double total_elm_bel_thres_x_fr, total_elm_bel_thres_y_fr;



	x = (float *) malloc(x_n * x_n * sizeof(float));
	y = (float *) malloc(n * n * sizeof(float));


	
	MPI_Datatype init_array; // l_in_col, l_out_col, r_in_col, r_out_col;	
	MPI_Type_contiguous(x_n*x_n, MPI_FLOAT, &init_array);
	MPI_Type_commit(&init_array);


	
//	initialize_arr(x,x_n);
	if (rank == 0){
		init_send_arr(x, x_n, nranks, comm_old);
	}
	else {
		MPI_Recv(&x[0], 1, init_array, 0, 1, comm_old, &status);
		cout << "RECEIVED INIT ARRAY! Rank: " << rank << endl;
	}
	


// Sendiing rows and columns to all Neighbors.
/* =============================================================================================================================== */

	//Creating Cartesion Topology
	//ivdim is an array that says how many  nodes reside along each dimension ex (3,3) or (2,2)
	//ivper is for periodicity, The value 0 means no periodicity. 
	int ivdim[2] = {NP, NP}, ivper[2] = {0,0};
	ierr = MPI_Cart_create(comm_old, 2, ivdim, ivper, 0,&comm_cart);	
	
	//Check the location of each rank within the new topology
	int coordinates[2]; // This is where you store the coordinates of each rank in the new topology
	ierr = MPI_Cart_coords(comm_cart, rank, 2, coordinates);
	//cout << "rank: " << rank << " dimension" << coordinates[0] << "," << coordinates[1] << endl;
	//cout << "IERR: " << ierr << endl;



	//Create Derived type in order to pass border columns
	float col_recv[n];	
	MPI_Datatype n_column; // l_in_col, l_out_col, r_in_col, r_out_col;	
	MPI_Type_vector(n, 1, x_n, MPI_FLOAT, &n_column);
	MPI_Type_commit(&n_column);


	//CARTESIAN SHIFT
	//create all the source and receive ranks for the 4 cart_shifts (up, down, left, right)
	int src_rank_r, dest_rank_r;
	int src_rank_l, dest_rank_l;
	int src_rank_u, dest_rank_u;
	int src_rank_d, dest_rank_d;
	
	int a1 = 10, a2, b1 = 20, b2, cnt = 1;
	float arr_up[n], arr_down[n];
	MPI_Request req;
	MPI_Status stat;


	//obtain all the ranks for each cart shift for each node.

	//shift right coordinates
	ierr = MPI_Cart_shift(comm_cart, 1, 1, &src_rank_r, &dest_rank_r);
	//cout << "Shift right - My rank: " << rank << " Receiving from: " << src_rank_r << " Sending to: " << dest_rank_r << endl;


	//shift left coordinates
	ierr = MPI_Cart_shift(comm_cart, 1, -1, &src_rank_l, &dest_rank_l);
	//cout << "Shift left - My rank: " << rank << " Receiving from: " << src_rank_l << " Sending to: " << dest_rank_l << endl;


	//shift up coordinates
	ierr = MPI_Cart_shift(comm_cart, 0, -1, &src_rank_u, &dest_rank_u);
	//cout << "Shift up - My rank: " << rank << " Receiving from: " << src_rank_u << " Sending to: " << dest_rank_u << endl;

	
	//shift down 
	ierr = MPI_Cart_shift(comm_cart, 0, 1, &src_rank_d, &dest_rank_d);
	//cout << "Shift down - My rank: " << rank << " Receiving from: " << src_rank_d << " Sending to: " << dest_rank_d << endl;

/* ##################################################################################################################################### */

//Send all the data to the appropriate nodes.

 	//Sending column vector right
	if (dest_rank_r != -1){
		//cout << "dest_rank_r does not equal -1" << endl;	
		ierr = MPI_Isend(&x[(2*x_n)-2], 1, n_column, dest_rank_r, 1, comm_old, &req);
	//	cout << "SENT RIGHT" << " My rank: " << rank << " To: " << dest_rank_r << 
	//	" Sending right column  with first value : " << x[(2*x_n)-2] << endl;
		//MPI_Wait;
	}

	if (src_rank_r != -1){
		//cout << "src_rank_r does not equal -1" << endl;
		ierr = MPI_Irecv(&col_recv[0], n, MPI_FLOAT, src_rank_r, 1, comm_old, &req);
		MPI_Wait(&req, &stat);
		for (int i = x_n, j = 0; i < x_n*(x_n -1); i += x_n, j++){
                       x[i] = col_recv[j];
                }
	//	cout << "RECEIVED RIGHT! My rank: " << rank << " Receiving from: " << src_rank_r << 
	//	" first value of array: " << col_recv[0] << endl;
	}
		
//	MPI_Barrier(comm_old);

/*	//Check that column vectors are successfully being shifted right.		
	if (rank == 0){
		cout << "Send rank: " << rank << endl;
		for (int i = (2*x_n)-2; i < x_n*(x_n -1); i += x_n){
			cout << x[i] << ", ";
		}
		cout << endl << endl;
	}

	
	if (rank == 1){
		cout << "Receive rank: " << rank << endl;
		for (int i = 0; i < n; i++){
			cout << col_recv[i] << ", ";
		}
		cout << endl << endl;
	}
*/	
	

/* #################################################################################################### */

	//Sending a column vector left
	if (dest_rank_l != -1){
		//cout << "dest_rank_r does not equal -1" << endl;	
		ierr = MPI_Isend(&x[x_n+1], 1, n_column, dest_rank_l, 1, comm_old, &req);
	//	cout << "SENT LEFT "<< "My rank: " << rank << " To: " << dest_rank_l << 
	//	" Sending left column  with first value : " << x[x_n+1] << endl;
		//MPI_Wait(&req, &stat);
	}

	if (src_rank_l != -1){
		//cout << "src_rank_r does not equal -1" << endl;
		ierr = MPI_Irecv(&col_recv[0], n, MPI_FLOAT, src_rank_l, 1, comm_old, &req);
		MPI_Wait(&req, &stat);
		for (int i = (2*x_n)-1, j = 0; i < x_n*(x_n -1); i += x_n, j++){
                       	x[i] = col_recv[j];
			//cout << x[i] << ", ";
                }
	//	cout << "RECEIVED LEFT! My rank: " << rank << " Receiving from: " << src_rank_l << 
	//	" first value of array: " << col_recv[0] << endl;
	}
		
//	MPI_Barrier(comm_old);

	//Check that column vectors are successfully being shifted LEFT		
/*
	if (rank == 1){
		cout << "Send rank: " << rank << endl;
		for (int i = x_n+1; i < x_n*(x_n -1); i += x_n){
			cout << x[i] << ", ";
		}
		cout << endl << endl;
	}


	if (rank == 0){
		cout << "Receive rank: " << rank << endl;
		for (int i = 0; i < n; i++){
			cout << col_recv[i] << ", ";
		}
		cout << endl << endl;
	}
*/




/* #################################################################################################### */
	
	//Sending up a row	
	if (dest_rank_u != -1){
		//cout << "dest_rank_r does not equal -1" << endl;	
		ierr = MPI_Isend(&x[n+1], n, MPI_FLOAT, dest_rank_u, 1, comm_old, &req);
	//	cout << "SENDING UP! Rank: " << rank << " sent up To: " << dest_rank_u << 
	//	" first value of top row: " << x[n+1] <<endl;
		MPI_Wait(&req,&stat);
	}

	if (src_rank_u != -1){
		//cout << "src_rank_r does not equal -1" << endl;
		ierr = MPI_Irecv(&x[(n+1)*(n+2)+1], n, MPI_FLOAT, src_rank_u, 1, comm_old, &req);
		MPI_Wait(&req, &stat);
		//MPI_Barrier;
	//	cout << "RECEIVED UP! My rank: " << rank << " receiving from: " << src_rank_u << 
	//	" Value received: " << x[(n+1)*(n+2)+1] << endl;
	}

/* #################################################################################################### */

	//Sending down a row
	if (dest_rank_d != -1){
		//cout << "dest_rank_r does not equal -1" << endl;	
		ierr = MPI_Isend(&x[(n)*(n+2)+1], n, MPI_FLOAT, dest_rank_d, 1, comm_old, &req);
	//	cout << "SENDING DOWN! Rank: " << rank << " Sent down to: " << dest_rank_d <<  
	//	" first value of bottom row: " << x[(n)*(n+2)+1] << endl;
		MPI_Wait(&req,&stat);
	}

	if (src_rank_d != -1){
		//cout << "src_rank_r does not equal -1" << endl;
		ierr = MPI_Irecv(&x[1], n, MPI_FLOAT, src_rank_d, 1, comm_old, &req);
		MPI_Wait(&req, &stat);
		//MPI_Barrier;
	//	cout << "RECEIVED DOWN! My rank: " << rank << " Receiving from: " << src_rank_d << 
	//	" Value received: " << x[1] << endl;
	}

/* ================================================================================================================================= */



	smooth(x, y ,x_n);


 	//Prints out the x matrix of a node
//	print_x_array(x, x_n, rank, 0);

 	//Prints out the y matrix of a node.
//	print_y_array(y, n, rank, 0);


	
	ierr = MPI_Barrier(comm_old);
	
	count(x, x_n, &elm_bel_thres_x_ct, t);
	if (rank == 0){ cout << "Element count below threshold in X: " << elm_bel_thres_x_ct << endl;}

	count(y, n, &elm_bel_thres_y_ct, t);
	if (rank == 0){ cout << "Element count below threshold in Y: " << elm_bel_thres_y_ct << endl;}

	

	ierr= MPI_Barrier(comm_old);
	//Send X array below threshold count
	MPI_Reduce(&elm_bel_thres_x_ct, &total_elm_bel_thres_x_ct, nranks, MPI_INT, MPI_SUM, 0, comm_old);
	MPI_Reduce(&elm_bel_thres_y_ct, &total_elm_bel_thres_y_ct, nranks, MPI_INT, MPI_SUM, 0, comm_old);
	ierr = MPI_Barrier(comm_old);
	
/*	if (rank == 0){
		cout << "Total Elements in X below threshold: " << total_elem_bel_thres_x_ct << endl;
		cout << "Total Elements in Y below threshold: " << total_elem_bel_thres_y_ct << endl;
		cout << "Total tasks: " << nranks << endl;
	}
*/

	if (rank == 0){
	double x_divisor =( x_n * x_n) * nranks;	
	double y_divisor =(n * n) * nranks;	
	total_elm_bel_thres_x_fr = total_elm_bel_thres_x_ct/(x_n*x_n*nranks);
	total_elm_bel_thres_y_fr = total_elm_bel_thres_y_ct/(n*n*nranks);
	

	cout << "Summary" << endl;
	cout << "---------------------------------------------------------------------" << endl;
	cout << "Total number of tasks				::	" << nranks << endl;
	cout << "Number of elements in a row/column /task	::	" << x_n << endl;
	cout << "Number of inner elements in a row/column /task	::	" << n << endl;
	cout << "Total number of inner elements			::	" << n*n*nranks << endl;
	cout << "Number of inner elements/task			::	" << n*n << endl;
	cout << "Memory (GB) used per array			::	" << array_size << endl;
	cout << "Threshold					::	" << t << endl;
	cout << "Smoothing constants (b,c)			::	" << b << ", "<< c << endl;
	cout << "Total number of elements below threshold (X)	::	" << total_elm_bel_thres_x_ct << endl;
	cout << "Total fraction of elements below threshold(X)	::	" <<( total_elm_bel_thres_x_ct/(x_divisor))  << endl;
	cout << "Total number of elements below threshold (Y)	::	" <<total_elm_bel_thres_y_ct << endl;
	cout << "Total fraction of elements below threshold(Y)	::	" <<( total_elm_bel_thres_y_ct/(y_divisor)) << endl;
	cout << "---------------------------------------------------------------------" << endl;


	}




	
	ierr = MPI_Finalize();	
}























/*
	clock_t t_alloc_x, t_alloc_y, t_init_x, t_smooth_y, t_count_x, t_count_y;	
	float T_main, T_alloc_x, T_alloc_y, T_init_x, T_smooth_y, T_count_x, T_count_y;	


	double r_t_main = omp_get_wtime();

	clock_t t_main = clock();

	t_alloc_x = clock();
	x = (float *) malloc(n * n * sizeof(size_t));
	T_alloc_x = (clock() - t_alloc_x)/(float)CLOCKS_PER_SEC;

		
	t_alloc_y = clock();
	y = (float *) malloc(n * n * sizeof(size_t));
	T_alloc_y = (clock() - t_alloc_y)/(float)CLOCKS_PER_SEC;

	
	t_init_x = clock();
	initialize_arr(x,n);
	T_init_x = (clock() - t_init_x)/(float)CLOCKS_PER_SEC;

	
	t_smooth_y = clock();
	smooth(x, y, n);
	T_smooth_y = (clock() - t_smooth_y)/(float)CLOCKS_PER_SEC;


	

	t_count_x = clock();
	count(x, n, &elm_bel_thres_x_ct, t);
	T_count_x = (clock() - t_count_x)/(float)CLOCKS_PER_SEC;	

	elm_bel_thres_x_fr = (float)elm_bel_thres_x_ct/(n*n);
	

	t_count_y = clock();
	count(y, n, &elm_bel_thres_y_ct, t);
	T_count_y = (clock() - t_count_y)/(float)CLOCKS_PER_SEC;
	

	elm_bel_thres_y_fr = (float)elm_bel_thres_y_ct/(n*n);


	T_main = (clock() - t_main)/(float)CLOCKS_PER_SEC;
*/	
	
/*
	cout << "Summary" << endl;
	cout << "---------------------------------------------------------------------" << endl;
	cout << "Number of elements in a row/column		::	" << n << endl;
	cout << "Number of inner elements in a row/column	::	" << n-2 << endl;
	cout << "Total number of elements			::	" << n*n << endl;
	cout << "Number of inner elements			::	" << (n-2)*(n-2) << endl;
	cout << "Memory (GB) used per array			::	" << array_size << endl;
	cout << "Threshold					::	" << t << endl;
	cout << "Smoothing constants (a,b,c)			::	" << a << " " << b << " "<< c << endl;
	cout << "Number of elements below threshold (X)		::	" << elm_bel_thres_x_ct << endl;
	cout << "Fraction of elements below threshold		::	" << elm_bel_thres_x_fr << endl;
	cout << "Number of elements below threshold (Y)		::	" << elm_bel_thres_y_ct << endl;
	cout << "Fraction of elements below threshold		::	" << elm_bel_thres_y_fr << endl;
	cout << "---------------------------------------------------------------------" << endl;
	cout << "Action		::	time/s" << endl;
	cout << "-----------------------------" << endl;
	cout << "CPU: Alloc-X	::	" << T_alloc_x << endl;
	cout << "CPU: Alloc-Y	::	" << T_alloc_y << endl;
	cout << "CPU: Init-X	::	" << T_init_x << endl;
	cout << "CPU: Smooth	::	" << T_smooth_y << endl;
	
	
	cout << "CPU: Count-X	::	" << T_count_x << endl;
	
	
	cout << "CPU: Count-Y	::	" << T_count_y << endl;
	
	
	cout << "CPU:Main Time	::	" << T_main << endl;
	

}

*/












