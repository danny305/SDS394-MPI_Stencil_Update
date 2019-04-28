#include <iostream>
#include <random>
#include "time.h"
#include <ctime>
#include <mpi.h>

using namespace std;



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




void smooth (float* x, float* y, uint64_t dim_size=98306, float a=0.05, float b=0.1, float c=0.4) {
//	cout << "Entered smooth function" << endl;
	uint64_t n = dim_size;
	int ub_i, lb_i, ub_j, lb_j;
	ub_i = 0; ub_j = 0;
	lb_i = n - 1; lb_j = n - 1;
	float new_val;
	for (uint64_t i = 0; i < n*n; i += n){
		for (uint64_t j = i; j < i + n; j++){
			if (j < n || j >= n*(n-1) || j == i || j == i + n -1){ 
				continue;
			}
			
			y[j] = 	  a * (x[j-n-1] + x[j-n+1] + x[j+n-1] + x[j+n+1]) + 
				  b * (x[j-n+0] + x[j+n+0] + x[j-1] + x[j+1]) +
				  c * (x[j]);
//			cout << y[j] << " ";
		}
//		cout << endl;	
	}
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


void main (int argc, char* argv[]) {

	//Initialize MPI
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Status status;
	MPI_Request request;


	int nranks, rank = -1, ierr, recv_buff= -1;
        ierr = MPI_Init(&argc, &argv);
        ierr = MPI_Comm_size(comm, &nranks);
        ierr = MPI_Comm_rank(comm, &rank);
	
	//Begin doing MPI

	srand(time(NULL));
	float * x1;
	float * x; float * y;
	uint64_t n = 100;	
//	uint64_t n = 98306;	
	double array_size = (double)sizeof(float)*n*n/1073741824;

	float a = 0.05, b = 0.1, c = 0.4, t = 0.1;
	uint64_t elm_bel_thres_x_ct, elm_bel_thres_y_ct;
	float elm_bel_thres_x_fr, elm_bel_thres_y_fr;




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












