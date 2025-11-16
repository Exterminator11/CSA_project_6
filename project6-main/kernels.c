/********************************************************
 * Kernels to be optimized for the  Performance Project
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include <sys/time.h>
#include <sys/resource.h>

#define BLOCK 32
#define BLOCK_SMOOTH 8

/* Below are statements to set up the performance measurement utilities */
/* we use rdtsc, clock, and getusage utilities to measure performance */

//#define rdtscll(val) __asm__ __volatile__("rdtsc" : "=A" (val))
#if defined(__i386__)

static __inline__ unsigned long long rdtsc(void)
{
  unsigned long long int x;
     __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
     return x;
}
#elif defined(__x86_64__)


static __inline__ unsigned long long rdtsc(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}
#endif

/* end of definitions to set up measurement utilities */


/* 
 * Please fill in the following team struct 
 */
team_t team = {
    "5697",              /* Team name */

    "Rachit R. Das",     /* First member full name */
    "rachit.das@gwu.edu",  /* First member email address */

    "",                   /* Second member full name (leave blank if none) */
    ""                    /* Second member email addr (leave blank if none) */
};

/***************
 * ROTATE KERNEL
 ***************/

/******************************************************
 * Your different versions of the rotate kernel go here
 ******************************************************/

/* the getUserTime function is used for measurement, you should not change the code for this function */

long int getUserTime()
{
	int who= RUSAGE_SELF;
	int ret;
	struct rusage usage;
	struct rusage *p=&usage;
	//long int current_time;

	ret=getrusage(who,p);
	if(ret == -1)
	{
		printf("Could not get GETRUSAGE to work in function %s at line %d in file %s\n",
				__PRETTY_FUNCTION__, __LINE__, __FILE__);
		exit(1);
	}
	return (p->ru_utime.tv_sec * 1000000 + p->ru_utime.tv_usec);
}

/* 
 * naive_rotate - The naive baseline version of rotate 
 */
 /* The parameters, pointers, rusage_time, rdtsc_time, and cpu_time_used are used to measure performance and return values to caller. */
 /* You should not change the code that uses these parameters and variables. */
 
char naive_rotate_descr[] = "naive_rotate: Naive baseline implementation";
void naive_rotate(int dim, pixel *src, pixel *dst, int *rusage_time, unsigned long long *rdtsc_time) 
{
	int i, j;
	/* the variables below are used for performance measurement and not for computing the results of the algorithm */
	long int rusage_start_time, rusage_end_time = 0;
	unsigned long long rdtsc_start_time, rdtsc_end_time = 0;
	/* call system functions to start measuring performance. you should not bother to change these. */
	
	rusage_start_time = getUserTime();
	rdtsc_start_time = rdtsc();

/* below is the main computations for the rotate function */

	for (j = 0; j < dim; j++)
		for (i = 0; i < dim; i++)
			dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];

/* the remaining lines in this function stop the measurement and set the values before returning. */

	rusage_end_time = getUserTime();
	rdtsc_end_time = rdtsc();

	*rusage_time = rusage_end_time - rusage_start_time;
	*rdtsc_time = rdtsc_end_time - rdtsc_start_time;
}

 /* The parameters, pointers, rusage_time, rdtsc_time, and cpu_time_used are used to measure performance and return values to caller. */
 /* You should not change the code that uses these parameters and variables. */
char my_rotate_descr[] = "my_rotate: Naive baseline implementation";
void my_rotate(int dim, pixel *src, pixel *dst, int *rusage_time, unsigned long long *rdtsc_time) 
{
	int i, j;
		/* the variables below are used for performance measurement and not for computing the results of the algorithm */
	long int rusage_start_time, rusage_end_time = 0;
        unsigned long long rdtsc_start_time, rdtsc_end_time = 0;
	/* call system functions to start measuring performance. you should not bother to change these. */
        rusage_start_time = getUserTime();
	rdtsc_start_time = rdtsc();

/* ANY CHANGES ARE MADE HERE */
/* below are the main computations for your implementation of the rotate. Any changes in implementation will go here or the other functions it may call */
	// for (j = 0; j < dim; j++)
	// 	for (i = 0; i < dim; i++)
	// 		dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];

	// loop interchange does not work
	// for (i = 0; i < dim; i++)
	// 	for (j = 0; j < dim; j++)
	// 		dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];

	//Blocking
	// for (j = 0; j < dim; j+=BLOCK)
	// 	for (i = 0; i < dim; i+=BLOCK)
	// 		for(int jj=j;jj<j+BLOCK && jj<dim;jj++)
	// 			for(int ii=i;ii<i+BLOCK && ii<dim;ii++)
	// 				dst[RIDX(dim-1-jj, ii, dim)] = src[RIDX(ii, jj, dim)];

	//Blocking + Loop interchange
	// for (i = 0; i < dim; i+=BLOCK)
	// 	for (j = 0; j < dim; j+=BLOCK)
	// 		for(int ii=i;ii<i+BLOCK && ii<dim;ii++)
	// 			for(int jj=j;jj<j+BLOCK && jj<dim;jj++)
	// 				dst[RIDX(dim-1-jj, ii, dim)] = src[RIDX(ii, jj, dim)];

	//Blocking + unrolling(4 times)
	// for (j = 0; j < dim; j+=BLOCK)
	// 	for (i = 0; i < dim; i+=BLOCK)
	// 		for(int jj=j;jj<j+BLOCK && jj<dim;jj++)
	// 			for(int ii=i;ii<i+BLOCK && ii<dim;ii+=4){
	// 				dst[RIDX(dim-1-jj, ii, dim)] = src[RIDX(ii, jj, dim)];
	// 				if(ii+1<dim) dst[RIDX(dim-1-jj, ii+1, dim)] = src[RIDX(ii+1, jj, dim)];
	// 				if(ii+2<dim) dst[RIDX(dim-1-jj, ii+2, dim)] = src[RIDX(ii+2, jj, dim)];
	// 				if(ii+3<dim) dst[RIDX(dim-1-jj, ii+3, dim)] = src[RIDX(ii+3, jj, dim)];
	// 			}

	//Blocking + unrolling(8 times)
	// for (j = 0; j < dim; j+=BLOCK)
	// 	for (i = 0; i < dim; i+=BLOCK)
	// 		for(int jj=j;jj<j+BLOCK && jj<dim;jj++)
	// 			for(int ii=i;ii<i+BLOCK && ii<dim;ii+=8){
	// 				dst[RIDX(dim-1-jj, ii, dim)] = src[RIDX(ii, jj, dim)];
	// 				if(ii+1<dim) dst[RIDX(dim-1-jj, ii+1, dim)] = src[RIDX(ii+1, jj, dim)];
	// 				if(ii+2<dim) dst[RIDX(dim-1-jj, ii+2, dim)] = src[RIDX(ii+2, jj, dim)];
	// 				if(ii+3<dim) dst[RIDX(dim-1-jj, ii+3, dim)] = src[RIDX(ii+3, jj, dim)];
	// 				if(ii+4<dim) dst[RIDX(dim-1-jj, ii+4, dim)] = src[RIDX(ii+4, jj, dim)];
	// 				if(ii+5<dim) dst[RIDX(dim-1-jj, ii+5, dim)] = src[RIDX(ii+5, jj, dim)];
	// 				if(ii+6<dim) dst[RIDX(dim-1-jj, ii+6, dim)] = src[RIDX(ii+6, jj, dim)];
	// 				if(ii+7<dim) dst[RIDX(dim-1-jj, ii+7, dim)] = src[RIDX(ii+7, jj, dim)];

	// 			}

	//Blocking + unrolling(8 times) + common subexpression elimination
	// for (j = 0; j < dim; j+=BLOCK)
	// 	for (i = 0; i < dim; i+=BLOCK)
	// 		for(int jj=j;jj<j+BLOCK && jj<dim;jj++){
	// 			int dIdx=dim-1-jj;
	// 			for(int ii=i;ii<i+BLOCK && ii<dim;ii+=8){
	// 				dst[RIDX(dIdx, ii, dim)] = src[RIDX(ii, jj, dim)];
	// 				if(ii+1<dim) dst[RIDX(dIdx, ii+1, dim)] = src[RIDX(ii+1, jj, dim)];
	// 				if(ii+2<dim) dst[RIDX(dIdx, ii+2, dim)] = src[RIDX(ii+2, jj, dim)];
	// 				if(ii+3<dim) dst[RIDX(dIdx, ii+3, dim)] = src[RIDX(ii+3, jj, dim)];
	// 				if(ii+4<dim) dst[RIDX(dIdx, ii+4, dim)] = src[RIDX(ii+4, jj, dim)];
	// 				if(ii+5<dim) dst[RIDX(dIdx, ii+5, dim)] = src[RIDX(ii+5, jj, dim)];
	// 				if(ii+6<dim) dst[RIDX(dIdx, ii+6, dim)] = src[RIDX(ii+6, jj, dim)];
	// 				if(ii+7<dim) dst[RIDX(dIdx, ii+7, dim)] = src[RIDX(ii+7, jj, dim)];
	// 			}
	// 		}

	// Blocking + unrolling(8 times) + common subexpression elimination + eliminating branches using a singular limit instead if
	// for (j = 0; j < dim; j+=BLOCK)
	// 	for (i = 0; i < dim; i+=BLOCK)
	// 		for(int jj=j;jj<j+BLOCK && jj<dim;jj++){
	// 			int dIdx=dim-1-jj;
	// 			int max_i = (i + BLOCK < dim) ? (i + BLOCK) : dim;
	// 			/* round max_i down to multiple of 8 relative to i */
	// 			int limit = max_i - ((max_i - i) % 8);

	// 			int ii;
	// 			for(ii=i;ii<limit;ii+=8){
	// 				dst[RIDX(dIdx, ii, dim)] = src[RIDX(ii, jj, dim)];
	// 				dst[RIDX(dIdx, ii+1, dim)] = src[RIDX(ii+1, jj, dim)];
	// 				dst[RIDX(dIdx, ii+2, dim)] = src[RIDX(ii+2, jj, dim)];
	// 				dst[RIDX(dIdx, ii+3, dim)] = src[RIDX(ii+3, jj, dim)];
	// 				dst[RIDX(dIdx, ii+4, dim)] = src[RIDX(ii+4, jj, dim)];
	// 				dst[RIDX(dIdx, ii+5, dim)] = src[RIDX(ii+5, jj, dim)];
	// 				dst[RIDX(dIdx, ii+6, dim)] = src[RIDX(ii+6, jj, dim)];
	// 				dst[RIDX(dIdx, ii+7, dim)] = src[RIDX(ii+7, jj, dim)];
	// 			}
	// 			for (; ii < max_i; ii++)
    //             	dst[RIDX(dIdx, ii, dim)] = src[RIDX(ii, jj, dim)];
	// 		}

	//Blocking + unrolling(8 times) + common subexpression elimination + eliminating RIDX to reduce multiplications
	// for (j = 0; j < dim; j+=BLOCK)
	// 	for (i = 0; i < dim; i+=BLOCK)
	// 		for(int jj=j;jj<j+BLOCK && jj<dim;jj++){
	// 			int dIdx=dim-1-jj;
	// 			int dstBase = dIdx * dim;
	// 			for(int ii=i;ii<i+BLOCK && ii<dim;ii+=8){
	// 				dst[dstBase+ii] = src[(ii)*(dim)+jj];
	// 				if(ii+1<dim) dst[dstBase+ii+1] = src[(ii+1)*(dim)+jj];
	// 				if(ii+2<dim) dst[dstBase+ii+2] = src[(ii+2)*(dim)+jj];
	// 				if(ii+3<dim) dst[dstBase+ii+3] = src[(ii+3)*(dim)+jj];
	// 				if(ii+4<dim) dst[dstBase+ii+4] = src[(ii+4)*(dim)+jj];
	// 				if(ii+5<dim) dst[dstBase+ii+5] = src[(ii+5)*(dim)+jj];
	// 				if(ii+6<dim) dst[dstBase+ii+6] = src[(ii+6)*(dim)+jj];
	// 				if(ii+7<dim) dst[dstBase+ii+7] = src[(ii+7)*(dim)+jj];
	// 			}
	// 		}

	// Blocking + loop unrolling(8 times) + Block size 32 + common subexpression elimination + eliminating branches + eliminating RIDX to reduce multiplications + eliminating braches
	for (j = 0; j < dim; j+=BLOCK)
		for (i = 0; i < dim; i+=BLOCK)
			for(int jj=j;jj<j+BLOCK && jj<dim;jj++){
				int max_i = (i + BLOCK < dim) ? (i + BLOCK) : dim;
				/* round max_i down to multiple of 8 relative to i */
				int limit = max_i - ((max_i - i) % 8);
				int dstBase = (dim-1-jj) * dim;
				int ii;
				pixel *srcCol = &src[jj];
				for(ii=i;ii<limit;ii+=8){
					dst[dstBase+ii] = srcCol[ii*dim];
					dst[dstBase+ii+1] = srcCol[(ii+1)*dim];
					dst[dstBase+ii+2] = srcCol[(ii+2)*dim];
					dst[dstBase+ii+3] = srcCol[(ii+3)*dim];
					dst[dstBase+ii+4] = srcCol[(ii+4)*dim];
					dst[dstBase+ii+5] = srcCol[(ii+5)*dim];
					dst[dstBase+ii+6] = srcCol[(ii+6)*dim];
					dst[dstBase+ii+7] = srcCol[(ii+7)*dim];
				}
			}




/* end of computation for rotate function. any changes you make should be made above this line. */
/* END OF CHANGES in this function */

/* the remaining lines in this function stop the measurement and set the values before returning. */
	rusage_end_time = getUserTime();
        rdtsc_end_time = rdtsc();

	*rusage_time = rusage_end_time - rusage_start_time;
	*rdtsc_time = rdtsc_end_time - rdtsc_start_time;
}





/***************
 * SMOOTH KERNEL
 **************/

/***************************************************************
 * Various typedefs and helper functions for the smooth function
 * You may modify these any way you like.
 **************************************************************/

/* A struct used to compute averaged pixel value */
typedef struct {
	int red;
	int green;
	int blue;
	int num;
} pixel_sum;

/* Compute min and max of two integers, respectively */
static int minimum(int a, int b) 
{ return (a < b ? a : b); }
static int maximum(int a, int b) 
{ return (a > b ? a : b); }

/* 
 * initialize_pixel_sum - Initializes all fields of sum to 0 
 */
static void initialize_pixel_sum(pixel_sum *sum) 
{
	sum->red = sum->green = sum->blue = 0;
	sum->num = 0;
	return;
}

/* 
 * accumulate_sum - Accumulates field values of p in corresponding 
 * fields of sum 
 */
static void accumulate_sum(pixel_sum *sum, pixel p) 
{
	sum->red += (int) p.red;
	sum->green += (int) p.green;
	sum->blue += (int) p.blue;
	sum->num++;
	return;
}

/* 
 * assign_sum_to_pixel - Computes averaged pixel value in current_pixel 
 */
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum) 
{
	current_pixel->red = (unsigned short) (sum.red/sum.num);
	current_pixel->green = (unsigned short) (sum.green/sum.num);
	current_pixel->blue = (unsigned short) (sum.blue/sum.num);
	return;
}

/* 
 * avg - Returns averaged pixel value at (i,j) 
 */
static pixel avg(int dim, int i, int j, pixel *src) 
{
	int ii, jj;
	pixel_sum sum;
	pixel current_pixel;

	initialize_pixel_sum(&sum);

	for(ii = maximum(i-1, 0); ii <= minimum(i+1, dim-1); ii++) 
		for(jj = maximum(j-1, 0); jj <= minimum(j+1, dim-1); jj++)
			accumulate_sum(&sum, src[RIDX(ii, jj, dim)]);

	assign_sum_to_pixel(&current_pixel, sum);
	return current_pixel;
}

static pixel my_avg(int dim, int i, int j, pixel *src) 
{
	int ii, jj;
	pixel_sum sum;
	pixel current_pixel;

	// initialize_pixel_sum(&sum);

	//function inlining
	sum.red = sum.green = sum.blue = 0;
	sum.num = 0;

	// for(ii = maximum(i-1, 0); ii <= minimum(i+1, dim-1); ii++) 
	// 	for(jj = maximum(j-1, 0); jj <= minimum(j+1, dim-1); jj++) {
	// 		// accumulate_sum(&sum, src[RIDX(ii, jj, dim)]);
			
	// 		//function inlining
	// 		pixel p=src[RIDX(ii, jj, dim)];
	// 		sum.red += (int) p.red;
	// 		sum.green += (int) p.green;
	// 		sum.blue += (int) p.blue;
	// 		sum.num++;
	// 	}

	//unroll the above loop as only 9 cases 
	// Handle all 9 cases based on position in the grid
	pixel p;

	// Top-left corner (i=0, j=0)
	if (i == 0 && j == 0) {
		// Only 4 pixels: (0,0), (0,1), (1,0), (1,1)
		p = src[RIDX(0, 0, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(0, 1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(1, 0, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(1, 1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
	}
	// Top-right corner (i=0, j=dim-1)
	else if (i == 0 && j == dim-1) {
		// Only 4 pixels: (0,dim-2), (0,dim-1), (1,dim-2), (1,dim-1)
		p = src[RIDX(0, dim-2, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(0, dim-1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(1, dim-2, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(1, dim-1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
	}
	// Bottom-left corner (i=dim-1, j=0)
	else if (i == dim-1 && j == 0) {
		// Only 4 pixels: (dim-2,0), (dim-2,1), (dim-1,0), (dim-1,1)
		p = src[RIDX(dim-2, 0, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(dim-2, 1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(dim-1, 0, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(dim-1, 1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
	}
	// Bottom-right corner (i=dim-1, j=dim-1)
	else if (i == dim-1 && j == dim-1) {
		// Only 4 pixels: (dim-2,dim-2), (dim-2,dim-1), (dim-1,dim-2), (dim-1,dim-1)
		p = src[RIDX(dim-2, dim-2, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(dim-2, dim-1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(dim-1, dim-2, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(dim-1, dim-1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
	}
	// Top edge (i=0, 0<j<dim-1)
	else if (i == 0) {
		// 6 pixels: (0,j-1), (0,j), (0,j+1), (1,j-1), (1,j), (1,j+1)
		p = src[RIDX(0, j-1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(0, j, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(0, j+1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(1, j-1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(1, j, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(1, j+1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
	}
	// Bottom edge (i=dim-1, 0<j<dim-1)
	else if (i == dim-1) {
		// 6 pixels: (dim-2,j-1), (dim-2,j), (dim-2,j+1), (dim-1,j-1), (dim-1,j), (dim-1,j+1)
		p = src[RIDX(dim-2, j-1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(dim-2, j, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(dim-2, j+1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(dim-1, j-1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(dim-1, j, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(dim-1, j+1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
	}
	// Left edge (0<i<dim-1, j=0)
	else if (j == 0) {
		// 6 pixels: (i-1,0), (i-1,1), (i,0), (i,1), (i+1,0), (i+1,1)
		p = src[RIDX(i-1, 0, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(i-1, 1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(i, 0, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(i, 1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(i+1, 0, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(i+1, 1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
	}
	// Right edge (0<i<dim-1, j=dim-1)
	else if (j == dim-1) {
		// 6 pixels: (i-1,dim-2), (i-1,dim-1), (i,dim-2), (i,dim-1), (i+1,dim-2), (i+1,dim-1)
		p = src[RIDX(i-1, dim-2, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(i-1, dim-1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(i, dim-2, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(i, dim-1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(i+1, dim-2, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(i+1, dim-1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
	}
	// Interior (all 9 pixels)
	else {
		p = src[RIDX(i-1, j-1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(i-1, j, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(i-1, j+1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(i, j-1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(i, j, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(i, j+1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(i+1, j-1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(i+1, j, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
		
		p = src[RIDX(i+1, j+1, dim)];
		sum.red += (int) p.red; sum.green += (int) p.green; sum.blue += (int) p.blue; sum.num++;
	}

	// assign_sum_to_pixel(&current_pixel, sum);

	//function inlining
	current_pixel.red = (unsigned short) (sum.red/sum.num);
	current_pixel.green = (unsigned short) (sum.green/sum.num);
	current_pixel.blue = (unsigned short) (sum.blue/sum.num);

	return current_pixel;
}

/******************************************************
 * Your different versions of the smooth kernel go here
 ******************************************************/

/*
 * naive_smooth - The naive baseline version of smooth 
 */
  /* The parameters, pointers, rusage_time, rdtsc_time, and cpu_time_used are used to measure performance and return values to caller. */
 /* You should not change the code that uses these parameters and variables. */
 
char naive_smooth_descr[] = "naive_smooth: Naive baseline implementation";
void naive_smooth(int dim, pixel *src, pixel *dst, int *rusage_time, unsigned long long *rdtsc_time) 
{
	int i, j;
	
	/* the variables below are used for performance measurement and not for computing the results of the algorithm */
	long int rusage_start_time, rusage_end_time = 0;
        unsigned long long rdtsc_start_time, rdtsc_end_time = 0;

	/* call system functions to start measuring performance. you should not bother to change these. */
        rusage_start_time = getUserTime();
	rdtsc_start_time = rdtsc();

/* below are the main computations for the smooth function */
	for (j = 0; j < dim; j++)
		for (i = 0; i < dim; i++)
			dst[RIDX(i, j, dim)] = avg(dim, i, j, src);

/* the remaining lines in this function stop the measurement and set the values before returning. */
	rusage_end_time = getUserTime();
        rdtsc_end_time = rdtsc();

	*rusage_time = rusage_end_time - rusage_start_time;
	*rdtsc_time = rdtsc_end_time - rdtsc_start_time;
}

 /* The parameters, pointers, rusage_time, rdtsc_time, and cpu_time_used are used to measure performance and return values to caller. */
 /* You should not change the code that uses these parameters and variables. */
 
char my_smooth_descr[] = "my_smooth: Naive baseline implementation";
void my_smooth(int dim, pixel *src, pixel *dst, int *rusage_time, unsigned long long *rdtsc_time) 
{
	int i, j;
	
	/* the variables below are used for performance measurement and not for computing the results of the algorithm */
	long int rusage_start_time, rusage_end_time = 0;
        unsigned long long rdtsc_start_time, rdtsc_end_time = 0;
	/* call system functions to start measuring performance. you should not bother to change these. */
        rusage_start_time = getUserTime();
	rdtsc_start_time = rdtsc();
	
/* ANY CHANGES TO BE MADE SHOULD BE BELOW HERE */
/* below are the main computations for your implementation of the smooth function. Any changes in implementation will go here or the other functiosn it calls */

	// for (j = 0; j < dim; j++)
	// 	for (i = 0; i < dim; i++)
	// 		dst[RIDX(i, j, dim)] = avg(dim, i, j, src);

	//loop interchange
	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++)
			dst[RIDX(i, j, dim)] = my_avg(dim, i, j, src);

	//loop interchange + loop unrolling
	// for (i = 0; i < dim; i++)
	// 	for (j = 0; j < dim; j+=8){
	// 		dst[RIDX(i, j, dim)] = avg(dim, i, j, src);
	// 		if(j+1<dim) dst[RIDX(i, j+1, dim)] = avg(dim, i, j+1, src);
	// 		if(j+2<dim) dst[RIDX(i, j+2, dim)] = avg(dim, i, j+2, src);
	// 		if(j+3<dim) dst[RIDX(i, j+3, dim)] = avg(dim, i, j+3, src);
	// 		if(j+4<dim) dst[RIDX(i, j+4, dim)] = avg(dim, i, j+4, src);
	// 		if(j+5<dim) dst[RIDX(i, j+5, dim)] = avg(dim, i, j+5, src);
	// 		if(j+6<dim) dst[RIDX(i, j+6, dim)] = avg(dim, i, j+6, src);
	// 		if(j+7<dim) dst[RIDX(i, j+7, dim)] = avg(dim, i, j+7, src);
	// 	}

	//loop interchange + loop unrolling + remove RIDX
	// for (i = 0; i < dim; i++)
	// 	for (j = 0; j < dim; j+=8){
	// 		dst[(i)*(dim)+(j)] = avg(dim, i, j, src);
	// 		if(j+1<dim) dst[(i)*(dim)+(j+1)] = avg(dim, i, j+1, src);
	// 		if(j+2<dim) dst[(i)*(dim)+(j+2)] = avg(dim, i, j+2, src);
	// 		if(j+3<dim) dst[(i)*(dim)+(j+3)] = avg(dim, i, j+3, src);
	// 		if(j+4<dim) dst[(i)*(dim)+(j+4)] = avg(dim, i, j+4, src);
	// 		if(j+5<dim) dst[(i)*(dim)+(j+5)] = avg(dim, i, j+5, src);
	// 		if(j+6<dim) dst[(i)*(dim)+(j+6)] = avg(dim, i, j+6, src);
	// 		if(j+7<dim) dst[(i)*(dim)+(j+7)] = avg(dim, i, j+7, src);
	// 	}

	//loop interchange + loop unrolling + blocking(small block size)
	// for (i = 0; i < dim; i+=BLOCK_SMOOTH)
	// 	for (j = 0; j < dim; j+=BLOCK_SMOOTH)
	// 		for(int ii=i; ii<i+BLOCK_SMOOTH && ii<dim;ii++)
	// 			for(int jj=j;jj<j+BLOCK_SMOOTH && jj<dim;jj+=8){
	// 				dst[RIDX(ii, jj, dim)] = avg(dim, ii, jj, src);
	// 				if(jj+1<dim) dst[RIDX(ii, jj+1, dim)] = avg(dim, ii, jj+1, src);
	// 				if(jj+2<dim) dst[RIDX(ii, jj+2, dim)] = avg(dim, ii, jj+2, src);
	// 				if(jj+3<dim) dst[RIDX(ii, jj+3, dim)] = avg(dim, ii, jj+3, src);
	// 				if(jj+4<dim) dst[RIDX(ii, jj+4, dim)] = avg(dim, ii, jj+4, src);
	// 				if(jj+5<dim) dst[RIDX(ii, jj+5, dim)] = avg(dim, ii, jj+5, src);
	// 				if(jj+6<dim) dst[RIDX(ii, jj+6, dim)] = avg(dim, ii, jj+6, src);
	// 				if(jj+7<dim) dst[RIDX(ii, jj+7, dim)] = avg(dim, ii, jj+7, src);
	// 			}
		




/* end of computation for smooth function. so don't change anything after this in this function. */
/* END OF CHANGES */

/* the remaining lines in this function stop the measurement and set the values before returning. */
	rusage_end_time = getUserTime();
        rdtsc_end_time = rdtsc();

	*rusage_time = rusage_end_time - rusage_start_time;
	*rdtsc_time = rdtsc_end_time - rdtsc_start_time;
}
