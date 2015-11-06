#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <complex.h>
#define PI 3.14159265358979323846
#define debug 0
#define SIZE_MAX ((size_t)-1)
typedef double complex cplx;
void _fft(cplx buf[], cplx out[], int n, int step); 
void fft(cplx buf[], int n);
void show(const char * s, cplx buf[],int N);
int isPowerOfTwo(unsigned int x);
int check(int rc);
double log_2(double n);
void swapPoint(cplx* a,cplx*b);




// Private function prototypes
static size_t reverse_bits(size_t x, unsigned int n);
static void *memdup(const void *src, size_t n);




int transform(double real[], double imag[], size_t n) {
	if (n == 0)
		return 1;
	else if ((n & (n - 1)) == 0)  // Is power of 2
		return transform_radix2(real, imag, n);
	else  // More complicated algorithm for arbitrary sizes
		return transform_bluestein(real, imag, n);
}


int inverse_transform(double real[], double imag[], size_t n) {
	return transform(imag, real, n);
}


int transform_radix2(double real[], double imag[], size_t n) {
	// Variables
	int status = 0;
	unsigned int levels;
	double *cos_table, *sin_table;
	size_t size;
	size_t i;
	
	// Compute levels = floor(log2(n))
	{
		size_t temp = n;
		levels = 0;
		while (temp > 1) {
			levels++;
			temp >>= 1;
		}
		if (1u << levels != n)
			return 0;  // n is not a power of 2
	}
	
	// Trignometric tables
	if (SIZE_MAX / sizeof(double) < n / 2)
		return 0;
	size = (n / 2) * sizeof(double);
	cos_table = malloc(size);
	sin_table = malloc(size);
	if (cos_table == NULL || sin_table == NULL)
		goto cleanup;
	for (i = 0; i < n / 2; i++) {
		cos_table[i] = cos(2 * M_PI * i / n);
		sin_table[i] = sin(2 * M_PI * i / n);
	}
	
	// Bit-reversed addressing permutation
	for (i = 0; i < n; i++) {
		size_t j = reverse_bits(i, levels);
		if (j > i) {
			double temp = real[i];
			real[i] = real[j];
			real[j] = temp;
			temp = imag[i];
			imag[i] = imag[j];
			imag[j] = temp;
		}
	}
	
	// Cooley-Tukey decimation-in-time radix-2 FFT
	for (size = 2; size <= n; size *= 2) {
		size_t halfsize = size / 2;
		size_t tablestep = n / size;
		for (i = 0; i < n; i += size) {
			size_t j;
			size_t k;
			for (j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
				double tpre =  real[j+halfsize] * cos_table[k] + imag[j+halfsize] * sin_table[k];
				double tpim = -real[j+halfsize] * sin_table[k] + imag[j+halfsize] * cos_table[k];
				real[j + halfsize] = real[j] - tpre;
				imag[j + halfsize] = imag[j] - tpim;
				real[j] += tpre;
				imag[j] += tpim;
			}
		}
		if (size == n)  // Prevent overflow in 'size *= 2'
			break;
	}
	status = 1;
	
cleanup:
	free(cos_table);
	free(sin_table);
	return status;
}


int transform_bluestein(double real[], double imag[], size_t n) {
	// Variables
	int status = 0;
	double *cos_table, *sin_table;
	double *areal, *aimag;
	double *breal, *bimag;
	double *creal, *cimag;
	size_t m;
	size_t size_n, size_m;
	size_t i;
	
	// Find a power-of-2 convolution length m such that m >= n * 2 + 1
	{
		size_t target;
		if (n > (SIZE_MAX - 1) / 2)
			return 0;
		target = n * 2 + 1;
		for (m = 1; m < target; m *= 2) {
			if (SIZE_MAX / 2 < m)
				return 0;
		}
	}
	
	// Allocate memory
	if (SIZE_MAX / sizeof(double) < n || SIZE_MAX / sizeof(double) < m)
		return 0;
	size_n = n * sizeof(double);
	size_m = m * sizeof(double);
	cos_table = malloc(size_n);
	sin_table = malloc(size_n);
	areal = calloc(m, sizeof(double));
	aimag = calloc(m, sizeof(double));
	breal = calloc(m, sizeof(double));
	bimag = calloc(m, sizeof(double));
	creal = malloc(size_m);
	cimag = malloc(size_m);
	if (cos_table == NULL || sin_table == NULL
			|| areal == NULL || aimag == NULL
			|| breal == NULL || bimag == NULL
			|| creal == NULL || cimag == NULL)
		goto cleanup;
	
	// Trignometric tables
	for (i = 0; i < n; i++) {
		double temp = M_PI * (size_t)((unsigned long long)i * i % ((unsigned long long)n * 2)) / n;
		// Less accurate version if long long is unavailable: double temp = M_PI * i * i / n;
		cos_table[i] = cos(temp);
		sin_table[i] = sin(temp);
	}
	
	// Temporary vectors and preprocessing
	for (i = 0; i < n; i++) {
		areal[i] =  real[i] * cos_table[i] + imag[i] * sin_table[i];
		aimag[i] = -real[i] * sin_table[i] + imag[i] * cos_table[i];
	}
	breal[0] = cos_table[0];
	bimag[0] = sin_table[0];
	for (i = 1; i < n; i++) {
		breal[i] = breal[m - i] = cos_table[i];
		bimag[i] = bimag[m - i] = sin_table[i];
	}
	
	// Convolution
	if (!convolve_complex(areal, aimag, breal, bimag, creal, cimag, m))
		goto cleanup;
	
	// Postprocessing
	for (i = 0; i < n; i++) {
		real[i] =  creal[i] * cos_table[i] + cimag[i] * sin_table[i];
		imag[i] = -creal[i] * sin_table[i] + cimag[i] * cos_table[i];
	}
	status = 1;
	
	// Deallocation
cleanup:
	free(cimag);
	free(creal);
	free(bimag);
	free(breal);
	free(aimag);
	free(areal);
	free(sin_table);
	free(cos_table);
	return status;
}


int convolve_real(const double x[], const double y[], double out[], size_t n) {
	double *ximag, *yimag, *zimag;
	int status = 0;
	ximag = calloc(n, sizeof(double));
	yimag = calloc(n, sizeof(double));
	zimag = calloc(n, sizeof(double));
	if (ximag == NULL || yimag == NULL || zimag == NULL)
		goto cleanup;
	
	status = convolve_complex(x, ximag, y, yimag, out, zimag, n);
cleanup:
	free(zimag);
	free(yimag);
	free(ximag);
	return status;
}


int convolve_complex(const double xreal[], const double ximag[], const double yreal[], const double yimag[], double outreal[], double outimag[], size_t n) {
	int status = 0;
	size_t size;
	size_t i;
	double *xr, *xi, *yr, *yi;
	if (SIZE_MAX / sizeof(double) < n)
		return 0;
	size = n * sizeof(double);
	xr = memdup(xreal, size);
	xi = memdup(ximag, size);
	yr = memdup(yreal, size);
	yi = memdup(yimag, size);
	if (xr == NULL || xi == NULL || yr == NULL || yi == NULL)
		goto cleanup;
	
	if (!transform(xr, xi, n))
		goto cleanup;
	if (!transform(yr, yi, n))
		goto cleanup;
	for (i = 0; i < n; i++) {
		double temp = xr[i] * yr[i] - xi[i] * yi[i];
		xi[i] = xi[i] * yr[i] + xr[i] * yi[i];
		xr[i] = temp;
	}
	if (!inverse_transform(xr, xi, n))
		goto cleanup;
	for (i = 0; i < n; i++) {  // Scaling (because this FFT implementation omits it)
		outreal[i] = xr[i] / n;
		outimag[i] = xi[i] / n;
	}
	status = 1;
	
cleanup:
	free(yi);
	free(yr);
	free(xi);
	free(xr);
	return status;
}


static size_t reverse_bits(size_t x, unsigned int n) {
	size_t result = 0;
	unsigned int i;
	for (i = 0; i < n; i++, x >>= 1)
		result = (result << 1) | (x & 1);
	return result;
}


static void *memdup(const void *src, size_t n) {
	void *dest = malloc(n);
	if (dest != NULL)
		memcpy(dest, src, n);
	return dest;
}


















void swapPoint(cplx* a,cplx*b){
    cplx temp = *a;
    *a = *b;
    *b = temp;
}

double log_2(double n){
return log(n)/log(2);
}

int isPowerOfTwo(unsigned int x){
    return((x!=0) && ((x&(~x+1))==x));
}

void _fft(cplx buf[], cplx out[], int n, int step)
{    
    
	if (step < n) {
		_fft(out, buf, n, step * 2);
		_fft(out + step, buf + step, n, step * 2);
        int i;
		for (i = 0; i < n; i += 2 * step) {
			cplx t = cexp(-I * PI * i / n) * out[i + step];
			buf[i / 2]     = out[i] + t;
			buf[(i + n)/2] = out[i] - t;

		}
	}
}
 
void fft(cplx buf[], int n)
{
	//cplx out[n];
    
    cplx *out;
    out = (cplx *) malloc(sizeof(cplx)*n);
    if(out==NULL){fprintf(stderr,"Cannot allocate");exit(2);}

    int i;
	for (i = 0; i < n; i++) out[i] = buf[i];
 
	_fft(buf, out, n, 1);
    free(out);
}
 
 
void show(const char * s, cplx buf[],int N) {
	printf("%s", s);
    int i;
	for (i = 0; i < N; i++)
		if (!cimag(buf[i]))
			printf("%g ", creal(buf[i]));
		else
			printf("(%g, %g) ", creal(buf[i]), cimag(buf[i]));
    printf("\n");
}
 
int main(int argc, char **argv)
{
    int rank, size, l, rc;
    rc = MPI_Init(&argc, &argv);
    MPI_Status statSend,statRecv;
    MPI_Request reqSend, reqRecv;
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &size);  
    check(rc);
    if (isPowerOfTwo((unsigned int) size) == 0 && size!=1){
        fprintf(stderr,"Proc not power of 2");
        exit(2);
    }
    if (argc != 2 && argc !=3){
        fprintf(stderr,"Wrong input N is, 2^N inputs per proc, ./Filename N");
        exit(2);    
    }
    
    clock_t begin, end;
    double time_spent;
    if(rank == 0){    
    begin = clock();}
    
	int elementsPowProc = atoi(argv[1]);
    int addN = 0;
    if (argc ==3){
        addN = atoi(argv[2]);
        printf("checkPoint\n");
    }
    int N = (int) pow(2,elementsPowProc);
    printf("n value before %d\n",N);
    N = N + addN;
    printf("n value after %d\n",N);
    //printf("every calculate N %d \n",N);
    //printf("Malloc size  %f \n",(double) sizeof(cplx)*N);   
    cplx *fCoefficient;
    cplx *recvbuffer;
    cplx *temp;
    fCoefficient = (cplx *) malloc(sizeof(cplx)*N);
    recvbuffer   = (cplx *) malloc(sizeof(cplx)*N);    
    if(fCoefficient == NULL || recvbuffer == NULL){
        fprintf(stderr,"cannot alloc memory");
        exit(2);
        };     
    if(debug) printf("allocated memory done");
    srandom(rank+1);
    
    //for (l = 0; l < N;l++){
    //    fCoefficient[l] = ((cplx) random()/(RAND_MAX) -1/2)*100;
    //}
    //if(debug) printf("begin with seq FFT");
    //fft(fCoefficient, N);//Sequencial FFT
    
    double *refoutreal, *refoutimag;
    double valueR, valueI;
    refoutreal = malloc(N * sizeof(double));
	refoutimag = malloc(N * sizeof(double));
    for (l = 0; l < N;l++){
        refoutreal[l] = ((double) random()/(RAND_MAX) -1/2)*100;
        refoutimag[l] = 0;
    }
    transform(refoutreal, refoutimag, N);
    for (l = 0; l < N;l++){
        valueR = refoutreal[l];
        valueI = refoutimag[l];
        fCoefficient[l] = valueR + valueI * I;
    }
    free(refoutreal);
	free(refoutimag);

    if(debug) printf("done with seq fft");
    if (debug == 3){
        show("\nFFT : ", fCoefficient,N);
    }
    int numOfComm = (int)(log_2(size)-0.00001);
    int target = 0;
    int bit,tag,j,n,k;
    target ^= 1 << 4;
    //printf("target new number %d\n",target);
    if(debug) printf("beginning with par FFT");
    for(l=0;numOfComm-l>0;l++){
        target = rank;    
        target ^=1<< l;
        //target ^=1<< numOfComm-1-l; //bitfliped target
        bit = (rank >> numOfComm-1-l)&1; //Determine who sends first
        if(0 == bit){
            rc = MPI_Isend(fCoefficient,N,MPI_DOUBLE_COMPLEX,target,tag,MPI_COMM_WORLD,&reqSend);
            check(rc);
            rc = MPI_Irecv(recvbuffer,N,MPI_DOUBLE_COMPLEX,target,tag,MPI_COMM_WORLD,&reqRecv);
            check(rc);
        }
        else{//BIT = 1 IS ODD SHOULD MULTIPLY WITH TWINDLE FACTOR
            rc = MPI_Irecv(recvbuffer,N,MPI_DOUBLE_COMPLEX,target,tag,MPI_COMM_WORLD,&reqRecv);
            check(rc);
            //int newVectSize = N*size/(numOfComm-l);
            int newVectSize = N*pow(2,l+1);    
            for(n=0;n<N;n++){
                k = n + rank%pow(2,l);
                fCoefficient[n] = fCoefficient[n]*cexp(-2*PI*I*k/(newVectSize));
            }
          
            rc = MPI_Isend(fCoefficient,N,MPI_DOUBLE_COMPLEX,target,tag,MPI_COMM_WORLD,&reqSend);
            check(rc);            
        }

        rc = MPI_Wait(&reqRecv,&statRecv);
        if(bit){//SUM NEGATIVELY
            //swapPoint(recvbuffer,fCoefficient);
            temp = recvbuffer;
            recvbuffer = fCoefficient;
            fCoefficient = temp;
            for(j=0;j<N;j++){
                fCoefficient[j] -= recvbuffer[j];
            }
        }
        else{ //BIT 0  POSITIVE SUM
            for(j=0;j<N;j++){
                fCoefficient[j] += recvbuffer[j];
            }
        }        
        //MAKE CORRECT SUMMATION
        rc = MPI_Wait(&reqSend,&statSend);
        
    }
        
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(rank ==0){    
        end = clock();
        time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("\n TIME SPENT %f np value %d and input value %d \n",time_spent,size,elementsPowProc);
        if(argc==3)
            printf("With added data %d\n",addN);
        
    }
    rc = MPI_Finalize();
    return 0;
}

int check(int rc)
{
    if(rc == MPI_SUCCESS)
    {
        return 0;
    }
    else
    {
        fprintf(stderr,"MPI failure code: %d",rc);
        exit(2);
    }
}
