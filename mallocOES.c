#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <time.h> 
#define PI 3.14159265359

int check(int rc);
int cmpfunc(const void * a, const void * b);
// input:  N,Nsteps,
int main(int argc, char **argv)
{
    //int n; Initializing MPI
    int rc, N,P,myrank, i;
    clock_t tic = 0;
    rc = MPI_Init(&argc,&argv);
    check(rc);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    check(rc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (argc<2){printf("Missing number of elements");}
    MPI_Status status;
    if(myrank==0){
        tic = clock();}
    //Divides the element 
    N = atoi(argv[1]);
    int NofUneven = N%P;
    double subL = N/P;
    int iSL = floor(subL);
    //Keep track of number of elements at the ranks neighbour
    int prevLength;
    int nextLength;
    int subLength;
    if(myrank < NofUneven){
        subLength = iSL+1;}
    else{subLength = iSL;
        }
    if(myrank < NofUneven+1){
        prevLength = iSL+1;}
    else{
        prevLength = iSL;}
    if(myrank+2>NofUneven){
        nextLength = iSL;}
    else{
        nextLength = iSL+1;}

    srandom(myrank+1);
    //Allocate the memory
    double *xList;
    xList = (double *)malloc(sizeof(double)*iSL);
    double *tempList;
    tempList = (double *)malloc(sizeof(double)*iSL);
    double *anotherList;
    anotherList = (double *)malloc(sizeof(double)*(1+iSL));
    //Makes sure we get the allocated memory otherwise exit
    if ((xList == NULL || tempList == NULL) || anotherList == NULL){
        return -1;}
    for (i = 0; i < subLength;i++){
        xList[i] = ((double) random())/(RAND_MAX);
    }
    //The initial sort, quicksort usualy quickest
    qsort(xList,subLength,sizeof(double),cmpfunc);
    
    int sortStepNum;
    int evenPhase;
    int evenProcess = (myrank%2)+1;
    int tag = 0;
    //Begin with Odd/even trans. Sort
    for(sortStepNum = 1; sortStepNum < P;sortStepNum++){
        MPI_Barrier(MPI_COMM_WORLD);
        evenPhase = (sortStepNum%2)+1;
        anotherList[iSL+1] = 2;//as all values is limited from 0 to 1
        if (evenPhase==evenProcess){//Processors belonging to right phase sends to next first then recieves
            if (myrank != P-1){ 
                rc = MPI_Send(xList,subLength,MPI_DOUBLE, myrank+1, tag, MPI_COMM_WORLD);
                check(rc);
                rc = MPI_Recv(anotherList,nextLength,MPI_DOUBLE,myrank+1,tag,MPI_COMM_WORLD,&status);
                check(rc);
                int xLcout = 0;
                int anLcou = 0;
                for(i = 0;i<subLength;i++){ //Only partial merge sort the first lowest values
                    if ((xList[xLcout] < anotherList[anLcou])|| anLcou ==nextLength){
                        tempList[i] = xList[xLcout];
                        xLcout++;}
                    else{
                        tempList[i] = anotherList[anLcou];
                        anLcou++;}
                }
                for (i = 0;i<subLength;i++){
                xList[i] = tempList[i];}                
             }
        }
        else{
            if(myrank!=0){ //Processors belonging to 'wrong' phase recieves first then sends
                rc = MPI_Recv(anotherList,prevLength,MPI_DOUBLE,myrank-1,tag,MPI_COMM_WORLD,&status);
                check(rc);
                rc = MPI_Send(xList,subLength,MPI_DOUBLE,myrank-1,tag,MPI_COMM_WORLD);
                check(rc);
                int xLcout = subLength-1;//Begin backwards to find bigest value
                int anLcou = prevLength-1;
                for(i = subLength-1;i>-1;i--){ //Partial mergesort the biggest relevant onces
                    if((xList[xLcout] > anotherList[anLcou])||anLcou == -1){ //Or statement to handle when list has not equal length
                        tempList[i] = xList[xLcout];
                        xLcout--;}
                    else{
                        tempList[i] = anotherList[anLcou];
                        anLcou--;}
                }
                for (i = 0;i<subLength;i++){
                xList[i] = tempList[i];}
             }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //printf("finsihing computation free memory"); //Control that the end list is actually sorted
    //OUTPUT THE VALUE TO CONTROL SORTING IS WORKING
    //for( n = 0 ; n < subLength; n++ ) {
    //   printf("\nvalue %f and pos %d for processor %d\n", xList[n],n,myrank);  
    //}
    

    //Free up the memory and check the elapsed time
    xList = NULL;
    free(xList);
    anotherList = NULL;
    free(anotherList);
    tempList = NULL;
    free(tempList);
    if(myrank==0){
        clock_t toc = clock();
        printf("Elapsed time with P %d processors, N %d elements, is %f seconds time\n",P,N,(double)(toc-tic)/CLOCKS_PER_SEC);
    }
    MPI_Finalize();
    return 0;

}


//Compare function for quicksort
int cmpfunc(const void * a, const void * b)
{
    if(*(double*)a > *(double*)b)
        return 1;
    else if (*(double*)a < *(double*)b)
        return -1;
    else
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
