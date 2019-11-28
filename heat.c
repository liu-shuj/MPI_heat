#define ENABLE_PRINT 0

#include "mpi.h"
#include<stdio.h>
#include<stdlib.h>
#include<float.h>
#include<math.h>

// All indexes start from 0
int ind2row(int ind, int rows, int cols)
{
    return ind/cols;
}

int ind2col(int ind, int rows, int cols)
{
    return ind-(ind/cols)*cols;
}

int rc2ind(int r, int c, int rows, int cols)
{
    return r*cols+c;
}

double** alloc_matrix(int rows, int cols)
{
    double** matrix;
    matrix = (double**) malloc(rows * sizeof(double *));
    matrix[0] = (double*) malloc(rows * cols * sizeof(double));
    for (int i = 1; i < rows; i++)
        matrix[i] = matrix[0] + i*cols;
    return matrix;
}

double maxdiff(double** m1, double** m2, int rows, int cols)
{
    double max_val=DBL_MIN;
    for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++)
            if(fabs(m1[i][j]-m2[i][j])>max_val)
                max_val=fabs(m1[i][j]-m2[i][j]);
    return max_val;
}

int main(int argc, char **argv)
{
    int np,myrank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&np);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Request sendl,sendr,sendu,sendd,recvl,recvr,recvu,recvd;
    MPI_Status statsl,statsr,statsu,statsd,statrl,statrr,statru,statrd;
    int rows,cols,r,c,vnp,hnp;
    double eps,temp;
    if(argc!=9)
    {
	if(myrank==0){
            printf("Usage: heat [rows] [cols] [row of heat source] [col of heat source] [temp of heat source] [vertical cores] [horizontal cores] [EPS]\n");
            printf("Using default settings: 10000*10000, (4,4), 50 degrees, 2*2, 1e-3...\n");
	}
        rows=10000;
        cols=10000;
        r=4;
        c=4;
        temp=50.0;
        vnp=2;
        hnp=2;
        eps=1e-3;
    }
    else
    {
        rows=atoi(argv[1]);
        cols=atoi(argv[2]);
        r=atoi(argv[3])-1;
        c=atoi(argv[4])-1;
        if((r<0) || (r>=rows) || (c<0) || (c>=cols))
        {
            printf("Error: Invalid heat source.\n");
            exit(1);
        }
        temp=atof(argv[5]);
        vnp=atoi(argv[6]);
        hnp=atoi(argv[7]);
        eps=atof(argv[8]);
        if(myrank==0)
            printf("Starting with:\n %d*%d grid,\n heat source at (%d,%d),\n temp=%f,\n splitted into %d*%d cores,\n EPS=%lf...\n", rows, cols, r, c, temp, vnp, hnp, eps);
    }
    if(myrank==0 && vnp*hnp!=np)
    {
        printf("Error: Num of cores not equal to number of MPI processes.\n");
        exit(1);
    }
 
    int myvrank=ind2row(myrank,vnp,hnp);
    int myhrank=ind2col(myrank,vnp,hnp);
    
    // distribute workload
    int myrows=rows/vnp;
    if(myvrank==vnp-1)
        myrows+=(rows%vnp);
    int mycols=cols/hnp;
    if(myhrank==hnp-1)
        mycols+=(cols%hnp);
    
    // calculate the block which heat source is in and its local coordinates
    int srcvrank=r/(rows/vnp);
    if(srcvrank==vnp)
        srcvrank-=1;
    int srchrank=c/(cols/hnp);
    if(srchrank==hnp)
        srchrank-=1;
    int srcrank=rc2ind(srcvrank,srchrank,vnp,hnp);
    int localr=r-((rows/vnp)*srcvrank);
    int localc=c-((cols/vnp)*srchrank);    
    
    // allocation and initialization
    double **my_old=alloc_matrix(myrows,mycols);
    double **my_new=alloc_matrix(myrows,mycols);
    
    for(int i=0;i<myrows;i++)
        for(int j=0;j<mycols;j++)
        {
            my_old[i][j]=0.0;
            my_new[i][j]=0.0;
        }
    if(myrank==srcrank){
        my_old[localr][localc]=temp;
        my_new[localr][localc]=temp;
    }
    
    double* sendbufu=(double*)malloc(mycols*sizeof(double));
    double* sendbufd=(double*)malloc(mycols*sizeof(double));
    double* sendbufl=(double*)malloc(myrows*sizeof(double));
    double* sendbufr=(double*)malloc(myrows*sizeof(double));
    
    // ghost cells
    double* recvbufu=(double*)malloc(mycols*sizeof(double));
    double* recvbufd=(double*)malloc(mycols*sizeof(double));
    double* recvbufl=(double*)malloc(myrows*sizeof(double));
    double* recvbufr=(double*)malloc(myrows*sizeof(double));
    
    while(1)
    {
        if(myvrank!=0) // non-blocking communication with block above
        {
            int urank=rc2ind(myvrank-1,myhrank,vnp,hnp);
            for(int i=0;i<mycols;i++)
                sendbufu[i]=my_old[0][i];
            MPI_Isend(sendbufu,mycols,MPI_DOUBLE,urank,123,MPI_COMM_WORLD,&sendu);
            MPI_Irecv(recvbufu,mycols,MPI_DOUBLE,urank,345,MPI_COMM_WORLD,&recvu);
        }
        if(myvrank!=vnp-1) //below
        {
            int drank=rc2ind(myvrank+1,myhrank,vnp,hnp);
            for(int i=0;i<mycols;i++)
                sendbufd[i]=my_old[myrows-1][i];
            MPI_Isend(sendbufd,mycols,MPI_DOUBLE,drank,345,MPI_COMM_WORLD,&sendd);
            MPI_Irecv(recvbufd,mycols,MPI_DOUBLE,drank,123,MPI_COMM_WORLD,&recvd);
        }
        if(myhrank!=0) //left
        {
            int lrank=rc2ind(myvrank,myhrank-1,vnp,hnp);
            for(int i=0;i<myrows;i++)
                sendbufl[i]=my_old[i][0];
            MPI_Isend(sendbufl,myrows,MPI_DOUBLE,lrank,567,MPI_COMM_WORLD,&sendl);
            MPI_Irecv(recvbufl,myrows,MPI_DOUBLE,lrank,789,MPI_COMM_WORLD,&recvl);
        }
        if(myhrank!=hnp-1) //right
        {
            int rrank=rc2ind(myvrank,myhrank+1,vnp,hnp);
            for(int i=0;i<myrows;i++)
                sendbufr[i]=my_old[i][mycols-1];
            MPI_Isend(sendbufr,myrows,MPI_DOUBLE,rrank,789,MPI_COMM_WORLD,&sendr);
            MPI_Irecv(recvbufr,myrows,MPI_DOUBLE,rrank,567,MPI_COMM_WORLD,&recvr);
        }
        
        // update interiors
        for(int i=1;i<myrows-1;i++)
            for(int j=1;j<mycols-1;j++)
                my_new[i][j]=0.25*(my_old[i-1][j]+my_old[i+1][j]+my_old[i][j-1]+my_old[i][j+1]);            
        if(myrank==srcrank){
            my_old[localr][localc]=temp;
            my_new[localr][localc]=temp;
        }
        
        //printf("#%d waiting for communication from above...\n",myrank); 
        if(myvrank!=0)
        {
            MPI_Wait(&sendu,&statsu);
            MPI_Wait(&recvu,&statru);
        }
        //printf("#%d waiting for communication from below...\n",myrank); 
        if(myvrank!=vnp-1)
        {
            MPI_Wait(&sendd,&statsd);
            MPI_Wait(&recvd,&statrd);
        }
        //printf("#%d waiting for communication from left...\n",myrank); 
        if(myhrank!=0)
        {
            MPI_Wait(&sendl,&statsl);
            MPI_Wait(&recvl,&statrl);
        }
        //printf("#%d waiting for communication from right...\n",myrank); 
        if(myhrank!=hnp-1)
        {
            MPI_Wait(&sendr,&statsr);
            MPI_Wait(&recvr,&statrr);
        }
        
        int lastc=mycols-1;
        int lastr=myrows-1;
        
        if(myvrank!=0)  // fill in borders: up
        {
            if(myhrank!=0)
                my_new[0][0]=0.25*(recvbufu[0]+recvbufl[0]+my_old[0][1]+my_old[1][0]);
            if(myhrank!=hnp-1)
                my_new[0][lastc]=0.25*(recvbufu[lastc]+recvbufr[0]+my_old[0][lastc-1]+my_old[1][lastc]);
            for(int i=1;i<lastc;i++)
                my_new[0][i]=0.25*(recvbufu[i]+my_old[0][i-1]+my_old[0][i+1]+my_old[1][i]);
        }
        if(myvrank!=vnp-1)  //down
        {
            if(myhrank!=0)
                my_new[lastr][0]=0.25*(recvbufd[0]+recvbufl[lastr]+my_old[lastr-1][0]+my_old[lastr][1]);
            if(myhrank!=hnp-1)
                my_new[lastr][lastc]=0.25*(recvbufd[lastc]+recvbufr[lastr]+my_old[lastr][lastc-1]+my_old[lastr-1][lastc]);
            for(int i=1;i<lastc;i++)
                my_new[lastr][i]=0.25*(recvbufd[i]+my_old[lastr][i-1]+my_old[lastr][i+1]+my_old[lastr-1][i]);
        }
        if(myhrank!=0)  //left
        {
            for(int i=1;i<lastr;i++)
                my_new[i][0]=0.25*(recvbufl[i]+my_old[i-1][0]+my_old[i+1][0]+my_old[i][1]);
        }
        if(myhrank!=hnp-1)  //right
        {
            for(int i=1;i<lastr;i++)
                my_new[i][lastc]=0.25*(recvbufr[i]+my_old[i-1][lastc]+my_old[i+1][lastc]+my_old[i][lastc-1]);
        }
        //if(myrank==0){
            //for(int i=0;i<myrows;i++){
        //	for(int j=0;j<mycols;j++){
        //		printf("%lf ",my_new[i][j]);
        //	}
        //	printf("\n");
        //} 
        //}
        double diff=maxdiff(my_old,my_new,myrows,mycols);
        double cur_maxdiff;
        // terminating criteria
        MPI_Allreduce(&diff,&cur_maxdiff,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        //if(myrank==0)
            //    printf("%lf\n",cur_maxdiff);
        if(cur_maxdiff<eps)
            break;
        
        double **c=NULL;
        c=my_new;
        my_new=my_old;
        my_old=c;
    }
    
    if(ENABLE_PRINT){
        FILE* f=fopen("heat_dist","w");
        fclose(f);
        for(int j=0;j<vnp;j++) 
        {
            int nlines=rows/vnp;
            if(j==vnp-1)
                nlines+=rows%vnp;
            for(int l=0;l<nlines;l++) 
            {
                for(int i=0;i<hnp;i++) 
                {
                    if(myrank==(rc2ind(j,i,vnp,hnp)))
                    {
                        f=fopen("heat_dist","a");
                        for(int k=0;k<mycols;k++)
                            fprintf(f,"%lf ",my_new[l][k]);
                        if(i==hnp-1){
                            printf("Printing line %d...\n",j*(rows/vnp)+l);
                            fprintf(f,"\n"); 
                        }
                        fclose(f);
                    }
                    MPI_Barrier(MPI_COMM_WORLD);
                }
            }
        }
        printf("Results saved to 'heat_dist'...\n");
    }
    MPI_Finalize();
    return 0;
}
        
        
        
        
