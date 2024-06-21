/*
 OBJECTIVE:
 
Do BD simulation for Rouse model under the uniaxial tension 

!!! with the correction output in this version   
  
% The uniaxial tension flow is applied in the z direction 
% chainIdx: Index of chain (simulation)
% chainLength: Number of beads per chain (N)
% delta_t: Time step
% rate: the applied Hencky strain rate in the z direction 
% startStep: starting step of the simulation 
% endStep: ending step of the simulation 
% startOfDump: Start of dump
% dumpPeriod: Number of steps between dumps

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mpi.h>

double myrandom();
double normalrandom();
void Rouse_model_BD_simulation(int chainIdx, int chainLength, float delta_t, double rate, int startStep, int endStep, int startOfDump,int dumpPeriod); 


int main(int argc, char** argv){
		
//====mpi variables
    int me=0,nprocs=1;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    int chain_per_mpi_thread = 0; //Number of chain per mpi thread
    int max_id = 0; //the cpu ids that have one more chain than the chain per mpi thread 
    
    clock_t start, end; //timing variables
    if (me==0) start = clock();
	
    int numberOfChains,chainLength;
    int startStep,endStep,startOfDump,dumpPeriod ; 
    float  delta_t; 
    double rate;
     
    char buffer[256],TempC[30];


// Reading the in put file 
if(me==0){    
    while(1){
        fgets(buffer,254,stdin);
        if(feof(stdin)!=0){
            break;
         }
        else{
            if(strstr(buffer,"NumberChains")){
                sscanf(buffer, "%s %i", &TempC,&numberOfChains);
                printf("# %s %i\n", TempC,numberOfChains);
            }
            if(strstr(buffer,"BeadsPerChain")){
                sscanf(buffer, "%s %i", &TempC,&chainLength);
                printf("# %s %i\n",TempC,chainLength);
            }
            if(strstr(buffer,"dt")){
                sscanf(buffer, "%s %f", &TempC,&delta_t);
                printf("# %s %f\n",TempC,delta_t);
            }
            if(strstr(buffer,"StartStep")){
                sscanf(buffer, "%s %i", &TempC,&startStep);
                printf("# %s %i\n",TempC,startStep);
            }
            if(strstr(buffer,"EndStep")){
                sscanf(buffer, "%s %i", &TempC,&endStep);
                printf("# %s %i\n",TempC,endStep);
            }
            if(strstr(buffer,"StartOfDump")){
                sscanf(buffer, "%s %i", &TempC,&startOfDump);
                printf("# %s %i\n",TempC,startOfDump);
            }
            
            if(strstr(buffer,"DumpPeriod")){
                sscanf(buffer, "%s %i", &TempC,&dumpPeriod);
                printf("# %s %i\n",TempC,dumpPeriod);
            }
            if(strstr(buffer,"ERate")){
                sscanf(buffer, "%s %lf", &TempC,&rate);
                printf("# %s %.16f\n",TempC,rate);
            }
         }
    }
    
    
    if(nprocs>numberOfChains){
        printf("!!!!\n");
        printf("ERROR: We excpect number of processor smaller than the chains!!!!\n");
        printf("!!!!\n");
        exit(EXIT_FAILURE);
    }
    chain_per_mpi_thread = (int) floorf(1.0*numberOfChains/nprocs);
    max_id=(numberOfChains%nprocs);
    
    printf("# MPI threads %i\n",nprocs);
    printf("# chains/mpi thread %i\n",chain_per_mpi_thread);
    
	
}

    MPI_Bcast(&numberOfChains,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&chainLength,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&delta_t,1,MPI_FLOAT,0,MPI_COMM_WORLD);
    MPI_Bcast(&startStep,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&endStep,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&startOfDump,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&dumpPeriod,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&rate,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&chain_per_mpi_thread,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&max_id,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
//====sperad the chain over the cpus.       
    int mychain;  
    if(me<max_id){
        mychain=chain_per_mpi_thread+1;
    }else{
        mychain=chain_per_mpi_thread; 
    } 
    
    int   *chain_index; 
    chain_index= (int *) malloc (mychain*sizeof(int));
    
    for(int i=0;i<mychain;i++)
    {
        if(me<max_id) chain_index[i]=me*mychain+i+1;  
        else chain_index[i] =me*(mychain+1)+i+1-(me-max_id);     
    }        
    MPI_Barrier(MPI_COMM_WORLD);    
    
//====output the chain spreading informatoin
    FILE *fpw0;
    sprintf(buffer,"Rank_%i_chain_info.dat",me);
    fpw0=fopen(buffer,"w");   
    fprintf(fpw0,"# rank %i: chain id: ",me);
    for(int i=0;i<mychain;i++){
    fprintf(fpw0," %i ",chain_index[i]);    
    }
    fprintf(fpw0,"#\n");
    fclose(fpw0);
    
//====start the BD simulation 
    if (me==0){
        printf("#---------------------\n");
        printf("#Start the BD simulation here \n");
    }
    
    srand((unsigned)time(NULL)+me*nprocs); /* initialize rand() */  


    for(int i=0;i<mychain;i++){
        Rouse_model_BD_simulation(chain_index[i],chainLength,delta_t,rate,startStep,endStep,startOfDump,dumpPeriod); 
        printf("# chain : %i is finished\n", chain_index[i]);        
    }


    if(me==0){
        end = clock();
        double total_time = ((double)(end-start))/CLOCKS_PER_SEC;
        printf("# total_time %f\n", total_time);
    }
    MPI_Finalize();

}

//=========================================================
void Rouse_model_BD_simulation(int chainIdx, int chainLength, float delta_t, double rate, int startStep, int endStep, int startOfDump,int dumpPeriod){
    
    
    FILE *fpw0;
    char buffer[256];
    
    float beadSize=1.0; 
    // By default b = 1.0
    // !!!! Please do not change beadSize in this version of program !!!!
    // !!!! The program is written with reduced units, for which b = 1.0 !!!!
 
//====initialize the configuration with the first bead position at (0,0,0)
    int numBonds = chainLength - 1;
    double      *configQ;
    long double *configR;
    configQ       = (double *) malloc (numBonds*3*sizeof(double));
    configR       = (long double *) malloc (chainLength*3*sizeof(long double));
    
    long double *Cumul_disp; 
    Cumul_disp    = (long double *) malloc (chainLength*3*sizeof(long double)); // store the cumulative flow induced displacemnt
    
    for(int i=0; i<numBonds; i++){
        for(int j=0; j<3;j++){
         // Generating a series of x with mean = 0 and variance b^2/3   
         configQ[i*3+j]= (beadSize/sqrt(3.0))*normalrandom();     
        }
    }
    
    if (startStep==0){  //set the inital configurations at the very beginning.   
        configR[0]=0.0;
        configR[1]=0.0;
        configR[2]=0.0; 
        for(int i=0; i<numBonds; i++){
            for(int j=0; j<3;j++){      
                configR[(i+1)*3+j]=configR[i*3+j]+configQ[i*3+j];
            }
        }
        
        for(int i=0; i<chainLength*3; i++) Cumul_disp[i]=0.0;
        
        sprintf(buffer,"./chain_%i/dump_chain_%i_%i.txt",chainIdx,chainIdx,startStep);   // out put the inital configurations 
        fpw0=fopen(buffer,"w");
        for(int i=0; i<chainLength; i++){    
            fprintf(fpw0,"%i %.16Lf %.16Lf %.16Lf %.16Lf %.16Lf %.16Lf\n",i+1, configR[i*3+0],configR[i*3+1],configR[i*3+2],Cumul_disp[i*3+0],Cumul_disp[i*3+1],Cumul_disp[i*3+2]); 
        }
        fclose(fpw0);
    } 
    else
    {   
        //restart the simulations, if required.
        sprintf(buffer,"./chain_%i/dump_chain_%i_%i.txt",chainIdx,chainIdx,startStep);        
        fpw0=fopen(buffer,"r");      
        int id;       
        for(int i=0; i<chainLength;i++)
        {
            fgets(buffer,254,fpw0);
            sscanf(buffer,"%i %Lf %Lf %Lf %Lf %Lf %Lf",&id,&configR[i*3+0],&configR[i*3+1],&configR[i*3+2], &Cumul_disp[i*3+0],&Cumul_disp[i*3+1],&Cumul_disp[i*3+2]);
        }
        fclose(fpw0);       
        
    } 
        
//====Generate the Rouse matrix    
    float *RouseMatrix;
    RouseMatrix       = (float *) malloc (chainLength*chainLength*sizeof(float));
    
    for(int i=0; i<chainLength; i++){
        for(int j=0; j<chainLength;j++){
            if(i==j) RouseMatrix[i*chainLength+j]=(-3.0/2.0)*delta_t*2.0;
            else if ( (i==j+1) || (i== j-1) ) RouseMatrix[i*chainLength+j]= (-3.0/2.0)*delta_t*(-1.0);
            else RouseMatrix[i*chainLength+j]= 0.0;
        }
    }
    
    RouseMatrix[0]=(-3.0/2.0)*delta_t*1.0;
    RouseMatrix[chainLength*chainLength-1]=(-3.0/2.0)*delta_t*1.0;
        
//====Prepare for the generation of random force
    double upperBound   = sqrt(3.0*delta_t); // ub for generating W_j
    double lowerBound   = -upperBound;       // lb for generating W_j
    double VelGrandient = rate*delta_t;      // the velocity grandient matrix 
    

    long double *increment;
    increment       = (long double*) malloc (chainLength*3*sizeof(long double));
       
//====Update the configurations     
    for(int stepIdx=startStep+1; stepIdx<=endStep; stepIdx++){
            
        // setup the incremnt matrix             
        for(int i=0; i<chainLength; i++){                    
            for(int j=0; j<3; j++){ 
                increment[i*3+j] = 0.0;
                    
                // random force                    
                increment[i*3+j] += lowerBound+(upperBound-lowerBound)*myrandom( );  
                // spring force                    
                for(int k=0; k<chainLength; k++) increment[i*3+j] += RouseMatrix[i*chainLength+k]*configR[k*3+j]; 
                // shear force                
                if (j==0 || j==1) // x, y direction 
                {
                    increment[i*3+j]     += -0.5*configR[i*3+j]*VelGrandient;
                    Cumul_disp[i*3+j]    += -0.5*configR[i*3+j]*VelGrandient;                   
                }
                else             // z direction 
                {
                    increment[i*3+j]     += configR[i*3+j]*VelGrandient;
                    Cumul_disp[i*3+j]    += configR[i*3+j]*VelGrandient;                         
                }                 
            }              
        }
            
        //update the configurations
        for(int i=0; i<chainLength; i++){
            for(int j=0; j<3; j++) configR[i*3+j]=configR[i*3+j]+increment[i*3+j];    
        }
        
        //output the configurations   
        if(stepIdx >= startOfDump && ((stepIdx-startOfDump) % dumpPeriod)==0)
        {
            sprintf(buffer,"./chain_%i/dump_chain_%i_%i.txt",chainIdx,chainIdx,stepIdx);
            fpw0=fopen(buffer,"w");
            for(int i=0; i<chainLength; i++){    
                fprintf(fpw0,"%i %.16Lf %.16Lf %.16Lf %.16Lf  %.16Lf %.16Lf\n",i+1, configR[i*3+0],configR[i*3+1],configR[i*3+2],Cumul_disp[i*3+0],Cumul_disp[i*3+1],Cumul_disp[i*3+2]); 
            }
            fclose(fpw0);
            
        }
            
    }  
}


//=============================
//==== random number functions
double myrandom() {
   return rand()/(double)RAND_MAX;
}

double normalrandom() {
   // return a normally distributed random value
   double v1=myrandom();
   double v2=myrandom();
   return cos(2.0*3.14*v2)*sqrt(-2.0*log(v1));
}
