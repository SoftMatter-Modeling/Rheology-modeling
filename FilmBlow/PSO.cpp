/*
 OBJECTIVE:
 07/26/2022: PSO C code based on the PSO matlab code 
 single core version 
 
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "./Rate_RDP_Cost.h"

//====random number generator function ======
//=========================================
double  alea(double a, double b){
    //Random real  number between a and b
      double r; 
      double uni_r; 
      
      uni_r= ((double) rand())/(double) RAND_MAX;
      r= a+uni_r*(b-a);
      return r; 
}

double alea_normal (double mean, double std_dev){
    // Use the polar form of the Box-Muller transformation to obtain a pseudo
    // random number from a Gaussian distribution
    double w=2.0;
    double x1, x2;
    double y1;
    
    while (w>=1.0){
        x1=2.0*alea(0,1)-1.0;
        x2=2.0*alea(0,1)-1.0;
        w=x1*x1+x2*x2;
    }
     w = sqrt(-2.0*log(w)/w);
     y1 = x1*w;
    
    if (alea(0.0,1.0)<0.5) y1=-y1;
    y1 = y1*std_dev+mean;
    return y1; 
}

void alea_sphere(int D, double radius, double *x){
//  ******* Random point in a hypersphere ********
// Put a random point inside the hypersphere S(0,radius) (center 0, radius 1).

    double l=0.0; 
    for(int j=0; j<D; j++){
        x[j]=alea_normal(0.0,1.0);
        l +=pow(x[j],2.0);
    }
    l=sqrt(l);
    
    double r; 
    r = alea(0.0,1.0);
    r = pow(r,1.0/D);
    for(int j=0; j<D; j++) x[j]=r*radius*x[j]/l;
}
//=========================================
int main(int argc, char** argv){
    
    srand((unsigned)time(NULL)); // for different random number 
    int    nVar=2; //Number of Decision ai  Variables 
    double VarMin[nVar], VarMax[nVar]; 
    
//====define the bounadry 
    for(int i=0; i<nVar; i++){
        VarMin[i]= -20.0;
        VarMax[i]=  40.0;
    }
    
// =====PSO Parameters
    int MaxIt =1000;      // Maximum Number of Iterations
    int nPop  =10;         // Population Size (Swarm Size)
    double w=1.0/(2.0*log(2.0));       // Inertia Weight
    double wdamp=1.0;                  // Inertia Weight Damping Ratio
    double c1=0.5*log(2.0);            // Personal Learning Coefficient
    double c2=0.5+log(2.0);            // Global Learning Coefficient
   
    double particle_Position[nPop][nVar], particle_Velocity[nPop][nVar], particle_Cost[nPop];
    double particle_Best_Pos[nPop][nVar], particle_Best_Cost[nPop]; 

    double GlobalBest_Cost = 9.0e60, GlobalBest_Position[nVar];
    int    IndexGlobal  =0; 
 
// =====initilization 
    for (int i=0; i<nPop;i++){
       for(int j=0; j<nVar; j++){
           //Initialize Position
           particle_Position[i][j]  =alea(VarMin[j],VarMax[j]);
           //Initialize Velocity
           particle_Velocity[i][j]  =alea(VarMin[j]-particle_Position[i][j],VarMax[j]-particle_Position[i][j]);
       }
       //Evaluation
       particle_Cost[i]      = CostFunction(particle_Position[i]); // program the CostFunction later 
       //printf("#Init, pop:%d, position=%le %le, cost=%le\n",i,particle_Position[i][0],particle_Position[i][1], particle_Cost[i]);
       //Update Personal Best
       for(int j=0; j<nVar; j++) particle_Best_Pos[i][j]  = particle_Position[i][j];
       particle_Best_Cost[i] = particle_Cost[i];
       
      //Update Global Best
      if (GlobalBest_Cost > particle_Best_Cost[i]){
          GlobalBest_Cost = particle_Best_Cost[i];
          IndexGlobal=i;
          for(int j=0; j<nVar; j++) GlobalBest_Position[j] = particle_Best_Pos[i][j];
      }
   }
   
   double BestCost[MaxIt], BestPosition[MaxIt][nVar];
    
   for(int it=0; it<MaxIt; it++){
        BestCost[it]  = 0.0; 
        for(int j=0; j<nVar; j++) BestPosition[it][j] = 0.0; 
    }
    
    double p_x_p[nVar], p_x_l[nVar];
    double Gi[nVar],x_p[nVar]; 
    double rad; 
    
    FILE *fpw0;
    char buffer[256]; 
    
    for (int it=0; it<MaxIt;it++){
        
        for(int i=0; i<nPop; i++){
            
                rad=0.0; 
                for(int j=0; j<nVar; j++){
                    //Define the centre gravity 
                    //define a point p' on x-p, beyond p
                    p_x_p[j] = particle_Position[i][j] + c1*(particle_Best_Pos[i][j]-particle_Position[i][j]);
                    // define a point g' on x-g, beyond g
                    p_x_l[j] = particle_Position[i][j] + c2*(GlobalBest_Position[j]-particle_Position[i][j]);
                
                    if  (IndexGlobal==i)  Gi[j]= 0.5*(particle_Position[i][j]+ p_x_p[j]);
                    else Gi[j] = (particle_Position[i][j]+p_x_p[j]+p_x_l[j])/3.0;
                    
                    rad += pow((Gi[j]-particle_Position[i][j]),2.0); 
               }
                rad = sqrt(rad); //radius = Euclidean norm of x-Gi
                
                double  sph_ran[nVar];
                alea_sphere(nVar,rad,sph_ran);
                
                for(int j=0; j<nVar; j++){
                    x_p[j]= sph_ran[j]+ Gi[j];     //Generate a random point in the hyper-sphere around G 
                    //Update Velocity
                    particle_Velocity[i][j]=w*particle_Velocity[i][j]+x_p[j]-particle_Position[i][j]; 
                    //Update Position
                    particle_Position[i][j]=particle_Position[i][j]+particle_Velocity[i][j]; 
                    
                    // Boundary confinment 
                    if(particle_Position[i][j] < VarMin[j]){
                        
                        particle_Position[i][j] = VarMin[j];
                        particle_Velocity[i][j] = -0.5*particle_Velocity[i][j];
                    }
                    
                   if(particle_Position[i][j] > VarMax[j]){
                        particle_Position[i][j] =VarMax[j];
                        particle_Velocity[i][j] =-0.5*particle_Velocity[i][j];
                    }
                }
                
                //Evaluation
                particle_Cost[i]  = CostFunction(particle_Position[i]); // program the CostFunction later 
                //printf("#step %d, pop:%d, position=%le %le, cost=%le\n",it,i,particle_Position[i][0],particle_Position[i][1], particle_Cost[i]);
                //Update Personal Best
                if(particle_Cost[i]<particle_Best_Cost[i]){
                      for(int j=0; j<nVar; j++) particle_Best_Pos[i][j]=particle_Position[i][j];
                      particle_Best_Cost[i] = particle_Cost[i];
                }

       
                //Update Global Best
                if (GlobalBest_Cost > particle_Best_Cost[i]){
                GlobalBest_Cost     = particle_Best_Cost[i];
                IndexGlobal=i;
                for(int j=0; j<nVar; j++) GlobalBest_Position[j] = particle_Best_Pos[i][j];
                }
        }

        w=w*wdamp; 
        BestCost[it]=GlobalBest_Cost;
        for(int j=0; j<nVar; j++) BestPosition[it][j]=GlobalBest_Position[j];
        printf("#step=%d, cost=%f, position=%f %f\n",it,GlobalBest_Cost,GlobalBest_Position[0],GlobalBest_Position[1]);
        sprintf(buffer,"Step_%i.dat",it);
        fpw0=fopen(buffer,"w");
        fprintf(fpw0,"%d  %le %le %le\n",it,GlobalBest_Cost,GlobalBest_Position[0],GlobalBest_Position[1]);
        fclose(fpw0);
    }

    


}
