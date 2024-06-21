/*
 OBJECTIVE:
 07/26/2022: Cost function  check the blance equation 
 based on the formula obtained from Campbell_Formula
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "./RDP_history.h"

double fun_R(double z,double zf, double R0,double a1,double a2,double a3,double a4){
    return R0*(a1-(a1-1.0+a2*(z/zf)+a3*(z/zf)*(z/zf))*exp(-a4*(z/zf))); 
}

double fun_U(double z,double zf, double u0,double b1,double b2,double b3,double b4){
    return u0*(b1-(b1-1.0+b2*(z/zf)+b3*(z/zf)*(z/zf))*exp(-b4*(z/zf))); 
}

double fun_dRdz(double z,double zf, double R0,double a1,double a2,double a3,double a4){
    return (exp(-a4*(z/zf))*(a4*(a1-1.0+a3*(z/zf)*(z/zf))+a2*(a4*(z/zf)-1.0)-2.0*a3*(z/zf)))*R0/zf; 
}

double fun_dudz(double z,double zf, double u0,double b1,double b2,double b3,double b4){
    return (exp(-b4*(z/zf))*(b4*(b1-1.0+b3*(z/zf)*(z/zf))+b2*(b4*(z/zf)-1.0)-2.0*b3*(z/zf)))*u0/zf; 
}

double fun_d2Rdz2(double z,double zf, double R0,double a1,double a2,double a3,double a4){
    return (exp(-a4*(z/zf))*(a4*(-a1*a4+a2*(2.0-a4*(z/zf))+a4)-a3*(a4*a4*(z/zf)*(z/zf)-4.0*a4*(z/zf)+2.0)))*R0/(zf*zf);    
}


double  CostFunction(double *X){
    double cost; 
//===parameters from the process line 
    double u0    =  0.0160;           //velocity at the die
    double u1    =  0.0833;           //velocity at the forst line 
    double R0    =  0.0192;           //radius of the die 
    double Rf    =  0.0384;           //radius of the bubble at the forst line  
    double zf    =  0.2;             //the length of the forst line  
    double Q     =  1.2654e-6;       //assuming the Q value does not change during the blowing 
    
 
    double a3 =   0.0; 
    double b3 =   0.0;
    
    double a4 = X[0];
    double b4 = X[1];
    
    if (abs(a4)<1e-20 || abs(b4)<1e-20){ //exclude the condition of zero 
        cost =1e20; 
        return cost; 
    }
    else{
        
    double a1, a2, b1, b2; 
    
    double R_ratio=Rf/R0;
    double U_ratio=u1/u0;
     
    a1=(a3+1.0+R_ratio*exp(a4)*(a4-1.0))/(a4+(exp(a4)-1.0)*(a4-1.0));
    a2=(a1-R_ratio)*exp(a4)-(a1-1.0)-a3;
     
    b1=(b3+1.0+U_ratio*exp(b4)*(b4-1.0))/(b4+(exp(b4)-1.0)*(b4-1.0));
    b2=(b1-U_ratio)*exp(b4)-(b1-1.0)-b3;
    
    int Nz=1000; 
    double Z[Nz],zbin;
    for(int i=0; i<Nz; i++) Z[i]= (double) 1.0*i*zf/Nz; 
    zbin=Z[1]-Z[0];
    
    //get the strain rate histories and the parameters for force balance check 
    double R_value[Nz], dR_dzvalue[Nz],cos_value[Nz];
    double U_value[Nz], du_dZ[Nz],Vz[Nz];
    double edot[Nz][3];
    double inv_Rm[Nz], inv_Rp[Nz], thickness[Nz]; 
    //FILE *fpw0, *fpw1;
    //fpw0=fopen("R_Check.dat","w");
    //fpw1=fopen("U_Check.dat","w");
    
    for(int i=0; i<Nz; i++){
        R_value[i]       = fun_R(Z[i],zf,R0,a1,a2,a3,a4);
        dR_dzvalue[i]    = fun_dRdz(Z[i],zf,R0,a1,a2,a3,a4);
        cos_value[i]     = cos(atan(dR_dzvalue[i]));
        U_value[i]       = fun_U(Z[i],zf,u0,b1,b2,b3,b4);
        du_dZ[i]         = fun_dudz(Z[i],zf,u0,b1,b2,b3,b4);
        edot[i][0]       = cos_value[i]*du_dZ[i]; // mm 
        edot[i][1]       = U_value[i]*cos_value[i]*dR_dzvalue[i]/R_value[i]; //pp 
        edot[i][2]       = -edot[i][0]-edot[i][1]; //nn 
        Vz[i]            = U_value[i]*cos_value[i];
        inv_Rm[i]        = -fun_d2Rdz2(Z[i],zf,R0,a1,a2,a3,a4)*pow(cos_value[i],3.0);
        inv_Rp[i]        = cos_value[i]/R_value[i]; 
        thickness[i]     = Q/R_value[i]/U_value[i]/(2.0*M_PI); 
        
        //fprintf(fpw0,"%10.5e %10.5e %10.5e  %10.5e %10.5e\n",Z[i],R_value[i],dR_dzvalue[i],inv_Rm[i],inv_Rp[i]);
        //fprintf(fpw1,"%10.5e %10.5e %10.5e %10.5e %10.5e %10.5e\n",Z[i],Vz[i],du_dZ[i],edot[i][0],edot[i][1],edot[i][2]);
    }
    //fclose(fpw0);
    //fclose(fpw1);
   // Run the RDP model
   int out_count; 
   out_count=RDP(Nz, Z, edot,Vz);
   
   double  stress[out_count][3];
   double  Balance[out_count][4];
   FILE *par  =fopen("stress.dat","r");
   for(int i=0; i<out_count; i++) fscanf(par, "%le %le %le", &stress[i][0], &stress[i][1], &stress[i][2]);
   // Evaluate the Cost function 
   fclose(par); 
   
   //printf("#Nstress %d\n", out_count);
   //par  =fopen("stress_check.dat","w");
   //for(int i=0; i<out_count; i++) fprintf( par, "%le %le %le %le\n",stress[i][0], stress[i][1], stress[i][2],stress[i][3]);
   //fclose(par); 
   
   double zz, c_R, c_cos,c_thick,c_inv_Rm,c_inv_Rp; 
   for(int i=0; i<out_count; i++){
        zz=stress[i][0];
        // do the linear interpolate
        int index = (int) floor(zz/zbin);
        if(index<Nz-1 && index>=0){
            c_R       = R_value[index]+(zz-Z[index])*(R_value[index+1]-R_value[index])/(Z[index+1]-Z[index]);
            c_cos     = cos_value[index]+(zz-Z[index])*(cos_value[index+1]-cos_value[index])/(Z[index+1]-Z[index]);
            c_thick   = thickness[index]+(zz-Z[index])*(thickness[index+1]-thickness[index])/(Z[index+1]-Z[index]);
            c_inv_Rm  = inv_Rm[index]+(zz-Z[index])*(inv_Rm[index+1]-inv_Rm[index])/(Z[index+1]-Z[index]);
            c_inv_Rp  = inv_Rp[index]+(zz-Z[index])*(inv_Rp[index+1]-inv_Rp[index])/(Z[index+1]-Z[index]);
            
        }else{
        printf("!!! indiex is not correct in cost function %f  %d\n",zz, index);
        }
       Balance[i][0] = zz;
       Balance[i][1] = (stress[i][1]*c_inv_Rm + stress[i][2]*c_inv_Rp)*c_thick;
       Balance[i][2] = 2.0*M_PI*c_R*c_thick*stress[i][1]*c_cos;
       Balance[i][3] = c_R*c_R-R0*R0; 
   }
   
   double Mean_DP, Mean_F0; 
   Mean_DP=0.0;
   Mean_F0=0.0; 
   
    for(int i=0; i<out_count; i++) {
        Mean_DP += Balance[i][1]/out_count; 
        Mean_F0 += (Balance[i][2]-Balance[i][1]*Balance[i][3]*M_PI)/out_count; 
    }
    
    double var1, var2;
    cost=0.0;
    for(int i=0; i<out_count; i++){
        var1 = pow(Balance[i][1]-Mean_DP,2.0); 
        var2 = pow(Balance[i][2]-Balance[i][1]*Balance[i][3]*M_PI-Mean_F0,2.0); 
        cost += pow(M_PI*R0*R0,2.0)*var1+var2; 
    }
    cost =cost/out_count; 
    return cost;
    }

}
