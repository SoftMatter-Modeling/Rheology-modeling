#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>                  // Library for file removal and rename
#include <chrono>
#include <iomanip>
#include <sstream>                  // For the int to string conversion
#include <cstdlib>

using namespace std; 

double repcr(double taud1, double taud2,  double atube)
{
	double da;
	da=-1.0/taud1*(atube-1)-1.0/taud2*(atube-1);
	return da;
}

double ret(double l, double taus, double fe, double atube, int q)
{
	double da;
	if (q>1)
	{
		da=-2.0*(1-1.0/l)/taus*atube*fe*exp(2.0/double(q)*(l-1));
	}
	else
		{da=-2.0*(1-1.0/l)/taus*fe*atube;}
	return da;
}


double ccr(double l1,double l2,  double taus, double fe, double atube, int q)
{
	double da;
	double beta_ccr=1.0, delta_ccr=-0.5;
	if (q>1)
		{da=-2.0*beta_ccr*(1.0-1.0/l1)/taus*fe*pow(l2,2.0*delta_ccr)*(atube-1)*exp(2.0/q*(l1-1));}
	else
		{da=-2.0*beta_ccr*(1.0-1.0/l1)/taus*fe*pow(l2, 2.0*delta_ccr)*(atube-1);}
	return da;
}


int  RDP(int N_filed, double *Z, double edot[][3], double *Vz)
{

	int num; // number of ensembles 

	FILE * par  =fopen("component.dat","r");
    FILE * hist_tem =fopen("L10Fitted_Tempzf0.2.dat","r");
	FILE * output=fopen("stress.dat","w");

	fscanf(par,"%d", &num);
	double phi[num], taud[num], taus[num];
	int    q[num];
	double l[num], fe[num];
	double lmax[num];   //may be changed per case
	double a1tube_ij[num][num], a2tube_ij[num][num], a3tube_ij[num][num];
	double G0;
    int    dumpPeriod=5000;

//====initialize the data: a1tube_ij[i][j], a2tube_ij[i][j], and  a3tube_ij[i][j] 
//====  l[i] fe[i] are uptating every time step 

    for (int i=0; i<num;i++)
	{ 
		fscanf(par, "%le %le %le %d %le", &phi[i], &taud[i], &taus[i], &q[i], &lmax[i]);
		l[i]=1.0; 
		fe[i]=(1.0-1.0/lmax[i]/lmax[i])/(1.0-l[i]*l[i]/lmax[i]/lmax[i]);

		for (int j=0;j<num;j++)
		{
			a1tube_ij[i][j]=1.0; 
			a2tube_ij[i][j]=1.0;
            a3tube_ij[i][j]=1.0;
		}
	}
	fscanf(par, "%le", &G0);

//====the history of the filed. 
    double Ea =  66.0; 
    double T0 =  150.0;
    double Rgas= 8.31446261815324/1000.0;

    double temp[N_filed];     

    for (int t=0; t<N_filed;t++) 
        {
        //fscanf(hist_str,"%le %le %le %le %le", &Z[t], &edot[t][0],&edot[t][1],&edot[t][2], &Vz[t]);
        fscanf(hist_tem,"%le", &temp[t]);
        }
    fclose (par);  
    fclose (hist_tem);
    
    //FILE *fpw0=fopen("hist.dat","w");
    //for (int t=0; t<N_filed;t++){
    //   fprintf(fpw0,"%10.5e %10.5e %10.5e %10.5e %10.5e %10.5e\n",Z[t],Vz[t],temp[t],edot[t][0],edot[t][1],edot[t][2]); 
        
    //}
    //fclose (fpw0); 
   
    double zbin = Z[2]-Z[1];
    double zmax = Z[N_filed-1]-2.0*zbin;
    
    double zz    = 1.0e-7;
    int    zcout = 0;
    int    out_count=0; 
    double zstep = 1.0e-7;
    
//====start the integration  
	while (zz<zmax)
	{
//====do the linear interpolate of history profile 
        double c_temp;
        double c_Vz; 
        double c_edot[3];
        int index = (int) floor(zz/zbin);  
        if(index< N_filed-1 && index>=0){
            c_temp    = temp[index]+(zz-Z[index])*(temp[index+1]-temp[index])/(Z[index+1]-Z[index]);
            c_Vz      = Vz[index]+(zz-Z[index])*(Vz[index+1]-Vz[index])/(Z[index+1]-Z[index]);
            c_edot[0] = edot[index][0]+(zz-Z[index])*(edot[index+1][0]-edot[index][0])/(Z[index+1]-Z[index]);
            c_edot[1] = edot[index][1]+(zz-Z[index])*(edot[index+1][1]-edot[index][1])/(Z[index+1]-Z[index]);
            c_edot[2] = edot[index][2]+(zz-Z[index])*(edot[index+1][2]-edot[index][2])/(Z[index+1]-Z[index]);
        }else{
        printf("!!! indiex is not correct in RDP %d\n",index);
        break;
        }
        
//=====calculate the relaxation time at the current step
        double taud_scaled[num], taus_scaled[num];
        double a1tube[num], a2tube[num],a3tube[num]; //, s1[num], s2[num], s3[num], eta;
        double factor; 
        
        for (int i=0;i<num;i++)
        {
            factor=exp(Ea*(1.0/(c_temp+273.15)-1.0/(T0+273.15))/Rgas);
            taud_scaled[i]=taud[i]*factor; 
            taus_scaled[i]=taus[i]*factor; 
            //printf("%le %le \n", taud_scaled[i],taus_scaled[i]);
          
            a1tube[i]=0;
            a2tube[i]=0;
            a3tube[i]=0;
        }
        
//=====calculate the conformational tensor 
        double da[3]; 
        double d_repcr[3];
        double d_ret[3];
        double d_ccr[3];
        
    for (int i=0;i<num;i++)
       {
       for (int j=0;j<num;j++)
        {
                d_repcr[0]=repcr(taud_scaled[i], taud_scaled[j], a1tube_ij[i][j]); 
                d_repcr[1]=repcr(taud_scaled[i], taud_scaled[j], a2tube_ij[i][j]);
                d_repcr[2]=repcr(taud_scaled[i], taud_scaled[j], a3tube_ij[i][j]);
                
                d_ret[0]=ret(l[i], taus_scaled[i], fe[i], a1tube_ij[i][j],q[i]);
                d_ret[1]=ret(l[i], taus_scaled[i], fe[i], a2tube_ij[i][j],q[i]);
                d_ret[2]=ret(l[i], taus_scaled[i], fe[i], a3tube_ij[i][j],q[i]);
                
                d_ccr[0]=ccr(l[j], l[i],taus_scaled[j], fe[j], a1tube_ij[i][j], q[j]);
                d_ccr[1]=ccr(l[j], l[i],taus_scaled[j], fe[j], a2tube_ij[i][j], q[j]);
                d_ccr[2]=ccr(l[j], l[i],taus_scaled[j], fe[j], a3tube_ij[i][j], q[j]);
                
                //printf("%le %le %le\n", d_repcr[0],d_ret[0], d_ccr[0]);
                da[0]  = 2.0*c_edot[0]*a1tube_ij[i][j]+d_repcr[0]+d_ret[0]+d_ccr[0];
                da[1]  = 2.0*c_edot[1]*a2tube_ij[i][j]+d_repcr[1]+d_ret[1]+d_ccr[1];
                da[2]  = 2.0*c_edot[2]*a3tube_ij[i][j]+d_repcr[2]+d_ret[2]+d_ccr[2];
                
                a1tube_ij[i][j] += da[0]*zstep/c_Vz;
                a2tube_ij[i][j] += da[1]*zstep/c_Vz;
                a3tube_ij[i][j] += da[2]*zstep/c_Vz;
                
                a1tube[i]+= phi[j]*a1tube_ij[i][j];
                a2tube[i]+= phi[j]*a2tube_ij[i][j];
                a3tube[i]+= phi[j]*a3tube_ij[i][j];
            }
       }
	    
 //=====update the stress, stretch value 
    double trA;
    double sig1,sig2,sig3;
    sig1=0.0;
    sig2=0.0;
    sig3=0.0;
   
    for (int i=0;i<num;i++)
    {
        trA    = a1tube[i]+a2tube[i]+a3tube[i];
        l[i]   = sqrt((trA/3.0));     // update the strate ratio 
        fe[i]  = (1.0-1.0/lmax[i]/lmax[i])/(1.0-l[i]*l[i]/lmax[i]/lmax[i]); //update the fe 
         //
        //s1[i]  =  1.5*(a1tube[i]/trA-1.0/3.0); //add the information later 
        //ssum   += phi[i]*s1[i];
        sig1  += G0*phi[i]*fe[i]*a1tube[i];
        sig2  += G0*phi[i]*fe[i]*a2tube[i];
        sig3  += G0*phi[i]*fe[i]*a3tube[i];
    }

    if( (zcout % dumpPeriod)==0){
        out_count +=1; 
        fprintf(output,"%le %le %le\n",zz, sig1-sig3, sig2-sig3);
    }
        
     zz += zstep;
     zcout +=1;

	}
    fclose (output); 
    
	//for (int i=0;i<num;i++) printf("%le \n", l[i]);
    return out_count; 
    
	//auto stop = high_resolution_clock::now(); 
	//auto duration = duration_cast<microseconds>(stop - start); 
  
    //cout << "Time taken by function: "
    //<< duration.count() << " microseconds" << endl; 

}
