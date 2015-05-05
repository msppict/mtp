#include <stdio.h>
#include <pthread.h>
#include <signal.h>
#include <stdlib.h>
#include <stdint.h>
#include "math.h"



void battery(double *Vload_prev, double *Vload, double *Iload_prev, double *Iload, double *temp, double *time, double *Vb_prev, double *Vc_prev, double *Vb, double *Vc)
{
    double alpha = 0,
        beta = 0,
        Rc = 10,
        Rt = 20,
        Re = 10,     //circuit parameters , values to be found out from proper source
        Cb = 600,
        Cc = 50,
        time_period = 50e-6,
        Ibat_limit = 90,
	Ibat_limit_rev = -90,
        x = 0,
        Ibat_maxlimit = 90,
        SOC = 1320;
        
    double k1 = 0,k2 = 0,k3 = 0,k4 = 0;     //Runge kutta method variables
    double l1 = 0,l2 = 0,l3 = 0,l4 = 0;
    
   double delta = 0;
   char buffer[BUFSIZ];


	char filename_Vb[] 	= "/home/mayur/RnDproject/ahirmaster/Documents/Tutorials/MTP/battery/src/Vb.m";
   	char filename_soc[] 	= "/home/mayur/RnDproject/ahirmaster/Documents/Tutorials/MTP/battery/src/soc.m";
	FILE *fp_Vb = NULL;

    alpha = 1/(Cb*(Re+Rc));
    beta = 1/(Cc*(Re+Rc));
    fprintf(stderr,"alpha = %20.18f, beta = %20.18f \n",alpha, beta);	
     if ( 0.0 == *time ) {  

     *Vb_prev = 0.0;
     *Vc_prev = 0.0;
     *Vb = 0.0;
     *Vc = 0.0;
     *Iload_prev = 0.0;
     *Iload = 0.0;  
     *Vload_prev = 0.0;
     *Vload = 0.0;
	*time = *time + delta;

   	
     	fp_Vb = fopen(filename_Vb, "w"); 
	if(fp_Vb == NULL) {
		printf("Failed to open file for writing\n");
	}
	fputs("Vb, Here it is! Herr Hitler! = [ ", fp_Vb); 
	fclose(fp_Vb);

/*
   	char filename_Vb[] 	= "/home/mayur/RnDproject/ahirmaster/Documents/Tutorials/MTP/battery/src/Vb.m";
   	char filename_soc[] 	= "/home/mayur/RnDproject/ahirmaster/Documents/Tutorials/MTP/battery/src/soc.m";
	FILE *fp_Vb = NULL;
     	fp_Vb = fopen(filename_Vb, "w"); 
	if(fp_Vb == NULL) {
		printf("Failed to open file for writing\n");
	}
	fputs("Vb, Here it is! Herr Hitler! = [ ", fp_Vb); 
	fclose(fp_Vb);
*/
	fprintf(stderr,"Initialization done!");
       
		
	
     }

     else {

     delta = time_period;
     *Vload_prev = *Vload;
     *Iload_prev = *Iload;
     *Vb_prev = *Vb;
     *Vc_prev = *Vc;

//  Vb' = -alpha(Vb - Vc + RcIload)
//  Vc' = -beta(-Vb +Vc + RcIload)

     k1 = -(alpha*(*Vb_prev - *Vc_prev + (Rc*(*Iload_prev))));

     l1 = -(beta*(-(*Vb_prev) + *Vc_prev + (Rc*(*Iload_prev))));

     k2 = -(alpha*(*Vb_prev + k1*(delta/2) -(*Vc_prev + l1*(delta/2)) + Rc*(*Iload_prev)));

     l2 = -(beta*(-(*Vb_prev) - k1*(delta/2) +(*Vc_prev + l1*(delta/2)) + Rc*(*Iload_prev)));
     
     k3 = -(alpha*(*Vb_prev + k2*(delta/2) -(*Vc_prev + l2*(delta/2)) + Rc*(*Iload_prev)));

     l3 = -(beta*(-(*Vb_prev) - k2*(delta/2) +(*Vc_prev + l2*(delta/2)) + Rc*(*Iload_prev)));

     k4 = -(alpha*(*Vb_prev + k1*(delta) -(*Vc_prev + l1*(delta)) + Rc*(*Iload_prev)));

     l4 = -(beta*(-(*Vb_prev) - k1*(delta) +(*Vc_prev + l1*(delta)) + Rc*(*Iload_prev)));

     

	 *Vb = (*Vb_prev) + delta*(k1 + 2*k2 + 2*k3 + k4)/6;
	 *Vc = (*Vc_prev) + delta*(l1 + 2*l2 + 2*l3 + l4)/6;
     
     if(*Iload_prev < 0){
	if(*Iload_prev > Ibat_limit_rev)
	{
		*Iload = *Iload_prev;
	}
	else
	*Iload = Ibat_limit_rev;
}
else{
	
     if(*Iload_prev < Ibat_maxlimit){

     *Iload = *Iload_prev;   
     }
     else {
     *Iload = Ibat_limit;     //current limiting to save the battery damage
     }       
}   
     x = *Iload;
     *Vload_prev = *Vload;
     *Vload = Cb*alpha*(Rc*(*Vb) + Re*(*Vc)) - ((Rt + (Rc*Re/(Rc+Re)))*(*Iload));
    
     *time = *time + delta;
     SOC = Cb*(*Vb) + Cc*(*Vc);
     printf("\nIn battery module,\n SOC = %20.18f, \t Vload =%20.18f, \n Vb = %20.18f, Vc = %20.18f ",SOC,*Vload,*Vb,*Vc); 

     fp_Vb = fopen(filename_Vb, "a"); // open file for appending !!! 
		fprintf(fp_Vb," %20.18f ",*Vb);
		fclose(fp_Vb);    

     }
     

}

void main()
{
printf("\nHello world!\n");
int i,j;
double Vload_prev =54 , Vload = 55, Iload_prev = -3, Iload = -4444631.1, temp = 25, time = 1, Vb_prev = 5, Vc_prev = 5, Vb=5, Vc=4;
/*
   	char filename_Vb[] 	= "/home/mayur/RnDproject/ahirmaster/Documents/Tutorials/MTP/battery/src/Vb.m";
   	char filename_soc[] 	= "/home/mayur/RnDproject/ahirmaster/Documents/Tutorials/MTP/battery/src/soc.m";
	FILE *fp_Vb = NULL;
     	fp_Vb = fopen(filename_Vb, "w"); 
	if(fp_Vb == NULL) {
		printf("Failed to open file for writing\n");
	}
	fputs("Vb, Here it is! Herr Hitler! = [ ", fp_Vb); 
	fclose(fp_Vb);
*/
for(j=0;j<100;j++)
{
for(i=0;i<200;i++)
{
//printf("\n In Main..hi!\n");
battery(&Vload_prev, &Vload, &Iload_prev, &Iload, &temp, &time, &Vb_prev, &Vc_prev, &Vb, &Vc);

/*     fp_Vb = fopen(filename_Vb, "a"); // open file for appending !!! 
		fprintf(fp_Vb," %20.18lf ",&Vb);
		fclose(fp_Vb);
*/
}
printf("\n Value after %dth iteration \n Iload = %20.18f",j*i,Iload);


}
printf("\nOver!\n");
    
/*
fp_Vb = fopen(filename_Vb, "a"); // open file for appending !!! 
	fprintf(fp_Vb," ]; ");
	fclose(fp_Vb);
*/
}
