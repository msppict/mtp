#include <stdio.h>
#include <pthread.h>
#include <signal.h>
#include <stdlib.h>
#include <stdint.h>
//#include "math.h"
#include "prog.h"


//void battery(double *Vload_prev, double *Vload, double *Iload_prev, double *Iload, double *temp, double *time, double *Vb_prev, double *Vc_prev, double *Vb, double *Vc)
void battery()
{
    float alpha = 0,
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
        SOC = 1320,
        Vload_prev = 0,
        Vload = 0,
        Iload_prev = 0,
        Iload = 0,
        temp = 0,
        time = 0,
        Vb_prev = 0,
        Vc_prev = 0,
        Vb = 0,
        Vc = 0;
        
	float k1 = 0,k2 = 0,k3 = 0,k4 = 0, k2_2 = 0, k3_2 = 0;     //Runge kutta method variables
	float l1 = 0,l2 = 0,l3 = 0,l4 = 0, l2_2 = 0, l3_2 = 0;
	float k1_temp0 = 0, k1_temp1 = 0, k1_temp2 = 0;
	float l1_temp0 = 0, l1_temp1 = 0, l1_temp2 = 0, l1_temp3 = 0;
	float k2_temp0 = 0, k2_temp1 = 0, k2_temp2 = 0,k2_temp3 = 0, k2_temp4 = 0, k2_temp5 = 0, k2_temp6 = 0;
	float l2_temp0 = 0, l2_temp1 = 0, l2_temp2 = 0,l2_temp3 = 0, l2_temp4 = 0, l2_temp5 = 0, l2_temp6 = 0;
	float k3_temp0 = 0, k3_temp1 = 0, k3_temp2 = 0,k3_temp3 = 0, k3_temp4 = 0, k3_temp5 = 0, k3_temp6 = 0;
	float l3_temp0 = 0, l3_temp1 = 0, l3_temp2 = 0,l3_temp3 = 0, l3_temp4 = 0, l3_temp5 = 0, l3_temp6 = 0;
	float k4_temp0 = 0, k4_temp1 = 0, k4_temp2 = 0,k4_temp3 = 0, k4_temp4 = 0, k4_temp5 = 0, k4_temp6 = 0;
	float l4_temp0 = 0, l4_temp1 = 0, l4_temp2 = 0,l4_temp3 = 0, l4_temp4 = 0, l4_temp5 = 0, l4_temp6 = 0;
	float alpha_temp1 = 0, alpha_temp = 0, beta_temp = 0, beta_temp1 = 0;
	float Vb_temp0 = 0, Vb_temp1 = 0, Vb_temp2 = 0, Vb_temp3 = 0, Vb_temp4 = 0;
	float Vc_temp0 = 0, Vc_temp1 = 0, Vc_temp2 = 0, Vc_temp3 = 0, Vc_temp4 = 0;
	

   float delta = 0, delta_2 = 0;

	// alpha = 1/(Cb*(Re+Rc));
	alpha_temp = fpadd32f(Re,Rc);
	alpha_temp1 = fpmul32f(Cb,alpha_temp);
	alpha = fdiv32(1,alpha_temp1);
	// beta = 1/(Cc*(Re+Rc));
	beta_temp = fpadd32f(Re,Rc);
	beta_temp1 = fpmul32f(Cc,beta_temp);
	beta = fdiv32(1,beta_temp1);
	
//    fprintf(stderr,"alpha = %20.18f, beta = %20.18f \n",alpha, beta);	

while(1)
{
     delta = time_period;

//  Vb' = -alpha(Vb - Vc + RcIload)
//  Vc' = -beta(-Vb +Vc + RcIload)

//     k1 = -(alpha*(Vb_prev - Vc_prev + (Rc*Iload_prev)));
	k1_temp0 = fpmul32f(Rc,Iload_prev);
	k1_temp1 = fpsub32f(Vb_prev, Vc_prev);
	k1_temp2 = fpadd32f(k1_temp0,k1_temp1);
	k1 = -(fpmul32f(alpha,k1_temp2));

	
//     l1 = -(beta*(-(*Vb_prev) + *Vc_prev + (Rc*(*Iload_prev))));


	l1_temp0 = fpmul32f(Rc,Iload_prev);
	l1_temp1 = fpsub32f(Vc_prev, Vb_prev);
	l1_temp2 = fpadd32f(l1_temp0,l1_temp1);
	l1_temp3 = -(fpmul32f(beta,l1_temp2));

//	delta_2 = delta/2
	delta_2 = fdiv32(delta,2);
	
//     k2 = -(alpha*(*Vb_prev + k1*(delta/2) -(*Vc_prev + l1*(delta/2)) + Rc*(*Iload_prev)));
	k2_temp0 = fpmul32f(Rc,Iload_prev);
	k2_temp1 = fpmul32f(l1,delta_2);
	k2_temp2 = fpadd32f(k2_temp0,k2_temp1);
	k2_temp3 = fpadd32f(Vc_prev,k2_temp2);
	k2_temp4 = fpmul32f(k1,delta_2);
	k2_temp5 = fpsub32f(k2_temp4,k2_temp3);
	k2_temp6 = fpadd32f(Vb_prev,k2_temp5);
	k2 = -(fpmul32f(alpha,k2_temp6));

	
//     l2 = -(beta*(-(*Vb_prev) - k1*(delta/2) +(*Vc_prev + l1*(delta/2)) + Rc*(*Iload_prev)));

	l2_temp0 = fpmul32f(Rc,Iload_prev);  //Rc*(*Iload_prev)
	l2_temp1 = fpmul32f(l1,delta_2);     //l1*(delta/2)
	l2_temp2 = fpadd32f(l2_temp0,l2_temp1);  //l1*(delta/2)) + Rc*(*Iload_prev)
	l2_temp3 = fpadd32f(Vc_prev,l2_temp2);  //(*Vc_prev + l1*(delta/2)) + Rc*(*Iload_prev)));
	l2_temp4 = fpmul32f(k1,delta_2);  // k1*(delta/2)
	l2_temp5 = fpsub32f(l2_temp3,l2_temp4);    // (*Vc_prev + l1*(delta/2)) + Rc*(*Iload_prev)) - k1*(delta/2)
	l2_temp6 = fpsub32f(l2_temp5,Vb_prev);    // (*Vc_prev + l1*(delta/2)) + Rc*(*Iload_prev)) - k1*(delta/2) - Vb_prev
	l2 = -(fpmul32f(beta,l2_temp6)); 

//     k2 = -(alpha*(*Vb_prev + k1*(delta/2) -(*Vc_prev + l1*(delta/2)) + Rc*(*Iload_prev)));
//     k3 = -(alpha*(*Vb_prev + k2*(delta/2) -(*Vc_prev + l2*(delta/2)) + Rc*(*Iload_prev)));
	k3_temp0 = fpmul32f(Rc,Iload_prev);
	k3_temp1 = fpmul32f(l2,delta_2);
	k3_temp2 = fpadd32f(k3_temp0,k3_temp1);
	k3_temp3 = fpadd32f(Vc_prev,k3_temp2);
	k3_temp4 = fpmul32f(k2,delta_2);
	k3_temp5 = fpsub32f(k3_temp4,k3_temp3);
	k3_temp6 = fpadd32f(Vb_prev,k3_temp5);
	k3 = -(fpmul32f(alpha,k3_temp6));

//      l2 = -(beta*(-(*Vb_prev) - k1*(delta/2) +(*Vc_prev + l1*(delta/2)) + Rc*(*Iload_prev)));
//      l3 = -(beta*(-(*Vb_prev) - k2*(delta/2) +(*Vc_prev + l2*(delta/2)) + Rc*(*Iload_prev)));


	l3_temp0 = fpmul32f(Rc,Iload_prev);
	l3_temp1 = fpmul32f(l2,delta_2);
	l3_temp2 = fpadd32f(l3_temp0,l3_temp1);
	l3_temp3 = fpadd32f(Vc_prev,l3_temp2);
	l3_temp4 = fpmul32f(k2,delta_2);
	l3_temp5 = fpsub32f(l3_temp3,l3_temp4);
	l3_temp6 = fpsub32f(l3_temp5,Vb_prev);
	l3 = -(fpmul32f(beta,l3_temp6));



//   k2 = -(alpha*(*Vb_prev + k1*(delta/2) -(*Vc_prev + l1*(delta/2)) + Rc*(*Iload_prev)));
//   k3 = -(alpha*(*Vb_prev + k2*(delta/2) -(*Vc_prev + l2*(delta/2)) + Rc*(*Iload_prev)));
//   k4 = -(alpha*(*Vb_prev + k1*(delta) -(*Vc_prev + l1*(delta)) + Rc*(*Iload_prev)));
	k4_temp0 = fpmul32f(Rc,Iload_prev);
	k4_temp1 = fpmul32f(l1,delta);
	k4_temp2 = fpadd32f(k4_temp0,k4_temp1);
	k4_temp3 = fpadd32f(Vc_prev,k4_temp2);
	k4_temp4 = fpmul32f(k1,delta);
	k4_temp5 = fpsub32f(k4_temp4,k4_temp3);
	k4_temp6 = fpadd32f(Vb_prev,k4_temp5);
	k4 = -(fpmul32f(alpha,k4_temp6));




//     l4 = -(beta*(-(*Vb_prev) - k1*(delta) +(*Vc_prev + l1*(delta)) + Rc*(*Iload_prev)));


	l4_temp0 = fpmul32f(Rc,Iload_prev);  //Rc*(*Iload_prev)
	l4_temp1 = fpmul32f(l1,delta);     //l1*(delta)
	l4_temp2 = fpadd32f(l4_temp0,l4_temp1);  //l1*(delta)) + Rc*(*Iload_prev)
	l4_temp3 = fpadd32f(Vc_prev,l4_temp2);  //(*Vc_prev + l1*(delta)) + Rc*(*Iload_prev)));
	l4_temp4 = fpmul32f(k1,delta);  // k1*(delta)
	l4_temp5 = fpsub32f(l4_temp3,l4_temp4);    // (*Vc_prev + l1*(delta)) + Rc*(*Iload_prev)) - k1*(delta)
	l4_temp6 = fpsub32f(l4_temp5,Vb_prev);    // (*Vc_prev + l1*(delta)) + Rc*(*Iload_prev)) - k1*(delta) - Vb_prev
	l4 = -(fpmul32f(beta,l4_temp6)); 
     



//      *Vb = (*Vb_prev) + delta*(k1 + 2*k2 + 2*k3 + k4)/6;
//	*Vc = (*Vc_prev) + delta*(l1 + 2*l2 + 2*l3 + l4)/6;
	Vb_temp0 = fpadd32f(k1,k4);
	k2_2 = fpmul32f(k2,2);
	k3_2 = fpmul32f(k3,2);
	Vb_temp1 = fpadd32f(k2_2,k3_2);
	Vb_temp2 = fpadd32f(Vb_temp0,Vb_temp1);
	Vb_temp3 = fdiv32f(Vb_temp2,3);
	Vb_temp4 = fpmul32f(Vb_temp3,delta_2);
	Vb = fpadd32f(Vb_prev, Vb_temp4);
	
//	*Vc = (*Vc_prev) + delta*(l1 + 2*l2 + 2*l3 + l4)/6;


	Vc_temp0 = fpadd32f(l1,l4);
	l2_2 = fpmul32f(l2,2);
	l3_2 = fpmul32f(l3,2);
	Vc_temp1 = fpadd32f(l2_2,l3_2);
	Vc_temp2 = fpadd32f(Vc_temp0,Vc_temp1);
	Vc_temp3 = fdiv32f(Vc_temp2,3);
	Vc_temp4 = fpmul32f(Vc_temp3,delta_2);
	Vc = fpadd32f(Vc_prev, Vc_temp4);
	 



//     *Vload_prev = *Vload;
//     *Iload_prev = *Iload;
     Vb_prev = Vb;
     Vc_prev = Vc;

     
     if(Iload_prev < 0){
	if(Iload_prev > Ibat_limit_rev)
	{
		Iload = Iload_prev;
	}
	else
	Iload = Ibat_limit_rev;}

	else{
	
     if(Iload_prev < Ibat_maxlimit){
     Iload = Iload_prev;   
     }
     else {
     Iload = Ibat_limit;     //current limiting to save the battery damage
     }       
}   
     x = Iload;
     Vload_prev = Vload;
     Vload = Cb*alpha*(Rc*(Vb) + Re*(Vc)) - ((Rt + (Rc*Re/(Rc+Re)))*(Iload));
    
     time = time + delta;
     SOC = Cb*(Vb) + Cc*(Vc);

	//Out puts from the module -- Pipes to other modules
	// Vload
	// I load
	// SOC
     }
     

}


