/*
float fpadd32f(float x, float y);
float fpsub32f(float x, float y);
uint32_t fpadd32fi(uint32_t x, uint32_t y);
uint32_t fpsub32fi(uint32_t x, uint32_t y);
float fpmul32f(float x, float y);
float fdiv32(float a, float b);
*/
#include "prog.h"

void motor()   
{
//float iq, iq_prev, id, id_prev, flq, flq_prev, fld, fld_prev, spd_prev, vd, vq, torque, load_torque, time, spd, theta, theta_prev, theta_sin,theta_cos;


	float	alpha = 0,
		beta = 0,
		sigma = 0,
		mu = 0,
		gamma = 0,
		omega = 314.15,
		inertia = 0.026,
		rs = 4.9,
		rr = 8.1,
		lls = 0.03215,
		llr = 0.03215,
		ls = 0,
		lr = 0,
		lm = 0.8096,
		poles = 4,
		cont = 0,
		time_period = 50e-6,
		time = 0;
                                                   
	float k1 = 0,k2 = 0,k3 = 0,k4 = 0;
	float l1 = 0,l2 = 0,l3 = 0,l4 = 0;
	float m1 = 0,m2 = 0,m3 = 0,m4 = 0;
	float n1 = 0,n2 = 0,n3 = 0,n4 = 0;
	float o1 = 0,o2 = 0,o3 = 0,o4 = 0;
	float delta = 0;
	float f_rotor = 0, omega_r = 0;
	float sigma1 = 0, sigma2 = 0;
	float beta1 = 0;
	float mu1 = 0, mu2 = 0, mu3 = 0, mu4 = 0;
	float gamma1 = 0, gamma2 = 0, gamma3 = 0, gamma4 = 0, gamma5 = 0, gamma6 = 0;
 	float cont1 = 0, cont2 = 0, cont3 = 0;
	float 	iq_prev = 0.0,
	id_prev = 0.0,
	flq_prev = 0.0,
	fld_prev = 0.0,
	spd_prev = 0.0,
	iq = 0.0,
	id = 0.0,
	flq = 0.0,
	fld = 0.0,
	spd = 0.0,
	torque = 0.0,
	vd = 0,
	vq = 0,
	load_torque = 0,
//	iq_prev =  iq,
//	id_prev =  id,
//	flq_prev =  flq,
//	fld_prev =  fld,
//	spd_prev =  spd,
//	delta = time_period,
//	time =  time + delta,
	theta = 0,
	theta_prev = 0,
	theta_sin = 0,
	theta_cos = 0,
	delta_2 = 0;
//	omega = 0;

	float k1_1 = 0, k1_2 = 0, k1_3 = 0, k1_4 = 0, k1_5 = 0, k1_6 = 0, k1_7 = 0, k1_8 = 0, k1_9 = 0, k1_10 = 0;
	float l1_1 = 0, l1_2 = 0, l1_3 = 0, l1_4 = 0, l1_5 = 0, l1_6 = 0, l1_7 = 0, l1_8 = 0, l1_9 = 0, l1_10 = 0;
	float m1_1 = 0, m1_2 = 0, m1_3 = 0, m1_4 = 0, m1_5 = 0, m1_6 = 0, m1_7 = 0, m1_8 = 0, m1_9 = 0, m1_10 = 0;
	float n1_1 = 0, n1_2 = 0, n1_3 = 0, n1_4 = 0, n1_5 = 0, n1_6 = 0, n1_7 = 0, n1_8 = 0, n1_9 = 0, n1_10 = 0;
	float o1_1 = 0, o1_2 = 0, o1_3 = 0, o1_4 = 0, o1_5 = 0, o1_6 = 0, o1_7 = 0, o1_8 = 0, o1_9 = 0, o1_10 = 0;
	float k2_1 = 0, k2_2 = 0, k2_3 = 0, k2_4 = 0, k2_5 = 0, k2_6 = 0, k2_7 = 0, k2_8 = 0, k2_9 = 0, k2_10 = 0;
	float k2_11 = 0, k2_12 = 0, k2_13 = 0, k2_14 = 0, k2_15 = 0, k2_16 = 0, k2_17 = 0, k2_18 = 0, k2_19 = 0;
	float l2_1 = 0, l2_2 = 0, l2_3 = 0, l2_4 = 0, l2_5 = 0, l2_6 = 0, l2_7 = 0, l2_8 = 0, l2_9 = 0, l2_10 = 0;
	float l2_11 = 0, l2_12 = 0, l2_13 = 0, l2_14 = 0, l2_15 = 0, l2_16 = 0, l2_17 = 0, l2_18 = 0, l2_19 = 0;
	float m2_1 = 0, m2_2 = 0, m2_3 = 0, m2_4 = 0, m2_5 = 0, m2_6 = 0, m2_7 = 0, m2_8 = 0, m2_9 = 0, m2_10 = 0;
	float n2_1 = 0, n2_2 = 0, n2_3 = 0, n2_4 = 0, n2_5 = 0, n2_6 = 0, n2_7 = 0, n2_8 = 0, n2_9 = 0, n2_10 = 0;
	float o2_1 = 0, o2_2 = 0, o2_3 = 0, o2_4 = 0, o2_5 = 0, o2_6 = 0, o2_7 = 0, o2_8 = 0, o2_9 = 0, o2_10 = 0;
	float k3_1 = 0, k3_2 = 0, k3_3 = 0, k3_4 = 0, k3_5 = 0, k3_6 = 0, k3_7 = 0, k3_8 = 0, k3_9 = 0, k3_10 = 0;
	float k3_11 = 0, k3_12 = 0, k3_13 = 0, k3_14 = 0, k3_15 = 0, k3_16 = 0, k3_17 = 0, k3_18 = 0, k3_19 = 0;
	float l3_1 = 0, l3_2 = 0, l3_3 = 0, l3_4 = 0, l3_5 = 0, l3_6 = 0, l3_7 = 0, l3_8 = 0, l3_9 = 0, l3_10 = 0;
	float l3_11 = 0, l3_12 = 0, l3_13 = 0, l3_14 = 0, l3_15 = 0, l3_16 = 0, l3_17 = 0, l3_18 = 0, l3_19 = 0;
	float m3_1 = 0, m3_2 = 0, m3_3 = 0, m3_4 = 0, m3_5 = 0, m3_6 = 0, m3_7 = 0, m3_8 = 0, m3_9 = 0, m3_10 = 0;
	float n3_1 = 0, n3_2 = 0, n3_3 = 0, n3_4 = 0, n3_5 = 0, n3_6 = 0, n3_7 = 0, n3_8 = 0, n3_9 = 0, n3_10 = 0;
	float o3_1 = 0, o3_2 = 0, o3_3 = 0, o3_4 = 0, o3_5 = 0, o3_6 = 0, o3_7 = 0, o3_8 = 0, o3_9 = 0, o3_10 = 0;
	float k4_1 = 0, k4_2 = 0, k4_3 = 0, k4_4 = 0, k4_5 = 0, k4_6 = 0, k4_7 = 0, k4_8 = 0, k4_9 = 0, k4_10 = 0;
	float k4_11 = 0, k4_12 = 0, k4_13 = 0, k4_14 = 0, k4_15 = 0, k4_16 = 0, k4_17 = 0, k4_18 = 0, k4_19 = 0;
	float l4_1 = 0, l4_2 = 0, l4_3 = 0, l4_4 = 0, l4_5 = 0, l4_6 = 0, l4_7 = 0, l4_8 = 0, l4_9 = 0, l4_10 = 0;
	float l4_11 = 0, l4_12 = 0, l4_13 = 0, l4_14 = 0, l4_15 = 0, l4_16 = 0, l4_17 = 0, l4_18 = 0, l4_19 = 0;
	float m4_1 = 0, m4_2 = 0, m4_3 = 0, m4_4 = 0, m4_5 = 0, m4_6 = 0, m4_7 = 0, m4_8 = 0, m4_9 = 0, m4_10 = 0;
	float n4_1 = 0, n4_2 = 0, n4_3 = 0, n4_4 = 0, n4_5 = 0, n4_6 = 0, n4_7 = 0, n4_8 = 0, n4_9 = 0, n4_10 = 0;
	float o4_1 = 0, o4_2 = 0, o4_3 = 0, o4_4 = 0, o4_5 = 0, o4_6 = 0, o4_7 = 0, o4_8 = 0, o4_9 = 0, o4_10 = 0;
	float k23 = 0, k32 = 0, kadd1 = 0, kadd2 = 0, kadd = 0, kdel = 0;
	float l23 = 0, l32 = 0, ladd1 = 0, ladd2 = 0, ladd = 0, ldel = 0;
	float m23 = 0, m32 = 0, madd1 = 0, madd2 = 0, madd = 0, mdel = 0;
	float n23 = 0, n32 = 0, nadd1 = 0, nadd2 = 0, nadd = 0, ndel = 0;
	float o23 = 0, o32 = 0, oadd1 = 0, oadd2 = 0, oadd = 0, odel = 0;
	float omega_r1 = 0, omega_r2 = 0, omega_r3 = 0;

	ls = fpadd32f(lm,lls); //lm + lls;
	lr = fpadd32f(lm,llr); //lm + llr;
	alpha = fdiv32(rr,lr); // rr/lr;
	sigma1 = fdiv32(lm,lr);  //lm/lr
	sigma2 = fpmul32f(lm,sigma1); // lm*lm/lr
	sigma = fpsub32f(ls, sigma2); // ls - lm*lm/lr;
	beta1 = fpmul32f(sigma,lr);  // sigma*lr
	beta = fdiv32(lm,beta1); //  lm/(sigma*lr);

 //	mu = (3/2)*(poles/2)*(lm/(inertia*lr));	
	mu1 = fpmul32f(inertia,lr); //  inertia * lr
	mu2 = fdiv32(lm,mu1); //   lm/(inertia * lr)
	mu3 = fdiv32(poles,2);
	mu4 = fpmul32f(mu2,mu3);   //  (poles/2) * (lm/(inertia*lr))
	mu = fpmul32f(1.5,mu4);   //	mu = (3/2)*(poles/2)*(lm/(inertia*lr));


//	gamma = lm*lm*rr/(sigma*lr*lr)+rs/sigma;
	gamma1 = fdiv32(rs,sigma);  //rs/sigma
	gamma2 = fpmul32f(lr,lr);
	gamma3 = fpmul32f(sigma,gamma2);  // sigma*lr*lr
	gamma4 = fpmul32f(lm,lm);
	gamma5 = fpmul32f(gamma4,rr);  //  lm*lm*lr
	gamma6 = fdiv32(gamma5, gamma3);
	gamma = fpadd32f(gamma6, gamma1);   //	gamma = lm*lm*rr/(sigma*lr*lr)+rs/sigma;
	
	

//	cont = (3*poles*lm/(4*lr));
	cont1 = fpmul32f(4,lr);    //  4*lr
	cont2 = fpmul32f(poles,lm);
	cont3 = fpmul32f(3,cont2);   // 3*poles*lm
	cont = fdiv32(cont3,cont1);   // cont = (3*poles*lm/(4*lr))


while(1)
{	
// Inputs to the motor module
	// load_torque, vd, vq, iq, id, flq, fld, spd 


//same will be the inputs to motor
		id  = read_float32("in_data");
		iq  = read_float32("in_data");
		vq  = read_float32("in_data");
		vd  = read_float32("in_data");
		iq  = read_float32("in_data");
		flq  = read_float32("in_data");
		fld  = read_float32("in_data");	
		spd  = read_float32("in_data");
		load_torque  = read_float32("in_data");	
		//omega_m  = read_float32("in_data");



/*Output Data from motor  --
		id  = read_float32("in_data");
		iq  = read_float32("in_data");
		speed  = read_float32("in_data");
		speed_ref_temp  = read_float32("in_data");	
		omega_m  = read_float32("in_data");
				
*/

// put these 2 lines at the end
//		*iq_prev = *iq;
//		*id_prev = *id;
		delta = time_period;

// REMOVE ALL POINTERS -- NO NEED OF POINTERS AS ALL THE DATA WILL BE IN THIS MODULE. INPUTS AND OUTPUTS WILL ONLY BE SENT THROUGH PIPES.


//		k1 = -gamma*(iq_prev) - omega*(id_prev) + alpha*beta*(flq_prev) - beta*(spd_prev)*(fld_prev) + (vq)/sigma;
		k1_1 = fdiv32(vq,sigma); //--
		k1_2 = fpmul32f(fld_prev,spd_prev);
		k1_3 = -(fpmul32f(k1_2,beta));  //--
		k1_4 = fpmul32f(alpha,beta);
		k1_5 = fpmul32f(k1_4,flq_prev); // --
		k1_6 = -(fpmul32f(omega,id_prev)); // --
		k1_7 = -(fpmul32f(gamma,iq_prev)); // --
		k1_8 = fpadd32f(k1_1,k1_3); // ++
 		k1_9 = fpadd32f(k1_5,k1_6);
		k1_10 = fpadd32f(k1_9,k1_7); // ++
		k1 = fpadd32f(k1_10,k1_8); //+

//		k1 = -gamma*(iq_prev) - omega*(id_prev) + alpha*beta*(flq_prev) - beta*(spd_prev)*(fld_prev) + (vq)/sigma;
//		l1 =  omega*(iq_prev) - gamma*(id_prev) + alpha*beta*(fld_prev) + beta*(spd_prev)*(flq_prev) + (vd)/sigma;

		l1_1 = fdiv32(vd,sigma); //--
		l1_2 = fpmul32f(flq_prev,spd_prev);
		l1_3 = fpmul32f(l1_2,beta);  //--
		l1_4 = fpmul32f(alpha,beta);
		l1_5 = fpmul32f(l1_4,fld_prev); // --
		l1_6 = -(fpmul32f(gamma,id_prev)); // --
		l1_7 = (fpmul32f(omega,iq_prev)); // --
		l1_8 = fpadd32f(l1_1,l1_3); // ++
 		l1_9 = fpadd32f(l1_5,l1_6);
		l1_10 = fpadd32f(l1_9,l1_7); // ++
		l1 = fpadd32f(l1_10,l1_8); //+









//		m1 = alpha*lm*(*iq_prev) - alpha*(*flq_prev) - (omega-(*spd_prev))*(*fld_prev);

//		m1 = alpha*lm*iq_prev - alpha*flq_prev - omega*fld_prev + spd_prev*fld_prev;
		m1_1 = fpmul32f(spd_prev,fld_prev); //--
		m1_2 = -( fpmul32f(omega,fld_prev)); // --
		m1_3 = -(fpmul32f(alpha,flq_prev)); // --
		m1_4 = fpmul32f(alpha,lm);
		m1_5 = fpmul32f(iq_prev,m1_4); // --
		m1_6 = fpadd32f(m1_1,m1_2);
		m1_7 = fpadd32f(m1_3,m1_5);
		m1 = fpadd32f(m1_6,m1_7);
//		m1 = alpha*lm*iq_prev - alpha*flq_prev - omega*fld_prev + spd_prev*fld_prev;		
//		n1 = alpha*lm*(id_prev) + (omega-(*spd_prev))*(*flq_prev) - alpha*(*fld_prev);  
//		n1 = alpha*lm*id_prev + omega*flq_prev - alpha*fld_prev - spd_prev*flq_prev;  // simplified version


		n1_1 = -(fpmul32f(spd_prev,flq_prev)); //--
		n1_2 = -( fpmul32f(omega,fld_prev)); // --
		n1_3 = fpmul32f(omega,flq_prev); // --
		n1_4 = fpmul32f(alpha,lm);
		n1_5 = fpmul32f(id_prev,n1_4); // --
		n1_6 = fpadd32f(n1_1,n1_2);
		n1_7 = fpadd32f(n1_3,n1_5);
		n1 = fpadd32f(n1_6,n1_7);


//		o1 = [     {    cont*(fld_prev*iq_prev - flq_prev*id_prev)   }    -load_torque   ]/inertia;

		o1_1 = -(fpmul32f(flq_prev,id_prev));
		o1_2 = fpmul32f(fld_prev,iq_prev);
		o1_3 = fpadd32f(o1_1,o1_2);
		o1_4 = fpmul32f(cont,o1_3);
		o1_5 = fpsub32f(o1_4,load_torque);
		o1 = fdiv32(o1_5,inertia);


		delta_2 = fdiv32(delta,2);
		
//		k1 = -gamma*(iq_prev) - omega*(id_prev) + alpha*beta*(flq_prev) - beta*(spd_prev)*(fld_prev) + (vq)/sigma;
//		k2 = -gamma*((*iq_prev) + delta/2*k1) - omega*((*id_prev)+delta/2*l1) + alpha*beta*((*flq_prev)+delta/2*m1) - beta*((*spd_prev)+delta/2*o1)*((*fld_prev)+delta/2*n1) + (vq)/sigma;


//		k2 = -gamma*(iq_prev + delta/2*k1) - omega*(id_prev+delta/2*l1) + alpha*beta*(flq_prev+delta/2*m1) 
//			- beta*(spd_prev+delta/2*o1)*(fld_prev+delta/2*n1) + (vq)/sigma;

		// k1_1 = fdiv32(vq,sigma); //--
		k2_1 = fpmul32f(delta_2,n1);
		k2_2 = fpadd32f(fld_prev,k2_1);
		k2_3 = fpmul32f(delta_2,o1);
		k2_4 = fpadd32f(spd_prev,k2_3);
		k2_5 = fpmul32f(k2_2,k2_4);
		k2_6 = -(fpmul32f(beta,k2_5)); //--
		k2_7 = fpmul32f(delta_2,m1);
		k2_8 = fpadd32f(flq_prev,k2_7);
		k2_9 = fpmul32f(alpha,beta);
		k2_10 = fpmul32f(k2_9,k2_8); //--
		k2_11 = fpmul32f(delta_2,l1);
		k2_12 = fpadd32f(id_prev,k2_11);
		k2_13 = -(fpmul32f(omega,k2_12));//--
		k2_14 = fpmul32f(delta_2,k1);
		k2_15 = fpadd32f(iq_prev,k2_14);
		k2_16 = -(fpmul32f(gamma,k2_15)); //--
		k2_17 = fpadd32f(k1_1,k2_6);  //++
		k2_18 = fpadd32f(k2_10,k2_13);  //++
		k2_19 = fpadd32f(k2_16,k2_17);  //-+
		k2 = fpadd32f(k2_19,k2_18);
		
		

//		l2 = omega*((*iq_prev) + delta/2*k1) - gamma*((*id_prev) + delta/2*l1) + beta*((*spd_prev) + delta/2*o1)*((*flq_prev) + delta/2*m1) + alpha*beta*((*fld_prev) + delta/2*n1) + (vd)/sigma;

//		k2 = -gamma*(iq_prev + delta/2*k1) - omega*(id_prev+delta/2*l1) + alpha*beta*(flq_prev+delta/2*m1) 
//			- beta*(spd_prev+delta/2*o1)*(fld_prev+delta/2*n1) + (vq)/sigma;

//		l2 = omega*(iq_prev + delta/2*k1) - gamma*(id_prev + delta/2*l1) + alpha*beta*(fld_prev + delta/2*n1) 
//  			+ beta*(spd_prev + delta/2*o1)*(flq_prev + delta/2*m1)+ (vd)/sigma;

		// l1_1 = fdiv32(vd,sigma); //--
		l2_1 = fpmul32f(delta_2,m1);
		l2_2 = fpadd32f(flq_prev,l2_1);
		l2_3 = fpmul32f(delta_2,o1);
		l2_4 = fpadd32f(spd_prev,l2_3);
		l2_5 = fpmul32f(l2_2,l2_4);
		l2_6 = fpmul32f(beta,l2_5); //--
		l2_7 = fpmul32f(delta_2,n1);
		l2_8 = fpadd32f(fld_prev,l2_7);
		l2_9 = fpmul32f(alpha,beta);
		l2_10 = fpmul32f(l2_9,l2_8); //--
		l2_11 = fpmul32f(delta_2,l1);
		l2_12 = fpadd32f(id_prev,l2_11);
		l2_13 = -(fpmul32f(gamma,l2_12));//--
		l2_14 = fpmul32f(delta_2,k1);
		l2_15 = fpadd32f(iq_prev,l2_14);
		l2_16 = fpmul32f(omega,l2_15); //--
		l2_17 = fpadd32f(l1_1,l2_6);  //++
		l2_18 = fpadd32f(l2_10,l2_13);  //++
		l2_19 = fpadd32f(l2_16,l2_17);  //-+
		l2 = fpadd32f(l2_19,l2_18);
k1_1 = 0, k1_2 = 0,

//		m2 = alpha*lm*((*iq_prev) + delta/2*k1) - alpha*((*flq_prev) + delta/2*m1) - (omega-((*spd_prev) + delta/2*o1))*((*fld_prev) + delta/2*n1);


//		l2 = omega*(iq_prev + delta/2*k1) - gamma*(id_prev + delta/2*l1) + alpha*beta*(fld_prev + delta/2*n1) 
//  			+ beta*(spd_prev + delta/2*o1)*(flq_prev + delta/2*m1)+ (vd)/sigma;

//		m2 = alpha*lm*(iq_prev + delta/2*k1) - alpha*(flq_prev + delta/2*m1) - (omega-(spd_prev + delta/2*o1))*(fld_prev + delta/2*n1);


//		l2_1 = fpmul32f(delta_2,m1);
//		l2_2 = fpadd32f(flq_prev,l2_1);   
//		l2_14 = fpmul32f(delta_2,k1);
//		l2_15 = fpadd32f(iq_prev,l2_14);
//		l2_16 = fpmul32f(omega,l2_15);
//		l2_7 = fpmul32f(delta_2,n1);
//		l2_8 = fpadd32f(fld_prev,l2_7);   //*
//		l2_3 = fpmul32f(delta_2,o1);
//		l2_4 = fpadd32f(spd_prev,l2_3);


		m2_1 = -(fpmul32f(alpha,l2_2)); //--    - alpha*(flq_prev + delta/2*m1)
		m2_2 = fpmul32f(alpha,lm);
		m2_3 = fpmul32f(m2_2,l2_15);  //--      alpha*lm*(iq_prev + delta/2*k1)	
		m2_4 = fpsub32f(omega,l2_4);  //*
		m2_5 = -(fpmul32f(m2_4,l2_8));  //--    - (omega-(spd_prev + delta/2*o1))*(fld_prev + delta/2*n1)
		m2_6 = fpadd32f(m2_1,m2_3);
		m2 = fpadd32f(m2_6,m2_5);
		


//		n2 = alpha*lm*((*id_prev) + delta/2*l1) + (omega-((*spd_prev) + delta/2*o1))*((*flq_prev) + delta/2*m1) - alpha*((*fld_prev) + delta/2*n1);

//		l2 = omega*(iq_prev + delta/2*k1) - gamma*(id_prev + delta/2*l1) + alpha*beta*(fld_prev + delta/2*n1) 
//  			+ beta*(spd_prev + delta/2*o1)*(flq_prev + delta/2*m1)+ (vd)/sigma;

//		n2 = alpha*lm*(id_prev + delta/2*l1) + (omega-(spd_prev + delta/2*o1))*(flq_prev + delta/2*m1) - alpha*(fld_prev + delta/2*n1);


//		l2_1 = fpmul32f(delta_2,m1);
//		l2_2 = fpadd32f(flq_prev,l2_1);  // *
//		l2_3 = fpmul32f(delta_2,o1);
//		l2_4 = fpadd32f(spd_prev,l2_3);
//		l2_7 = fpmul32f(delta_2,n1);
//		l2_8 = fpadd32f(fld_prev,l2_7);
//		l2_11 = fpmul32f(delta_2,l1);
//		l2_12 = fpadd32f(id_prev,l2_11); //**

		n2_1 = fpsub32f(omega,l2_4); //*
		n2_2 = fpmul32f(n2_1,l2_2);    // --     	 + (omega-(spd_prev + delta/2*o1))*(flq_prev + delta/2*m1)
		n2_3 = -(fpmul32f(alpha,l2_8)); //--    	- alpha*(fld_prev + delta/2*n1)
		n2_4 = fpmul32f(alpha,lm); // **
		n2_5 = fpmul32f(n2_4,l2_12);  // --    		alpha*lm*(id_prev + delta/2*l1)
		n2_6 = fpadd32f(n2_2,n2_3);
		n2 = fpadd32f(n2_6,n2_5);
		

//		l2 = omega*(iq_prev + delta/2*k1) - gamma*(id_prev + delta/2*l1) + alpha*beta*(fld_prev + delta/2*n1) 
//  			+ beta*(spd_prev + delta/2*o1)*(flq_prev + delta/2*m1)+ (vd)/sigma;


//		o2 = ((cont*(((*fld_prev) + delta/2*n1)*((*iq_prev) + delta/2*k1) - ((*flq_prev) + delta/2*m1)*((*id_prev) + delta/2*l1)))-load_torque)/inertia;
//		o2 = ((cont*((fld_prev + delta/2*n1)*(iq_prev + delta/2*k1) - (flq_prev + delta/2*m1)*(id_prev + delta/2*l1)))-load_torque)/inertia;

//		l2_7 = fpmul32f(delta_2,n1);
//		l2_8 = fpadd32f(fld_prev,l2_7);  //$$
//		l2_15 = fpadd32f(iq_prev,l2_14);
//		l2_16 = fpmul32f(omega,l2_15); //$$
//		l2_1 = fpmul32f(delta_2,m1);   
//		l2_2 = fpadd32f(flq_prev,l2_1);  //$$$$
//		l2_11 = fpmul32f(delta_2,l1);
//		l2_12 = fpadd32f(id_prev,l2_11); //$$$$

		o2_1 = fpmul32f(l2_8,l2_16);     //--     (fld_prev + delta/2*n1)*(iq_prev + delta/2*k1)
		o2_2 = -(fpmul32f(l2_2,l2_12));     //--	  (flq_prev + delta/2*m1)*(id_prev + delta/2*l1)
		o2_3 = fpadd32f(o2_1,o2_2);
		o2_4 = fpmul32f(cont,o2_3);
		o2_5 = fpsub32f(o2_4,load_torque);
		o2 = fdiv32f(o2_5,inertia);


//		k3 = -gamma*((*iq_prev) + delta/2*k2) - omega*((*id_prev)+delta/2*l2) + alpha*beta*((*flq_prev)+delta/2*m2) - beta*((*spd_prev)+delta/2*o2)*((*fld_prev)+delta/2*n2) + (vq)/sigma;

//		k2 = -gamma*(iq_prev + delta/2*k1) - omega*(id_prev+delta/2*l1) + alpha*beta*(flq_prev+delta/2*m1) 
//			- beta*(spd_prev+delta/2*o1)*(fld_prev+delta/2*n1) + (vq)/sigma;

//		k3 = -gamma*(iq_prev + delta/2*k2) - omega*(id_prev+delta/2*l2) + alpha*beta*(flq_prev+delta/2*m2) 
//			- beta*(spd_prev+delta/2*o2)*(fld_prev+delta/2*n2) + (vq)/sigma;




		// k1_1 = fdiv32(vq,sigma); //--
		k3_1 = fpmul32f(delta_2,n2);
		k3_2 = fpadd32f(fld_prev,k3_1);
		k3_3 = fpmul32f(delta_2,o2);
		k3_4 = fpadd32f(spd_prev,k3_3);
		k3_5 = fpmul32f(k3_2,k3_4);
		k3_6 = -(fpmul32f(beta,k3_5)); //--
		k3_7 = fpmul32f(delta_2,m2);
		k3_8 = fpadd32f(flq_prev,k3_7);
		k3_9 = fpmul32f(alpha,beta);
		k3_10 = fpmul32f(k3_9,k3_8); //--
		k3_11 = fpmul32f(delta_2,l2);
		k3_12 = fpadd32f(id_prev,k3_11);
		k3_13 = -(fpmul32f(omega,k3_12));//--
		k3_14 = fpmul32f(delta_2,k2);
		k3_15 = fpadd32f(iq_prev,k3_14);
		k3_16 = -(fpmul32f(gamma,k3_15)); //--
		k3_17 = fpadd32f(k1_1,k3_6);  //++
		k3_18 = fpadd32f(k3_10,k3_13);  //++
		k3_19 = fpadd32f(k3_16,k3_17);  //-+
		k3 = fpadd32f(k3_19,k3_18);



//		l3 = omega*((*iq_prev) + delta/2*k2) - gamma*((*id_prev) + delta/2*l2) + beta*((*spd_prev) + delta/2*o2)*((*flq_prev) + delta/2*m2) 
//  			+ alpha*beta*((*fld_prev) + delta/2*n2) + (vd)/sigma;


//		l2 = omega*(iq_prev + delta/2*k1) - gamma*(id_prev + delta/2*l1) + alpha*beta*(fld_prev + delta/2*n1) 
//  			+ beta*(spd_prev + delta/2*o1)*(flq_prev + delta/2*m1) + (vd)/sigma;
//		l3 = omega*(iq_prev + delta/2*k2) - gamma*(id_prev + delta/2*l2)  + alpha*beta*(fld_prev + delta/2*n2)
//			+ beta*(spd_prev + delta/2*o2)*(flq_prev + delta/2*m2) + (vd)/sigma;

		// l1_1 = fdiv32(vd,sigma); //--
		l3_1 = fpmul32f(delta_2,m2);
		l3_2 = fpadd32f(flq_prev,l3_1);
		l3_3 = fpmul32f(delta_2,o2);
		l3_4 = fpadd32f(spd_prev,l3_3);
		l3_5 = fpmul32f(l3_2,l3_4);
		l3_6 = fpmul32f(beta,l3_5); //--
		l3_7 = fpmul32f(delta_2,n2);
		l3_8 = fpadd32f(fld_prev,l3_7);
		l3_9 = fpmul32f(alpha,beta);
		l3_10 = fpmul32f(l3_9,l3_8); //--
		l3_11 = fpmul32f(delta_2,l2);
		l3_12 = fpadd32f(id_prev,l3_11);
		l3_13 = -(fpmul32f(gamma,l3_12));//--
		l3_14 = fpmul32f(delta_2,k2);
		l3_15 = fpadd32f(iq_prev,l3_14);
		l3_16 = fpmul32f(omega,l3_15); //--
		l3_17 = fpadd32f(l1_1,l3_6);  //++
		l3_18 = fpadd32f(l3_10,l3_13);  //++
		l3_19 = fpadd32f(l3_16,l3_17);  //-+
		l3 = fpadd32f(l3_19,l3_18);


//		m3 = alpha*lm*((*iq_prev) + delta/2*k2) - alpha*((*flq_prev) + delta/2*m2) - (omega-((*spd_prev) + delta/2*o2))*((*fld_prev) + delta/2*n2);

//		l3 = omega*(iq_prev + delta/2*k2) - gamma*(id_prev + delta/2*l2)  + alpha*beta*(fld_prev + delta/2*n2)
//			+ beta*(spd_prev + delta/2*o2)*(flq_prev + delta/2*m2) + (vd)/sigma;

//		m3 = alpha*lm*(iq_prev + delta/2*k2) - alpha*(flq_prev + delta/2*m2) - (omega-(spd_prev + delta/2*o2))*(fld_prev + delta/2*n2);




//		l3_1 = fpmul32f(delta_2,m2);
//		l3_2 = fpadd32f(flq_prev,l3_1);
//		l3_14 = fpmul32f(delta_2,k2);
//		l3_15 = fpadd32f(iq_prev,l3_14);
//		l3_3 = fpmul32f(delta_2,o2);
//		l3_4 = fpadd32f(spd_prev,l3_3);
//		l3_7 = fpmul32f(delta_2,n2);
//		l3_8 = fpadd32f(fld_prev,l3_7);


		m3_1 = -(fpmul32f(alpha,l3_2)); //--    - alpha*(flq_prev + delta/2*m2)
		//m3_2 = fpmul32f(alpha,lm);   == m2_2
		m3_3 = fpmul32f(m2_2,l3_15);  //--      alpha*lm*(iq_prev + delta/2*k2)	
		m3_4 = fpsub32f(omega,l3_4);  //*
		m3_5 = -(fpmul32f(m3_4,l3_8));  //--    - (omega-(spd_prev + delta/2*o2))*(fld_prev + delta/2*n2)
		m3_6 = fpadd32f(m3_1,m3_3);
		m3 = fpadd32f(m3_6,m3_5);
		




//		n3 = alpha*lm*((*id_prev) + delta/2*l2) + (omega-((*spd_prev) + delta/2*o2))*((*flq_prev) + delta/2*m2) - alpha*((*fld_prev) + delta/2*n2);

//		n3 = alpha*lm*(id_prev + delta/2*l2) + (omega-(spd_prev + delta/2*o2))*(flq_prev + delta/2*m2) - alpha*(fld_prev + delta/2*n2);
//		n2 = alpha*lm*(id_prev + delta/2*l1) + (omega-(spd_prev + delta/2*o1))*(flq_prev + delta/2*m1) - alpha*(fld_prev + delta/2*n1);

//		l3_3 = fpmul32f(delta_2,o2);
//		l3_4 = fpadd32f(spd_prev,l3_3);
//		l3_7 = fpmul32f(delta_2,n2);
//		l3_8 = fpadd32f(fld_prev,l3_7);

		n3_1 = fpsub32f(omega,l3_4); //*
		n3_2 = fpmul32f(n3_1,l3_2);    // --     	 + (omega-(spd_prev + delta/2*o1))*(flq_prev + delta/2*m1)
		n3_3 = -(fpmul32f(alpha,l3_8)); //--    	- alpha*(fld_prev + delta/2*n1)
		//n3_4 = fpmul32f(alpha,lm); //  == m2_2
		n3_5 = fpmul32f(m2_2,l3_12);  // --    		alpha*lm*(id_prev + delta/2*l1)
		n3_6 = fpadd32f(n3_2,n3_3);
		n3 = fpadd32f(n3_6,n3_5);
		

//		o3 = ((cont*(((*fld_prev) + delta/2*n2)*((*iq_prev) + delta/2*k2) - ((*flq_prev) + delta/2*m2)*((*id_prev) + delta/2*l2)))-load_torque)/inertia;
//		o3 = ((cont*((fld_prev + delta/2*n2)*(iq_prev + delta/2*k2) - (flq_prev + delta/2*m2)*(id_prev + delta/2*l2)))-load_torque)/inertia;
//		o2 = ((cont*((fld_prev + delta/2*n1)*(iq_prev + delta/2*k1) - (flq_prev + delta/2*m1)*(id_prev + delta/2*l1)))-load_torque)/inertia;


//		l3_7 = fpmul32f(delta_2,n2);
//		l3_8 = fpadd32f(fld_prev,l3_7); //$$
//		l3_15 = fpadd32f(iq_prev,l3_14);
//		l3_16 = fpmul32f(omega,l3_15); //$$
//		l3_1 = fpmul32f(delta_2,m2);
//		l3_2 = fpadd32f(flq_prev,l3_1);  //$$$$
//		l3_11 = fpmul32f(delta_2,l2);
//		l3_12 = fpadd32f(id_prev,l3_11); //$$$$

		
		o3_1 = fpmul32f(l3_8,l3_16);     //--     (fld_prev + delta/2*n1)*(iq_prev + delta/2*k1)
		o3_2 = -(fpmul32f(l3_2,l3_12));     //--	  (flq_prev + delta/2*m1)*(id_prev + delta/2*l1)
		o3_3 = fpadd32f(o3_1,o3_2);
		o3_4 = fpmul32f(cont,o3_3);
		o3_5 = fpsub32f(o3_4,load_torque);
		o3 = fdiv32f(o3_5,inertia);


//		k4 = -gamma*((*iq_prev) + delta*k3) - omega*((*id_prev)+delta*l3) + alpha*beta*((*flq_prev)+delta*m3) - beta*((*spd_prev)+delta*o3)*((*fld_prev)+delta*n3) + (vq)/sigma;

//		k2 = -gamma*(iq_prev + delta/2*k1) - omega*(id_prev+delta/2*l1) + alpha*beta*(flq_prev+delta/2*m1) 
//			- beta*(spd_prev+delta/2*o1)*(fld_prev+delta/2*n1) + (vq)/sigma;

//		k4 = -gamma*(iq_prev + delta*k3) - omega*(id_prev+delta*l3) + alpha*beta*(flq_prev+delta*m3) 
//			- beta*(spd_prev+delta*o3)*(fld_prev+delta*n3) + (vq)/sigma;


		// k1_1 = fdiv32(vq,sigma); //--
		k4_1 = fpmul32f(delta,n3);
		k4_2 = fpadd32f(fld_prev,k4_1);
		k4_3 = fpmul32f(delta,o3);
		k4_4 = fpadd32f(spd_prev,k4_3);
		k4_5 = fpmul32f(k4_2,k4_4);
		k4_6 = -(fpmul32f(beta,k4_5)); //--
		k4_7 = fpmul32f(delta,m3);
		k4_8 = fpadd32f(flq_prev,k4_7);
		k4_9 = fpmul32f(alpha,beta);
		k4_10 = fpmul32f(k4_9,k4_8); //--
		k4_11 = fpmul32f(delta,l3);
		k4_12 = fpadd32f(id_prev,k4_11);
		k4_13 = -(fpmul32f(omega,k4_12));//--
		k4_14 = fpmul32f(delta,k3);
		k4_15 = fpadd32f(iq_prev,k4_14);
		k4_16 = -(fpmul32f(gamma,k4_15)); //--
		k4_17 = fpadd32f(k1_1,k4_6);  //++
		k4_18 = fpadd32f(k4_10,k4_13);  //++
		k4_19 = fpadd32f(k4_16,k4_17);  //-+
		k4 = fpadd32f(k4_19,k4_18);



//		l4 = omega*((*iq_prev) + delta*k3) - gamma*((*id_prev) + delta*l3) + beta*((*spd_prev) + delta*o3)*((*flq_prev) + delta*m3) + alpha*beta*((*fld_prev) + delta*n3) + (vd)/sigma;
//		l4 = omega*(iq_prev + delta*k3) - gamma*(id_prev + delta*l3) + alpha*beta*(fld_prev + delta*n3)
//			+ beta*(spd_prev + delta*o3)*(flq_prev + delta*m3) + (vd)/sigma;

//		l2 = omega*(iq_prev + delta/2*k1) - gamma*(id_prev + delta/2*l1) + alpha*beta*(fld_prev + delta/2*n1) 
//  			+ beta*(spd_prev + delta/2*o1)*(flq_prev + delta/2*m1)+ (vd)/sigma;

		// l1_1 = fdiv32(vd,sigma); //--
		l4_1 = fpmul32f(delta,m3);
		l4_2 = fpadd32f(flq_prev,l4_1);
		l4_3 = fpmul32f(delta,o3);
		l4_4 = fpadd32f(spd_prev,l4_3);
		l4_5 = fpmul32f(l4_2,l4_4);
		l4_6 = fpmul32f(beta,l4_5); //--
		l4_7 = fpmul32f(delta,n3);
		l4_8 = fpadd32f(fld_prev,l4_7);
		l4_9 = fpmul32f(alpha,beta);
		l4_10 = fpmul32f(l4_9,l4_8); //--
		l4_11 = fpmul32f(delta,l3);
		l4_12 = fpadd32f(id_prev,l4_11);
		l4_13 = -(fpmul32f(gamma,l4_12));//--
		l4_14 = fpmul32f(delta,k3);
		l4_15 = fpadd32f(iq_prev,l4_14);
		l4_16 = fpmul32f(omega,l4_15); //--
		l4_17 = fpadd32f(l1_1,l4_6);  //++
		l4_18 = fpadd32f(l4_10,l4_13);  //++
		l4_19 = fpadd32f(l4_16,l4_17);  //-+
		l4 = fpadd32f(l4_19,l4_18);


//		m4 = alpha*lm*((*iq_prev) + delta*k3) - alpha*((*flq_prev) + delta*m3) - (omega-((*spd_prev) + delta*o3))*((*fld_prev) + delta*n3);
//		m2 = alpha*lm*(iq_prev + delta/2*k1) - alpha*(flq_prev + delta/2*m1) - (omega-(spd_prev + delta/2*o1))*(fld_prev + delta/2*n1);
//		m4 = alpha*lm*(iq_prev + delta*k3)   - alpha*(flq_prev + delta*m3)   - (omega-(spd_prev + delta*o3))*(fld_prev + delta*n3);

//		l4 = omega*(iq_prev + delta*k3) - gamma*(id_prev + delta*l3) + alpha*beta*(fld_prev + delta*n3)
//			+ beta*(spd_prev + delta*o3)*(flq_prev + delta*m3) + (vd)/sigma;






		m4_1 = -(fpmul32f(alpha,l4_2)); //--    - alpha*(flq_prev + delta*m3)
		//m3_2 = fpmul32f(alpha,lm);   == m2_2
		m4_3 = fpmul32f(m2_2,l4_15);  //--      alpha*lm*(iq_prev + delta*k3)	
		m4_4 = fpsub32f(omega,l4_4);  //*
		m4_5 = -(fpmul32f(m4_4,l4_8));  //--    - (omega-(spd_prev + delta*o3))*(fld_prev + delta*n3)
		m4_6 = fpadd32f(m4_1,m4_3);
		m4 = fpadd32f(m3_6,m4_5);
		




//		n4 = alpha*lm*((*id_prev) + delta*l3) + (omega-((*spd_prev) + delta*o3))*((*flq_prev) + delta*m3) - alpha*((*fld_prev) + delta*n3);
//		n3 = alpha*lm*(id_prev + delta/2*l2) + (omega-(spd_prev + delta/2*o2))*(flq_prev + delta/2*m2) - alpha*(fld_prev + delta/2*n2);
//		n4 = alpha*lm*(id_prev + delta*l3) + (omega-(spd_prev + delta*o3))*(flq_prev + delta*m3) - alpha*(fld_prev + delta*n3);


		n4_1 = fpsub32f(omega,l4_4); //*
		n4_2 = fpmul32f(n4_1,l4_2);    // --     	 + (omega-(spd_prev + delta*o3))*(flq_prev + delta*m3)
		n4_3 = -(fpmul32f(alpha,l4_8)); //--    	- alpha*(fld_prev + delta*n3)
		//n3_4 = fpmul32f(alpha,lm); //  == m2_2
		n4_5 = fpmul32f(m2_2,l4_12);  // --    		alpha*lm*(id_prev + delta*l3)
		n4_6 = fpadd32f(n4_2,n4_3);
		n4 = fpadd32f(n4_6,n4_5);


//		o4 = ((cont*(((*fld_prev) + delta*n3)*((*iq_prev) + delta*k3) - ((*flq_prev) + delta*m3)*((*id_prev) + delta*l3)))-load_torque)/inertia;
//		o3 = ((cont*((fld_prev + delta/2*n2)*(iq_prev + delta/2*k2) - (flq_prev + delta/2*m2)*(id_prev + delta/2*l2)))-load_torque)/inertia;
//		o4 = ((cont*((fld_prev + delta*n3)*(iq_prev + delta*k3) - (flq_prev + delta*m3)*(id_prev + delta*l3)))-load_torque)/inertia;

		o4_1 = fpmul32f(l4_8,l4_16);     //--     (fld_prev + delta*n3)*(iq_prev + delta*k3)
		o4_2 = -(fpmul32f(l4_2,l4_12));     //--	  (flq_prev + delta*m3)*(id_prev + delta*l3)
		o4_3 = fpadd32f(o4_1,o4_2);
		o4_4 = fpmul32f(cont,o4_3);
		o4_5 = fpsub32f(o4_4,load_torque);
		o4 = fdiv32f(o4_5,inertia);



//Outputs from motor module
	// iq, id, flq, fld, spd 

// outputs


		k23 = fpadd32f(k2,k3);
		k32 = fpmul32f(k23,2);
		kadd1 = fpadd32f(k32,k1);
		kadd2 = fpadd32f(kadd1,k4);
		kadd = fdiv32f(kadd2,6);
		kdel = fpmul32f(delta,kadd);
		iq = fpadd32f(iq_prev,kdel);
//		iq = iq_prev + delta*(k1 + 2*k2 + 2*k3 + k4)/6;



//		id = id_prev + delta*(l1 + 2*l2 + 2*l3 + l4)/6;


		l23 = fpadd32f(l2,l3);
		l32 = fpmul32f(l23,2);
		ladd1 = fpadd32f(l32,l1);
		ladd2 = fpadd32f(ladd1,l4);
		ladd = fdiv32f(ladd2,6);
		ldel = fpmul32f(delta,ladd);
		id = fpadd32f(id_prev,ldel);


//		flq = flq_prev + delta*(m1 + 2*m2 + 2*m3 + m4)/6;

		m23 = fpadd32f(m2,m3);
		m32 = fpmul32f(m23,2);
		madd1 = fpadd32f(m32,m1);
		madd2 = fpadd32f(madd1,m4);
		madd = fdiv32f(madd2,6);
		mdel = fpmul32f(delta,madd);
		flq = fpadd32f(flq_prev,mdel);






//		fld = fld_prev + delta*(n1 + 2*n2 + 2*n3 + n4)/6;

		n23 = fpadd32f(n2,n3);
		n32 = fpmul32f(n23,2);
		nadd1 = fpadd32f(n32,n1);
		nadd2 = fpadd32f(nadd1,n4);
		nadd = fdiv32f(nadd2,6);
		ndel = fpmul32f(delta,nadd);
		fld = fpadd32f(fld_prev,ndel);



//		spd = spd_prev + delta*(o1 + 2*o2 + 2*o3 + o4)/6;	// output to the vehicle module
		o23 = fpadd32f(o2,o3);
		o32 = fpmul32f(o23,2);
		oadd1 = fpadd32f(o32,o1);
		oadd2 = fpadd32f(oadd1,o4);
		oadd = fdiv32f(oadd2,6);
		odel = fpmul32f(delta,oadd);
		spd = fpadd32f(spd_prev,odel);





//		torque = cont*((*iq)*(*fld) - (*id)*(*flq));  // Output to the vehicle module
		torque = cont*(iq*fld - id*flq);  // Output to the vehicle module

//Pipe out_motor for motor -> vector control daemon

		
// Make a squre root function 		
		f_rotor = sqrt(  (fld*fld) + (flq*flq)  );
		
//		omega_r = (lm * (*iq))/((lr/rr)*f_rotor);
//		omega_r = (lm * iq)/((lr/rr)*f_rotor);
		omega_r1 = fpmul32f(lm,iq);
		omega_r2 = fdiv32f(lr,rr);
		omega_r3 = fpmul32f(omega_r2,f_rotor);
		omega_r = fdiv32f(omega_r1,omega_r3);
		
		if (iq == 0)
			omega_r = 0;
		else
			omega_r = omega_r; 
			
		//omega_r = 0;
		omega = omega_r;
		
	
		
//  *********** What about below line of code?		
//		theta = theta_prev + (+omega+*spd) * delta;
		theta = theta_prev + (omega*spd) * delta;
		
		time = time + delta;
		//*iq_prev = *iq;
		//*id_prev = *id;
		flq_prev = flq;
		fld_prev = fld;
		spd_prev = spd;
		//*theta_prev_prev = *theta_prev;
		theta_cos = theta;theta_sin = theta;
		theta_prev = theta;
	}
}

