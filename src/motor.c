/*
float fpadd32f(float x, float y);
float fpsub32f(float x, float y);
uint32_t fpadd32fi(uint32_t x, uint32_t y);
uint32_t fpsub32fi(uint32_t x, uint32_t y);
float fpmul32f(float x, float y);
float fdiv32(float a, float b);
*/


void im_zep()   
{
float iq, iq_prev, id, id_prev, flq, flq_prev, fld, fld_prev, spd_prev, vd, vq, torque, load_torque, time, spd, theta, theta_prev, theta_sin,theta_cos;


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
		time_period = 50e-6;
                                                   
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

// Inputs to the motor module
	// load_torque, vd, vq, iq, id, flq, fld, spd 

//Outputs from motor module
	// iq, id, flq, fld, spd 
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
	
        if ( 0.0 == time ) {           
                  
		iq_prev = 0.0;
		id_prev = 0.0;
		flq_prev = 0.0;
		fld_prev = 0.0;
		spd_prev = 0.0;
		iq = 0.0;
		id = 0.0;
		flq = 0.0;
		fld = 0.0;
		spd = 0.0;
		torque = 0.0;
		iq_prev =  iq;
		id_prev =  id;
		flq_prev =  flq;
		fld_prev =  fld;
		spd_prev =  spd;
		delta = time_period;
		time =  time + delta;
		theta = 0;
		theta_prev = 0;
		theta_sin = 0;
		theta_cos = 0;
		omega = 0;
        }
        else {       
		*iq_prev = *iq;
		*id_prev = *id;
		delta = time_period;

		k1 = -gamma*(*iq_prev) - omega*(*id_prev) + alpha*beta*(*flq_prev) - beta*(*spd_prev)*(*fld_prev) + (vq)/sigma;
		l1 = omega*(*iq_prev) - gamma*(*id_prev) + beta*(*spd_prev)*(*flq_prev) + alpha*beta*(*fld_prev) + (vd)/sigma;
		m1 = alpha*lm*(*iq_prev) - alpha*(*flq_prev) - (omega-(*spd_prev))*(*fld_prev);
		n1 = alpha*lm*(*id_prev) + (omega-(*spd_prev))*(*flq_prev) - alpha*(*fld_prev);
		o1 = ((cont*((*fld_prev)*(*iq_prev) - (*flq_prev)*(*id_prev)))-load_torque)/inertia;

		k2 = -gamma*((*iq_prev) + delta/2*k1) - omega*((*id_prev)+delta/2*l1) + alpha*beta*((*flq_prev)+delta/2*m1) - beta*((*spd_prev)+delta/2*o1)*((*fld_prev)+delta/2*n1) + (vq)/sigma;
		l2 = omega*((*iq_prev) + delta/2*k1) - gamma*((*id_prev) + delta/2*l1) + beta*((*spd_prev) + delta/2*o1)*((*flq_prev) + delta/2*m1) + alpha*beta*((*fld_prev) + delta/2*n1) + (vd)/sigma;
		m2 = alpha*lm*((*iq_prev) + delta/2*k1) - alpha*((*flq_prev) + delta/2*m1) - (omega-((*spd_prev) + delta/2*o1))*((*fld_prev) + delta/2*n1);
		n2 = alpha*lm*((*id_prev) + delta/2*l1) + (omega-((*spd_prev) + delta/2*o1))*((*flq_prev) + delta/2*m1) - alpha*((*fld_prev) + delta/2*n1);
		o2 = ((cont*(((*fld_prev) + delta/2*n1)*((*iq_prev) + delta/2*k1) - ((*flq_prev) + delta/2*m1)*((*id_prev) + delta/2*l1)))-load_torque)/inertia;

		k3 = -gamma*((*iq_prev) + delta/2*k2) - omega*((*id_prev)+delta/2*l2) + alpha*beta*((*flq_prev)+delta/2*m2) - beta*((*spd_prev)+delta/2*o2)*((*fld_prev)+delta/2*n2) + (vq)/sigma;
		l3 = omega*((*iq_prev) + delta/2*k2) - gamma*((*id_prev) + delta/2*l2) + beta*((*spd_prev) + delta/2*o2)*((*flq_prev) + delta/2*m2) + alpha*beta*((*fld_prev) + delta/2*n2) + (vd)/sigma;
		m3 = alpha*lm*((*iq_prev) + delta/2*k2) - alpha*((*flq_prev) + delta/2*m2) - (omega-((*spd_prev) + delta/2*o2))*((*fld_prev) + delta/2*n2);
		n3 = alpha*lm*((*id_prev) + delta/2*l2) + (omega-((*spd_prev) + delta/2*o2))*((*flq_prev) + delta/2*m2) - alpha*((*fld_prev) + delta/2*n2);
		o3 = ((cont*(((*fld_prev) + delta/2*n2)*((*iq_prev) + delta/2*k2) - ((*flq_prev) + delta/2*m2)*((*id_prev) + delta/2*l2)))-load_torque)/inertia;

		k4 = -gamma*((*iq_prev) + delta*k3) - omega*((*id_prev)+delta*l3) + alpha*beta*((*flq_prev)+delta*m3) - beta*((*spd_prev)+delta*o3)*((*fld_prev)+delta*n3) + (vq)/sigma;
		l4 = omega*((*iq_prev) + delta*k3) - gamma*((*id_prev) + delta*l3) + beta*((*spd_prev) + delta*o3)*((*flq_prev) + delta*m3) + alpha*beta*((*fld_prev) + delta*n3) + (vd)/sigma;
		m4 = alpha*lm*((*iq_prev) + delta*k3) - alpha*((*flq_prev) + delta*m3) - (omega-((*spd_prev) + delta*o3))*((*fld_prev) + delta*n3);
		n4 = alpha*lm*((*id_prev) + delta*l3) + (omega-((*spd_prev) + delta*o3))*((*flq_prev) + delta*m3) - alpha*((*fld_prev) + delta*n3);
		o4 = ((cont*(((*fld_prev) + delta*n3)*((*iq_prev) + delta*k3) - ((*flq_prev) + delta*m3)*((*id_prev) + delta*l3)))-load_torque)/inertia;


// outputs
		*iq = (*iq_prev) + delta*(k1 + 2*k2 + 2*k3 + k4)/6;
		*id = (*id_prev) + delta*(l1 + 2*l2 + 2*l3 + l4)/6;
		*flq = (*flq_prev) + delta*(m1 + 2*m2 + 2*m3 + m4)/6;
		*fld = (*fld_prev) + delta*(n1 + 2*n2 + 2*n3 + n4)/6;
		*spd = (*spd_prev) + delta*(o1 + 2*o2 + 2*o3 + o4)/6;	// output to the vehicle module
		*torque = cont*((*iq)*(*fld) - (*id)*(*flq));  // Output to the vehicle module

//Pipe out_motor for motor -> vector control daemon

		
		
		f_rotor = sqrt(((*fld)*(*fld)) + ((*flq)*(*flq)));
		
		omega_r = (lm * (*iq))/((lr/rr)*f_rotor);
		if (*iq == 0)
			omega_r = 0;
		else
			omega_r = omega_r; 
			
		//omega_r = 0;
		omega = omega_r;
		
	
		
		
		*theta = *theta_prev + (+omega+*spd) * delta;
		//*theta = omega * (*time);
		//*theta = atan((*flq)/((*fld)));
		
		//*theta = *theta + (omega_r + *spd) * delta;
		/*if(*flq == 0){
			*theta_sin = 0;
			*theta_cos = 0;
		//	*theta = 0;
		}
		else{
			*theta_cos = acos((*fld)/(f_rotor));
			*theta_sin = asin((*flq)/(f_rotor));
		//	*theta = atan((*flq)/((*fld)));
		}
		
		//if(*fld == 0)
		//	*theta_cos = 0;
		//else
			
		
		*/
		
		*time = *time + delta;
		//*iq_prev = *iq;
		//*id_prev = *id;
		*flq_prev = *flq;
		*fld_prev = *fld;
		*spd_prev = *spd;
		//*theta_prev_prev = *theta_prev;
		*theta_cos = *theta;*theta_sin = *theta;
		*theta_prev = *theta;
	}
}

