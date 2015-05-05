#include "prog.h"


/*
float fpadd32f(float x, float y);
float fpsub32f(float x, float y);
uint32_t fpadd32fi(uint32_t x, uint32_t y);
uint32_t fpsub32fi(uint32_t x, uint32_t y);
float fpmul32f(float x, float y);
float fdiv32(float a, float b);
*/


// Read this 

// http://scholar.lib.vt.edu/theses/available/etd-5440202339731121/unrestricted/CHAP3_DOC.pdf

// https://www.physicsforums.com/threads/how-to-calculate-max-velocity-of-a-vehicle-from-torque-and-rpm.355334/

// http://ev-bg.com/wordpress1/wp-content/uploads/2011/08/electric-vehicle-technology-explained-2003-j-larminie.pdf

void vehicle()
{
	float F_tr = 0,    // tractive force - to be found out in this function
	mu_r = 0.3, 	   // coeffiient of friction
	R = 0.5,	   // Radius of the wheel
	mass_chasis = 600,
	mass_battery = 400,
	mass = fpadd32f(mass_chasis,mass_battery), 	   // mass of the vehicle in Kg
	mass_linear = fpmul32f(mass,1.1),
	distance = 0,	
	g = 9.8, 	   // gravitational acceleration
	rho = 0.1, 	   // air drag coefficient
	A = 2, 		   // frontal area of car in m^2
	v_prev = 0, 	   // previous velocity of car
	v_new = 0, 	   // new velocity of car
	theta = 0, 	   // angle of elevation
	acc = 0, 	   // acceleration
	acc_prev = 0,	   // acceleration value in previous iteration
	Cd = 0.3, 	   // air drag coefficient
	brake = 0,
	brake_prev = 0,
	speed_max = 120,
	speed_req = 0,
	speed = 0,	   // current speed
 	speed_prev = 0,	   // previous iteration speed
	F_rr = 0,	   // rolling resistance force
	F_hill = 0,	   // Uphill/downhill force
	F_aero = 0,	   // Aerodynamic resistance force
	F_linear = 0,	   // Linear force = mass * acceleration
	F_brake = 0,	   // Braking force
	F_motor = 0,
	elev = 0,	   // sine of the elevation angle -- taken from testbench for now
	torque = 0,
	speed_sq = 0,
	
	delta = 10e-6;

	float temp0 = 0, temp1 = 0, temp2 = 0, temp3 = 0, tempx = 0, tempy = 0, tempz = 0;
	
while(1)
{


	torque = read_float32("in_data1");		//Get the  current torque from motor 	  -- motor module
	float acc = read_float32("in_data2");		//Get the acceleration from User  -- testbench
	float brake = read_float32("in_data2");		//Get the brake value from User	  -- testbench	
	float elev = read_float32("in_data2");		//Get sin of elevation from testbench -- testbench	**Need to add Sin and Cos in this module 
	

	// brake range 0-10. Max brake value = 10.
	if(brake < 0)
	brake = 0;
	else 
	{
	if(brake > 10)	brake = 10;
	else brake = brake;
	}


	// acc range 0-10. Max acc value = 10.
	if(acc < 0)
	acc = 0;
	else 
	{
	if(acc > 10)	acc = 10;
	else acc = acc;
	}


	//getting the max speed depending upon the brake position
	temp1 = fdiv32(brake,10);
	speed_max = fpmul32(temp1,120);
	
	

	//Rolling Resistance 
	temp1 = fpmul32(mu_r,mass); //mu= rolling resistance coefficient
	F_rr = fpmul32(temp1,g); 

	//Aerodynamic resistance , rho = air density A = frontal area, Cd = drag coefficient
	//F_aero = (0.5)*rho*A*Cd*speed^2; 

	temp1 = fpmul32(rho,A);
	temp2 = fpmul32(Cd,temp1);
	temp3 = fpmul32(temp2,0.5);
	speed_sq = fpmul32(speed_prev,speed_prev);
	//temp4 = fpmul32(temp3,speed_prev);
	F_aero = fpmul32(temp3,speed_sq) ; 		// v^2 term

	//Hill climbing force  psi = elevation angle	//F_hill = mass*g*(sin(psi));
	temp1 = fpmul32(mass,g);
	F_hill = fpmul32(temp1,elev);

	// Linear accelerating force
//	F_linear = fpmul32(mass,acceleration); 	// dv/dt

	// Braking force
	temp1 = fpmul32(brake, -0.1);
	F_brake = fpmul32(temp1, 555555.55); //check project book for calculations 

	// Force due to motor torque
	F_motor = fdiv32(torque,R);
	
	// F_motor = F_rr + F_aero + F_hill + F_linear + F_brake
	
	temp1 = fpsub32f(F_motor, F_hill);
	temp0 = fpsub32f(temp1,F_brake);	
	temp2 = fpsub32f(temp0, F_rr); 		
	tempx = fpsub32f(temp2, F_aero);
	temp1 = fdiv32(tempx,mass_linear);
	temp2 = fdiv32(F_aero,mass_linear);
	tempy = fpsub32f(temp1,temp2);
	tempz = fpmul32(tempy,delta);
	speed = fpadd32f(tempz,speed_prev);

	speed_req = speed + fpmul32(delta,acc);  // --acc = input from the user

	distance = distance + fpmul32(speed_prev,delta);

// Ouptputs to various modules
	write_float32("out_data",speed);
	write_float32("out_data",speed_req);
//	write_float32("out_data",theta);
//	write_float32("out_data",flux_rotor);


}
}
