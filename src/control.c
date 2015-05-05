#include "prog.h"
void vector_control_daemon(){

	float id = 0; float iq = 0; float torque_ref = 0; float flux_ref = 0; float speed = 0; 
	acc  = read_float32("in_data"); //acceleration from User
	iq  = read_float32("in_data");
	speed  = read_float32("in_data");
	speed_ref_temp  = read_float32("in_data");

	//Rolling Resistance 
	F_rolling = mu*mass*g; //mu= rolling resistance coefficient

	//Aerodynamic resistance , rho = air density A = frontal area, Cd = drag coefficient
	F_aero = (0.5)*rho*A*Cd*speed; 

	//Hill climbing force  psi = elevation angle
	F_hill = mass*g*(sin(psi));

	// Linear accelerating force
	F_linear = mass*acceleration;   //define acceleration in upper section

	// Angular acceleration force/torque
	
