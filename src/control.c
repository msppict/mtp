#include "prog.h"
void control()
{ 

	float id = 0; float iq = 0; float torque_ref = 0; float flux_ref = 0; float speed = 0; 
	speed  = read_float32("in_data"); //current speed from vehicle model
	speed_req  = read_float32("in_data"); //required speed from vehicle model
	speed  = read_float32("in_data");
	speed_ref_temp  = read_float32("in_data");



}	
