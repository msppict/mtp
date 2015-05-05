#include "prog.h"


/*
float fpadd32f(float x, float y);
float fpsub32f(float x, float y);
uint32_t fpadd32fi(uint32_t x, uint32_t y);
uint32_t fpsub32fi(uint32_t x, uint32_t y);
float fpmul32f(float x, float y);
float fdiv32(float a, float b);
*/



void control()
{ 
float	speed = 0,
	speed_req = 0,
	error = 0,
	prev_error = 0,
	Kp = 1,
	Ki = 0.5,
	Kd = 1,
	integral = 0,
	derivative = 0,
	proportional = 0,
	delta = 5e-6,
	output = 0;	

float temp1 = 0, temp2 = 0, temp_d = 0, temp_d1 = 0, temp_prop = 0, temp_i = 0, temp_output, temp = 0;



while(1)
{
	speed  = read_float32("in_data"); //current speed from vehicle model
	speed_req  = read_float32("in_data"); //required speed from vehicle model
	//speed  = read_float32("in_data");
	//speed_ref_temp  = read_float32("in_data");

	/*
	*
        * Pseudocode from Wikipedia
        * 
        previous_error = 0
        integral = 0 
      	start:
        error = setpoint - PV(actual_position)
        integral = integral + error*dt
        derivative = (error - previous_error)/dt
        output = Kp*error + Ki*integral + Kd*derivative
        previous_error = error
        wait(dt)
        goto start
        */


        // calculate the difference between the desired value and the actual value
        error = speed_req - speed;  
        // track error over time, scaled to the timer interval
	temp = fpmul32f(error, delta);
	temp2 = fpadd32(integral, temp);
//        integral = integral + (error * delta);

        // determin the amount of change from the last time checked
        temp_d = fpsub32f(error, prev_error); // / delta; 
	derivative = fdiv32f(temp_d,delta);
        // calculate how much drive the output in order to get to the 
        // desired setpoint. 

	temp_prop = fpmul32f(Kp,error);
	temp_i = fpmul32f(Ki,temp2);
	temp_d1 = fpmul32f(Kd,derivative);
       	
	temp_output = fpadd32f(temp_prop,temp_i);
	output = fpadd32f(temp_output,temp_d1);
	
        // remember the error for the next time around.
        prev_error = error;             
            

}
}	
