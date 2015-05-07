#include <pthread.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthreadUtils.h>
#include <Pipes.h>
#include <pipeHandler.h>
#ifndef SW
#include "vhdlCStubs.h"
#endif
void topModule();
void thread1();
void thread2();
#ifdef SW
DEFINE_THREAD(topModule);
DEFINE_THREAD(thread2)
DEFINE_THREAD(thread1)
#endif

int main(int argc, char* argv[])
{
#ifdef SW
	init_pipe_handler();
	PTHREAD_DECL(topModule);   // declare the Daemon thread.
	PTHREAD_CREATE(topModule); // start the Daemon..
    register_pipe("in_data",1,32,0);
	register_pipe("in_data1",1,32,0);
	register_pipe("in_data2",1,32,0);
	register_pipe("out_data1",1,32,0);
	register_pipe("synch_pipe12",1,32,0);
	register_pipe("synch_pipe21",1,32,0);

//	register_pipe("out_data2");
/*	PTHREAD_DECL(thread2);   // declare the Daemon thread.
	PTHREAD_CREATE(thread2); // start the Daemon..
	PTHREAD_DECL(thread1);   // declare the Daemon thread.
	PTHREAD_CREATE(thread1); // start the Daemon..*/
#endif
//	while(1)
	//	{
		uint32_t a, b;
		//char yn;
/*		fprintf(stdout,"Enter 2 input values\n");
		scanf("%d", &a);
		scanf("%d", &b);

		if(a == b)
		{
         fprintf(stderr,"Lost your chance!! Giving up.\n");
 //        break;   
		}
		else{
		fprintf(stderr, "Great! Good to go..\n \n");
		}
/*
#ifdef SW
    fprintf(stderr,"Creating software threads. \n");
    PTHREAD_DECL(thread2);   // declare the Daemon thread.
	PTHREAD_CREATE(thread2); // start the Daemon..
	PTHREAD_DECL(thread1);   // declare the Daemon thread.
	PTHREAD_CREATE(thread1); // start the Daemon..
#endif */

while(1)
{
#ifdef SW
    fprintf(stderr,"Creating software threads. \n");
    PTHREAD_DECL(thread2);   // declare the Daemon thread.
	PTHREAD_CREATE(thread2); // start the Daemon..
	PTHREAD_DECL(thread1);   // declare the Daemon thread.
	PTHREAD_CREATE(thread1); // start the Daemon..
#endif 

//		write_uint32("in_data",a);
//		write_uint32("in_data",b);
//		fprintf(stderr,"Great! Data sent! \n ");

//while(1)
//{

        
	    uint32_t c = read_uint32("out_data1");
       	uint32_t d = read_uint32("out_data1");

		fprintf(stderr,"Result\n c = %d.\n d = %d\n", c, d);


#ifdef SW
        fprintf(stderr,"Do you want to continue? Type y/n?\n");
        char yn = scanf("%c", &yn);

        if(yn == 'n')
        {
        break;
        }
#endif
		//if(a == 0)
		//	break;

//}	
#ifdef SW
    PTHREAD_CANCEL(topModule);
	PTHREAD_CANCEL(thread2);
	PTHREAD_CANCEL(thread1);
	close_pipe_handler();
#endif	

	}


	return(0);
}

