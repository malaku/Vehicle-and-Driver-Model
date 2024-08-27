/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: ert_main.c
 *
 * Code generated for Simulink model 'VandD'.
 *
 * Model version                  : 1.68
 * Simulink Coder version         : 9.7 (R2022a) 13-Nov-2021
 * C/C++ source code generated on : Sun Jan 28 01:50:34 2024
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: ARM Compatible->ARM Cortex
 * Emulation hardware selection:
 *    Differs from embedded hardware (MATLAB Host)
 * Code generation objectives:
 *    1. Execution efficiency
 *    2. RAM efficiency
 * Validation result: Not run
 */

#include <stdio.h>
#include <stdlib.h>
#include "VandD.h"
#include "rtwtypes.h"
#include "limits.h"
#include "linuxinitialize.h"
#include <string.h>
#include <sys/socket.h>
#include <sys/ioctl.h>
#include <net/if.h>
#include <unistd.h> // Add this line for write() and close()
#include <linux/can.h>
#include <linux/can/raw.h>
#include <arpa/inet.h>

#define UNUSED(x)                      x = x
#define NAMELEN                        16
#define SERVER_PORT 8080
int det=0;
/* Function prototype declaration*/
int openCANSocket();
void closeCANSocket(int socketDescriptor);
void sendCANMessage(int socketDescriptor, uint32_t messageID, float value) ;
void receiveCANMessage(int socketDescriptor, int targetID, float *tempvar);
void exitFcn(int sig);
void *terminateTask(void *arg);
void *baseRateTask(void *arg);
void *subrateTask(void *arg);
volatile boolean_T stopRequested = false;
volatile boolean_T runModel = true;
sem_t stopSem;
sem_t baserateTaskSem;
pthread_t schedulerThread;
pthread_t baseRateThread;
void *threadJoinStatus;
int terminatingmodel = 0;
int switchvar=1;
int openEthernetSocket() {
    int sockfd;
    struct sockaddr_in serverAddr;

    // Create socket
    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd == -1) {
        perror("socket");
        exit(EXIT_FAILURE);
    }

    // Set up server address
    memset(&serverAddr, 0, sizeof(serverAddr));
    serverAddr.sin_family = AF_INET;
    serverAddr.sin_port = htons(SERVER_PORT); // Example port number
    serverAddr.sin_addr.s_addr = inet_addr("127.0.0.1"); // Example IP address

    // Connect to the server
    if (connect(sockfd, (struct sockaddr*)&serverAddr, sizeof(serverAddr)) == -1) {
        perror("connect");
        close(sockfd);
        exit(EXIT_FAILURE);
    }

    return sockfd;
}

void closeEthernetSocket(int socketDescriptor) {
    close(socketDescriptor);
}

void sendEthernetMessage(int socketDescriptor, const char *message) {
    if (send(socketDescriptor, message, strlen(message), 0) == -1) {
        perror("send");
    }
}

void receiveEthernetMessage(int socketDescriptor, char *buffer, size_t bufferSize) {
    ssize_t bytesRead = recv(socketDescriptor, buffer, bufferSize - 1, 0);
    if (bytesRead == -1) {
        perror("recv");
    } else {
        buffer[bytesRead] = '\0'; // Null-terminate the received data
        printf("Received Ethernet message: %s\n", buffer);
    }
}

void *baseRateTask(void *arg)
{
int det=0;
int socketDescriptor = openCANSocket();
int socketDescriptorEth = openEthernetSocket();
int targetMessageID1 = 0x101; // Desired FW1
int targetMessageID2 = 0x201; // Desired FW1
int targetMessageID3 = 0x301; // Desired FW1
int targetMessageID4 = 0x401; // Desired FW1
uint32_t messageID1 = 0x011;
uint32_t messageID2 = 0x012;
uint32_t messageID3 = 0x022;
uint32_t messageID4 = 0x032;
uint32_t messageID5 = 0x042;
uint32_t messageID6 = 0x053;
uint32_t messageID7 = 0x054;
float ID1;
float ID2;
float ID3;
float ID4;
float t=0;
uint8_t receivedData[4]; // Buffer to store received data
float lamda[4];

  runModel = (rtmGetErrorStatus(rtM) == (NULL)) && !rtmGetStopRequested(rtM);
if (switchvar==1){
  while (runModel) {
    sem_wait(&baserateTaskSem);
   if (det==0){

	receiveCANMessage(socketDescriptor, 0x111, &ID1);
    	receiveCANMessage(socketDescriptor, targetMessageID1, &lamda[0]);
	receiveCANMessage(socketDescriptor, 0x222, &ID2);
    	receiveCANMessage(socketDescriptor, targetMessageID2, &lamda[1]);
	receiveCANMessage(socketDescriptor, 0x333, &ID3);
    	receiveCANMessage(socketDescriptor, targetMessageID3, &lamda[2]);
	receiveCANMessage(socketDescriptor, 0x444, &ID4);
    	receiveCANMessage(socketDescriptor, targetMessageID4, &lamda[3]);

	if (ID1==0 && ID2==0 && ID3==0 && ID4==0){

    		rtU.lamda[0] = lamda[0];
    		rtU.lamda[1] = lamda[1];
    		rtU.lamda[2] = lamda[2];
    		rtU.lamda[3] = lamda[3];
		sendCANMessage(socketDescriptor,0x555, 0);
		VandD_step();
		float tf= rtY.throttleforce[0];
		float Fy1=rtY.Fy1;
		float Fy2=rtY.Fy2;
		float Fy3=rtY.Fy3;
		float Fy4=rtY.Fy4;
		float swa=rtY.Steeringwheelangle;
		float av = rtY.ActualVelocity;
		printf("Lateral Acceleration: %f\n", rtY.LateralAcceleration);
		printf("Lateral Velocity: %f\n", rtY.LateralVelocity);
		printf("Angular Velocity: %f\n", rtY.AngularVelocity);
		printf("X position: %f\n", rtY.x_position);
		printf("Y position: %f\n", rtY.y_position);
		printf("Orientation: %f\n", rtY.orientation);
		sendCANMessage(socketDescriptor,messageID1, tf);
		usleep(100);

		sendCANMessage(socketDescriptor,messageID2, Fy1);
		usleep(100);

		sendCANMessage(socketDescriptor,messageID3, Fy2);
		usleep(100);

		sendCANMessage(socketDescriptor,messageID4, Fy3);
		usleep(100);

		sendCANMessage(socketDescriptor,messageID5, Fy4);
		usleep(100);
		sendCANMessage(socketDescriptor,messageID6, swa);
		usleep(100);

		sendCANMessage(socketDescriptor,messageID7, av);

	} else {
		sendCANMessage(socketDescriptor, 0x555, 1);
		det=1;
	}
} else {
	if (det==1) {
		det=2;

		VandD_step();
float swa=rtY.Steeringwheelangle;

		sendCANMessage(socketDescriptor,messageID6, swa);
	} else {
		receiveCANMessage(socketDescriptor, targetMessageID1, &lamda[0]);
		receiveCANMessage(socketDescriptor, targetMessageID2, &lamda[1]);
    		rtU.lamda[0] = lamda[0];
    		rtU.lamda[1] = lamda[1];
		rtU.lamda[2]=0;
		rtU.lamda[3]=0;
		VandD_step();
float swa=rtY.Steeringwheelangle;

		sendCANMessage(socketDescriptor,messageID6, swa);

	}
}
t=t+0.010000;
printf("time: %f\n",t);

    /* Get model outputs here */
    stopRequested = !((rtmGetErrorStatus(rtM) == (NULL)) && !rtmGetStopRequested
                      (rtM));
    runModel = !stopRequested;
  }
closeCANSocket(socketDescriptor);
closeEthernetSocket(socketDescriptorEth);

  runModel = 0;
  terminateTask(arg);
  pthread_exit((void *)0);
  return NULL;
}else{
while (runModel) {
    sem_wait(&baserateTaskSem);
receiveCANMessage(socketDescriptor, targetMessageID1, &lamda[0]);
receiveCANMessage(socketDescriptor, targetMessageID2, &lamda[1]);
receiveCANMessage(socketDescriptor, targetMessageID3, &lamda[2]);
receiveCANMessage(socketDescriptor, targetMessageID4, &lamda[3]);

		rtU.lamda[0] = lamda[0];
    		rtU.lamda[1] = lamda[1];
    		rtU.lamda[2] = lamda[2];
    		rtU.lamda[3] = lamda[3];
VandD_step();
		float tf= rtY.throttleforce[0];
		float Fy1=rtY.Fy1;
		float Fy2=rtY.Fy2;
		float Fy3=rtY.Fy3;
		float Fy4=rtY.Fy4;
		float swa=rtY.Steeringwheelangle;
		float av = rtY.ActualVelocity;
		printf("Lateral Acceleration: %f\n", rtY.LateralAcceleration);
		printf("Lateral Velocity: %f\n", rtY.LateralVelocity);
		printf("Angular Velocity: %f\n", rtY.AngularVelocity);
		printf("X position: %f\n", rtY.x_position);
		printf("Y position: %f\n", rtY.y_position);
		printf("Orientation: %f\n", rtY.orientation);
		sendCANMessage(socketDescriptor,messageID1, tf);
		usleep(100);

		sendCANMessage(socketDescriptor,messageID2, Fy1);
		usleep(100);

		sendCANMessage(socketDescriptor,messageID3, Fy2);
		usleep(100);

		sendCANMessage(socketDescriptor,messageID4, Fy3);
		usleep(100);

		sendCANMessage(socketDescriptor,messageID5, Fy4);
		usleep(100);
		sendCANMessage(socketDescriptor,messageID6, swa);
		usleep(100);

		sendCANMessage(socketDescriptor,messageID7, av);

    stopRequested = !((rtmGetErrorStatus(rtM) == (NULL)) && !rtmGetStopRequested
                      (rtM));
    runModel = !stopRequested;
  }
}
closeCANSocket(socketDescriptor);
closeEthernetSocket(socketDescriptorEth);

  runModel = 0;
  terminateTask(arg);
  pthread_exit((void *)0);
  return NULL;

}

void exitFcn(int sig)
{
  UNUSED(sig);
  rtmSetErrorStatus(rtM, "stopping the model");
}

void *terminateTask(void *arg)
{
  UNUSED(arg);
  terminatingmodel = 1;

  {
    runModel = 0;
  }

  sem_post(&stopSem);
  return NULL;
}

int main(int argc, char **argv)
{
  rtmSetErrorStatus(rtM, 0);

  /* Initialize model */
  VandD_initialize();

  /* Call RTOS Initialization function */
  myRTOSInit(0.01, 0);

  /* Wait for stop semaphore */
  sem_wait(&stopSem);

#if (MW_NUMBER_TIMER_DRIVEN_TASKS > 0)

  {
    int i;
    for (i=0; i < MW_NUMBER_TIMER_DRIVEN_TASKS; i++) {
      CHECK_STATUS(sem_destroy(&timerTaskSem[i]), 0, "sem_destroy");
    }
  }

#endif

  return 0;
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
