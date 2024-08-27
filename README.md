# Vehicle-and-Driver-Model Repository

## Overview

Welcome to the Vehicle-and-Driver-Model repository, a crucial part of a comprehensive Vehicle ECU and CAN Modeling project. This project simulates vehicle and driver dynamics, focusing on real-time communication between various ECUs through the CAN and Ethernet protocols.

This repository serves as the central hub of the project, integrating data from the following related repositories:

- **FW1**: Front left wheel ECU.
- **FW2**: Front right wheel ECU.
- **RW1**: Rear left wheel ECU.
- **RW2**: Rear right wheel ECU.
- **WAC1**: Wheel and axle controller 1.
- **WAC2**: Wheel and axle controller 2.
- **WAC3**: Wheel and axle controller 3.
- **WAC4**: Wheel and axle controller 4.

## Code Content

The Vehicle-and-Driver-Model repository contains the core simulation code, which models vehicle dynamics and driver behavior. It communicates with the wheel and axle ECUs to synchronize the entire vehicle system.

### Key Files:

- **`ert_main.c`**: Contains the main program for vehicle and driver simulation, including CAN and Ethernet communication.
- **`VandD.c`**: Core logic for the vehicle and driver model.
- **`VandD_data.c`**: Data management and configuration for the vehicle and driver model.
- **`linuxinitialize.c`**: Linux-specific initialization routines.

## Compilation

To compile the code, use the following GCC command:

```bash
gcc -o mycode /home/debian/VandD_ert_rtw/ert_main.c /home/debian/VandD_ert_rtw/VandD.c /home/debian/VandD_ert_rtw/VandD_data.c /home/debian/VandD_ert_rtw/linuxinitialize.c -I/home/debian/VandD_ert_rtw -lm -lpthread
```
This command will compile the necessary source files into an executable named mycode.

## CAN Setup
Before running the compiled code, ensure that the CAN interface is correctly set up. The setup is done using the following bash script:

```bash
bash can_setup.sh can1 1000000
```
This command configures the CAN interface (can1) with a baud rate of 1 Mbps.

## Execution
Once the CAN interface is set up and the code is compiled, you can execute the program using:

```bash
./mycode
```
This will run the Vehicle and Driver Model code, enabling it to communicate with other ECUs in the simulated vehicle network.
