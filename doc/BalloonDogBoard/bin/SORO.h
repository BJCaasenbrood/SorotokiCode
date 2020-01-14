/*  SOROTOKI: Toolkit for the soft robotics
    Copyright (c) 2018 Brandon Caasenbrood  <b.j.caasenbrood@tue.nl>
*/
#ifndef _SORO_H_
#define _SORO_H_

#include "Arduino.h"

#undef max
#undef min

#include "Motion.h"
#include "Control.h"
#include "Communication.h"
#include "Sensor.h"
#include "Timer.h"
#include "Wire.h"
#include "src/interface.hpp"
#include "src/rang.hpp"

using namespace std;
using namespace rang;

#define SORO_TIMESCALE 1e-6
#define SORO_LINECHAR '-'
#define SORO_STOPCHAR '~'
#define SORO_CMDLENGHT 62
#define RUNTIME 300
#define AUTOLOG true

#define UPDATE_RATE_HIGHLEVEL 1000  
#define UPDATE_RATE_LOWLEVEL 10000

enum ARGUMENTS
{
  MODE,
  WHO,
  WHAT,
  VALUE,
};

enum _MODES 
{ MODE_HELP, 
  MODE_INIT, 
  MODE_IDLE,
  MODE_RUN, 
  MODE_STOP, 
  MODE_TEST,
  MODE_GAINS,
  MODE_REFERENCE,
  MODE_IMU, 
  MODE_READ_SEN,
  MODE_MOT, 
  MODE_SOL,
  MODE_I2CSCAN, 
  MODE_SETPOINT,
  MODE_TCASCAN, 
  MODE_A,
  MODE_S, 
  MODE_D, 
  MODE_W,
  MODE_GRIPPER 
};

enum _WHAT{
  ON, // on
  OFF, // hold/off
  RESET, // reset
  GET, // get
  SET // set
};

enum _WHO{
  SYSTEM,
  IMU, 
  MOTOR_ID1,
  MOTOR_ID2,
  MOTOR_ID3,
  MOTOR_ID4,
  PRESS_ID1,
  PRESS_ID2,
  PRESS_ID3,
  PRESS_ID4,
};

//const int MODE_HELP = -2;
//const int MODE_INIT = -1;
//const int MODE_IDLE = 0;
//const int MODE_RUN  = 1;
//const int MODE_STOP = -3;
//const int MODE_TEST  = 2;
//const int MODE_SET_GAINS = 10;
//const int MODE_SET_REFERENCE = 11;
//const int MODE_SET_IMU = 13;
//const int MODE_SET_HX711 = 14;
//const int MODE_READ_SEN = 15;
//const int MODE_SET_MOT = 16;
//const int MODE_SET_SOL = 17;
//const int MODE_I2CSCAN = 18;
//const int MODE_SETPOINT = 19;
//const int MODE_TCASCAN = 20;
//const int MODE_A = 21;
//const int MODE_S = 22;
//const int MODE_D = 23;
//const int MODE_W = 24;
//const int MODE_GRIPPER = 25;

class SORO: 
public Sensor,
public Motion,
public Control,
public Communication{
public:

    int  mode;
    bool isPaused = false;
    bool runActive = false;
    bool printActive = false;
    bool GripperState = false;
    float t = 0.0;
    float dt = 0.0;
    
	  SORO();
    
    void begin();
    void run();

    void processMonitor();
    void setPause(bool tmp);
    void setRun(bool tmp);
    
    void showCycleTime();
    void setClockRate();
    void setLoopTime(int tmp);
    void updateTime();
 
    void parseDataProcess();

    void showData();
    void displayInfo();
   
};

#endif
