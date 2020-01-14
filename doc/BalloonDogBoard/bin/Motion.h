#ifndef _MOTION_H_
#define _MOTION_H_

#include "src/interface.hpp"
//#include <iostream>
//#include <Eigen.h>
//#include <Dense>
#define FREQLP 2
#define DEADZONE 300
#define AMPLITUDE 50
#define OMEGA .25

// M1
const int EN_A = 8;
const int EN_B = 13;
const int EN_C = 2;
const int EN_D = 7;
const int IN1 = 10;
const int IN2 = 9;
const int IN3 = 11;
const int IN4 = 12;
const int IN5 = 4;
const int IN6 = 3;
const int IN7 = 5;
const int IN8 = 6;
const int LEDPIN = 13;
const int PWMMAX = 4095;
const int PIDMAX = 9999;

struct MotorState{
  float pwmA = 0;
  float pwmB = 0;
  float pwmC = 0;
  float pwmD = 0;
};

class Motion 
{
  public:
  
    bool motionActive = false;
    bool checkActive = false;
    bool motionTrigger = true;

    const int ParseMsgPwmMotA = 0;
    const int ParseMsgPwmMotB = 1;
    const int ParseMsgTime = 2;
    const int ParseMsgErrorX = 3;
    const int ParseMsgErrorY = 4;
    const int ParseMsgErrorP1 = 5;
    const int ParseMsgErrorP2 = 6;
    const int ParseMsgErrorP3 = 7;

    float InputMotA = 0;
    float InputMotB = 0;
  
    const float VMAX = 12;
    const float VMIN = 0;

    float dt;
    float t;
  
    void beginMotion();
    void motionProcess();
    void setMotionActive(bool tmp);
    void setCheckActive(bool tmp);
    void setMotionTrigger(bool tmp);
    void setMotor(int id, float state);
    void updateTimeMotion(float tmp);
    void setUpdateRateMotion(float tmp);
    void setSetPoint(float X, float Y);
    void checkHardware();
    void showMotion();

    double FeedForwardMotor(double PWM);
    void OptimalTorqueSolutions(double X, double Y, double* v);

    float getTimeMotion();

    void parseDataMotion(int msg, float data);
    float getDataMotion(int msg);
};

#endif
