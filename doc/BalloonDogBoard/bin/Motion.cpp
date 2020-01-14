#include <iostream>
#include <Eigen.h>
#include <Dense>
#include "SORO.h"
#include "Motion.h"
#include "Wire.h"
#include "src/PID.h"
#include "src/PCA9685.h"
#include "src/FilterOnePole.h"
#include "src/eigquadprog.hpp"
#include "src/stepsignal.hpp"
using namespace Eigen;

#define PI 3.1415926535897

float ErrorAngleX = 0;
float ErrorAngleY = 0;

double Xd, Xin, Xout;
double Yd, Yin, Yout;

double Pd1, Pin1, Pout1;
double Pd2, Pin2, Pout2;
double Pd3, Pin3, Pout3;

double Kp = 1, Ki = 0, Kd = 0;
double Kpm = 1000, Kim = 100, Kdm = 5;

//PID PIDX(&Xin, &Xout, &Xd, Kp, Ki, Kd, DIRECT);
//PID PIDY(&Yin, &Yout, &Yd, Kp, Ki, Kd, DIRECT);
//
PID PIDM1(&Pin1, &Pout1, &Pd1, Kpm, Kim, Kdm, DIRECT);

FilterOnePole lp1(LOWPASS, FREQLP);

MotorState Motor;
PCA9685 PwmIC;

/*Constructor (...)*********************************************************
 ***************************************************************************/
void Motion::beginMotion() {
  Wire.begin();
  PwmIC.resetDevices();
  PwmIC.init();
  PwmIC.setPWMFrequency(1500);

  Xd = 0;
  Yd = 0;

  Pout1 = 0;
  Pout2 = 0;
  Pout3 = 0;
    

//  ePD1 << 0, 0, 0, 0;
//  uPD1 << 0, 0, 0, 0;
//  ePD2 << 0, 0, 0, 0, 0, 0;
//  uPD2 << 0, 0, 0, 0, 0, 0;

 
 //    PIDX.SetMode(AUTOMATIC);
  //    PIDY.SetMode(AUTOMATIC);
//      PIDM1.SetMode(AUTOMATIC);
//      PIDM2.SetMode(AUTOMATIC);
//      PIDM3.SetMode(AUTOMATIC);
  //    PIDX.SetOutputLimits(-PIDMAX, PIDMAX);
  //    PIDY.SetOutputLimits(-PIDMAX, PIDMAX);
//      PIDM1.SetOutputLimits(-PIDMAX, PIDMAX);
//      PIDM2.SetOutputLimits(-PIDMAX, PIDMAX);
//      PIDM3.SetOutputLimits(-PIDMAX, PIDMAX);
}

/*Constructor (...)*********************************************************
 ***************************************************************************/
void Motion::motionProcess() {

//  double e1, e2, u1, u2;
//  double E1, E2, E3, U1, U2, U3;
//
//  if (motionActive) {
//
//    float gainPWM = 1;
//    float a = (1 - exp(-0.2 * t));
//    Xd = a * AMPLITUDE * cos(OMEGA * t);
//    Yd = 0 * AMPLITUDE * sin(OMEGA * t);
//
//    e1 = Xd - Xin;  
//    e2 = Yd - Yin;  
//
//    VectorXd tau(2); 
//
//    tau(0) = Kp*e1;
//    tau(1) = Kp*e2;
//
//    // compute optimal pressure values
//    double U[3] = {0, 0, 0};
//    OptimalTorqueSolutions(tau(0), tau(1)*0, U);
//
//    //Pd1 = U[0];
//    //Pd2 = U[1];
//    //Pd3 = U[2];
//
//    Pd1 = a*((20 + 10*sin(t)));
//
//    PIDM1.Compute();
//
//    setMotor(1,Pout1);
    //setMotor(2,Pout2);
    //setMotor(3,Pout3);

  if (motionActive) {
    motionTrigger = false;
    setMotor(2, stepsignal(t));
    setMotor(3, stepsignal(t-10));
    setMotor(1, stepsignal(t-20));
  }

  if (motionTrigger == false && motionActive == false) {

    Pout1 = 0;
    Pout2 = 0;
    Pout3 = 0;
    
    motionTrigger = true;
    setMotor(1, 0);
    setMotor(2, 0);
    setMotor(3, 0);
    setMotor(4, 0);
  }
}

/*Constructor (...)*********************************************************
 ***************************************************************************/
double Motion::FeedForwardMotor(double PWM) {

  float u = 0;

  if (PWM <= 0) {
    return u;
  }
  else {
    u = DEADZONE + PWM;
    return u;
  }
}

////////////////////////////////////////////////////////////////////////
void Motion::OptimalTorqueSolutions(double X, double Y, double *v) {

  MatrixXd G(3, 3);
  VectorXd g0(3);
  MatrixXd CE(3, 2);
  VectorXd ce0(2);
  MatrixXd CI(3, 3);
  VectorXd ci0(3);
  VectorXd x(3);

  G = MatrixXd::Identity(3, 3);

  g0  << 0, 0, 0;
  CE  << -1, 0, 0.5,  sqrt(3) / 2, 0.5, -sqrt(3) / 2;
  ce0 << -X, -Y;
  CI  << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  ci0 << 0, 0, 0;

  solve_quadprog(G, g0,  CE, ce0,  CI, ci0, x);

  v[0] = x(0);
  v[1] = x(1);
  v[2] = x(2);
}

////////////////////////////////////////////////////////////////////////
void Motion::setMotionActive(bool tmp) {
  motionActive = tmp;
}

////////////////////////////////////////////////////////////////////////
void Motion::setCheckActive(bool tmp) {
  checkActive = tmp;
}

////////////////////////////////////////////////////////////////////////
void Motion::setMotionTrigger(bool tmp) {
  motionTrigger = tmp;
}

////////////////////////////////////////////////////////////////////////
void Motion::setMotor(int id, float state) {

  if (id == 1 && state != Motor.pwmA) {
    if (state > 0 && state < PWMMAX) {
      PwmIC.setChannelPWM(IN1, 4095);
      PwmIC.setChannelPWM(IN2, 0);
      PwmIC.setChannelPWM(EN_A, state);
      Motor.pwmA = state;
    }
    else if (state >= PWMMAX) {
      PwmIC.setChannelPWM(IN1, 4095);
      PwmIC.setChannelPWM(IN2, 0);
      PwmIC.setChannelPWM(EN_A, state);
      Motor.pwmA = PWMMAX;
    }
    else {
      PwmIC.setChannelPWM(IN1, 0);
      PwmIC.setChannelPWM(IN2, 0);
      PwmIC.setChannelPWM(EN_A, 0);
      Motor.pwmA = 0;
    }
  }

  if (id == 2 && state != Motor.pwmB) {
    if (state > 0 && state < PWMMAX) {
      PwmIC.setChannelPWM(IN3, 4095);
      PwmIC.setChannelPWM(IN4, 0);
      PwmIC.setChannelPWM(EN_B, state);
      Motor.pwmB = state;
    }
    else if (state >= PWMMAX) {
      PwmIC.setChannelPWM(IN3, 4095);
      PwmIC.setChannelPWM(IN4, 0);
      PwmIC.setChannelPWM(EN_B, state);
      Motor.pwmB = PWMMAX;
    }
    else {
      PwmIC.setChannelPWM(IN3, 0);
      PwmIC.setChannelPWM(IN4, 0);
      PwmIC.setChannelPWM(EN_B, 0);
      Motor.pwmB = 0;
    }
  }

  if (id == 3 && state != Motor.pwmC) {
    if (state > 0 && state < PWMMAX) {
      PwmIC.setChannelPWM(IN5, 4095);
      PwmIC.setChannelPWM(IN6, 0);
      PwmIC.setChannelPWM(EN_C, state);
      Motor.pwmC = state;
    }
    else if (state >= PWMMAX) {
      PwmIC.setChannelPWM(IN5, 4095);
      PwmIC.setChannelPWM(IN6, 0);
      PwmIC.setChannelPWM(EN_C, state);
      Motor.pwmC = PWMMAX;
    }
    else {
      PwmIC.setChannelPWM(IN5, 0);
      PwmIC.setChannelPWM(IN6, 0);
      PwmIC.setChannelPWM(EN_C, 0);
      Motor.pwmC = 0;
    }
  }

  if (id == 4 && state != Motor.pwmD) {
    if (state > 0 && state < PWMMAX) {
      PwmIC.setChannelPWM(IN7, 4095);
      PwmIC.setChannelPWM(IN8, 0);
      PwmIC.setChannelPWM(EN_D, state);
      Motor.pwmD = state;
    }
    else if (state >= PWMMAX) {
      PwmIC.setChannelPWM(IN7, 4095);
      PwmIC.setChannelPWM(IN8, 0);
      PwmIC.setChannelPWM(EN_D, state);
      Motor.pwmD = PWMMAX;
    }
    else {
      PwmIC.setChannelPWM(IN7, 0);
      PwmIC.setChannelPWM(IN8, 0);
      PwmIC.setChannelPWM(EN_D, 0);
      Motor.pwmD = 0;
    }
  }
}

////////////////////////////////////////////////////////////////////////
void Motion::updateTimeMotion(float tmp) {
  t = tmp;
}

////////////////////////////////////////////////////////////////////////
float Motion::getTimeMotion() {
  return t;
}

////////////////////////////////////////////////////////////////////////
void Motion::setUpdateRateMotion(float tmp) {
  dt = tmp*1e-6;
}

////////////////////////////////////////////////////////////////////////
void Motion::setSetPoint(float X, float Y) {
  //dt = 1/tmp;
  Xd = X;
  Yd = Y;
  Console.print("* Setpoint IMU | x: ");
  Console.print(Xd);
  Console.print(" | y: ");
  Console.print(Yd);
  Console.println(" | (deg)");
}

////////////////////////////////////////////////////////////////////////
void Motion::checkHardware() {
  setCheckActive(false);
}

////////////////////////////////////////////////////////////////////////
void Motion::showMotion() {

  if (motionActive) {
    Console.print(Pd1);
    Console.print("\t");
    Console.print(Pin1);
    Console.print("\t");
    Console.print(Pout1);
    Console.print("\t");
//    Console.print(Pd2);
//    Console.print("\t");
//    Console.print(Pin2);
//    Console.print("\t");
//    Console.print(Pout2);
//    Console.print("\t");
//    Console.print(Pd3);
//    Console.print("\t");
//    Console.print(Pin3);
//    Console.print("\t");
//    Console.print(Pout3);
//    Console.print("\t");
  }
}

////////////////////////////////////////////////////////////////////////
void Motion::parseDataMotion(int msg, float data) {
  if (msg == ParseMsgPwmMotA) {
    Motor.pwmA = data;
  }
  if (msg == ParseMsgPwmMotB) {
    Motor.pwmB = data;
  }
  if (msg == ParseMsgErrorX) {
    Xin = data;
  }
  if (msg == ParseMsgErrorY) {
    Yin = data;
  }
  if (msg == ParseMsgErrorP1) {
    Pin1 = data;
  }
  if (msg == ParseMsgErrorP2) {
    Pin2 = data;
  }
  if (msg == ParseMsgErrorP3) {
    Pin3 = data;
  }
  if (msg == ParseMsgTime) {
    t = data;
  }
}

////////////////////////////////////////////////////////////////////////
float Motion::getDataMotion(int msg) {
  if (msg == ParseMsgPwmMotA) {
    return Motor.pwmA;
  }
  if (msg == ParseMsgPwmMotB) {
    return Motor.pwmB;
  }
}
