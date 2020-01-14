/***************************************************************************
 * Communication SORO-board - Version 1.0.1
 * by Brandon Caasenbrood: - b.j.caasenbrood@tue.nl
 *
 * This code is licensed under the TU/e license
 ***************************************************************************/
#include <iostream>
#include <Eigen.h>
#include <Dense>
#include "SORO.h"
#include "Control.h"
#include "src/eigstatespace.hpp"

using namespace Eigen;

VectorXd Y(STATEDIM); 
VectorXd U(STATEDIM);

/*Constructor (...)*********************************************************
 *    Setup the control process. Not required to be called.
 ***************************************************************************/
void Control::beginControl(){
  initStateSpace();
}

/*Constructor (...)*********************************************************
 *    Setup the control process. Not required to be called.
 ***************************************************************************/
void Control::controlProcess(){
  if(runActive){
    updateTime();         // update time t =+ dt
    update_statespace(t,dt,Y, U);
  }
  else{
    initStateSpace();     // reset the state space
  }
}

/*Constructor (...)*********************************************************
 *    Setup the communication process. Not required to be called.
 ***************************************************************************/
void Control::initStateSpace(){
  Y << 0,2.5;
  U << 0,1;
  t = 0;
}

/*Constructor (...)*********************************************************
 *    Setup the control process. Not required to be called.
 ***************************************************************************/
bool Control::getControlProcessStatus(){
  return runActive;
}

/*Constructor (...)*********************************************************
 *
 ***************************************************************************/
void Control::setActive(bool tmp){
  runActive = tmp;
}

/*Constructor (...)*********************************************************
 *
 ***************************************************************************/
void Control::updateTime(){
  t += dt;
}

/*Constructor (...)*********************************************************
 *
 ***************************************************************************/
void Control::setUpdateRateControl(float tmp){
  dt = tmp*1e-6;
}

/*Constructor (...)*********************************************************
 *
 ***************************************************************************/
void Control::setControlGain(float tmp, int id){
    if(id == 1){Kp = tmp; }
    if(id == 2){Kd = tmp; }
    if(id == 3){Ki = tmp; }
}

////////////////////////////////////////////////////////////////////////
void Control::showStates(){
  for(int i = 0; i < STATEDIM; i++){
    Console.print(Y(i),DECIMEL);
    Console.print("\t");
  }
}

////////////////////////////////////////////////////////////////////////
void Control::showControllerTime(){
  Console.print(t,2);
  Console.print("\t");
}

////////////////////////////////////////////////////////////////////////
void Control::parseDataControl(int msg, float data){
  if(msg == ParseMsgY0){Y(0) = data;}
  if(msg == ParseMsgY1){Y(1) = data;}
  if(msg == ParseMsgY2){Y(2) = data;}
  if(msg == ParseMsgU0){U(0) = data;}
  if(msg == ParseMsgU1){U(1) = data;}
  if(msg == ParseMsgU2){U(2) = data;}
}

////////////////////////////////////////////////////////////////////////
float Control::getDataControl(int msg){
  if(msg == ParseMsgY0){return Y[0];}
  if(msg == ParseMsgY1){return Y[1];}
  if(msg == ParseMsgY2){return Y[2];}
  if(msg == ParseMsgU0){return U[0];}
  if(msg == ParseMsgU1){return U[1];}
  if(msg == ParseMsgU2){return U[2];}
  if(msg == ParseMsgT){return t;}
}
