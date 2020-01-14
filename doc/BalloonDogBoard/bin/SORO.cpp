/***************************************************************************
   SOROTOKI: Toolkit for the soft robotics - version 1.0.1
   by Brandon Caasenbrood: - b.j.caasenbrood@tue.nl

   This c++ file is licensed under the TU/e license
 ***************************************************************************/
#include <iostream>
#include <fstream>
#include "SORO.h"

ofstream autolog;
Timer timeHighLevel;
Timer timeLowLevel;

double POSX = 0;
double POSY = 0;

SORO::SORO() {}

//***************************************************************************
void SORO::begin() {
  Console.begin();
  DisplayLogo();
  setClockRate();
  beginMotion();
  beginControl();
  beginCommunication();
  beginSensor();
  newLine();
  Prompt();
}

//***************************************************************************
void SORO::run() {
  if (timeHighLevel >= UPDATE_RATE_HIGHLEVEL) {
    processMonitor();
    communicationProcess();
    controlProcess();
    parseDataProcess();
    updateTime();
    timeHighLevel = timeHighLevel - UPDATE_RATE_HIGHLEVEL;
  }

  if (timeLowLevel >= UPDATE_RATE_LOWLEVEL) {
    sensorProcess();
    motionProcess();
    showData();
    timeLowLevel = timeLowLevel - UPDATE_RATE_LOWLEVEL;
  }

  if (getDataControl(-1) > RUNTIME) {
    setProcessMode(MODE_STOP);
    setProcessStatus(true);
  }

}

//***************************************************************************
void SORO::processMonitor() {

  awaitDataFromPC();

  if (getProcessStatus() && isPaused != true && getProcessMode() == MODE_GAINS) {
    //setGain(getProcessData(0), getProcessData(1), getProcessData(2));
    //if (getControlProcessStatus()) {
    //  setProcessMode(MODE_IDLE);
    //}
  }

  if (getProcessStatus() && isPaused != true && getProcessMode() == MODE_SETPOINT) {
    setSetPoint(getProcessData(0), getProcessData(1));
  }

  if (getProcessStatus() && isPaused != true && getProcessMode() == MODE_IMU) {
    setImu(getProcessData(0));
  }

  if (getProcessStatus() && isPaused != true && getProcessMode() == MODE_SOL) {
    //setSolenoidPWM(getProcessData(0),getProcessData(1));
  }

  if (getProcessStatus() && isPaused != true && getProcessMode() == MODE_MOT) {
    Console.println(getProcessData(1));
    setMotor(getProcessData(0), getProcessData(1));
  }

  if (getProcessStatus() && isPaused != true && getProcessMode() == MODE_A) {
    POSX += 500;
    Console.print(" POSX: "); Console.print(POSX);
    Console.print("\t | POSY: "); Console.println(POSY);
    if (POSX > 0) {
      setMotor(1, POSX);
      setMotor(3, 0);
    }
    else if (POSX < 0) {
      setMotor(1, 0);
      setMotor(3, -POSX);
    }
  }

  if (getProcessStatus() && isPaused != true && getProcessMode() == MODE_D) {
    POSX -= 500;
    Console.print(" POSX: "); Console.print(POSX);
    Console.print("\t | POSY: "); Console.println(POSY);
    if (POSX > 0) {
      setMotor(1, POSX);
      setMotor(3, 0);
    }
    else if (POSX < 0) {
      setMotor(1, 0);
      setMotor(3, -POSX);
    }
  }

  if (getProcessStatus() && isPaused != true && getProcessMode() == MODE_W) {
    POSY += 500;
    Console.print(" POSX: "); Console.print(POSX);
    Console.print("\t | POSY: "); Console.println(POSY);
    setMotor(2, POSY);
  }

  if (getProcessStatus() && isPaused != true && getProcessMode() == MODE_S) {
    POSX = 0; POSY = 0;
    Console.print(" POSX: "); Console.print(POSX);
    Console.print("\t | POSY: "); Console.println(POSY);
    setMotor(1, 0);
    setMotor(2, 0);
    setMotor(3, 0);
    setMotor(4, 0);
  }

  if (getProcessStatus() && isPaused != true && getProcessMode() == MODE_GRIPPER) {
    if (GripperState)
    {
      Console.println(" opening soft gripper ");
      setMotor(4, 0);
      GripperState = false;
    }
    else {
      Console.println(" closing soft gripper ");
      setMotor(4, 3000);
      GripperState = true;
    }
  }

  if (getProcessStatus() && isPaused != true && getProcessMode() == MODE_TEST) {
    double p[3];

    OptimalTorqueSolutions(getProcessData(0), getProcessData(1), p);

    Console.print(p[0]); Console.print(" ");
    Console.print(p[1]); Console.print(" ");
    Console.print(p[2]); Console.println(" ");

  }

  if (getProcessStatus() && isPaused != true && getProcessMode() == MODE_I2CSCAN) {
    setI2CScanActive(true);
  }

  if (getProcessStatus() && isPaused != true && getProcessMode() == MODE_TCASCAN) {
    setTCAScanActive(true);
  }

  if (getProcessStatus() && isPaused != true && getProcessMode() == MODE_READ_SEN) {
    readSensor(true);
    Console.println("");
  }

  if (getProcessStatus() && isPaused != true && getProcessMode() == MODE_RUN) {
    setPause(false);
    setRun(true);
    setActive(true);
    setCommunicationActive(true);
    setMotionActive(true);
    setProcessStatus(true);
  }

  if (getProcessStatus() && isPaused != true && getProcessMode() == MODE_STOP) {
    setRun(false);
    //setPause(true);
    setActive(false);
    setCommunicationActive(false);
    //setSensorActive(false);
    setMotionActive(false);
    setMotionTrigger(false);
    setProcessStatus(true);
  }

  if (getProcessStatus() && isPaused != true && getProcessMode() == MODE_IDLE) {
    setRun(false);
    //setPause(true);
    setActive(false);
    setCommunicationActive(false);
    //setSensorActive(false);
    setMotionActive(false);
    setMotionTrigger(false);
    setProcessMode(MODE_HELP);
    setProcessStatus(true);
  }

  if (getProcessStatus() && getProcessMode() == MODE_HELP) {
    setRun(false);
    setActive(false);
    displayInfo();
  }

  processCompleted();
}

/*Constructor (...)*********************************************************
 ***************************************************************************/
void SORO::setPause(bool tmp) {
  isPaused = tmp;
  if (isPaused) {Console.println("");}
}

/*Constructor (...)*********************************************************
 ***************************************************************************/
void SORO::setRun(bool tmp) {
  runActive = tmp;
}

/*Constructor (...)*********************************************************
 ***************************************************************************/
void SORO::updateTime() {
  if (runActive){t = getDataControl(-1);}
  else {t = 0;}
}

/*Constructor (...)*********************************************************

 ***************************************************************************/
void SORO::showCycleTime() {
  Console.print(t, 3);
  Console.print("\t");
}

/*Constructor (...)*********************************************************
 ***************************************************************************/
void SORO::setClockRate() {
  setUpdateRateControl(UPDATE_RATE_HIGHLEVEL);
  setLoopTime(UPDATE_RATE_LOWLEVEL);
  setUpdateRateMotion(UPDATE_RATE_LOWLEVEL);
  breakLine();
  Console.print("* Sample freq. = ");
  Console.print(1 / (UPDATE_RATE_HIGHLEVEL * 1e-6), 1);
  Console.print("Hz | ");
  Console.print(UPDATE_RATE_HIGHLEVEL, 1);
  Console.println("ms");

  Console.print("* Output freq. = ");
  Console.print(1 / (UPDATE_RATE_LOWLEVEL * 1e-6), 1);
  Console.print("Hz | ");
  Console.print(UPDATE_RATE_LOWLEVEL, 1);
  Console.println("ms");
}

/*Constructor (...)*********************************************************
 ***************************************************************************/
void SORO::setLoopTime(int tmp) {
  dt = tmp * SORO_TIMESCALE;
}

/*Constructor (...)*********************************************************
 ***************************************************************************/
void SORO::displayInfo() {
  breakLine();
  titleLine("Help menu (list of built-in commands)");
  breakLine();
  printout("help");
  printout("run");
  printout("stop");
  printout("i2c");
  printout("tca");
  printout("imu (bool)");
  printout("mot (id) (pwm)");
  printout("test");
  breakLine();
}

/*Constructor (...)*********************************************************
 ***************************************************************************/
void SORO::showData() {

  if (runActive) {
    //showCycleTime();
    //showControllerTime();
    //showStates();
    //showSensor();
    //showMotion();
    //Console.println("");
  }

  if (AUTOLOG) {
    if (!printActive && runActive) {
      //autolog.precision(4);
      autolog.open("soro.log");
      printActive = true;
    }
    if (printActive && runActive) {
      autolog << getDataControl(-1) << ", " << getDataSensor(0) <<  ", " << getDataSensor(1) << ", " << getDataSensor(5) << endl;
      //autolog << getDataControl(-1) << ", " << getDataSensor(3) << ", " << getDataMotion(0) << endl;
    }
    if (printActive && !runActive) {
      printActive = false;
      autolog.close();
    }
  }
}

/*Constructor (...)*********************************************************

 ***************************************************************************/
void SORO::parseDataProcess() {
  parseDataMotion(2, getDataControl(-1)); // parse time
  parseDataMotion(3, getDataSensor(0)); // parse roll angle
  parseDataMotion(4, getDataSensor(1)); // parse pitch angle
  parseDataMotion(5, getDataSensor(3)); // parse pressure diff 1
  parseDataMotion(6, getDataSensor(4)); // parse pressure diff 2
  parseDataMotion(7, getDataSensor(5)); // parse pressure diff 3
}
