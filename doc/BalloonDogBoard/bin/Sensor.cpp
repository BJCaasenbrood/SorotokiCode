#include "SORO.h"
#include "Sensor.h"
#include "src/MahonyAHRS.h"
#include "src/MPU6050_tockn.h"
#include "src/i2cmultiplex.hpp"
#include "src/Adafruit_MPRLS.h"
#include "src/FilterOnePole.h"

extern "C" {
#include "utility/twi.h"
}

#include <math.h>
#include <Wire.h>

#define PI 3.141592653589793
#define RAD2DEG 57.29577951308200
#define DEG2RAD 0.017453292519943
#define DEC 3
#define RESET_PIN -1
#define EOC_PIN -1
#define TCAADDR 0x70

Adafruit_MPRLS mpr = Adafruit_MPRLS(RESET_PIN, EOC_PIN);
MPU6050 mpu(Wire);
Mahony filter;

FilterOnePole lp_p1(LOWPASS, FREQLP);
FilterOnePole lp_p2(LOWPASS, FREQLP);
FilterOnePole lp_p3(LOWPASS, FREQLP);

SensorState sen;

/**************************************************************************/
// Prompt symbol
/**************************************************************************/
void Sensor::beginSensor() {
  initSensors();
}

//***************************************************************************
void Sensor::initSensors() {
  Wire.begin();
  if (ImuConnection) {
    startIMUConnection();
  }
  else {
    ImuConnection = false;
  }

  if (MprConnection) {
    startMPRConnection(0);
    startMPRConnection(1);
    startMPRConnection(2);
  }
  else {
    MprConnection = false;
  }
}

//---------------------------------------------------------------------
void Sensor::setSensorActive(bool tmp) {
  sensorActive = tmp;
}

//---------------------------------------------------------------------
void Sensor::setImu(bool tmp) {
  if (tmp) {
    Console.println("* searching IC2 connection BN055...");
  }
  else {
    Console.println("* BN055 connection terminated! ");
  };
  ImuConnection = tmp;
  startIMUConnection();
}

//---------------------------------------------------------------------
void Sensor::startIMUConnection() {
  if (!mpu.begin())
  {
    Console.println("* IMU error! no MPU650 detected ... ");
    Console.println("* check your wiring or I2C ADDRESS!");
    ImuConnection = false;
  }
  else {
    Console.println("* MPU650 IMU connected!");
    mpu.calcGyroOffsets(true);
    breakLine();
    ImuConnection = true;
  }
}

//---------------------------------------------------------------------
void Sensor::startMPRConnection(int tcaID) {
  //Wire.begin();
  tcaSelect(tcaID);
  if (!mpr.begin())
  {
    Console.println("* MPR error! no MPRLS sensor detected ... ");
    Console.println("* check your wiring or I2C ADDRESS!");
    MprConnection = false;
  }
  else {
    Console.print("* @tca~");
    Console.print(tcaID);
    Console.print(" | ");
    Console.println("MPRLS sensor connected! ");
    //mpr.calcPressureOffset(true);
    MprConnection = true;
  }
}

//---------------------------------------------------------------------
void Sensor::setHX711(bool tmp) {

}

//---------------------------------------------------------------------
bool Sensor::getSensorProcessStatus() {
  return sensorActive;
}

//---------------------------------------------------------------------
void Sensor::sensorProcess() {
  
  if (ImuConnection) {
    mpu.update();
    //    float ax = mpu.getAccX();
    //    float ay = mpu.getAccY();
    //    float az = mpu.getAccZ();
    //    float gx = mpu.getGyroX();
    //    float gy = mpu.getGyroY();
    //    float gz = mpu.getGyroZ();
    //    filter.updateIMU(gx, gy, gz, ax, ay, az);
  }

  readSensor(false);

  if (tcaScanActive) {
    scanTCAAdress();
  }

}

//---------------------------------------------------------------------
void Sensor::readSensor(bool enablePrint) {

  if (ImuConnection) {

    sen.IMUx = mpu.getAngleX();
    sen.IMUy = mpu.getAngleY();
    sen.IMUax = mpu.getAccX();
    sen.IMUay = mpu.getAccY();
    sen.IMUaz = mpu.getAccZ();
    sen.IMUax = mpu.getGyroX();
    sen.IMUay = mpu.getGyroX();
    sen.IMUaz = mpu.getGyroX();

    if (enablePrint) {
      Console.print(sen.IMUx, DEC); Console.print("\t");
      Console.print(sen.IMUy, DEC); Console.print("\t");
    }
  }
  if (MprConnection) {
    readPressure();
    if (enablePrint) {
      Console.print(sen.P1, DEC); Console.print("\t");
      Console.print(sen.P2, DEC); Console.print("\t");
      Console.print(sen.P3, DEC); Console.print("\t");
    }
  }
}

//---------------------------------------------------------------------
void Sensor::readPressure() {

  for (int i = 0; i <= 2; i++) {
    tcaSelect(i);
    mpr.prepareRead();
  }
  for (int i = 0; i <= 2; i++) {
    tcaSelect(i);
    if (i == 0) {
      sen.P1 = lp_p1.input(mpr.readPressure()) - 100;
    }
    if (i == 1) {
      sen.P2 = lp_p2.input(mpr.readPressure()) - 100;
    }
    if (i == 2) {
      sen.P3 = lp_p3.input(mpr.readPressure()) - 100;
    }
  }
}

//---------------------------------------------------------------------
void Sensor::tarePressure() {
}

//---------------------------------------------------------------------
void Sensor::readIMU() {

}

//---------------------------------------------------------------------
void Sensor::showSensor() {
  Console.print(sen.IMUx, DEC);
  Console.print("\t");
  Console.print(sen.IMUy, DEC);
  Console.print("\t");
  Console.print(sen.P1, DEC);
  Console.print("\t");
  Console.print(sen.P2, DEC);
  Console.print("\t");
  Console.print(sen.P3, DEC);
  Console.print("\t");
}

//---------------------------------------------------------------------
void Sensor::setTCAScanActive(bool tmp) {
  tcaScanActive = tmp;
}

//---------------------------------------------------------------------
void Sensor::scanTCAAdress() {
  for (uint8_t t = 0; t < 8; t++) {
    tcaSelect(t);
    Console.print("TCA Port #");
    Console.print(t); Console.println(": ");

    for (uint8_t addr = 0; addr <= 127; addr++) {
      if (addr == TCAADDR) continue;

      uint8_t data;
      if (! twi_writeTo(addr, &data, 0, 1, 1)) {
        Console.print("* Found i2c device at address: 0x");
        Console.println(addr, HEX);
      }
    }
    Console.println("");
  }
  setTCAScanActive(false);
}

//---------------------------------------------------------------------
float Sensor::angleWrap(float x) {
  x = fmod(x, 360);
  if (x < 0) {
    x += 360;
  }
  return x;
}

//---------------------------------------------------------------------
void Sensor::parseDataSensor(int msg, float data) {
  if (msg == ParseImuX) {}
  if (msg == ParseImuY) {}
  if (msg == ParseImuZ) {}
}

//---------------------------------------------------------------------
float Sensor::getDataSensor(int msg) {
  if (msg == ParseImuX) {
    return sen.IMUx;
  }
  if (msg == ParseImuY) {
    return sen.IMUy;
  }
  if (msg == ParseImuZ) {
    return sen.IMUz;
  }
  if (msg == ParseP1) {
    return sen.P1;
  }
  if (msg == ParseP2) {
    return sen.P2;
  }
  if (msg == ParseP3) {
    return sen.P3;
  }
}
