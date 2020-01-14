#include "MPU6050_tockn.h"
#include "Arduino.h"

MPU6050::MPU6050(TwoWire &w){
  wire = &w;
  accCoef = 0.01f; // 0.02f
  gyroCoef = 0.99f; // 0.98f
}

////////////////////////////////////////////////////////////////////////
MPU6050::MPU6050(TwoWire &w, float aC, float gC){
  wire = &w;
  accCoef = aC;
  gyroCoef = gC;
}

////////////////////////////////////////////////////////////////////////
bool MPU6050::begin(){
  writeMPU6050(MPU6050_SMPLRT_DIV, 0x00);
  writeMPU6050(MPU6050_CONFIG, 0x00);
  writeMPU6050(MPU6050_GYRO_CONFIG, 0x08);
  writeMPU6050(MPU6050_ACCEL_CONFIG, 0x00);
  writeMPU6050(MPU6050_PWR_MGMT_1, 0x01);
  this->update();
  angleGyroX = 0;
  angleGyroY = 0;
  angleX = this->getAccAngleX();
  angleY = this->getAccAngleY();
  preInterval = millis();

  if (!isConnected())
  {
        return false;
  }
}
////////////////////////////////////////////////////////////////////////
void MPU6050::setGyroOffsets(float x, float y, float z){
  gyroXoffset = x;
  gyroYoffset = y;
  gyroZoffset = z;
}

////////////////////////////////////////////////////////////////////////
void MPU6050::calcGyroOffsets(bool console, uint16_t delayBefore, uint16_t delayAfter){
	float x = 0, y = 0, z = 0;
	int16_t rx, ry, rz;

  //delay(delayBefore);
	if(console){
    Console.println();
    Console.print("* IMU: calibration in progress");
  }
  for(int i = 0; i < 3000; i++){
    if(console && i % 600 == 0){
      Console.print(".");
    }
    wire->beginTransmission(MPU6050_ADDR);
    wire->write(0x43);
    wire->endTransmission(false);
    wire->requestFrom((int)MPU6050_ADDR, 6);

    rx = wire->read() << 8 | wire->read();
    ry = wire->read() << 8 | wire->read();
    rz = wire->read() << 8 | wire->read();

    x += ((float)rx) / 65.5;
    y += ((float)ry) / 65.5;
    z += ((float)rz) / 65.5;
  }
  gyroXoffset = x / 3000;
  gyroYoffset = y / 3000;
  gyroZoffset = z / 3000;

  if(console){
    Console.println();
    Console.print("* Gyro offsets | x: ");Console.print(gyroXoffset,3);
    Console.print("| y: ");Console.print(gyroYoffset,4);
    Console.print("| z: ");Console.println(gyroZoffset,4);
		//delay(delayAfter);
	}
}

////////////////////////////////////////////////////////////////////////
void MPU6050::update(){
	wire->beginTransmission(MPU6050_ADDR);
	wire->write(0x3B);
	wire->endTransmission(false);
	wire->requestFrom((int)MPU6050_ADDR, 14);

  rawAccX = wire->read() << 8 | wire->read();
  rawAccY = wire->read() << 8 | wire->read();
  rawAccZ = wire->read() << 8 | wire->read();
  rawTemp = wire->read() << 8 | wire->read();
  rawGyroX = wire->read() << 8 | wire->read();
  rawGyroY = wire->read() << 8 | wire->read();
  rawGyroZ = wire->read() << 8 | wire->read();

  temp = (rawTemp + 12412.0) / 340.0;

  accX = -((float)rawAccX) / 16384.0;
  accY = ((float)rawAccY) / 16384.0;
  accZ = -((float)rawAccZ) / 16384.0;

  angleAccX = atan2(accY, accZ + abs(accX)) * 360 / 2.0 / PI;
  angleAccY = atan2(accX, accZ + abs(accY)) * 360 / -2.0 / PI;

  gyroX = -((float)rawGyroX) / 65.5;
  gyroY = ((float)rawGyroY) / 65.5;
  gyroZ = -((float)rawGyroZ) / 65.5;

  gyroX -= gyroXoffset;
  gyroY -= gyroYoffset;
  gyroZ -= gyroZoffset;

  interval = (millis() - preInterval) * 0.001;

  angleGyroX += gyroX * interval;
  angleGyroY += gyroY * interval;
  angleGyroZ += gyroZ * interval;

  angleX = (gyroCoef * (angleX + gyroX * interval)) + (accCoef * angleAccX);
  angleY = (gyroCoef * (angleY + gyroY * interval)) + (accCoef * angleAccY);
  angleZ = angleGyroZ;

  preInterval = millis();

}

////////////////////////////////////////////////////////////////////////
void MPU6050::writeMPU6050(byte reg, byte data){
  wire->beginTransmission(MPU6050_ADDR);
  wire->write(reg);
  wire->write(data);
  wire->endTransmission();
}

////////////////////////////////////////////////////////////////////////
byte MPU6050::readMPU6050(byte reg) {
  wire->beginTransmission(MPU6050_ADDR);
  wire->write(reg);
  wire->endTransmission(true);
  wire->requestFrom(MPU6050_ADDR, 1);
  byte data = wire->read();
  return data;
}

////////////////////////////////////////////////////////////////////////
bool MPU6050::isConnected()
{
    return read8(MPU6050_WHO_AM_I) == MPU6050_ADDR;
}

////////////////////////////////////////////////////////////////////////
byte MPU6050::read8(byte registerAddr)
{
    wire->beginTransmission(MPU6050_ADDR);
    wire->write(registerAddr);
    wire->endTransmission(false);
    wire->requestFrom(MPU6050_ADDR,1,true);

    return wire->read();
}

////////////////////////////////////////////////////////////////////////
int16_t MPU6050::read16(byte registerAddr)
{
    wire->beginTransmission(MPU6050_ADDR);
    wire->write(registerAddr);
    wire->endTransmission(false);
    wire->requestFrom(MPU6050_ADDR,2,true);

    // Concatenate the two bytes
    return wire->read()<<8 | wire->read();
}

////////////////////////////////////////////////////////////////////////
byte MPU6050::write8(byte registerAddr, byte value)
{
    wire->beginTransmission(MPU6050_ADDR);
    wire->write(registerAddr);
    wire->write(value);
    wire->endTransmission(true);
}

////////////////////////////////////////////////////////////////////////
void MPU6050::breakLine() {
  for (int i = 0; i <= SORO_CMDLENGHT ; i++) {
    Console.print(SORO_LINECHAR);
  }
  Console.println("");
}
