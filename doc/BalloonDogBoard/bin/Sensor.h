#ifndef _SENSOR_H_
#define _SENSOR_H_

#define FREQLP 1

struct SensorState{
  float IMUx = 0;
  float IMUy = 0;
  float IMUz = 0;
  float IMUax = 0;
  float IMUay = 0;
  float IMUaz = 0;
  float IMUgx = 0;
  float IMUgy = 0;
  float IMUgz = 0;
  float P1 = 0;
  float P2 = 0;
  float P3 = 0;
};

class Sensor 
{
  public:

    const int ParseImuX = 0;
    const int ParseImuY = 1;
    const int ParseImuZ = 2;
    const int ParseP1   = 3;
    const int ParseP2   = 4;
    const int ParseP3   = 5;

    bool sensorActive = false;
    bool ImuConnection = true;
    bool MprConnection = true;
    bool tcaScanActive = false;

    //Vector LastReading;

    void beginSensor();
    void sensorProcess();
    void setSensorActive(bool tmp);
    void setImu(bool tmp);
    void setHX711(bool tmp);
    void startIMUConnection();
    void startMPRConnection(int tcaID = 0);
    bool getSensorProcessStatus();
    void initSensors();
    void showSensor();
    void readSensor(bool enablePrint);
    void readPressure();
    void tarePressure();
    void showQuat();
    
    void readIMU();
    void scanTCAAdress();
    void setTCAScanActive(bool tmp);

    float angleWrap(float x);

    void parseDataSensor(int msg, float data);
    float getDataSensor(int msg);
  
};

#endif
