#ifndef _COMMUNICATION_H_
#define _COMMUNICATION_H_
#define LIBRARY_VERSION	1.0.1

#include "src/interface.hpp"

#define ENDMARKER '\n'
#define buffSize 80

struct instruction
{
  int mode = -2;
  float N1 = -1;
  float N2 = -1;
  float N3 = -1;
  char *data[buffSize];
  bool active = false;
};

class Communication
{
  public:

      //Constants used in some of the functions below
    char *argv[20];
    int  argc;
    char inputBuffer[buffSize];
    byte bytesRecvd = 0;
    char messageFromPC[buffSize];
    bool readInProgress = false;
    bool newDataFromPC = false;
    bool dataProcessed = false;
    bool RunActive = false;
    bool i2cScanActive = false;
    float data1 = 0.0;
    float data2 = 0.0;
    float data3 = 0.0;

    void beginCommunication();
    void communicationProcess();
    void awaitDataFromPC();
    void getDataFromPC();
    void parseData();
    void processData();

    void processCompleted();
    void setProcessData(float tmp1, float tmp2, float tmp3);
    void setProcessStatus(bool tmp);
    void setProcessMode(int tmp);
    void setCommunicationActive(bool tmp);
    void setI2CScanActive(bool tmp);

    bool getProcessStatus();
    int getProcessMode();
    int getProcessData(int address);
    void getArguments(char* outStr[20]);

    void scanI2Cadress();

    instruction getProcess();

};

#endif
