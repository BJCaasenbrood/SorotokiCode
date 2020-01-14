/***************************************************************************
   Communication SORO-board - Version 1.0.1
   by Brandon Caasenbrood: - b.j.caasenbrood@tue.nl

   This c++ file is licensed under the TU/e license
 ***************************************************************************/
#include "SORO.h"
#include "Communication.h"

char TemporaryChar;
instruction process;

/*Constructor (...)*********************************************************
 ***************************************************************************/
static inline char *getWord(char *buf) {
  static char *ptr = NULL;
  char *start, *scan;
  char term = ' ';

  if (buf != NULL) {ptr = buf;}

  while ((*ptr == ' ' || *ptr == '\t') && *ptr != '\0') {
    ptr++;
  }
  if (*ptr == '\0') {
    return NULL;
  }

  if (*ptr == '"' || *ptr == '\'') {
    term = *ptr;
    ptr++;
  }
  start = ptr;

  while (*ptr != '\0') {
    if (*ptr == '\\') {
      for (scan = ptr; *scan != '\0'; scan++) {
        *scan = *(scan + 1);
      }
      ptr++;
      continue;
    }
    if (*ptr == term || (term == ' ' && *ptr == '\t')) {
      *ptr = '\0';
      ptr++;
      return start;
    }
    ptr++;
  }
  if (ptr == start) {
    return NULL;
  }
  return start;
}

/*Constructor (...)*********************************************************
      Setup the communication process. Not required to be called.
 ***************************************************************************/
void Communication::beginCommunication() {

}

/*Constructor (...)*********************************************************
      Starting the communication process.
 ***************************************************************************/
void Communication::communicationProcess() {
  if (i2cScanActive){scanI2Cadress();}
}

/*Constructor (...)*********************************************************
      void process to wait for the pc data
 ***************************************************************************/
void Communication::awaitDataFromPC() {
  getDataFromPC();
  processData();
}

/***************************************************************************
     Process awaits data from the pc until the ENDMARKER (i.e., enter key) is
     detected. The process fills the inputBuffer[] with characters untill a
     full command is constructed. This command is parsed by void parseData().
 ***************************************************************************/
void Communication::getDataFromPC() {
  while (Console.available() > 0) {
    char tmp = Console.read();
    if (tmp == ENDMARKER) {
      if (!RunActive) {Console.print(tmp);}
      readInProgress = false;
      newDataFromPC = true;
      inputBuffer[bytesRecvd] = 0;
      parseData();
    }
    else if (!readInProgress) {
      bytesRecvd = 0;
      readInProgress = true;
    }

    if (readInProgress) {
      if (tmp == 8 || tmp == 127) {
        if (bytesRecvd > 0) {
          Console.print("\b \b");
          bytesRecvd--;
        }
      }
      else {
        if (!RunActive) {
          Console.print(tmp);
        }
        inputBuffer[bytesRecvd] = tmp;
        bytesRecvd ++;
      }
      if (bytesRecvd == buffSize) {
        bytesRecvd = buffSize - 1;
      }
    }

    TemporaryChar = tmp;

  }
}

/*Constructor (...)*********************************************************
     After the void awaitDataFromPC(), the bytes are converted to a String,
     constructing a command which can be processed by the void processData().
 ***************************************************************************/
void Communication::parseData() {

  char* strtokIndx;
  int ii = -1;
  char *w;

  argc = 0;
  w = getWord(inputBuffer);
  while ((argc < 20) && (w != NULL)) {
    argv[argc++] = w;
    w = getWord(NULL);
  }

//  messageFromPC = argv[0];
//  strtokIndx = strtok(inputBuffer, " ");
//  strcpy(messageFromPC, strtokIndx);

  while (strtokIndx != NULL) {
    if (ii == 0) {data1 = atof(strtokIndx);}
    if (ii == 1) {data2 = atof(strtokIndx);}
    if (ii == 2) {data3 = atof(strtokIndx);}

    strtokIndx = strtok(NULL, " ");
    ii++;
  }
}

/*Constructor (...)*********************************************************
     After the void awaitDataFromPC(), the bytes are converted to a String,
     constructing a command which can be processed by the void processData().
 ***************************************************************************/
void Communication::processData() {

  String str = String(argv[0]);
  float cmd[argc];
  float type[argc];

  for(int i = 0; i < argc; i++){
    cmd[i] = -1;
    type[i] = -1;
  }

  for(int i = 0; i < argc; i++){
  string arg = argv[i];
    if(isalpha(arg.at(0))){
      type[i] = MODE;
      if (arg == "help"){cmd[i] = MODE_HELP;};
      if (arg == "run"){cmd[i] = MODE_RUN;};
      if (arg == "stop"){cmd[i] = MODE_STOP;};
      if (arg == "gain"){cmd[i] = MODE_GAINS;};
      if (arg == "imu"){cmd[i] = MODE_IMU;};
    }
    else if(!isalpha(arg.at(0)) && !isdigit(arg.at(0))) {
      type[i] = WHAT;
      if (arg == "-o"){cmd[i] = ON;};
      if (arg == "-h"){cmd[i] = OFF;};
    }
    else if(isdigit(arg.at(0))){
      type[i] = VALUE;
      cmd[i] = stof(arg);
    }
  }

  if (newDataFromPC) {

  for(int i = 0; i < argc; i++){
    cout << " " << cmd[i] << " " << type[i] << endl;
  }
    
    newDataFromPC = false;

    if (str == "run") {
      setProcessMode(MODE_RUN);
      setProcessStatus(true);
    }
    else if (str == "stop") {
      setProcessMode(MODE_STOP);
      setProcessStatus(true);
    }
    else if (str == "gain") {
      setProcessMode(MODE_GAINS);
      setProcessStatus(true);
      setProcessData(data1, data2, data3);
    }
    else if (str == "imu") {
      setProcessMode(MODE_IMU);
      setProcessStatus(true);
      setProcessData(data1, data2, data3);
    }
    else if (str == "help") {
      setProcessMode(MODE_HELP);
      setProcessStatus(true);
    }
    else if (str == "test") {
      setProcessMode(MODE_TEST);
      setProcessStatus(true);
      setProcessData(data1, data2, data3);
    }
    else if (str == "read") {
      setProcessMode(MODE_READ_SEN);
      setProcessStatus(true);
    }
    else if (str == "mot") {
      setProcessMode(MODE_MOT);
      setProcessStatus(true);
      setProcessData(data1, data2, data3);
    }
    else if (str == "sol") {
      setProcessMode(MODE_SOL);
      setProcessStatus(true);
      setProcessData(data1, data2, data3);
    }
    else if (str == "i2c") {
      setProcessMode(MODE_I2CSCAN);
      setProcessStatus(true);
    }
    else if (str == "tca") {
      setProcessMode(MODE_TCASCAN);
      setProcessStatus(true);
    }
    else if (str == "setpoint") {
      setProcessMode(MODE_SETPOINT);
      setProcessStatus(true);
      setProcessData(data1, data2, data3);
    }
    else if (str == "a") {
      setProcessMode(MODE_A);
      setProcessStatus(true);
    }
    else if (str == "s") {
      setProcessMode(MODE_S);
      setProcessStatus(true);
    }
    else if (str == "d") {
      setProcessMode(MODE_D);
      setProcessStatus(true);
    }
    else if (str == "w") {
      setProcessMode(MODE_W);
      setProcessStatus(true);
    }
    else if (str == "p") {
      setProcessMode(MODE_GRIPPER);
      setProcessStatus(true);
    }
    else if (str == "-1") {
      Console.println("");
      Prompt();
    }
    else {
      Console.print("");
      cout << fg::red << argv[0] << " = unknown command" << style::reset 
      << " - check 'help' for built-in commands" << endl;
      argv[0] = "-1";
      Prompt();
    }
  }
}

/*Constructor (...)*********************************************************
     Set 'processes completed' confirmation
 ***************************************************************************/
void Communication::processCompleted() {
  if (getProcessStatus() && !RunActive) {
    //breakLine();
    Prompt();
  }
  setProcessData(-1, -1, -1);
  setProcessStatus(false);
  setProcessMode(-1);
}

/*Constructor (...)*********************************************************
     Prefroms backspace operation on the virtual command line with SOROTOKI.
 ***************************************************************************/
void Communication::setProcessData(float x, float y, float z) {
  process.N1 = x;
  process.N2 = y;
  process.N3 = z;
}

/*Constructor (...)*********************************************************
 ***************************************************************************/
void Communication::setProcessStatus(bool tmp) {
  process.active = tmp;
}

/*Constructor (...)*********************************************************
 ***************************************************************************/
void Communication::setCommunicationActive(bool tmp) {
  RunActive = tmp;
}

/*Constructor (...)*********************************************************
 ***************************************************************************/
void Communication::setI2CScanActive(bool tmp) {
  i2cScanActive = tmp;
}

/*Constructor (...)*********************************************************
 ***************************************************************************/
void Communication::setProcessMode(int tmp) {
  process.mode = tmp;
}

/*Constructor (...)*********************************************************
 ***************************************************************************/
bool Communication::getProcessStatus() {
  return process.active;
}

/*Constructor (...)*********************************************************
 ***************************************************************************/
void Communication::getArguments(char* outStr[20]) {
  for(int i=0; i < argc; ++i){
    outStr[i] = argv[i];
  }
}


/*Constructor (...)*********************************************************
 ***************************************************************************/
int Communication::getProcessMode() {
  return process.mode;
}

/*Constructor (...)*********************************************************
 ***************************************************************************/
int Communication::getProcessData(int address) {
  if (address == 0) {return process.N1;}
  else if (address == 1){return process.N2;}
  else {return process.N3;}
}

/*Constructor (...)*********************************************************
    Preform I2C scan to check all active devices
 ***************************************************************************/
void Communication::scanI2Cadress() {
  byte error, address;
  int nDevices;

  Console.println("");
  Console.println("Start scanning i2c ports...");

  nDevices = 0;
  for (address = 1; address < 127; address++ ) {
    Wire.beginTransmission(address);
    error = Wire.endTransmission();

    if (error == 0)
    {
      Console.print("* Found i2c device at address: 0x");
      if (address < 16)
        Console.print("0");
      Console.println(address, HEX);
      nDevices++;
    }
    else if (error == 4)
    {
      Console.print("* Unknown error at address: 0x");
      if (address < 16)
        Console.print("0");
      Console.println(address, HEX);
    }
  }
  if (nDevices == 0) {
    Console.println("* No i2c devices found...\n");
  }
  else {
    Console.print("* scan completed with (");
    Console.print(nDevices);
    Console.println(") devices found");
  }

  setI2CScanActive(false);
}
