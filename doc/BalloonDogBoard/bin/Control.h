#ifndef _CONTROL_H_
#define _CONTROL_H_

#define STATEDIM 2
#define DECIMEL 4
#define ParseMsgT -1
#define ParseMsgY0 0
#define ParseMsgY1 1
#define ParseMsgY2 2
#define ParseMsgU0 3
#define ParseMsgU1 4
#define ParseMsgU2 5

class Control
{
  public:
    float eps = 0.0001;
    bool runActive = false;
    float dt = 0.0;
    float t = 0.0;
    float Kp = 20;
    float Ki = 1;
    float Kd = 0.0001;

    void beginControl();
    void controlProcess();
    bool getControlProcessStatus();
    void updateTime();
    void setUpdateRateControl(float tmp);
    void setActive(bool tmp);

    void setControlGain(float tmp, int id);
    
    void initStateSpace();
    void showStates();
    void showControllerTime();

    void parseDataControl(int msg, float data);
    float getDataControl(int msg);
};

#endif
