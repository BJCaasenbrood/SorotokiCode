#---------------------------------------------------------------------
#- base imports (PLEASE DONT MODIFY) ---------------------------------
import os.path
import math
import sys
import time
import multiprocessing
import numpy as np

# base soft robot class
from baseSoftRobot import baseSoftRobot, baseSensor

# base multiple i2c busses on the RPi 
from adafruit_extended_bus import ExtendedI2C as I2C

# base i2c libaries for Analog-Digital and Digital-Analog converters
import adafruit_mcp4725
import adafruit_ads1x15.ads1015 as ADS
from adafruit_ads1x15.analog_in import AnalogIn
from adafruit_ads1x15.ads1x15 import Mode

#---------------------------------------------------------------------
#- adding additional libaries for custom sensors (please refer to the 
#  extensive libaries of Adafruit or Sparkfun!).

# examples:
# from adafruit_mprls import MPRLS
# from adafruit_mcp9808 import MCP

#---------------------------------------------------------------------

class VeabSensor(baseSensor):
    def __init__(self, i2c=1, address=(0x48), RATE = 3300):
        ads = ADS.ADS1015(I2C(i2c),address=address)
        # ADC Configuration
        ads.mode      = Mode.CONTINUOUS	 # read last measurement
        ads.data_rate = RATE             # set measurement rate
        ads.gain      = 1                # set adc gain
        
        chan = AnalogIn(ads, ADS.P0)     # set channel read
        super().__init__(chan)           # init adc
    def readSensor(self):
        return self.instance.voltage*3
        
        
#---------------------------------------------------------------------
#- main SoftRC class -------------------------------------------------
class SoftRobot(baseSoftRobot):
    def __init__(self,i2c=[1], port = 12345, frequency = 400):

        self.port      = port
        self.nSensors  = len(i2c)
        self.channels  = i2c
        self.frequency = frequency
        self.period    = 1.0/frequency
        
        ## Set-up multiprocessing for different sensors
        self.sensorsValues = multiprocessing.Array('d',[0.0]*(self.nSensors))
        self.sensors = []
        for i in range(self.nSensors):
            self.sensors.append(VeabSensor(i2c=self.channels[i]))

        ## Set up actuators
        self.nMotors = len(i2c)
        self.motors = []
        for i in range(self.nMotors):
            self.motors.append(adafruit_mcp4725.MCP4725(I2C(self.channels[i]),address=0x60))
        self.motorsValues = multiprocessing.Array('d',[0.5]*self.nMotors)
        #self.setActuators(0.9)
        ## Call __init__ of the parent class (Set up multi-processes and TCP comm)
        super().__init__(self.nSensors, self.port)
        
#---------------------------------------------------------------------        
#    def addMPR(self,i2c = 1):
#        if self.IMUadded == False:
#            self.nSensors = self.nSensors+1
#            self.sensors.append(PressureSensor())
#            self.sensorsValues = multiprocessing.Array('d',[0.0]*(self.nSensors))
#            super().__init__(self.nSensors, self.port)
#            print("Adding one MPR sensor on I2C channel ", i2c)
#        else:
#            raise Exception("IMU added, cannot add pressure sensors.")

#    def addIMU(self):
#        self.sensorsValues = multiprocessing.Array('d',[0.0]*(self.nSensors+3))
#        self.nSensors = self.nSensors+1
#        super().__init__(self.nSensors, self.port)
#        self.processes.append(multiprocessing.Process(target=self.readIMU))

#---------------------------------------------------------------------
    def readSensors(self, index):
        while not self.stopFlag.value:
            try:
                self.sensorsValues[index]  = self.sensors[index].readSensor()
                self.sensorsUpdated[index] = True
                
                time.sleep(self.period - time.time() * self.frequency % 1 / self.frequency)
            except Exception as e:
                print('Error in readSensors:',e)
                self.stopFlag.value = True

#---------------------------------------------------------------------
    def setActuators(self,value = 0.5):
        for i,p in enumerate(range(self.nMotors)):
            self.motors[p].normalized_value = value

    def controlActuators(self):
        while not self.stopFlag.value:
            try:           
                for i,p in enumerate(range(self.nMotors)):
                    self.motors[p].normalized_value = self.motorsValues[i]

                time.sleep(self.period - time.time() * self.frequency % 1 / self.frequency)
            except Exception as e:
                print('Error in control Actuators:',e)
                self.stopFlag.value = True
        self.setActuators(0.5)  
#---------------------------------------------------------------------        
