/*!
 * @file Adafruit_MPRLS.cpp
 *
 * @mainpage Adafruit MPRLS Pressure sensor
 *
 * @section intro_sec Introduction
 *
 * Designed specifically to work with the MPRLS sensor from Adafruit
 * ----> https://www.adafruit.com/products/3965
 *
 * These sensors use I2C to communicate, 2 pins (SCL+SDA) are required
 * to interface with the breakout.
 *
 * Adafruit invests time and resources providing this open source code,
 * please support Adafruit and open-source hardware by purchasing
 * products from Adafruit!
 *
 * @section dependencies Dependencies
 *
 *
 * @section author Author
 *
 * Written by Limor Fried/Ladyada for Adafruit Industries.
 *
 * @section license License
 *
 * MIT license, all text here must be included in any redistribution.
 *
 */

#if (ARDUINO >= 100)
 #include "Arduino.h"
#else
 #include "WProgram.h"
#endif
#include "Wire.h"

#include "Adafruit_MPRLS.h"

/**************************************************************************/
/*! 
    @brief constructor initializes default configuration value
    @param reset_pin Optional hardware reset pin, default set to -1 to skip
    @param EOC_pin Optional End-of-Convert indication pin, default set to -1 to skip
    @param PSI_min The minimum PSI measurement range of the sensor, default 0
    @param PSI_max The maximum PSI measurement range of the sensor, default 25
*/
/**************************************************************************/
Adafruit_MPRLS::Adafruit_MPRLS(int8_t reset_pin, int8_t EOC_pin, 
			       uint8_t PSI_min, uint8_t PSI_max) {

    _reset = reset_pin;
    _eoc = EOC_pin;
    _PSI_min = PSI_min;
    _PSI_max = PSI_max;
}

/**************************************************************************/
/*! 
    @brief  setup and initialize communication with the hardware
    @param i2c_addr The I2C address for the sensor (default is 0x18)
    @param twoWire Optional pointer to the desired TwoWire I2C object. Defaults to &Wire
    @returns True on success, False if sensor not found
*/
/**************************************************************************/
boolean Adafruit_MPRLS::begin(uint8_t i2c_addr, TwoWire *twoWire) {
  _i2c_addr = i2c_addr;
  _i2c = twoWire;

  _i2c->begin();

  if (_reset != -1) {
    pinMode(_reset, OUTPUT);
    digitalWrite(_reset, HIGH);
    digitalWrite(_reset, LOW);
    delay(10);
    digitalWrite(_reset, HIGH);
  }
  if (_eoc != -1) {
    pinMode(_eoc, INPUT);
  }

  //delay(10); // startup timing

  //Serial.print("Status: ");
  uint8_t stat = readStatus();
  //Serial.println(stat);
  return !(stat & 0b10011110);
}

/**************************************************************************/
/*! 
    @brief Read and calculate the pressure
    @returns The measured pressure, in hPa on success, NAN on failure
*/
/**************************************************************************/
float Adafruit_MPRLS::readPressure(void) {
  uint32_t raw_psi = readData();
  if (raw_psi == 0xFFFFFFFF) {
    return NAN;
  }

  // All is good, calculate the PSI and convert to kPA
  // use the 10-90 calibration curve
  float psi = (raw_psi - 0x19999A) * (_PSI_max - _PSI_min);
  psi /= (float) (0xE66666 - 0x19999A);
  psi += _PSI_min;
    // convert PSI to kPA
  return psi * 6.8947572932;
}


/**************************************************************************/
/*! 
    @brief Read 24 bits of measurement data from the device
    @returns -1 on failure (check status) or 24 bits of raw ADC reading
*/
/**************************************************************************/
uint32_t Adafruit_MPRLS::readData(void) {
/*   _i2c->beginTransmission(_i2c_addr);
  _i2c->write(0xAA);   // command to read pressure
  _i2c->write(0x0);
  _i2c->write(0x0);
  _i2c->endTransmission(); */
  
  // Use the gpio to tell end of conversion
  if (_eoc != -1) {
    while (!digitalRead(_eoc)) {
      //delay(10);
    }
  } else {
    // check the status byte
    uint8_t stat;
    while ((stat = readStatus()) & MPRLS_STATUS_BUSY) {
      //Serial.print("Status: "); Serial.println(stat, HEX);
      //delay(10);
    }
  }
  _i2c->requestFrom(_i2c_addr, (uint8_t)4);
  
  uint8_t status = _i2c->read();
  if (status & MPRLS_STATUS_MATHSAT) {
    return 0xFFFFFFFF;
  }
  if (status & MPRLS_STATUS_FAILED) {
    return 0xFFFFFFFF;
  }

  uint32_t ret;
  ret = _i2c->read(); ret <<= 8;
  ret |= _i2c->read(); ret <<= 8;
  ret |= _i2c->read();
  
  return ret;
}
/**************************************************************************/
/*! 
    @brief Read just the status byte, see datasheet for bit definitions
    @returns 8 bits of status data
*/
/**************************************************************************/
void Adafruit_MPRLS::prepareRead(void) {
_i2c->beginTransmission(_i2c_addr);
_i2c->write(0xAA);   // command to read pressure
_i2c->write(0x0);
_i2c->write(0x0);
_i2c->endTransmission();
}

/**************************************************************************/
/*! 
    @brief Read just the status byte, see datasheet for bit definitions
    @returns 8 bits of status data
*/
/**************************************************************************/
uint8_t Adafruit_MPRLS::readStatus(void) {
  _i2c->requestFrom(_i2c_addr, (uint8_t)1);

  return _i2c->read();
}


/**************************************************************************/
/*! 
    @brief Read just the status byte, see datasheet for bit definitions
    @returns 8 bits of status data
*/
/**************************************************************************/
void Adafruit_MPRLS::calcPressureOffset(bool console) {
	float x = 0;

  //delay(delayBefore);
	if(console){
    Console.println();
    Console.print("* MPRLS: calibration in progress");
  }
  
  for(int i = 0; i < 3000; i++){
    if(console && i % 600 == 0){
      Console.print(".");
    }
    x += readData();
  }
  
  PressureOffset = x / 3000;

  if(console){
    Console.println();
    Console.print("* Pressure offset | p: ");Console.print(PressureOffset,3);
		//delay(delayAfter);
	}
}


