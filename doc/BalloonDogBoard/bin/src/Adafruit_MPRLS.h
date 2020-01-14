/*!
 * @file Adafruit_MPRLS.h
 *
 * Designed specifically to work with the MPRLS sensors from Adafruit
 * ----> https://www.adafruit.com/products/3965
 *
 * These sensors use I2C to communicate, 2 pins (SCL+SDA) are required
 * to interface with the breakout.
 *
 * Adafruit invests time and resources providing this open source code,
 * please support Adafruit and open-source hardware by purchasing
 * products from Adafruit!
 *
 * Written by Limor Fried/Ladyada for Adafruit Industries.  
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

#define MPRLS_DEFAULT_ADDR         (0x18) ///< Most common I2C address
#define MPRLS_STATUS_POWERED       (0x40) ///< Status SPI powered bit
#define MPRLS_STATUS_BUSY          (0x20) ///< Status busy bit
#define MPRLS_STATUS_FAILED        (0x04) ///< Status bit for integrity fail
#define MPRLS_STATUS_MATHSAT       (0x01) ///< Status bit for math saturation


/**************************************************************************/
/*! 
    @brief  Class that stores state and functions for interacting with MPRLS sensor IC
*/
/**************************************************************************/
class Adafruit_MPRLS {
 public:
  Adafruit_MPRLS(int8_t reset_pin = -1, int8_t EOC_pin = -1, 
		 uint8_t PSI_min = 0, uint8_t PSI_max = 25);

  boolean begin(uint8_t i2c_addr = MPRLS_DEFAULT_ADDR,
	     TwoWire *twoWire = &Wire);

  uint8_t readStatus(void);
  float   readPressure(void);
  
  void calcPressureOffset(bool console = false);
  
  float getPressureoffset(){ return PressureOffset; };
  
  void prepareRead(void);



 private:
  uint32_t readData(void);

  uint8_t _i2c_addr;
  int8_t _reset, _eoc;
  uint8_t _PSI_min, _PSI_max;
  float PressureOffset;

  TwoWire *_i2c;
};
