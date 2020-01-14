#ifndef _I2CMULTIPLEX_HPP_
#define _I2CMULTIPLEX_HPP_

/*

*/
#include <Arduino.h>
#include <Wire.h>

#define TCAADDR 0x70

uint8_t tcaChannel = -1;

inline void tcaSelect(uint8_t i) {
  if (i > 7) return;
  tcaChannel = i;
  Wire.beginTransmission(TCAADDR);
  Wire.write(1 << i);
  Wire.endTransmission();  
}

inline uint8_t tcaGetChannel(){
	return tcaChannel;
}

#endif
