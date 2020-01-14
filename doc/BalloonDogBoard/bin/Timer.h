#ifndef Timer_h
#define Timer_h
#ifdef __cplusplus

#if ARDUINO >= 100
#include "Arduino.h"
#else
#include "WProgram.h"
#endif

class Timer
{
private:
	unsigned long us;
public:
	Timer(void) { us = micros(); }
	Timer(unsigned long val) { us = micros() - val; }
	Timer(const Timer &orig) { us = orig.us; }
	operator unsigned long () const { return micros() - us; }
	Timer & operator = (const Timer &rhs) { us = rhs.us; return *this; }
	Timer & operator = (unsigned long val) { us = micros() - val; return *this; }
	Timer & operator -= (unsigned long val)      { us += val ; return *this; }
	Timer & operator += (unsigned long val)      { us -= val ; return *this; }
	Timer operator - (int val) const           { Timer r(*this); r.us += val; return r; }
	Timer operator - (unsigned int val) const  { Timer r(*this); r.us += val; return r; }
	Timer operator - (long val) const          { Timer r(*this); r.us += val; return r; }
	Timer operator - (unsigned long val) const { Timer r(*this); r.us += val; return r; }
	Timer operator + (int val) const           { Timer r(*this); r.us -= val; return r; }
	Timer operator + (unsigned int val) const  { Timer r(*this); r.us -= val; return r; }
	Timer operator + (long val) const          { Timer r(*this); r.us -= val; return r; }
	Timer operator + (unsigned long val) const { Timer r(*this); r.us -= val; return r; }
};

#endif 
#endif
