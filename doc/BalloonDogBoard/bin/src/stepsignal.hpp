#ifndef _STEPSIGNAL_HPP_
#define _STEPSIGNAL_HPP_

/**************************************************************************/
// New line
/**************************************************************************/
inline float stepsignal(float t){
	if(t <= 9 && t > 0){
		return 2000 + 1000*sin(0.6*t-1.1);
	}
	else{
		return 0;
	}
}

#endif
