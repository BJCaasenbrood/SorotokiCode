#ifndef DEBUGGER_H
#define DEBUGGER_H
#include <iostream>

//#define DEBUG

void debug(const char *str){
	#ifdef DEBUG
		std::cout << "[ok]: " << str <<  << std::endl;
	#endif
}

#endif