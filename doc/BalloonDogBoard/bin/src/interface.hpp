#ifndef _INTERFACE_HPP_
#define _INTERFACE_HPP_

#define SORO_LINECHAR '-'
#define SORO_STOPCHAR '~'
#define SORO_CMDLENGHT 62
#define BASE_PROMPT " "

#include "rang.hpp"

using namespace rang;
using namespace std;

/**************************************************************************/
// New line
/**************************************************************************/
inline void newLine(){
	cout << " " << endl;
}

/**************************************************************************/
// Break line
/**************************************************************************/
inline void breakLine() {
  for (int i = 0; i <= SORO_CMDLENGHT ; i++) {
    Console.print(SORO_LINECHAR);
  }
  newLine();
}

/**************************************************************************/
// Print out
/**************************************************************************/
inline void printout(String title) {
	Console.print(BASE_PROMPT);
	Console.print(title);
	newLine();
}

/**************************************************************************/
// Print title
/**************************************************************************/
inline void titleLine(String title) {
  cout << " " << fg::yellow << flush;
  Console.print(title);
  cout << " " << style::reset << flush;

  int N = title.length();
  for (int i = 0; i <= SORO_CMDLENGHT - 2 - N ; i++) {
    Console.print(SORO_LINECHAR);
  }
  newLine();
}

/**************************************************************************/
// Double space 
/**************************************************************************/
inline void DoubleSpace() {
  newLine();
}

/**************************************************************************/
// Prompt symbol
/**************************************************************************/
inline void Prompt() {
  cout << fg::yellow << "@" << style::bold << "sorotoki"
  << style::reset << ":" << fg::blue << "~ $ " 
  << style::reset << flush;
}

/**************************************************************************/
// LOGO SOROTOKIs
/**************************************************************************/
inline void DisplayLogo() {
  DoubleSpace();
  breakLine();
  breakLine();
  DoubleSpace();
  Console.println("Sorotoki v.2.1.1 (2018.01.05 14:20+0000) Built-in shell (ash)");
  Console.println("Enter 'help' for a list of built-in commands");
  DoubleSpace();
#ifndef SORO_NOLOGO
  Console.print("███████╗ ██████╗ ██████╗  ██████╗ ████████╗ ██████╗ ██╗  ██╗██╗\n");
  Console.print("██╔════╝██╔═══██╗██╔══██╗██╔═══██╗╚══██╔══╝██╔═══██╗██║ ██╔╝██║\n");
  Console.print("███████╗██║   ██║██████╔╝██║   ██║   ██║   ██║   ██║█████╔╝ ██║\n");
  Console.print("╚════██║██║   ██║██╔══██╗██║   ██║   ██║   ██║   ██║██╔═██╗ ██║\n");
  Console.print("███████║╚██████╔╝██║  ██║╚██████╔╝   ██║   ╚██████╔╝██║  ██╗██║\n");
  Console.print("╚══════╝ ╚═════╝ ╚═╝  ╚═╝ ╚═════╝    ╚═╝    ╚═════╝ ╚═╝  ╚═╝╚═╝\n");
  Console.println("");
#endif
  titleLine(" Sorotoki Firmware (Version 0.0.1 rev.1)");
}

#endif
