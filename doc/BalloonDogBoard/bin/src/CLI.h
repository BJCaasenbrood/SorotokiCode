/*
 * Copyright (c) 2013, Majenko Technologies
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 * 
 * * Neither the name of Majenko Technologies nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _CLI_H
#define _CLI_H

#if ARDUINO >= 100
# include <Arduino.h>
#else
# include <WProgram.h>
#endif

#include <stdlib.h>

#define CLI_BUFFER 80 

class CLIServer;

class CLIClient : public Stream {
    private:
        Stream *dev;
        char input[CLI_BUFFER];
        int pos;
        char *prompt;
        boolean connected;

        boolean willEcho;

        uint8_t testConnected();

        void *_sessionData;

        void (*_redirect)(CLIClient *, char *, int);


    public:
        static const uint8_t IDLE = 0;
        static const uint8_t CONNECTED = 1;
        static const uint8_t DISCONNECTED = 2;

        CLIClient(Stream *d);
        virtual ~CLIClient();
        void echo(boolean e) { willEcho = e; }
        void printPrompt();
        int readline();
        int parseCommand();
        void setPrompt(const char *p);
        void setSessionData(void *data);
        void *getSessionData();

        void redirectStart(void (*func)(CLIClient *, char *, int)) { _redirect = func; }
        void redirectEnd() { _redirect = NULL; }

        size_t write(uint8_t);
        int available() { return dev->available(); }
        int read() { return dev->read(); }
        void flush() { dev->flush(); }
        int peek() { return dev->peek(); }
        
    friend class CLIServer;
};

#define CLI_COMMAND(X) int X(CLIClient *dev, int argc, char **argv)

typedef struct _CLIClientList {
    CLIClient *client;
    struct _CLIClientList *next;
} CLIClientList;

#define CLI_IS_PREFIX 0x01

typedef struct _CLICommand {
    char *command;
    uint8_t flags;
    int (*function)(CLIClient *, int, char **);
    struct _CLICommand *next;
} CLICommand;

class CLIServer : public Print {
    private:
        CLIClientList *clients;
        CLICommand *commands;
        char *prompt;
        boolean _caseSensitive;

        int (*_onConnect)(CLIClient *, int, char **);
        int (*_onDisconnect)(CLIClient *, int, char **);

    public:
        CLIServer();
        void setCaseInsensitive();
        void setCaseSensitive();
        boolean isCaseSensitive();
        void addCommand(const char *command, int (*function)(CLIClient *, int, char **));
        void addPrefix(const char *command, int (*function)(CLIClient *, int, char **));
        CLIClient *addClient(Stream *dev, void *data);
        CLIClient *addClient(Stream &dev, void *data);
        CLIClient *addClient(Stream *dev);
        CLIClient *addClient(Stream &dev);
        void removeClient(CLIClient *c);
        void removeClient(CLIClient &c);
        void removeClient(Stream *dev);
        void removeClient(Stream &dev);
        void onConnect(int (*function)(CLIClient *, int, char **));
        void onDisconnect(int (*function)(CLIClient *, int, char **));
        void process();
        void setDefaultPrompt(const char *p);
        void broadcast(char *);
#if (ARDUINO >= 100) 
        size_t write(uint8_t);
#else
        void write(uint8_t);
#endif

    friend class CLIClient;
};

extern CLIServer CLI;

#endif
