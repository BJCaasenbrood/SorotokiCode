import multiprocessing
import socket
import struct
from ctypes import c_bool
import time

#---------------------------------------------------------------------
#- main SoftRC class -------------------------------------------------
class baseSoftRobot:
    def __init__(self,nSensors,port):
        self.nSensors = nSensors

        ## Set up TCP server
        # Set up socket
        host = ''
        server_address = (host, port) 
        self.clients = []
        self.control_processes = []
        self.comm_processes = []
        self.clients_addresses = []
        self.buffersize = 8*self.nMotors
        self.socket_TCP = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket_TCP.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.socket_TCP.bind(server_address)

        ## Multiprocessing
        if self.nSensors > 0:
            self.sensorsUpdated = multiprocessing.Array(c_bool,[False]*self.nSensors)
        self.stopFlag = multiprocessing.Value(c_bool,False)
        
        print("port: ",self.port)
        print("i2c:  ",self.port)

#---------------------------------------------------------------------
#- main SoftRC clas -------------------------------------------------
    def waitForClient(self):
        """Automatically accept one incoming client"""
        self.socket_TCP.listen(1)
        print("TCP server is running. Waiting for client...")
        client, client_address = self.socket_TCP.accept()
        print("Client ",client_address," connected. Comm's ready!")
        client.settimeout(1)
        time.sleep(0.9)
        self.clients.append(client)
        self.clients_addresses.append(client_address)

    def repeatedlySend(self):
        """Send data to all clients until stopFlag == True"""
        while not self.stopFlag.value:
            if all(self.sensorsUpdated):
                try:
                    data = struct.pack('d'*len(self.sensorsValues),*self.sensorsValues)
                    for client in self.clients:
                        client.sendall(data)

                    for i in range(self.nSensors):
                        self.sensorsUpdated[i] = False
                except Exception as e:
                    print('Send failed')
                    print(e)
                    self.stopFlag.value = True
                    break
        self.socket_TCP.close()

    def receive(self):
        while not self.stopFlag.value:
            try:
                raw = self.clients[0].recv(self.buffersize)
                 # print('Data received')
                unpackedData = struct.unpack('d'*int(self.buffersize/8),raw)
                
                for i in range(int(self.buffersize/8)):
                    self.motorsValues[i]=unpackedData[i]
            except Exception as e:
                print('Error in receive: ',e)
                self.stopFlag.value = True
                
        self.socket_TCP.close()

    def readSensors(self):
        raise NotImplementedError("readSensors must be implemented in the child class")
        
    def controlActuators(self):
        raise NotImplementedError("controlActuators must be implemented in the child class")

    def createProcesses(self):
        ## Processes to read sensors
        for i in range(self.nSensors):
            self.control_processes.append(multiprocessing.Process(target=self.readSensors, args=(i,)))

        ## Process to control the motors
        self.control_processes.append(multiprocessing.Process(target=self.controlActuators))

        ## Processes for TCP/IP comm
        self.comm_processes.append(multiprocessing.Process(target=self.repeatedlySend))
        self.comm_processes.append(multiprocessing.Process(target=self.receive))

    def run_control(self):
        for p in self.control_processes:
            p.start()
            
    def run_comm(self):
        for p in self.comm_processes:
            p.start()
            
    def run(self):
    	self.run_control()
    	self.waitForClient()
    	self.run_comm()

    def waitForProcesses(self):
        for p in self.processes:
            p.join()
            
            
#---------------------------------------------------------------------
#- base template for sensor class ---------------------------------------
class baseSensor:
    def __init__(self, sensorInstance):
        self.instance = sensorInstance

    def readSensor(self):
        raise NotImplementedError("readSensor must be implemented in the child class")
