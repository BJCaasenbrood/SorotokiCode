import multiprocessing
import socket
import struct
from ctypes import c_bool

class baseSoftRobot:
    def __init__(self,nSensors,port):
        self.nSensors = nSensors

        ## Set up TCP server
        # Set up socket
        host = ''
        server_address = (host, port) 
        self.buffersize = 8*self.nMotors
        self.socket_TCP = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket_TCP.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.socket_TCP.bind(server_address)
        self.clients = []
        self.clients_addresses = []

        ## Multiprocessing
        if self.nSensors > 0:
            self.sensorsUpdated = multiprocessing.Array(c_bool,[False]*self.nSensors)
        self.stopFlag = multiprocessing.Value(c_bool,False)
        self.processes = []

    def waitForClient(self):
        """Automatically accept one incoming client"""
        self.socket_TCP.listen(1)
        print("TCP server is running. Waiting for client...")
        client, client_address = self.socket_TCP.accept()
        print("Client ",client_address," connected. Comm's ready!")
        client.settimeout(0.5)
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
            self.processes.append(multiprocessing.Process(target=self.readSensors, args=(i,)))

        ## Process to control the motors
        self.processes.append(multiprocessing.Process(target=self.controlActuators))

        ## Processes for TCP/IP comm
        self.processes.append(multiprocessing.Process(target=self.repeatedlySend))
        self.processes.append(multiprocessing.Process(target=self.receive))

    def run(self):
        for p in self.processes:
            p.start()

    def waitForProcesses(self):
        for p in self.processes:
            p.join()

class baseSensor:
    def __init__(self, sensorInstance):
        self.instance = sensorInstance

    def readSensor(self):
        raise NotImplementedError("readSensor must be implemented in the child class")
