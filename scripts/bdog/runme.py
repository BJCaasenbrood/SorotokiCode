import sys, getopt
from SoftRobot  import SoftRobot

i2c = [1]

def main(argv):

    # setting up main class
    robot = SoftRobot(i2c=i2c) # port = 8888, frequency = 400 by default

    #robot.waitForClient()      # start TCP when we recieve start request from client
    robot.createProcesses()    # Initialize all the processes needed for I2C sensors, motors, TCP/IP comm
    robot.run()                # Start the processes
    robot.waitForProcesses()   # Wait for the processes to end/exit python scripts
    robot.socket_TCP.close()   # close TCP connections
    
if __name__ == "__main__":
    main(sys.argv[1:])
