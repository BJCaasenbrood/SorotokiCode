import sys
import smbus

id   = int(sys.argv[1])
bus  = smbus.SMBus(id) # id indicates /dev/i2c-id
Ndev = 0

for device in range(128):

      try:
         bus.read_byte(device)
         print(hex(device))
         Ndev = Ndev + 1
      except: # exception if read_byte fails
         pass

#if Ndev == 0:
      #print("0")
