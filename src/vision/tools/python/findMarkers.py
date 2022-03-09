import numpy as np
import cv2, PIL
from cv2 import aruco
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

#matplotlib nbagg

aruco_dict = aruco.Dictionary_get(aruco.DICT_6X6_250)

frame = cv2.imread("tmp.png")
#plt.figure()
#plt.imshow(frame)
#plt.show()

## time
gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
aruco_dict = aruco.Dictionary_get(aruco.DICT_6X6_250)

parameters =  aruco.DetectorParameters_create()

corners, ids, rejectedImgPoints = aruco.detectMarkers(gray, aruco_dict, parameters=parameters)

frame_markers = aruco.drawDetectedMarkers(frame.copy(), corners, ids)

#plt.figure()
#plt.imshow(frame_markers)

with open('readme.txt', 'w') as f:

    for i in range(len(ids)):
        c = corners[i][0]
        for j in range(0, 4):
            f.write(str(c[j,0]))
            f.write(', ')
            f.write(str(c[j,1]))
            f.write('\n')
	    
    #plt.plot([c[:, 0].mean()], [c[:, 1].mean()], "o", label = "id={0}".format(ids[i]))

#plt.legend()
#plt.show()
