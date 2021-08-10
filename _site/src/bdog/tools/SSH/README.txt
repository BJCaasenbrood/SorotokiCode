
Ganymed SSH-2 (Java) for for Matlab - m1
===========================================

See the example folder for many examples on how one can 
use the SSH-2 Matlab library. 

Known Issues: Transferred files to the PC do not seem to be released 
from Matlab until the application is closed.

- Please read the included LICENCE.txt
- Latest changes can be found in HISTORY.txt

Please feel free to contact me. I welcome feedback of any kind!
Contact: David S. Freedman, dfreedma@bu.edu

David Freedman 
Jan 1st 2013

(An unsupported ReadMatlab function was added to the SFTPv3Client to 
support SFTP-get files in Matlab. The original file doesn't allow 
pass by reference for byte arrays, 
see http://www.mathworks.com/help/matlab/matlab_external/handling-data-returned-from-a-java-method.html 
and http://www.mathworks.com/matlabcentral/newsreader/view_thread/71084 
for more info. Custom (and compiled) java source is included but not enabled by default).