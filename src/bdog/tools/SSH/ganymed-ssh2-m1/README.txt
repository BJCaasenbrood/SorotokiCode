
Ganymed SSH-2 for Java - build 250
==================================

http://www.cleondris.ch/ssh2

Ganymed SSH-2 for Java is a library which implements the SSH-2 protocol in pure Java
(tested on J2SE 1.4.2, 5 and 6). It allows one to connect to SSH servers from within
Java programs. It supports SSH sessions (remote command execution and shell access),
local and remote port forwarding, local stream forwarding, X11 forwarding, SCP and SFTP.
There are no dependencies on any JCE provider, as all crypto functionality is included.

Ganymed SSH-2 for Java was originally developed for the Ganymed replication project
and a couple of other projects at the IKS group at ETH Zurich (Switzerland).

This distribution contains the source code, examples, javadoc and the FAQ.
It also includes a pre-compiled jar version of the library which is ready to use.

- Please read the included LICENCE.txt
- Latest changes can be found in HISTORY.txt

The latest version of the FAQ is available on the website.

Please feel free to contact us. We welcome feedback of any kind!
Contact: Christian Plattner, christian.plattner@cleondris.ch

Christian Plattner
March 2010

(An unsupported ReadMatlab function was added to the SFTPv3Client to support SFTP-get files in Matlab. The original file doesn't allow pass by reference for byte arrays, see http://www.mathworks.com/matlabcentral/newsreader/view_thread/71084 for more info).
David S. Freedman
January 2013