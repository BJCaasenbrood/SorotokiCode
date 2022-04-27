clr;
%%
port = 8888;
ip = 'localhost';

t = tcpip(ip,port,'NetworkRole', 'client');
fopen(t)

A = fread(t, 10);
%%
fclose(t)