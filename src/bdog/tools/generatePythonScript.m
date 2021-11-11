function FID = generatePythonScript(board,FID)

numProcesses = board.get('Ndev');

fprintf(FID,'import busio \n');
fprintf(FID,'from adafruit_extended_bus import ExtendedI2C as I2C \n');
fprintf(FID,'from adafruit_mprls import MPRLS \n');
fprintf(FID,'import multiprocessing \n');
fprintf(FID,'import time \n');
fprintf(FID,'import board \n');
fprintf(FID,'import socket \n');
fprintf(FID,'import sys \n');
fprintf(FID,'import struct \n');
fprintf(FID,'from ctypes import c_bool \n');
fprintf(FID,'\n');

%% setting up main
fprintf(FID,"if __name__ == '__main__': \n");
fprintf(FID,"\t#set up host \n");
fprintf(FID,"\thost = '' \n");
fprintf(FID,"\tport = 8889 \n");
fprintf(FID,"\tbuffersize = 8\n");
fprintf(FID,"\tserver_address = (host, port)\n");
fprintf(FID,"\tsocket_TCP = socket.socket(socket.AF_INET, socket.SOCK_STREAM)\n");
fprintf(FID,"\tsocket_TCP.bind(server_address)\n");
fprintf(FID,"\tsocket_TCP.listen(1)\n\n");

fprintf(FID,"\tprint('TCP server is running. Waiting for client...')\n");
fprintf(FID,"\tsocket_TCP.bind(server_address)\n");
fprintf(FID,"\tsocket_TCP.listen(1)\n\n");

fprintf(FID,"\t#set up I2C busses\n");

dev = [board.Sensor;board.Actuator];
i2cbuslist = unique(vertcat(dev{:,2}));
for ii = 1:length(i2cbuslist)
    [SCL,SDA] = i2cbusPi(i2cbuslist(ii));
    fprintf(FID,['\ti2c',num2str(i2cbuslist(ii)),' = busio.I2C(',...
        num2str(SCL),',',num2str(SDA),',frequency= 400000)\n']);
end

fprintf(FID,"\n\t#set up sensor libaries\n");

for ii = 1:length(board.Sensor)
    dev = board.Sensor(ii,:);
    if strcmp(dev{1},'0x18')
        fprintf(FID,['\tmpr',num2str(ii),' = MPRLS(','i2c',num2str(dev{2}),...
            ',psi_min=0, psi_max=25)\n']);
    end
    %fprintf(FID,['\t\tp',num2str(ii),'.start()\n']);
end

fprintf(FID,"\n\t#set up sensor libaries\n");

for ii = 1:length(board.Actuator)
    dev = board.Actuator(ii,:);
    if strcmp(dev{1},'0x56')
        %fprintf(FID,['\tmpr',num2str(ii),' = MPRLS(','i2c',num2str(dev{2}),...
        %    ',psi_min=0, psi_max=25)\n']);
    end
    %fprintf(FID,['\t\tp',num2str(ii),'.start()\n']);
end

fprintf(FID,"\n");
fprintf(FID,"\ttry: \n");
for ii = 1:numProcesses
    fprintf(FID,['\t\tp',num2str(ii),'.start()\n']);
end

fprintf(FID,"\n");
fprintf(FID,"\texcept: \n");
for ii = 1:numProcesses
    fprintf(FID,['\t\tp',num2str(ii),'.terminate()\n']);
end
% 
%     print("TCP server is running. Waiting for client...")
%     client, client_address = socket_TCP.accept()
%     print(client_address," has connected.")
% 
%     # Start I2C
%     i2c = busio.I2C(1,0,frequency= 400000)
%     mpr = MPRLS(i2c, psi_min=0, psi_max=25)
% 
%     # Multiprocessing
%     sensorUpdated = multiprocessing.Value(c_bool,False)
%     sensorValue = multiprocessing.Value('d',0.0)
%     p1 = multiprocessing.Process(target=readsensor, args=(mpr,sensorValue,sensorUpdated))
%     p2 = multiprocessing.Process(target=senddata, args=(client,sensorValue,sensorUpdated))
%     
%     try:
%         p1.start()
%         p2.start()
% 
%     except:
%         print('Interrupted')
%         p1.terminate()
%         p2.terminate()
%         sys.exit(1)
%         pass # Handle error here. 
end

