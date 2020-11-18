clr;
%% assign free DOF
ip = '169.254.168.25';
usr = 'pi';
pwd = 'raspberry';

board = Bdog(ip,usr,pwd);

%% generate dynamic model
board = board.connect();
board = board.run();

%% get data
A = board.read('soro.log');

figure;
hold on;
plot(A(:,1),A(:,2));
plot(A(:,1),A(:,3));
plot(A(:,1),A(:,4));

figure;
hold on;
plot(A(:,1),A(:,5));
plot(A(:,1),A(:,6));
plot(A(:,1),A(:,7));

board = board.disconnect();