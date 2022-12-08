classdef Bdog < handle

    properties (Access = public)
        t;
        U;      % control signal
        Y;      % measurement signal
        P0;     %
        Ip;
        Usr;
        Pwd;
        Port;
        Log;
        Sensor;
        Actuator;
        Frequency;
        NVeab;
    end
    
    properties (Access = private)
        out;
        Ndev;
        pyFile;
        outFile;
        conFile;
        autoConnect;
        autoCompile;
        autoTransfer;
        autoi2cCheck;
        isConnected;
        i2cBusses;
        SSH;
        TCPclient;
        LoopCounter;
        LoopIndex;
    end
    
%--------------------------------------------------------------------------
methods  
%-------------------------------------------------------- ballong-dog class 
function obj = Bdog(usr,ip,pwd,varargin)
    
    obj.Port = 8888;
    obj.Ip  = ip;
    obj.Usr = usr;
    obj.Pwd = pwd;
    
    obj.t       = 0;
    obj.Ndev    = 0;
    obj.pyFile  = 'runme.py';
    obj.outFile = 'soro.log';
    obj.conFile = 'config.txt';
    obj.autoCompile = false;
    obj.autoConnect = false;
    obj.i2cBusses   = [1,3,4,5,6];
    obj.isConnected = false;
    obj.autoConnect = true;
    obj.LoopCounter = 0;
    obj.LoopIndex = uint8(1e6);
    obj.Frequency = 400;
    obj.Log = struct('y',[],'t',[]);
    obj.NVeab = 1;
    
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
    
    if isempty(obj.P0)
        obj.P0 = zeros(2*obj.NVeab,1);
    end
    
    obj.U = pdac(obj.P0);
    
    % setup TCP client
    %client = tcpip(obj.Ip,obj.Port,'NetworkRole', 'client');
    client = tcpclient(obj.Ip,obj.Port);
    client.ByteOrder = 'little-endian';
    obj.TCPclient = client;
      
    if obj.autoConnect
        obj = connect(obj); 
    end
end
%---------------------------------------------------------------------- get     
function varargout = get(Bdog,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = Bdog.(varargin{ii});
        end
    else
        varargout = Bdog.(varargin);
    end
end     
%---------------------------------------------------------------------- set
function Bdog = set(Bdog,varargin)
    for ii = 1:2:length(varargin)
        Bdog.(varargin{ii}) = varargin{ii+1};
    end
end
%---------------------------------------------------------------------- set
function Bdog = setInput(Bdog,x,varargin)
    
    if (numel(x) == 1) && ~isempty(varargin)
       X = Bdog.U*0;
       x(round(x)) = varargin{1};
    end

    Bdog.U = clamp(pdac(x(:) + Bdog.P0(:)),0,1);
end
%---------------------------------------------------------------------- set
function y = getOutput(Bdog)
   y = padc(Bdog.Log.y(end,:)); 
end
%---------------------------------------------------------------------- set
function resetInput(Bdog)
    Bdog.U = clamp(pdac(zeros(Bdog.NVeab*2,1) + Bdog.P0(:)),0,1);
    Bdog.tcpSendData(Bdog.U);
end

%--------------------------------------------------------- connect to board
function Bdog = connect(Bdog)

    if ~Bdog.autoConnect
        str = action({'(y)es, continue with SSH connection',...
            '(n)o, stop connection'},'s');
        
        if strcmp(str,'n'), return;
        elseif strcmp(str,'y'), pause(0.1);
        else, error('terminated');
        end
    end
    
    ip_   = Bdog.Ip;
    usr_  = Bdog.Usr;
    pwd_  = Bdog.Pwd;
    fprintf('* Starting SSH server \n');
    fprintf('* Establishing connection to target-pc: \n    dev -- ');
    cout('hyper',['',Bdog.Usr, '@' , Bdog.Ip,': ...']); pause(0.25);
    
    SSH_ = ssh2_config(ip_,usr_,pwd_);
    Bdog.SSH = ssh2(SSH_);
    
    cout('\b\b\b');
    cout('green','Connected!\n');
    Bdog.isConnected = true;

    %if checkexist(Bdog,Bdog.cppFile) && Bdog.autoCompile
    %    str = action({'(y)es, recompile executable file',...
    %            '(n)o'},'s');
    %     CallExecuted([FILENAME,' found!']);
    %     request = CallRequest('recompile executable file?','y/n');
    %     if strcmp(request,'y')
    %         brd = CommandShell(brd,'chmod 755 Soro*',0);
    %         CallExecuted(['compiled ',FILENAME,'!']);
    %     end
    %     if existsFile(brd,'soro.log')
    %         CallExecuted('cleared log file');
    %         brd = CommandShell(brd,'rm soro.log',0);
    %     end
    % else
    %     CallExecuted([FILENAME,' not found! File is required!'])
    % end
    
end
%--------------------------------------------------------------- disconnect   
function Bdog = disconnect(Bdog)
    fprintf('* Closing connection with target-pc: \n    dev -- '); 
    cout('hyper',['',Bdog.Usr, '@' , Bdog.Ip,': ...']); pause(0.25);
    
    Bdog.SSH = ssh2_close(Bdog.SSH);
    %Bdog.SSH = [];

    cout('\b\b\b');
    cout('key','Closed!\n');
end
%------------------------------ add sensor or actuator to device on i2c bus   
function Bdog = addDevice(Bdog,varargin)
   
    for ii = 1:3:length(varargin)
        
        i2cbus = dev2i2c(varargin{ii+1});
        ID = {vertcat(i2cbus{:}),(varargin{ii+2})};
        
        if Bdog.isConnected && Bdog.autoi2cCheck
            Bdog = Bdog.shell(['python3 i2cscan.py ',...
                num2str(varargin{ii+2})]);
            
            if ismember(Bdog.out,horzcat(i2cbus{:}))
               cout('key',['    dev -- ',varargin{ii+1},...
                   ' found on i2c bus: ',num2str(varargin{ii+2}),'\n']);
            else
               cout('err',['    dev -- ',varargin{ii+1}, ...
                   ' not found on i2c bus: ',num2str(varargin{ii+2}),...
                   '\n \t   Perhaps check wiring, or wrong i2c port?\n']); 
            end
        end
        
        Bdog.(varargin{ii}) = [Bdog.(varargin{ii});ID];
        Bdog.Ndev = Bdog.Ndev  + 1;
    end
end
%------------------------------ send terminal command (with printed output)    
function Bdog = command(Bdog,Arg)
    [Bdog.SSH, Bdog.out] = ssh2_command(Bdog.SSH, Arg, 1);
end
%---------------------------------------- send terminal command (no output) 
function Bdog = shell(Bdog,Arg)
    [Bdog.SSH, Bdog.out] = ssh2_command(Bdog.SSH, Arg, 0);
end
%------------------------------------------------------ transfer file to Pi   
function Bdog = transfer(Bdog,Arg)
    Bdog.SSH = scp_simple_put(Bdog.Ip,Bdog.Usr,Bdog.Pwd,Arg);
end
%--------------------------------------------- transfer python script to Pi   
function Bdog = transfertoPi(Bdog)
    
    N = 1;
    % while loop for matlab to finish file generations
    while ~exist(Bdog.pyFile,'file') && (N < 500)
       N = N + 1;
       pause(0.1);
    end
    
    Bdog.SSH = scp_simple_put(Bdog.ip,Bdog.usr,Bdog.pwd,Bdog.pyFile);
end
%-------------------------- Send data as a TCP client (formatted as double)
function Bdog = tcpSendData(Bdog, data)
    x = clamp(data,0,1);
    %fwrite(Bdog.TCPclient,clamp(x(:),0,1),"double");
    Bdog.TCPclient.write(clamp(x(:),0,1),"double");
end
%-------------------------- Recv data as a TCP client (formatted as double)
function data = tcpRecvData(Bdog,data_length)
    
    data = zeros(1,Bdog.NVeab*2);
    for ii = 1:(Bdog.NVeab*2)
        %data(ii) = fread(Bdog.TCPclient,data_length,"double");
        data(ii) = Bdog.TCPclient.read(data_length,"double");
    end
    
    %data1 = fread(Bdog.TCPclient,data_length,"double");
    %data2 = fread(Bdog.TCPclient,data_length,"double");
    
    %data = [data1,data2];
    Yn = padc(data(:).') - Bdog.P0(:).';
    Bdog.Log.y = vappend(Bdog.Log.y,Yn);
    
    Bdog.t = (1/Bdog.Frequency)*Bdog.LoopCounter;
    Bdog.Log.t = vappend(Bdog.Log.t,Bdog.t);
end
%----------------------------------------------------------- run executable 
function flag = loop(Bdog,Ts)
    
    if Bdog.LoopCounter == 0
        %fopen(Bdog.TCPclient);

        tic;
        fprintf('Waiting for VEABs to calibrate, please wait...\n');
        while toc<3
            tcpRecvData(Bdog,1);
            Bdog.tcpSendData(Bdog.U);
            pause(0.01);
        end

        tic;
        Bdog.LoopCounter = 1;
        Bdog.LoopIndex   = Ts/(1/Bdog.Frequency);
        tcpRecvData(Bdog,1);
        Bdog.tcpSendData(Bdog.U);
    else
        
        tcpRecvData(Bdog,1);
        Bdog.tcpSendData(Bdog.U);
        
%         while toc < 1/Bdog.Frequency
%             continue
%         end
%         
%         tic;
        Bdog.LoopCounter = Bdog.LoopCounter + 1;
        
    end
    
    if Bdog.LoopCounter <= Bdog.LoopIndex
        flag = true;
    else
        flag = false;
    end
    
end
%----------------------------------------------------------- run executable 
function reset(Bdog)
    tcpSendData(Bdog,0);
end
%--------------------------------------------- reads log file, return array
function A = read(Bdog,filename)
    if checkexist(Bdog,filename)
        Bdog = shell(Bdog,['tail ',Bdog.outFile,' -n 50000']);
        A = (cellfun(@str2num,Bdog.out(1:end),'un',0));
        A = vertcat(A{:});
    else
        A = [];
        cout('error','* log file does not exist!\n');
    end
end
%----------------------------------------------------- check if file exists  
function bool = checkexist(Bdog,name)
    Bdog = shell(Bdog,['[ -f ',name,' ] && echo "1" || echo "0"']);
    bool = str2double(Bdog.out{end});
end
end

methods (Access = private)

end
end
