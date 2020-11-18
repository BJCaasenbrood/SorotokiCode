classdef Bdog < handle

    properties (Access = public)
        ip;
        usr;
        pwd;
        prt;
        out;
    end
    
    properties (Access = private)
        cppFile;
        outFile;
        conFile;
        autoConnect;
        autoCompile;
        SSH;
    end
    
%--------------------------------------------------------------------------
methods  
%---------------------------------------------------------------------- get   
function obj = Bdog(Ip,User,Password,varargin)
    obj.prt = '22';
    obj.ip  = Ip;
    obj.usr = User;
    obj.pwd = Password;
    
    obj.cppFile = 'main.cpp';
    obj.outFile = 'soro.log';
    obj.conFile = 'config.txt';
    obj.autoCompile = false;
    obj.autoConnect = false;
    
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
    
    if obj.autoConnect, obj = connect(obj); end
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

%------------------------------------------------------------ CONNECT BOARD
function Bdog = connect(Bdog)

    if ~Bdog.autoConnect
        str = action({'(y)es, continue with SSH connection',...
            '(n)o, stop connection'},'s');
        
        if strcmp(str,'n'), return;
        elseif strcmp(str,'y'), pause(0.1);
        else, error('terminated');
        end
    end
    
    fprintf('* establishing board connection... \n');
    
    Bdog.SSH = connectBoard(Bdog);

    if checkexist(Bdog,Bdog.cppFile) && Bdog.autoCompile
        str = action({'(y)es, recompile executable file',...
                '(n)o'},'s');
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
    end
end

%---------------------------------------------------------------------- get   
function Bdog = disconnect(Bdog)
    fprintf('* disconnecting...\n ');
    Bdog.SSH = ssh2_close(Bdog.SSH);
    Bdog.SSH = [];
    cout('green','* boards disconnected!\n');
end

%----------------------------------------------------------------- command    
function Bdog = command(Bdog,Arg)
    [Bdog.SSH, Bdog.out] = ssh2_command(Bdog.SSH, Arg, 1);
end

%----------------------------------------------------------------- command    
function Bdog = shell(Bdog,Arg)
    [Bdog.SSH, Bdog.out] = ssh2_command(Bdog.SSH, Arg, 0);
end

%---------------------------------------------------------------------- get  
function bool = checkexist(Bdog,name)
    Bdog = shell(Bdog,['[ -f ',name,' ] && echo "1" || echo "0"']);
    bool = str2double(Bdog.out{end});
end

%---------------------------------------------------------------------- get  
function A = read(Bdog,filename)
    if checkexist(Bdog,filename)
        Bdog = shell(Bdog,'tail soro.log -n 50000');
        A = (cellfun(@str2num,Bdog.out(1:end),'un',0));
        A = vertcat(A{:});
    else
        A = [];
        cout('error','* log file does not exist!\n');
    end
end

%---------------------------------------------------------------------- get   
function Bdog = run(Bdog)
    if checkexist(Bdog,'soro.log')
        Bdog = shell(Bdog,'rm soro.log');
    end
    pause(0.1);
    %Bdog = command(Bdog,['sudo nohup ./main config.txt >soro.log',...
    %    '</dev/null 2>/dev/null& && \n']);
    cprintf('error', '>> running executable on real-time target... \n');
    efile = ['./',erase(Bdog.cppFile,'.cpp')];
    %CMD = ['sudo nohup ',efile,' ',Bdog.conFile, ' >',Bdog.outFile];
    CMD = ['sudo ',efile,' ',Bdog.conFile];

    Bdog = command(Bdog,CMD);
    cprintf('green', '>> done!\n');
    cout(['* log file - ', Bdog.outFile, '\n']);
end

end
end
