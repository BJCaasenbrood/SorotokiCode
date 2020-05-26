classdef Process

    properties (Access = public)
    Name;
    Batch;
    BatchData;
    end
    
    properties (Access = private)
    Time0;
    cmdlength;
    end
    
%--------------------------------------------------------------------------
methods  
    
%------------------------------------------------------------ Process Class
function obj = Process(str,varargin) 
     if nargin < 1, str = 'Process'; end
     obj.Name = str;
     obj.Time0 = clock;
     obj.cmdlength = 80;
     obj.Batch= {};
     obj.BatchData = {};
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
end

%---------------------------------------------------------------------- get     
function varargout = get(Process,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = Model.(varargin{ii});
        end
    else
        varargout = Process.(varargin);
    end
end
        
%---------------------------------------------------------------------- set
function Process = set(Process,varargin)
    for ii = 1:2:length(varargin)
        Process.(varargin{ii}) = varargin{ii+1};
    end
end

%---------------------------------------------------------------------- set
function Process = add(Process,varargin)
    for ii = 1:length(varargin)
        Process.Batch{end+1,1} = {varargin{ii}};
        Process.Batch{length(Process.Batch),2} = -1;
    end
    t0 = vertcat(Process.Batch{:,1});
    t1 = vertcat(Process.Batch{:,2});
    
    Process.Batch = {t0,t1};
end

%------------------------------------------------------------- time request
function t = time(Process,Request)
    
    if nargin < 2, Request = false; end
    
    switch(Request)
    case('now');    t = clock;
    case('active'); t = etime(clock,Process.Time0);
    case('init');   t = Process.Time0;
    otherwise;      t = etime(clock,Process.Time0);
    end
    
end

function line(Process)
    n = Process.cmdlength;
    fprintf('+');
    for ii = 2:n-1, fprintf('-'); end
    fprintf('+\n');
end

end
end