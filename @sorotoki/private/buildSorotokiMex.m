function buildSorotokiMex(installPath)

global auto_approve

log = Log();
log.options.isDebug = true;
library = '/lib/';

if ~auto_approve
    disp('Build Sorotoki .mex executables?');
    reply = input(i18n('confirm'), 's');
    if isempty(reply)
        reply = i18n('confirm_yes');
    end

    if ~strcmpi(reply(1), i18n('confirm_yes'))
        disp(i18n('confirm_nvm'));
        return;
    end
end

log.info('Starting Sorotoki Mex Compiler');
log.info('Please be patient, this make take a while...');
log.hline()

%% buildMexLNH
log.debug('Nagivating towards libary directory');
cd([installPath, library]);

log.debug('Nagivating towards ..mex/localsNH/');
cd('sorotokifem/mex/localsNH/');

log.debug('Building mex localsNH');
fprintf('./sorotokifem/mex/localsNH build: ...');
flag = buildMexLNH;
fprintf(repmat('\b',1,60));
log.bool('./sorotokifem/mex/localsNH build:',flag,{'Build succesfull','Failed'});

%% buildMexLNH
log.debug('Nagivating towards libary directory');
cd([installPath, library]);

log.debug('Nagivating towards ..mex/localsYH/');
cd('sorotokifem/mex/localsYH/');

log.debug('Building mex localsYH');
fprintf('./sorotokifem/mex/localsYH build: ...');
flag = buildMexLYH;
fprintf(repmat('\b',1,60));
log.bool('./sorotokifem/mex/localsYH build:',flag,{'Build succesfull','Failed'});

%% buildCurveSweep
log.debug('Nagivating towards libary directory');
cd([installPath, library]);

log.debug('Nagivating towards ..mex/curvesweep/');
cd('sorotokimodel/mex/curvesweep/');

log.debug('Building mex buildCurveSweep');
fprintf('./sorotokimodel/mex/curvesweep build: ...');
flag = buildCurveSweep;
fprintf(repmat('\b',1,60));
log.bool('./sorotokimodel/mex/curvesweep build:',flag,{'Build succesfull','Failed'});

%% buildForwardKin
log.debug('Nagivating towards libary directory');
cd([installPath, library]);

log.debug('Nagivating towards ..mex/curvesweep/');
cd('sorotokimodel/mex/forwardkin/');

log.debug('Building mex buildForwardKin');
fprintf('./sorotokimodel/mex/forwardkin build: ...');
flag = buildForwardKin;
fprintf(repmat('\b',1,60));
log.bool('./sorotokimodel/mex/forwardkin build:',flag,{'Build succesfull','Failed'});

%% buildLagrangian
log.debug('Nagivating towards libary directory');
cd([installPath, library]);

log.debug('Nagivating towards ..mex/curvesweep/');
cd('sorotokimodel/mex/lagrangian/');

log.debug('Building mex buildLagrangian');
fprintf('./sorotokimodel/mex/lagrangian build: ...');
flag = buildLagrangian;
fprintf(repmat('\b',1,60));
log.bool('./sorotokimodel/mex/lagrangian build:',flag,{'Build succesfull','Failed'});

%% return to main folder
cd(installPath);
end


%% buildMexLNH
function flag = buildMexLNH
    flag = 0;
    try
    cfg = coder.config('mex');
    cfg.IntegrityChecks = false;
    cfg.ExtrinsicCalls = false;
    cfg.IntegrityChecks = false;
    cfg.SaturateOnIntegerOverflow = false;
    cfg.ResponsivenessChecks = false;
    cfg.NumberOfCpuThreads = 16;

    warning off
    codegen -config cfg LocalsNHFast -silent ...
        -args {coder.typeof(10,[inf,inf]), coder.typeof(10,[inf,inf]), coder.typeof(10,[1,1]), coder.typeof(10,[1,1]), coder.typeof(10,[1,1]), coder.typeof(10,[inf,inf]), coder.typeof(10,[inf,1,inf]), coder.typeof(10,[inf,inf,inf]), coder.typeof(10,[inf,1]), coder.typeof(10,[inf,inf]), coder.typeof(10,[inf,inf]), coder.typeof(10,[1,1]), coder.typeof(10,[1,1]), coder.typeof(10,[inf,1]), coder.typeof(10,[1,1]), coder.typeof(10,[1,1])}
        
    warning on    
    flag = 1;
    end
end

%% buildMexLYH
function flag = buildMexLYH
    flag = 0;
    try
    cfg = coder.config('mex');
    cfg.IntegrityChecks = false;
    cfg.ExtrinsicCalls = false;
    cfg.IntegrityChecks = false;
    cfg.SaturateOnIntegerOverflow = false;
    cfg.ResponsivenessChecks = false;
    cfg.NumberOfCpuThreads = 16;
    
    warning off
    codegen -config cfg LocalsYHFast -silent ...
        -args {coder.typeof(10,[inf,inf]), coder.typeof(10,[inf,inf]), coder.typeof(10,[1,1]), coder.typeof(10,[1,1]), coder.typeof(10,[1,1]), coder.typeof(10,[inf,inf]), coder.typeof(10,[inf,1,inf]), coder.typeof(10,[inf,inf,inf]), coder.typeof(10,[inf,1]), coder.typeof(10,[inf,inf]), coder.typeof(10,[inf,inf]), coder.typeof(10,[1,1]), coder.typeof(10,[1,1]), coder.typeof(10,[inf,1]), coder.typeof(10,[3,1]), coder.typeof(10,[3,1])}
        
    warning on    
    flag = 1;
    end
end

%% buildCurveSweep
function flag = buildCurveSweep
    flag = 0;
    try
    cfg = coder.config('mex');
    cfg.IntegrityChecks = false;
    cfg.ExtrinsicCalls = false;
    cfg.IntegrityChecks = false;
    cfg.SaturateOnIntegerOverflow = false;
    cfg.ResponsivenessChecks = false;
    cfg.NumberOfCpuThreads = 16;
    
    warning off
    codegen -config cfg curveSweepModifierFast -silent ...
        -args {coder.typeof(10,[inf,3]), coder.typeof(10,[inf,3]), coder.typeof(10,[inf,1]), coder.typeof(10,[4,4,inf])}
        
    warning on    
    flag = 1;
    end
end

%% buildForwardKin
function flag = buildForwardKin
    flag = 0;
    try
    cfg = coder.config('mex');
    cfg.IntegrityChecks = false;
    cfg.ExtrinsicCalls = false;
    cfg.IntegrityChecks = false;
    cfg.SaturateOnIntegerOverflow = false;
    cfg.ResponsivenessChecks = false;
    cfg.NumberOfCpuThreads = 16;
    
    codegen -config cfg computeForwardKinematicsGaussFast ...
        -args {coder.typeof(10,[inf,1]), coder.typeof(10,[inf,1]), coder.typeof(10,[3,1]), coder.typeof(10,[3,3]), coder.typeof(10,[6,1,inf,inf]), coder.typeof(10,[inf,inf,inf,inf]), coder.typeof(10,[6,inf]), coder.typeof(10,[inf,inf]), coder.typeof(10,[inf,inf])}
        
    warning on    
    flag = 1;
    end
end

%% buildLagrangian
function flag = buildLagrangian
    flag = 0;
    try
    cfg = coder.config('mex');
    cfg.IntegrityChecks = false;
    cfg.ExtrinsicCalls = false;
    cfg.IntegrityChecks = false;
    cfg.SaturateOnIntegerOverflow = false;
    cfg.ResponsivenessChecks = false;
    cfg.NumberOfCpuThreads = 16;
    
    codegen -config cfg computeLagrangianGaussFast ...
        -args {coder.typeof(10,[inf,1]), coder.typeof(10,[inf,1]), coder.typeof(10,[4,4,inf,inf]), coder.typeof(10,[6,inf,inf,inf]), coder.typeof(10,[6,inf,inf,inf]), coder.typeof(10,[inf,inf,inf,inf]), coder.typeof(10,[6,inf]), coder.typeof(10,[6,6,inf,inf]), coder.typeof(10,[6,6,inf,inf]), coder.typeof(10,[6,6,inf,inf]), coder.typeof(10,[3,1]), coder.typeof(10,[inf,inf]), coder.typeof(10,[inf,inf])}
    
    warning on    
    flag = 1;
    end
end