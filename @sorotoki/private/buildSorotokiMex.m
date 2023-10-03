function buildSorotokiMex(installPath)

    log = Log();
    log.options.isDebug = true;
    library = '/src/';
    path = installPath(1:end-10);

    disp('Build Sorotoki .mex executables?');
    reply = input(i18n('confirm'), 's');
    if isempty(reply)
        reply = i18n('confirm_yes');
    end

    if ~strcmpi(reply(1), i18n('confirm_yes'))
        disp(i18n('confirm_nvm'));
        return;
    end

log.info('Starting Sorotoki Mex Compiler');
log.info('Please be patient, this make take a while...');
log.hline()

%% buildMexLNH
log.debug('Nagivating towards libary directory');
cd([path, library]);

log.debug('Nagivating towards ..mex/localsNH/');
cd('SorotokiFem/mex/localsNH/');

log.debug('Building mex localsNH');
fprintf('./SorotokiFem/mex/localsNH build: ...');
flag = buildMexLNH;
fprintf(repmat('\b',1,60));
log.bool('./SorotokiFem/mex/localsNH build:',flag,{'Build succesfull','Failed'});

%% buildMexLNH
log.debug('Nagivating towards libary directory');
cd([path, library]);

log.debug('Nagivating towards ..mex/localsYH/');
cd('SorotokiFem/mex/localsYH/');

log.debug('Building mex localsYH');
fprintf('./SorotokiFem/mex/localsYH build: ...');
flag = buildMexLYH;
fprintf(repmat('\b',1,60));
log.bool('./SorotokiFem/mex/localsYH build:',flag,{'Build succesfull','Failed'});

end


%%
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