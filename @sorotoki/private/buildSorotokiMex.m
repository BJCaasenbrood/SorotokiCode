function buildSorotokiMex(installPath, base)

global auto_approve

log = Log();
log.options.isDebug = true;
library = base;

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
log.bool('./sorotokifem/mex/localsNH build:',flag,{'Build successful','Failed'});

%% buildMexLNH
log.debug('Nagivating towards libary directory');
cd([installPath, library]);

log.debug('Nagivating towards ..mex/localsYH/');
cd('sorotokifem/mex/localsYH/');

log.debug('Building mex localsYH');
fprintf('./sorotokifem/mex/localsYH build: ...');
flag = buildMexLYH;
fprintf(repmat('\b',1,60));
log.bool('./sorotokifem/mex/localsYH build:',flag,{'Build successful','Failed'});

%% buildCurveSweep
log.debug('Nagivating towards libary directory');
cd([installPath, library]);

log.debug('Nagivating towards ..mex/curvesweep/');
cd('sorotokimodel/mex/curvesweep/');

log.debug('Building mex buildCurveSweep');
fprintf('./sorotokimodel/mex/curvesweep build: ...');
flag = buildMexCSMF;
fprintf(repmat('\b',1,60));
log.bool('./sorotokimodel/mex/curvesweep build:',flag,{'Build successful','Failed'});

%% buildForwardKin
log.debug('Nagivating towards libary directory');
cd([installPath, library]);

log.debug('Nagivating towards ..mex/curvesweep/');
cd('sorotokimodel/mex/forwardkin/');

log.debug('Building mex buildForwardKin');
fprintf('./sorotokimodel/mex/forwardkin build: ...');
flag = buildMexFKFG;
fprintf(repmat('\b',1,60));
log.bool('./sorotokimodel/mex/forwardkin build:',flag,{'Build successful','Failed'});

%% buildLagrangian
log.debug('Nagivating towards libary directory');
cd([installPath, library]);

log.debug('Nagivating towards ..mex/curvesweep/');
cd('sorotokimodel/mex/lagrangian/');

log.debug('Building mex buildLagrangian');
fprintf('./sorotokimodel/mex/lagrangian build: ...');
flag = buildMexCLGF;
fprintf(repmat('\b',1,60));
log.bool('./sorotokimodel/mex/lagrangian build:',flag,{'Build successful','Failed'});

%% return to main folder
cd(installPath);
end