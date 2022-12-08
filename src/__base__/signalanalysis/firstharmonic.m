function [F,P,A] = firstharmonic(X,t,varargin)
%Fs = 1000;            % Sampling frequency                    
%T = 1/Fs;             % Sampling period       
%L = 1500;             % Length of signal
%t = (0:L-1)*T;        % Time vector

for ii = 1:2:length(varargin)
    if strcmp(varargin{ii},'window')
        W = varargin{ii+1};
        [~,id1] = min(abs(t - W(1)));
        [~,id2] = min(abs(t - W(2)));
        t = t(id1:id2);
    end
end

dt = mean(diff(t));
Fs = 1/dt;
L  = numel(t);

Y = fft(X(:));
P2 = abs(Y/L);
Z2 = angle(Y/L);

P1 = P2(1:(L/2+1));
P1(2:end-1) = 2*P1(2:end-1);

Z1 = Z2(1:(L/2+1));

% Z1 = Z2(1:L/2+1);
% Z1(2:end-1) = 2*Z1(2:end-1);

f = Fs*(0:(L/2))/L;
%[~,I] = maxk(P1,10);
Ys = P1; Ys(1) = 0;
[A,loc] = findresonance(Ys,5);

F = interp1(0:(L/2),f,loc);
P = interp1(0:(L/2),Z1,loc);
end

