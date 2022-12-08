function [phase, filtered_sig] = hilberttransform(signal,Fs,bp)
% de-trend the signal with a band-pass filter
% you may need the signal processing toolbox ....
bp = bp * 2 / Fs; % convert Hz to radians/S
[N, Wn] = buttord( bp, bp .* [.5 1.5], 3, 20);
[B,A] = butter(N,Wn);
filtered_sig = filtfilt(B,A,signal); % zero-phase filtering

% remove negative frequency component of fourier transform
X = fft(filtered_sig);
halfway = 1 + ceil(length(X)/2);
X(halfway:end) = 0;
ht_signal = ifft(X);

% keep phase
phase = angle( ht_signal );
end

