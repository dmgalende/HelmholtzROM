function Y = getSignalAmplituteFreq(fval)

dt = 4e-3;
nt = 1001;
wvlt_trace=ricker_wavelet(dt,nt,0.1,10);

Y = zeros(length(fval),1);
for i = 1:length(fval)
    freq = fval(i);
    Y(i) = extract_single_freq(wvlt_trace,dt,freq);
end

% %Final frequency              
% Fs = 1/dt;
% 
% %Wavelet
% L = length(y);
% 
% %Fast fourier tranform (spectra) and interpolation
% NFFT = 2^nextpow2(L);   %Next power of 2 from length of y
% Y = fft(y,NFFT)/L;
% f = Fs/2*linspace(0,1,NFFT/2+1);
% Y = Y(1:NFFT/2+1);
% Ys = interp1(f,Y,f0,'linear');