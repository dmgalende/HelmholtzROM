% extract_single_freq(s,dt,fval);
%  This function obtains sample of the Fourier transform at a target frequency
%    s:   temporal signal
%    dt:  time discretization
%    fval: the frequency value where we want the sample


function [outval]=extract_single_freq(s,dt,fval);

nt=length(s);
if(mod(nt,2)==1),s=s(1:end-1);nt=length(s);end
nt_2 = ceil(nt/2);
fs = 1/dt;  
bin_vals = [0 : nt-1];
faxis = (bin_vals-nt_2)*fs/nt;
S=fftshift(fft(s));

myval_re=interp1(faxis,real(S),fval);
myval_im=interp1(faxis,imag(S),fval);
outval=myval_re+sqrt(-1)*myval_im;

end