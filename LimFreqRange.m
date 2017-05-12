function fidsLimFreq = LimFreqRange(MRS_struct)

fids = MRS_struct.fids.data;
N = size(fids,1);
freqrange = MRS_struct.p.sw/MRS_struct.p.LarmorFreq;
freq = (N+1-(1:1:N))/N*freqrange+4.68-freqrange/2.0;
fidsFT = fftshift(fft(fids,[],1),1);

metab_range = freq >= 1.9 & freq <= 3.5;
fidsFTLimFreq = fidsFT(metab_range,:);
% fidsLimFreq = fftshift(ifft(fidsFTLimFreq,[],1));
if mod(size(fidsFTLimFreq,1),1)==0
    %disp('Length of vector is even.  Doing normal conversion');
    fidsLimFreq = fftshift(ifft(fidsFTLimFreq,[],1));
else
    %disp('Length of vector is odd.  Doing circshift by 1');
    fidsLimFreq = ifft(circshift(fftshift(fidsFTLimFreq)),[],1);
end