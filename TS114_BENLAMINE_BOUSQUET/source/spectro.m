function[Sx, f, t] = spectro(x, w, d,Nfft, Fs)
N = length(w);
[X, f, t] = stft(x, w, d, Nfft, Fs);
Sx = (abs(X).*abs(X))/N;


