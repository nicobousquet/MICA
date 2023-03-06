clear;
close all;

%Partie spectrogramme ----------------------------------------------------------------------------------------

load ../../../MICA_project/Sujet/data/ecg_PVC.mat;
ecg = ecg;

Nfft = 4000;
N = 1600; %largeur fenêtre
d = 800; %oveerlap=N-d

[X, f, t] = stft(ecg, hamming(N), d, Nfft, Fs);

[Sx, f, t] = spectro(ecg, hamming(N), d, Nfft, Fs);

figure,
imagesc(t, f, 10*log10((Sx)));
h = colorbar;
h.Label.String = "Power/Frequency (dB/Hz)";
caxis([-100 60]);
xlabel("Time(s)");
ylabel("Frequency (Hz)")

%Partie Pathologies ------------------------------------------------------------------------------------------

Ts=1/Fs;
B_derivative =[-1 -2 0 2 1];
A_derivative= [8*Ts];
B_low = [1 0 0 0 0 0 -1];
A_low = [1 -1];
B_high = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
A_high = [1 -1];
delay_filter = dfilt.delay(2);

%Pan and Thompkins algorithm

signal_low_pass=filter(B_low,A_low,filter(B_low,A_low,ecg));
signal_band_pass=filter(B_high,A_high,signal_low_pass);
signal_band_pass = filter(delay_filter, signal_band_pass);
signal_derivative=filter(B_derivative, A_derivative, signal_band_pass);
signal_squared = signal_derivative.^2;

%window integration

N = 100;
integral = zeros(size(signal_squared));
for k = 1:1:length(signal_squared)
    for i = 1:1:N
        if k-i>0
            integral(k) = integral(k)+ 1/N * signal_squared(k-i);
        end
    end
end

coeff_dir = []; 
for k=1:1:length(integral)-1
    coeff_dir(k)=(integral(k+1)-integral(k));
end

X_rwave_tresh = [1];
X=[];
bool = 0;
for k = 1:1:length(coeff_dir)
    if(k+1000<length(coeff_dir))
        if(coeff_dir(k)>max(coeff_dir(k:k+1000)/10))  %threshold = max(coeff_dir(k:k+1000)/10))
            if(not(bool))
                X_rwave_tresh = [X_rwave_tresh k];
            end
            bool = 1; 
        elseif (coeff_dir(k)<min(coeff_dir(k:k+1000)/10)) %threshold = min(coeff_dir(k:k+1000)/10))
            bool = 0;
        end
    end
end


X_rwave = [];
Y_rwave = [];
for k = 1:1:length(X_rwave_tresh)-1
   [y,x] = max(   ecg(X_rwave_tresh(k):X_rwave_tresh(k+1))   );
   x = x + X_rwave_tresh(k);
   X_rwave = [X_rwave x];
   Y_rwave = [Y_rwave y];
end
t_rwave = X_rwave/Fs;

bpm = [];
abcisse_bpm = t_rwave(1:end-10);
for k = 1:1:length(X_rwave)-10
    delta = X_rwave(k+10)-X_rwave(k); %on prend la moyenne de tous les 10 battements pour avoir des changements de fréquence cardiaque un peu plus cohérents
    bpm = [bpm 10*60*Fs/delta]; 
end

figure, %figure du nombre de battements par minute
plot(abcisse_bpm, bpm);
xlabel("Time(s)");
ylabel("Heartbeat (beat/min)");
title("Heart rate versus time");

mean_bpm = ((length(X_rwave)-1)*60)/(length(ecg)*Ts); %nombre de bpm moyen par seconde

%DETECTION SWAVE

X_swave = [];
Y_swave=[];
for i=1:1:length(X_rwave)-1
    [yS,xS] = min(ecg(X_rwave(i):X_rwave(i)+floor(0.25*(X_rwave(i+1)-X_rwave(i)))));
    x=xS+X_rwave(i);
    X_swave = [X_swave x];
    Y_swave = [Y_swave yS];
end

%DETECTION QWAVE

X_qwave = [];
Y_qwave=[];
for i=2:1:length(X_rwave)
    [yQ,xQ] = min(ecg(X_rwave(i)- floor(0.1*(X_rwave(i)-X_rwave(i-1))):X_rwave(i)));
    x=xQ+(X_rwave(i)-floor(0.1*(X_rwave(i)-X_rwave(i-1))));
    X_qwave = [X_qwave x];
    Y_qwave = [Y_qwave yQ];
end

%DETECTION PWAVE

x=0;
X_pwave = [];
Y_pwave=[];
for i=2:1:length(X_qwave)
    [yp,xp] = max(ecg(X_qwave(i)- floor(0.15*(X_qwave(i)-X_qwave(i-1))):X_qwave(i)-1));
    x=xp+(X_qwave(i)-floor(0.15*(X_qwave(i)-X_qwave(i-1))));
    X_pwave = [X_pwave x];
    Y_pwave = [Y_pwave yp];
end

%DETECTION TWAVE

x=0;
X_twave = [];
Y_twave=[];
for i=1:1:length(X_swave)-1
    [yt,xt] = max(ecg(X_swave(i)+1:X_swave(i) + floor(0.70*(X_rwave(i+1)-X_rwave(i))-(X_swave(i)-X_rwave(i)))));
    x=xt+X_swave(i);
    X_twave = [X_twave x];
    Y_twave = [Y_twave yt];
end

%BPM

bpm = [];
abcisse_bpm = t_rwave(1:end-1);
for k = 1:1:length(X_rwave)-1
    delta = X_rwave(k+1)-X_rwave(k);
    bpm = [bpm 60*Fs/delta];
end

%ECTOPIC BEAT
delta = [];
for k = 1:1:length(X_rwave)-1
    delta(k) = X_rwave(k+1)-X_rwave(k);
end

X_ectopic_begin = [];
X_ectopic_finish = [];
for k=1:1:length(delta)-1
    if delta(k+1)>=2*delta(k)
        X_ectopic_begin = [X_ectopic_begin  X_rwave(k)];
        X_ectopic_finish = [X_ectopic_finish X_rwave(k+2)];
    end
end

X_ectopic = [X_ectopic_begin' X_ectopic_finish']; 
X_ectopic = X_ectopic.*Ts; %la matrice X_ectopic montre en chaque ligne le début et la fin d'un battement ectopique

%Atrial Fibrilation

N = length(delta);
mean_delta = 60/mean_bpm;
gamma = zeros(1, N);
delta=delta.*Ts;

for k=1:1:N
    for n=1:1:N-k
        gamma(k) = gamma(k) + (1/(N-k-1))*(delta(n+k)-mean_delta)*(delta(n)-mean_delta);
    end
end

%PLOT
figure, %figure de l'ecg avec PQRST
plot(Ts:Ts:length(ecg)*Ts, ecg);
hold on;
plot(X_rwave.*Ts, Y_rwave,"+r");
text(X_rwave.*Ts, Y_rwave," R",'Color','red');
plot(X_swave.*Ts, Y_swave,"+b");
text(X_swave.*Ts, Y_swave," S",'Color','blue');
plot(X_qwave.*Ts, Y_qwave,"+r");
text(X_qwave.*Ts, Y_qwave," Q",'Color','red');
plot(X_pwave.*Ts, Y_pwave,"+r");
text(X_pwave.*Ts, Y_pwave," P",'Color','red');
plot(X_twave.*Ts, Y_twave, "+b");
text(X_twave.*Ts, Y_twave," T",'Color','blue');
xlabel("Time (s)");
ylabel("Cardiac intensity (mV)");
title("PQRST waves detection");

figure, %figure de l'autocovariance pour l'atrial fibrillation
plot((1:1:N)*mean_delta, gamma);
xlabel("Time(s)");
ylabel("Autocovariance");
title("Autocovariance of an electrocardiogram");



