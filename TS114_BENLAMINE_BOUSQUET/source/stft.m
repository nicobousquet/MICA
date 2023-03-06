function[X, f, t] = stft(x, w, d, Nfft, Fs)
M = floor(length(x)/d); %nombre de colonnes de la matrice
N = length(w); %longueur de la fenêtre
X = zeros(N, M+1); %Matrice X
t = 0:N/Fs:(M+1)*(N/Fs) - (N/Fs); %temps au début de chaque fenêtre
f = 0:Fs/Nfft:Fs-(Fs/Nfft); %fréquences en Hz

for i=1:1:M+1
    if d*(i-1)+N < length(x) %on ne dépasse pas la longueur du vecteur x
        X(:, i) = x(1+d*(i-1):d*(i-1)+N)';
        X(:, i) = X(:, i).*w;  
    else %pour la dernière colonne, problème car fin de x plus petit que colonne de X, donc problème de dimension
        indice = 1;
        for j=d*(i-1)+1:1:length(x)
            X(indice, i) = x(j); %on décide de remplir la dernière colonne de X terme à terme
            indice = indice+1;
        end
        X(:, i) = X(:, i).*w; %on multiplie les dernières colonnes par la fenêtre
    end
end
X = fft(X,Nfft, 1); %on fait la TF de chaque colonne



