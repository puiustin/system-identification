clc; clear;

% configurare parametri simulare
N = 100;    % numarul de esantioane (puncte de date)
Ts = 0.1;   % perioada de esantionare (secunde)
t = (0:N-1)' * Ts;  % vectorul de timp generat pe coloana

% generarea datelor de iesire (y)
% prima coloana contine timpul (t), urmatoarele 3 coloane sunt date random
y = [t, randn(N,4)];    % avem in total 3 canale de iesire

% generarea datelor de intrare (u)
% generam o matrice cu  coloane (2 canale de intrare)
u = randn(N,3);      

% apelam functia creata pentru a genera obiectul iddata mimo
% y contine timpul si iesirile, u contine intrarile
DATA = make_IDDATA_MIMO(y,u);

% afisare rezultat in consola pentru verificare
disp('test mimo cu 3 iesiri si 2 intrari:');
disp(DATA);
