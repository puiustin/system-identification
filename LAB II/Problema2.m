
clc;
close all;
N = 10000; % orizont de masurare
tau_max = 100; %pivoti autocovarianta
nr = 1; % numar de realizari

c = 0.4;
a = - 0.95;
C = [1 c]; %polinom asociat zgomotului e[n]
A = [1 a]; %polinom asociat iesirii y[n]
ISLAB_2A(C,A,N,tau_max,nr);
%Problema 2.2
%Ca in oricare sistem discret trebuie impusa conditia de stailitate, astfel
%valoarea polului trebuie sa aiba |a| < 1, adica sa fie in cadrul
%cercului unitate.
% a) Analiza functiilor de covarianta: variaza N si tau_max in
%functie de diferite locatii ale polilor:

% Pentru cazul in care avem a = 0.4,  N = 100, tau_max = 90
%cu cat ne apropiem mai tare cu tau de N...covarianta devine instabila..
%folosind Ipoteza ergodica...aproxima autocovarianta cu media, iar fiindca
%tau devine mare -> aceasta medie are putini termeni, deci are valori mari.
% valoarea polului este una stabila (si nu este la limita stabilitatii)
% atunci autocovarianta nu are un caracter oscilant
%Pentru cazul in care avem  a = 0.4, N = 150, tau_max = 50
%grafiul estimat se aproprie cat mai mult de cel simulat, iar conform IE
% cu cat am mari mai mult N (orizontul de masura), cu atat are trebui ca
%eroarea dintre cele doua grafice sa fie din ce in ce mai mica...voi atasa 
%poze cu graficele pentru un pivot fix si N variabil: N = 200, N = 500,
% N = 1000.
% Despre alegerea lui tau se vorbeste ca ar trebui sa fie maxim pana in N/4
% (Anexa A)


% alegem acuma un pol care sa fie aproape de instabilitate (iese din cercul
% unitate) a = 0.95...avand in vedere limita de stabilitate pe care o
% provoaca amplasarea polului aproape in exterior...autocovarianta capata
% un caracter oscilant...!! care incepe sa dispara cu cat tau este mai
% mare, astfel ca autocovrarianta se duce spre 0, cu cate pivotul dintre
% elemente este mai mare.
% dar aceleasi remarci sunt valide...N mare si tau
% destul de mic (raportat la N) produc ca graficul estimat si cel real 
% (~) coincid.

%Pentru   a = -0.95 se observa o diferenta a faptului ca primul element al
%autocovariantei nu mai este 1 ci este mult mai mare, iar scaderea este
%nu mai este la fel de brusca ca atunci cand aveam pol negativ... Se mai
%observa ca scaderea spre 0 devine mai lenta cu cat modulul lui a devine
%mai mare |a| ~= 1.

%%

clc;
close all;
clear;
N = 10000; % orizont de masurare
tau_max = 100; %pivoti autocovarianta
nr = 1; % numar de realizari

%AR[1] si MA[1]...dupa cum zice IE -> cu cat mai mare orizontul de masura
%cu atat valoarea este mai precisa, mai exact, mai aproape de
%realitate...astfel ca alegerea unui N foarte mare aduce o eroare foarte
%mica intre cele doua grafice (cel estimat vs cel real).
c = 0.5;
a =  0.5;
C = [1 c]; %polinom asociat zgomotului e[n]
%A = [1 a]; %polinom asociat iesirii y[n]
ISLAB_2A(C,1,N,tau_max,nr);
