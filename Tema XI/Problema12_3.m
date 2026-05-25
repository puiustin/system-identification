% Problema 12.3: Identificarea parametrilor paraziti si diagnoza defecte

clear all; close all; clc;
global FIG
FIG = 1;

% parametri de baza (Motor sanatos)
K0 = 4;
T0 = 0.5;
Tp0 = 0.05;  % constanta parazita (10 ori mai mică decât T0)
Tmax = 80;
Ts = 0.1;
U = 0.5;

% vector de intensitati zgomot pentru analiză (SNR descrescător)
lambda_vals = [0.05, 0.2, 0.5]; 

% rulare simulator diagnoză
% acesta va testa 3 tipuri de semnale de intrare (Square, Random, SPAB)
% si va afisa precizia identificarii celor 3 parametri (K, T, Tp).
ISLAB_12G(K0, T0, Tp0, Tmax, Ts, U, lambda_vals);

% nota privind Diagnoza (punctul d):
% in ISLAB_12G, s-a implementat un sistem cu 5 niveluri de uzura:
% 1. Incipient (Normal)
% 2. Mica
% 3. Moderata
% 4. Avansata
% 5. Severa (se recomanda oprirea motorului pentru reparatii)
% decizia se bazeaza pe cresterea Tp identificat fata de valoarea nominala.
