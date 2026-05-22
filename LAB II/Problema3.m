clc; 
clear; 
close all; 
SNR =  0.1;
ISLAB_2B([],[],SNR);

%Problema 2.3
%ISLAB_2B(parte reala poli, parte imaginara, SNR);


%Pentru SNR =3 si Parte reala poli = 0.5 si parte imaginara = 0.5 ->
%observam ca avem o pereche de poli conjugati..si ca cele doua modele nu 
%au denistatea spectrala egala -> deci cele doua nu sunt echivalente inca

%SNR = 100 -> se observa ca polii raman in aceleasi locuri, dar zerourile
%se indreapta catre centru, densitatile spectrale par sa tinda si astfel
%cele doua modele sunt echivalente
 
%SNR = 1 -> zerourile se indreapta catre valoarea polilor, iar densitatile
%spectrale ale modelelor tind spre cea a zgomotului alb

%SNR = 0.1 -> zerourile si polii se suprapun (aproape)


%Concluzie...SNR influenteaza doar pozitia zerourilor..SNR mare -> pozitia
%zerourilor se duce spre centru cercului unitate...SNR mic -> zerourile se
%suprapun peste poli

%2 %SNR = 100 -> zerourile aroape se suprapun in 0, iar cele doua densitati
%spectrale sunt aproape identice..cu cat SNR creste cu atat influenta
%zgomotului este mai redusa..si astfel cele doua modele coincid

%SNR = 10000 -> aceiasi situatie, valoarea densitatii spectrale a
%zgomotului este mult sub cea a modelului.
%SNR = 0.01 -> zgomotul foarte puternic -> polii si zerourile coincid ->
%densitatea spectrala o  urmareste in totalitate pe cea a  zgomotului

%Concluzie: SNR mare -> zgmotul nu mai influenteaza densitatea spectralaa
% SNR mic -> denistatea spectrala a modelului coincide cu cea a zgomotului,
% astfel ca zgomotul influenteaza in totalitate raspunsul.
