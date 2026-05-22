clc;
clear;
close all;

%% ARX[2,2] - filtrata 4D;
[a4D,b4D,lambda4D,magi4D,phii4D,mag4D,phi4D,f4D] = ISLAB_4D();

%% ARX[2,2] - nefiltrata 4c;

[a4c,b4c,lambda4c,magi4c,phii4c,mag4c,phi4c,f4c] = ISLAB_4C();

%% ARX[1,1] - filtrata 4b
[a4b,b4b,lambda4b,magi4b,phii4b,mag4b,phi4b,f4b] = ISLAB_4B();


%% ARX[1,1] - nefiltrata 4a 
[a4a,b4a,lambda4a,magi4a,phii4a,mag4a,phi4a,f4a] = ISLAB_4A();
