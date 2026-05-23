function [Mid,Did,Dva] = ISLAB_6C_bj(B,C,F,D,nk,N,sigma,lambda) 
%
% ISLAB_6C_BJ   Modul care estimeaza un model Box-Jenkins folosind MCMMPE.
%               (Adaptare dupa ISLAB_6B pentru Problema 6.3 - BJ Extended)
%
% Inputs:	B, C, F, D - Polinoamele reale ale procesului BJ
%           nk, N, sigma, lambda - Parametrii simularii
%
% Outputs:  Mid - Modelul identificat (obiect IDPOLY)
%           Did, Dva - Datele generate
%
% Explanation:	Utilizeaza metoda MCMMPE (implementata in bj_e) 
%               pentru a identifica un proces Box-Jenkins.
%

%
% BEGIN
% 
global FIG ;			% Figure number handler 
FIG = 1;                 

%
% Constants
% ~~~~~~~~~
alpha = 3 ;			% Weighting factor of confidence disks radius. 
pf = 1 ; 			% Plot flag: 0=no, 1=yes. (Modificat pe 1 pt vizualizare)
Nb = 3 ;			% Maximum index of X part. 
Nc = 3 ;			% Maximum index of MA part. 
Nf = 3 ;			% Maximum index of Input Denom part.
Nd = 3 ;			% Maximum index of Noise Denom part. 
Ts = 1 ; 			% Sampling period. 

% 
% Faults preventing & Input Validation (Identic cu 6B)
% ~~~~~~~~~~~~~~~~~
if (nargin < 8), lambda = 1 ; end 
if (isempty(lambda)), lambda = 1 ; end 
lambda = abs(lambda(1)) ; 
if (~lambda), lambda = 1 ; end 

if (nargin < 7), sigma = 1 ; end 
if (isempty(sigma)), sigma = 1 ; end 
sigma = abs(sigma(1)) ; 

if (nargin < 6), N = 250 ; end
if (isempty(N)), N = 250 ; end 
N = abs(fix(N(1))) ; 
if (~N), N = 250 ; end 

if (nargin < 5), nk = 1 ; end
if (isempty(nk)), nk = 1 ; end
nk = abs(fix(nk(1))) ; 
if (~nk), nk = 1 ; end 

if (nargin < 4), C = [1 -1 0.2] ; end
if (isempty(C)), C = [1 -1 0.2] ; end 

if (nargin < 3), B = [1 0.5] ; end
if (isempty(B)), B = [1 0.5] ; end 

if (nargin < 2), F = [1 -1.5 0.7] ; end
if (isempty(F)), F = [1 -1.5 0.7] ; end

if (nargin < 1), D = [1 1.5 0.7] ; end
if (isempty(D)), D = [1 1.5 0.7]; end

% Stabilizare polinoame numitor
F = roots(F) ; 
F(abs(F)>=1) = 1./F(abs(F)>=1) ; 
F = poly(F) ;
D = roots(D) ; 
D(abs(D)>=1) = 1./D(abs(D)>=1) ; 
D = poly(D) ; 

% 
% Constructing the data provider process 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
P = idpoly(1,[zeros(1,nk) B],C,D,F,lambda,Ts) ; 

% 
% Generating the identification data 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Did = gen_data(P,N,sigma,lambda) ; 
Dva = gen_data(P,N,sigma,lambda) ; 

% 
% Estimating all models via MCMMPE (Utilizand bj_e)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Nb = Nb+1 ; Nc = Nc+1 ; Nf = Nf+1; Nd = Nd+1;

M = cell(Nb,Nc,Nd,Nf) ; 			
Lambda = 1000*lambda*ones(Nb,Nc,Nd,Nf) ; 	
Yid = cell(Nb,Nc,Nd) ; 			
Yva = cell(Nb,Nc,Nd,Nf) ; 			
PEid = cell(Nb,Nc,Nd,Nf) ; 		
PEva = cell(Nb,Nc,Nd,Nf) ; 		
Eid = zeros(Nb,Nc,Nd,Nf) ; 		
Eva = zeros(Nb,Nc,Nd,Nf) ; 		
yNid = Did.y-mean(Did.y) ;		
yNva = Dva.y-mean(Dva.y) ;		
Viid = zeros(Nb,Nc,Nd,Nf) ; 		
Viva = zeros(Nb,Nc,Nd,Nf) ; 		
optidx = zeros(3,4) ; 			

if (~pf)
   war_err([blanks(4) '* MCMMPE (BJ) Models estimation started. ' ... 
                      'This may take few moments. Please wait...']) ; 
end 

for nf=2:Nf
  if (~pf), disp([blanks(10) 'nf = ' int2str(nf-1)]) ; end 
  for nd=2:Nd
    if (~pf), disp([blanks(16) 'nd = ' int2str(nd-1)]) ; end 
    for nb=2:Nb
      if (~pf), disp([blanks(22) 'nb = ' int2str(nb-1)]) ; end 
      for nc=2:Nc
          if (~pf), disp([blanks(22) 'nc = ' int2str(nc-1)]) ; end 
      
          % Conditia de estimare
          if ((nb>1) | (nc>1) | (nd>1) | (nf>1))
             
             % >>>>>>>>>>>> ZONA MODIFICATA: APEL bj_e <<<<<<<<<<<<<<
             % Construim vectorul de structura: [nb nc nd nf nk]
             % Scadem 1 pentru ca bucla incepe de la 2 (index Matlab) dar ordinul e de la 1
             si_curr = [nb-1, nc-1, nd-1, nf-1, nk];
             
             Mid = bj_e(Did, si_curr); 
             % >>>>>>>>>>>>>>> SFARSIT MODIFICARE <<<<<<<<<<<<<<<<<<<
     
             M{nb,nc,nd,nf} = Mid ; 		
             Lambda(nb,nc,nd,nf) = Mid.NoiseVariance ; 
             
             % Calcule performanta
             PEid{nb,nc,nd,nf} = resid(Mid,Did) ;
             PEva{nb,nc,nd,nf} = resid(Mid,Dva) ; 
             ys = compare(Mid,Did) ;	
             Yid{nb,nc,nd,nf} = ys.y ;	% Modificat .OutputData -> .y pentru siguranta
             ys = compare(Mid,Dva) ;        
             Yva{nb,nc,nd,nf} = ys.y ;    
    					
             Eid(nb,nc,nd,nf) = 100*(1-norm(PEid{nb,nc,nd,nf}.y)/norm(yNid)) ; 
             Eva(nb,nc,nd,nf) = 100*(1-norm(PEva{nb,nc,nd,nf}.y)/norm(yNva)) ; 
             
             Viid(nb,nc,nd,nf) = valid_LS(Mid,Did) ; 
             Viva(nb,nc,nd,nf) = valid_LS(Mid,Dva) ; 
             
             % Selectia indicilor optimi (Testele F)
             if ((nb>1) & (nc>1) & (nd>1) & (nf>1))	
                if (~sum(optidx(1,:)))	% F-test prediction error. 
                   ys = [Lambda(nb-1,nc,nd,nf) Lambda(nb,nc-1,nd,nf) ...
                         Lambda(nb,nc,nd-1,nf) Lambda(nb,nc,nd,nf-1)] ; 
                   ys = ys/Lambda(nb,nc,nd,nf) - 1 ; 
                   if (sum(ys<(4/N))==3) | ((nb==Nb) & (nc==Nc) & (nd==Nd) & (nf==Nf))
                      optidx(1,:) = [nb nc nd nf] ;
                   end 
                end 
                if (~sum(optidx(2,:)))	% F-test fitness (ID). 
                   ys = [Eid(nb-1,nc,nd,nf) Eid(nb,nc-1,nd,nf) ...
                         Eid(nb,nc,nd-1,nf) Eid(nb,nc,nd,nf-1)]; 
                   ys = 1-ys/Eid(nb,nc,nd,nf) ; 
                   if (sum(ys<(4/N))==3) | ((nb==Nb) & (nc==Nc) & (nd==Nd) & (nf==Nf))
                      optidx(2,:) = [nb nc nd nf] ;
                   end 
                end 
                if (~sum(optidx(3,:)))	% F-test fitness (VA). 
                   ys = [Eva(nb-1,nc,nd,nf) Eva(nb,nc-1,nd,nf) ...
                         Eva(nb,nc,nd-1,nf) Eva(nb,nc,nd,nf-1)]; 
                   ys = 1-ys/Eva(nb,nc,nd,nf) ; 
                   if (sum(ys<(4/N))==3) | ((nb==Nb) & (nc==Nc) & (nd==Nd) & (nf==Nf))
                      optidx(3,:) = [nb nc nd nf] ;
                   end 
                end 
             end 
             
             % Afisare Grafica
             if (pf)			
                figure(FIG),clf ;
                   fig_look(FIG,1.5) ; 
                   subplot(321)
                      plot(1:N,Did.y,'-b',1:N,Yid{nb,nc,nd,nf},'-r') ; 
                      title('Identification data') ; 
                      ylabel('Outputs') ; 
                      ys = [min(min(Did.y),min(Yid{nb,nc,nd,nf})) ... 
                            max(max(Did.y),max(Yid{nb,nc,nd,nf}))] ; 
                      dy = ys(2)-ys(1) ; 
                      axis([0 N+1 ys(1)-0.05*dy ys(2)+0.2*dy]) ;
                      text(1.22*N,ys(2)+0.45*dy,['MCMMPE(BJ) nb=' int2str(nb-1) ' nc=' int2str(nc-1) ...
                             ' nd=' int2str(nd-1) ' nf=' int2str(nf-1)]) ; 
                      text(N/2,ys(2)+0.05*dy, ... 
                           ['Fitness E_N = ' sprintf('%g',Eid(nb,nc,nd,nf)) ' %']) ;
                   subplot(322)
                      plot(1:N,Dva.y,'-b',1:N,Yva{nb,nc,nd,nf},'-r') ; 
                      title('Validation data') ; 
                      ylabel('Outputs') ; 
                      ys = [min(min(Dva.y),min(Yva{nb,nc,nd,nf})) ... 
                            max(max(Dva.y),max(Yva{nb,nc,nd,nf}))] ; 
                      dy = ys(2)-ys(1) ; 
                      axis([0 N+1 ys(1)-0.05*dy ys(2)+0.2*dy]) ;
                      text(N/2,ys(2)+0.05*dy, ... 
                           ['Fitness E_N = ' sprintf('%g',Eva(nb,nc,nd,nf)) ' %']) ;
                      legend('y','ym') ; 
                   subplot(323)
                      plot(1:N,PEid{nb,nc,nd,nf}.y,'-m') ; 
                      ylabel('Prediction error') ; 
                      title(['Lambda^2 = ' sprintf('%g',std(PEid{nb,nc,nd,nf}.y,1)^2)]);
                   subplot(324)
                      plot(1:N,PEva{nb,nc,nd,nf}.y,'-m') ; 
                      ylabel('Prediction error') ; 
                      title(['Lambda^2 = ' sprintf('%g',std(PEva{nb,nc,nd,nf}.y,1)^2)]);
                   subplot(325)
                      [r,K] = xcov(PEid{nb,nc,nd,nf}.y,'unbiased') ; 
                      r = r(K>=0) ; K = ceil(length(r)/2) ; r = r(1:K) ; 
                      stem(1:K,r,'-g','filled') ; 
                      ylabel('Auto-covariance') ; 
                      title(['Val. Idx (ID) = ' num2str(Viid(nb,nc,nd,nf))]);
                   subplot(326)
                      [r,K] = xcov(PEva{nb,nc,nd,nf}.y,'unbiased') ; 
                      r = r(K>=0) ; K = ceil(length(r)/2) ; r = r(1:K) ; 
                      stem(1:K,r,'-g','filled') ; 
                      ylabel('Auto-covariance') ; 
                      title(['Val. Idx (VA) = ' num2str(Viva(nb,nc,nd,nf))]);
                      text(1.3*K,min(r),'<Press a key>') ;
                
                FIG = FIG+1 ; 
                % pause ; % Decomenteaza pentru pas-cu-pas
                
                % Poli-Zerouri SISTEM (Input/Output)
                figure(FIG),clf ;
                   fig_look(FIG,2) ; 
                   % idpoly returneaza F si B. Pentru iopzplot/map, mapam la a si b
                   P_sys = Mid; 
                   P_sys.a = Mid.f; % Polinomul F este numitorul sistemului
                   P_sys.b = Mid.b; % Polinomul B este numaratorul sistemului
                   P_sys.c = 1; P_sys.d = 1; % Ignoram zgomotul aici
                   
                   iopzmap(P_sys,'b','SD',alpha) ; 
                   title('Poles-Zeros representation (system F, B)') ; 
                   xlabel('Real axis') ; ylabel('Imaginary axis') ; 
                   text(0,0,'<Press a key>') ; 
                FIG = FIG+1 ; 
                % pause ;
                
                % Poli-Zerouri ZGOMOT
                   figure(FIG),clf ;
                      fig_look(FIG,2) ; 
                      P_noise = Mid;
                      P_noise.a = Mid.d; % Polinomul D este numitorul zgomotului
                      P_noise.b = Mid.c; % Polinomul C este numaratorul zgomotului
                      P_noise.c = 1; P_noise.d = 1; P_noise.f = 1; P_noise.nk = 0;
                      
                      iopzmap(P_noise,'r','SD',alpha) ; 
                      title('Poles-Zeros representation (noise D, C)') ; 
                      xlabel('Real axis') ; ylabel('Imaginary axis') ; 
                      text(0,0,'<Press a key>') ;
                   % pause ;
                FIG = FIG-2 ; 
             end 
          end 
      end 
    end 
  end 
end 

if (~pf)
   war_err([blanks(6) '... Done.']) ; 
end 

% 
% Proposing the optimal structure
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Nota: GAIC_R4 trebuie sa suporte 4 dimensiuni. Daca nu ai GAIC_R4, foloseste GAIC_R3 adaptat
[nb,nc,nd,nf] = GAIC_R4(Lambda,N) ; 
optidx = [optidx ; [nb nc,nd,nf]] ; 
war_err('* Proposed optimal indices (MCMMPE - BJ):') ; 
disp(['<F-test pred err>: [nb nc nd nf] = [' sprintf(' %d',optidx(1,:)-1) ']']) ; 
disp(['<F-test fit ID>:   [nb nc nd nf] = [' sprintf(' %d',optidx(2,:)-1) ']']) ; 
disp(['<F-test fit VA>:   [nb nc nd nf] = [' sprintf(' %d',optidx(3,:)-1) ']']) ; 
disp(['<GAIC-Rissanen>:   [nb nc nd nf] = [' sprintf(' %d',optidx(4,:)-1) ']']) ; 

% 
% Selecting the optimal structure
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('# Insert optimal indices: ') ;
nb = input('nb: ');
nc = input('nc: ');
nd = input('nd: ');
nf = input('nf: ');
% Fallback
if isempty(nb), nb=2; end; if isempty(nc), nc=2; end
if isempty(nd), nd=2; end; if isempty(nf), nf=2; end

nf = min(Nf,nf)+1 ; 
nd = min(Nd,nd)+1 ;
nc = min(Nc,nc)+1 ; 
nb = min(Nb,nb)+1 ;
% optidx = [optidx ; [nb nc nd nf]-1] ; 

Mid = M{nb,nc,nd,nf} ; 
Yid = Yid{nb,nc,nd,nf} ; 
Yva = Yva{nb,nc,nd,nf} ; 
PEid = PEid{nb,nc,nd,nf}.y ; 
PEva = PEva{nb,nc,nd,nf}.y ; 
Eid = Eid(nb,nc,nd,nf) ; 
Eva = Eva(nb,nc,nd,nf) ; 
Viid = Viid(nb,nc,nd,nf) ; 
Viva = Viva(nb,nc,nd,nf) ; % Corectat indexarea aici (era nb,nc)

war_err('o Optimum model (MCMMPE - BJ): ') ; 
Mid
war_err([blanks(25) '<Press a key>']) ; 
pause

% 
% Plotting model performances (FINAL)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
% (Blocul de afisare finala - identic cu cel din bucla)
 figure(FIG),clf ;
fig_look(FIG,1.5) ; 
subplot(321)
  plot(1:N,Did.y,'-b',1:N,Yid,'-r') ; 
  title('Identification data') ; 
  ylabel('Outputs') ; 
  ys = [min(min(Did.y),min(Yid)) max(max(Did.y),max(Yid))] ; 
  dy = ys(2)-ys(1) ; 
  axis([0 N+1 ys(1)-0.05*dy ys(2)+0.2*dy]) ;
  text(N/2,ys(2)+0.05*dy, ['Fitness E_N = ' sprintf('%g',Eid) ' %']) ;
subplot(322)
  plot(1:N,Dva.y,'-b',1:N,Yva,'-r') ; 
  title('Validation data') ; 
  ylabel('Outputs') ; 
  ys = [min(min(Dva.y),min(Yva)) max(max(Dva.y),max(Yva))] ; 
  dy = ys(2)-ys(1) ; 
  axis([0 N+1 ys(1)-0.05*dy ys(2)+0.2*dy]) ;
  text(N/2,ys(2)+0.05*dy, ['Fitness E_N = ' sprintf('%g',Eva) ' %']) ;
  legend('y','ym') ; 
subplot(323)
  plot(1:N,PEid,'-m') ; ylabel('Prediction error');
  title(['Lambda^2 = ' sprintf('%g',std(PEid,1)^2)]) ;
subplot(324)
  plot(1:N,PEva,'-m') ; ylabel('Prediction error');
  title(['Lambda^2 = ' sprintf('%g',std(PEva,1)^2)]) ; 
subplot(325)
  [r,K] = xcov(PEid,'unbiased') ; r = r(K>=0) ; K = ceil(length(r)/2) ; r = r(1:K) ; 
  stem(1:K,r,'-g','filled') ; ylabel('Auto-covariance') ; 
  title(['Val. Idx = ' num2str(Viid)]) ; 
subplot(326)
  [r,K] = xcov(PEva,'unbiased') ; r = r(K>=0) ; K = ceil(length(r)/2) ; r = r(1:K) ; 
  stem(1:K,r,'-g','filled') ; ylabel('Auto-covariance') ; 
  title(['Val. Idx = ' num2str(Viva)]) ; 

FIG = FIG+1 ; 
pause ; 

% Plot final Poli-Zerouri
figure(FIG),clf ;
fig_look(FIG,2) ; 
P_sys = Mid; P_sys.a = Mid.f; P_sys.b = Mid.b; P_sys.c = 1; P_sys.d = 1;
iopzmap(P_sys,'b','SD',alpha) ; 
title('Poles-Zeros representation (system)') ; 

FIG = FIG+1 ; 
pause ;

figure(FIG),clf ;
fig_look(FIG,2) ; 
P_noise = Mid; P_noise.a = Mid.d; P_noise.b = Mid.c; P_noise.c = 1; P_noise.d = 1; P_noise.f=1; P_noise.nk=0;
iopzmap(P_noise,'r','SD',alpha) ; 
title('Poles-Zeros representation (noise)') ; 

FIG = FIG-2 ;
% END