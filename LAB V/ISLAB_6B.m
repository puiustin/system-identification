function [Mid,Did,Dva] = ISLAB_6B(B,C,F,D,nk,N,sigma,lambda) 
%
% ISLAB_6B   Simulator pentru identificarea Box-Jenkins (MMEP sau MCMMPE).
%            Include grafica corectata (stil ISLAB_6A) si GAIC_R4.
%
% Inputs:    B, C, F, D - Polinoamele procesului BJ
%            nk, N, sigma, lambda - Parametri simulare
%
% Outputs:   Mid, Did, Dva
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
pf = 1 ; 			% Plot flag: 0=no, 1=yes. 
Nb = 3 ;			% Maximum index of B part. 
Nc = 3 ;			% Maximum index of C part. 
Nf = 3 ;			% Maximum index of F part.
Nd = 3 ;			% Maximum index of D part. 
Ts = 1 ; 			% Sampling period. 

% 
% Faults preventing (Validare Parametri)
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
% Constructing the data provider process (Box-Jenkins)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% P = idpoly(A, B, C, D, F, ...) unde A=1
P = idpoly(1,[zeros(1,nk) B],C,D,F,lambda,Ts) ; 

% 
% Generating the identification data 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Did = gen_data(P,N,sigma,lambda) ; 

% 
% Generating the validation data
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Dva = gen_data(P,N,sigma,lambda) ; 

% 
% Estimating all models
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Crestem limitele pentru bucle (indexare MATLAB 1...N+1)
Nb = Nb+1 ; Nc = Nc+1 ; Nf = Nf+1; Nd = Nd+1;

M = cell(Nb,Nc,Nd,Nf) ; 			
Lambda = 1000*lambda*ones(Nb,Nc,Nd,Nf) ; 	
Yid = cell(Nb,Nc,Nd,Nf) ; 			
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
   war_err([blanks(4) '* Box-Jenkins Estimation started. Please wait...']) ; 
end 

% --- BUCLELE PRINCIPALE ---
for nf=2:Nf
  if (~pf)
      % AICI ERA EROAREA: Afisam 'nf', nu 'nb'
      disp([blanks(10) 'nf = ' int2str(nf-1)]) ; 
  end 
  for nd=2:Nd
    if (~pf)
        % AICI ERA EROAREA: Afisam 'nd', nu 'nc'
        disp([blanks(16) 'nd = ' int2str(nd-1)]) ;
    end 
    for nb=2:Nb
      if (~pf)
         disp([blanks(22) 'nb = ' int2str(nb-1)]) ;
      end 
      for nc=2:Nc
          if (~pf)
             disp([blanks(22) 'nc = ' int2str(nc-1)]) ;
          end 
      
          % Conditia de estimare (macar un parametru > 0)
          if ((nb>1) || (nc>1) || (nd>1) || (nf>1))
             
             % ----------------------------------------------------
             % ESTIMAREA MODELULUI
             % ----------------------------------------------------
             % Folosim 'bj' standard. 
             % Daca vrei MCMMPE, schimba in: Mid = bj_e(Did, [nb-1 nc-1 nd-1 nf-1 nk]);
             Mid = bj(Did, [nb-1 nc-1 nd-1 nf-1 nk]) ; 
             
             % Salvare model si performante
             M{nb,nc,nd,nf} = Mid ; 		
             Lambda(nb,nc,nd,nf) = Mid.NoiseVariance ; 
             
             PEid{nb,nc,nd,nf} = resid(Mid,Did) ;
             PEva{nb,nc,nd,nf} = resid(Mid,Dva) ; 
             ys = compare(Mid,Did) ;	
             Yid{nb,nc,nd,nf} = ys.y ; 
             ys = compare(Mid,Dva) ;        
             Yva{nb,nc,nd,nf} = ys.y ;    
    					
             Eid(nb,nc,nd,nf) = 100*(1-norm(PEid{nb,nc,nd,nf}.y)/norm(yNid)) ; 
             Eva(nb,nc,nd,nf) = 100*(1-norm(PEva{nb,nc,nd,nf}.y)/norm(yNva)) ; 
             
             Viid(nb,nc,nd,nf) = valid_LS(Mid,Did) ; 
             Viva(nb,nc,nd,nf) = valid_LS(Mid,Dva) ; 
             
             % Actualizare indici optimi (Testele F)
             if ((nb>1) && (nc>1) && (nd>1) && (nf>1))	
                if (~sum(optidx(1,:)))	% F-test pred error
                   ys = [Lambda(nb-1,nc,nd,nf) Lambda(nb,nc-1,nd,nf) ...
                         Lambda(nb,nc,nd-1,nf) Lambda(nb,nc,nd,nf-1)] ; 
                   ys = ys/Lambda(nb,nc,nd,nf) - 1 ; 
                   if (sum(ys<(4/N))==3) || ((nb==Nb) && (nc==Nc) && (nd==Nd) && (nf==Nf))
                      optidx(1,:) = [nb nc nd nf] ;
                   end 
                end 
                if (~sum(optidx(2,:)))	% F-test fitness (ID) 
                   ys = [Eid(nb-1,nc,nd,nf) Eid(nb,nc-1,nd,nf) ...
                         Eid(nb,nc,nd-1,nf) Eid(nb,nc,nd,nf-1)]; 
                   ys = 1-ys/Eid(nb,nc,nd,nf) ; 
                   if (sum(ys<(4/N))==3) || ((nb==Nb) && (nc==Nc) && (nd==Nd) && (nf==Nf))
                      optidx(2,:) = [nb nc nd nf] ;
                   end 
                end 
                if (~sum(optidx(3,:)))	% F-test fitness (VA)
                   ys = [Eva(nb-1,nc,nd,nf) Eva(nb,nc-1,nd,nf) ...
                         Eva(nb,nc,nd-1,nf) Eva(nb,nc,nd,nf-1)]; 
                   ys = 1-ys/Eva(nb,nc,nd,nf) ; 
                   if (sum(ys<(4/N))==3) || ((nb==Nb) && (nc==Nc) && (nd==Nd) && (nf==Nf))
                      optidx(3,:) = [nb nc nd nf] ;
                   end 
                end 
             end 
             
             % ----------------------------------------------------
             % AFISARE GRAFICA (Stil ISLAB_6A adaptat BJ)
             % ----------------------------------------------------
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
                      text(1.22*N,ys(2)+0.45*dy,['nb=' int2str(nb-1) ' nc=' int2str(nc-1) ...
                             ' nd=' int2str(nd-1) ' nf=' int2str(nf-1)]) ; 
                      text(N/2,ys(2)+0.05*dy, ... 
                           ['Fitness E_N = ' sprintf('%g',Eid(nb,nc,nd,nf)) ' %']) ;
                           
                   % Subplot 2: Validare
                   subplot(322)
                      plot(1:N,Dva.y,'-b',1:N,Yva{nb,nc,nd,nf},'-r') ; 
                      title('Validation data') ; 
                      ylabel('Outputs') ; 
                      ys = [min(min(Dva.y),min(Yva{nb,nc,nd,nf})) ... 
                            max(max(Dva.y),max(Yva{nb,nc,nd,nf}))] ; 
                      dy = ys(2)-ys(1) ; 
                      % CORECTIE: Fortam axa X sa fie exact 1:N
                      axis([1 N ys(1)-0.05*dy ys(2)+0.2*dy]) ; 
                      text(N/2,ys(2)+0.05*dy, ... 
                           ['Fitness E_N = ' sprintf('%g',Eva(nb,nc,nd,nf)) ' %']) ;
                      legend('y','ym') ; 
                      
                   % Subplot 3: Eroare predictie (ID)
                   subplot(323)
                      plot(1:N,PEid{nb,nc,nd,nf}.y,'-m') ; 
                      ylabel('PE (ID)') ; 
                      ys = [min(PEid{nb,nc,nd,nf}.y) max(PEid{nb,nc,nd,nf}.y)];
                      dy = ys(2)-ys(1);
                      % CORECTIE: Fortam axa X sa fie exact 1:N
                      axis([1 N ys(1)-0.05*dy ys(2)+0.2*dy]) ;
                      title(['\lambda^2 = ' sprintf('%g',std(PEid{nb,nc,nd,nf}.y,1)^2)]);
                      
                   % Subplot 4: Eroare predictie (VA)
                   subplot(324)
                      plot(1:N,PEva{nb,nc,nd,nf}.y,'-m') ; 
                      ylabel('PE (VA)') ; 
                      ys = [min(PEva{nb,nc,nd,nf}.y) max(PEva{nb,nc,nd,nf}.y)];
                      dy = ys(2)-ys(1);
                      % CORECTIE: Fortam axa X sa fie exact 1:N
                      axis([1 N ys(1)-0.05*dy ys(2)+0.2*dy]) ;
                      title(['\lambda^2 = ' sprintf('%g',std(PEva{nb,nc,nd,nf}.y,1)^2)]);
                      
                   % Subplot 5: Autocovarianta (ID)
                   subplot(325)
                      [r,K] = xcov(PEid{nb,nc,nd,nf}.y,'unbiased') ; 
                      r = r(K>=0) ; K = ceil(length(r)/2) ; r = r(1:K) ; 
                      stem(1:K,r,'-g','filled') ; 
                      ylabel('ACF (ID)') ; 
                      ys = [min(r) max(r)]; dy = ys(2)-ys(1);
                      % CORECTIE: Fortam axa X sa fie exact 1:K
                      axis([1 K ys(1)-0.05*dy ys(2)+0.2*dy]) ;
                      title(['Val. Idx = ' num2str(Viid(nb,nc,nd,nf))]);
                      
                   % Subplot 6: Autocovarianta (VA)
                   subplot(326)
                      [r,K] = xcov(PEva{nb,nc,nd,nf}.y,'unbiased') ; 
                      r = r(K>=0) ; K = ceil(length(r)/2) ; r = r(1:K) ; 
                      stem(1:K,r,'-g','filled') ; 
                      ylabel('ACF (VA)') ; 
                      ys = [min(r) max(r)]; dy = ys(2)-ys(1);
                      % CORECTIE: Fortam axa X sa fie exact 1:K
                      axis([1 K ys(1)-0.05*dy ys(2)+0.2*dy]) ;
                      title(['Val. Idx = ' num2str(Viva(nb,nc,nd,nf))]);
                      
                      % CORECTIE: Pozitionarea textului 'Press a key' mai spre dreapta
                      % Folosim coordonate relative la axa curenta
                      text(K*1.1, min(r), '<Press a key>', 'HorizontalAlignment', 'left') ;
                FIG = FIG+1 ; 
                % pause ; % Decomenteaza daca vrei pauza
                
                % Poli-Zerouri SISTEM (B/F)
                figure(FIG),clf ;
                   fig_look(FIG,2) ; 
                   % Mapam la obiect ARMA(Num=B, Den=F) pentru afisare
                   P_sys = idpoly(Mid.f, Mid.b, [], [], [], 1, 1);
                   iopzmap(P_sys,'b','SD',alpha) ; 
                   title('Poles-Zeros: System (B/F)') ; 
                   xlabel('Real axis') ; ylabel('Imaginary axis') ; 
                   
                FIG = FIG+1 ; 
                
                % Poli-Zerouri ZGOMOT (C/D)
               if ((nd>1) || (nc>1))
                   figure(FIG),clf ;
                      fig_look(FIG,2) ; 
                      P_noise = idpoly(Mid.d, Mid.c, [], [], [], 1, 1);
                      P_noise.nk = 0 ; 
                      iopzmap(P_noise,'r','SD',alpha) ; 
                      title('Poles-Zeros: Noise (C/D)') ; 
                      xlabel('Real axis') ; ylabel('Imaginary axis') ; 
                      
                      % CORECTIE: Luam limitele actuale ale graficului
                      ys = axis; 
                      % ys(2) este limita dreapta (X max), ys(3) este limita jos (Y min)
                      % Punem textul putin in afara cercului, dreapta-jos
                      text(ys(2)*0.8, ys(3)*0.9, '<Press a key>', ...
                           'HorizontalAlignment', 'left', 'Color', 'k');
                   % pause ;
                end
                FIG = FIG-2 ; 
             end  % [if (pf)]
          end 
      end 
    end 
  end 
end 

if (~pf)
   war_err([blanks(6) '... Done.']) ; 
end  

% 
% Propunerea structurii optime (GAIC_R4)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[nb,nc,nd,nf] = GAIC_R4(Lambda,N) ; 
optidx = [optidx ; [nb nc nd nf]] ; 

war_err('* Proposed optimal indices:') ; 
disp(['<F-test Pred Err>: [nb nc nd nf] = [' sprintf(' %d',optidx(1,:)-1) ']']) ; 
disp(['<F-test Fit ID>:   [nb nc nd nf] = [' sprintf(' %d',optidx(2,:)-1) ']']) ; 
disp(['<F-test Fit VA>:   [nb nc nd nf] = [' sprintf(' %d',optidx(3,:)-1) ']']) ; 
disp(['<GAIC-Rissanen>:   [nb nc nd nf] = [' sprintf(' %d',optidx(4,:)-1) ']']) ; 

% 
% Selectia finala
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp(' ') ;
disp('# Insert optimal indices for BJ: ') ; 
nb_in = input('nb: '); nc_in = input('nc: ');
nd_in = input('nd: '); nf_in = input('nf: ');

% Defaults
if isempty(nb_in), nb_in=2; end; if isempty(nc_in), nc_in=2; end
if isempty(nd_in), nd_in=2; end; if isempty(nf_in), nf_in=2; end

nf = min(Nf, abs(round(nf_in)) + 1) ; 
nd = min(Nd, abs(round(nd_in)) + 1) ; 
nc = min(Nc, abs(round(nc_in)) + 1) ; 
nb = min(Nb, abs(round(nb_in)) + 1) ; 

Mid = M{nb,nc,nd,nf} ; 
Yid = Yid{nb,nc,nd,nf} ; 
Yva = Yva{nb,nc,nd,nf} ; 
PEid = PEid{nb,nc,nd,nf}.y ; 
PEva = PEva{nb,nc,nd,nf}.y ; 
Eid = Eid(nb,nc,nd,nf) ; 
Eva = Eva(nb,nc,nd,nf) ; 
Viid = Viid(nb,nc,nd,nf) ; 
Viva = Viva(nb,nc,nd,nf) ; 

war_err('o Optimum model: ') ; 
Mid
war_err([blanks(25) '<Press a key>']) ; 
pause

% 
% Afisare finala (Identica cu cea din bucla, dar statica)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% Plotting model performances (FINAL - OPTIMUM MODEL)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure(FIG),clf
   fig_look(FIG,1.5) ; 
   
   % --- Subplot 321: Date Identificare ---
   subplot(321)
      plot(1:N,Did.y,'-b',1:N,Yid,'-r') ; 
      title('Identification data') ; ylabel('Outputs') ; 
      ys = [min(min(Did.y),min(Yid)) max(max(Did.y),max(Yid))] ; 
      dy = ys(2)-ys(1) ; 
      % CORECTIE: Limite exacte [1 N]
      axis([1 N ys(1)-0.05*dy ys(2)+0.2*dy]) ; 
      % Afisam structura (nb, nc...) in titlu sau text. 
      % Deoarece axa e fixa, textul din dreapta (1.22*N) nu se mai vede. 
      % Il punem sub titlu sau folosim titlul.
      text(1.22*N,ys(2)+0.45*dy, ['nb=' int2str(nb-1) ' | nc=' int2str(nc-1) ' | nd=' int2str(nd-1) ' | nf=' int2str(nf-1)], 'HorizontalAlignment', 'center') ; 
      text(N/2,ys(2)+0.05*dy, ['Fitness E_N = ' sprintf('%g',Eid) ' %'], 'HorizontalAlignment', 'center') ;

   % --- Subplot 322: Date Validare ---
   subplot(322)
      plot(1:N,Dva.y,'-b',1:N,Yva,'-r') ; 
      title('Validation data') ; ylabel('Outputs') ; 
      ys = [min(min(Dva.y),min(Yva)) max(max(Dva.y),max(Yva))] ; 
      dy = ys(2)-ys(1) ; 
      % CORECTIE: Limite exacte [1 N]
      axis([1 N ys(1)-0.05*dy ys(2)+0.2*dy]) ; 
      text(N/2,ys(2)+0.05*dy, ['Fitness E_N = ' sprintf('%g',Eva) ' %'], 'HorizontalAlignment', 'center') ;

   % --- Subplot 323: Eroare Predictie (ID) ---
   subplot(323)
      plot(1:N,PEid,'-m') ; ylabel('Prediction error') ; 
      title(['\lambda^2 = ' sprintf('%g',std(PEid,1)^2)]) ;
      % CORECTIE: Adaugam calcul limite si axis [1 N]
      ys = [min(PEid) max(PEid)]; dy = ys(2)-ys(1);
      axis([1 N ys(1)-0.05*dy ys(2)+0.2*dy]) ;

   % --- Subplot 324: Eroare Predictie (VA) ---
   subplot(324)
      plot(1:N,PEva,'-m') ; ylabel('Prediction error') ; 
      title(['\lambda^2 = ' sprintf('%g',std(PEva,1)^2)]) ; 
      % CORECTIE: Adaugam calcul limite si axis [1 N]
      ys = [min(PEva) max(PEva)]; dy = ys(2)-ys(1);
      axis([1 N ys(1)-0.05*dy ys(2)+0.2*dy]) ;

   % --- Subplot 325: Autocovarianta (ID) ---
   subplot(325)
      [r,K] = xcov(PEid,'unbiased') ; r = r(K>=0) ; K = ceil(length(r)/2) ; r = r(1:K) ; 
      stem(1:K,r,'-g','filled') ; ylabel('ACF (ID)') ; 
      title(['Val. Idx = ' num2str(Viid)]) ; 
      % CORECTIE: Adaugam calcul limite si axis [1 K]
      ys = [min(r) max(r)]; dy = ys(2)-ys(1);
      axis([1 K ys(1)-0.05*dy ys(2)+0.2*dy]) ;

   % --- Subplot 326: Autocovarianta (VA) ---
   subplot(326)
      [r,K] = xcov(PEva,'unbiased') ; r = r(K>=0) ; K = ceil(length(r)/2) ; r = r(1:K) ; 
      stem(1:K,r,'-g','filled') ; ylabel('ACF (VA)') ; 
      title(['Val. Idx = ' num2str(Viva)]) ; 
      % CORECTIE: Adaugam calcul limite si axis [1 K]
      ys = [min(r) max(r)]; dy = ys(2)-ys(1);
      axis([1 K ys(1)-0.05*dy ys(2)+0.2*dy]) ;
FIG = FIG+1 ; 
pause ;
figure(FIG),clf
   fig_look(FIG,2) ; 
   P_sys = idpoly(Mid.f, Mid.b, [], [], [], 1, 1);
   iopzmap(P_sys,'b','SD',alpha) ; 
   title('System Poles-Zeros (B/F)') ; 
   
FIG = FIG+1 ; 
pause ;
if ((nd>1) || (nc>1))
   P_noise = idpoly(Mid.d, Mid.c, [], [], [], 1, 1); P_noise.nk = 0 ; 
   figure(FIG),clf ;
      fig_look(FIG,2) ; 
      iopzmap(P_noise,'r','SD',alpha) ; 
      title('Noise Poles-Zeros (C/D)') ; 
end  
% END