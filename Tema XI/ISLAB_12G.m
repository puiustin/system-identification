function [R] = ISLAB_12G(mt,K0,T0,Tp0,Tmax,Ts,U)
%
% ISLAB_12G   Mini-simulator pentru identificarea parametrilor K, T si Tp
%             ai motorului de curent continuu cu parametrul parazit Tp.
%
% Model continuu:
%
%                  K
% H(s) = ------------------------
%        s*(1+T*s)*(1+Tp*s)
%
% Cerinta:
%   - se testeaza toate cele 3 tipuri de intrari;
%   - se testeaza valori crescatoare ale lui lambda;
%   - se observa conditiile in care Tp este identificat satisfacator.
%
% Intrari:
%   mt   # tip model de identificare:
%          0 -> OE
%          1 -> ARX
%   K0   # amplificare reala
%   T0   # constanta principala reala
%   Tp0  # constanta parazita reala
%   Tmax # durata simularii
%   Ts   # perioada de esantionare
%   U    # amplitudinea intrarii
%
% Iesire:
%   R    # tabel/structura cu rezultatele obtinute
%

%
% BEGIN
%

global FIG ;
FIG = 1 ;

%
% Faults preventing
% ~~~~~~~~~~~~~~~~~
if (nargin < 7)
   U = 0.5 ;
end
if (isempty(U))
   U = 0.5 ;
end
U = U(1) ;
if (abs(U)<eps)
   U = 0.5 ;
end

if (nargin < 6)
   Ts = 0.1 ;
end
if (isempty(Ts))
   Ts = 0.1 ;
end
Ts = abs(Ts(1)) ;
if (Ts<eps)
   Ts = 0.1 ;
end

if (nargin < 5)
   Tmax = 80 ;
end
if (isempty(Tmax))
   Tmax = 80 ;
end
Tmax = abs(Tmax(1)) ;
if (Tmax<eps)
   Tmax = 80 ;
end
Tmax = max(250*Ts,Tmax) ;

if (nargin < 4)
   Tp0 = 0.05 ;
end
if (isempty(Tp0))
   Tp0 = 0.05 ;
end
Tp0 = abs(Tp0(1)) ;
if (Tp0<eps)
   Tp0 = 0.05 ;
end

if (nargin < 3)
   T0 = 0.5 ;
end
if (isempty(T0))
   T0 = 0.5 ;
end
T0 = abs(T0(1)) ;
if (T0<eps)
   T0 = 0.5 ;
end

if (nargin < 2)
   K0 = 4 ;
end
if (isempty(K0))
   K0 = 4 ;
end
K0 = K0(1) ;
if (abs(K0)<eps)
   K0 = 4 ;
end

if (nargin < 1)
   mt = 0 ;
end
if (isempty(mt))
   mt = 0 ;
end
mt = abs(round(mt(1))) ;

%
% Seturi de test impuse de cerinta
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ip_values = [0 1 2] ;                 % cele 3 tipuri de intrari
lambda_values = [0.1 0.5 1 2 3] ;     % lambda crescator => SNR descrescator

nIP = length(ip_values) ;
nL  = length(lambda_values) ;

%
% Initializare rezultate
% ~~~~~~~~~~~~~~~~~~~~~~
R = struct ;
cnt = 1 ;

fprintf('\n<ISLAB_12G>: Identificarea parametrilor paraziti K, T si Tp\n') ;
fprintf('<ISLAB_12G>: Model real: K = %.4f, T = %.4f, Tp = %.4f\n',K0,T0,Tp0) ;

if (~mt)
   model_name = 'OE' ;
else
   model_name = 'ARX' ;
end

fprintf('<ISLAB_12G>: Model de identificare: %s\n\n',model_name) ;

fprintf(' ip   lambda      SNR        K_est      T_est      Tp_est     err_Tp      sigma_y^2\n') ;
fprintf('-----------------------------------------------------------------------------------\n') ;

%
% Simulari pentru toate intrarile si toate valorile lambda
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for i = 1:nIP

   ip = ip_values(i) ;

   for j = 1:nL

      lambda = lambda_values(j) ;

      %
      % Generare date I/O cu parametrul parazit Tp
      % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      [D,V,P] = IOdata_DCeng(ip,0,K0,T0,Tp0,Tmax,Ts,U,lambda) ;

      Tmax = Ts*round(Tmax/Ts) ;
      t = 0:Ts:Tmax ;

      %
      % Calcul SNR
      % ~~~~~~~~~~
      y0 = D.y ;
      e0 = V.y ;

      SNR = 10*log10(var(y0)/var(e0)) ;

      %
      % Identificare model discret de ordin 3
      % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (~mt)

         M = oe(D,[3 3 1]) ;

         A = M.f ;
         B = M.b ;

      else

         M = arx(D,[3 3 1]) ;

         A = M.a ;
         B = M.b ;

      end

      %
      % Pregatire coeficienti pentru Hz2Hs_DCeng
      % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      A = A(:).' ;
      B = B(:).' ;

      if (length(A) < 4)
         A = [zeros(1,4-length(A)) A] ;
      end

      if (length(B) < 3)
         B = [zeros(1,3-length(B)) B] ;
      end

      A4 = A(end-3:end) ;
      B3 = B(end-2:end) ;

      %
      % Conversie model discret -> parametri fizici
      % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      try

         [Kest,Test,Tpest] = Hz2Hs_DCeng(B3,A4,Ts) ;

      catch

         Kest = NaN ;
         Test = NaN ;
         Tpest = NaN ;

      end

      %
      % Simulare iesire model identificat
      % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ys = sim(M,D.u) ;

      if (isa(ys,'iddata'))
         ysim = ys.y ;
      else
         ysim = ys ;
      end

      vy = D.y - ysim ;
      sigma2_y = var(vy) ;

      errTp = abs(Tpest - Tp0) ;

      %
      % Salvare rezultate
      % ~~~~~~~~~~~~~~~~~
      R(cnt).ip = ip ;
      R(cnt).lambda = lambda ;
      R(cnt).SNR = SNR ;
      R(cnt).K_real = K0 ;
      R(cnt).T_real = T0 ;
      R(cnt).Tp_real = Tp0 ;
      R(cnt).K_est = Kest ;
      R(cnt).T_est = Test ;
      R(cnt).Tp_est = Tpest ;
      R(cnt).err_Tp = errTp ;
      R(cnt).sigma2_y = sigma2_y ;
      R(cnt).model = model_name ;

      fprintf(' %1d   %6.2f   %8.4f   %8.4f   %8.4f   %8.4f   %8.4f   %10.4f\n', ...
              ip,lambda,SNR,Kest,Test,Tpest,errTp,sigma2_y) ;

      cnt = cnt + 1 ;

   end

end

fprintf('-----------------------------------------------------------------------------------\n') ;

%
% Alegerea celui mai bun caz pentru Tp
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
err = [R.err_Tp] ;
[~,best_idx] = min(err) ;

best = R(best_idx) ;

fprintf('\n<ISLAB_12G>: Cel mai bun caz pentru identificarea lui Tp:\n') ;
fprintf('ip = %d, lambda = %.2f, SNR = %.4f dB\n',best.ip,best.lambda,best.SNR) ;
fprintf('K_est = %.4f, T_est = %.4f, Tp_est = %.4f\n',best.K_est,best.T_est,best.Tp_est) ;
fprintf('eroare Tp = %.4f\n',best.err_Tp) ;
fprintf('sigma_y^2 = %.4f\n\n',best.sigma2_y) ;

%
% Refacerea celui mai bun caz pentru grafice
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[D,V,P] = IOdata_DCeng(best.ip,0,K0,T0,Tp0,Tmax,Ts,U,best.lambda) ;

if (~mt)

   M = oe(D,[3 3 1]) ;
   A = M.f ;
   B = M.b ;

else

   M = arx(D,[3 3 1]) ;
   A = M.a ;
   B = M.b ;

end

A = A(:).' ;
B = B(:).' ;

if (length(A) < 4)
   A = [zeros(1,4-length(A)) A] ;
end

if (length(B) < 3)
   B = [zeros(1,3-length(B)) B] ;
end

A4 = A(end-3:end) ;
B3 = B(end-2:end) ;

[Kest,Test,Tpest] = Hz2Hs_DCeng(B3,A4,Ts) ;

ys = sim(M,D.u) ;

if (isa(ys,'iddata'))
   ysim = ys.y ;
else
   ysim = ys ;
end

vy = D.y - ysim ;
sigma2_y = var(vy) ;

t = 0:Ts:Tmax ;

%
% Grafic 1: intrare si iesire masurata
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure(FIG),clf
   fig_look(FIG,1.5) ;
   plot(t,D.u,'-b',t,D.y,'-r') ;
   FN = scaling([D.u D.y]) ;
   axis([0 Tmax FN]) ;
   title('Date intrare-iesire pentru motorul cu parametrul parazit Tp') ;
   xlabel('Time [s]') ;
   ylabel('Magnitude') ;
   legend('intrare','iesire masurata') ;
FIG = FIG+1 ;

%
% Grafic 2: iesire masurata vs iesire simulata
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure(FIG),clf
   fig_look(FIG,1.5) ;
   plot(t,ysim,'-b',t,D.y,'-r',t,ysim,'-b') ;
   FN = scaling([ysim D.y]) ;
   axis([0 Tmax FN]) ;
   title(['Iesire simulata versus iesire masurata - ' model_name]) ;
   xlabel('Time [s]') ;
   ylabel('Magnitude') ;
   legend('iesire simulata','iesire masurata') ;
FIG = FIG+1 ;

%
% Grafic 3: eroare pe iesire
% ~~~~~~~~~~~~~~~~~~~~~~~~~~
figure(FIG),clf
   fig_look(FIG,1.5) ;
   plot(t,vy,'-k') ;
   FN = scaling(vy) ;
   axis([0 Tmax FN]) ;
   title(['Eroare pe iesire, sigma_y^2 = ' num2str(sigma2_y)]) ;
   xlabel('Time [s]') ;
   ylabel('v_y') ;
   legend(['sigma_y^2 = ' num2str(sigma2_y)]) ;
FIG = FIG+1 ;

%
% Grafic 4: parametri reali si estimati
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure(FIG),clf
   fig_look(FIG,1.5) ;
   bar([K0 Kest ; T0 Test ; Tp0 Tpest]) ;
   grid ;
   title('Parametri fizici reali si estimati') ;
   ylabel('Valoare') ;
   set(gca,'XTickLabel',{'K','T','Tp'}) ;
   legend('real','estimat') ;
FIG = FIG+1 ;

%
% Observatii automate
% ~~~~~~~~~~~~~~~~~~~
fprintf('<ISLAB_12G>: Observatie:\n') ;
fprintf('Intrarea cu ip = 2 reprezinta SPAB si are ordin de persistenta mai mare.\n') ;
fprintf('In general, estimarea lui Tp este mai buna pentru intrari persistente si lambda mic.\n') ;
fprintf('Cand lambda creste, SNR scade si eroarea de identificare a lui Tp creste.\n') ;

%
% END
%