function [F,ID,SD] = ISLAB_12F(K0,T0,Tmax,Ts,U,lambda)
%
% ISLAB_12F   Identificarea recursiva a parametrilor variabili K si T
%             folosind rpem si reprezentarea pe stare.
%

%
% BEGIN
%

global FIG ;
FIG = 1 ;

cv = 1 ;

alpha = 0.1 ;
beta  = 0.1 ;

FN = '<ISLAB_12F>: ' ;
WB = [FN 'Recursive state-space estimation using rpem. Please wait ...'] ;
WE = [blanks(length(FN)) '... Done.'] ;

if (nargin < 6)
   lambda = 1 ;
end
if (isempty(lambda))
   lambda = 1 ;
end
lambda = abs(lambda(1)) ;

if (nargin < 5)
   U = 0.5 ;
end
if (isempty(U))
   U = 0.5 ;
end

if (nargin < 4)
   Ts = 0.1 ;
end
if (isempty(Ts))
   Ts = 0.1 ;
end

if (nargin < 3)
   Tmax = 80 ;
end
if (isempty(Tmax))
   Tmax = 80 ;
end

if (nargin < 2)
   T0 = 0.5 ;
end
if (isempty(T0))
   T0 = 0.5 ;
end

if (nargin < 1)
   K0 = 4 ;
end
if (isempty(K0))
   K0 = 4 ;
end

war_err(WB) ;

[ID,V,P] = gdata_DCeng(cv,K0,T0,Tmax,Ts,U,lambda) ;

Tmax = Ts*round(Tmax/Ts) ;
t = 0:Ts:Tmax ;
N = length(ID.y) ;

%
% Identificare recursiva cu rpem
% Model discret de tip OE/BJ:
% B(q)/F(q), cu model de zgomot C(q)
%

na = 0 ;
nb = 2 ;
nc = 1 ;
nd = 0 ;
nf = 2 ;
nk = 1 ;

theta = rpem(ID,[na nb nc nd nf nk],'ff',0.999) ;

war_err(WE) ;

%
% Coeficientii estimati recursiv:
% theta = [b1 b2 c1 f1 f2]
%

b1 = theta(:,1) ;
b2 = theta(:,2) ;
c1 = theta(:,3) ;
f1 = theta(:,4) ;
f2 = theta(:,5) ;

%
% Parametrii fizici K si T
%

F.K = (b1+b2)./(1-f2)/Ts ;
F.T = Ts*(f2.*b1+b2)./(b1+b2)./(1-f2) ;

F.K(~isfinite(F.K)) = K0 ;
F.T(~isfinite(F.T)) = T0 ;
F.T(F.T<=0) = T0 ;

%
% Simulare iesire folosind modelul discret recursiv identificat
%

ys = zeros(N,1) ;

for k = 3:N

   ys(k) = -f1(k)*ys(k-1) - f2(k)*ys(k-2) + ...
            b1(k)*ID.u(k-1) + b2(k)*ID.u(k-2) ;

   if (~isfinite(ys(k)))
      ys(k) = ys(k-1) ;
   end

end

SD = iddata(ys,V.u,Ts) ;

%
% Grafice date I/O
%

figure(FIG),clf
   fig_look(FIG,1.5) ;
   plot(t,ID.u,'-b',t,ID.y,'-r') ;
   FN2 = scaling([ID.u ID.y]) ;
   axis([0 Tmax FN2]) ;
   title('Input-output data provided by a DC engine.') ;
   xlabel('Time [s]') ;
   ylabel('Magnitude') ;
   legend('input (square wave)','output') ;
FIG = FIG+1 ;

%
% Iesire simulata vs iesire masurata
%

figure(FIG),clf
   fig_look(FIG,1.5) ;
   plot(t,SD.y,'-b',t,ID.y,'-r',t,SD.y,'-b') ;
   FN2 = scaling([SD.y ID.y]) ;
   axis([0 Tmax FN2]) ;
   title('Output data provided by a DC engine and its state-space model.') ;
   xlabel('Time [s]') ;
   ylabel('Magnitude') ;
   legend('simulated output','measured output') ;
FIG = FIG+1 ;

%
% Variatia parametrilor fizici
%

figure(FIG),clf
   fig_look(FIG,1.5) ;

   subplot(211)
      plot(t,F.K,'-b',t,P.num{1},'-r') ;
      FN2 = scaling([F.K P.num{1}']) ;
      axis([0 Tmax FN2]) ;
      title('Physical parameters variation (DC engine).') ;
      xlabel('Time [s]') ;
      ylabel('Gain K') ;
      legend('estimated','true') ;

   subplot(212)
      plot(t,F.T,'-b',t,P.den{1},'-r') ;
      FN2 = scaling([F.T P.den{1}']) ;
      axis([0 Tmax FN2]) ;
      xlabel('Time [s]') ;
      ylabel('Time constant T [s]') ;

FIG = FIG+1 ;

%
% Variatia parametrilor in regim stationar
%

WEz = (t>(0.5*Tmax)) ;
tz = t(WEz) ;

figure(FIG),clf
   fig_look(FIG,1.5) ;

   subplot(211)
      plot(tz,F.K(WEz),'-b',tz,P.num{1}(WEz),'-r') ;
      FN2 = scaling([F.K(WEz) P.num{1}(WEz)']) ;
      axis([tz(1) Tmax FN2]) ;
      title('Physical parameters variation - steady-state (DC engine).') ;
      xlabel('Time [s]') ;
      ylabel('Gain K') ;
      legend('estimated','true') ;

   subplot(212)
      plot(tz,F.T(WEz),'-b',tz,P.den{1}(WEz),'-r') ;
      FN2 = scaling([F.T(WEz) P.den{1}(WEz)']) ;
      axis([tz(1) Tmax FN2]) ;
      xlabel('Time [s]') ;
      ylabel('Time constant T [s]') ;

FIG = FIG+1 ;

%
% Zgomotul pe iesire
%

vy = ID.y - SD.y ;
sigma2_y = var(vy) ;

figure(FIG),clf
   fig_look(FIG,1.5) ;
   plot(t,vy,'-r') ;
   FN2 = scaling(vy) ;
   axis([0 Tmax FN2]) ;
   title(['Output colored noise v_y, sigma_y^2 = ' num2str(sigma2_y)]) ;
   xlabel('Time [s]') ;
   ylabel('v_y') ;
   legend(['sigma_y^2 = ' num2str(sigma2_y)]) ;
FIG = FIG+1 ;

%
% Zgomotul parametrului K
%

trueK = P.num{1}(:) ;
vK = trueK - F.K ;
sigma2_K = var(vK) ;

figure(FIG),clf
   fig_look(FIG,1.5) ;
   plot(t,vK,'-b') ;
   FN2 = scaling(vK) ;
   axis([0 Tmax FN2]) ;
   title(['Parameter noise v_K, sigma_K^2 = ' num2str(sigma2_K)]) ;
   xlabel('Time [s]') ;
   ylabel('v_K') ;
   legend(['sigma_K^2 = ' num2str(sigma2_K)]) ;
FIG = FIG+1 ;

%
% Zgomotul parametrului T
%

trueT = P.den{1}(:) ;
vT = trueT - F.T ;
sigma2_T = var(vT) ;

figure(FIG),clf
   fig_look(FIG,1.5) ;
   plot(t,vT,'-k') ;
   FN2 = scaling(vT) ;
   axis([0 Tmax FN2]) ;
   title(['Parameter noise v_T, sigma_T^2 = ' num2str(sigma2_T)]) ;
   xlabel('Time [s]') ;
   ylabel('v_T') ;
   legend(['sigma_T^2 = ' num2str(sigma2_T)]) ;
FIG = FIG+1 ;

disp(' ') ;
disp([FN 'alpha = ' num2str(alpha)]) ;
disp([FN 'beta = ' num2str(beta)]) ;
disp([FN 'Noise polynomial C(q): 1 + c1*q^-1']) ;
disp([FN 'Mean c1 = ' num2str(mean(c1))]) ;
disp([FN 'Noise dispersions:']) ;
disp(['sigma_y^2 = ' num2str(sigma2_y)]) ;
disp(['sigma_K^2 = ' num2str(sigma2_K)]) ;
disp(['sigma_T^2 = ' num2str(sigma2_T)]) ;

%
% Identificarea zgomotului colorat v_y cu modele ARMA/AR/MA
%

best_lambda2 = inf ;
best_e = vy ;
best_vhat = vy ;
best_name = 'none' ;

rng('shuffle') ;
NID = iddata(vy,[],Ts) ;

for k = 1:10

   na_n = randi([1 20]) ;
   nc_n = randi([1 20]) ;

   try
      NM = armax(NID,[na_n nc_n]) ;
      eID = resid(NM,NID) ;
      e_n = eID.OutputData ;
      lambda2_n = var(e_n) ;

      if (lambda2_n < best_lambda2)
         best_lambda2 = lambda2_n ;
         best_e = e_n ;
         best_vhat = vy - e_n ;
         best_name = ['ARMA(' num2str(na_n) ',' num2str(nc_n) ')'] ;
      end

   catch
   end

end

for k = 1:5

   na_n = randi([1 20]) ;

   try
      rv = xcorr(vy,na_n,'biased') ;
      rv = rv(na_n+1:end) ;

      [Aar,lambda2_n] = levinson(rv,na_n) ;
      e_n = filter(Aar,1,vy) ;

      if (lambda2_n < best_lambda2)
         best_lambda2 = lambda2_n ;
         best_e = e_n ;
         best_vhat = vy - e_n ;
         best_name = ['AR(' num2str(na_n) ')'] ;
      end

   catch
   end

end

for k = 1:5

   nc_n = randi([1 20]) ;

   try
      NM = armax(NID,[0 nc_n]) ;
      eID = resid(NM,NID) ;
      e_n = eID.OutputData ;
      lambda2_n = var(e_n) ;

      if (lambda2_n < best_lambda2)
         best_lambda2 = lambda2_n ;
         best_e = e_n ;
         best_vhat = vy - e_n ;
         best_name = ['MA(' num2str(nc_n) ')'] ;
      end

   catch
   end

end

disp(' ') ;
disp([FN 'Best output noise model: ' best_name]) ;
disp([FN 'lambda^2 = ' num2str(best_lambda2)]) ;

%
% Iesire masurata vs iesire cu zgomot colorat estimat
%

yhat_noise = SD.y + best_vhat ;

figure(FIG),clf
   fig_look(FIG,1.5) ;
   plot(t,ID.y,'-r',t,yhat_noise,'-b') ;
   FN2 = scaling([ID.y yhat_noise]) ;
   axis([0 Tmax FN2]) ;
   title('Measured output and output with estimated colored noise') ;
   xlabel('Time [s]') ;
   ylabel('Magnitude') ;
   legend('measured output',best_name) ;
FIG = FIG+1 ;

%
% Zgomot alb estimat
%

figure(FIG),clf
   fig_look(FIG,1.5) ;
   plot(t,best_e,'-k') ;
   FN2 = scaling(best_e) ;
   axis([0 Tmax FN2]) ;
   title(['Estimated white noise e, lambda^2 = ' num2str(best_lambda2)]) ;
   xlabel('Time [s]') ;
   ylabel('e') ;
   legend(['lambda^2 = ' num2str(best_lambda2)]) ;
FIG = FIG+1 ;

%
% END
%