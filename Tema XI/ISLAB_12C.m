function [F,D,M] = ISLAB_12C(met,K0,T0,Tmax,Ts,U,lambda)
%
% ISLAB_12C   Mini-simulator that performs off-line identification of
%             physical parameters of a DC engine using state-space models.
%
% Inputs: met    # identification method:
%                    0 -> n4sid (default)
%                    1 -> pem
%         K0     # constant gain (4, by default)
%         T0     # constant time constant (0.5 s, by default)
%         Tmax   # simulation duration (80 s, by default)
%         Ts     # sampling period (0.1 s, by default)
%         U      # amplitude of input square wave (0.5, by default)
%         lambda # standard deviation of white noise (1, by default)
%
% Outputs: F     # estimated physical parameters:
%                    F.K -> gain
%                    F.T -> time constant
%         D      # IDDATA object representing the I/O data
%         M      # estimated state-space model
%

%
% BEGIN
%

global FIG ;
FIG = 1 ;

%
% Constants
% ~~~~~~~~~
% .........................................
% Freely chosen constants for state representation
% These values are changed here when testing alpha and beta
% .........................................
alpha = 0.5 ;
beta  = 0.5 ;
% .........................................

%
% Messages
% ~~~~~~~~
FN = '<ISLAB_12C>: ' ;
NFN = length(FN) ;
PK = [blanks(70) '<Press a key>'] ;
M1 = [FN 'Physical parameters:'] ;
M2 = [blanks(NFN+7) 'True' blanks(8) 'Estimated'] ;
M3 = [blanks(NFN) 'K: %8.4f' blanks(6) '%8.4f \n'] ;
M4 = [blanks(NFN) 'T: %8.4f' blanks(6) '%8.4f \n'] ;

%
% Faults preventing
% ~~~~~~~~~~~~~~~~~
if (nargin < 7)
   lambda = 1 ;
end
if (isempty(lambda))
   lambda = 1 ;
end
lambda = abs(lambda(1)) ;
if (~lambda)
   lambda = 1 ;
end

if (nargin < 6)
   U = 0.5 ;
end
if (isempty(U))
   U = 0.5 ;
end
U = U(1) ;
if (abs(U)<eps)
   U = 0.5 ;
end

if (nargin < 5)
   Ts = 0.1 ;
end
if (isempty(Ts))
   Ts = 0.1 ;
end
Ts = abs(Ts(1)) ;
if (Ts<eps)
   Ts = 0.1 ;
end

if (nargin < 4)
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

if (nargin < 3)
   T0 = 0.5 ;
end
if (isempty(T0))
   T0 = 0.5 ;
end
T0 = T0(1) ;
if (abs(T0)<eps)
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
   met = 0 ;
end
if (isempty(met))
   met = 0 ;
end
met = abs(round(met(1))) ;

%
% Generate identification data
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[D,V] = gdata_DCeng(0,K0,T0,Tmax,Ts,U,lambda) ;

%
% State-space identification
% ~~~~~~~~~~~~~~~~~~~~~~~~~~
% .........................................
% met = 0 -> n4sid
% met = 1 -> pem
% .........................................
if (~met)

   try
      opt = n4sidOptions ;
      opt.Focus = 'simulation' ;
      opt.EnforceStability = true ;
      M = n4sid(D,2,opt) ;
   catch
      M = n4sid(D,2) ;
   end

else

   try
      opt1 = n4sidOptions ;
      opt1.Focus = 'simulation' ;
      opt1.EnforceStability = true ;
      M0 = n4sid(D,2,opt1) ;

      opt2 = pemOptions ;
      opt2.Focus = 'simulation' ;
      opt2.EnforceStability = true ;
      M = pem(D,M0,opt2) ;
   catch
      M0 = n4sid(D,2) ;
      M = pem(D,M0) ;
   end

end

%
% Convert state-space model to transfer function
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% .........................................
% The physical parameters K and T are recovered from
% the discrete transfer function obtained from the
% identified state-space model.
% .........................................
G = tf(M) ;
[num,den] = tfdata(G,'v') ;

if (iscell(num))
   num = num{1} ;
end
if (iscell(den))
   den = den{1} ;
end

num = num(:).' ;
den = den(:).' ;

if (abs(den(1)) > eps)
   num = num/den(1) ;
   den = den/den(1) ;
end

if (length(num) < 3)
   num = [zeros(1,3-length(num)) num] ;
end

if (length(den) < 3)
   den = [zeros(1,3-length(den)) den] ;
end

b1 = num(end-1) ;
b2 = num(end) ;
a2 = den(end) ;

if (abs(b1+b2)<eps)
   F.K = NaN ;
   F.T = NaN ;
else
   F.K = (b1+b2)/(1-a2)/Ts ;
   F.T = Ts*(a2*b1+b2)/(b1+b2)/(1-a2) ;
end

%
% Build physical state-space representation
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% .........................................
% Continuous physical state-space form:
%
% x1_dot = theta2*x1 + alpha*u
% x2_dot = theta1*x1
% y      = beta*x2
%
% T = -1/theta2
% K = -alpha*beta*theta1/theta2
% .........................................
if (~isnan(F.K) && ~isnan(F.T) && F.T > eps)

   theta2 = -1/F.T ;
   theta1 = F.K/(alpha*beta*F.T) ;

   Ac = [theta2 0 ; theta1 0] ;
   Bc = [alpha ; 0] ;
   Cc = [0 beta] ;
   Dc = 0 ;

   Pc = ss(Ac,Bc,Cc,Dc) ;
   Pd = c2d(Pc,Ts,'zoh') ;

   Tmax = Ts*round(Tmax/Ts) ;
   t = 0:Ts:Tmax ;
   ysim = lsim(Pd,D.u,t') ;

else

   Tmax = Ts*round(Tmax/Ts) ;
   t = 0:Ts:Tmax ;

   ytmp = sim(M,D.u) ;

   if (isa(ytmp,'iddata'))
      ysim = ytmp.y ;
   else
      ysim = ytmp ;
   end

   Ac = NaN ;
   Bc = NaN ;
   Cc = NaN ;

end

%
% Display physical parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
war_err(M1) ;
disp(M2) ;
fprintf(1,M3,K0,F.K) ;
fprintf(1,M4,T0,F.T) ;

disp(' ') ;
if (~met)
   disp([FN 'State-space identification method: n4sid']) ;
else
   disp([FN 'State-space identification method: pem']) ;
end

disp([FN 'alpha = ' num2str(alpha)]) ;
disp([FN 'beta  = ' num2str(beta)]) ;

disp(' ') ;
disp([FN 'Physical state-space model:']) ;
disp('A = ') ;
disp(Ac) ;
disp('B = ') ;
disp(Bc) ;
disp('C = ') ;
disp(Cc) ;

try
   disp(' ') ;
   disp([FN 'Endogenous noise matrix K from state-space model:']) ;
   disp(M.K) ;
catch
   disp(' ') ;
   disp([FN 'No endogenous noise matrix was available.']) ;
end

war_err(PK) ;
pause ;

%
% Plot I/O data
% ~~~~~~~~~~~~~
figure(FIG),clf
   fig_look(FIG,1.5) ;
   plot(t,D.u,'-b',t,D.y,'-r') ;
   FN2 = scaling([D.u D.y]) ;
   axis([0 Tmax FN2]) ;
   title('Input-output data provided by a DC engine.') ;
   xlabel('Time [s]') ;
   ylabel('Magnitude') ;
   legend('input','output') ;
FIG = FIG+1 ;

%
% Plot measured output versus simulated output
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure(FIG),clf
   fig_look(FIG,1.5) ;
   plot(t,ysim,'-b',t,D.y,'-r',t,ysim,'-b') ;
   FN2 = scaling([ysim D.y]) ;
   axis([0 Tmax FN2]) ;
   title('Measured output versus simulated output') ;
   xlabel('Time [s]') ;
   ylabel('Magnitude') ;
   legend('simulated output','measured output') ;
FIG = FIG+1 ;

%
% Colored noise
% ~~~~~~~~~~~~~
% .........................................
% Difference between measured output and simulated output
% .........................................
v = D.y - ysim ;
sigma2 = var(v) ;

figure(FIG),clf
   fig_look(FIG,1.5) ;
   plot(t,v,'-r') ;
   FN2 = scaling(v) ;
   axis([0 Tmax FN2]) ;
   title(['Colored noise v, sigma^2 = ' num2str(sigma2)]) ;
   xlabel('Time [s]') ;
   ylabel('v') ;
   legend(['sigma^2 = ' num2str(sigma2)]) ;
FIG = FIG+1 ;

%
% Exogenous noise identification
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% .........................................
% pem and n4sid provide the endogenous noise matrix.
% The exogenous output noise is identified here
% with AR, MA and ARMA models.
% .........................................
best_lambda2 = inf ;
best_e = v ;
best_vhat = v ;
best_name = 'none' ;

NID = iddata(v,[],Ts) ;

rng('shuffle') ;

% .........................................
% 10 ARMA models
% .........................................
for k = 1:10

   na_n = randi([1 20]) ;
   nc_n = randi([1 20]) ;

   try

      NM = armax(NID,[na_n nc_n]) ;

      e_n = resid(NM,NID) ;
      e_n = e_n.OutputData ;

      lambda2_n = var(e_n) ;

      if (lambda2_n < best_lambda2)

         best_lambda2 = lambda2_n ;
         best_e = e_n ;
         best_vhat = v - e_n ;
         best_name = ['ARMA(' num2str(na_n) ',' num2str(nc_n) ')'] ;

      end

   catch

   end

end

% .........................................
% 5 MA models
% .........................................
for k = 1:5

   nc_n = randi([1 20]) ;

   try

      NM = armax(NID,[0 nc_n]) ;

      e_n = resid(NM,NID) ;
      e_n = e_n.OutputData ;

      lambda2_n = var(e_n) ;

      if (lambda2_n < best_lambda2)

         best_lambda2 = lambda2_n ;
         best_e = e_n ;
         best_vhat = v - e_n ;
         best_name = ['MA(' num2str(nc_n) ')'] ;

      end

   catch

   end

end

% .........................................
% 5 AR models
% .........................................
for k = 1:5

   na_n = randi([1 20]) ;

   try

      rv = xcorr(v,na_n,'biased') ;
      rv = rv(na_n+1:end) ;

      [Aar,lambda2_n] = levinson(rv,na_n) ;

      e_n = filter(Aar,1,v) ;

      if (lambda2_n < best_lambda2)

         best_lambda2 = lambda2_n ;
         best_e = e_n ;
         best_vhat = v - e_n ;
         best_name = ['AR(' num2str(na_n) ')'] ;

      end

   catch

   end

end

disp(' ') ;
disp([FN 'Best exogenous noise model: ' best_name]) ;
disp([FN 'lambda^2 = ' num2str(best_lambda2)]) ;

%
% Output with estimated colored noise
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% .........................................
% Measured output versus simulated output corrected
% with the estimated colored noise
% .........................................
yhat_noise = ysim + best_vhat ;

figure(FIG),clf
   fig_look(FIG,1.5) ;
   plot(t,D.y,'-r',t,yhat_noise,'-b') ;
   FN2 = scaling([D.y yhat_noise]) ;
   axis([0 Tmax FN2]) ;
   title('Measured output and output with estimated noise') ;
   xlabel('Time [s]') ;
   ylabel('Magnitude') ;
   legend('measured output',best_name) ;
FIG = FIG+1 ;

%
% Estimated white noise
% ~~~~~~~~~~~~~~~~~~~~~
% .........................................
% Estimated white noise obtained after filtering
% the colored noise
% .........................................
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