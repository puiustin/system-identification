function [F,D,M,MODEL_BEST] = ISLAB_12B(mt,K0,T0,Tmax,Ts,U,lambda) 
%
% BEGIN
%

% Messages
FN = '<ISLAB_12B>: ' ;
NFN = length(FN) ;
PK = [blanks(70) '<Press a key>'] ;
M1 = [FN 'Physical parameters:'] ;
M2 = [blanks(NFN+7) 'True' blanks(8) 'Estimated'] ;
M3 = [blanks(NFN) 'K: %8.4f' blanks(6) '%8.4f'] ;
M4 = [blanks(NFN) 'T: %8.4f' blanks(6) '%8.4f'] ;

% Faults preventing
if (nargin < 7), lambda = 1 ; end
if (isempty(lambda)), lambda = 1 ; end
lambda = abs(lambda(1)) ;
if (~lambda), lambda = 1 ; end

if (nargin < 6), U = 0.5 ; end
if (isempty(U)), U = 0.5 ; end
U = U(1) ;
if (abs(U) < eps), U = 0.5 ; end

if (nargin < 5), Ts = 0.1 ; end
if (isempty(Ts)), Ts = 0.1 ; end
Ts = abs(Ts(1)) ;
if (Ts < eps), Ts = 0.1 ; end

if (nargin < 4), Tmax = 80 ; end
if (isempty(Tmax)), Tmax = 80 ; end
Tmax = abs(Tmax(1)) ;
if (Tmax < eps), Tmax = 80 ; end
Tmax = max(250*Ts, Tmax) ;

if (nargin < 3), T0 = 0.5 ; end
if (isempty(T0)), T0 = 0.5 ; end
T0 = T0(1) ;
if (abs(T0) < eps), T0 = 0.5 ; end

if (nargin < 2), K0 = 4 ; end
if (isempty(K0)), K0 = 4 ; end
K0 = K0(1) ;
if (abs(K0) < eps), K0 = 4 ; end

if (nargin < 1), mt = 0 ; end
if (isempty(mt)), mt = 0 ; end
mt = abs(round(mt(1))) ;

% Generate identification data
[D,V] = gdata_DCeng(0,K0,T0,Tmax,Ts,U,lambda) ;

% Estimate discrete time parameters
if (~mt)
   M = oe(D,[2 2 1]) ;          % OE model
else
   M = arx(D,[2 2 1]) ;         % ARX model
end

% Estimate physical parameters
% K = (b1+b2) / (Ts*(1-a2)),  T = Ts*(a2*b1+b2) / ((b1+b2)*(1-a2))
if (~mt)                        % OE model
   b1 = M.b(2) ; b2 = M.b(3) ;
   a2 = M.f(3) ;
else                            % ARX model
   b1 = M.b(2) ; b2 = M.b(3) ;
   a2 = M.a(3) ;
end
F.K = (b1 + b2) / (Ts * (1 - a2)) ;
F.T = Ts * (a2*b1 + b2) / ((b1 + b2) * (1 - a2)) ;

% Display physical parameters
disp(M1) ;
disp(M2) ;
disp(sprintf(M3, K0, F.K)) ;
disp(sprintf(M4, T0, F.T)) ;
disp(PK) ;
pause ;

% Time axis
Tmax = Ts*round(Tmax/Ts) ;
t = (0:Ts:Tmax)' ;
N = length(t) ;

% Simulated output
y_sim = sim(M, [D.u zeros(size(D.u))]) ;

% Colored noise: v = measured - simulated  (correct sign: noise corrupts measured)
v = D.y - y_sim ;

% Running variance of v: sigma2(n) = (1/n) * sum_{k=1}^{n} v(k)^2
sigma2_v = cumsum(v.^2) ./ (1:N)' ;

% ---- Figure 1: I/O data ----
FIG = 10 ;
figure(FIG), clf
   plot(t, D.u, '-b', t, D.y, '-r') ;
   title('Input-output data provided by a DC engine.') ;
   xlabel('Time [s]') ; ylabel('Magnitude') ;
   legend('input (square wave)', 'output', 0) ;
FIG = FIG + 1 ;

% ---- Figure 2: simulated vs measured + noise + variance ----
figure(FIG), clf
   subplot(3,1,1)
      plot(t, y_sim, '-b', t, D.y, '-r') ;
      title('Simulated output vs measured output.') ;
      xlabel('Time [s]') ; ylabel('Magnitude') ;
      legend('simulated output', 'measured output', 0) ;
   subplot(3,1,2)
      plot(t, v, '-k') ;
      title(sprintf('Colored noise v = y_{meas} - y_{sim}  (\\sigma^2 = %.4f)', var(v))) ;
      xlabel('Time [s]') ; ylabel('v') ;
   subplot(3,1,3)
      plot(t, sigma2_v, '-m') ;
      yline(var(v), '--r', sprintf('\\sigma^2 = %.4f', var(v))) ;
      title('Running variance \sigma^2(n) of colored noise v') ;
      xlabel('Time [s]') ; ylabel('\sigma^2') ;
FIG = FIG + 1 ;

disp(PK) ;
pause ;

% ================================================================
% Noise model identification
% 20 models: 10 ARMA(na,nc), 5 AR(na) via levinson, 5 MA(nc)
% Structural indices chosen pseudo-randomly in [1,20]
% ================================================================
disp(' ') ;
disp('Identificarea zgomotului colorat v') ;
disp('===================================') ;

ERR_DATA = iddata(v, [], Ts) ;   % noise-only iddata (no input)

% -- 10 ARMA models (armax with random na, nc in [1,20]) --
for i = 1:10
   na = randi(20) ; nc = randi(20) ;
   Mnoise(i) = armax(ERR_DATA, [na nc]) ;
   fprintf('ARMA model %2d: na=%d, nc=%d\n', i, na, nc) ;
end

% -- 5 AR models via levinson (pure AR: nc=0) --
% levinson(r, p) fits AR(p) from autocorrelation r
r_v = xcorr(v, 'biased') ;       % autocorrelation of v
r_v = r_v(N:end) ;               % keep non-negative lags
for i = 1:5
   p = randi(20) ;
   [a_lev, err_lev] = levinson(r_v, p) ;
   % Build idmodel from AR coefficients
   Mnoise(10+i) = armax(ERR_DATA, [p 0]) ;  % armax with nc=0 = pure AR
   fprintf('AR  model %2d: na=%d (levinson)\n', 10+i, p) ;
end

% -- 5 MA models (armax with na=0, nc random in [1,20]) --
for i = 1:5
   nc = randi(20) ;
   Mnoise(15+i) = armax(ERR_DATA, [0 nc]) ;
   fprintf('MA  model %2d: nc=%d\n', 15+i, nc) ;
end

% ================================================================
% Evaluate all 20 models
% ================================================================
best_lambda2 = Inf ;
best_idx = 1 ;

for i = 1:20
   % White noise residuals
   e = resid(ERR_DATA, Mnoise(i)) ;
   lambda2 = var(e.y) ;          % variance of white residual
   
   if lambda2 < best_lambda2
      best_lambda2 = lambda2 ;
      best_idx = i ;
   end
   
   fprintf('\n--- Model %2d | lambda^2 = %.6f ---\n', i, lambda2) ;
end

fprintf('\nAutomatic best model: M%d (lambda^2 = %.6f)\n', best_idx, best_lambda2) ;

% ---- Figure 3: variance comparison across models ----
figure(FIG), clf
   lambda2_all = zeros(1,20) ;
   for i = 1:20
      e = resid(ERR_DATA, Mnoise(i)) ;
      lambda2_all(i) = var(e.y) ;
   end
   bar(lambda2_all) ;
   hold on
   bar(best_idx, lambda2_all(best_idx), 'r') ;
   title('White noise variance \lambda^2 for each noise model') ;
   xlabel('Model index') ; ylabel('\lambda^2') ;
   legend('Models', 'Best') ;
FIG = FIG + 1 ;

disp(' ') ;
option = input(sprintf('Alegeti cel mai bun model (1->20) [automat: %d]: ', best_idx)) ;
if isempty(option), option = best_idx ; end
MODEL_BEST = Mnoise(option) ;

% ---- Figure 4: best model residuals ----
e_best = resid(ERR_DATA, MODEL_BEST) ;
sigma2_e = var(e_best.y) ;
figure(FIG), clf
   subplot(2,1,1)
      plot(t, v, '-k', t, filter(-MODEL_BEST.c, MODEL_BEST.a, v), '-r') ;
      title(sprintf('Colored noise v and model fit (model M%d)', option)) ;
      xlabel('Time [s]') ; ylabel('Amplitude') ;
      legend('v (colored)', 'model output', 0) ;
   subplot(2,1,2)
      plot(t, e_best.y, '-b') ;
      title(sprintf('White noise e (residual) | \\lambda^2 = %.6f', sigma2_e)) ;
      xlabel('Time [s]') ; ylabel('e') ;
FIG = FIG + 1 ;

disp(' ') ;
fprintf('Dispersia zgomotului alb rezultat: lambda^2 = %.6f\n', sigma2_e) ;
disp(PK) ;
pause ;

%
% END
%
