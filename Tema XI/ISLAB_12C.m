function [F,D,M,MODEL_BEST] = ISLAB_12C(mt,K0,T0,Tmax,Ts,U,lambda) 
%
% BEGIN
%

% Messages
FN_str = '<ISLAB_12C>: ' ;
NFN = length(FN_str) ;
PK = [blanks(70) '<Press a key>'] ;
M1 = [FN_str 'Physical parameters:'] ;
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

% State-space identification (order 2)
if (~mt)
   disp([FN_str 'Using n4sid for state-space identification...']) ;
   M = n4sid(D, 2) ;             % order 2 forced
else
   disp([FN_str 'Using pem for state-space identification...']) ;
   M0 = n4sid(D, 2) ;            % initial model for pem
   M  = pem(D, M0) ;
end

% ================================================================
% Extract physical parameters from discrete state-space matrices
% Discrete A (Euler): A_d = I + Ts*A_c
%   A_c = [theta1  0 ]   =>  A_d(1,1) = 1 + Ts*theta1
%         [alpha   0 ]       A_d(2,1) = Ts*alpha
%   B_c = [theta2]          B_d(1,1) = Ts*theta2
%         [0     ]
%   C_c = [0  beta]         C_d = C_c (unchanged by Euler)
%
%   theta1 = (A_d(1,1) - 1) / Ts
%   alpha  = A_d(2,1) / Ts
%   theta2 = B_d(1,1) / Ts
%   beta   = C_d(1,2)
%
%   T = -1/theta1
%   K = -(alpha * beta * theta2) / theta1
% ================================================================
Ad = M.A ;
Bd = M.B ;
Cd = M.C ;

theta1 = (Ad(1,1) - 1) / Ts ;
alpha  = Ad(2,1) / Ts ;
theta2 = Bd(1,1) / Ts ;
beta   = Cd(1,2) ;

if abs(theta1) < eps
   warning([FN_str 'theta1 ~ 0, T cannot be estimated reliably.']) ;
   F.T = Inf ;
else
   F.T = -1 / theta1 ;
end

F.K = -(alpha * beta * theta2) / theta1 ;

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
y_sim = sim(M, D) ;
if isobject(y_sim)
   y_sim = y_sim.y ;
end

% Colored noise v = measured - simulated
v = D.y - y_sim ;
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
      title('Simulated output (state-space) vs measured output.') ;
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
% Noise model identification: 10 ARMA + 5 AR + 5 MA = 20 models
% Structural indices pseudo-randomly in [1,20]
% Note: pem/n4sid provide endogenous noise matrix (K in state-space).
%       Exogenous noise v must still be identified separately via ARMA.
% ================================================================
disp(' ') ;
disp('Identificarea zgomotului colorat v (exogen)') ;
disp('=============================================') ;
disp('Nota: pem/n4sid furnizeaza zgomotul endogen (matricea K).') ;
disp('      Zgomotul exogen v este identificat separat cu ARMA.') ;
disp(' ') ;

ERR_DATA = iddata(v, [], Ts) ;

% 10 ARMA models
for i = 1:10
   na = randi(20) ; nc = randi(20) ;
   Mnoise(i) = armax(ERR_DATA, [na nc]) ;
   fprintf('ARMA model %2d: na=%d, nc=%d\n', i, na, nc) ;
end

% 5 AR models via levinson
r_v = xcorr(v, 'biased') ;
r_v = r_v(N:end) ;
for i = 1:5
   p = randi(20) ;
   [a_lev, ~] = levinson(r_v, p) ;
   Mnoise(10+i) = armax(ERR_DATA, [p 0]) ;
   fprintf('AR  model %2d: na=%d (levinson)\n', 10+i, p) ;
end

% 5 MA models
for i = 1:5
   nc = randi(20) ;
   Mnoise(15+i) = armax(ERR_DATA, [0 nc]) ;
   fprintf('MA  model %2d: nc=%d\n', 15+i, nc) ;
end

% Find best model (minimum white noise variance lambda^2)
best_lambda2 = Inf ;
best_idx = 1 ;
lambda2_all = zeros(1,20) ;

for i = 1:20
   e = resid(ERR_DATA, Mnoise(i)) ;
   lambda2_all(i) = var(e.y) ;
   if lambda2_all(i) < best_lambda2
      best_lambda2 = lambda2_all(i) ;
      best_idx = i ;
   end
   fprintf('Model %2d: lambda^2 = %.6f\n', i, lambda2_all(i)) ;
end

fprintf('\nAutomatic best model: M%d (lambda^2 = %.6f)\n', best_idx, best_lambda2) ;

% ---- Figure 3: lambda^2 comparison ----
figure(FIG), clf
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
figure(FIG), clf
   subplot(2,1,1)
      plot(t, v, '-k') ;
      title(sprintf('Colored noise v | best model: M%d', option)) ;
      xlabel('Time [s]') ; ylabel('v') ;
   subplot(2,1,2)
      plot(t, e_best.y, '-b') ;
      title(sprintf('White noise e | \\lambda^2 = %.6f', var(e_best.y))) ;
      xlabel('Time [s]') ; ylabel('e') ;
FIG = FIG + 1 ;

disp(PK) ;
pause ;

%
% END
%
