function [F,ID,SD,MODEL_BEST] = ISLAB_12F(alpha_c,beta_c,K0,T0,Tmax,Ts,U,lambda) 
%
% BEGIN
%

% Constants
cv = 1 ;    % variable parameters

% Messages
FN_str = '<ISLAB_12F>: ' ;
WB = [FN_str 'Recursive state-space estimation. This may take a minute. Please wait ...'] ;
WE = [blanks(length(FN_str)) '... Done.'] ;
PK = [blanks(70) '<Press a key>'] ;

% Faults preventing
if (nargin < 8), lambda = 1 ; end
if (isempty(lambda)), lambda = 1 ; end
lambda = abs(lambda(1)) ;
if (~lambda), lambda = 1 ; end

if (nargin < 7), U = 0.5 ; end
if (isempty(U)), U = 0.5 ; end
U = U(1) ;
if (abs(U) < eps), U = 0.5 ; end

if (nargin < 6), Ts = 0.1 ; end
if (isempty(Ts)), Ts = 0.1 ; end
Ts = abs(Ts(1)) ;
if (Ts < eps), Ts = 0.1 ; end

if (nargin < 5), Tmax = 80 ; end
if (isempty(Tmax)), Tmax = 80 ; end
Tmax = abs(Tmax(1)) ;
if (Tmax < eps), Tmax = 80 ; end
Tmax = max(250*Ts, Tmax) ;

if (nargin < 4), T0 = 0.5 ; end
if (isempty(T0)), T0 = 0.5 ; end
T0 = T0(1) ;
if (abs(T0) < eps), T0 = 0.5 ; end

if (nargin < 3), K0 = 4 ; end
if (isempty(K0)), K0 = 4 ; end
K0 = K0(1) ;
if (abs(K0) < eps), K0 = 4 ; end

if (nargin < 2), beta_c = 1 ; end
if (isempty(beta_c)), beta_c = 1 ; end
beta_c = beta_c(1) ;
if (abs(beta_c) < eps)
   warning([FN_str 'beta_c=0 would cause K=0. Reset to 1.']) ;
   beta_c = 1 ;
end

if (nargin < 1), alpha_c = 1 ; end
if (isempty(alpha_c)), alpha_c = 1 ; end
alpha_c = alpha_c(1) ;
if (abs(alpha_c) < eps)
   warning([FN_str 'alpha_c=0 would cause K=0. Reset to 1.']) ;
   alpha_c = 1 ;
end

% Generate identification data
disp(WB) ;
[ID,V,P] = gdata_DCeng(cv,K0,T0,Tmax,Ts,U,lambda) ;

if (~cv)
   P.num{1} = K0 * ones(1,length(ID.y)) ;
   P.den{1} = T0 * ones(1,length(ID.y)) ;
end

% ================================================================
% Recursive state-space identification using rpem
% Model order 2 (DC engine has 2 states)
% rpem structure for state-space: use n4sid as initial model,
% then rpem recursively updates.
%
% rpem([na nb nc nd nf nk]) for Box-Jenkins / OE-like:
%   For a state-space model of order 2 with 1 input, 1 output:
%   use rpem(ID, 2, 'ff', 0.999) with integer order argument.
% ================================================================
[theta, SD] = rpem(ID, 2, 'ff', 0.999) ;

disp(WE) ;
SD = iddata(SD, V.u, Ts) ;

% ================================================================
% Extract physical parameters from rpem output
% rpem with integer order returns a structure array of idss models.
% For each time step n, theta(n) is an idss model.
% We extract A(1,1), A(2,1), B(1,1) from each step.
%
% theta1[n] = (A_d[n](1,1) - 1) / Ts
% theta2[n] = B_d[n](1,1) / Ts
% T[n] = -1/theta1[n]
% K[n] = -(alpha_c * beta_c * theta2[n]) / theta1[n]
% ================================================================
N_steps = length(theta) ;
F.K = zeros(N_steps,1) ;
F.T = zeros(N_steps,1) ;

for n = 1:N_steps
   An = theta(n).A ;
   Bn = theta(n).B ;
   theta1_n = (An(1,1) - 1) / Ts ;
   theta2_n = Bn(1,1) / Ts ;
   if abs(theta1_n) < eps
      F.T(n) = F.T(max(1,n-1)) ;   % hold previous value
      F.K(n) = F.K(max(1,n-1)) ;
   else
      F.T(n) = -1 / theta1_n ;
      F.K(n) = -(alpha_c * beta_c * theta2_n) / theta1_n ;
   end
end

% Time axis
Tmax = Ts*round(Tmax/Ts) ;
t = (0:Ts:Tmax)' ;
N = length(t) ;

% ================================================================
% Colored noise signals
% ================================================================
v_y = ID.y - SD.y ;
v_K = P.num{1}' - F.K ;
v_T = P.den{1}' - F.T ;

sigma2_vy = cumsum(v_y.^2) ./ (1:N)' ;
sigma2_vK = cumsum(v_K.^2) ./ (1:N)' ;
sigma2_vT = cumsum(v_T.^2) ./ (1:N)' ;

% ---- Figure 1: I/O data ----
FIG = 10 ;
figure(FIG), clf
   plot(t, ID.u, '-b', t, ID.y, '-r') ;
   title('Input-output data provided by a DC engine.') ;
   xlabel('Time [s]') ; ylabel('Magnitude') ;
   legend('input (square wave)', 'output', 0) ;
FIG = FIG + 1 ;

% ---- Figure 2: simulated vs measured + v_y + sigma^2_vy ----
figure(FIG), clf
   subplot(3,1,1)
      plot(t, SD.y, '-b', t, ID.y, '-r') ;
      title('Simulated (state-space rpem) vs measured output.') ;
      xlabel('Time [s]') ; ylabel('Magnitude') ;
      legend('simulated output', 'measured output', 0) ;
   subplot(3,1,2)
      plot(t, v_y, '-k') ;
      title(sprintf('Colored noise v_y  (\\sigma^2 = %.4f)', var(v_y))) ;
      xlabel('Time [s]') ; ylabel('v_y') ;
   subplot(3,1,3)
      plot(t, sigma2_vy, '-m') ;
      yline(var(v_y), '--r', sprintf('\\sigma^2_v = %.4f', var(v_y))) ;
      title('Running variance \sigma^2(n) of v_y') ;
      xlabel('Time [s]') ; ylabel('\sigma^2') ;
FIG = FIG + 1 ;

disp(PK) ;
pause ;

% ---- Figure 3: parameter tracking + tracking errors ----
figure(FIG), clf
   subplot(2,2,1)
      plot(t, F.K, '-b', t, P.num{1}', '-r') ;
      title('Gain K: estimated vs true') ;
      xlabel('Time [s]') ; ylabel('Gain K') ;
      legend('estimated', 'true', 0) ;
   subplot(2,2,2)
      plot(t, F.T, '-b', t, P.den{1}', '-r') ;
      title('Time constant T: estimated vs true') ;
      xlabel('Time [s]') ; ylabel('T [s]') ;
      legend('estimated', 'true', 0) ;
   subplot(2,2,3)
      plot(t, v_K, '-k') ;
      title(sprintf('v_K = K_{true} - K_{est}  (\\sigma^2_K = %.4f)', var(v_K))) ;
      xlabel('Time [s]') ; ylabel('v_K') ;
   subplot(2,2,4)
      plot(t, v_T, '-k') ;
      title(sprintf('v_T = T_{true} - T_{est}  (\\sigma^2_T = %.4f)', var(v_T))) ;
      xlabel('Time [s]') ; ylabel('v_T') ;
FIG = FIG + 1 ;

% ---- Figure 4: steady-state zoom ----
WE_idx = (t > 0.5*Tmax) ;
t_ss = t(WE_idx) ;
figure(FIG), clf
   subplot(2,1,1)
      plot(t_ss, F.K(WE_idx), '-b', t_ss, P.num{1}(WE_idx)', '-r') ;
      title('Steady-state tracking: Gain K') ;
      xlabel('Time [s]') ; ylabel('Gain K') ;
      legend('estimated', 'true', 0) ;
   subplot(2,1,2)
      plot(t_ss, F.T(WE_idx), '-b', t_ss, P.den{1}(WE_idx)', '-r') ;
      xlabel('Time [s]') ; ylabel('Time constant T [s]') ;
      legend('estimated', 'true', 0) ;
FIG = FIG + 1 ;

disp(PK) ;
pause ;

% ================================================================
% Adaptive ARMA noise identification for v_y
% Note: rpem provides endogenous noise (K_ss matrix).
%       Exogenous v_y must be identified separately.
% 10 ARMA + 5 AR + 5 MA = 20 adaptive models (rarmax)
% ================================================================
disp(' ') ;
disp('Identificarea adaptiva a zgomotului colorat v_y (exogen)') ;
disp('Note: rpem ofera zgomot endogen (matricea K_ss).') ;
disp('      Zgomotul exogen v_y este identificat separat cu rarmax.') ;
disp(' ') ;

ERR_DATA = iddata(v_y, [], Ts) ;

adm = 'ff' ;
adg = 0.999 ;    % numeric forgetting factor

% 10 ARMA adaptive models
for i = 1:10
   na = randi(20) ; nc = randi(20) ;
   Mnoise(i) = rarmax(ERR_DATA, [na nc], adm, adg) ;
   fprintf('ARMA model %2d: na=%d, nc=%d\n', i, na, nc) ;
end

% 5 AR adaptive models
for i = 1:5
   na = randi(20) ;
   Mnoise(10+i) = rarmax(ERR_DATA, [na 0], adm, adg) ;
   fprintf('AR   model %2d: na=%d\n', 10+i, na) ;
end

% 5 MA adaptive models
for i = 1:5
   nc = randi(20) ;
   Mnoise(15+i) = rarmax(ERR_DATA, [0 nc], adm, adg) ;
   fprintf('MA   model %2d: nc=%d\n', 15+i, nc) ;
end

% Select best model
best_lambda2 = Inf ;
best_idx = 1 ;
lambda2_all = zeros(1,20) ;

for i = 1:20
   th = Mnoise(i) ;
   th_last = th(end,:) ;
   try
      lambda2_all(i) = var(v_y) / (1 + norm(th_last)^2) ;
   catch
      lambda2_all(i) = Inf ;
   end
   if lambda2_all(i) < best_lambda2
      best_lambda2 = lambda2_all(i) ;
      best_idx = i ;
   end
   fprintf('Model %2d: lambda^2 ~ %.6f\n', i, lambda2_all(i)) ;
end

fprintf('\nAutomatic best model: M%d\n', best_idx) ;

% ---- Figure 5: lambda^2 bar chart ----
figure(FIG), clf
   bar(lambda2_all) ;
   hold on
   bar(best_idx, lambda2_all(best_idx), 'r') ;
   title('White noise variance \lambda^2 per adaptive noise model') ;
   xlabel('Model index') ; ylabel('\lambda^2') ;
   legend('Models', 'Best') ;
FIG = FIG + 1 ;

disp(' ') ;
option = input(sprintf('Alegeti cel mai bun model (1->20) [automat: %d]: ', best_idx)) ;
if isempty(option), option = best_idx ; end
MODEL_BEST = Mnoise(option) ;

fprintf('Modelul ales: M%d\n', option) ;
disp(PK) ;
pause ;

%
% END
%
