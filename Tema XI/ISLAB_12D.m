function [F,ID,SD] = ISLAB_12D(mt,K0,T0,Tmax,Ts,U,lambda) 
%
% BEGIN
%

% Constants
cv = 1 ;    % variable parameters

% Messages
FN = '<ISLAB_12D>: ' ;
WB = [FN 'Recursive estimation of parameters. This may take a minute. Please wait ...'] ;
WE = [blanks(length(FN)) '... Done.'] ;
PK = [blanks(70) '<Press a key>'] ;

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
disp(WB) ;
[ID,V,P] = gdata_DCeng(cv,K0,T0,Tmax,Ts,U,lambda) ;

% If constant parameters, build time-invariant parameter vectors
if (~cv)
   P.num{1} = K0 * ones(1,length(ID.y)) ;
   P.den{1} = T0 * ones(1,length(ID.y)) ;
end

% Recursive identification with forgetting factor ff=0.999
if (~mt)
   % OE model: roe returns theta = [B | F], reorder to [F | B]
   [theta,SD] = roe(ID,[2 2 1],'ff',0.999) ;
   theta = theta(:,[3 4 1 2]) ;   % now: theta = [f1 f2 b1 b2]
   % Columns: 1=f1, 2=f2, 3=b1, 4=b2
   a2_col = 2 ; b1_col = 3 ; b2_col = 4 ;
else
   % ARX model: rarx returns theta = [A | B] = [a1 a2 b1 b2]
   [theta,SD] = rarx(ID,[2 2 1],'ff',0.999) ;
   % Columns: 1=a1, 2=a2, 3=b1, 4=b2
   a2_col = 2 ; b1_col = 3 ; b2_col = 4 ;
end

disp(WE) ;
SD = iddata(SD, V.u, Ts) ;

% Extract physical parameters at each time step
% K[n] = (b1[n]+b2[n]) / (Ts*(1-a2[n]))
% T[n] = Ts*(a2[n]*b1[n]+b2[n]) / ((b1[n]+b2[n])*(1-a2[n]))
a2 = theta(:,a2_col) ;
b1 = theta(:,b1_col) ;
b2 = theta(:,b2_col) ;

denom_K = Ts * (1 - a2) ;
denom_T = (b1 + b2) .* (1 - a2) ;

F.K = (b1 + b2) ./ denom_K ;
F.T = Ts * (a2.*b1 + b2) ./ denom_T ;

% Time axis
Tmax = Ts*round(Tmax/Ts) ;
t = (0:Ts:Tmax)' ;

% ---- Figure 1: I/O data ----
FIG = 10 ;
figure(FIG), clf
   plot(t, ID.u, '-b', t, ID.y, '-r') ;
   title('Input-output data provided by a DC engine.') ;
   xlabel('Time [s]') ; ylabel('Magnitude') ;
   legend('input (square wave)', 'output', 0) ;
FIG = FIG + 1 ;

% ---- Figure 2: simulated vs measured ----
figure(FIG), clf
   plot(t, SD.y, '-b', t, ID.y, '-r') ;
   title('Output: simulated vs measured (DC engine).') ;
   xlabel('Time [s]') ; ylabel('Magnitude') ;
   legend('simulated output', 'measured output', 0) ;
FIG = FIG + 1 ;

% ---- Figure 3: parameter tracking (full time) ----
figure(FIG), clf
   subplot(2,1,1)
      plot(t, F.K, '-b', t, P.num{1}', '-r') ;
      title('Physical parameters variation (DC engine).') ;
      xlabel('Time [s]') ; ylabel('Gain K') ;
      legend('estimated', 'true', 0) ;
   subplot(2,1,2)
      plot(t, F.T, '-b', t, P.den{1}', '-r') ;
      xlabel('Time [s]') ; ylabel('Time constant T [s]') ;
      legend('estimated', 'true', 0) ;
FIG = FIG + 1 ;

% ---- Figure 4: parameter tracking (steady-state zoom) ----
WE_idx = (t > 0.5*Tmax) ;
t_ss = t(WE_idx) ;
figure(FIG), clf
   subplot(2,1,1)
      plot(t_ss, F.K(WE_idx), '-b', t_ss, P.num{1}(WE_idx)', '-r') ;
      title('Physical parameters variation - steady-state (DC engine).') ;
      xlabel('Time [s]') ; ylabel('Gain K') ;
      legend('estimated', 'true', 0) ;
   subplot(2,1,2)
      plot(t_ss, F.T(WE_idx), '-b', t_ss, P.den{1}(WE_idx)', '-r') ;
      xlabel('Time [s]') ; ylabel('Time constant T [s]') ;
      legend('estimated', 'true', 0) ;
FIG = FIG + 1 ;

disp(PK) ;
pause ;

%
% END
%
