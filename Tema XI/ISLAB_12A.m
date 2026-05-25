function [F,D,M] = ISLAB_12A(mt,K0,T0,Tmax,Ts,U,lambda) 
%
% BEGIN
%

% Messages
FN = '<ISLAB_12A>: ' ; 
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

% Estimate physical parameters from discrete model coefficients
% Based on Euler discretization of H(s) = K / (s*(1+T*s))
% Hd(z) = (b1*z^-1 + b2*z^-2) / (1 + a1*z^-1 + a2*z^-2)
% K = (b1+b2) / (Ts*(1-a2))
% T = Ts * (a2*b1 + b2) / ((b1+b2)*(1-a2))
if (~mt)                        % OE model: theta = [f1 f2 b1 b2]
   b1 = M.b(2) ; b2 = M.b(3) ;
   a2 = M.f(3) ;                % f2 coefficient (= -a2 in standard form)
   F.K = (b1 + b2) / (Ts * (1 - a2)) ;
   F.T = Ts * (a2*b1 + b2) / ((b1 + b2) * (1 - a2)) ;
else                            % ARX model: theta = [a1 a2 b1 b2]
   b1 = M.b(2) ; b2 = M.b(3) ;
   a2 = M.a(3) ;
   F.K = (b1 + b2) / (Ts * (1 - a2)) ;
   F.T = Ts * (a2*b1 + b2) / ((b1 + b2) * (1 - a2)) ;
end

% Display physical parameters
disp(M1) ; 
disp(M2) ; 
disp(sprintf(M3, K0, F.K)) ; 
disp(sprintf(M4, T0, F.T)) ; 
disp(PK) ; 
pause ;

% Plot I/O data
Tmax = Ts*round(Tmax/Ts) ;
t = 0:Ts:Tmax ;
FIG = 10 ;
figure(FIG), clf
   plot(t, D.u, '-b', t, D.y, '-r') ; 
   FN2 = [min([D.u; D.y])-0.5, max([D.u; D.y])+0.5] ;
   axis([0 Tmax FN2]) ; 
   title('Input-output data provided by a DC engine.') ; 
   xlabel('Time [s]') ; 
   ylabel('Magnitude') ; 
   legend('input (square wave)', 'output', 0) ; 
FIG = FIG + 1 ;

% Plot simulated vs measured output
V.y = sim(M, [D.u zeros(size(D.u))]) ; 
figure(FIG), clf
   plot(t, V.y, '-b', t, D.y, '-r') ; 
   FN2 = [min([V.y; D.y])-0.5, max([V.y; D.y])+0.5] ;
   axis([0 Tmax FN2]) ; 
   title('Output data provided by a DC engine and its discrete model.') ; 
   xlabel('Time [s]') ; 
   ylabel('Magnitude') ; 
   legend('simulated output', 'measured output', 0) ; 
FIG = FIG + 1 ;

disp(PK) ;
pause ;

%
% END
%
