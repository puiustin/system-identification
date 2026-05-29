function [F,D,M] = ISLAB_12A(mt,K0,T0,Tmax,Ts,U,lambda) 
%
% ISLAB_12A   Module that performs off-line identification of 
%             physical parameters of a DC engine (gain and 
%             time constant). The parameters are CONSTANT here. 
%
% Inputs:	mt     # model type: 
%                         0 -> OE (default)
%                         1 -> ARX
%               K0     # constant gain (4, by default)
%               T0     # constant time constant (0.5 s, by default)
%               Tmax   # simulation duration 
%                        (80 s, by default)
%               Ts     # sampling period
%                        (0.1 s, by default) 
%               U      # amplitude of input square wave 
%                        (0.5, by default) 
%               lambda # standard deviation of white noise 
%                        (1, by default)
%
% Outputs:      F      # structure representing the estimated 
%                        physical parameters: 
%                          F.K  -> gain
%                          F.T  -> time constant
%               D      # IDDATA object representing the I/O data 
%                        employed in identification
%               M      # IDMODEL object representing the estimated 
%                        discrete time model 
%
% Explanation:	A second order continuous transfer function 
%               with a null pole, gain K and time constant T 
%               (the model of a DC engine) is stimulated with 
%               a square wave in order to provide identification 
%               data. Constant parameters K and T are identified 
%               by discretizing the transfer function.  
%               (See the function GDATA_DCENG.)
%
% Author:   Dan Stefanoiu (*)
% Revised:  Dan Stefanoiu (*)
%           Lavinius Ioan Gliga (*)
%
% Created: April 29, 2004
% Revised: January 30, 2012
%          August 9, 2018
%
% Copyright: (*) "Politehnica" University of Bucharest, ROMANIA
%                Department of Automatic Control & Computer Science
%

%
% BEGIN
% 

global FIG ;			% Figure number handler 
FIG = 1;            

% 
% Messages
% ~~~~~~~~
FN = '<ISLAB_12A>: ' ; 
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
if (isempty(Ts))
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
   mt = 0 ;
end
if (isempty(mt))
   mt = 0 ;
end 
mt = abs(round(mt(1))) ; 
% 
% Generate identification data 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[D,V] = gdata_DCeng(0,K0,T0,Tmax,Ts,U,lambda) ; 
% 
% Estimate discrete time parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (~mt)
   M = oe(D,[2 2 1]) ; 	        % Model of tipe OE. 
else
   M = arx(D,[2 2 1]) ;         % Model of type ARX. 
end 
% 
% Estimate physical parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (~mt)                        % OE model. 
   F.K = sum(M.b)/(1-M.f(3))/Ts ; 
   F.T = Ts*(M.f(3)*M.b(2)+M.b(3))/sum(M.b)/(1-M.f(3)) ; 
%   F.T = -Ts/log(M.f(3)) ; 
else                            % ARX model.
   F.K = sum(M.b)/(1-M.a(3))/Ts ; 
   F.T = Ts*(M.a(3)*M.b(2)+M.b(3))/sum(M.b)/(1-M.a(3)) ; 
%   F.T = -Ts/log(M.a(3)) ; 
end 
% 
% Display physical parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
war_err(M1) ; 
disp(M2) ; 
fprintf(1, M3, K0, F.K);
fprintf(1, M4, T0, F.T);
war_err(PK) ; 
pause ;
% 
% Plot I/O data
% ~~~~~~~~~~~~~
Tmax = Ts*round(Tmax/Ts) ; 	% Correct Tmax. 
t = 0:Ts:Tmax ; 		% Set time axis. 
figure(FIG),clf
   fig_look(FIG,1.5) ; 
   plot(t,D.u,'-b',t,D.y,'-r') ; 
   FN = scaling([D.u D.y]) ;     % Re-scale the axes. 
   axis([0 Tmax FN]) ; 
   title('Input-output data provided by a DC engine.') ; 
   xlabel('Time [s]') ; 
   ylabel('Magnitude') ; 
   set(FIG,'DefaultTextHorizontalAlignment','left') ; 
   legend('input (square wave)','output') ; 
FIG = FIG+1 ;
% 
% Plot I/O simulated data
% ~~~~~~~~~~~~~~~~~~~~~~~
V.y = sim(M,[D.u zeros(size(D.u))]) ; 
figure(FIG),clf
   fig_look(FIG,1.5) ; 
   plot(t,V.y,'-b',t,D.y,'-r',t,V.y,'-b') ; 
   FN = scaling([V.y D.y]) ;     % Re-scale the axes. 
   axis([0 Tmax FN]) ; 
   title(['Output data provided by a DC engine ' ... 
          'and its discrete model.']) ; 
   xlabel('Time [s]') ; 
   ylabel('Magnitude') ; 
   set(FIG,'DefaultTextHorizontalAlignment','left') ; 
   legend('simulated output','measured output') ; 
FIG = FIG+1 ;
%
% END
%