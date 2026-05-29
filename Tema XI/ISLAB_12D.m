function [F,ID,SD] = ISLAB_12D(mt,K0,T0,Tmax,Ts,U,lambda) 
%
% ISLAB_12D   Module that performs on-line identification of 
%             physical parameters of a DC engine (gain and 
%             time constant). The parameters are VARIABLE here. 
%
% Inputs:	mt     # model type: 
%                         0 -> OE (default)
%                         1 -> ARX
%               K0     # basic gain (4, by default)
%               T0     # basic time constant (0.5 s, by default)
%               Tmax   # simulagtion duration 
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
%                          F.K  -> gain variation
%                          F.T  -> time constant variation
%               ID     # IDDATA object representing the I/O data 
%                        employed in identification
%               SD     # IDDATA object representing the simulated 
%                        data provided by the discrete transfer 
%                        functions, within SD.y (output) and 
%                        SD.u (noise); the input is saved 
%                        within ID.u. 
%
% Explanation:	A second order continuous transfer function 
%               with a null pole, gain K and time constant T 
%               (the model of a DC engine) is stimulated with 
%               a square wave in order to provide identification 
%               data. Variable parameters K and T are 
%               identified and pursued by discretizing the 
%               transfer function. 
%               (See the function GDATA_DCENG.)
%
% Author:   Dan Stefanoiu (*)
% Revised:  Dan Stefanoiu (*)
%           Lavinius Ioan Gliga (*)
%
% Created: April 29, 2004
% Revized: January 31, 2012
%          August 9, 2018
%
% Copyright: (*) "Politehnica" University of Bucharest, ROMANIA
%                Department of Automatic Control & Computer Science
%

%
% BEGIN
% 

global FIG ;			% Figure number handler 
FIG = 1;                                % (to be set before running the routine). 
% 
% Constants
% ~~~~~~~~~
cv = 1 ; 			% Flag indicating the type of 
                                % physical parameters: 
                                %  cv = 0 -> constant parameters
                                %  cv = 1 -> variable parameters
% 
% Messages
% ~~~~~~~~
FN = '<ISLAB_12D>: ' ; 
WB = [FN 'Recursive estimation of parameters. ' ...
         'This may take a minute. Please wait ...'] ; 
WE = [blanks(length(FN)) '... Done.'] ;  
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
war_err(WB) ; 
[ID,V,P] = gdata_DCeng(cv,K0,T0,Tmax,Ts,U,lambda) ; 
if (~cv)
   FN = ones(1,length(ID.y)) ; 
   P.num{1} = K0*FN ; 
   P.den{1} = T0*FN ; 
end 
% 
Tmax = Ts*round(Tmax/Ts) ; 	% Correct Tmax. 
% Estimate and pursue discrete time parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (~mt)			% Model of type OE.
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEW STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v = version('-release'); % get release of Matlab
    if (strcmp(v, '2014b') || strcmp(v, '2015a')) 
        [theta,SD] = roe(ID,[2 2 1],'ff',0.999) ; 
        theta = theta(:,[3 4 1 2]) ; % theta = [B|F]->[F|B].
    else % if the release of Matlab is at least 2015b
        EstOE = recursiveOE([2 2 1]); % initialize the ARX estimator
        EstOE.ForgettingFactor = 0.999;  % choose the estimation method and 
                                   % forgetting factor
        NrSam = Tmax * 1 / Ts + 1;
        theta = zeros(NrSam, 4);
        SD = zeros(NrSam, 1);
        for i = 1 : NrSam   % for each entry in the dataset
            [A, B, SD(i)] = step(EstOE, ID.y(i), ID.u(i)); % identify the 
                                                         % parameters
            theta(i, 1:2) = A(2:3); % save the new parameters, to work with
            theta(i, 3:4) = B(2:3); % the rest of the program
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEW STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   theta = theta(:,[3 4 1 2]) ; % theta = [B|F]->[F|B].
else				% Model of type ARX. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEW STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v = version('-release'); % get release of Matlab
    if (strcmp(v, '2014b') || strcmp(v, '2015a')) 
        [theta,SD] = rarx(ID,[2 2 1],'ff',0.999) ;
				% theta = [A|B]. 
    else % if the release of Matlab is at least 2015b
        EstARX = recursiveARX([2 2 1]); % initialize the ARX estimator
        EstARX.ForgettingFactor = 0.999;  % choose the estimation method
                                   % and the forgetting factor
        NrSam = Tmax * 1 / Ts + 1; % number of samples
        theta = zeros(NrSam, 4);
        SD = zeros(NrSam, 1);
        for i = 1 : NrSam   % for each entry in the dataset
            [A, B, SD(i)] = step(EstARX, ID.y(i), ID.u(i)); % identify the 
                                                         % parameters
            theta(i, 1:2) = A(2:3); % save the new parameters, to work
            theta(i, 3:4) = B(2:3); % with the rest of the program
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEW STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end 
war_err(WE) ; 
SD = iddata(SD,V.u,Ts) ; 	% Pack simulation data. 
% 
% Estimate physical parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   F.K = (theta(:,3)+theta(:,4))./(1-theta(:,2))/Ts ; 
   F.T = Ts*(theta(:,2).*theta(:,3)+theta(:,4))./ ... 
         (theta(:,3)+theta(:,4))./(1-theta(:,2)) ; 
%   F.T = -Ts./log(theta(:,2)) ; 
% 
% Plot I/O data
% ~~~~~~~~~~~~~
t = 0:Ts:Tmax ; 		% Set time axis. 
figure(FIG),clf
   fig_look(FIG,1.5) ; 
   plot(t,ID.u,'-b',t,ID.y,'-r') ; 
   FN = scaling([ID.u ID.y]) ;	% Re-scale the axes. 
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
figure(FIG),clf
   fig_look(FIG,1.5) ; 
   plot(t,SD.y,'-b',t,ID.y,'-r',t,SD.y,'-b') ; 
   FN = scaling([SD.y ID.y]) ;	% Re-scale the axes. 
   axis([0 Tmax FN]) ; 
   title(['Output data provided by a DC engine ' ... 
          'and its discrete model.']) ; 
   xlabel('Time [s]') ; 
   ylabel('Magnitude') ; 
   set(FIG,'DefaultTextHorizontalAlignment','left') ; 
   legend('simulated output','measured output') ; 
FIG = FIG+1 ;
% 
% Plot physical parameters variation
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure(FIG),clf
   fig_look(FIG,1.5) ; 
   subplot(211)
      plot(t,F.K,'-b',t,P.num{1},'-r') ; 
      FN = scaling([F.K P.num{1}']) ;	% Re-scale the axes. 
      axis([0 Tmax FN]) ; 
      title(['Physical parameters variation ' ... 
             '(DC engine).']) ; 
      xlabel('Time [s]') ; 
      ylabel('Gain K') ; 
      set(FIG,'DefaultTextHorizontalAlignment','left') ; 
      legend('estimated','true') ; 
      set(FIG,'DefaultTextHorizontalAlignment','center') ; 
   subplot(212)
      plot(t,F.T,'-b',t,P.den{1},'-r') ; 
      FN = scaling([F.T P.den{1}']) ;	% Re-scale the axes. 
      axis([0 Tmax FN]) ; 
      ylabel('Time constant T [s]') ; 
      xlabel('Time [s]') ; 
FIG = FIG+1 ;
% 
% Zoom on physical parameters variation (steady-state)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%WE = (t>(0.1*Tmax/T0)) ; 
WE = (t>(0.5*Tmax)) ; 
t = t(WE) ;
figure(FIG),clf
   fig_look(FIG,1.5) ; 
   subplot(211)
      plot(t,F.K(WE),'-b',t,P.num{1}(WE),'-r') ; 
					% Re-scale the axes. 
      FN = scaling([F.K(WE) P.num{1}(WE)']) ;	
      axis([t(1) Tmax FN]) ; 
      title(['Physical parameters variation - steady-state ' ... 
             '(DC engine).']) ; 
      xlabel('Time [s]') ; 
      ylabel('Gain K') ; 
      set(FIG,'DefaultTextHorizontalAlignment','left') ; 
      legend('estimated','true') ; 
      set(FIG,'DefaultTextHorizontalAlignment','center') ; 
   subplot(212)
      plot(t,F.T(WE),'-b',t,P.den{1}(WE),'-r') ; 
					% Re-scale the axes. 
      FN = scaling([F.T(WE) P.den{1}(WE)']) ; 
      axis([t(1) Tmax FN]) ; 
      ylabel('Time constant T [s]') ; 
      xlabel('Time [s]') ; 
FIG = FIG+1 ;
%
% END
%