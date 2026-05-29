function [F,D,M] = ISLAB_12BC(mt,K0,T0,Tmax,Ts,U,lambda) 
%
% ISLAB_12A   Module that performs off-line identification of 
%             physical parameters of a DC engine (gain and 
%             time constant). The parameters are CONSTANT here. 


%              OE armax recursiv, AR levinson




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

%
% BEGIN
% 

global FIG ;			
FIG = 1;            

FN = '<ISLAB_12C>: ' ; 
NFN = length(FN) ; 
PK = [blanks(70) '<Press a key>'] ; 
M1 = [FN 'Physical parameters:'] ; 
M2 = [blanks(NFN+7) 'True' blanks(8) 'Estimated'] ; 
M3 = [blanks(NFN) 'K: %8.4f' blanks(6) '%8.4f \n'] ; 
M4 = [blanks(NFN) 'T: %8.4f' blanks(6) '%8.4f \n'] ; 

if (nargin < 7)
   lambda = 1 ;
end 
if (isempty(lambda))
   lambda = 1 ;
end 
lambda = abs(lambda(1)) ; 

if (nargin < 6)
   U = 0.5 ;
end

if (isempty(U))
   U = 0.5 ;
end 

if (nargin < 5)
   Ts = 0.1 ;
end 

if (isempty(Ts))
   Ts = 0.1 ;
end 

if (nargin < 4)
   Tmax = 80 ;
end 

if (isempty(Tmax))
   Tmax = 80 ;
end 

if (nargin < 3) 
   T0 = 0.5 ;
end 

if (isempty(T0))
   T0 = 0.5 ;
end 

if (nargin < 2) 
   K0 = 4 ;
end

if (isempty(K0))
   K0 = 4 ;
end  

if (nargin < 1)
   mt = 0 ;
end

if (isempty(mt))
   mt = 0 ;
end 

mt = abs(round(mt(1))) ; 

[D,V] = gdata_DCeng(0,K0,T0,Tmax,Ts,U,lambda) ; 

if (~mt)
   M = oe(D,[2 2 1]) ; 	
else
   M = arx(D,[2 2 1]) ; 
end 

if (~mt)
   F.K = sum(M.b)/(1-M.f(3))/Ts ; 
   F.T = Ts*(M.f(3)*M.b(2)+M.b(3))/sum(M.b)/(1-M.f(3)) ; 
else
   F.K = sum(M.b)/(1-M.a(3))/Ts ; 
   F.T = Ts*(M.a(3)*M.b(2)+M.b(3))/sum(M.b)/(1-M.a(3)) ; 
end 

war_err(M1) ; 
disp(M2) ; 
fprintf(1, M3, K0, F.K);
fprintf(1, M4, T0, F.T);
war_err(PK) ; 
pause ;

Tmax = Ts*round(Tmax/Ts) ; 	
t = 0:Ts:Tmax ; 		

figure(FIG),clf
   fig_look(FIG,1.5) ; 
   plot(t,D.u,'-b',t,D.y,'-r') ; 
   FN = scaling([D.u D.y]) ;     
   axis([0 Tmax FN]) ; 
   title('Input-output data provided by a DC engine.') ; 
   xlabel('Time [s]') ; 
   ylabel('Magnitude') ; 
   legend('input','output') ; 
FIG = FIG+1 ;

V.y = sim(M,[D.u zeros(size(D.u))]) ; 

figure(FIG),clf
   fig_look(FIG,1.5) ; 
   plot(t,V.y,'-b',t,D.y,'-r',t,V.y,'-b') ; 
   FN = scaling([V.y D.y]) ;     
   axis([0 Tmax FN]) ; 
   title('Measured output versus simulated output') ; 
   xlabel('Time [s]') ; 
   ylabel('Magnitude') ; 
   legend('simulated output','measured output') ; 
FIG = FIG+1 ;

% .........................................
% Colored noise
% .........................................

v = D.y - V.y ;
sigma2 = var(v) ;

figure(FIG),clf
   fig_look(FIG,1.5) ;
   plot(t,v,'-r') ;
   FN = scaling(v) ;
   axis([0 Tmax FN]) ;
   title(['Colored noise v, sigma^2 = ' num2str(sigma2)]) ;
   xlabel('Time [s]') ;
   ylabel('v') ;
   legend(['sigma^2 = ' num2str(sigma2)]) ;
FIG = FIG+1 ;

% .........................................
% Noise model identification
% .........................................

best_lambda2 = inf ;
best_e = v ;
best_vhat = v ;
best_name = 'none' ;

if (~mt)

   rng('shuffle') ;

   for k = 1:10

      na_n = randi([1 20]) ;
      nc_n = randi([1 20]) ;

      try

         NM = armax(iddata(v,[],Ts),[na_n nc_n]) ;

         e_n = resid(NM,iddata(v,[],Ts)) ;
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

else

   A_arx = M.a ;

   best_e = filter(A_arx,1,v) ;

   best_lambda2 = var(best_e) ;

   best_vhat = v - best_e ;

   best_name = 'AR model from polynomial A' ;

end

% .........................................
% Output with estimated colored noise
% .........................................

yhat_noise = V.y + best_vhat ;

figure(FIG),clf
   fig_look(FIG,1.5) ;

   plot(t,D.y,'-r',t,yhat_noise,'-b') ;

   FN = scaling([D.y yhat_noise]) ;

   axis([0 Tmax FN]) ;

   title(['Measured output and output with estimated noise']) ;

   xlabel('Time [s]') ;

   ylabel('Magnitude') ;

   legend('measured output',best_name) ;

FIG = FIG+1 ;

% .........................................
% Estimated white noise
% .........................................

figure(FIG),clf

   fig_look(FIG,1.5) ;

   plot(t,best_e,'-k') ;

   FN = scaling(best_e) ;

   axis([0 Tmax FN]) ;

   title(['Estimated white noise e, lambda^2 = ' num2str(best_lambda2)]) ;

   xlabel('Time [s]') ;

   ylabel('e') ;

   legend(['lambda^2 = ' num2str(best_lambda2)]) ;

FIG = FIG+1 ;

%
% END
%