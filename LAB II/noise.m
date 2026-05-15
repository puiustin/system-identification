function noise(operation)
%
% NOISE		Module illustrating stochastic processes. 
%
% Inputs:	operation
%                 # character string in range: 
%                        'close_noise'
%                        'close_noise_def'
%                        'init_noise'
%                        'move_p'
%                        'move_z'
%                        'moved_p'
%                        'moved_z'
%                        'moving_p'
%                        'moving_z'
%                        'noiseclear'
%                        'show'               (by default)
%                        'system'
%                        'winit_noise'
%
% Outputs:	-----------
%
% Author:   Helena Haglund (*)
% Revised:  Bjorn Wittenmark (*)
%           Dan Stefanoiu (**)
%           Lavinius Ioan Gliga (**)
%
% Last upgrade: (*)  January 3, 1997
%               (**) March   8, 2004
%                    February 26, 2018
%
% Copyright: (*)  Lund Institute of Technology, SWEDEN
%                 Department of Automatic Control
%            (**) "Politehnica" Unversity of Bucharest, ROMANIA
%                 Department of Automatic Control & Computer Science
%
%

%
% BEGIN
%
% Global variables
% ~~~~~~~~~~~~~~
global fig_noise ccs_col fig_ccs
global system_noise error_noise
global Bd Ad h
global pole zero x y
global disc_axes_noise
global tau t
global ry C e
global a1 a2 b1 Bd1 K
global cov_handle spectr_handle realiz_handle
global cov_axes spectrum_axes realization_axes
%
% Initial constants
% ~~~~~~~~~~~~~~
ccs_col=1;
%
% Faults preventing
% ~~~~~~~~~~~~~~~
if (nargin<1)
   operation = 'show'; 
end
if (isempty(operation))
   operation = 'show'; 
end
%
% Operation SHOW
% ~~~~~~~~~~~~~~
%
% - checks if window already exists
%
if strcmp(operation,'show')
   existFlag =  fig_exist('Noise');
   if ~existFlag
      noise('winit_noise');
      noise('init_noise');	
   else
      clf;
      noise('init_noise');
   end
%
% Operation SYSTEM
% ~~~~~~~~~~~~~~~~
%
% - otherwise, draw the window 
%
elseif strcmp(operation,'system') 
   watchon;
   figure(fig_noise);
   set(error_noise,'Visible','off');
   h = 1;
   %
   % - make plots go clear after next updating
   %
   %set(cov_handle,'EraseMode','XOR');
   %set(spectr_handle,'EraseMode','XOR');
   %set(realiz_handle,'EraseMode','XOR');
   if get(system_noise,'value')==1
      axes(disc_axes_noise);
      cla;
      %
      % - plot unit circle
      %
      t=0:.1:6.3;			
      plot(sin(t),cos(t),'k-');
      axes(cov_axes);
      cla;
      axes(spectrum_axes);
      cla;
      axes(realization_axes);
      cla;
   elseif get(system_noise,'value')==2
      Ad = [1 -0.5];
      a1 = -0.5;
      %
      % - gives var(y)=1
      %
      Bd1 = 1;
      Bd = sqrt(1-a1^2);
      [phi,gam,C,~] = tf2ss(Bd,Ad);     
   elseif get(system_noise,'Value')==3
      Ad = [1 -0.3 0.1];
      Bd1 = [1 -0.5];
      a1 = -0.3;
      a2 = 0.1;
      b1 = -0.5;
      %
      % - gives var(y)=1
      %
      vary = ((1+b1^2)*(1+a2)-2*b1*a1)/...
             ((1-a2^2)*(1+a2)-(a1-a1*a2)*a1);
      K = 1/sqrt(vary);
      Bd = K*Bd1;
      [phi,gam,C,~] = tf2ss(Bd,Ad); 
   end
   if get(system_noise,'Value')==2 || ... 
      get(system_noise,'Value')==3 || ...
      get(system_noise,'Value')==4
      axes(disc_axes_noise);
      cla;
      %
      % - plot unit circle
      %
      t=0:.1:6.3;			
      plot(sin(t),cos(t),'k-');
      %
      % - plot poles and zeros of sampled system
      %
      if ccs_col==1
         pole = plot(real(roots(Ad)),imag(roots(Ad)),'rx');
         set(pole,'Linewidth',2,  ...
                  'Markersize',10);
         zero = plot(real(roots(Bd)),imag(roots(Bd)),'ro');
         set(zero,'Linewidth',2, ...
                  'Markersize',7);
      else
         pole = plot(real(roots(Ad)),imag(roots(Ad)),'kx');
         set(pole,'Linewidth',2,  ...
                  'Markersize',10);
         zero = plot(real(roots(Bd)),imag(roots(Bd)),'ko');
         set(zero,'Linewidth',2, ...
                  'Markersize',7);
      end
      %
      % - makes poles movable
      %
      set(pole,'ButtonDownFcn','noise(''move_p'')');
      set(zero,'ButtonDownFcn','noise(''move_z'')');
      %
      % - calculates covariance function
      %
      ry = zeros(20,1);
      tau = 20;
      rx = dlyap(phi,gam*gam');
      ry(1) = C*rx*C';
      for k=1:tau
         rx = phi*rx;
         ry(k+1) = C*rx*C'; 
      end
      axes(cov_axes);
      %
      % - plot covariance function
      %
      set(cov_handle,'XData',0:tau,'YData',ry);
      set(cov_handle,'LineWidth',2 );
      %
      % - calculate spectrum
      %
      [w_ab,fi_ab]=d_spektr(Ad,Bd,1);

      %
      % - plot spectrum
      %
      axes(spectrum_axes);
      set(spectr_handle,'XData',w_ab,'YData',fi_ab);
      set(spectr_handle,'LineWidth',2 );
      %
      % - calculate realization
      %
      axes(realization_axes);
      %
      % - white noise			
      %
      e = randn(1,50);
      y = filter(Bd,Ad,e);
      t = 0:49;
      %
      % - plot realization
      %
      set(realiz_handle,'XData',t,'YData',y);
      set(realiz_handle,'LineWidth',2 );
   end
   axes(disc_axes_noise);
   figure(fig_noise);
   watchoff;
%
% Operation MOVE_P, MOVED_P, MOVING_P (poles)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif strcmp(operation,'move_p')
   set(fig_noise,'WindowButtonMotionFcn', ...
                 'noise(''moving_p'');', ...
                 'WindowButtonUpFcn', ...
                 'noise(''moved_p'');');
elseif strcmp(operation,'moving_p')
   currpoint = get(disc_axes_noise,'CurrentPoint');
   x = currpoint(1,1);
   y = currpoint(1,2);
   if get(system_noise,'Value')==2
      if ccs_col==1
         set(pole,'XData',x,'YData',0,'Color','r');
      else
         set(pole,'XData',x,'YData',0,'Color','k');
      end
   elseif get(system_noise,'Value')==3
      if ccs_col==1
         set(pole,'XData',[x x],'YData',[y -y],'Color','r');
      else
         set(pole,'XData',[x x],'YData',[y -y],'Color','k');
      end
   end
elseif strcmp(operation,'moved_p')  	
   %
   % - update plots after completed move
   %
   if get(system_noise,'Value')==2
      if ccs_col==1
         set(pole,'XData',x,'YData',0,'Color','r');
      else
         set(pole,'XData',x,'YData',0,'Color','k');
      end
      Ad = [1 -x];
      a1 = -x;
      %
      % - test stability
      %
      if (abs(x)>=1) 
         Bd_test = 1; 
         Bd = Bd1;
      else
         %
         % - gives var(y)=1
         %
         Bd_test = 0; 
         Bd=sqrt(3)/2;
      end
   elseif get(system_noise,'Value')==3
      if ccs_col==1
         set(pole,'XData',[x x],'YData',[y -y],'Color','r');
      else
         set(pole,'XData',[x x],'YData',[y -y],'Color','k');
      end
      Ad = [1 -2*x x^2+y^2];
      a1 = -2*x;
      a2 = x^2+y^2;
      %
      % - test stability
      %
      if (max(abs(roots(Ad)))>=1)
         Bd_test = 1; 
         Bd = Bd1;
      else
         Bd_test = 0; 
         %
         % - gives var(y)=1
         %
         vary = ((1+b1^2)*(1+a2)-2*b1*a1)/ ...
                ((1-a2^2)*(1+a2)-(a1-a1*a2)*a1);
         K = 1/sqrt(vary);
         Bd = K*Bd1;
      end
   end
   %
   % - warning message if system is unstable
   %
   if (Bd_test==1 && get(system_noise,'Value')==2) || ...
      (Bd_test==1 && get(system_noise,'Value')==3)
      set(error_noise,'Visible','on');
   else
      set(error_noise,'Visible','off');		
   end
   [phi,gam,C,~] = tf2ss(Bd,Ad);
   set(fig_noise,'WindowButtonMotionFcn','', ...
                 'WindowButtonUpFcn','');
   %
   % - covariance function and spectrum if system is stable
   %
   if Bd_test==0
      %
      % - calculate covariance function
      %
      rx = dlyap(phi,gam*gam');
      ry(1) = C*rx*C';
      for k=1:tau
         rx = phi*rx;
         ry(k+1) = C*rx*C'; 
      end
      %
      % - plot covariance function
      %
      axes(cov_axes);
      set(cov_handle,'XData',0:tau,'YData',ry);
      %
      % - calculate spectrum
      %
      [w_ab,fi_ab]=d_spektr(Ad,Bd,1);

      %      
      % - plot spectrum
      %
      axes(spectrum_axes);
      set(spectr_handle,'XData',w_ab,'YData',fi_ab);
      set(spectr_handle,'LineWidth',2 );
   end
   %
   % - calculate and plot realization
   %
   axes(realization_axes);
   y = filter(Bd,Ad,e);
   t=0:49;
   set(realiz_handle,'XData',t,'YData',y);
   axes(disc_axes_noise);
   figure(fig_noise);
   watchoff;
%
% Operation MOVE_Z, MOVED_Z, MOVING_Z (zeros)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif strcmp(operation,'move_z')
   set(fig_noise,'WindowButtonMotionFcn', ...
                 'noise(''moving_z'');', ...
                 'WindowButtonUpFcn', ...
                 'noise(''moved_z'');');
elseif strcmp(operation,'moving_z')
   currpoint = get(disc_axes_noise,'CurrentPoint');
   x = currpoint(1,1);
   y = currpoint(1,2);
   set(zero,'XData',x,'YData',0);
elseif strcmp(operation,'moved_z')  	
   %
   % - update plots after completed move
   %
   watchon;
   set(zero,'XData',x,'YData',0);
   Bd1 = [1 -x];
   b1 = -x;
   %
   % - test stability
   %
   if (max(abs(roots(Ad)))>=1)
      Bd_test = 1; 
      Bd = Bd1;
   else
      Bd_test = 0;
      %
      % - gives var(y)=1
      %
      vary = ((1+b1^2)*(1+a2)-2*b1*a1)/ ... 
             ((1-a2^2)*(1+a2)-(a1-a1*a2)*a1);
      K = 1/sqrt(vary);
      Bd = K*Bd1;
   end
   [phi,gam,C,~] = tf2ss(Bd,Ad);
   set(fig_noise,'WindowButtonMotionFcn','', ...
                 'WindowButtonUpFcn','');
   %
   % - covariance function and spectrum if system is stable
   %
   if Bd_test==0
      %
      % - calculate covariance function
      %
      rx = dlyap(phi,gam*gam');
      ry(1) = C*rx*C';
      for k=1:tau
         rx = phi*rx;
         ry(k+1) = C*rx*C'; 
      end
      %
      % - plot covariance function
      %
      axes(cov_axes);
      set(cov_handle,'XData',0:tau,'YData',ry);
      %
      % - calculate spectrum
      %
      [w_ab,fi_ab]=d_spektr(Ad,Bd,1);
                
      %
      % - plot spectrum
      %
      axes(spectrum_axes);
      set(spectr_handle,'XData',w_ab,'YData',fi_ab);
      set(spectr_handle,'LineWidth',2 );      
   end
   %
   % - calculate and plot realization
   %
   axes(realization_axes);
   y = filter(Bd,Ad,e);
   t=0:49;
   set(realiz_handle,'XData',t,'YData',y);
   axes(disc_axes_noise);
   figure(fig_noise);
   watchoff;
%
% Operation NOISECLEAR
% ~~~~~~~~~~~~~~~~~~~~
elseif strcmp(operation,'noiseclear')
   if get(system_noise,'Value')==2 || ... 
      get(system_noise,'Value')==3 || ...
      get(system_noise,'Value')==4
      %
      % - plot the most recent covariance function
      %
      axes(cov_axes);
      cla;
      if ccs_col==1
	 cov_handle = plot(0:tau, ry,'r');
      else
	 cov_handle = plot(0:tau, ry,'k');
      end
      set(cov_handle,'LineWidth',2 );
      %
      % - plot the most recent spectrum
      %
      axes(spectrum_axes);
      cla;
      if ccs_col==1
         spectr_handle = loglog(w_ab,fi_ab,'r');
      else
	 spectr_handle = loglog(w_ab,fi_ab,'k');
      end
      set(spectr_handle,'LineWidth',2 );
      %
      % - plot the most recent realization
      %
      axes(realization_axes);
      cla;
      if ccs_col==1				
         realiz_handle = plot(t,y,'r');
      else
         realiz_handle = plot(t,y,'k');
      end
      set(realiz_handle,'LineWidth',2 );
      figure(fig_noise);
      axes(disc_axes_noise);
   end
%
% Operation WINIT_NOISE
% ~~~~~~~~~~~~~~~~~~~~~
elseif strcmp(operation,'winit_noise')
   %
   % - create main window
   %
   fig_noise = figure('Name','Noise','NumberTitle','off', ...
	              'Units','Normalized','Position', ...
                       [0.2561 0.4400 0.4861 0.4667 ], ... 
                      'BackingStore','Off',...
                      'DefaultUicontrolFontSize',11);
   set(fig_noise,'Color',[0.8 0.8 0.8]);
%
% Operation INIT_NOISE
% ~~~~~~~~~~~~~~~~~~~
elseif strcmp(operation,'init_noise')
   watchon;
   figure(fig_noise);
   %
   % - frame left
   %
   close_noise = uicontrol(fig_noise,'Style','Push', ... 
                                     'String','Quit', ...
                                     'Units','Normalized', ...
                                     'Position', ...
                                      [0.0339 0.0690 0.1429 0.0595 ], ...
                                     'BackgroundColor',[1 0.4 0.4], ...
                                     'Callback', ... 
                                     'noise(''close_noise_def'');');
   %
   % - frame middle
   %
   frame_middle = uicontrol(fig_noise,'Style','Frame', ...
                                      'Units','Normalized', ...
                                      'Position', ... 
                                       [0.2036 0.7119 0.3214 0.2619 ]);
   system_noise = uicontrol(fig_noise,'Style','popup', ...
                                      'Units','Normalized', ...
                                      'Position', ...
                                       [0.2214 0.8667 0.2857 0.0595 ], ...
                                      'string',...
                                       ['Select system | ' ...
                                        'b/(z+a) | ' ... 
                                        '(b0z+b1)/(z^2+a1*z+a2)']);
   set(system_noise,'Callback','noise(''system'');');
   %
   % - create pole/zero axes
   %
   disc_axes_noise = axes('position',[0.23 0.15 0.28 0.28]);
   %
   % - plot unit circle
   %
   t=0:.1:6.3;
   plot(sin(t),cos(t),'k-');
   grid on;
   axis('equal');
   title('Poles/Zeros','Color','k', ...
	 'FontName','Times','Fontsize',11);
   set(disc_axes_noise,'XLim',[-1.05 1.05],'YLim',[-1.05 1.05],...
                       'Clipping','Off','XLimMode','Manual', ... 
                       'YLimMode','Manual','YTick',[-1 0 1], ... 
                       'SortMethod','childorder','Xcolor','k', ... 
                       'Ycolor','k','FontName','Times','Fontsize',11);
   hold on;
   %
   % - create covariance diagram
   %
   cov_axes = axes('Position',[0.6 0.71 0.35 0.22]);
   grid on;
   set(cov_axes, 'XLim',[0 20],'YLim', [-1 2], 'SortMethod', ... 
                 'childorder','Clipping','Off','XLimMode','Manual', ... 
                 'YLimMode','Manual','XColor','k','YColor','k', ...
                 'FontName','Times','Fontsize',11);
   title('Covariance function','Color','k', ...
         'FontName','Times','Fontsize',11);
   hold on;
   %
   % - create spectrum diagram
   %
   spectrum_axes = axes('Position',[0.6 0.39 0.35 0.22]);
   grid on;
   set(spectrum_axes,'XLim',[0.01 4],'YLim',[0.01 5], ...
                     'XScale','log','YScale','log', ... 
                     'XColor','k','YColor','k', ...
	             'FontName','Times','Fontsize',11);
   title('Spectrum','Color','k',...
         'FontName','Times','Fontsize',11);	
   hold on;
   %
   % - create realization diagram
   %
   realization_axes = axes('Position',[0.6 0.05 0.35 0.22]);
   grid on;
   set(realization_axes,'XLim',[0 50],'YLim', [-5 5], ...
                        'SortMethod','childorder','Clipping','Off', ... 
                        'XLimMode','Manual','YLimMode','Manual', ...
                        'XColor','k','YColor','k', ...
                        'FontName','Times','Fontsize',11);
   title('Realization','Color','k', ...
         'FontName','Times','Fontsize',11);
   hold on;
   %
   % - create handles
   %
   axes(cov_axes);
   if ccs_col==1
      cov_handle = plot(NaN,NaN,'r');
      axes(spectrum_axes);
      spectr_handle = loglog(NaN,NaN,'r');
      axes(realization_axes);
      realiz_handle = plot(NaN,NaN,'r');
   else
      cov_handle = plot(NaN,NaN,'k');
      axes(spectrum_axes);
      spectr_handle = loglog(NaN,NaN,'k');
      axes(realization_axes);
      realiz_handle = plot(NaN,NaN,'k');
   end
   watchoff;
   %
   % - error mesage
   %
   error_noise = uicontrol(fig_noise,'Style','text',...
                                     'Units', 'Normalized', ...
                                     'Position', ...
                                      [0.2304 0.5202 0.2679 0.0950 ], ...
                                     'String','Unstable system!', ... 
                                     'Fontsize',12, ...
                                     'BackgroundColor','r');
   set(error_noise,'Visible','off');
%
% Operation CLOSE_NOISE
% ~~~~~~~~~~~~~~~~~~~~~
elseif strcmp(operation, 'close_noise')
   existFlag =  fig_exist('Noise');
   if existFlag
      close(fig_noise);	
   end
%
% Operation CLOSE_NOISE_DEF
% ~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif strcmp(operation, 'close_noise_def')
   existFlag =  fig_exist('Noise');
   if existFlag
      close(fig_noise);	
   end
%
% Operation unknown
% ~~~~~~~~~~~~~~~~~
else
    error(['### Unknown argument <<' operation '>>. Nothing to do.'])  ;
end
%
% END
%