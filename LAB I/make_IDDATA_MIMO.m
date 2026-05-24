function DATA = make_IDDATA_MIMO(y, u)

%
%  BEGIN
%
% Messages 
% ~~~~~~~~
        warning('off','MATLAB:dispatcher:InexactMatch') ; 
	FN = '<MAKE_DATA>: ' ; 
	E1 = [FN 'Missing or empty input data. Empty output. Exit.'] ; 
	BL = [blanks(3) '* Insert '] ; 
	EN = ' (ENTER means none): ' ;
	I1 = [BL 'the data name block [ENTER means ''DATA'']: '] ;
	I2 = [BL 'data notes' EN] ; 
	I3 = [BL 'the experiment name' EN] ; 
	I4 = [BL 'the time unit' EN] ; 
	I5 = [BL 'the starting date in format <dd-mmm-yyyy HH:MM:SS> (ENTER means NOW): '] ; 
	I6 = [BL 'user info (such as the starting date) as a string: '] ; 
	I7 = [BL 'data name on channel %d' EN] ; 
	I8 = [BL 'unit on channel %d' EN] ; 
	S  = [FN 'Data saved in file <%s.MAT>.'] ;
%
% Faults preventing 
% ~~~~~~~~~~~~~~~~~
	DATA = iddata ; 
	if (nargin < 2) % <2 pentru ca nargin verifica cate argumente primeste,inainte primea doar unul si acum primeste 2
	   war_err(E1) ; 
	   return ; 
	end  
	if (isempty(y))
	   war_err(E1) ; 
	   return ; 
	end 
%
% Building the DATA object
% ~~~~~~~~~~~~~~~~~~~~~~~~
%
	if (isscalar(y) || isvector(y))		% Storing the data ...
	   DATA.y = vectorize(y).' ; 
	   DATA.Ts = 1 ; 			% the sampling period ...
	else
	   DATA.y = y(:,2:end) ;
	   DATA.SamplingInstants = y(:,1) ; 	% and the sampling instants (if any). 
    end 
    
    % adaugare matrice de intrari u (specific pentru modelul mimo)
    if(~isempty(u))
        DATA.u = u;
    end

	war_err(FN) ;
	DATA.Name = input(I1,'s') ; 		% Setting the name of data block. 
	if (isempty(DATA.Name))
	   DATA.Name = 'DATA' ; 
	end 
	DATA.Notes = input(I2,'s') ; 		% Setting the notes on data (what they mean).
	DATA.ExperimentName = {input(I3,'s')} ;	% Setting the name of experiment or supplementary information. 
	DATA.TimeUnit = input(I4,'s') ; 	% Setting the time unit (e.g. ms, s, hours, days, etc.)
	if (~isempty(DATA.Ts)) 			% Setting the starting date and/or time 
	   FN = input(I5,'s') ;			% (only allowed for uniform sampling). 
	   if (isempty(FN))
	      DATA.Tstart = now ; 
	   else
	      DATA.Tstart = datenum(FN) ; 
	   end 
	else					% Here the starting date can be specified as a string 
	   DATA.UserData = input(I6,'s') ; 	% in a preferred format (such as 'dd-Mmm-yyyy'). 
	end  
	EN = size(DATA.y,2) ; 
	FN = input(sprintf(I7,1),'s') ; 	% Setting the name of each output channel. 
	BL = input(sprintf(I8,1),'s') ; 	% Setting the unit of each output channel. 
	for (n=2:EN)
	   FN = [FN ; {input(sprintf(I7,n),'s')}] ; 
	   BL = [BL ; {input(sprintf(I8,n),'s')}] ; 
	end 
	DATA.OutputName = FN ; 
	DATA.OutputUnit = BL ; 

    EN = size(DATA.u,2) ; % aici am luat ce a fost folosit pentru output si am folosit pentru input
	FN = input(sprintf(I7,1),'s') ; 	% Setting the name of each input channel. 
	BL = input(sprintf(I8,1),'s') ; 	% Setting the unit of each input channel. 
	for (n=2:EN)
	   FN = [FN ; {input(sprintf(I7,n),'s')}] ; 
	   BL = [BL ; {input(sprintf(I8,n),'s')}] ; 
	end 
	DATA.InputName = FN ; 
	DATA.InputUnit = BL ; 
%	DATA.Domain = 'Time' ; 			% Data are in time domain. 
	Y = DATA ; 
	eval(['save ' DATA.Name '.mat Y']) ;	% Save the DATA object. 
	war_err(sprintf(S,DATA.Name)) ; 
%
%  END
%