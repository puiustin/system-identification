%
%	File MAKE_DATA.M
%
%	Function: MAKE_DATA
%
%	Synopsis: DATA = make_DATA(y) ; 
%
%	Creates the IDDATA object from acquired data y. The first column of 
%	y sould be the sampling instants vector. The remaining columns are then
%	understood as measured data on different channels. If y is a vector 
%	(row or column), then it is understood as measured data on one channel. 
%	Once DATA being created, it is saved on disk in the current directory, 
%	by means of MATLAB command SAVE. The DATA object can be loaded in MATLAB 
%	workspace by means of dual command LOAD. 
%
%	The user is also invited to input some auxiliary data (if any), 
%	such as: 
%
%	   DATA.Name           = name of data block; the same name is used for the  
%	                         file saved on disk; (char string); 
%	   DATA.Notes          = what the data set means (char string); 
%	   DATA.ExperimentName = name of the measuring experiment (char string); 
%	   DATA.TimeUnit       = unit of time (e.g. seconds, minutes, hours, days, etc.); 
%	                         (char string); 
%	   DATA.Tstart         = integer that encodes the starting date of measurements; 
%	                         can be created by: 
%	                                     datenum('dd-mmm-yyyy HH:MM:SS') ; 
% 	                         to retrieve the date, call:  
%	                                     datestr(DATA.Tstart) ; 
%	                         the user is invited to insert the string: 'dd-mmm-yyyy HH:MM:SS'; 
%	   DATA.Ts             = sampling period (scalar) ; 
%	   DATA.OutputName     = name of each measuring channel (cell array of char strings);
%	   DATA.OutputUnit     = unit of measured data (cell array of char strings). 
%
%	Example: 
%	   DATA.Name = 'T_Bucharest' ; 
%	   DATA.Notes = 'Systematic measurements on 2 channels' ; 
%	   DATA.ExperimentName = {'Daily min and max average temperatures in Bucharest'} ; 
%	   DATA.TimeUnit = 'day' ; 
%	   DATA.Tstart = datenum('29-Nov-2007 12:00:00') ; 
%	   DATA.Ts = 1 ; 
%	   DATA.OutputName = {'Minimum temperature' ; 
%			      'Maximum temperature'} ; 
%	   DATA.OutputUnit = {'ºC' ; 
%			      'ºC'} ; 
%	Note: The field DATA.Tstart cannot be set unless the sampling is uniform. 
%	      Non uniform sampling leaves the field empty. 
%	      Shall the starting date of measurements has to be specified, 
%	      the fiels DATA.UserData can be used. For example: 
%	      DATA.UserData = '29-Nov-2007' ; 
%
%	Empty input enforces empty output. 
%
%	Uses:	 VECTORIZE
%		 WAR_ERR
%
%	Author:  Dan STEFANOIU
%   2nd Revisor: Lavinius Ioan GLIGA
%	Created: March    08, 2009
%	Revised: November 01, 2010
%            February 26, 2018

function DATA = make_DATA(y)

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
	if (nargin < 1)
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
%	DATA.Domain = 'Time' ; 			% Data are in time domain. 
	Y = DATA ; 
	eval(['save ' DATA.Name '.mat Y']) ;	% Save the DATA object. 
	war_err(sprintf(S,DATA.Name)) ; 
%
%  END
%