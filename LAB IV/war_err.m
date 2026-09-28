%
%	File WAR_ERR.M
%
%	Function: WAR_ERR
%
%	Call: war_err(msg)
%
%	Displays a warning/error message. If msg is 
%	present and non void, then it is displayed. 
%	Otherwise (missing or void), a standard error 
%	message is displayed. 
%
%	Author:  Dan Stefanoiu
%   2nd Revisor: Lavinius Ioan Gliga
%	Created: June       2, 1999
%	Revised: September 22, 2007
%            February  26, 2018
%

function war_err(msg)

%
%	BEGIN
%
    if ((nargin > 0) && (~isempty(msg)))
	   disp(' ') ;
	   disp(msg) ;
	   disp(' ') ;  
    else
	   disp(' ') ;
	   disp('<WAR_ERR>: An error occured.') ;
	   disp(' ') ;
    end
%
%
%	END
%

