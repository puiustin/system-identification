function [ existFlag ] = fig_exist( LookFor)
%
% fig_exist     Checks to see if a figure, with a given name, exists 
%
% Inputs:	LookFor  # the name of the figure whose existance should be
%                      checked
%
% Outputs:	existFlag # 0 if the figure with the specifier name doesn't 
%                       exist
%                     # 1 if the figure with the specified name exists
%
% Author:   Lavinius Ioan Gliga (**)
%
% Last upgrade: (*) February 26, 2018
%
% Copyright: (*)  "Politehnica" Unversity of Bucharest, ROMANIA
%                 Department of Automatic Control & Computer Science
%

% BEGIN
% 
% check to see if the user input an argument
if (nargin ~= 1)
    disp('No figure name has been given')
end

    % figflag
    existFlag = 0; %% the figure might not exist, we are pesimistic :)
    h = findobj(); %% get all object handles
    for i = 1 : length(h)
        %% check to see if each of them is a figure and if it the figure
        %% of interest
        if (strcmp(h(i).Type, 'figure') == 1 && ...
                strcmp( h(i).Name, LookFor) == 1) 
            existFlag = 1;
            break;
        end
    end
    % /figflag

end

