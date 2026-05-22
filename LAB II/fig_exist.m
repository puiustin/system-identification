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
% verificarea daca utilizatorul a introdus un argument
if (nargin ~= 1)
    disp('No figure name has been given')
end

    % initializare flag existenta (presupunem initial ca nu exista)
    existFlag = 0; 
    h = findobj(); %% obtinerea tuturor obiectelor grafice
    for i = 1 : length(h)
        %% verificam daca obiectul este de tip figura si daca are numele cautat
        if (strcmp(h(i).Type, 'figure') == 1 && ...
                strcmp( h(i).Name, LookFor) == 1) 
            existFlag = 1; % figura a fost gasita
            break;
        end
    end
    % /figflag

end

