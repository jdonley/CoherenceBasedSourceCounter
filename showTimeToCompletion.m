function num_curr_char = showTimeToCompletion( percent_complete, num_prev_char )
%SHOWTIMETOCOMPLETION Prints the time to completion and expected finish of a looped simulation based on linear extrapolation.
% 
% Syntax:	[ num_curr_char ] = showTimeToCompletion( percent_complete, num_prev_char )
%   Note that before using this function in a loop the in-built MATLAB
%   function tic should be called.
% 
% Inputs: 
% 	percent_complete - A decimal number between 0 and 1 representing the
% 	percentage completion.
% 	num_prev_char - Number of previous characters printed to the screen
% 	(Usually ok to begin with 0 and then reuse num_curr_char)
% 
% Outputs: 
% 	num_curr_char - Number of characters printed to the screen. Usually
% 	feed this number back into this function on the next iteration or
% 	increment appropriately if other characters have been printed between
% 	function calls.
%
% Example: 
%       fprintf('\t Completion: ');
%       n=0; tic;
%       len=1e2;
%       for i = 1:len
%           pause(1);
%           n = showTimeToCompletion( i/len, n);
%       end
% 
% See also: tic, toc

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2017
% Date: 25 August 2017
% Version: 1.0 (Simplified)
% 
% Original Source URL: https://github.com/JacobD10/SoundZone_Tools
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
tElapsed = toc;
ratio = percent_complete;
tRem = (1-ratio) / ratio * tElapsed;
tTot = tElapsed + tRem;

fprintf(repmat('\b',1,num_prev_char));
num_curr_char=fprintf( ...
    ['%.2f%%\n' ...
    '      Remaining: %s\n' ...
    '          Total: %s\n' ...
    'Expected Finish: %s\n'], ...
    ratio * 100, ...
    datestr(seconds(tRem),'hh:MM:SS'), ... floor(tRem/60), ...    rem(tRem,60), ...
    datestr(seconds(tTot),'hh:MM:SS'), ... floor(tTot/60), ...    rem(tTot,60), ...
    [strrep(datestr(datetime + seconds(tRem),'hh:MM:SS AM'),' ',''),'  ', ...
     datestr(datetime + seconds(tRem),'dd-mmm-yyyy')]);

end

