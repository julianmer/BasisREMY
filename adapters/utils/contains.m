% contains.m
% Adapter function for MATLAB's "contains" to ensure compatibility with Octave.
%
% AU
%
% USAGE:
% tf = contains(str, pattern)
%
% DESCRIPTION:
% Checks if the input string or cell array of strings contains the specified pattern.
% This function mimics MATLAB's "contains" using "strfind", which is available in Octave.
%
% INPUTS:
% str     = Input string or cell array of strings
% pattern = Pattern to search for (string)
%
% OUTPUTS:
% tf      = Logical array indicating if pattern is found in each element of str
%
% EXAMPLE:
% tf = contains({'Hello', 'World'}, 'or'); % returns [false true]
%
% NOTE:
% Place this file in your project directory or adapters folder to override missing "contains" in Octave.

function tf = contains(str, pattern)
    if ischar(str)
        tf = ~isempty(strfind(str, pattern));
    elseif iscell(str)
        tf = cellfun(@(s) ~isempty(strfind(s, pattern)), str);
    else
        error('Input must be a string or cell array of strings.');
    end
end