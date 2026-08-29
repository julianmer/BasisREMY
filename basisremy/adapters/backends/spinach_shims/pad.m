function s=pad(s,n,side,padchar)
% BasisREMY Octave shim for MATLAB's pad(str,n[,side[,padchar]]) on char rows
% (Spinach's report() pads every log line with it).
if nargin<2, n=numel(s); end
if nargin<3, side='right'; end
if nargin<4, padchar=' '; end
m=n-numel(s); if m<=0, return; end
switch lower(side)
    case 'left',  s=[repmat(padchar,1,m) s];
    case 'both',  s=[repmat(padchar,1,floor(m/2)) s repmat(padchar,1,ceil(m/2))];
    otherwise,    s=[s repmat(padchar,1,m)];
end
end
