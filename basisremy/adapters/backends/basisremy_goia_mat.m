function ok = basisremy_goia_mat(src, dst, var)
% BASISREMY_GOIA_MAT  Write a generated GOIA waveform as the RF struct MRSCloud loads.
%
%   ok = basisremy_goia_mat(src, dst, var)
%
%   MRSCloud's load_parameters.m runs `load Philips_GOIA_WURST_100pts.mat` (variable
%   Sweep2) or `load GE_GOIA_WURST_100pts.mat` (Sweep_GE_100) for sLASER; neither file
%   is in the public repository. BasisREMY generates a GOIA-WURST (core/pulse_library.py)
%   as a 4-column FID-A .txt (src); this function turns it into the FID-A RF struct and
%   saves it as `var` in dst.
%
%   The waveform is read by BasisREMY's headless io_loadRFwaveform (this folder), not
%   MRSCloud's bundled copy, which asks the user to type w1max for phase-modulated pulses.

here = fileparts(mfilename('fullpath'));
p = path();
addpath(here, '-begin');
rf = io_loadRFwaveform(src, 'ref', 0);
path(p);
eval([var ' = rf;']);
save('-v7', dst, var);
ok = 1;
end
