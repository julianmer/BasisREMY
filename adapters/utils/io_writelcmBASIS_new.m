% io_writelcmBASIS_new
% -------------------------------------------------------------------------
%
% Writes an LCModel .BASIS file from a cell array of FID-A compatible metabolite outputs.
%
% USAGE:
%   BASIS = io_writelcmBASIS_new(outputs, outfile, vendor, SEQ, B0, linewidth, Npts, sw, te, centreFreq)
%
% INPUTS:
%   outputs     = Cell array of metabolite structs (each containing fids, t, name, centerFreq, w1max)
%   outfile     = Full path and filename without extension for the .BASIS file
%   vendor      = System label (e.g., 'Siemens', 'Philips', 'GE')
%   SEQ         = Sequence name (e.g., 'sLASER', 'PRESS', 'MEGA', 'HERMES')
%   B0          = Magnetic field strength in Tesla
%   linewidth   = Linewidth in Hz
%   Npts        = Number of spectral points
%   sw          = Spectral width in Hz
%   te          = Echo time in ms
%   centreFreq  = Center frequency in ppm
%
% DESCRIPTION:
%   This function converts a set of simulated metabolite FIDs into an LCModel-
%   compatible .BASIS file. Each metabolite is written with its real and imaginary
%   components, using the acquisition parameters provided. No zero-filling is applied.
%
% AUTHOR:
%   J. P. Merkofer
%
% CREDITS:
%   Based on FID-A and Osprey toolbox functions and conventions
% ------------------------------------------------------------------------

function BASIS = io_writelcmBASIS_new(outputs, outfile, vendor, SEQ, B0, linewidth, Npts, sw, te, centreFreq)

% Ensure output directory exists
[outDir, baseName, ext] = fileparts(outfile);
disp(['Parsed outfile path - Directory: "' outDir '", Base: "' baseName '", Ext: "' ext '"']);

% Handle empty directory (use current directory)
if isempty(outDir)
    outDir = pwd;
    disp(['WARNING: No directory in outfile path. Using current directory: ' outDir]);
    outfile = fullfile(outDir, [baseName ext]);
end

% Check and create directory
if ~exist(outDir, 'dir')
    disp(['Creating directory: ' outDir]);
    [status, msg, msgid] = mkdir(outDir);
    if status == 1
        disp(['Successfully created directory: ' outDir]);
    else
        error('Failed to create directory: %s. Error: %s', outDir, msg);
    end
else
    disp(['Directory already exists: ' outDir]);
end

% Ensure .basis extension
if isempty(ext)
    outfile = [outfile '.basis'];
    disp(['Added .basis extension. Final path: ' outfile]);
end

% Open file for writing
fid = fopen(outfile, 'w+');
if fid < 0
    [errmsg, errnum] = ferror(fid);
    error('FOPEN FAILED! Could not open file for writing: %s\nError number: %d\nError message: %s\nCurrent dir: %s\nAttempted path: %s', ...
          outfile, errnum, errmsg, pwd, outfile);
end

% Write header (like io_writelcmBASIS)
fprintf(fid,' $SEQPAR\n');
fprintf(fid,' FWHMBA = %5.6f,\n', linewidth/(42.577*B0));
fprintf(fid,' HZPPPM = %5.6f,\n', 42.577*B0);
fprintf(fid,' ECHOT = %5.2f,\n', te);
fprintf(fid,' SEQ = ''%s''\n', SEQ);
fprintf(fid,' $END\n');

fprintf(fid,' $BASIS1\n');
fprintf(fid,' IDBASI = ''%s %s %d ms'',\n', vendor, SEQ, round(te));
fprintf(fid,' FMTBAS = ''(2E15.6)'',\n');
fprintf(fid,' BADELT = %5.6f,\n', 1/sw);
fprintf(fid,' NDATAB = %i\n', Npts);
fprintf(fid,' $END\n');

% Loop through outputs and write each metabolite
for m = 1:numel(outputs)
    out = outputs{m};
    specs = ifft(out.fids(:));
    RF = [real(specs), imag(specs)];

    fprintf(fid,' $NMUSED\n XTRASH = 0\n $END\n');
    fprintf(fid,' $BASIS\n ID = ''%s'',\n METABO = ''%s'',\n CONC = 1.,\n TRAMP = 1.,\n VOLUME = 1.,\n ISHIFT = 0\n AUTOPH = F\n AUTOSC = F\n NOSHIF = T\n $END\n', out.name, out.name);
    fprintf(fid,' %7.6e  %7.6e\n', RF');
end

fclose(fid);

% Return the output path
BASIS = outfile;
