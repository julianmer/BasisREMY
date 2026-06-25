function [fid_re, fid_im, npts, sw_out, cf_mhz] = mrscloud_run_metab( ...
        metab, vendor, sequence, localization, te, field_str, ...
        edit_target, edit_on, edit_off, edit_tp, spatial_points, save_dir, ...
        samples, bandwidth)
% MRSCLOUD_RUN_METAB  Adapter that runs the MRSCloud workflow for ONE metabolite.
%
%   [fid_re, fid_im, npts, sw_out, cf_mhz] = ...
%       mrscloud_run_metab(metab, vendor, sequence, localization, te, ...
%                          field_str, edit_target, edit_on, edit_off, ...
%                          edit_tp, spatial_points, save_dir, ...
%                          samples, bandwidth)
%
%   This function mirrors the per-metabolite portion of MRSCloud's
%   externals/mrscloud/run/run_simulations_cloud.m but exposes a clean
%   numeric interface so the Python backend (BasisREMY) can drive it via
%   oct2py without dealing with structs, JSON files, or zip output.
%
%   Inputs (additions vs. upstream)
%     samples         double  number of complex spectral points (overrides
%                             MRSCloud's hard-coded Npts = 8192)
%     bandwidth       double  spectral width in Hz (overrides MRSCloud's
%                             hard-coded sw = 4000)
%   The other inputs are unchanged — see file header above.
    if nargin < 13 || isempty(samples);   samples   = 0; end
    if nargin < 14 || isempty(bandwidth); bandwidth = 0; end
% MRSCLOUD_RUN_METAB  Adapter that runs the MRSCloud workflow for ONE metabolite.
%
%   [fid_re, fid_im, npts, sw_out, cf_mhz] = ...
%       mrscloud_run_metab(metab, vendor, sequence, localization, te, ...
%                          field_str, edit_target, edit_on, edit_off, ...
%                          edit_tp, spatial_points, save_dir)
%
%   This function mirrors the per-metabolite portion of MRSCloud's
%   externals/mrscloud/run/run_simulations_cloud.m but exposes a clean
%   numeric interface so the Python backend (BasisREMY) can drive it via
%   oct2py without dealing with structs, JSON files, or zip output.
%
%   Inputs
%     metab           char  metabolite name, e.g. 'NAA'
%     vendor          char  'Philips' | 'Philips_universal' | 'Siemens' | 'GE'
%     sequence        char  'UnEdited' | 'MEGA' | 'HERMES' | 'HERCULES'
%     localization    char  'PRESS' | 'sLASER' | 'STEAM_7T'
%     te              double  echo time [ms]
%     field_str       char  '1.5T' | '3T' | '7T'
%     edit_target     char  'GABA' | 'GSH' | 'Lac' | 'PE'  ('' for unEdited)
%     edit_on         double  editing-on offset [ppm] (MEGA only)
%     edit_off        double  editing-off offset [ppm] (MEGA only)
%     edit_tp         double  editing-pulse duration [ms]
%     spatial_points  double  spatial grid (41 acceptable, 101 ideal)
%     save_dir        char  scratch directory MRSCloud needs for intermediate .mat
%
%   Outputs
%     fid_re, fid_im  Nx1 real/imag parts of the simulated FID
%     npts            number of complex points in the FID
%     sw_out          spectral width [Hz]
%     cf_mhz          centre / Larmor frequency [MHz]
%
%   For edited sequences (MEGA, HERMES, HERCULES) the returned FID is the
%   first sub-spectrum (off / 'a'); higher-order sub-spectra remain on disk
%   in `save_dir` and can be loaded explicitly if needed.
%
%   See also: externals/mrscloud/run/run_simulations_cloud.m
%             externals/mrscloud/functions/load_parameters.m

    if ~exist(save_dir, 'dir'); mkdir(save_dir); end

    % ---------- build MRS_temp like run_simulations_cloud.m does ----------
    MRS_temp                = struct();
    MRS_temp.metab          = metab;
    MRS_temp.Nmetab         = 1;
    MRS_temp.flipAngle      = 180;
    MRS_temp.centreFreq     = 4.65;  % carrier at water (4.65 ppm) — must match plot offset
    MRS_temp.edit_flipAngle = 180;
    MRS_temp.excflipAngle   = 90;
    MRS_temp.nX             = spatial_points;
    MRS_temp.nY             = spatial_points;
    MRS_temp.nZ             = spatial_points;
    MRS_temp.localization   = {localization};
    MRS_temp.FieldStr       = {field_str};
    MRS_temp.vendor         = {vendor};
    MRS_temp.seq            = {sequence};
    MRS_temp.TEs            = {te};
    MRS_temp.save_dir       = save_dir;

    % Editing parameters. We always set them — even for UnEdited — so that a
    % stray sequence value never produces the cryptic Octave error
    %   "structure has no member 'editON'"
    % from load_parameters.m line ~206. The values are ignored for UnEdited.
    MRS_temp.editTp = edit_tp;
    switch sequence
        case 'MEGA'
            MRS_temp.editON = num2cell([edit_on edit_off]);
        case {'HERMES', 'HERCULES', 'HERMES_GABA_GSH_EtOH'}
            % MRSCloud overrides TE=80 internally for HERMES/HERCULES.
            % These four offsets are the canonical HERMES/HERCULES pattern.
            MRS_temp.editON = num2cell([4.56 1.90 (4.56+1.9)/2 7.5]);
        otherwise
            % UnEdited (or anything unrecognised). Provide a benign stub so
            % load_parameters cannot crash if it ever inspects the field.
            MRS_temp.editON = num2cell([1.90 7.50 4.56 4.18]);
    end

    % Hard guard: only the four sequences below are supported by MRSCloud.
    if ~any(strcmp(sequence, {'UnEdited','MEGA','HERMES','HERCULES','HERMES_GABA_GSH_EtOH'}))
        error('mrscloud_run_metab:badSequence', ...
              'Sequence must be UnEdited / MEGA / HERMES / HERCULES (got "%s")', sequence);
    end

    % ---------- load defaults (pulse waveforms, geometry, etc.) ----------
    if strcmp(field_str, '1.5T')
        MRS_opt = load_parameters_1_5T(MRS_temp);
    else
        MRS_opt = load_parameters(MRS_temp);
    end

    % MRSCloud's load_parameters hard-codes Npts = 8192 and sw = 4000
    % regardless of the user's acquisition parameters. Override here so
    % the returned FID matches the (Samples, Bandwidth) the GUI was
    % configured with — otherwise the basis set is silently 4× too long
    % at the wrong spectral width.
    if samples > 0
        MRS_opt.Npts = double(samples);
    end
    if bandwidth > 0
        MRS_opt.sw = double(bandwidth);
    end

    % ---------- dispatch to the right simulator ----------
    switch localization
        case 'PRESS'
            switch sequence
                case {'HERMES','HERCULES','HERMES_GABA_GSH_EtOH'}
                    [~, outA, ~, ~, ~] = sim_signals(MRS_opt);
                case 'MEGA'
                    [~, outA, ~]       = sim_signals(MRS_opt);
                otherwise
                    [~, outA]          = sim_signals(MRS_opt);
            end
        case 'sLASER'
            switch sequence
                case {'HERMES','HERCULES','HERMES_GABA_GSH_EtOH'}
                    if strcmp(vendor, 'GE')
                        [~, outA, ~, ~, ~] = sim_signals_GE_sLASER_HERMES(MRS_opt);
                    else
                        [~, outA, ~, ~, ~] = sim_signals_sLASER(MRS_opt);
                    end
                case 'MEGA'
                    [~, outA, ~] = sim_signals_sLASER(MRS_opt);
                otherwise
                    [~, outA]    = sim_signals_sLASER(MRS_opt);
            end
        case 'STEAM_7T'
            [~, outA] = sim_signals_STEAM(MRS_opt);
        otherwise
            error('mrscloud_run_metab: unsupported localization "%s"', localization);
    end

    % ---------- pull the complex FID out of the FID-A-style struct ----------
    if isstruct(outA) && isfield(outA, 'fids')
        fid = outA.fids(:);
    elseif isnumeric(outA)
        fid = outA(:);
    else
        error('mrscloud_run_metab: unexpected output type for %s', metab);
    end

    fid_re = real(fid);
    fid_im = imag(fid);
    npts   = numel(fid);

    if isstruct(outA) && isfield(outA, 'spectralwidth')
        sw_out = outA.spectralwidth;
    else
        sw_out = NaN;
    end
    if isstruct(outA) && isfield(outA, 'txfrq')
        cf_mhz = outA.txfrq / 1e6;
    elseif isstruct(outA) && isfield(outA, 'Bo')
        cf_mhz = outA.Bo * 42.577;
    else
        cf_mhz = NaN;
    end
end






