function [fid_re, fid_im, npts, sw_out, cf_mhz, sub_re, sub_im] = mrscloud_run_metab( ...
        metab, vendor, sequence, localization, te, field_str, ...
        edit_target, edit_on, edit_off, edit_tp, spatial_points, save_dir, ...
        samples, bandwidth, b0, design_json)
% MRSCLOUD_RUN_METAB  Adapter that runs the MRSCloud workflow for ONE metabolite.
%
%   [fid_re, fid_im, npts, sw_out, cf_mhz, sub_re, sub_im] = ...
%       mrscloud_run_metab(metab, vendor, sequence, localization, te, ...
%                          field_str, edit_target, edit_on, edit_off, ...
%                          edit_tp, spatial_points, save_dir, ...
%                          samples, bandwidth, b0)
%
%   Outputs (additions vs. upstream)
%     sub_re, sub_im  the further sub-spectra as columns: MEGA — the edit-OFF
%                     scan (fid_re/fid_im is the edit-ON one); HERMES /
%                     HERCULES — sub-experiments B, C, D (fid is A). Empty
%                     for UnEdited.
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
%     b0              double  field strength in T (overrides the per-vendor
%                             value load_parameters hard-codes; 0 = keep it)
%   A BasisREMY sequence design (PRESS / sLASER and their MEGA / HERMES / HERCULES
%   variants) as one JSON string (basisremy/backends/mrscloud_backend.py _design_args);
%   omitted or '' = MRSCloud's own pulses and timing:
%     layout      'PRESS' | 'sLASER'
%     ref_file    refocusing waveform (.pta; .txt with the gradient column for GOIA / FOCI)
%     ref_tp      its duration [ms];  ref_ideal 1: a hard pulse at one spatial point, no
%                 gradient;  ref_gm 1: gradient-modulated (the waveform's own gradient)
%     thk         refocusing slab [cm]
%     te1         first echo time [ms] (unedited PRESS; 0 = not set)
%     taus        pulse-to-pulse gaps [ms] in MRSCloud's order (edited PRESS, sLASER);
%                 empty = MRSCloud's own (needs the BasisREMY patch of sim_signals*.m)
%     edit_file   single-lobe editing waveform; edit_tp its duration [ms]
%     edit_files  one waveform per sub-experiment (dual-lobe ones generated), edit_dual
%                 1 where dual-lobe, edit_on the editing centres [ppm] (MRSCloud's editON)
%   The other inputs are unchanged — see file header above.
    if nargin < 13 || isempty(samples);   samples   = 0; end
    if nargin < 14 || isempty(bandwidth); bandwidth = 0; end
    if nargin < 15 || isempty(b0);        b0        = 0; end
    design = struct();
    if nargin >= 16 && ~isempty(design_json)
        if exist('jsondecode', 'builtin') || exist('jsondecode', 'file')
            design = jsondecode(design_json);
        else
            design = loadjson(design_json);          % jsonlab, bundled with MRSCloud
        end
    end
    sub_re = []; sub_im = [];
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
%     field_str       char  '1.5T' | '3T' — which MRSCloud parameter set
%                           (pulses, timings) to load; the field itself is b0
%     edit_target     char  accepted for interface stability, unused
%                           (MRSCloud fixes the editing offsets per scheme)
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
%   MEGA returns both sub-spectra (fid = edit-ON at edit_on, sub = edit-OFF
%   at edit_off). HERMES / HERCULES return sub-experiment A in fid and B, C,
%   D as the columns of sub (editing pulses at the scheme's four offsets, in
%   the order of MRS_temp.editON).
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
    if isfield(design, 'edit_on') && ~isempty(design.edit_on)
        sequence_editON = num2cell(double(design.edit_on(:))');   % the design's editing centres
    else
        sequence_editON = {};
    end
    switch sequence
        case 'MEGA'
            MRS_temp.editON = num2cell([edit_on edit_off]);
        % Fixed editing schemes, mirrored from run_simulations_cloud.m
        % (A / B single-lobe, C / D dual-lobe pulses).
        case 'HERMES'
            MRS_temp.editON = num2cell([4.56 1.90 (4.56+1.9)/2 7.50]);
        case 'HERMES_GABA_GSH_EtOH'
            MRS_temp.editON = num2cell([(4.56+1.9)/2 (3.67+1.9)/2 (3.67+4.56)/2 7.50]);
        case 'HERCULES'
            MRS_temp.editON = num2cell([4.58 4.18 (4.58+1.9)/2 (4.18+1.9)/2]);
        otherwise
            % UnEdited (or anything unrecognised). Provide a benign stub so
            % load_parameters cannot crash if it ever inspects the field.
            MRS_temp.editON = num2cell([1.90 7.50 4.56 4.18]);
    end

    if ~isempty(sequence_editON)
        MRS_temp.editON = sequence_editON;
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

    % ---------- simulate at the acquisition's own field strength ----------
    % load_parameters hard-codes B0 per vendor (Siemens 2.89 T, otherwise
    % 3 T; 7 T only for STEAM_7T). Rebuild what it computed at that field:
    % the Hamiltonian, and the editing pulses' frequency shifts (the offset
    % in Hz scales with B0; rf_freqshift composes additively).
    if b0 > 0 && abs(b0 - MRS_opt.Bfield) > 1e-6
        b0_old = MRS_opt.Bfield;
        MRS_opt.Bfield = double(b0);
        [MRS_opt.H, MRS_opt.d] = sim_Hamiltonian_mgs(MRS_opt.sys, MRS_opt.Bfield);
        % the refocusing propagators were built with the old Hamiltonian
        % (load_parameters.m, "Creating propagators for Refoc pulse") — redo
        if isfield(MRS_opt, 'Qrefoc')
            for X = 1:length(MRS_opt.x)
                MRS_opt.Qrefoc{X} = calc_shapedRF_propagator_refoc( ...
                    MRS_opt.H, MRS_opt.refRF, MRS_opt.refTp, MRS_opt.flipAngle, ...
                    0, MRS_opt.y(X), MRS_opt.Gx);
            end
        end
        names = {'editRFonA', 'editRFonB', 'editRFonC', 'editRFonD'};
        for k = 1:numel(names)
            if isfield(MRS_opt, names{k}) && k <= numel(MRS_opt.editON)
                dF = (MRS_opt.centreFreq - MRS_opt.editON{k}) ...
                     * (MRS_opt.Bfield - b0_old) * MRS_opt.gamma / 1e6;
                MRS_opt.(names{k}) = rf_freqshift(MRS_opt.(names{k}), MRS_opt.editTp, dF);
            end
        end
    end

    % ---------- a BasisREMY sequence design: its pulses and TE1 ----------
    if ~isempty(fieldnames(design))
        MRS_opt = apply_design(MRS_opt, design, te);
    end

    % ---------- dispatch to the right simulator ----------
    switch localization
        case 'PRESS'
            switch sequence
                case {'HERMES','HERCULES','HERMES_GABA_GSH_EtOH'}
                    [~, outA, outB, outC, outD] = sim_signals(MRS_opt);
                case 'MEGA'
                    [~, outA, outB]    = sim_signals(MRS_opt);
                otherwise
                    [~, outA]          = sim_signals(MRS_opt);
            end
        case 'sLASER'
            switch sequence
                case {'HERMES','HERCULES','HERMES_GABA_GSH_EtOH'}
                    if strcmp(vendor, 'GE')
                        [~, outA, outB, outC, outD] = sim_signals_GE_sLASER_HERMES(MRS_opt);
                    else
                        [~, outA, outB, outC, outD] = sim_signals_sLASER(MRS_opt);
                    end
                case 'MEGA'
                    [~, outA, outB] = sim_signals_sLASER(MRS_opt);
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
    switch sequence
        case 'MEGA'                                          % edit-OFF
            sub = outB.fids(:);
        case {'HERMES', 'HERCULES', 'HERMES_GABA_GSH_EtOH'}  % B, C, D
            sub = [outB.fids(:), outC.fids(:), outD.fids(:)];
        otherwise
            sub = [];
    end
    sub_re = real(sub);
    sub_im = imag(sub);

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


function MRS_opt = apply_design(MRS_opt, design, te)
% Replace MRSCloud's refocusing / editing pulses and timing by the design's and rebuild what
% load_parameters derived from them (slice gradient, grid, refocusing propagators).
    if isfield(design, 'te1') && design.te1 > 0
        MRS_opt.TE1 = design.te1;
        MRS_opt.TE2 = te - design.te1;
    end
    if isfield(design, 'taus') && ~isempty(design.taus)
        MRS_opt.taus_design = double(design.taus(:))';   % read by the patched sim_signals*.m
    end
    if isfield(design, 'thk') && design.thk > 0
        MRS_opt.thkX = design.thk;
        MRS_opt.thkY = design.thk;
    end
    if isfield(design, 'edit_file') && ~isempty(design.edit_file)
        single = load_rf_headless(design.edit_file, 'inv');
        files = {design.edit_file, design.edit_file};
        dual = [0 0];
        if isfield(design, 'edit_files') && ~isempty(design.edit_files)
            files = cellstr(design.edit_files);
            dual = double(design.edit_dual(:))';
        end
        for k = 1:numel(files)
            if dual(k)
                rf = load_rf_headless(files{k}, 'inv');
                rf.tw1 = single.tw1 * 2;              % two lobes: MRSCloud's dual-lobe scaling
            else
                rf = single;
            end
            MRS_opt.(sprintf('editRF%d', k)) = rf;
        end
    end
    if ~isfield(design, 'ref_file') || isempty(design.ref_file)
        return
    end
    MRS_opt.refRF = load_rf_headless(design.ref_file, 'ref');
    MRS_opt.refTp = design.ref_tp;
    gm = isfield(design, 'ref_gm') && design.ref_gm;
    if isfield(design, 'ref_ideal') && design.ref_ideal
        % a hard pulse acts on every position alike: one point, no gradient
        MRS_opt.x = 0; MRS_opt.y = 0; MRS_opt.nX = 1; MRS_opt.nY = 1;
        MRS_opt.Gx = 0;
    else
        if gm
            % gradient-modulated: the waveform's own gradient column (G/cm at the design slab)
            MRS_opt.refRF.f0 = 0;
            MRS_opt.refRF.tthk = MRS_opt.thkX * (MRS_opt.refTp / 1000);
            MRS_opt.Gx = 0;
        else
            MRS_opt.Gx = (MRS_opt.refRF.tbw / (MRS_opt.refTp / 1000)) ...
                         / (MRS_opt.gamma * MRS_opt.thkX / 10000);        % [G/cm]
        end
        MRS_opt.fovX = 1.5 * MRS_opt.thkX;                                % slab + transition bands
        MRS_opt.fovY = MRS_opt.fovX;
        MRS_opt.x = linspace(-MRS_opt.fovX / 2, MRS_opt.fovX / 2, MRS_opt.nX);
        MRS_opt.y = linspace(-MRS_opt.fovY / 2, MRS_opt.fovY / 2, MRS_opt.nY);
    end
    Qrefoc = cell(1, numel(MRS_opt.y));
    for X = 1:numel(MRS_opt.y)
        Qrefoc{X} = calc_shapedRF_propagator_refoc(MRS_opt.H, MRS_opt.refRF, MRS_opt.refTp, ...
                                                   MRS_opt.flipAngle, 0, MRS_opt.y(X), MRS_opt.Gx);
    end
    MRS_opt.Qrefoc = Qrefoc;
end


function rf = load_rf_headless(file, type)
% BasisREMY's io_loadRFwaveform (this folder): reads .pta / .RF / .txt (with a gradient
% column) without the interactive w1max prompt MRSCloud's copy makes for phase-modulated
% pulses.
    here = fileparts(mfilename('fullpath'));
    p = path();
    addpath(here, '-begin');
    try
        rf = io_loadRFwaveform(file, type, 0);
    catch err
        path(p);
        rethrow(err);
    end
    path(p);
end
