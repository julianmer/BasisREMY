function [fid_re, fid_im, npts, sw_out, cf_mhz] = fida_run(metab, kind, varargin)
% FIDA_RUN  Single-entry dispatcher for the BasisREMY FID-A backend family.
%
%   [fid_re, fid_im, npts, sw_out, cf_mhz] = fida_run(metab, kind, ...)
%
%   This adapter is the ONLY MATLAB entry point for the FID-A backend family
%   in BasisREMY. The Python side picks a `kind` (e.g. 'ideal' / 'press_shaped')
%   and the trailing positional arguments; this function loads the requested
%   spin system from `metabolites/spinSystems.mat`, dispatches to the right
%   FID-A simulator, and returns the FID as plain numeric vectors so oct2py
%   can ferry them to Python without dealing with structs.
%
%   To add a new FID-A simulator: append a `case` to the switch below and
%   add a Python subclass of FidaBackend that forwards the right
%   parameters via its `_build_args()` method.
%
%   Inputs (always)
%     metab  : char  metabolite name (must match a `sys<NAME>` field in
%              externals/fidA/simulationTools/metabolites/spinSystems.mat)
%     kind   : char  dispatch key — see switch below
%
%   Outputs
%     fid_re, fid_im : Nx1 real / imag parts of the FID
%     npts           : length of the FID
%     sw_out         : spectral width [Hz]
%     cf_mhz         : carrier frequency [MHz]

    % ---------- spin system --------------------------------------------
    S = load('metabolites/spinSystems.mat');
    sysFieldName = ['sys' metab];
    if ~isfield(S, sysFieldName)
        error('fida_run: unknown metabolite "%s" (no %s in spinSystems.mat)', ...
              metab, sysFieldName);
    end
    sys = S.(sysFieldName);

    % ---------- dispatch -----------------------------------------------
    sw_out = NaN;   % most simulators set this; defaults for safety
    cf_mhz = NaN;
    switch lower(kind)

        % ----- IDEAL  (Spin Echo / PRESS / STEAM / LASER) ----------
        % Use 2 outputs: [RF, out] = sim_lcmrawbasis(...).  We ignore RF
        % and take out.fids directly so we get the raw FID-A signal WITHOUT
        % the op_complexConj that io_writelcmraw would apply (which would
        % invert the spectrum if we later process with fft).
        % makeraw='n' suppresses the .RAW file write inside sim_lcmrawbasis;
        % a separate export path handles that if needed.
        %   args: n, sw, Bfield, lw, tau1, tau2, seq, out_path
        case 'ideal'
            [n, sw, Bfield, lw, tau1, tau2, seq, out_path] = ...
                deal(varargin{1:8});
            [~, out] = sim_lcmrawbasis(n, sw, Bfield, lw, metab, ...
                                       tau1, tau2, 'n', 'n', seq, out_path);
            fid    = out.fids(:);
            fid_re = real(fid); fid_im = imag(fid);
            npts   = numel(fid);
            sw_out = sw;
            cf_mhz = Bfield * 42.577;

        % ----- PRESS shaped --------------------------------------------
        %   args: n, sw, Bfield, lw, tau1, tau2, pulse_path, tp,
        %         thkX, thkY, fovX, fovY, nX, nY, flipAngle, centreFreq
        case 'press_shaped'
            [n, sw, Bfield, lw, tau1, tau2, pulse_path, tp, ...
             thkX, thkY, fovX, fovY, nX, nY, flipAngle, centreFreq] = ...
                deal(varargin{1:16});
            if ~exist(pulse_path, 'file')
                error('fida_run/press_shaped: pulse waveform not found: "%s"', pulse_path);
            end
            % NOTE: io_loadRFwaveform expects type ∈ {'exc','ref','inv'} (or
            % a numeric flip angle). Passing 'refoc' triggers a length-mismatch
            % crash inside its `type=='exc'` test ("mx_el_eq: nonconformant
            % arguments (op1 is 1x5, op2 is 1x3)").
            RF = io_loadRFwaveform(pulse_path, 'ref', 0);
            % For gradient-modulated (GM / adiabatic) pulses such as GOIA the
            % gradient waveform is already stored in column 4 of RF.waveform.
            % sim_shapedRF will ERROR if you ALSO supply an explicit Gx/Gy
            % ("You cannot supply GM pulse AND separately specify the Gradient
            % strength"). For non-GM pulses we derive Gx/Gy analytically from
            % the time-bandwidth product so the slice thickness matches thkX/Y.
            if RF.isGM
                Gx = 0;
                Gy = 0;
            else
                Gx = (RF.tbw / (tp/1000)) / (4258 * thkX);
                Gy = (RF.tbw / (tp/1000)) / (4258 * thkY);
            end
            if nX < 2; nX = 2; end
            if nY < 2; nY = 2; end
            x = linspace(-fovX/2, fovX/2, nX);
            y = linspace(-fovY/2, fovY/2, nY);
            accumFid = [];
            for ix = 1:nX
                for iy = 1:nY
                    out = sim_press_shaped(n, sw, Bfield, lw, sys, tau1, tau2, ...
                                           RF, tp, x(ix), y(iy), Gx, Gy, ...
                                           flipAngle, centreFreq);
                    if isempty(accumFid)
                        accumFid = out.fids(:);
                    else
                        accumFid = accumFid + out.fids(:);
                    end
                end
            end
            accumFid = accumFid / (nX * nY);
            fid_re = real(accumFid); fid_im = imag(accumFid);
            npts   = numel(accumFid);
            sw_out = sw;
            cf_mhz = Bfield * 42.577;

        % ----- LASER (ideal AFP, six equally spaced echoes) ------------
        %   args: n, sw, Bfield, lw, te
        case 'laser'
            [n, sw, Bfield, lw, te] = deal(varargin{1:5});
            out = sim_laser(n, sw, Bfield, lw, sys, te);
            fid    = out.fids(:);
            fid_re = real(fid); fid_im = imag(fid);
            npts   = numel(fid);
            sw_out = sw;
            cf_mhz = Bfield * 42.577;

        % ----- Spin Echo xN (multi-echo train) -------------------------
        %   args: n, sw, Bfield, lw, tau, nechoes
        case 'spinecho_xn'
            [n, sw, Bfield, lw, tau, nechoes] = deal(varargin{1:6});
            out = sim_spinecho_xN(n, sw, Bfield, lw, sys, tau, nechoes);
            fid    = out.fids(:);
            fid_re = real(fid); fid_im = imag(fid);
            npts   = numel(fid);
            sw_out = sw;
            cf_mhz = Bfield * 42.577;

        % ----- One pulse (ideal pulse-acquire FID) ---------------------
        %   args: n, sw, Bfield, lw
        case 'onepulse'
            [n, sw, Bfield, lw] = deal(varargin{1:4});
            out = sim_onepulse(n, sw, Bfield, lw, sys);
            fid    = out.fids(:);
            fid_re = real(fid); fid_im = imag(fid);
            npts   = numel(fid);
            sw_out = sw;
            cf_mhz = Bfield * 42.577;

        % ----- semi-LASER shaped (Oz et al. 2018; one AFP waveform) ----
        %   args: n, sw, Bfield, lw, te, pulse_path, tp,
        %         thkX, thkY, fovX, fovY, nX, nY, flipAngle, centreFreq
        case 'semilaser_shaped'
            [n, sw, Bfield, lw, te, pulse_path, tp, ...
             thkX, thkY, fovX, fovY, nX, nY, flipAngle, centreFreq] = ...
                deal(varargin{1:15});
            if ~exist(pulse_path, 'file')
                error('fida_run/semilaser_shaped: pulse waveform not found: "%s"', pulse_path);
            end
            RF = io_loadRFwaveform(pulse_path, 'ref', 0);
            gamma = 42577000;
            % Gradient conventions from jbss run_mysLASERShaped_fast.m:
            % GM (e.g. GOIA) pulses scale via the time*thickness product,
            % conventional pulses via the time-bandwidth product.
            if RF.isGM
                Gx = (RF.tthk / (tp/1000)) / thkX;
                Gy = (RF.tthk / (tp/1000)) / thkY;
            else
                Gx = (RF.tbw / (tp/1000)) / (gamma * thkX / 10000);
                Gy = (RF.tbw / (tp/1000)) / (gamma * thkY / 10000);
            end
            if nX < 2; nX = 2; end
            if nY < 2; nY = 2; end
            x = linspace(-fovX/2, fovX/2, nX);
            y = linspace(-fovY/2, fovY/2, nY);
            accumFid = [];
            for ix = 1:nX
                for iy = 1:nY
                    out = sim_semiLASER_shaped(n, sw, Bfield, lw, sys, te, ...
                                               RF, tp, x(ix), y(iy), Gx, Gy, ...
                                               flipAngle, centreFreq);
                    if isempty(accumFid)
                        accumFid = out.fids(:);
                    else
                        accumFid = accumFid + out.fids(:);
                    end
                end
            end
            accumFid = accumFid / (nX * nY);
            fid_re = real(accumFid); fid_im = imag(accumFid);
            npts   = numel(accumFid);
            sw_out = sw;
            cf_mhz = Bfield * 42.577;

        % ----- STEAM shaped (shaped 90s, per run_simSteamShaped.m) -----
        %   args: n, sw, Bfield, lw, te, tm, pulse_path, tp,
        %         thkX, thkY, fovX, fovY, nX, nY, flipAngle, centreFreq
        case 'steam_shaped'
            [n, sw, Bfield, lw, te, tm, pulse_path, tp, ...
             thkX, thkY, fovX, fovY, nX, nY, flipAngle, centreFreq] = ...
                deal(varargin{1:16});
            if ~exist(pulse_path, 'file')
                error('fida_run/steam_shaped: pulse waveform not found: "%s"', pulse_path);
            end
            RF = io_loadRFwaveform(pulse_path, 'exc', 0);
            gamma = 42577000;
            if RF.isGM
                Gx = (RF.tthk / (tp/1000)) / thkX;
                Gy = (RF.tthk / (tp/1000)) / thkY;
            else
                Gx = (RF.tbw / (tp/1000)) / (gamma * thkX / 10000);
                Gy = (RF.tbw / (tp/1000)) / (gamma * thkY / 10000);
            end
            if nX < 2; nX = 2; end
            if nY < 2; nY = 2; end
            x = linspace(-fovX/2, fovX/2, nX);
            y = linspace(-fovY/2, fovY/2, nY);
            accumFid = [];
            for ix = 1:nX
                for iy = 1:nY
                    out = sim_steam_shaped(n, sw, Bfield, lw, sys, te, tm, ...
                                           RF, tp, x(ix), y(iy), Gx, Gy, ...
                                           flipAngle, centreFreq);
                    if isempty(accumFid)
                        accumFid = out.fids(:);
                    else
                        accumFid = accumFid + out.fids(:);
                    end
                end
            end
            accumFid = accumFid / (nX * nY);
            fid_re = real(accumFid); fid_im = imag(accumFid);
            npts   = numel(accumFid);
            sw_out = sw;
            cf_mhz = Bfield * 42.577;

        % ----- Spin Echo shaped (1-D; subtractive [0,90] phase cycle) --
        %   args: n, sw, Bfield, lw, te, pulse_path, tp, thk, fov, npos
        case 'spinecho_shaped'
            [n, sw, Bfield, lw, te, pulse_path, tp, thk, fov, npos] = ...
                deal(varargin{1:10});
            if ~exist(pulse_path, 'file')
                error('fida_run/spinecho_shaped: pulse waveform not found: "%s"', pulse_path);
            end
            RF = io_loadRFwaveform(pulse_path, 'ref', 0);
            gamma = 42577000;
            if RF.isGM
                G = (RF.tthk / (tp/1000)) / thk;
            else
                G = (RF.tbw / (tp/1000)) / (gamma * thk / 10000);
            end
            if npos < 2; npos = 2; end
            pos = linspace(-fov/2, fov/2, npos);
            phCyc = [0, 90];   % per run_simSpinEchoShaped.m: RP1 - RP2
            accumFid = [];
            for ip = 1:npos
                out1 = sim_spinecho_shaped(n, sw, Bfield, lw, sys, te, ...
                                           RF, tp, G, pos(ip), phCyc(1));
                out2 = sim_spinecho_shaped(n, sw, Bfield, lw, sys, te, ...
                                           RF, tp, G, pos(ip), phCyc(2));
                posFid = out1.fids(:) - out2.fids(:);
                if isempty(accumFid)
                    accumFid = posFid;
                else
                    accumFid = accumFid + posFid;
                end
            end
            accumFid = accumFid / npos;
            fid_re = real(accumFid); fid_im = imag(accumFid);
            npts   = numel(accumFid);
            sw_out = sw;
            cf_mhz = Bfield * 42.577;

        % ----- MEGA-PRESS ideal (per-spin instantaneous editing) -------
        %   args: n, sw, Bfield, lw, t1..t5 [ms], editPpm, editBand, editOn
        %   The ideal editing pulse inverts every spin within
        %   editPpm +/- editBand/2 (in the metabolite's own ppm frame);
        %   editOn = 0 leaves all spins untouched (the OFF sub-spectrum).
        case 'megapress_ideal'
            [n, sw, Bfield, lw, t1, t2, t3, t4, t5, ...
             editPpm, editBand, editOn] = deal(varargin{1:12});
            taus = [t1 t2 t3 t4 t5];
            nsub = numel(sys);
            refoc1Flip = cell(1, nsub);
            refoc2Flip = cell(1, nsub);
            editFlip   = cell(1, nsub);
            for k = 1:nsub
                shifts = sys(k).shifts(:)';
                refoc1Flip{k} = 180 * ones(size(shifts));
                refoc2Flip{k} = 180 * ones(size(shifts));
                if editOn
                    editFlip{k} = 180 * double(abs(shifts - editPpm) <= editBand/2);
                else
                    editFlip{k} = zeros(size(shifts));
                end
            end
            out = sim_megapress(n, sw, Bfield, lw, sys, taus, ...
                                refoc1Flip, refoc2Flip, editFlip);
            fid    = out.fids(:);
            fid_re = real(fid); fid_im = imag(fid);
            npts   = numel(fid);
            sw_out = sw;
            cf_mhz = Bfield * 42.577;

        % ----- MEGA-PRESS shaped-edit (ideal refoc; shaped editing) ----
        %   Conventions from run_simMegaPressShapedEdit.m: edit pulse loaded
        %   as 'inv', frequency-shifted to the ON/OFF target, [0,90]x[0,90,
        %   180,270] edit phase cycle summed and averaged.
        %   args: n, sw, Bfield, lw, t1..t5 [ms], edit_path, editTp,
        %         editOnFreq, editOffFreq, centreFreq, edit_on_flag
        case 'megapress_shapededit'
            [n, sw, Bfield, lw, t1, t2, t3, t4, t5, edit_path, editTp, ...
             editOnFreq, editOffFreq, centreFreq, editOn] = ...
                deal(varargin{1:15});
            if ~exist(edit_path, 'file')
                error('fida_run/megapress_shapededit: edit pulse waveform not found: "%s"', edit_path);
            end
            taus = [t1 t2 t3 t4 t5];
            editRF = io_loadRFwaveform(edit_path, 'inv', 0);
            gamma = 42577000;
            if editOn
                targetFreq = editOnFreq;
            else
                targetFreq = editOffFreq;
            end
            editRFshift = rf_freqshift(editRF, editTp, ...
                (centreFreq - targetFreq) * Bfield * gamma / 1e6);
            editPhCyc1 = [0 90];
            editPhCyc2 = [0 90 180 270];
            accumFid = [];
            for EP1 = 1:numel(editPhCyc1)
                for EP2 = 1:numel(editPhCyc2)
                    out = sim_megapress_shapedEdit(n, sw, Bfield, lw, taus, ...
                        sys, editRFshift, editTp, ...
                        editPhCyc1(EP1), editPhCyc2(EP2), centreFreq);
                    if isempty(accumFid)
                        accumFid = out.fids(:);
                    else
                        accumFid = accumFid + out.fids(:);
                    end
                end
            end
            accumFid = accumFid / (numel(editPhCyc1) * numel(editPhCyc2));
            fid_re = real(accumFid); fid_im = imag(accumFid);
            npts   = numel(accumFid);
            sw_out = sw;
            cf_mhz = Bfield * 42.577;

        % ----- stubs ---------------------------------------------------
        % The Python side already raises NotImplementedError for these
        % kinds before reaching Octave, but we guard here as well so a
        % buggy caller gets a clear MATLAB-side error too.
        case {'megapress_shaped','megaspecial_shaped'}
            error('fida_run: kind "%s" is a registered FID-A wrapper but the Octave-side branch is not implemented yet.', kind);

        otherwise
            error('fida_run: unknown kind "%s"', kind);
    end
end

