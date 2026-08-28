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
    sim_centre = 4.65;   % rotating-frame centre [ppm] of the FID the simulator returns
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
            % Gx/Gy [G/cm] from the time-bandwidth product so the slice
            % thickness matches thkX/Y; for gradient-modulated (GOIA-type)
            % waveforms the same numbers are the scale factors that
            % gm_prepare (end of this file) applies to the stored gradient.
            if RF.isGM
                Gx = (RF.tthk / (tp/1000)) / thkX;
                Gy = (RF.tthk / (tp/1000)) / thkY;
            else
                Gx = (RF.tbw / (tp/1000)) / (4258 * thkX);
                Gy = (RF.tbw / (tp/1000)) / (4258 * thkY);
            end
            [RF, Gx, Gy] = gm_prepare(RF, 'press_shaped', Gx, Gy);
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
            [RF, Gx, Gy] = gm_prepare(RF, 'steam_shaped', Gx, Gy);
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
            [RF, G] = gm_prepare(RF, 'spinecho_shaped', G);
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

        % ----- One pulse with ADC-onset delay (sim_onepulse_delay) --------
        %   args: n, sw, Bfield, lw, delay [ms]
        case 'onepulse_delay'
            [n, sw, Bfield, lw, delay] = deal(varargin{1:5});
            out = sim_onepulse_delay(n, sw, Bfield, lw, sys, delay);
            fid    = out.fids(:);
            fid_re = real(fid); fid_im = imag(fid);
            npts   = numel(fid);
            sw_out = sw;
            cf_mhz = Bfield * 42.577;

        % ----- One pulse with arbitrary excitation phase (sim_onepulse_arbPh)
        %   args: n, sw, Bfield, lw, phase [deg]
        case 'onepulse_arbph'
            [n, sw, Bfield, lw, ph] = deal(varargin{1:5});
            out = sim_onepulse_arbPh(n, sw, Bfield, lw, sys, ph);
            fid    = out.fids(:);
            fid_re = real(fid); fid_im = imag(fid);
            npts   = numel(fid);
            sw_out = sw;
            cf_mhz = Bfield * 42.577;

        % ----- One pulse, shaped frequency-selective excitation ----------
        %   sim_onepulse_shaped without a gradient (frequency selective,
        %   pulse centred on the simulation centre). The waveform is loaded
        %   for the requested flip angle (numeric type -> flipCyc = angle/360).
        %   args: n, sw, Bfield, lw, pulse_path, tp [ms], flipAngle [deg]
        case 'onepulse_shaped'
            [n, sw, Bfield, lw, pulse_path, tp, flipAngle] = deal(varargin{1:7});
            if ~exist(pulse_path, 'file')
                error('fida_run/onepulse_shaped: pulse waveform not found: "%s"', pulse_path);
            end
            RF  = io_loadRFwaveform(pulse_path, flipAngle, 0);
            out = sim_onepulse_shaped(n, sw, Bfield, lw, sys, RF, tp, 0, 0);
            fid    = out.fids(:);
            fid_re = real(fid); fid_im = imag(fid);
            npts   = numel(fid);
            sw_out = sw;
            cf_mhz = Bfield * 42.577;

        % ----- semi-LASER shaped, 4-step phase cycle -----------------------
        %   sim_semiLASER_shaped_phCyc with FID-A's run_simSemiLASERShaped_phCyc
        %   scheme: [ph1 ph2, ph3 ph4] = [0 0,0 0] - [0 0,0 90] - [0 90,0 0]
        %   + [0 90,0 90], averaged over the spatial grid.
        %   args: n, sw, Bfield, lw, te, pulse_path, tp,
        %         thkX, thkY, fovX, fovY, nX, nY, flipAngle, centreFreq
        case 'semilaser_shaped_phcyc'
            [n, sw, Bfield, lw, te, pulse_path, tp, ...
             thkX, thkY, fovX, fovY, nX, nY, flipAngle, centreFreq] = ...
                deal(varargin{1:15});
            if ~exist(pulse_path, 'file')
                error('fida_run/semilaser_shaped_phcyc: pulse waveform not found: "%s"', pulse_path);
            end
            RF = io_loadRFwaveform(pulse_path, 'ref', 0);
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
            ph1 = [0 0 0 0]; ph2 = [0 0 90 90];
            ph3 = [0 0 0 0]; ph4 = [0 90 0 90];
            sgn = [1 -1 -1 1];
            accumFid = [];
            for ix = 1:nX
                for iy = 1:nY
                    for m = 1:4
                        out = sim_semiLASER_shaped_phCyc(n, sw, Bfield, lw, sys, te, ...
                                  RF, tp, x(ix), y(iy), Gx, Gy, ...
                                  ph1(m), ph2(m), ph3(m), ph4(m), flipAngle, centreFreq);
                        if isempty(accumFid)
                            accumFid = sgn(m) * out.fids(:);
                        else
                            accumFid = accumFid + sgn(m) * out.fids(:);
                        end
                    end
                end
            end
            accumFid = accumFid / (nX * nY * 4);
            fid_re = real(accumFid); fid_im = imag(accumFid);
            npts   = numel(accumFid);
            sw_out = sw;
            cf_mhz = Bfield * 42.577;

        % ----- MEGA-PRESS, shaped refocusing + ideal editing ---------------
        %   sim_megapress_shapedRefoc: refoc phase cycle [0,90]x[0,90] combined
        %   as in run_simMegaPressShaped.m (subtract when exactly one pulse is
        %   at 90), averaged over the spatial grid. Ideal editing flips as in
        %   the 'megapress_ideal' kind (180 within editBand of editPpm).
        %   args: n, sw, Bfield, lw, t1..t5 [ms], editPpm, editBand,
        %         refoc_path, refTp, thkX, thkY, fovX, fovY, nX, nY, edit_on_flag
        case 'megapress_shapedrefoc'
            [n, sw, Bfield, lw, t1, t2, t3, t4, t5, editPpm, editBand, ...
             refoc_path, refTp, thkX, thkY, fovX, fovY, nX, nY, editOn] = ...
                deal(varargin{1:20});
            if ~exist(refoc_path, 'file')
                error('fida_run/megapress_shapedrefoc: refocusing waveform not found: "%s"', refoc_path);
            end
            taus = [t1 t2 t3 t4 t5];
            nsub = numel(sys);
            editFlip = cell(1, nsub);
            for k = 1:nsub
                shifts = sys(k).shifts(:)';
                if editOn
                    editFlip{k} = 180 * double(abs(shifts - editPpm) <= editBand/2);
                else
                    editFlip{k} = zeros(size(shifts));
                end
            end
            refRF = io_loadRFwaveform(refoc_path, 'ref', 0);
            gamma = 42577000;
            if refRF.isGM
                Gx = (refRF.tthk / (refTp/1000)) / thkX;
                Gy = (refRF.tthk / (refTp/1000)) / thkY;
            else
                Gx = (refRF.tbw / (refTp/1000)) / (gamma * thkX / 10000);
                Gy = (refRF.tbw / (refTp/1000)) / (gamma * thkY / 10000);
            end
            [refRF, Gx, Gy] = gm_prepare(refRF, 'megapress_shapedrefoc', Gx, Gy);
            if nX < 2; nX = 2; end
            if nY < 2; nY = 2; end
            x = linspace(-fovX/2, fovX/2, nX);
            y = linspace(-fovY/2, fovY/2, nY);
            refPh = [0 90];
            accumFid = [];
            for ix = 1:nX
                for iy = 1:nY
                    for RP1 = 1:2
                        for RP2 = 1:2
                            out = sim_megapress_shapedRefoc(n, sw, Bfield, lw, taus, sys, ...
                                      editFlip, refRF, refTp, Gx, Gy, x(ix), y(iy), ...
                                      refPh(RP1), refPh(RP2));
                            if xor(RP1 == 2, RP2 == 2); sgn = -1; else; sgn = 1; end
                            if isempty(accumFid)
                                accumFid = sgn * out.fids(:);
                            else
                                accumFid = accumFid + sgn * out.fids(:);
                            end
                        end
                    end
                end
            end
            accumFid = accumFid / (nX * nY * 4);
            fid_re = real(accumFid); fid_im = imag(accumFid);
            npts   = numel(accumFid);
            sw_out = sw;
            cf_mhz = Bfield * 42.577;

        % ----- MEGA-PRESS, fully shaped (refocusing + editing) -------------
        %   sim_megapress_shaped: edit phase cycle [0,90]x[0,90,180,270]
        %   summed, refoc cycle [0,90]x[0,90] combined with signs, averaged
        %   over the spatial grid (8 x 4 x nX x nY simulations per call).
        %   args: n, sw, Bfield, lw, t1..t5 [ms], edit_path, editTp,
        %         editOnFreq, editOffFreq, refoc_path, refTp,
        %         thkX, thkY, fovX, fovY, nX, nY, centreFreq, edit_on_flag
        %   centreFreq is accepted for interface stability but has no effect:
        %   sim_megapress_shaped simulates in a 3 ppm frame (re-referenced below).
        case 'megapress_shaped'
            [n, sw, Bfield, lw, t1, t2, t3, t4, t5, edit_path, editTp, ...
             editOnFreq, editOffFreq, refoc_path, refTp, ...
             thkX, thkY, fovX, fovY, nX, nY, centreFreq, editOn] = ...
                deal(varargin{1:23});
            if ~exist(edit_path, 'file')
                error('fida_run/megapress_shaped: edit pulse waveform not found: "%s"', edit_path);
            end
            if ~exist(refoc_path, 'file')
                error('fida_run/megapress_shaped: refocusing waveform not found: "%s"', refoc_path);
            end
            taus = [t1 t2 t3 t4 t5];
            gamma = 42577000;
            editRF = io_loadRFwaveform(edit_path, 'inv', 0);
            if editOn
                targetFreq = editOnFreq;
            else
                targetFreq = editOffFreq;
            end
            editRFshift = rf_freqshift(editRF, editTp, ...
                (3 - targetFreq) * Bfield * gamma / 1e6);   % this FID-A function's frame is 3 ppm
            refRF = io_loadRFwaveform(refoc_path, 'ref', 0);
            if refRF.isGM
                Gx = (refRF.tthk / (refTp/1000)) / thkX;
                Gy = (refRF.tthk / (refTp/1000)) / thkY;
            else
                Gx = (refRF.tbw / (refTp/1000)) / (gamma * thkX / 10000);
                Gy = (refRF.tbw / (refTp/1000)) / (gamma * thkY / 10000);
            end
            [refRF, Gx, Gy] = gm_prepare(refRF, 'megapress_shaped', Gx, Gy);
            if nX < 2; nX = 2; end
            if nY < 2; nY = 2; end
            x = linspace(-fovX/2, fovX/2, nX);
            y = linspace(-fovY/2, fovY/2, nY);
            editPhCyc1 = [0 90];
            editPhCyc2 = [0 90 180 270];
            refPh = [0 90];
            accumFid = [];
            for ix = 1:nX
                for iy = 1:nY
                    for EP1 = 1:numel(editPhCyc1)
                        for EP2 = 1:numel(editPhCyc2)
                            for RP1 = 1:2
                                for RP2 = 1:2
                                    out = sim_megapress_shaped(n, sw, Bfield, lw, taus, sys, ...
                                              editRFshift, editTp, editPhCyc1(EP1), editPhCyc2(EP2), ...
                                              refRF, refTp, Gx, Gy, x(ix), y(iy), ...
                                              refPh(RP1), refPh(RP2));
                                    if xor(RP1 == 2, RP2 == 2); sgn = -1; else; sgn = 1; end
                                    if isempty(accumFid)
                                        accumFid = sgn * out.fids(:);
                                    else
                                        accumFid = accumFid + sgn * out.fids(:);
                                    end
                                end
                            end
                        end
                    end
                end
            end
            accumFid = accumFid / (nX * nY * numel(editPhCyc1) * numel(editPhCyc2) * 4);
            fid_re = real(accumFid); fid_im = imag(accumFid);
            npts   = numel(accumFid);
            sw_out = sw;
            cf_mhz = Bfield * 42.577;

        % ----- MEGA-SPECIAL shaped (1-D; per run_simMegaSpecialShaped) -
        %   Shaped editing (loaded 'inv', frequency-shifted) + one shaped
        %   1-D refocusing pulse. Edit phase cycles [0,90]x[0,90,180,270]
        %   are summed; the refoc cycle [0,90] is combined subtractively.
        %   args: n, sw, Bfield, lw, t1..t4 [ms], edit_path, editTp,
        %         editOnFreq, editOffFreq, refoc_path, refTp,
        %         thk, fov, npos, centreFreq, edit_on_flag
        %   centreFreq is accepted for interface stability but has no effect:
        %   sim_megaspecial_shaped simulates in a 3 ppm frame (re-referenced below).
        case 'megaspecial_shaped'
            [n, sw, Bfield, lw, t1, t2, t3, t4, edit_path, editTp, ...
             editOnFreq, editOffFreq, refoc_path, refTp, thk, fov, npos, ...
             centreFreq, editOn] = deal(varargin{1:19});
            if ~exist(edit_path, 'file')
                error('fida_run/megaspecial_shaped: edit pulse waveform not found: "%s"', edit_path);
            end
            if ~exist(refoc_path, 'file')
                error('fida_run/megaspecial_shaped: refoc pulse waveform not found: "%s"', refoc_path);
            end
            taus = [t1 t2 t3 t4];
            refRF  = io_loadRFwaveform(refoc_path, 'ref', 0);
            editRF = io_loadRFwaveform(edit_path, 'inv', 0);
            gamma = 42577000;
            if editOn
                targetFreq = editOnFreq;
            else
                targetFreq = editOffFreq;
            end
            editRFshift = rf_freqshift(editRF, editTp, ...
                (3 - targetFreq) * Bfield * gamma / 1e6);   % this FID-A function's frame is 3 ppm
            if refRF.isGM
                G = (refRF.tthk / (refTp/1000)) / thk;
            else
                G = (refRF.tbw / (refTp/1000)) / (gamma * thk / 10000);
            end
            [refRF, G] = gm_prepare(refRF, 'megaspecial_shaped', G);
            if npos < 2; npos = 2; end
            pos = linspace(-fov/2, fov/2, npos);
            editPhCyc1 = [0 90];
            editPhCyc2 = [0 90 180 270];
            refPhCyc   = [0 90];
            accumFid = [];
            for ip = 1:npos
                for EP1 = 1:numel(editPhCyc1)
                    for EP2 = 1:numel(editPhCyc2)
                        out1 = sim_megaspecial_shaped(n, sw, Bfield, lw, ...
                            taus, sys, editRFshift, editTp, ...
                            editPhCyc1(EP1), editPhCyc2(EP2), ...
                            refRF, refTp, G, pos(ip), refPhCyc(1));
                        out2 = sim_megaspecial_shaped(n, sw, Bfield, lw, ...
                            taus, sys, editRFshift, editTp, ...
                            editPhCyc1(EP1), editPhCyc2(EP2), ...
                            refRF, refTp, G, pos(ip), refPhCyc(2));
                        stepFid = out1.fids(:) - out2.fids(:);
                        if isempty(accumFid)
                            accumFid = stepFid;
                        else
                            accumFid = accumFid + stepFid;
                        end
                    end
                end
            end
            accumFid = accumFid / (npos * numel(editPhCyc1) * numel(editPhCyc2));
            fid_re = real(accumFid); fid_im = imag(accumFid);
            npts   = numel(accumFid);
            sw_out = sw;
            cf_mhz = Bfield * 42.577;

        otherwise
            error('fida_run: unknown kind "%s"', kind);
    end

    % ---------- reference the FID to the 4.65 ppm frame -------------------
    % BasisREMY places the returned FID on a ppm axis centred at 4.65. The
    % shaped simulators run in the frame given by 'Sim Centre (ppm)', and
    % three FID-A functions hard-code a 3 ppm frame (sim_megapress_shaped,
    % sim_megapress_shapedRefoc, sim_megaspecial_shaped), so demodulate here.
    switch lower(kind)
        case {'press_shaped', 'semilaser_shaped', 'semilaser_shaped_phcyc', ...
              'steam_shaped', 'megapress_shapededit'}
            sim_centre = centreFreq;
        case {'megaspecial_shaped', 'megapress_shaped', 'megapress_shapedrefoc'}
            sim_centre = 3;
    end
    if abs(sim_centre - 4.65) > 1e-9
        t   = (0:npts-1)' / sw_out;
        fid = (fid_re(:) + 1i * fid_im(:)) .* exp(-1i * 2 * pi * (4.65 - sim_centre) * cf_mhz * t);
        fid_re = real(fid); fid_im = imag(fid);
    end
end

% ---------- helper: gradient-modulated waveforms ----------------------------
function [RF, Gx, Gy] = gm_prepare(RF, kind, Gx, Gy)
% Gradient-modulated (GOIA-type) waveforms carry their own gradient shape in
% column 4, and sim_shapedRF refuses an explicit gradient for them. Scale the
% stored shape to the requested slice thickness (what sim_semiLASER_shaped
% does internally) and hand the simulator Gx = Gy = 0. FID-A's PRESS, STEAM,
% spin-echo and MEGA simulators apply ONE waveform in both directions, so a
% GM pulse needs thkX == thkY there. Conventional pulses pass through.
    if nargin < 4; Gy = Gx; end
    if ~RF.isGM; return; end
    if abs(Gx - Gy) > 1e-9 * max(abs([Gx Gy 1]))
        error(['fida_run/%s: gradient-modulated waveform with thkX ~= thkY: ' ...
               'this simulator applies the same waveform along X and Y, so ' ...
               'use equal thicknesses (or the semi-LASER backend)'], kind);
    end
    RF = rf_scaleGrad(RF, Gx);
    Gx = 0;
    Gy = 0;
end
