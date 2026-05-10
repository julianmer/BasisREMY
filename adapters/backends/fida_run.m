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

        % ----- stubs ---------------------------------------------------
        % The Python side already raises NotImplementedError for these
        % kinds before reaching Octave, but we guard here as well so a
        % buggy caller gets a clear MATLAB-side error too.
        case {'semilaser_shaped','steam_shaped','spinecho_shaped', ...
              'megapress_shaped','megaspecial_shaped','laser', ...
              'megapress_ideal','spinecho_xn','onepulse'}
            error('fida_run: kind "%s" is a registered FID-A wrapper but the Octave-side branch is not implemented yet.', kind);

        otherwise
            error('fida_run: unknown kind "%s"', kind);
    end
end

