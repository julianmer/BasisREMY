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
        % Shells out to the BasisREMY-patched sim_lcmrawbasis (which is
        % itself an override on the Octave path). We replicate just the
        % per-metab core here to avoid file I/O — sim_lcmrawbasis already
        % loads spinSystems.mat internally, so we delegate fully.
        %   args: n, sw, Bfield, lw, tau1, tau2, seq, out_path
        case 'ideal'
            [n, sw, Bfield, lw, tau1, tau2, seq, out_path] = ...
                deal(varargin{1:8});
            results = sim_lcmrawbasis(n, sw, Bfield, lw, metab, ...
                                      tau1, tau2, 'n', 'y', seq, out_path);
            fid    = results(:,1) + 1i*results(:,2);
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
            RF = io_loadRFwaveform(pulse_path, 'refoc', 0);
            % gradients chosen so RF bandwidth ≈ slice thickness
            Gx = (RF.tbw / (tp/1000)) / (4258 * thkX);
            Gy = (RF.tbw / (tp/1000)) / (4258 * thkY);
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

