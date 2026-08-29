function [fid_re, fid_im, npts, sw_out, cf_mhz] = spinach_run(metab, kind, spinach_root, ...
                                                              shims_dir, fida_root, varargin)
% SPINACH_RUN  BasisREMY adapter: one FID-A spin system simulated by Spinach under Octave.
%
%   [fid_re, fid_im, npts, sw_out, cf_mhz] = spinach_run(metab, kind, spinach_root, ...
%                                                        shims_dir, fida_root, args...)
%
%   kind / args:
%     'press' | 'spinecho' | 'laser' | 'steam'   (ideal pulses)
%         n, sw, Bfield, lw, te, tm
%     'press_shaped'      (shaped refocusing pulses on a spatial grid, fast method)
%         n, sw, Bfield, lw, tau1, tau2, pulse_path, tp, thkX, thkY, fovX, fovY, nX, nY, ...
%         flipAngle, centreFreq
%     'semilaser_shaped'  (two AFP pairs on a spatial grid, fast method)
%         n, sw, Bfield, lw, te, pulse_path, tp, thkX, thkY, fovX, fovY, nX, nY, ...
%         flipAngle, centreFreq
%   n points; sw, lw in Hz; Bfield in T; te, tm, tau, tp in ms; thk/fov in cm; centreFreq in ppm.
%   *_root: Spinach checkout (its kernel/ goes on the path), the BasisREMY shim directory (in
%   front of it) and the FID-A checkout (metabolite .mat definitions, rf_scaleGrad).
%
%   Spinach does the physics (create / basis / hamiltonian / shaped_pulse_xy / evolution); the
%   ideal pulses are rotations about Spinach's operators. The FID is returned in FID-A's
%   conventions so the two backends compare point by point: frame referenced to water
%   (4.65 ppm), amplitude 1 per proton (sim_readout's 2^(2-nspins) scaling, sys.scaleFactor
%   honoured), FID-A's sequence timings (PRESS 90x-TE/4-180y-TE/2-180y-TE/4; spin echo
%   90x-TE/2-180y-TE/2; LASER 90x + six 180y at TE/6; STEAM 90x-TE/2-(-90x)-TM-90x-TE/2 with
%   the coherence-order filters +1 / 0 / -1) and sim_readout's Lorentzian decay exp(-pi*lw*t).
%   The shaped kinds follow fida_run.m: FID-A's waveform loader and gradient conventions, the
%   fast spatial method of Zhang et al. 2017 (the X pulse(s) over the x grid, the Y pulse(s)
%   over the y grid, coherence-order filtered), simulated in the frame of 'Sim Centre (ppm)'
%   and demodulated to 4.65 ppm at the end. The shaped pulses themselves are Spinach's
%   piecewise-constant shaped_pulse_xy with the slice gradient as a third (z) control.
%   Spinach runs with sys.disable={'hygiene'} (skips its MATLAB-only start-up checks) on top of
%   spinach_shims/ and spinach_octave.patch.

    persistent path_ready
    if isempty(path_ready)
        % FID-A (io_readpta, rf_scaleGrad, bes, ...) below Spinach's kernel; BasisREMY's own
        % adapter directory above both (its headless io_loadRFwaveform must win over FID-A's),
        % and the shims in front of everything (they shadow Spinach's MATLAB-only helpers)
        addpath(genpath(fida_root));
        addpath(genpath(fullfile(spinach_root, 'kernel')));
        addpath(fileparts(mfilename('fullpath')));
        addpath(shims_dir, '-begin');
        warning('off', 'all');
        path_ready = true;
    end

    % FID-A's combined spin-system file (v7): the per-metabolite .mat files are partly
    % MATLAB v7.3 (GSH, EtOH, Ref0ppm), which Octave cannot read
    metab_dir = fullfile(fida_root, 'simulationTools', 'metabolites');
    S = load(fullfile(metab_dir, 'spinSystems.mat'));
    field = ['sys' metab];
    if isfield(S, field)
        parts = S.(field);
    elseif exist(fullfile(metab_dir, [metab '.mat']), 'file')
        % e.g. Lip: only shipped as its own (v7) file
        S = load(fullfile(metab_dir, [metab '.mat']));
        f = fieldnames(S);
        parts = S.(f{1});
    else
        error('spinach_run: unknown metabolite "%s" (no %s in spinSystems.mat)', metab, field);
    end

    opt = struct('kind', lower(kind));
    switch opt.kind
        case {'press', 'spinecho', 'laser', 'steam'}
            [n, sw, Bfield, lw, opt.te, opt.tm] = deal(varargin{1:6});
            opt.centre = 4.65;
        case 'press_shaped'
            [n, sw, Bfield, lw, opt.tau1, opt.tau2, pulse_path, opt.tp, thkX, thkY, ...
             fovX, fovY, nX, nY, opt.flip, opt.centre] = deal(varargin{1:16});
            opt.RF = load_pulse(pulse_path, opt.kind);
            [opt.RF, opt.Gx, opt.Gy] = pulse_gradients(opt.RF, opt.tp, thkX, thkY, false);
            [opt.x, opt.y] = spatial_grid(fovX, fovY, nX, nY);
            opt.phase = 90;
        case 'semilaser_shaped'
            [n, sw, Bfield, lw, opt.te, pulse_path, opt.tp, thkX, thkY, ...
             fovX, fovY, nX, nY, opt.flip, opt.centre] = deal(varargin{1:15});
            opt.RF = load_pulse(pulse_path, opt.kind);
            [opt.RF, opt.Gx, opt.Gy] = pulse_gradients(opt.RF, opt.tp, thkX, thkY, true);
            [opt.x, opt.y] = spatial_grid(fovX, fovY, nX, nY);
            opt.phase = 0;
        otherwise
            error('spinach_run: unknown kind "%s"', kind);
    end

    n = round(n);
    t = (0:n-1)' / sw;
    fid = zeros(n, 1);
    for p = 1:numel(parts)
        fid = fid + simulate_part(parts(p), opt, n, sw, Bfield);
    end
    fid = fid .* exp(-t * pi * lw);            % FID-A: T2 = 1/(pi*linewidth)

    cf_mhz = Bfield * 42.577;
    if abs(opt.centre - 4.65) > 1e-9          % demodulate to BasisREMY's 4.65 ppm axis
        fid = fid .* exp(-1i * 2 * pi * (4.65 - opt.centre) * cf_mhz * t);
    end

    fid_re = real(fid);
    fid_im = imag(fid);
    npts   = n;
    sw_out = sw;
end


function fid = simulate_part(part, opt, n, sw, Bfield)
    shifts = part.shifts(:)';
    nsp    = numel(shifts);
    J = part.J;
    if isempty(J), J = zeros(nsp); end
    scale = 1;
    if isfield(part, 'scaleFactor') && ~isempty(part.scaleFactor), scale = part.scaleFactor; end

    % ----- Spinach spin system in the rotating frame of opt.centre (ppm)
    sys = struct();
    sys.magnet   = Bfield;
    sys.isotopes = repmat({'1H'}, 1, nsp);
    sys.output   = 'hush';
    sys.disable  = {'hygiene'};
    inter = struct();
    inter.zeeman.scalar   = num2cell(shifts - opt.centre);
    inter.coupling.scalar = num2cell(triu(J, 1));      % FID-A reads J(i,j) for i<j
    bas.formalism     = 'zeeman-hilb';
    bas.approximation = 'none';
    ss = create(sys, inter);
    ss = basis(ss, bas);
    ss = assume(ss, 'nmr');
    op.H  = full(hamiltonian(ss));
    op.Lx = full(operator(ss, 'Lx', '1H'));
    op.Ly = full(operator(ss, 'Ly', '1H'));
    op.Lz = full(operator(ss, 'Lz', '1H'));
    op.ss = ss;

    % coherence order of each element, FID-A's sim_coherenceOrder convention
    M = real(diag(op.Lz));
    op.order = M' - M;

    rho0 = op.Lz * scale;                       % FID-A starts from Fz * scaleFactor
    switch opt.kind
        case 'press'
            te_s = opt.te / 1000;
            rho = pulse(rho0, op.Lx, pi/2);
            rho = delay(rho, op.H, te_s/4);
            rho = pulse(rho, op.Ly, pi);
            rho = delay(rho, op.H, te_s/2);
            rho = pulse(rho, op.Ly, pi);
            rho = delay(rho, op.H, te_s/4);
        case 'spinecho'
            te_s = opt.te / 1000;
            rho = pulse(rho0, op.Lx, pi/2);
            rho = delay(rho, op.H, te_s/2);
            rho = pulse(rho, op.Ly, pi);
            rho = delay(rho, op.H, te_s/2);
        case 'laser'
            tau = opt.te / 1000 / 6;
            rho = pulse(rho0, op.Lx, pi/2);
            rho = delay(rho, op.H, tau/2);
            for k = 1:6
                rho = pulse(rho, op.Ly, pi);
                if k < 6
                    rho = delay(rho, op.H, tau);
                else
                    rho = delay(rho, op.H, tau/2);
                end
            end
        case 'steam'
            te_s = opt.te / 1000;
            rho = pulse(rho0, op.Lx, pi/2);
            rho = rho .* (op.order == 1);
            rho = delay(rho, op.H, te_s/2);
            rho = pulse(rho, op.Lx, -pi/2);
            rho = rho .* (op.order == 0);
            rho = delay(rho, op.H, opt.tm / 1000);
            rho = pulse(rho, op.Lx, pi/2);
            rho = rho .* (op.order == -1);
            rho = delay(rho, op.H, te_s/2);
        case 'press_shaped'
            % fast method: X pulse over the x grid (density matrices summed), Y pulse over
            % the y grid; the average over the nX*nY grid is taken at the end
            d1 = opt.tau1 - opt.tp;
            d2 = opt.tau2 - opt.tp;
            if d1 < 0 || d2 < 0
                error(['spinach_run/press_shaped: the echo times must exceed the ' ...
                       'refocusing pulse duration (Tau 1 = %g, Tau 2 = %g, RefTp = %g ms)'], ...
                      opt.tau1, opt.tau2, opt.tp);
            end
            rho_x = zeros(size(rho0));
            for ix = 1:numel(opt.x)
                rho = pulse(rho0, op.Lx, pi/2);
                rho = rho .* (op.order == -1);
                rho = delay(rho, op.H, d1 / 2000);
                rho = shaped_rf(rho, op, opt.RF, opt.tp, opt.flip, opt.phase, opt.x(ix), opt.Gx);
                rho = rho .* (op.order == 1);
                rho = delay(rho, op.H, (d1 + d2) / 2000);
                rho_x = rho_x + rho;
            end
            fid = zeros(n, 1);
            for iy = 1:numel(opt.y)
                rho = shaped_rf(rho_x, op, opt.RF, opt.tp, opt.flip, opt.phase, opt.y(iy), opt.Gy);
                rho = rho .* (op.order == -1);
                rho = delay(rho, op.H, d2 / 2000);
                fid = fid + readout(rho, op, n, sw, nsp);
            end
            fid = fid / (numel(opt.x) * numel(opt.y));
            cleanup_scratch(ss);
            return
        case 'semilaser_shaped'
            if opt.te / 4 < opt.tp
                error(['spinach_run/semilaser_shaped: the refocusing pulse (%g ms) cannot ' ...
                       'be longer than a quarter of TE (%g ms)'], opt.tp, opt.te);
            end
            tau1 = (opt.te / 4 - opt.tp) / 2 / 1000;
            tau2 = (opt.te / 4 - opt.tp) / 1000;
            rho_x = zeros(size(rho0));
            for ix = 1:numel(opt.x)
                rho = pulse(rho0, op.Lx, pi/2);
                rho = rho .* (op.order == -1);
                rho = delay(rho, op.H, tau1);
                rho = shaped_rf(rho, op, opt.RF, opt.tp, opt.flip, opt.phase, opt.x(ix), opt.Gx);
                rho = rho .* (op.order == 1);
                rho = delay(rho, op.H, tau2);
                rho = shaped_rf(rho, op, opt.RF, opt.tp, opt.flip, opt.phase, opt.x(ix), opt.Gx);
                rho = rho .* (op.order == -1);
                rho = delay(rho, op.H, tau2);
                rho_x = rho_x + rho;
            end
            fid = zeros(n, 1);
            for iy = 1:numel(opt.y)
                rho = shaped_rf(rho_x, op, opt.RF, opt.tp, opt.flip, opt.phase, opt.y(iy), opt.Gy);
                rho = rho .* (op.order == 1);
                rho = delay(rho, op.H, tau2);
                rho = shaped_rf(rho, op, opt.RF, opt.tp, opt.flip, opt.phase, opt.y(iy), opt.Gy);
                rho = rho .* (op.order == -1);
                rho = delay(rho, op.H, tau1);
                fid = fid + readout(rho, op, n, sw, nsp);
            end
            fid = fid / (numel(opt.x) * numel(opt.y));
            cleanup_scratch(ss);
            return
    end
    fid = readout(rho, op, n, sw, nsp);
    cleanup_scratch(ss);
end


function cleanup_scratch(ss)
    % Spinach creates a scratch directory per create() call under its own tree and relies
    % on its 'hygiene' pass (disabled here) to remove it — do it ourselves.
    if isfield(ss.sys, 'scratch') && exist(ss.sys.scratch, 'dir')
        rmdir(ss.sys.scratch, 's');
    end
end


function fid = readout(rho, op, n, sw, nsp)
    % acquisition by Spinach: trace(coil' * rho(t)) at every dwell, coil = L+. The factor -1i
    % is FID-A's readout convention (trace(rho*F+) with a 90 degree receiver phase); measured
    % against sim_press: identical amplitudes, constant ratio -1i.
    Lp   = operator(op.ss, 'L+', '1H');
    fidp = evolution(op.ss, op.H, Lp, rho, 1/sw, n-1, 'observable');
    fid  = -1i * 2^(2-nsp) * fidp(:);
end


function rho = shaped_rf(rho, op, RF, tp, flipAngle, phase, x, G)
    % One shaped pulse at position x [cm] under slice gradient G [G/cm] (or the waveform's own
    % gradient shape for gradient-modulated pulses), as Spinach's piecewise-constant
    % shaped_pulse_xy. Amplitude calibration follows FID-A's sim_shapedRF: w1max [kHz] =
    % tw1 * flipAngle / (tp * 180), the waveform scaled to that maximum; the slice gradient
    % enters as a z control with the sign of the chemical-shift term of the Hamiltonian.
    gamma = 42577000;                                   % Hz/T, as in sim_shapedRF
    wf    = RF.waveform;
    w1max = RF.tw1 * flipAngle / (tp * 180) * 1000;     % Hz
    amp   = wf(:, 2) / max(wf(:, 2)) * w1max;           % Hz
    zeta  = (wf(:, 1) + phase) * pi / 180;
    dts   = wf(:, 3) * (tp / 1000) / sum(wf(:, 3));     % s
    if RF.isGM
        grad = wf(:, 4);                                % G/cm per step (scaled by rf_scaleGrad)
    else
        grad = G * ones(size(amp));
    end
    Cx = 2 * pi * amp .* cos(zeta);
    Cy = 2 * pi * amp .* sin(zeta);
    Cz = -2 * pi * gamma * grad * x * 1e-4;             % G/cm * cm -> T, Larmor sign convention
    rho = shaped_pulse_xy(op.ss, op.H, {op.Lx, op.Ly, op.Lz}, {Cx, Cy, Cz}, dts, rho, 'expv-pwc');
end


function RF = load_pulse(pulse_path, kind)
    if ~exist(pulse_path, 'file')
        error('spinach_run/%s: pulse waveform not found: "%s"', kind, pulse_path);
    end
    RF = io_loadRFwaveform(pulse_path, 'ref', 0);      % BasisREMY's headless FID-A loader
end


function [RF, Gx, Gy] = pulse_gradients(RF, tp, thkX, thkY, per_axis)
    % Slice gradients [G/cm] from the pulse's time-bandwidth (or time*thickness) product,
    % fida_run.m conventions. Gradient-modulated waveforms carry their own gradient shape,
    % which is scaled to the requested thickness; per_axis allows different X/Y thicknesses
    % (semi-LASER keeps two copies), otherwise thkX == thkY is required.
    if RF.isGM
        Gx = (RF.tthk / (tp/1000)) / thkX;
        Gy = (RF.tthk / (tp/1000)) / thkY;
        if ~per_axis && abs(Gx - Gy) > 1e-9 * max(abs([Gx Gy 1]))
            error(['spinach_run: gradient-modulated waveform with thkX ~= thkY: this ' ...
                   'kind applies the same waveform along X and Y, use equal thicknesses']);
        end
        RF = rf_scaleGrad(RF, Gx);
        Gx = 0;
        Gy = 0;
    else
        Gx = (RF.tbw / (tp/1000)) / (4258 * thkX);
        Gy = (RF.tbw / (tp/1000)) / (4258 * thkY);
    end
end


function [x, y] = spatial_grid(fovX, fovY, nX, nY)
    % FID-A convention: linspace(-fov/2, fov/2, n) — with n = 2 the points sit on the edges
    nX = max(2, round(nX));
    nY = max(2, round(nY));
    x = linspace(-fovX/2, fovX/2, nX);
    y = linspace(-fovY/2, fovY/2, nY);
end


function rho = pulse(rho, L, angle)
    P = expm(-1i * angle * L);
    rho = P * rho * P';
end


function rho = delay(rho, H, t)
    P = expm(-1i * H * t);
    rho = P * rho * P';
end
