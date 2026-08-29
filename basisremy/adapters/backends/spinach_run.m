function [fid_re, fid_im, npts, sw_out, cf_mhz] = spinach_run(metab, kind, n, sw, Bfield, lw, ...
                                                              te, tm, spinach_root, shims_dir, ...
                                                              fida_root)
% SPINACH_RUN  BasisREMY adapter: one FID-A spin system simulated by Spinach under Octave.
%
%   [fid_re, fid_im, npts, sw_out, cf_mhz] = spinach_run(metab, kind, n, sw, Bfield, lw, ...
%                                                        te, tm, spinach_root, shims_dir, fida_root)
%
%   kind   : 'press' | 'spinecho' | 'laser' | 'steam'   (ideal pulses)
%   n, sw  : points and spectral width [Hz]; Bfield [T]; lw [Hz]; te, tm [ms]
%   *_root : Spinach checkout (its kernel/ is added to the path), the BasisREMY shim directory
%            (added in front of it) and the FID-A checkout (metabolite .mat definitions)
%
%   Spinach does the physics (create / basis / hamiltonian / evolution); the pulses are ideal
%   rotations about Spinach's operators. The FID is returned in FID-A's conventions so the two
%   backends compare point by point: frame centred on water (4.65 ppm), amplitude 1 per proton
%   (sim_readout's 2^(2-nspins) scaling, sys.scaleFactor honoured), FID-A's sequence timings
%   (PRESS 90x-TE/4-180y-TE/2-180y-TE/4; spin echo 90x-TE/2-180y-TE/2; LASER 90x + six 180y at
%   TE/6; STEAM 90x-TE/2-(-90x)-TM-90x-TE/2 with the coherence-order filters +1 / 0 / -1), and
%   sim_readout's Lorentzian decay exp(-pi*lw*t). Spinach runs with sys.disable={'hygiene'}
%   (skips its MATLAB-only start-up checks) on top of spinach_shims/ and spinach_octave.patch.

    persistent path_ready
    if isempty(path_ready)
        addpath(genpath(fullfile(spinach_root, 'kernel')));
        addpath(shims_dir, '-begin');          % shims shadow Spinach's MATLAB-only helpers
        warning('off', 'all');
        path_ready = true;
    end

    S = load(fullfile(fida_root, 'simulationTools', 'metabolites', [metab '.mat']));
    f = fieldnames(S);
    parts = S.(f{1});

    n = round(n);
    t = (0:n-1)' / sw;
    fid = zeros(n, 1);
    for p = 1:numel(parts)
        fid = fid + simulate_part(parts(p), kind, n, sw, Bfield, te, tm);
    end
    fid = fid .* exp(-t * pi * lw);            % FID-A: T2 = 1/(pi*linewidth)

    fid_re = real(fid);
    fid_im = imag(fid);
    npts   = n;
    sw_out = sw;
    cf_mhz = Bfield * 42.577;
end


function fid = simulate_part(part, kind, n, sw, Bfield, te, tm)
    shifts = part.shifts(:)';
    nsp    = numel(shifts);
    J = part.J;
    if isempty(J), J = zeros(nsp); end
    scale = 1;
    if isfield(part, 'scaleFactor') && ~isempty(part.scaleFactor), scale = part.scaleFactor; end

    % ----- Spinach spin system: shifts relative to water so the frame sits at 4.65 ppm
    sys = struct();
    sys.magnet   = Bfield;
    sys.isotopes = repmat({'1H'}, 1, nsp);
    sys.output   = 'hush';
    sys.disable  = {'hygiene'};
    inter = struct();
    inter.zeeman.scalar   = num2cell(shifts - 4.65);
    inter.coupling.scalar = num2cell(triu(J, 1));      % FID-A reads J(i,j) for i<j
    bas.formalism     = 'zeeman-hilb';
    bas.approximation = 'none';
    ss = create(sys, inter);
    ss = basis(ss, bas);
    ss = assume(ss, 'nmr');
    H  = full(hamiltonian(ss));
    Lx = full(operator(ss, 'Lx', '1H'));
    Ly = full(operator(ss, 'Ly', '1H'));
    Lz = full(operator(ss, 'Lz', '1H'));
    Lp = operator(ss, 'L+', '1H');

    % coherence order of each element, FID-A's sim_coherenceOrder convention
    M = real(diag(Lz));
    order = M' - M;

    rho  = Lz * scale;                          % FID-A starts from Fz * scaleFactor
    te_s = te / 1000;
    tm_s = tm / 1000;
    switch lower(kind)
        case 'press'
            rho = pulse(rho, Lx, pi/2);
            rho = delay(rho, H, te_s/4);
            rho = pulse(rho, Ly, pi);
            rho = delay(rho, H, te_s/2);
            rho = pulse(rho, Ly, pi);
            rho = delay(rho, H, te_s/4);
        case 'spinecho'
            rho = pulse(rho, Lx, pi/2);
            rho = delay(rho, H, te_s/2);
            rho = pulse(rho, Ly, pi);
            rho = delay(rho, H, te_s/2);
        case 'laser'
            tau = te_s / 6;
            rho = pulse(rho, Lx, pi/2);
            rho = delay(rho, H, tau/2);
            for k = 1:6
                rho = pulse(rho, Ly, pi);
                if k < 6
                    rho = delay(rho, H, tau);
                else
                    rho = delay(rho, H, tau/2);
                end
            end
        case 'steam'
            rho = pulse(rho, Lx, pi/2);
            rho = rho .* (order == 1);
            rho = delay(rho, H, te_s/2);
            rho = pulse(rho, Lx, -pi/2);
            rho = rho .* (order == 0);
            rho = delay(rho, H, tm_s);
            rho = pulse(rho, Lx, pi/2);
            rho = rho .* (order == -1);
            rho = delay(rho, H, te_s/2);
        otherwise
            error('spinach_run: unknown kind "%s"', kind);
    end

    % ----- acquisition by Spinach: trace(coil' * rho(t)) at every dwell. The factor -1i is
    % FID-A's readout convention (trace(rho*F+) with a 90 degree receiver phase); measured
    % against sim_press: identical amplitudes, constant ratio -1i.
    fidp = evolution(ss, H, Lp, rho, 1/sw, n-1, 'observable');
    fid  = -1i * 2^(2-nsp) * fidp(:);
end


function rho = pulse(rho, L, angle)
    P = expm(-1i * angle * L);
    rho = P * rho * P';
end


function rho = delay(rho, H, t)
    P = expm(-1i * H * t);
    rho = P * rho * P';
end
