% sLASER_makebasisset_function.m
% -------------------------------------------------------------------------
%
% Generate a basis set for sLASER MRS simulations using custom RF pulses and
% metabolite spin systems. Simulates reference and metabolite signals, applies
% frequency shifts, saves intermediate results, and creates a basis set for fitting.
%
% AUTHORS:
%   Mahrshi Jani
%   Julian P. Merkofer (j.p.merkofer@tue.nl)
%
% USAGE:
% outputs = sLASER_makebasisset_function(curfolder, pathtofida, system, seq_name, basis_name, ...
%     B1max, flip_angle, refTp, Npts, sw, lw, Bfield, thkX, thkY, fovX, fovY, ...
%     nX, nY, te, centreFreq, spinSysList, tau1, tau2, path_to_pulse, path_to_save, ...
%     path_to_spin_system, display, make_basis, make_raw);
%
% INPUTS:
%   curfolder           - Current working directory
%   pathtofida          - Path to fidA toolbox
%   system              - System name (string)
%   seq_name            - Sequence name (string)
%   basis_name          - Output basis set name (string)
%   B1max               - Maximum B1 field [uT]
%   flip_angle          - Flip angle [deg]
%   refTp               - Refocusing pulse duration [ms]
%   Npts                - Number of spectral points
%   sw                  - Spectral width [Hz]
%   lw                  - Linewidth [Hz]
%   Bfield              - Magnetic field strength [T]
%   thkX, thkY          - Slice thickness [cm]
%   fovX, fovY          - Field of view [cm]
%   nX, nY              - Grid points in x/y
%   te                  - Echo time [ms]
%   centreFreq          - Center frequency [ppm]
%   spinSysList         - Cell array of metabolite spin system names
%   tau1, tau2          - PRESS sequence timings [ms]
%   path_to_pulse       - Path to RF pulse file
%   path_to_save        - Output directory
%   path_to_spin_system - Path to spin system .mat file
%   display             - Display flag ('y'/'n')
%   make_basis          - Flag to create .basis file ('y'/'n')
%   make_raw            - Flag to save .raw data ('y'/'n')
%
% OUTPUTS:
%   outputs:            - Cell array of simulated metabolite outputs
%
% -------------------------------------------------------------------------

function outputs = sLASER_makebasisset_function(curfolder, pathtofida, system, seq_name, ...
    basis_name, B1max, flip_angle, refTp, Npts, sw, lw, Bfield, thkX, thkY, fovX, fovY, ...
    nX, nY, te, centreFreq, spinSysList, tau1, tau2, path_to_pulse, path_to_save, ...
    path_to_spin_system, display, make_basis, make_raw)

    folder_to_save = path_to_save;
    path_to_spinsystem = path_to_spin_system;
    save_result = true;
    Waveform = path_to_pulse;

    % Simulation grid
    x = linspace(-fovX/2, fovX/2, nX);
    y = linspace(-fovY/2, fovY/2, nY);
    fovX = -x(1) + x(end);
    fovY = -y(1) + y(end);

    shift_in_ppm = (4.65 - centreFreq);

    % Add required paths
    addpath(genpath(pathtofida), '-begin');
    addpath(fullfile(curfolder, 'dependencies'), '-begin');

    % Load RF waveform
    rfPulse = io_loadRFwaveform(Waveform, 'inv', 0, B1max);
    w1max = rfPulse.w1max;

    % Reference system setup
    sysRef.J = 0;
    sysRef.shifts = 0;
    sysRef.scaleFactor = 1;
    sysRef.name = 'Ref_0ppm';
    sysRef.centreFreq = centreFreq;

    ref = run_mysLASERShaped_fast(rfPulse, refTp, Npts, sw, lw, Bfield, thkX, thkY, x, y, te, sysRef, flip_angle);
    refjustforppmrange = sim_press(Npts, sw, Bfield, lw, sysRef, tau1, tau2);

    % Frequency shift for reference
    freqShift_hz = shift_in_ppm * (Bfield * 42.577478);
    ref.fids = ref.fids .* exp(-(1i * 2 * pi * freqShift_hz) .* ref.t).';

    % Load spin systems
    load(path_to_spinsystem);

    outputs = cell(1, size(spinSysList, 2));
    for met_nr = 1:size(spinSysList, 2)
        spinSys = spinSysList{met_nr};
        sys = eval(['sys' spinSys]);
        sys(1).centreFreq = centreFreq;

        % Simulate metabolite
        out = run_mysLASERShaped_fast(rfPulse, refTp, Npts, sw, lw, Bfield, thkX, thkY, x, y, te, sys, flip_angle);

        % Save before shift
        save_out_mat = [folder_to_save, 'matfiles_pre'];
        if exist(save_out_mat, 'dir') == 0
            mkdir(save_out_mat);
        end
        save([save_out_mat, '/', spinSys], 'out');

        % Apply frequency shift
        out.fids = out.fids .* exp(-(1i * 2 * pi * freqShift_hz) .* ref.t).';

        % Add reference scan
        out = op_addScans(out, ref);

        % Create folder for figures
        if display == 1
            save_figure=[folder_to_save, 'figures'];

            if (exist(save_figure, 'dir')==0)
                mkdir(save_figure);
            end
        end

        % Display metabolite with reference
        if display == 'y' || display == 'Y'
            figure; plot(refjustforppmrange.ppm, real(ifftshift(ifft(out.fids))), 'b')
            set(gca, 'xdir', 'reverse')
            colormap; set(gcf, 'color', 'w');
            xlim([-1 5])
            xlabel('ppm');
            title(['met with ref ', spinSys])
            print('-dpng', '-r300', [save_figure, '\', spinSys])
        end

        % Add metadata
        out.name = io_sysname({sys.name});
        out.centerFreq = centreFreq;
        out.w1max = w1max;

        % Save raw data
        save_raw = [folder_to_save, 'raw'];
        if save_result
            if exist(save_raw, 'dir') == 0
                mkdir(save_raw);
            end
            RF = io_writelcmraw(out, [save_raw, '/', spinSys, '.RAW'], spinSys);
        end

        % Save after shift
        save_out_mat_end = [folder_to_save, 'matfiles_post'];
        if exist(save_out_mat_end, 'dir') == 0
            mkdir(save_out_mat_end);
        end
        save([save_out_mat_end, '/', spinSys], 'out');

        % Create .raw files
        if make_raw == 'y'|| make_raw == 'Y'
            metab = spinSys;
            RF = io_writelcmraw(out, [path metab '.RAW'], metab);
        elseif make_raw == 'n'||make_raw == 'N'
            RF = []
        end

        outputs{met_nr} = out;
    end

    % Create .basis file
    if make_basis == 'y' || make_basis == 'Y'
        disp('Running fit_makeLCMBasis...');
        BASIS = fit_makeLCMBasis(save_out_mat_end, false, [folder_to_save, '/', basis_name], system, seq_name);
    end
end