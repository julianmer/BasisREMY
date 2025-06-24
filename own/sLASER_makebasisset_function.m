function OUT =  sLASER_makebasisset_function(curfolder,pathtofida,system,seq_name,basis_name,B1max,flip_angle,refTp,Npts,sw,lw,Bfield,thkX,thkY,fovX,fovY,nX,nY,te,centreFreq,spinSysList,tau1,tau2,path_to_pulse,path_to_save,path_to_spin_system,display)


% ToolboxCheck

% Where everything is saved
folder_to_save = path_to_save;

%--------------------------------------------------------------------------
% Metabolite Information
%--------------------------------------------------------------------------
path_to_spinsystem = path_to_spin_system;
%--------------------------------------------------------------------------

save_result=true;
%
Waveform=path_to_pulse;
% B1max;% in mikro Tesla // set to [] to be prompted with figure
% flip_angle flip angle;
% refTp %duration of refocusing pulses[ms]
% Npts %number of spectral points
% sw= %spectral width [Hz]
% lw= %linewidth of the output spectrum [Hz]
% Bfield= %Magnetic field strength in [T]
% thkX= %slice thickness of x refocusing pulse [cm]
% thkY= %slice thickness of y refocusing pulse [cm]

% fovX= %size of the full simulation Field of View in the x-direction [cm]
% fovY= %size of the full simulation Field of View in the y-direction [cm]
%
% nX= %Number of grid points to simulate in the x-direction
% nY= %Number of grid points to simulate in the y-direction
%
% full voxel
x=linspace(-fovX/2,fovX/2,nX); %X positions to simulate [cm]
y=linspace(-fovY/2,fovY/2,nY); %Y positions to simulate [cm]
% 
% te %timing of the pulse sequence [ms]
% centreFreq %Centre frequency of MR spectrum [ppm]
%
fovX=-x(1)+x(end);
fovY=-y(1)+y(end);
% select the metabolites you want to simulate 
% the basis set is created with all metabolites in the matfiles_post

% spinSysList= mypar.spinSysList;
% shift
shift_in_ppm=(4.65-centreFreq);

% addpath 
addpath(genpath(pathtofida),'-begin');

% Add the scripts adapted from fidA 
% Will check here first to find the function
addpath(fullfile(curfolder, 'dependencies'), '-begin');

%--------------------------------------------------------------------------
%Load RF waveform
%--------------------------------------------------------------------------
%Niklaus : set inv / macht keinen unterschied oder?
rfPulse=io_loadRFwaveform(Waveform,'inv',0,B1max);
w1max=rfPulse.w1max; % to save the input value 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
sysRef.J=0;
sysRef.shifts=0;
sysRef.scaleFactor=1;
sysRef.name='Ref_0ppm';
%
sysRef.centreFreq=centreFreq;
%
[ ref] = run_mysLASERShaped_fast(rfPulse,refTp,Npts,sw,lw,Bfield,thkX,thkY,x,y,te,sysRef,flip_angle);
%
% tau1=15; tau2=13;%fake timing
refjustforppmrange=sim_press(Npts,sw,Bfield,lw,sysRef,tau1,tau2);
%-------------------------------------------------------------------------
% shift for the ref
%-------------------------------------------------------------------------
% https://ch.mathworks.com/matlabcentral/newsreader/view_thread/243061 

freqShift_hz=shift_in_ppm*(Bfield*42.577478); % in Hz
%--------------------------------------------------------------------------------------
ref.fids=ref.fids.*exp(-(1i*2*pi*freqShift_hz).*ref.t).';


%-------------------------------------------------------------------------
%Load spin systems 
%-------------------------------------------------------------------------
load(path_to_spinsystem)
%-------------------------------------------------------------------------

for met_nr=1:size(spinSysList,2)
    %
    spinSys=spinSysList{met_nr}; %spin system to simulate
    sys=eval(['sys' spinSys]);
    % Schreibe die einfach im ersten rein
    sys(1).centreFreq=centreFreq;
    
    %-------------------------------------------------------------------------
    % Simulation
    %-------------------------------------------------------------------------
    [ out] = run_mysLASERShaped_fast(rfPulse,refTp,Npts,sw,lw,Bfield,thkX,thkY,x,y,te,sys,flip_angle);
    
    % Save before the shift -
    save_out_mat=[folder_to_save,'matfiles_pre'];
    if (exist(save_out_mat,'dir')==0)
                 mkdir(save_out_mat);
    end
    save([save_out_mat,'/',spinSys],'out')

        %-------------------------------------------------------------------------
    % Add shift here and later for every simulated 
    %-------------------------------------------------------------------------
    out.fids=out.fids.*exp(-(1i*2*pi*freqShift_hz).*ref.t).';

    %-------------------------------------------------------------------------
    % add tms ref
    %-------------------------------------------------------------------------
    out=op_addScans(out,ref);
 
    if display == 0
    save_figure=[folder_to_save,'figures'];
   

    if (exist(save_figure,'dir')==0)
                 mkdir(save_figure);
    end
    
    end
    % figure
    if display == 0
    figure;plot(refjustforppmrange.ppm,real(ifftshift(ifft(out.fids))),'b')
    set(gca,'xdir','reverse')
    colormap;set(gcf,'color','w');
    xlim([-1 5])
    xlabel('ppm');
    title(['met with ref ',spinSys])
    print('-dpng','-r300',[save_figure,'\',spinSys])
    end
    %---------------------------------------------------------------------
    % add inforation that then can be stored in BASIS 
    %--------------------------------------------------------------------- 
    out.name=io_sysname({sys.name});
    out.centerFreq=centreFreq; % This is needed for the check within fit_LCMmakeBasis
    out.w1max=w1max;
    
    save_raw=[folder_to_save,'raw'];
    if save_result
        if (exist(save_raw,'dir')==0)
                 mkdir(save_raw);
        end

        RF=io_writelcmraw(out,[save_raw,'/',spinSys,'.RAW'],spinSys);
       
    end
    
    % Saving after shift   
    save_out_mat_end=[folder_to_save,'matfiles_post'];
    if (exist(save_out_mat_end,'dir')==0)
                 mkdir(save_out_mat_end);
    end
    save([save_out_mat_end,'/',spinSys],'out')

    
end

disp('Running fit_makeLCMBasis...');

BASIS=fit_makeLCMBasis(save_out_mat_end, false, [folder_to_save,'/', basis_name],system,seq_name);

OUT = BASIS;

end