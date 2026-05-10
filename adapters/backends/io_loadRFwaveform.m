% io_loadRFwaveform.m  —  BasisREMY headless-Octave override.
%
% This file SHADOWS externals/fidA/inputOutput/io_loadRFwaveform.m. The
% upstream version, when handed a phase-modulated (e.g. adiabatic / GOIA)
% pulse, asks the user to look at a Mz-vs-w1 plot and TYPE in the desired
% w1max:
%
%     [mv,sc]=bes(rf,Tp*1000,'b',0,0,5,40000);
%     plot(sc,mv(3,:));
%     w1max=input('Input desired w1max in kHz (for 5.00 ms pulse):  ');
%
% That is impossible inside our headless Docker Octave — `plot()` aborts
% with "ft_text_renderer: invalid bounding box, cannot render, unable to
% create graphics handle" and `input()` would hang forever anyway.
%
% This patched version replaces ONLY that interactive block with a
% deterministic search for w1max:
%   * 'exc'  / 90°  →  smallest w1 where Mz first crosses 0   (excitation)
%   * 'ref'  / 'inv' / 180°  →  smallest w1 where Mz first reaches the
%                               adiabatic plateau near -1 (or its argmin
%                               if the plateau is never reached)
%   * numeric flip α [deg]   →  smallest w1 where Mz first crosses cos(α)
%
% Everything else is byte-for-byte identical to the upstream file so we
% do not silently change pulse-loading semantics for the non-phase-modulated
% path.
%
% Keep this file on the Octave path AHEAD of externals/fidA/inputOutput/
% (BasisREMY's setup_octave_paths() adds adapters/backends/ AFTER the
% genpath() of externals/fidA/, which under Octave's default prepend
% behaviour means the adapter wins).

function RF_struct=io_loadRFwaveform(filename,type,off_res);

if nargin<3
    off_res=0;
end

if isnumeric(type)
    if type==90;  type='exc'; end
    if type==180; type='inv'; end
end

%Now read in the waveform:
if isstr(filename)
    if exist(filename)
        if filename(end-3:end)=='.pta'
            disp('Siemens format .pta RF pulse file detected!! Loading waveform now.');
            rf=io_readpta(filename);
        elseif filename(end-2:end)=='.RF'
            disp('Varian/Agilent format .RF RF pulse file detected!! Loading waveform now.');
            rf=io_readRF(filename);
        elseif filename(end-3:end)=='.inv'
            disp('Bruker format .inv RF pulse file detected!! Loading waveform now.');
            rf=io_readRFBruk(filename);
        elseif filename(end-3:end)=='.rfc'
            disp('Bruker format .pta RF pulse file detected!! Loading waveform now.');
            rf=io_readRFBruk(filename);
        elseif filename(end-3:end)=='.exc'
            disp('Bruker format .exc RF pulse file detected!! Loading waveform now.');
            rf=io_readRFBruk(filename);
        elseif filename(end-3:end)=='.txt'
            disp('Basic .txt format RF pulse file detected!! Loading waveform now.');
            rf=io_readRFtxt(filename);
        else
            error('ERROR:  RF Pulse file not recognized.  Aborting!');
        end
    else
        error('ERROR:  File not found!  Aborting!');
    end
else
    if ismatrix(filename) && ndims(filename)==2 && ( 2 <= size(filename,2) <= 4 )
        disp('Input is an RF waveform already in the matlab workspace.  Loading waveform now.');
        if size(filename,2)==2
            rf=filename;
            rf(:,3)=ones(length(filename(:,1)),1);
        elseif size(filename,2)==3 || size(filename,2)==4
            rf=filename;
        end
    end
end


Tp=0.005;  %assume a 5 ms rf pulse;

if off_res
    tmp_w1=0.03;
    tmp_min=1;
    switch type
        case 'exc'
            min_floor = 0.2;
        case {'inv','ref'}
            min_floor = -0.98;
    end
    while tmp_min > min_floor
        [mv,sc]=bes(rf,Tp*1000,'f',tmp_w1,-5,5,1000);
        [tmp_min,tmp_min_pos]=min(mv(3,:));
        freq_shift=sc(tmp_min_pos)*1000;
        tmp_w1=tmp_w1+0.02;
    end
    N=size(rf,1);
    dt=Tp/N;
    t=[0:dt:Tp-dt];
    phaseRamp=t*-freq_shift*360;
    rf(:,1)=rf(:,1)+phaseRamp';
end

a=(round(rf(:,1))==180)|(round(rf(:,1))==0);
if sum(a)<length(rf(:,1))
    isPhsMod=true;
else
    isPhsMod=false;
end

jumps=diff(rf(:,1));
jumpsAbs=(abs(jumps)>355 & abs(jumps)<365);
jumpIndex=find(jumpsAbs);
for n=1:length(jumpIndex)
    rf(jumpIndex(n)+1:end,1)=rf(jumpIndex(n)+1:end,1)-(360*(jumps(jumpIndex(n))/abs(jumps(jumpIndex(n)))));
end

%scale amplitude function so that maximum value is 1:
rf(:,2)=rf(:,2)./max(rf(:,2));

if ~isPhsMod
    if isstr(type)
        if type=='exc'
            flipCyc=0.25;
        elseif type=='ref'
            flipCyc=0.5;
        elseif type=='inv'
            flipCyc=0.5;
        end
    elseif isnumeric(type)
        flipCyc=type/360;
    end
    intRF=sum(rf(:,2).*((-2*(rf(:,1)>179))+1))/length(rf(:,2));
    if intRF~=0
        w1max=flipCyc/(intRF*Tp);
    else
        w1max=0;
    end
    tw1=Tp*w1max;
else
    % ----- BasisREMY headless replacement of the interactive block -----
    % Sweep B1max from 0..5 kHz, find the lowest w1 that reaches the
    % desired Mz target. Same Mz-vs-w1 curve the upstream code plots
    % for the user to inspect.
    [mv,sc]=bes(rf,Tp*1000,'b',0,0,5,40000);
    mz = mv(3,:);

    if isstr(type)
        if type=='exc'
            target = 0.0;     % 90°  → Mz crosses zero
            tol    = 0.02;
        else                  % 'ref' / 'inv'
            target = -1.0;    % 180° → Mz dips to / plateaus at -1
            tol    = 0.02;    % accept first arrival in adiabatic plateau
        end
    elseif isnumeric(type)
        target = cos(type*pi/180);
        tol    = 0.02;
    end

    idx = find(mz <= target + tol, 1, 'first');
    if isempty(idx)
        % pulse never reaches the target in [0, 5 kHz] → fall back to
        % the closest match (typical for soft / under-driven pulses).
        [~, idx] = min(abs(mz - target));
    end
    w1max = sc(idx) * 1000;   % sc is in kHz → convert to Hz
    tw1   = Tp * w1max;
    fprintf('io_loadRFwaveform (BasisREMY headless): phase-modulated pulse → auto w1max = %.4f kHz (target Mz = %.3f)\n', ...
            w1max/1000, target);
    % --------------------------------------------------------------------
end

%now it's time to find out the time-bandwidth product:
[mv,sc]=bes(rf,Tp*1000,'f',w1max/1000,-5,5,100000);
if isstr(type)
    if type=='exc'
        index=find(abs(mv(1,:)+1i*mv(2,:))>0.5);
        bw=sc(index(end))-sc(index(1));
    elseif type=='ref'
        index=find(mv(3,:)<0);
        bw=sc(index(end))-sc(index(1));
    elseif type=='inv'
        index=find(mv(3,:)<0);
        bw=sc(index(end))-sc(index(1));
    end
elseif isnumeric(type)
    mz=cos(type);
    thr=(1+mz)/2;
    index=find(mv(3,:)<thr);
    bw=sc(index(end))-sc(index(1));
end

[mv,sc]=bes(rf,Tp*1000,'f',w1max/1000,-bw,bw,100000);
if isstr(type)
    if type=='exc'
        index=find(abs(mv(1,:)+1i*mv(2,:))>0.5);
        bw=sc(index(end))-sc(index(1));
    elseif type=='ref'
        index=find(mv(3,:)<0);
        bw=sc(index(end))-sc(index(1));
    elseif type=='inv'
        index=find(mv(3,:)<0);
        bw=sc(index(end))-sc(index(1));
    end
elseif isnumeric(type)
    mz=cos(type);
    thr=(1+mz)/2;
    index=find(mv(3,:)<thr);
    bw=sc(index(end))-sc(index(1));
end

maxIndex=find(rf(:,2)==max(rf(:,2)));
rfCentre=mean(maxIndex)/length(rf(:,2));

RF_struct.waveform=rf;
RF_struct.type=type;
RF_struct.f0=0;
if size(rf,2)>3 && any(rf(:,4))
    RF_struct.tbw='N/A - gradient modulated pulse';
    RF_struct.isGM=true;
    RF_struct.tthk=bw*Tp;
else
    RF_struct.tbw=bw*Tp*1000;
    RF_struct.isGM=false;
    RF_struct.tthk='N/A - frequency selective pulse';
end
RF_struct.tw1=tw1;
RF_struct.rfCentre=rfCentre;

end

