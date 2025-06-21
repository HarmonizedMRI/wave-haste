% -------------------------------------------------------------------------
%% Create a HASTE sequence and export for execution
%% relies on the original version created by Maxim Zaitsev & Juergen Hennig
% -------------------------------------------------------------------------
% 
% The |Sequence| class provides functionality to create magnetic
% resonance sequences (MRI or NMR) from basic building blocks.
%
% This provides an implementation of the open file format for MR sequences
% described here: http://pulseq.github.io/specification.pdf
%
% This example performs the following steps:
% 
% # Create slice selective RF pulse for imaging.
% # Create readout gradient and phase encode strategy.
% # Loop through phase encoding and generate sequence blocks.
% # Write the sequence to an open file format suitable for execution on a
% scanner.
% 
%   Juergen Hennig <juergen.hennig@uniklinik-freiburg.de>
%   Maxim Zaitsev  <maxim.zaitsev@uniklinik-freiburg.de>
 
% -------------------------------------------------------------------------
%% Instantiation and gradient limits
% The system gradient limits can be specified in various units _mT/m_,
% _Hz/cm_, or _Hz/m_. However the limits will be stored internally in units
% of _Hz/m_ for amplitude and _Hz/m/s_ for slew. Unspecificied hardware
% parameters will be assigned default values.
% -------------------------------------------------------------------------

clear

addpath utils_seq_haste/

dG = 250e-6;  % 'standard' ramp time - makes sequence structure much simpler

gMax = 45;
sMax = 140;  % to avoid PNS

system = mr.opts('MaxGrad', gMax, 'GradUnit', 'mT/m', ...
    'MaxSlew', sMax, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 100e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

seq = mr.Sequence(system);


% -------------------------------------------------------------------------
%% Sequence events
% Some sequence parameters are defined using standard MATLAB variables
% -------------------------------------------------------------------------

use_wave_y = 1;         % set to 1 to use sine wave in Gy

fov = 240e-3;
Ny_pre = 36;            % lines of k-space prior to TE

Nx = 120; 
Ny = 120; 

Ry = 3;                 % ky acceleration factor

necho = Ny/2 + Ny_pre;  % ETL
necho = necho / Ry;     % ETL shortened due to acceleration

disp(['ETL = ', num2str(necho)])

Nslices = 1;

rflip=180;

if (numel(rflip)==1)
    rflip=rflip+zeros([1 necho]); 
end

sliceThickness=5e-3;
TE = 6.6e-3;            % controls echo spacing -> we want this to be as short as possible to reduce T2 blurring
TR = 2000e-3;
% TEeff=60e-3;


deltak = 1/fov;
kWidth = Nx*deltak;
nex = 1;

PEorder = ((-Ny_pre):Ry:Ny/2-1)';

disp(['PE order = ', num2str(PEorder.')])

pf_factor = (length(find(PEorder<0)) + length(find(PEorder>0))) / 2 / length(find(PEorder>0))

phaseAreas = PEorder*deltak;

idx = find(PEorder == 0);
TEeff = TE * idx;

disp(['TE eff = ', num2str(TEeff*1e3), ' ms'])

k0 = round(TEeff/TE);
PEtype = 'linear';

assert(length(PEorder) == necho)


% -------------------------------------------------------------------------
%% Define wave params
% -------------------------------------------------------------------------


os_factor = 5;
dwell           =   5e-6;
Tread           =   dwell * Nx * os_factor;


num_read_points = Nx * os_factor;

time_per_adc_point = Tread / num_read_points;


gro = mr.makeTrapezoid('x',system,'FlatArea',kWidth,'FlatTime',Tread,'riseTime',dG);


gYmax_wave = 5;    % mT/m
sYmax_wave = sMax; % T/m/s

num_cycles = 10;    % number of wave cycles
                    % make sure this divides the number of read points

assert(rem(num_read_points,num_cycles)==0)

waveParams = [];
waveParams.ADC_duration = Tread * 1e6;  % ADC duration in micro sec

waveParams.NcyclesGy = num_cycles;    % no of sinusiodal cycles in Gy

waveParams.sYmax = sYmax_wave*100;     % max slew rate in ky axis: Gauss / cm / sec

waveParams.gYmax = gYmax_wave/10;       % max gradient amplitude in kz axis: Gauss / cm


%--------------------------------------------------------------------------
%% make Gy wave gradient: sine
% consider oversampling  
%--------------------------------------------------------------------------

% Tread is equal to gx.flatTime. add delay so that Gy wave starts at the
% beginning of Gx flat top

wavepoints = round(Tread / system.gradRasterTime);     % number of readout points
T_wavepoints = system.gradRasterTime;                  % grad raster time in seconds

TimePerSineY = (Tread / system.gradRasterTime) * T_wavepoints / waveParams.NcyclesGy;     % time per sine wave in seconds

wY = 2*pi / TimePerSineY;

if waveParams.sYmax >= wY * waveParams.gYmax
    disp('wave amplitude is not slew limited')
    G0_Y = waveParams.gYmax;
else
    disp('wave amplitude is slew limited')
    G0_Y = waveParams.sYmax / wY;
end

SO_Y = G0_Y * wY;

% add one last point to go back to zero Gy (wavepoints+1 total)
% this seems to be done in the siemens code as well
GradTimePointsY = [0:wavepoints] * T_wavepoints;        

scaling_factor = system.gamma * 1e-2;                          % convert from G/cm to T/m by multiplying by 1e-2, then to Hz/m by multiplying by gamma

GwaveY = G0_Y * sin(wY * GradTimePointsY) * scaling_factor; % Gy amplitude in Hz/m


disp(['Gx amplitude: ', num2str(max(abs(gro.amplitude)*1e-3)), ' kHz/m'])
disp(['Gy amplitude: ', num2str(max(abs(GwaveY)*1e-3)), ' kHz/m'])


%--------------------------------------------------------------------------
% pad Gy waveform at the beginning and end with zeroes to account for rise and fall time of trapezoid
% SF: This is not needed here because of G5 and G7 
%--------------------------------------------------------------------------

%num_pad_pre = round(gro.riseTime / sys.gradRasterTime);
num_pad_pre = round(system.adcDeadTime / system.gradRasterTime);

% -1 to account for added point at the end to make sure sine wave goes back
% to zero due to [0:wavepoints] instead of [0:wavepoints-1]
% num_pad_post = round(gro.fallTime / sys.gradRasterTime) - 1;     
num_pad_post = round(system.adcDeadTime / system.gradRasterTime) - 1;     


GwaveY = cat(2, zeros(1,num_pad_pre), GwaveY, zeros(1,num_pad_post));

disp(['num points Gy:', num2str(length(GwaveY))])


% form the Gx waveform for display
GwaveX = gro.amplitude * ones(1, round(gro.flatTime / system.gradRasterTime));

GX_rampup = linspace(0, gro.amplitude, round(gro.riseTime / system.gradRasterTime));
GX_rampdown = linspace(gro.amplitude, 0, round(gro.fallTime / system.gradRasterTime));

GwaveX = cat(2, GX_rampup, GwaveX, GX_rampdown);

disp(['num points Gx:', num2str(length(GwaveX))])



% form the Gy object -> play on physical y axis
gy_wave = mr.makeArbitraryGrad('y', GwaveY);


%--------------------------------------------------------------------------
%%
%--------------------------------------------------------------------------

samplingTime= 6.4e-3;           % readout duration

readoutTime = Tread + 2*system.adcDeadTime;

disp(num2str(readoutTime*1e3))

tEx=2.5e-3;     % duration of excitation pulse
tExwd=tEx+system.rfRingdownTime+system.rfDeadTime;

tRef=2e-3;      % duration of refocusing pulse
tRefwd=tRef+system.rfRingdownTime+system.rfDeadTime;

tSp=0.5*(TE-readoutTime-tRefwd);
tSpex=0.5*(TE-tExwd-tRefwd);

fspR=1.0;
fspS=0.5;

rfex_phase=pi/2; % MZ: we need to maintain these as variables because we will overwrtite phase offsets for multiple slice positions
rfref_phase=0;


%--------------------------------------------------------------------------
%%
%--------------------------------------------------------------------------

%%% Base gradients
%%% Slice selection
% Key concepts in the sequence description are *blocks* and *events*.
% Blocks describe a group of events that are executed simultaneously. This
% hierarchical structure means that one event can be used in multiple
% blocks, a common occurrence in MR sequences, particularly in imaging
% sequences. 
%
% First, the slice selective RF pulses (and corresponding slice gradient)
% are generated using the |makeSincPulse| function.
% Gradients are recalculated such that their flattime covers the pulse plus
% the rfdead- and rfringdown- times.
%
flipex=90*pi/180;
[rfex, gz] = mr.makeSincPulse(flipex,system,'Duration',tEx,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,'PhaseOffset',rfex_phase,...
    'use', 'excitation');

GSex = mr.makeTrapezoid('z',system,'amplitude',gz.amplitude,'FlatTime',tExwd,'riseTime',dG);
% plotPulse(rfex,GSex);

flipref=rflip(1)*pi/180;

[rfref, gz] = mr.makeSincPulse(flipref,system,'Duration',tRef,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,'PhaseOffset',rfref_phase,'use','refocusing');

GSref = mr.makeTrapezoid('z',system,'amplitude',GSex.amplitude,'FlatTime',tRefwd,'riseTime',dG);
% plotPulse(rfref,GSref);

AGSex=GSex.area/2;
GSspr = mr.makeTrapezoid('z',system,'area',AGSex*(1+fspS),'duration',tSp,'riseTime',dG);
GSspex = mr.makeTrapezoid('z',system,'area',AGSex*fspS,'duration',tSpex,'riseTime',dG);


%--------------------------------------------------------------------------
%%
%--------------------------------------------------------------------------

%%% Readout gradient
% To define the remaining encoding gradients we need to calculate the
% $k$-space sampling. The Fourier relationship
%
% $$\Delta k = \frac{1}{FOV}$$
% 
% Therefore the area of the readout gradient is $n\Delta k$.


% GRacq = mr.makeTrapezoid('x',system,'FlatArea',kWidth,'FlatTime',readoutTime,'riseTime',dG);
GRacq = gro;

% adc = mr.makeAdc(Nx,'Duration',samplingTime, 'Delay', system.adcDeadTime);%,'Delay',GRacq.riseTime);
adc = mr.makeAdc(Nx*os_factor,'Duration',gro.flatTime, 'Delay', system.adcDeadTime);%,'Delay',GRacq.riseTime);

GRspr = mr.makeTrapezoid('x',system,'area',GRacq.area*fspR,'duration',tSp,'riseTime',dG);
GRspex = mr.makeTrapezoid('x',system,'area',GRacq.area*(1+fspR),'duration',tSpex,'riseTime',dG);


AGRspr=GRspr.area;%GRacq.area/2*fspR;
AGRpreph = GRacq.area/2+AGRspr;%GRacq.area*(1+fspR)/2;
GRpreph = mr.makeTrapezoid('x',system,'Area',AGRpreph,'duration',tSpex,'riseTime',dG);


%--------------------------------------------------------------------------
%%
%--------------------------------------------------------------------------

%%% Phase encoding
% To move the $k$-space trajectory away from 0 prior to the readout a
% prephasing gradient must be used. Furthermore rephasing of the slice
% select gradient is required.

%[PEorder,Ny] = myTSE_PEorder(Ny,necho,k0,PEtype);
%nex=floor(Ny/necho);
%pe_steps=(1:(necho*nex))-0.5*necho*nex-1;
%if 0==mod(necho,2)
%    pe_steps=circshift(pe_steps,[0,-round(nex/2)]); % for odd number of echoes we have to apply a shift to avoid a contrast jump at k=0
%end
%PEorder=reshape(pe_steps,[nex,necho])';
% nex=1;
% 
% % PEorder=((-Ny_pre):Ny)';
% PEorder=((-Ny_pre):Ry:Ny/2-1)';
% 
% phaseAreas = PEorder*deltak;


%--------------------------------------------------------------------------
% split gradients and recombine into blocks
%--------------------------------------------------------------------------

% lets start with slice selection....
GS1times=[0 GSex.riseTime];
GS1amp=[0 GSex.amplitude];
GS1 = mr.makeExtendedTrapezoid('z','times',GS1times,'amplitudes',GS1amp);

GS2times=[0 GSex.flatTime];
GS2amp=[GSex.amplitude GSex.amplitude];
GS2 = mr.makeExtendedTrapezoid('z','times',GS2times,'amplitudes',GS2amp);

GS3times=[0 GSspex.riseTime GSspex.riseTime+GSspex.flatTime GSspex.riseTime+GSspex.flatTime+GSspex.fallTime];
GS3amp=[GSex.amplitude GSspex.amplitude GSspex.amplitude GSref.amplitude];
GS3 = mr.makeExtendedTrapezoid('z','times',GS3times,'amplitudes',GS3amp);

GS4times=[0 GSref.flatTime];
GS4amp=[GSref.amplitude GSref.amplitude];
GS4 = mr.makeExtendedTrapezoid('z','times',GS4times,'amplitudes',GS4amp);

GS5times=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
GS5amp=[GSref.amplitude GSspr.amplitude GSspr.amplitude 0];
GS5 = mr.makeExtendedTrapezoid('z','times',GS5times,'amplitudes',GS5amp);

GS7times=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
GS7amp=[0 GSspr.amplitude GSspr.amplitude GSref.amplitude];
GS7 = mr.makeExtendedTrapezoid('z','times',GS7times,'amplitudes',GS7amp);

% and now the readout gradient....

GR3=GRpreph;%GRspex;

GR5times=[0 GRspr.riseTime GRspr.riseTime+GRspr.flatTime GRspr.riseTime+GRspr.flatTime+GRspr.fallTime];
GR5amp=[0 GRspr.amplitude GRspr.amplitude GRacq.amplitude];
GR5 = mr.makeExtendedTrapezoid('x','times',GR5times,'amplitudes',GR5amp);

GR6times=[0 readoutTime];
GR6amp=[GRacq.amplitude GRacq.amplitude];
GR6 = mr.makeExtendedTrapezoid('x','times',GR6times,'amplitudes',GR6amp);

GR7times=[0 GRspr.riseTime GRspr.riseTime+GRspr.flatTime GRspr.riseTime+GRspr.flatTime+GRspr.fallTime];
GR7amp=[GRacq.amplitude GRspr.amplitude GRspr.amplitude 0];
GR7 = mr.makeExtendedTrapezoid('x','times',GR7times,'amplitudes',GR7amp);


% and filltimes
tex=GS1.shape_dur+GS2.shape_dur+GS3.shape_dur;
tref=GS4.shape_dur+GS5.shape_dur+GS7.shape_dur+readoutTime;
tend=GS4.shape_dur+GS5.shape_dur;

tETrain=tex+necho*tref+tend;
TRfill=(TR-Nslices*tETrain)/Nslices;

% round to gradient raster
TRfill=system.gradRasterTime * round(TRfill / system.gradRasterTime);

if TRfill<0, TRfill=1e-3; 
    disp(strcat('Warning!!! TR too short, adapted to include all slices to : ',num2str(1000*Nslices*(tETrain+TRfill)),' ms')); 
else
    disp(strcat('TRfill : ',num2str(1000*TRfill),' ms')); 
end

delayTR = mr.makeDelay(TRfill);

% extra delay not present in TSE:
delayEnd = mr.makeDelay(5);


%--------------------------------------------------------------------------
%% include calibration lines after the HASTE acquisition using GRE for coil
%% sensitivity estimation
%% calculate timing: ACS
%--------------------------------------------------------------------------

Tpre = 3e-3;            % duration of Gz/Gy gradient blips
Ndummy = 10;

params.Ny_acs = 32;                 % external gre acs size ky

params.TE_acs = 7e-3;               % ACS TE
params.TR_acs = 14e-3;              % ACS TR

Ny_acs = params.Ny_acs;

TE_acs = params.TE_acs;
TR_acs = params.TR_acs;

flipex_acs=15*pi/180;

[rf_acs, gz_acs] = mr.makeSincPulse(flipex,system,'Duration',tEx,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,'PhaseOffset',rfex_phase,...
    'use', 'excitation');

gro_acs = mr.makeTrapezoid('x',system,'FlatArea',kWidth,'FlatTime',Tread,'riseTime',dG);

adc_acs = mr.makeAdc(Nx*os_factor,'Duration',gro.flatTime, 'Delay', gro_acs.riseTime); 

gxPre = mr.makeTrapezoid('x',system,'Area',-gro_acs.area/2,'Duration',Tpre);        % Gx prewinder

gxSpoil = mr.makeTrapezoid('x',system,'Area',gro_acs.area,'Duration',Tpre);         % Gx spoiler

gzSpoil = mr.makeTrapezoid('z',system,'Area',4/sliceThickness,'Duration',Tpre);

gzReph = mr.makeTrapezoid('z',system,'Area',-gz_acs.area/2,'Duration',Tpre);


areaY = ((0:Ny-1)-Ny/2)*deltak;


temp = 1:Ny;
iY_acs_indices = temp(1+end/2-Ny_acs/2:end/2+Ny_acs/2);

% Make trapezoids for inner loop to save computation
gyPre_acs = [];
gyReph_acs = [];
for iY = iY_acs_indices
    t1 = mr.makeTrapezoid('y','Area',areaY(iY),'Duration',Tpre);
    t2 = mr.makeTrapezoid('y','Area',-areaY(iY),'Duration',Tpre);

    gyPre_acs{iY} = t1;
    gyReph_acs{iY} = t2;    
end


delayTE_acs = ceil( (TE_acs - mr.calcDuration(rf_acs) + mr.calcRfCenter(rf_acs) + rf_acs.delay - mr.calcDuration(gxPre)  ...
    - mr.calcDuration(gro_acs)/2)/seq.gradRasterTime)*seq.gradRasterTime;

assert(delayTE_acs > 0)



delayTR_acs = ceil((TR_acs - mr.calcDuration(rf_acs) - mr.calcDuration(gxPre) ...
    - mr.calcDuration(gro_acs) - mr.calcDuration(gxSpoil) - delayTE_acs)/seq.gradRasterTime)*seq.gradRasterTime;

assert(delayTR_acs > 0)



%--------------------------------------------------------------------------
%% HASTE + ACS: Define sequence blocks
%--------------------------------------------------------------------------

% Next, the blocks are put together to form the sequence
for kex=1:nex % MZ: we start at 0 to have one dummy
    for s=1:Nslices
        rfex.freqOffset=GSex.amplitude*sliceThickness*(s-1-(Nslices-1)/2);
        rfref.freqOffset=GSref.amplitude*sliceThickness*(s-1-(Nslices-1)/2);
        rfex.phaseOffset=rfex_phase-2*pi*rfex.freqOffset*mr.calcRfCenter(rfex); % align the phase for off-center slices
        rfref.phaseOffset=rfref_phase-2*pi*rfref.freqOffset*mr.calcRfCenter(rfref); % dito
    
        seq.addBlock(GS1);
        seq.addBlock(GS2,rfex);
        seq.addBlock(GS3,GR3);
        %GS4.first=GS4f;
        %GS4.first=GS3.last;
        for kech=1:necho
            if (kex>0)
                disp(['PE line: ', num2str(PEorder(kech))])

                phaseArea=phaseAreas(kech,kex);
            else
                phaseArea=0;
            end

            GPpre = mr.makeTrapezoid('y',system,'Area',phaseArea,'Duration',tSp,'riseTime',dG);
            GPrew = mr.makeTrapezoid('y',system,'Area',-phaseArea,'Duration',tSp,'riseTime',dG);
            
            seq.addBlock(GS4,rfref);
            seq.addBlock(GS5,GR5,GPpre);

            if (kex>0)
                if use_wave_y == 1
                    seq.addBlock(GR6,gy_wave,adc);
                else
                    seq.addBlock(GR6,adc);
                end
            else
                seq.addBlock(GR6);
            end
            
            seq.addBlock(GS7,GR7,GPrew);
            %GS4.first=GS7.last;
        end
        seq.addBlock(GS4);
        seq.addBlock(GS5);
        seq.addBlock(delayTR);

        % add delay
        seq.addBlock(delayEnd);

        % acs scan after the delay
        rf_acs.freqOffset=GSex.amplitude*sliceThickness*(s-1-(Nslices-1)/2);

        for iY = 1:Ndummy
            % dummies for acs
            % RF
            rf_acs.phaseOffset = mod(117*(iY^2+iY+2)*pi/180,2*pi);  % add RF phase spoiling with 117 degrees
            seq.addBlock(rf_acs, gz_acs);
        
            % Gradients    
            seq.addBlock(gxPre,gyPre_acs{iY_acs_indices(1)},gzReph);                      % add Gx pre-winder, go to desired ky
            seq.addBlock(mr.makeDelay(delayTE_acs));            % add delay needed before the start of readout
            
            seq.addBlock(gro_acs);                                   % add readout Gx
            
            seq.addBlock(gyReph_acs{iY_acs_indices(1)},gxSpoil,gzSpoil);                   % add Gx spoiler, and go back to DC in ky
            seq.addBlock(mr.makeDelay(delayTR_acs))             % add delay to the end of TR
        end

        % now actual lines for acs
        for iY = iY_acs_indices
            % RF spoiling
            rf_acs.phaseOffset = mod(117*(iY^2+iY+2)*pi/180,2*pi);
            adc_acs.phaseOffset = rf_acs.phaseOffset;

            % Excitation
            seq.addBlock(rf_acs, gz_acs);

            % Encoding
            seq.addBlock(gxPre,gyPre_acs{iY},gzReph);        % Gz, Gy blips, Gx pre-winder
            seq.addBlock(mr.makeDelay(delayTE_acs));    % delay until readout

            seq.addBlock(gro_acs,adc_acs);                       % readout
        
            seq.addBlock(gyReph_acs{iY},gxSpoil,gzSpoil);% -Gz, -Gy blips, Gx spoiler
            seq.addBlock(mr.makeDelay(delayTR_acs))     % wait until end of TR
        end
    end
end

% extra delay at the end:
% seq.addBlock(delayEnd);

   

%--------------------------------------------------------------------------
%% check whether the timing of the sequence is correct
%--------------------------------------------------------------------------

[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end



%--------------------------------------------------------------------------
%% Plot and write sequence
%--------------------------------------------------------------------------

% Set definitions for recon
seq.setDefinition('FOV', [fov fov Nslices*sliceThickness]);
seq.setDefinition('Matrix', [Nx Ny Nslices]);
seq.setDefinition('necho', necho);
seq.setDefinition('Nslices', Nslices);
seq.setDefinition('os_factor', os_factor);
seq.setDefinition('num_cycles', num_cycles);
seq.setDefinition('Name', 'wave_tse');

write_seq = 1;
plot_ktraj = 1;
plot_seq = 1;

save_dir = [pwd, '/haste/'];

if ~isfolder(save_dir)
    mkdir(save_dir)
end

if write_seq
    seq.write([save_dir, 'wave_haste_acs_wY_', num2str(use_wave_y),'_Ry_', num2str(Ry),'.seq']) ;      % Write to pulseq file
end


if plot_seq
    seq.plot('timeRange', [0 10], 'TimeDisp', 'ms');
end


if plot_ktraj
    [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

    % plot k-spaces
    figure; plot(t_ktraj,ktraj'); title('k-space components as functions of time'); % plot the entire k-space trajectory
    figure; plot(ktraj(1,:),ktraj(2,:),'b',...
        ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
    axis('equal'); % enforce aspect ratio for the correct trajectory display
    title('2D k-space');
end


%--------------------------------------------------------------------------
%% load data
%--------------------------------------------------------------------------

data_path = '/autofs/cluster/berkin/berkin/Matlab_Code_New/pulseq/haste/2025_06_19_v1/';

tic
    dt1 = mapVBVD2([data_path, 'meas_MID00086_FID140246_pulseq_mgh_haste_wave0_R3_opt.dat']);
    dt2 = mapVBVD2([data_path, 'meas_MID00087_FID140247_pulseq_mgh_haste_wave1_R3_opt.dat']);
    
    k1 = dt1{end}.image.unsorted();
    k2 = dt2{end}.image.unsorted();
toc


%--------------------------------------------------------------------------
%% parse data
%--------------------------------------------------------------------------

Nx = 120; 
Ny = 120; 

% Ry = 2;
% Ny_pre = 12;               % lines of k-space prior to TE
Ry = 3;
Ny_pre = 36;               % lines of k-space prior to TE

PEorder = ((-Ny_pre):Ry:Ny/2-1)';

PEorder2use = PEorder + Ny/2+1;

PF_factor = 1 - length(find(PEorder<0)) / length(PEorder)

disp(['PE order = ', num2str(PEorder.')])

os_factor = 5;
Nslices = 1;

% necho = Ny/2+Ny_pre;      % ETL
necho = length(PEorder);

disp(['ETL = ', num2str(necho)])

Ny_acs = 32;

num_chan = size(k1,2);

% extract acs (last 24 lines)
kspace_acs1 = zeross([Nx*os_factor,num_chan,Ny]);
kspace_acs2 = zeross([Nx*os_factor,num_chan,Ny]);

kspace_acs1(:,:,1+end/2-Ny_acs/2:end/2+Ny_acs/2) = k1(:,:,end-Ny_acs+1:end);
kspace_acs2(:,:,1+end/2-Ny_acs/2:end/2+Ny_acs/2) = k2(:,:,end-Ny_acs+1:end);

kspace_acs1 = permute(kspace_acs1, [1,3,2]);
kspace_acs2 = permute(kspace_acs2, [1,3,2]);


% remove acs lines
k1_remAcs = k1(:,:,1:end-Ny_acs);
k2_remAcs = k2(:,:,1:end-Ny_acs);

kspace1 = zeross([Nx*os_factor,num_chan,Ny]);
kspace2 = zeross([Nx*os_factor,num_chan,Ny]);

kspace1(:,:,PEorder2use) = k1_remAcs;
kspace2(:,:,PEorder2use) = k2_remAcs;

kspace1 = permute(kspace1, [1,3,2]);
kspace2 = permute(kspace2, [1,3,2]);


mosaic(rsos(kspace1,3),1,1,1,'',[0,1e-3]), setGcf(.5)
mosaic(rsos(kspace2,3),1,1,2,'',[0,1e-3]), setGcf(.5)

mosaic(rsos(kspace_acs1,3),1,1,11,'',[0,1e-3])
mosaic(rsos(kspace_acs2,3),1,1,12,'',[0,1e-3])


%--------------------------------------------------------------------------
%% ifft
%--------------------------------------------------------------------------

img1 = ifft2call(kspace1);
img2 = ifft2call(kspace2);

img_acs1 = ifft2call(kspace_acs1);
img_acs2 = ifft2call(kspace_acs2);

mosaic(rsos(img1,3),1,1,11,'',[0,2e-3],-90), setGcf(.5)
mosaic(rsos(img2,3),1,1,12,'',[0,2e-3],-90), setGcf(.5)

mosaic(rsos(img_acs1,3),1,1,21,'acs1',[0,2e-3],-90), setGcf(.5)
mosaic(rsos(img_acs2,3),1,1,22,'acs2',[0,2e-3],-90), setGcf(.5)


%--------------------------------------------------------------------------
%% estimate psf-y 
%--------------------------------------------------------------------------

psf_y_est = exp(1i * angle( mean(mean(fftc(img2, 1) .* conj(fftc(img1, 1)), 3), 4) ));

mosaic(angle(psf_y_est),1,1,13,'',[-pi,pi], 90), colormap jet


img_deconv = ifftc(fftc(img2,1) .* repmat(conj(psf_y_est),[1,1,num_chan]),1);

mosaic(rsos(img_deconv,3),1,1,14,'',[0,2e-3],-90)


%--------------------------------------------------------------------------
%% load psf-y from separate TSE scans
%--------------------------------------------------------------------------

load([data_path, 'psf_y_tse.mat']);
psf_y_tse = psf_y;

img_deconv_tse = ifftc(fftc(img2,1) .* repmat(conj(psf_y_tse),[1,1,num_chan]),1);

mosaic(rsos(img_deconv_tse,3),1,1,15,'',[0,2e-3],-90)


%--------------------------------------------------------------------------
%% espirit on acs data
%--------------------------------------------------------------------------

addpath /autofs/cluster/berkin/berkin/Matlab_Code_New/TOOLBOXES/Espirit_matlab_only_toolbox/utils/


num_acs = Ny_acs;
kernel_size = [6,6];
eigen_thresh = 0.6;

img_acs = img_acs1(1+end/2-Nx/2:end/2+Nx/2,:,:);

receive = zeross(size(img_acs));

    
tic
for slc_select = 1:s(img_acs,4)     
    disp(num2str(slc_select))
    
    [maps, weights] = ecalib_soft( fft2c( sq(img_acs(:,:,:,slc_select)) ), num_acs, kernel_size, eigen_thresh );

    receive(:,:,:,slc_select) = dot_mult(maps, weights >= eigen_thresh );
end 
toc

mosaic(receive,4,8,1,'',[0,.5],-90)
mosaic(angle(receive),4,8,2,'',[-pi,pi],-90)

save([data_path, 'receive_R3.mat'], 'receive')


%--------------------------------------------------------------------------
%% sense recon
%--------------------------------------------------------------------------

addpath /autofs/cluster/berkin/berkin/Matlab_Code/TOOLBOXES/SENSE_LSQR_Toolbox

load([data_path,'receive_R3.mat'])

img1_crop = img1(1+end/2-Nx/2:end/2+Nx/2,:,:);

kspace1_crop = fft2call(img1_crop);

kspace_accl = zeros(size(kspace1_crop));
kspace_accl(:,1:Ry:end,:) = kspace1_crop(:,1:Ry:end,:);

m2d = zeros(size(kspace1_crop));
m2d(:,1:Ry:end,:) = 1;
 

% LSQR recon
lsqr_iter = 100;
lsqr_tol = 1e-3;

param = [];
param.m2d = m2d;
param.sens = receive;
param.N = [Nx,Ny];
param.num_chan = num_chan;
param.lambda = 1e-3;        % L2 reg

tic
    res = lsqr(@apply_sense_tikc, cat(1, kspace_accl(:), zeros(prod(param.N),1)), lsqr_tol, lsqr_iter, [], [], [], param);  
toc

Res = reshape(res, param.N);


mosaic(coil_combine(img_acs, receive, 3), 1, 1, 12, 'acs', [0,1e-3],-90), setGcf(.5)
mosaic(Res, 1, 1, 13, 'sense recon', [0,4e-3],-90), setGcf(.5)
mosaic(rsos(kspace_accl,3),1,1,14,'kspace accl',[0,3e-3],-90), setGcf(.5)
mosaic(fft2call(Res), 1, 1, 15, ' sense recon', [0,4e-3],-90), setGcf(.5)


save([data_path,'meas_MID00086_FID140246_pulseq_mgh_haste_wave0_R3_opt.mat'], 'Res')


%--------------------------------------------------------------------------
%% sense recon with psf-y
%--------------------------------------------------------------------------

kspace_accl = zeros(size(kspace2));
kspace_accl(:,1:Ry:end,:) = kspace2(:,1:Ry:end,:);

m2d = zeros(size(kspace_accl));
m2d(:,1:Ry:end,:) = 1;
 

% LSQR recon
lsqr_iter = 100;
lsqr_tol = 1e-3;

param = [];
param.m2d = m2d;
param.sens = receive;
param.N = [Nx,Ny];

param.psf_length = size(kspace2,1);

% use psf estimated from TSE scan
param.psfY = repmat(psf_y_tse, [1,1,num_chan]);
param.psfYconj = conj(param.psfY);

param.num_chan = num_chan;
param.lambda = 1e-3;        % L2 reg

tic
    res = lsqr(@apply_sense_tikc_waveY, cat(1, kspace_accl(:), zeros(prod(param.N),1)), lsqr_tol, lsqr_iter, [], [], [], param);  
toc

ResY = reshape(res, param.N);

mosaic(ResY, 1, 1, 23, 'wave sense recon', [0,4e-3],-90), setGcf(.5)
mosaic(rsos(kspace_accl,3),1,1,24,'kspace accl',[0,3e-3],-90), setGcf(.5)
mosaic(fft2call(ResY), 1, 1, 25, 'wave sense recon', [0,4e-3],-90), setGcf(.5)

save([data_path,'meas_MID00087_FID140247_pulseq_mgh_haste_wave1_R3_opt_ResY.mat'], 'ResY')



