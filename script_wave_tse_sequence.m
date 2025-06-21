%--------------------------------------------------------------------------
%% create wave-tse pulse sequence to esimate point spread function
%% which will be later used in wave-haste recon 
%--------------------------------------------------------------------------

addpath utils_seq_tse/

clear;clc; 

plot_seq = 1;
write_seq = 1;
plot_ktraj = 1;

save_PEorder = 1;           % set to 1 to save phase encoding order: needed for k-space recon    
use_wave_y = 1;             % set to 1 to enable wave encoding gradients

sys_type = 'prisma';

gMax = 45;
sMax = 140;                 % derate slew rate to avoid PNS

sys = mr.opts( ...
    'MaxGrad', gMax, 'GradUnit', 'mT/m', ...
    'MaxSlew', sMax, 'SlewUnit', 'T/m/s', ...
    'rfRingdownTime', 100e-6, ...
    'rfDeadTime', 100e-6, ...
    'adcDeadTime', 10e-6);


seq = mr.Sequence(sys);

fov = 240e-3;
Nx = 120;
Ny = 120;                   % 2 mm resolution

necho = 2;                  % echo train length (ETL)
Nslices = 1;
sliceThickness = 5e-3;

TE1 = 6.6e-3;               % echo time of the first echo in the train -> also controls the echo spacing (ESP)
TR = 2000e-3;

% TEeff = 100e-3; % the desired echo time (can only be achieved approximately)
TEeff = TE1;                % set to TE1 since contrast is not important for wave psf estimation


deltak = 1/fov;
kWidth = Nx*deltak;
dG=250e-6;                  % 'standard' ramp time - makes sequence structure much simpler


rflip=180;

if (numel(rflip)==1), rflip=rflip+zeros([1 necho]); end


%-------------------------------------------------------------------------
%% Define wave params
%-------------------------------------------------------------------------
        
os_factor = 5;              % readout oversampling factor, needed to avoid circular convolution artifacts
dwell           =   5e-6;
Tread           =   dwell*Nx * os_factor;


num_read_points = Nx * os_factor;

time_per_adc_point = Tread / num_read_points;

gro = mr.makeTrapezoid('x',sys,'FlatArea',kWidth,'FlatTime',Tread,'riseTime',dG);


gYmax_wave = 5;             % mT/m, can be increased for greater amount of spreading

sYmax_wave = sMax; % T/m/s

num_cycles = 10;            % number of wave cycles

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

wavepoints = round(Tread / sys.gradRasterTime);     % number of readout points
T_wavepoints = sys.gradRasterTime;                  % grad raster time in seconds

TimePerSineY = (Tread / sys.gradRasterTime) * T_wavepoints / waveParams.NcyclesGy;     % time per sine wave in seconds

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

scaling_factor = sys.gamma * 1e-2;                          % convert from G/cm to T/m by multiplying by 1e-2, then to Hz/m by multiplying by gamma

GwaveY = G0_Y * sin(wY * GradTimePointsY) * scaling_factor; % Gy amplitude in Hz/m


disp(['Gx amplitude: ', num2str(max(abs(gro.amplitude)*1e-3)), ' kHz/m'])
disp(['Gy amplitude: ', num2str(max(abs(GwaveY)*1e-3)), ' kHz/m'])


% pad Gy waveform at the beginning and end with zeroes to account for rise and fall time of trapezoid
% SF: This is not needed here because of G5 and G7 
num_pad_pre = round(sys.adcDeadTime / sys.gradRasterTime);

% -1 to account for added point at the end to make sure sine wave goes back
% to zero due to [0:wavepoints] instead of [0:wavepoints-1]
% num_pad_post = round(gro.fallTime / sys.gradRasterTime) - 1;     
num_pad_post = round(sys.adcDeadTime / sys.gradRasterTime) - 1;     


GwaveY = cat(2, zeros(1,num_pad_pre), GwaveY, zeros(1,num_pad_post));

disp(['num points Gy:', num2str(length(GwaveY))])


% form the Gx waveform for display
GwaveX = gro.amplitude * ones(1, round(gro.flatTime / sys.gradRasterTime));

GX_rampup = linspace(0, gro.amplitude, round(gro.riseTime / sys.gradRasterTime));
GX_rampdown = linspace(gro.amplitude, 0, round(gro.fallTime / sys.gradRasterTime));

GwaveX = cat(2, GX_rampup, GwaveX, GX_rampdown);

disp(['num points Gx:', num2str(length(GwaveX))])


% form the Gy object -> play on physical y axis
gy_wave = mr.makeArbitraryGrad('y', GwaveY);


readoutTime = Tread+ 2*sys.adcDeadTime;

disp(num2str(readoutTime*1e3))

tEx=2.5e-3;         % excite duration
tExwd=tEx+sys.rfRingdownTime+sys.rfDeadTime;
    
tRef=2.0e-3;        % refocus duration
tRefwd=tRef+sys.rfRingdownTime+sys.rfDeadTime;

tSp=0.5*(TE1-readoutTime-tRefwd);
tSpex=0.5*(TE1-tExwd-tRefwd);
fspR=1.0;
fspS=0.5;

rfex_phase=pi/2;
rfref_phase=0;


% create excitation pulse
flipex=90*pi/180;
[rfex, gz] = mr.makeSincPulse(flipex,sys,'Duration',tEx,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,'PhaseOffset',rfex_phase);
GSex = mr.makeTrapezoid('z',sys,'amplitude',gz.amplitude,'FlatTime',tExwd,'riseTime',dG);

% and refocusing pulse
flipref=rflip(1)*pi/180;
[rfref, gz2] = mr.makeSincPulse(flipref,sys,'Duration',tRef,... % it was a bug as 'gz' was owerwritten
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,'PhaseOffset',rfref_phase,'use','refocusing');
GSref = mr.makeTrapezoid('z',sys,'amplitude',GSex.amplitude,'FlatTime',tRefwd,'riseTime',dG);

AGSex=GSex.area/2;
GSspr = mr.makeTrapezoid('z',sys,'area',AGSex*(1+fspS),'duration',tSp,'riseTime',dG);
GSspex = mr.makeTrapezoid('z',sys,'area',AGSex*fspS,'duration',tSpex,'riseTime',dG);


% adc, readout and excitation gradients

GRacq = gro;
adc = mr.makeAdc(Nx*os_factor,'Duration',gro.flatTime, 'Delay',sys.adcDeadTime);     

GRspr = mr.makeTrapezoid('x',sys,'area',GRacq.area*fspR,'duration',tSp,'riseTime',dG);
GRspex = mr.makeTrapezoid('x',sys,'area',GRacq.area*(1+fspR),'duration',tSpex,'riseTime',dG);


AGRspr=GRspr.area;%GRacq.area/2*fspR;
AGRpreph = GRacq.area/2+AGRspr;%GRacq.area*(1+fspR)/2;
GRpreph = mr.makeTrapezoid('x',sys,'Area',AGRpreph,'duration',tSpex,'riseTime',dG);


nex = floor(Ny/necho);
pe_steps=(1:(necho*nex))-0.5*necho*nex-1;
if 0==mod(necho,2)
    pe_steps=circshift(pe_steps,[0,-round(nex/2)]); % for odd number of echoes we have to apply a shift to avoid a contrast jump at k=0
end

% TSE echo time magic??
[~,iPEmin]=min(abs(pe_steps));
k0curr=floor((iPEmin-1)/nex)+1;             % calculate the 'native' central echo index
k0prescr=max(round(TEeff/TE1),1);           % echo to be aligned to the k-space center
PEorder=circshift(reshape(pe_steps,[nex,necho])',k0prescr-k0curr);

phaseAreas = PEorder*deltak;


% rising edge of excitation gradient
GS1times=[0 GSex.riseTime];
GS1amp=[0 GSex.amplitude];
GS1 = mr.makeExtendedTrapezoid('z','times',GS1times,'amplitudes',GS1amp);

% flat top of the excitation gradient
GS2times=[0 GSex.flatTime];
GS2amp=[GSex.amplitude GSex.amplitude];
GS2 = mr.makeExtendedTrapezoid('z','times',GS2times,'amplitudes',GS2amp);

% rephaser + left spoiler combined (this is only after the excitation)
GS3times=[0 GSspex.riseTime GSspex.riseTime+GSspex.flatTime GSspex.riseTime+GSspex.flatTime+GSspex.fallTime];
GS3amp=[GSex.amplitude GSspex.amplitude GSspex.amplitude GSref.amplitude];
GS3 = mr.makeExtendedTrapezoid('z','times',GS3times,'amplitudes',GS3amp);

% slice select gradient during refocusing pulse
GS4times=[0 GSref.flatTime];
GS4amp=[GSref.amplitude GSref.amplitude];
GS4 = mr.makeExtendedTrapezoid('z','times',GS4times,'amplitudes',GS4amp);

% right lobe of spoiler (this is on the right side of the refocusing pulse)
GS5times=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
GS5amp=[GSref.amplitude GSspr.amplitude GSspr.amplitude 0];
GS5 = mr.makeExtendedTrapezoid('z','times',GS5times,'amplitudes',GS5amp);

% and here goes the ADC, gradient GR6

% left lobe of spoiler (this is on the left side of the next refocusing pulse)
GS7times=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
GS7amp=[0 GSspr.amplitude GSspr.amplitude GSref.amplitude];
GS7 = mr.makeExtendedTrapezoid('z','times',GS7times,'amplitudes',GS7amp);


% and now the readout gradient....

GR3=GRpreph;%GRspex;

GR5times=[0 GRspr.riseTime GRspr.riseTime+GRspr.flatTime GRspr.riseTime+GRspr.flatTime+GRspr.fallTime];
GR5amp=[0 GRspr.amplitude GRspr.amplitude GRacq.amplitude];
GR5 = mr.makeExtendedTrapezoid('x','times',GR5times,'amplitudes',GR5amp);

% read gradient on x axis
GR6times=[0 readoutTime];
GR6amp=[GRacq.amplitude GRacq.amplitude];
GR6 = mr.makeExtendedTrapezoid('x','times',GR6times,'amplitudes',GR6amp);

GR7times=[0 GRspr.riseTime GRspr.riseTime+GRspr.flatTime GRspr.riseTime+GRspr.flatTime+GRspr.fallTime];
GR7amp=[GRacq.amplitude GRspr.amplitude GRspr.amplitude 0];
GR7 = mr.makeExtendedTrapezoid('x','times',GR7times,'amplitudes',GR7amp);


% and filltimes
tex=mr.calcDuration(GS1)+mr.calcDuration(GS2)+mr.calcDuration(GS3);
tref=mr.calcDuration(GS4)+mr.calcDuration(GS5)+mr.calcDuration(GS7)+readoutTime;
tend=mr.calcDuration(GS4)+mr.calcDuration(GS5);


tETrain=tex+necho*tref+tend;
TRfill=(TR-Nslices*tETrain)/Nslices;
% round to gradient raster
TRfill=sys.gradRasterTime * round(TRfill / sys.gradRasterTime);
if TRfill<0, TRfill=1e-3;
    disp(strcat('Warning!!! TR too short, adapted to include all slices to : ',num2str(1000*Nslices*(tETrain+TRfill)),' ms'));
else
    disp(strcat('TRfill : ',num2str(1000*TRfill),' ms'));
end
delayTR = mr.makeDelay(TRfill);


% to do wave-Z, we need to modify GS5 and GS7
% for wave-Y only, we can start from putting sine wave on Gy so no need to 
% worry about rising edge of cosine



%
for kex=0:nex % MZ: we start at 0 to have one dummy
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
        for kech=1:necho,
            if (kex>0)
                phaseArea=phaseAreas(kech,kex);
            else
                phaseArea=0;
            end

            % these are the Gy phase encoding blips:
            GPpre = mr.makeTrapezoid('y',sys,'Area',phaseArea,'Duration',tSp,'riseTime',dG);
            GPrew = mr.makeTrapezoid('y',sys,'Area',-phaseArea,'Duration',tSp,'riseTime',dG);

            seq.addBlock(GS4,rfref);
            seq.addBlock(GS5,GR5,GPpre);
            
            if (kex>0)
                % here we insert Gy sine gradient during the ADC
                if (use_wave_y==1)
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
    end
end


% check whether the timing of the sequence is correct
[ok, error_report] = seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end
   


%--------------------------------------------------------------------------
% Plot and write sequence
%--------------------------------------------------------------------------

% Set definitions for recon
seq.setDefinition('FOV', [fov fov Nslices*sliceThickness]);
seq.setDefinition('Matrix', [Nx Ny Nslices]);
seq.setDefinition('necho', necho);
seq.setDefinition('Nslices', Nslices);
seq.setDefinition('os_factor', os_factor);
seq.setDefinition('num_cycles', num_cycles);
seq.setDefinition('Name', 'wave_tse');

save_dir = [pwd, '/tse/'];

if ~isfolder(save_dir)
    mkdir(save_dir)
end

if write_seq
    seq.write([save_dir, 'wave_tse_wY_', num2str(use_wave_y),'.seq']) ;      % Write to pulseq file
end


if plot_seq
    seq.plot('timeRange', [0 TR]*1.1, 'TimeDisp', 'ms');
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

if save_PEorder
    save([save_dir, 'PEorder_tse_wY_', num2str(use_wave_y)], 'PEorder')
end


