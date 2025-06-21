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
 

%--------------------------------------------------------------------------
%% load data: R=3 accelerated haste w/ and w/o wave encoding
%--------------------------------------------------------------------------

addpath 

tic
    dt1 = mapVBVD2('meas_MID00086_FID140246_pulseq_mgh_haste_wave0_R3_opt.dat');
    dt2 = mapVBVD2('meas_MID00087_FID140247_pulseq_mgh_haste_wave1_R3_opt.dat');
    
    k1 = dt1{end}.image.unsorted();
    k2 = dt2{end}.image.unsorted();
toc


%--------------------------------------------------------------------------
%% parse data
%--------------------------------------------------------------------------

Nx = 120; 
Ny = 120; 

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
%% load psf-y from separate TSE scans -> run script_wave_tse_recon.m
%--------------------------------------------------------------------------

load('psf_y_tse.mat');
psf_y_tse = psf_y;

img_deconv_tse = ifftc(fftc(img2,1) .* repmat(conj(psf_y_tse),[1,1,num_chan]),1);

mosaic(rsos(img_deconv_tse,3),1,1,15,'',[0,2e-3],-90)


%--------------------------------------------------------------------------
%% espirit on acs data
%--------------------------------------------------------------------------

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


%--------------------------------------------------------------------------
%% sense recon: for Haste without wave encoding
%--------------------------------------------------------------------------


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



%--------------------------------------------------------------------------
%% sense recon with psf-y: for Haste with wave encoding enabled
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




