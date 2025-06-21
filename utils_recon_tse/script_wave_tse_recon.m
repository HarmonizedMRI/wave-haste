%--------------------------------------------------------------------------
%% load data acquired with and without wave encoding using the TSE sequence
%--------------------------------------------------------------------------

tic
    dt1 = mapVBVD2('meas_MID00085_FID140245_pulseq_mgh_tse_wave0_opt.dat');
    dt2 = mapVBVD2('meas_MID00083_FID140243_pulseq_mgh_tse_wave1_opt.dat');
    
    k1 = dt1{end}.image.unsorted();
    k2 = dt2{end}.image.unsorted();
toc


%--------------------------------------------------------------------------
%% parse data
%--------------------------------------------------------------------------

% load phase encoding order per each TR
load('PEorder_tse_wY_0.mat')

Nx = 120; 
Ny = 120; 

os_factor = 5;
Nslices = 1;

num_chan = size(k1,2);

kspace1 = zeross([Nx*os_factor,num_chan,Ny]);
kspace2 = zeross([Nx*os_factor,num_chan,Ny]);

PEorder2use = PEorder + Ny/2+1;

necho = size(PEorder,1);

disp(['ETL ', num2str(necho)])
disp(['num TR ', num2str(size(PEorder,2))])

for t = 1:size(PEorder,2)
    temp1 = k1(:,:,1+(t-1)*necho:t*necho,:);
    temp2 = k2(:,:,1+(t-1)*necho:t*necho,:);

    kspace1(:,:,PEorder2use(:,t)) = temp1;
    kspace2(:,:,PEorder2use(:,t)) = temp2;
end


kspace1 = permute(kspace1, [1,3,2]);
kspace2 = permute(kspace2, [1,3,2]);

mosaic(rsos(kspace1,3),1,1,1,'',[0,1e-3]),setGcf(.5)
mosaic(rsos(kspace2,3),1,1,2,'',[0,1e-3]),setGcf(.5)


%--------------------------------------------------------------------------
%% ifft
%--------------------------------------------------------------------------

img1 = ifft2call(kspace1);
img2 = ifft2call(kspace2);

mosaic(rsos(img1,3),1,1,11,'',[0,6e-3],-90), setGcf(.5)
mosaic(rsos(img2,3),1,1,12,'',[0,6e-3],-90), setGcf(.5) 


%--------------------------------------------------------------------------
%% psf estimation
%--------------------------------------------------------------------------

psf_y = exp(1i * angle( mean(mean(fftc(img2, 1) .* conj(fftc(img1, 1)), 3), 4) ));

mosaic(angle(psf_y),1,1,12,'',[-pi,pi], -90), colormap jet


img_deconv = ifftc(fftc(img2,1) .* repmat(conj(psf_y),[1,1,num_chan]),1);

mosaic(rsos(img_deconv,3),1,1,13,'',[0,6e-3],-90)

% save estimated psf to be used in haste recon
save('psf_y_tse.mat', 'psf_y')

 