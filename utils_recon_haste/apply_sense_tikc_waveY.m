function [ res, tflag ] = apply_sense_tikc_waveY( in, params, tflag )
%APPLY_COLLAPSE Summary of this function goes here
%   Detailed explanation goes here

if strcmp(tflag,'transp')
    
    % Transposed SENSE operator:
    % IFFT coil k-space, multiply by conjugate of coil sensitivities, then
    % sum across channels

    b = in(1+params.num_chan*params.N(2)*params.psf_length:end);    

    kspace_coils = reshape(in(1:params.num_chan*params.N(2)*params.psf_length), [params.psf_length, params.N(2), params.num_chan]);    

    hybrid_coils = ifftc( kspace_coils .* params.m2d, 2 ) .* params.psfYconj;
    
    img_coils = ifftc(hybrid_coils,1);

    Res = sum(img_coils(1+end/2-params.N(1)/2:end/2+params.N(1)/2,:,:) .* conj(params.sens), 3);
    
    res = Res(:) + sqrt(params.lambda) * b;
    
else
    
    % Forward SENSE operator:
    % multiply by coil sensitivities, take undersampled FFT
    
    img_coils = repmat(reshape(in, params.N), [1,1,params.num_chan]) .* params.sens;

    hybrid_coils_pad = fftc(padarray(img_coils, params.psf_length/2 - params.N(1)/2), 1);

    % kspace_coils = fft2call(img_coils) .* params.m2d;
    kspace_coils = fftc(hybrid_coils_pad .* params.psfY, 2) .* params.m2d;

    res = cat(1, kspace_coils(:), sqrt(params.lambda)*in);
    
end

end
