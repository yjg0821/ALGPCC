close all;
clear;

path_src = 'images';
Dir = dir(path_src);
Nimgs = length(Dir) - 2;

Time_ALGPCC = zeros(Nimgs, 1);
RMSE_ALGPCC = zeros(Nimgs, 9);
PSNR_ALGPCC = zeros(Nimgs, 9);
SSIM_ALGPCC = zeros(Nimgs, 9);
APMR_ALGPCC = zeros(Nimgs, 1);

for i =  1 : Nimgs
    RGB_0 = double(imread(fullfile(path_src, Dir(i + 2).name, '0.png')));
    RGB_45 = double(imread(fullfile(path_src, Dir(i + 2).name, '45.png')));
    RGB_90 = double(imread(fullfile(path_src, Dir(i + 2).name, '90.png')));
    RGB_135 = double(imread(fullfile(path_src, Dir(i + 2).name, '135.png')));

    I_truth(:,:,1) = RGB_0(:,:,2);
    I_truth(:,:,2) = RGB_45(:,:,2);
    I_truth(:,:,3) = RGB_90(:,:,2);
    I_truth(:,:,4) = RGB_135(:,:,2);
    
    Stokes_truth = StokesCalculation(I_truth);
    
    subplot(2,3,1), imshow(Stokes_truth(10:end-10,10:end-10,5),[])
    subplot(2,3,2), imshow(Stokes_truth(10:end-10,10:end-10,6),[])
    subplot(2,3,3), imshow(Stokes_truth(10:end-10,10:end-10,7),[])
    subplot(2,3,4), imshow(Stokes_truth(10:end-10,10:end-10,8),[0,1])
    subplot(2,3,5), imshow(Stokes_truth(10:end-10,10:end-10,9),[0,180])
    
    [I_raw, I_mosaic, I_mask] = mosaic_polarized(I_truth);
    
    tic
    I_ALGPCC = ALGPCC(I_raw, I_mosaic, I_mask);
    Time_ALGPCC(i) = toc;
    
    Stokes_ALGPCC = StokesCalculation(I_ALGPCC);
    [RMSE_ALGPCC(i, :), PSNR_ALGPCC(i, :), SSIM_ALGPCC(i, :)] = RMSE_PSNR_SSIM(Stokes_ALGPCC, Stokes_truth, [65535, 1, 180], 10);
    
    figure,
    subplot(2,3,1), imshow(Stokes_ALGPCC(10:end-10,10:end-10,5),[])
    subplot(2,3,2), imshow(Stokes_ALGPCC(10:end-10,10:end-10,6),[])
    subplot(2,3,3), imshow(Stokes_ALGPCC(10:end-10,10:end-10,7),[])
    subplot(2,3,4), imshow(Stokes_ALGPCC(10:end-10,10:end-10,8),[0,1])
    subplot(2,3,5), imshow(Stokes_ALGPCC(10:end-10,10:end-10,9),[0,180])
end


