% Lecture: Photogrammetric Computer Vision
% Exercise 3: Background subtraction
% Author: M.Sc. Jakob Unger (unger@ipi.uni-hannover.de)
% Lecturer: Dipl.-Ing. Tobias Klinger (klinger@ipi.uni-hannover.de)
% Lecturer: M.Sc. Lin Chen (chen@ipi.uni-hannover.de)
% Group: <group number>
% Authors: <Chawda Vimal>

close all
clear all

%% Read image sequence
path2sequence = 'sequence';
search_string = fullfile(path2sequence, '*.jpeg');
file_list = dir(search_string);

%% Sequencial estimation of the gaussians parameters

% TODO: Learning rate
% alpha = ?;

% Initial values
im_RGB = imread(fullfile(path2sequence, file_list(1).name));
[m,n] = size(rgb2gray(im_RGB));
n1=50;
n2 = 1400;
alpha = 1/n2;
mu = single(rgb2gray(im_RGB)); 
sigma_square = ones(m,n)*100;

% Structure element for morphological operation
se = strel('square',3);

h1 = figure(1);
colormap(gray);

% Iterate over the sequence
for i = 2:length(file_list)
    colormap(gray);
    im_RGB = imread(fullfile(path2sequence, file_list(i).name));
    im = single(rgb2gray(im_RGB));
    
    % TODO: Thresholding for background subtraction
    %delta_g = ?
    %return_mask = ? -> background pixels are 1, foreground pixels are set to 0
    %return_mask = logical(return_mask); %convert to type "logical"
     del_g =abs(im-mu);                                                                                                                                               variable =2.5*sqrt(sigma_square);                                                                                                                              
     B=del_g<variable;                                                                                                                                        
     return_mask = logical(B);                                                                                                                                      
     before_mask =return_mask; 
     
    % Eliminate noise
    % TODO: Closing using structure element "se"
    return_mask = imdilate(return_mask,se);    
    return_mask = imerode(return_mask,se); 
    
    % Update Gaussian parameters 
    % TODO: Mean
    % TODO: Variance
    mu = (1-alpha)*mu+alpha*im;   
    sigma_square = (1-alpha)*sigma_square + alpha*((mu-im).^2); 
    % Output
    subplot(231); 
    imagesc(im_RGB, [0 255]);title('Input Image') % TODO: Show input image   we can write or im or im_RGB
    subplot(232); 
    imagesc(mu, [0 255]); title('Mean')% TODO: Show mean values
    subplot(233);  
    imagesc(sigma_square, [0 255]);title('Variance') % TODO: Show variances
    subplot(234);  
    imagesc(del_g, [0 255]);title('Difference DeltaG') % TODO: Show difference between current image and mean values
    subplot(235);  
    imagesc(before_mask, [0 1]);title('Unfiltered Background') % TODO: Show background mask (unfiltered)
    subplot(236);  
    imagesc(return_mask, [0 1]); title('Filtered Background')% TODO: Show background mask after eliminating noise
    drawnow; 
end