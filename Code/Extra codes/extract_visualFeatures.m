function [R, G, B, L, C, H, S, V ] = extract_visualFeatures(inputpatient)

% [R, G, B, L, C, H, S, V ] = extract_visualFeatures(inputpatient)
% Input is a string with the PId of the patient
% Output are the visual features of each used stimulus.




%% This code extracts visual features of the images; only the selected patient

folder= '/media/raid/Alex/Stimuli/Used Stimuli/';

folder=  [ folder , inputpatient,'sat1/'];
patdir=  [ folder];
folder=  [ folder, inputpatient,'sat1/'];

pp=1;


listdir = dir([folder, '/*.jpg']);
if numel(listdir) ~= 20
    error('not 20 stimuli in folder')
    
end

% check if strictly increasing the name
clear ff
for ll=1:numel(listdir)
    ff(ll,:) = listdir(ll).name;
    folder= '/media/raid/Alex/Stimuli/Used Stimuli/';
    
    folder=  [ folder , inputpatient,'sat1/'];
    patdir=  [ folder];
    folder=  [ folder, inputpatient,'sat1/'];
    
    pp=1;
    
    
    listdir = dir([folder, '/*.jpg']);
    if numel(listdir) ~= 20
        error('not 20 stimuli in folder')
        
    end
    
    % check if strictly increasing the name
    clear ff
    for ll=1:numel(listdir)
        ff(ll,:) = listdir(ll).name;
    end
    
    
    if all(diff(str2num(ff(:,1:3)))>0) % check if strictly increasing the name; otherwise error
        
        for stim = 1:numel(listdir)
            image = [folder, listdir(stim).name];
            
            im1= imread(image);
            
            pR = (im1(:,:,1)); % p means pixel-wise
            pG = (im1(:,:,2));
            pB = (im1(:,:,3));
            
            pL = 0.2126 * pR + 0.7152 * pG  +  0.0722 * pB;
            grayImage = rgb2gray(im1); % pixel intensity
            
            R(stim, pp) = mean2(pR);
            G(stim, pp) = mean2(pG);
            B(stim, pp) = mean2(pB);
            L(stim, pp) = mean2(pL);
            C(stim, pp) = std2(pL); % local contrast; Suzuki et al define it as the std of the luminance, others define it as the std of intensity (in gray color).
            I(stim, pp) = mean2(grayImage); % pixel intensity
            C2(stim, pp) = std2(grayImage); % local contrast
        
            HSV = rgb2hsv(im1);
            
            H(stim, pp) = mean2(HSV(:,:,1));
            S(stim, pp) = mean2(HSV(:,:,2));
            V(stim, pp) = mean2(HSV(:,:,3));
        end
    else
        disp([inputpatient, '  not strictly increasing stimulus ID'])
    end
    
end


if all(diff(str2num(ff(:,1:3)))>0) % check if strictly increasing the name; otherwise error
    
    for stim = 1:numel(listdir)
        image = [folder, listdir(stim).name];
        
        im1= imread(image);
        
        pR = (im1(:,:,1)); % p means pixel-wise
        pG = (im1(:,:,2));
        pB = (im1(:,:,3));
        
        pL = 0.2126 * pR + 0.7152 * pG  +  0.0722 * pB;
        grayImage = rgb2gray(im1); % pixel intensity
        
        R(stim, pp) = mean2(pR);
        G(stim, pp) = mean2(pG);
        B(stim, pp) = mean2(pB);
        L(stim, pp) = mean2(pL);
        C(stim, pp) = std2(pL); % pixel intensity
        I(stim, pp) = mean2(grayImage); % pixel intensity
        C2(stim, pp) = std2(grayImage); % local contrast
        
        HSV = rgb2hsv(im1);
        
        H(stim, pp) = mean2(HSV(:,:,1));
        S(stim, pp) = mean2(HSV(:,:,2));
        V(stim, pp) = mean2(HSV(:,:,3));
    end
else
    disp([inputpatient, '  not strictly increasing stimulus ID'])
end

% end



