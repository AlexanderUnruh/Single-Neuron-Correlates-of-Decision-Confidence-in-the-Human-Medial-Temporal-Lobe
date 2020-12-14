function generateplot1 

%% Figure 3
% This Script plots the Figure 3 of the paper, including a raster plot of
% a single cell and a distribution of the neuron when trials split by the
% presentation mode of the stimuli



dbstop if error

set(groot,'defaultAxesLinewidth', 1)
windowsize = 200; resolution = 1; stepsize=10;  stdkernel= 5;
newwinwosreg = [-5000,5000];  
rangePlot = [-5000:stepsize:7400-1];

possible_colors = [[0 0.4470 0.7410];...
    [0.9290 0.5940 0.1250];... 
    [0.9290 0.3250 0.0980];...
    [0.4940 0.1840 0.5560]];

x_ax_lim = {[-2350, 3601+2350]}; % plot from 2350 ms before trial onset to 2350 after trial onset


inputRegion2 = {'EC'}; % region
memberEC = [69]; % unit ID

%% folder and subfolders to save results and load data

satiationfolder = '/media/Projects/Alex/Reclustered analysis All/';
folder_to_save =  '/media/Projects/Alex/Reclustered analysis All/PaperPlots'; mkdir(folder_to_save)
cd(satiationfolder);     


for l=1:length(inputRegion2)
    mkdir([folder_to_save, '/' char(inputRegion2{l})])
    folder_to_save2{l}=[folder_to_save, '/' char(inputRegion2{l})];
end



%% start analysis
for region=1:length(inputRegion2)
    
    %load neural data
    inputRegion = inputRegion2{region};
    inputRegion=char(inputRegion);
    behavDir=[];
    load ([satiationfolder,'/Rasts matfiles/',char(inputRegion),'_rasts6.mat'], 'rasts', 'outfiles')
    
    for o=1:size(outfiles,1)
        behavDir{end+1}= outfiles{o,2};
    end
    
    
    currentdirectory = pwd;
       

    lastnumber=1;
    for sess=1:size(outfiles,1)
        
        unitind=ismember(cell2mat(rasts(:,2:6)),cell2mat(outfiles(sess,2:6)),'rows'); % check if the unit is part of the channel, compares the name of the patient, etc and region
        n_ses_units = sum(unitind & cell2mat(rasts(:,9)) == 7 & cell2mat(rasts(:,10)) == 3);
        
        for a=1:n_ses_units
            unitnum{sess}(a,1)= lastnumber;
            lastnumber = lastnumber+1;
        end
    end
    
    
    
    for sess=1:size(outfiles,1) % This loops selects the data for the region by session and subject
        
        cd(currentdirectory)
            
        unitind=ismember(cell2mat(rasts(:,2:6)),cell2mat(outfiles(sess,2:6)),'rows'); % check if the unit is part of the channel, compares the name of the patient, etc and region
        nombre_region= [ ' ', char(outfiles(sess,6))];      
        rast_mrg_sess = rasts(unitind ,:); % Units in rast_mrg recorded in this session and this region        
        ttl_n_cond = cell2mat(rast_mrg_sess(:,9:10)); %icludes when the event (7 or 43 or 44) and the criteria (0-20)
        
        pic_rat1 = find((ttl_n_cond(:,1) == 7) & (ttl_n_cond(:,2) == 3));     %  First block judge pics

        pic_A11 = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 11)));   % Pic A1
        pic_A12 = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 12)));   % Pic A1
        pic_B1  = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 13)));   % Pic B
        pic_B2  = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 14)));   % Pic B
        pic_A21  = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 15)));
        pic_A22  = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 16)));
        

        % extract behavioural data
        [satcue,ranking,events]=load_behavioural_data(behavDir{sess});
        [RT,endResponse] = extract_reaction_times(events);
        [Red, Green, Blue, Luminance, Contrast, Hue, Saturation, ~ ] = extract_visualFeatures(outfiles{1,2});
        [Kalorien, Protein, Kohlenhydrate, Zucker, Fett, Ballaststoffe, Salz] = naehwerte(outfiles{1,2});
        
        cd ..
        
  
        
        % define both experimental sets, before satiation (3) and after
        % satiation (6)
        for jj=[3,6]
            
            
            % Output values (subject's reports from likert scale)
            
            val_rat(:,jj/3)= cell2mat(satcue{jj-1}(:,14));
            val_2AFC(:,jj/3)= cell2mat(satcue{jj}(:,14));
            
            
            % Input values
            
            % Rating paradigm
            
            val_rat_Rat(:,jj/3)=nan(length(satcue{jj-1}(:,14)),1);          
            val_rat_Rank(:,jj/3)=nan(length(satcue{jj-1}(:,14)),1);
            
            
            
            for i=1:20
                
                ind3=find(cell2mat(satcue{jj-1}(:,2))==i);
                valRat=  ranking{jj-1}(i,4); % ranking              
                valRank= ranking{jj}(i,3);

                
                val_rat_Rat(ind3,jj/3)= valRat;
                val_rat_Rank(ind3,jj/3)= valRank;
 
                stimOrderRat(ind3,jj/3) = i;
                flavourRat(ind3, jj/3) = i<11;
                
                % visual features
                Red_Rat (ind3,jj/3) = Red(i,1);
                Green_Rat (ind3,jj/3) = Green(i,1);
                Blue_Rat (ind3,jj/3) = Blue(i,1);
                Luminance_Rat (ind3,jj/3) = Luminance(i,1);
                Contrast_Rat (ind3,jj/3) = Contrast(i,1);
                Hue_Rat (ind3,jj/3) = Hue(i,1);
                Saturation_Rat (ind3,jj/3) = Saturation(i,1);
                % nutritional  features
                
                Kalorien_Rat (ind3,jj/3) = Kalorien(i,1);
                Protein_Rat (ind3,jj/3) = Protein(i,1);
                Kohlenhydrate_Rat (ind3,jj/3) = Kohlenhydrate(i,1);
                Zucker_Rat (ind3,jj/3) = Zucker(i,1);
                Fett_Rat (ind3,jj/3) = Fett(i,1);
                Ballaststoffe_Rat (ind3,jj/3) = Ballaststoffe(i,1);
                Salz_Rat (ind3,jj/3) = Salz(i,1);
                
            end
            
            
            % 2AFC paradigm
            
            val_A_Rat(:,jj/3) = nan(length(satcue{jj}(:,14)),1);             % value of image on A given in previous rating trials (-300 to 300)
            val_B_Rat(:,jj/3) = nan(length(satcue{jj}(:,14)),1);             % value of image on B given in previous rating trials (-300 to 300)
            val_A_Rank(:,jj/3) = nan(length(satcue{jj}(:,14)),1);         % rank of image on A, calculated after all the choice trials (0 to 20)
            val_B_Rank(:,jj/3) = nan(length(satcue{jj}(:,14)),1);         % rank of image on B, calculated after all the choice trials (0 to 20)
            
            val_L_Rat(:,jj/3) = nan(length(satcue{jj}(:,14)),1);             % value of image on A given in previous rating trials (-300 to 300)
            val_R_Rat(:,jj/3) = nan(length(satcue{jj}(:,14)),1);             % value of image on B given in previous rating trials (-300 to 300)
            val_L_Rank(:,jj/3) = nan(length(satcue{jj}(:,14)),1);         % rank of image on A, calculated after all the choice trials (0 to 20)
            val_R_Rank(:,jj/3) = nan(length(satcue{jj}(:,14)),1);         % rank of image on B, calculated after all the choice trials (0 to 20)
            
            val_C_Rat(:,jj/3) = nan(length(satcue{jj}(:,14)),1);             % value of image on A given in previous rating trials (-300 to 300)
            val_U_Rat(:,jj/3) = nan(length(satcue{jj}(:,14)),1);             % value of image on B given in previous rating trials (-300 to 300)
            val_C_Rank(:,jj/3) = nan(length(satcue{jj}(:,14)),1);         % rank of image on A, calculated after all the choice trials (0 to 20)
            val_U_Rank(:,jj/3) = nan(length(satcue{jj}(:,14)),1);         % rank of image on B, calculated after all the choice trials (0 to 20)
            
            
            
            % transform the output of the likert and convert it so that Signed Confidence is defined for A s negative values and B Stimulus for positve Values 
            
            lr_to_ab3(:,jj/3) = cell2mat(satcue{jj}(:,8));  % conversion factor LR to AB for whole experiment
            val_SignedConfidenceAB(:,jj/3)=   val_2AFC(:,jj/3);
            binary= ((val_SignedConfidenceAB(:,jj/3)>0) ==lr_to_ab3(:,jj/3));    % compute when Right chosen (positive answer). Thus, ig cnversion factor 1 and oout pos, then R chosen, and R=B, so B chosen. Whenever binary==1, B chosen.
            val_SignedConfidenceAB(binary,jj/3) = abs(val_SignedConfidenceAB(binary,jj/3));  % convert B is positive.                        if negative (L chosen) then 0 in binary, and factor 0, then L=B, thus B chosen. Thus whenever Conv factor == to find(val_2AFC>0), then B chosen.
            val_SignedConfidenceAB(~binary,jj/3)= -abs(val_SignedConfidenceAB(~binary,jj/3)); %  A in negative.                              The opposite as the previous, if pos out, R chosen, and factor 0, then L=B, A chosen, and negative as 1 ≃0
               
            
            for i=1:20
                
                indL = find(cell2mat(satcue{jj}(:,2))==i);            % when this stimulus was on left
                indR = find(cell2mat(satcue{jj}(:,3))==i);           % when this stimulus was on right
                lr_to_ab = cell2mat(satcue{jj}(indL,8));              % conversion factor from left_rigth to A_B, this factor saved on the column 8
                lr_to_ab2 = cell2mat(satcue{jj}(indR,8));            % conversion factor from left_rigth to A_B
                
                lr_to_cu = cell2mat(satcue{jj}(:,14))>0;              % conversion factor from left_rigth to A_B, just if positive or negative the output, 1 if right chosen, zero if left chosen
                
                
                indA = cat(1,((lr_to_ab) .* indL), (~lr_to_ab2).* indR);  indA(indA==0)= [];  %indices for A
                indB = cat(1,((~lr_to_ab) .* indL), (lr_to_ab2).* indR);  indB(indB==0)= [];  %indices for B
                indC = cat(1,indL(~lr_to_cu(indL)), indR(lr_to_cu(indR)));   % indeces for chosen, 1 if right chosen, so 1*indr, and zero if left chosen, so ~0 * indl
                indU = cat(1,indL(lr_to_cu(indL)), indR(~lr_to_cu(indR)));
                
                valRat = ranking{jj-1}(i,4);                        % rating value (-300 to 300)
                valRank = ranking{jj}(i,3);                         % ranking (0 to 20)

                
                val_A_Rat(indA, jj/3) = valRat;
                val_B_Rat(indB, jj/3) = valRat;
                val_L_Rat(indL, jj/3)= valRat;
                val_R_Rat(indR, jj/3)= valRat;
                val_A_Rank(indA, jj/3) = valRank;
                val_B_Rank(indB, jj/3) = valRank;
                val_L_Rank(indL, jj/3) = valRank;
                val_R_Rank(indR, jj/3) = valRank;
                
                val_C_Rat(indC, jj/3) = valRat;
                val_U_Rat(indU, jj/3) = valRat;
                val_C_Rank(indC, jj/3) = valRank;
                val_U_Rank(indU, jj/3) = valRank;
                
 
                
                
                flavourA(indA,jj/3) = i<11;
                flavourB(indB,jj/3) = i<11;
                flavourC(indC,jj/3) = i<11;
                flavourU(indU,jj/3) = i<11;
                flavourL(indL,jj/3) = i<11;
                flavourR(indR,jj/3) = i<11;
                
                stimOrderA(indA,jj/3) = i;
                stimOrderB(indB,jj/3) = i;
                stimOrderC(indC,jj/3) = i;
                stimOrderU(indU,jj/3) = i;
                stimOrderL(indL,jj/3) = i;
                stimOrderR(indR,jj/3) = i;
                
                
                % visual features
                Red_A (indA,jj/3) = Red(i,1);
                Green_A (indA,jj/3) = Green(i,1);
                Blue_A (indA,jj/3) = Blue(i,1);
                Luminance_A (indA,jj/3) = Luminance(i,1);
                Contrast_A (indA,jj/3) = Contrast(i,1);
                Hue_A (indA,jj/3) = Hue(i,1);
                Saturation_A (indA,jj/3) = Saturation(i,1);
                
                Red_B (indB,jj/3) = Red(i,1);
                Green_B (indB,jj/3) = Green(i,1);
                Blue_B (indB,jj/3) = Blue(i,1);
                Luminance_B (indB,jj/3) = Luminance(i,1);
                Contrast_B (indB,jj/3) = Contrast(i,1);
                Hue_B (indB,jj/3) = Hue(i,1);
                Saturation_B (indB,jj/3) = Saturation(i,1);
                
                
                % nutritional features
                
                Kalorien_A (indA,jj/3) = Kalorien(i,1);
                Protein_A (indA,jj/3) = Protein(i,1);
                Kohlenhydrate_A (indA,jj/3) = Kohlenhydrate(i,1);
                Zucker_A (indA,jj/3) = Zucker(i,1);
                Fett_A (indA,jj/3) = Fett(i,1);
                Ballaststoffe_A (indA,jj/3) = Ballaststoffe(i,1);
                Salz_A (indA,jj/3) = Salz(i,1);
                
                
                Kalorien_B (indB,jj/3) = Kalorien(i,1);
                Protein_B (indB,jj/3) = Protein(i,1);
                Kohlenhydrate_B (indB,jj/3) = Kohlenhydrate(i,1);
                Zucker_B (indB,jj/3) = Zucker(i,1);
                Fett_B (indB,jj/3) = Fett(i,1);
                Ballaststoffe_B (indB,jj/3) = Ballaststoffe(i,1);
                Salz_B (indB,jj/3) = Salz(i,1);
            end
            
            
            
            % Extract RT
            
            shortRT_Rat(:,jj/3) = RT{jj-1}; % this is te RT variable
            shortRT_AFC(:,jj/3) = RT{jj};
                        
        end
        

       
        response_windowEnd={{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}};
        endline2= {{2, 2},{0, 0},{0, 0}};
        xlabelplot = {{'Time relative to stimulus A1 onset (s)','Time relative to response onset (s)'},{'Time relative to stimulus onset (s)','Time relative to stimulus onset (s)'} ,{'Time relative to response onset (s)','Time relative to response onset (s)'}};
            
        for unit= memberEC
                A21  = rast_mrg_sess{pic_A21(unit),12};
                A22  = rast_mrg_sess{pic_A22(unit),12};
                A11  = rast_mrg_sess{pic_A11(unit),12};
                A12  = rast_mrg_sess{pic_A12(unit),12};
                B1  = rast_mrg_sess{pic_B1(unit),12};
                B2  = rast_mrg_sess{pic_B2(unit),12};

                
                
                clear All1 All2
                for trialnum = 1: numel(A11)
                    All1{trialnum,1} = [A11{trialnum}(A11{trialnum} < 1200), (B1{trialnum}(B1{trialnum} >  0 & B1{trialnum} < 1200) + 1200) ,  (A21{trialnum}(A21{trialnum} >  0) + 2400)];
                    All2{trialnum,1} = [A12{trialnum}(A12{trialnum} < 1200), (B2{trialnum}(B2{trialnum} >  0 & B2{trialnum} < 1200) + 1200) ,  (A22{trialnum}(A22{trialnum} >  0) + 2400)];
                    assert(~any(diff(All1{trialnum,1})<0)) % check strictly increasing
                    assert(~any(diff(All2{trialnum,1})<0)) % check strictly increasing524
                    
                    if isequal(All1{trialnum,1} , zeros(0,0))
                        All1{trialnum,1} = zeros(1,0);
                    end
                    if isequal(All2{trialnum,1} , zeros(0,0))
                        All2{trialnum,1} = zeros(1,0);
                    end                    
                end
                
                
                
                
                % Do a gaussian kernel for instantaneous firing rate
                [mfrA11] = kernelgauss(A11, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                [mfrA12] = kernelgauss(A12, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                [mfrB1] = kernelgauss(B1, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                [mfrB2] = kernelgauss(B2, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                [mfrA21] = kernelgauss(A21, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                [mfrA22] = kernelgauss(A22, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                
        
                % Stack data to correct for timestamps error, just in case present
                steps= [newwinwosreg(1) :stepsize: newwinwosreg(2)-1];
                laststep= size(steps,2);
                presentation1 = find(steps==0);
                endpresentation1 = (presentation1)+(1200/stepsize)-1;
                presentation2 = (presentation1)+(1200/stepsize);
                endpresentation2= endpresentation1 +(1200/stepsize);
                presentation3= presentation2 +(1200/stepsize);
                totalfinal =  endpresentation2 + size(presentation1:laststep,2);
                mfrA11(:,presentation2:endpresentation2)= mfrB1(:,presentation1:endpresentation1);
                mfrA12(:,presentation2:endpresentation2)= mfrB2(:,presentation1:endpresentation1);                              
                mfrA11(:,presentation3:totalfinal)= mfrA21(:,presentation1:laststep);
                mfrA12(:,presentation3:totalfinal)= mfrA22(:,presentation1:laststep);

                
                % Plot the Figure 3          
                spikeshape{unit}= rast_mrg_sess{pic_rat1(unit),14};                                     
                plotFigure3('Confidence ',  unit, sess, mfrA11, mfrA12, All1, All2, abs(val_2AFC), abs(val_2AFC), val_2AFC, val_SignedConfidenceAB, nombre_region, x_ax_lim, response_windowEnd{3}, endline2{1},  spikeshape, folder_to_save2{region}, unitnum, xlabelplot{1}, {'', ''}, rangePlot, possible_colors, 1, 'Confidence')
        
                close all;            
        end
        
        
    end
    
end




function [histoconv2] = kernelgauss (data, kern_width, resolution, stepsize, analyseswindow, stdkernel)

t = analyseswindow(2)- analyseswindow(1);
histo = zeros(size(data,1),t/resolution) * nan;

for a = 1:size(data,1)
    histo(a,:) = histc(data{a,1},(-((t/2)-(resolution/2)):resolution:((t/2)-(resolution/2))));
end
mfr=1;


kern_length = kern_width * stdkernel;
kernel = normpdf(-kern_length:resolution:kern_length, 0, kern_width);
y= histo * ((1000/resolution)/mfr);

for kkk=1:size(histo, 1)
    histoconv(kkk,:)= conv(y(kkk,:), kernel, 'same');
end

histoconv2= histoconv(:,1:stepsize:end);






function plotFigure3(nombre, unit, sess, Stim1, Stim2,  Rast1, Rast2,  valueplot1, valueplot2, valueplot3, valueplot4,  nombre_region, x_ax_lim, response_windowEnd, endline2,   spikeshape, folder_to_save2, unitnum, xlabelplot, letra, rangePlot, possible_colors, legendtype ,  xlabelin)



dbstop if error


%% relativize in cm

figpos = [0,0, 30.6800, 28.4300];
pos2x = 16.6116; pos1x = 2.2860;
heightTicks = 0.5715; heightmfr = 7.2009; heightRaster = 7.2009;   heightVshape= 5.9150;
width =  12.1920;

posvshapey = 2; posmfry = posvshapey + heightVshape + 4;  
posrastery = posmfry + heightmfr + 0; 


spikeshapex = pos2x + width - 3.2918;
spikeshapey = posmfry + heightmfr -  2.0631;
spikepos = [spikeshapex  spikeshapey  3.6576    2.0];


legposx=  pos1x +0.8961; legposy = posmfry +  + heightmfr - 2.3214;
legendpos =  [ legposx    legposy    2.1702    2.0505];



posticksy = posrastery + heightRaster + 0.05; 

pos21    =  [pos2x    posmfry    width    heightmfr];
pos22    =  [pos2x    posticksy  width    heightTicks];
pos11    =  [pos1x    posmfry    width    heightmfr];
pos12    =  [pos1x    posticksy  width    heightTicks];


posRaster2 = [pos2x    posrastery    width    heightRaster];
posRaster1 = [pos1x    posrastery    width    heightRaster];
posVshape2 = [pos2x    posvshapey    width    heightVshape];
posVshape1 = [pos1x    posvshapey   width    heightVshape];


heightanot = posvshapey + heightVshape + 0.2401;
anot1pos = [3.6576-0.8 heightanot+0.8 6.5720 0.4286];
anot2pos = [9.7536-0.8 heightanot+0.8 6.5720 0.4286] ;
anot3pos = [18.2880-0.8 heightanot+0.8 6.5720 0.4286];
anot4pos =[24.3840-0.8 heightanot+0.8 6.5720 0.4286];




sizechosentype = 18;

close all
fig1 =  figure('color', 'w', 'units', 'centimeters', 'position', figpos, 'visible', 'off');


    markersize =6;
    
    sb1 = plot_Raster_only(Rast1, Rast2, valueplot1(:,1), valueplot1(:,2), x_ax_lim, response_windowEnd{2}, endline2{2},     xlabelplot{1},  letra{1},   fig1, posRaster1, rangePlot, 1, 0, possible_colors, markersize);
    sb2 = plot_Raster_only(Rast1(2:end), Rast2(2:end), valueplot1(1:end-1,1), valueplot1(1:end-1,2), x_ax_lim, response_windowEnd{2}, endline2{2},  xlabelplot{1},  letra{1},   fig1, posRaster2, rangePlot, 0, 0, possible_colors, markersize);
    
    set(sb1, 'units', 'centimeters', 'position', posRaster1);
    set(sb2, 'units', 'centimeters', 'position', posRaster2);
    [ylimtop] = max([sb1.YLim(2) ,sb2.YLim(2)]);
    
    sb1.YLim = [0, ylimtop];
    sb2.YLim = [0, ylimtop];
    sb1.Children(end-2).YData = [ylimtop, ylimtop];
    sb1.Children(end-2).YData = [ylimtop, ylimtop];
    sb1.Children(end-1).YData = [ylimtop, ylimtop];
    sb2.Children(end-1).YData = [ylimtop, ylimtop];
    sb2.Children(end).YData = [ylimtop, ylimtop];
    sb2.Children(end).YData = [ylimtop, ylimtop];
    

    
    
    
    [sb1, tickspos, leg1] = plot_sPSTH_only(Stim1, Stim2, valueplot2(:,1), valueplot2(:,2), x_ax_lim, response_windowEnd{2}, endline2{2},     xlabelplot{1},  letra{1},   fig1, pos11, rangePlot, 1,2, possible_colors, legendtype);
    [sb2, tickspos2] = plot_sPSTH_only(Stim1(2:end,:), Stim2(2:end,:), valueplot2(1:end-1,1), valueplot2(1:end-1,2), x_ax_lim, response_windowEnd{2}, endline2{2},  xlabelplot{1},  letra{1},   fig1, pos21, rangePlot, 0,2, possible_colors, legendtype);
    set(sb1, 'units', 'centimeters', 'position', pos11);
    set(sb2, 'units', 'centimeters', 'position', pos21);
    set(leg1, 'units', 'centimeters', 'position', legendpos);
    set(tickspos, 'units', 'centimeters', 'position', pos12);
    set(tickspos2, 'units', 'centimeters', 'position', pos22);
    
    [ylimtop] = max([sb1.YLim(2) ,sb2.YLim(2)]);
    
    sb1.YLim = [0, ylimtop];
    sb2.YLim = [0, ylimtop];
    sb1.Children(9).YData = [ylimtop, ylimtop];
    sb1.Children(10).YData = [ylimtop, ylimtop];
    sb1.Children(11).YData = [ylimtop, ylimtop];
    sb2.Children(9).YData = [ylimtop, ylimtop];
    sb2.Children(10).YData = [ylimtop, ylimtop];
    sb2.Children(11).YData = [ylimtop, ylimtop];
    
    
    denseaxes=axes( 'units', 'centimeters', 'Position', spikepos);
    
    density_plot(spikeshape{unit});
    set(denseaxes, 'YAxisLocation', 'right')
    
    denseaxes.FontSize = 9;
    denseaxes.FontWeight = 'bold';
    
    
    
    
    sb1 = plot_Confidence_Vshape(Stim1, Stim2, valueplot3(:,1),  valueplot3(:,2),  xlabelin,  fig1, pos11,  1, possible_colors,  legendtype);
    sb2 = plot_Confidence_Vshape(Stim1, Stim2, valueplot4(:,1),  valueplot4(:,2),  xlabelin,  fig1, pos21, 0, possible_colors, legendtype);
    set(sb1, 'units', 'centimeters', 'position', posVshape1);
    set(sb2, 'units', 'centimeters', 'position', posVshape2);
    
    
    [ylimtop] = max([sb1.YLim(2) ,sb2.YLim(2)]);
    
    sb1.YLim= [0,ylimtop];
    sb2.YLim= [0,ylimtop];
    
    for gg= 1:numel (sb1.Children)
        if isequal(sb1.Children(gg).XData, [4.5, 4.5])
            if isequal(sb1.Children(gg).Color , [0.6000    0.6000    0.6000])
                sb1.Children(gg).YData=[0,ylimtop];
            end
        end
    end
    
    
    for gg= 1:numel (sb2.Children)
        
        if isequal(sb2.Children(gg).XData, [4.5, 4.5])
            if isequal(sb2.Children(gg).Color , [0.6000    0.6000    0.6000])
                sb2.Children(gg).YData=[0,ylimtop];
            end
        end
    end
    
    
    annotation1 = annotation('textbox',[0.12 0.858+0.03 0.15 0.05],'string','Left chosen','Units','normalized','FontName','Helvetica','FontSize',sizechosentype,'FontWeight','bold','EdgeColor','none');
    annotation2 = annotation('textbox',[0.32 0.858+0.03 0.15 0.05],'string','Right chosen','Units','normalized','FontName','Helvetica','FontSize',sizechosentype,'FontWeight','bold','EdgeColor','none');
    annotation3 = annotation('textbox',[0.6 0.858+0.03 0.15 0.05],'string','A chosen','Units','normalized','FontName','Helvetica','FontSize',sizechosentype,'FontWeight','bold','EdgeColor','none');
    annotation4 = annotation('textbox',[0.8 0.858+0.03 0.15 0.05],'string','B chosen','Units','normalized','FontName','Helvetica','FontSize',sizechosentype,'FontWeight','bold','EdgeColor','none');
    
    
    
    set(annotation1, 'units', 'centimeters', 'position', anot1pos)
    set(annotation2, 'units', 'centimeters', 'position', anot2pos)
    set(annotation3, 'units', 'centimeters', 'position', anot3pos)
    set(annotation4, 'units', 'centimeters', 'position', anot4pos)

    
    % Save Figure 3
    name_region
    print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ',  Whole plot 300'], '-dpng', '-r300', '-opengl')  
    print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ',  Whole plot2 300'], '-depsc', '-r300', '-opengl')
    set(fig1,'Units','centimeters');
    pos = get(fig1,'Position');
    set(fig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ',  Whole plot 300'], '-dpdf', '-r300', '-opengl')
    
  

