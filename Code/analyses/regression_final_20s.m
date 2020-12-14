function regression_final_20s(windowtype, modeltypeIn, typeCorrIn, extraname, windowsize, stepsize, inputRegion2, autoregresive_corr, a, lignedtoDecision)

% This funciton does the main analysis but for a window of 20 seconds
% example  regression_final_20s(10, 1, 'Spearman', '_InversesNEWVARIABLES', 200, 10, {'PHC', 'A', 'EC', 'H'}, 0)

tic
dbstop if error

% settings

%% extra settings
global modeltype typeCorr KNNparameter
KNNparameter = 5;

if modeltypeIn<=4
    modeltype=modeltypeIn; % if ==1 then duplicate corr, if ==2 then mixed-model, if ==3 is a pooled-regression method, ==4 is mutualinformation
else
    error('modeltype must be smaller or equal to 3 and')
end

if ismember(typeCorrIn, {'Spearman', 'Pearson', 'Kendall',  'MutualInfo'})
    typeCorr=typeCorrIn;
else
    error('Type of correlation, must be a string and a valid type of correlation, i.e, Spearman, Pearson or Kendall')
end

if nargin < 4
    extraname = '';
end

if windowtype ==10
    resolution = 1;
    
    newwinwosreg = [-20000,20000];
    newwinwosregRating = [-20000,20000];
    newwindRespAFC = [-20000,20000];
    newwindRespRat = [-20000,20000];
    
    newwinwosreg = [-20000,20000];
    newwinwosregRating = [-20000,20000];
    newwindRespAFC = [-20000,20000];
    newwindRespRat = [-20000,20000];
    stdkernel= 5;
    
elseif windowtype == 20
    windowsize = 0;
    windowsizeRat = 0;
    
    wRat = [0,2000];    wRat2 = wRat + 200;  wRat3=[0,5000];
    
    wA1 = [0,1000];
    wB  = [1200,2200];
    wA2 = [2400,3400];
    wA12 = wA1 + 200;
    wB2  = wB + 200;
    wA22 = wA2+ 200;
    
    wafterB= [1200,5950];
    wafterA2=[3400,5950];
    wafterA2_200= [3600,5950];
    
    newwinwosreg =[wA1; wB; wA2; wA12; wB2; wA22; wafterB; wafterA2; wafterA2_200];
    newwinwosregRating =[wRat;wRat2];
    
elseif windowtype == 30
    resolution = 1;
    stdkernel= 5;
    newwinwosreg = [-2350,3600+2350];
    newwinwosregRating = [-2350,2000+2350];
   
end




%% Consumed products y patients
Consumedindex1 = [    9,    11,     20,     2,    14,     2,    16,   nan,     16,     4,     7,   17];
Consumedindex2 = {'001', '002',  '003', '004', '005', '006', '007', '008',  '009', '010', '011','012'};
Consumedindex3 = [nan(1,27), 9,nan(1,2)  11, 20,  2, 14, 2, nan(1,1), 16, nan(1,1), nan, 16,  nan(1,1), 4, nan(1,1), 7,nan(1,1),  17];
Consumedindex4 = {{'001',9}, {'002', 11},  {'003', 20}, {'004', 2}, {'005', 14}, {'006', 2}, {'007', 16}, ...
    {'008', nan}, {'009', 16}, {'010', 4}, {'011', 7}, {'012', 17}};

% folder and subfolders to save results and load data
satiationfolder = '/media/raid/Alex/Reclustered analysis All 2/';
cd(satiationfolder);


if autoregresive_corr == 0
    ana_title= ['20s_',num2str(modeltypeIn), '_',typeCorrIn,'_window',num2str(windowtype),'_kernel', num2str(windowsize), 'ms_step' , num2str(stepsize), 'ms_',extraname];
elseif autoregresive_corr == 1
    ana_title= ['20s_',num2str(modeltypeIn), '_',typeCorrIn,'_window',num2str(windowtype),'_kernel', num2str(windowsize), 'ms_step' , num2str(stepsize), 'ms_',extraname, '_autocorr'];
end


% alignedtoDecision ==0  means aligend to start of the trial, ==1 aligned
% to binary decision, ==2 aligned to OK pus for confidence report
if alignedtoDecision == 1
    ana_title= [ana_title, '_LockettoBinaryDecision'];
elseif alignedtoDecision == 2
    ana_title= [ana_title, '_LockettoOKDecision'];
end

mkdir([satiationfolder, ana_title]);
folder_to_save = [satiationfolder, ana_title];
disp('Folder to save the dat is:')
folder_to_save
bslwind = [-400, 0];


namefolderrasts= 'Rasts matfiles';
for l=1:length(inputRegion2)
    mkdir([folder_to_save, '/' char(inputRegion2{l})])
    folder_to_save2{l}=[folder_to_save, '/' char(inputRegion2{l})];
end



%% start analysis
for hhh=1:length(inputRegion2)
    clear unitnum
    inputRegion=inputRegion2{hhh};
    inputRegion=char(inputRegion);
    behavDir=[];
    load ([satiationfolder,'/', namefolderrasts, '/',char(inputRegion),'_rasts6.mat'], 'rasts', 'outfiles')
    for o=1:size(outfiles,1)
        behavDir{end+1}= outfiles{o,2};
    end
    
    %% regional analysis
    currentdirectory = pwd;
    %some changes for parallel processing
    lastnumber=1;
    for sess=1:size(outfiles,1)
        
        unitind=ismember(cell2mat(rasts(:,2:6)),cell2mat(outfiles(sess,2:6)),'rows'); % check if the unit is part of the channel, compares the name of the patient, etc and region
        n_ses_units = sum(unitind & cell2mat(rasts(:,9)) == 7 & cell2mat(rasts(:,10)) == 3);
        
        for a=1:n_ses_units
            unitnum{sess}(a,1)= lastnumber;
            lastnumber=lastnumber+1;
        end
    end
    
    
    for sess = 1:size(outfiles,1)
        
        cd(currentdirectory)
        
        
        unitind=ismember(cell2mat(rasts(:,2:6)),cell2mat(outfiles(sess,2:6)),'rows'); % check if the unit is part of the channel, compares the name of the patient, etc and region
        actualregion=char(outfiles(sess,6));
        
        PIDstr = (outfiles{sess,2});
        PID = str2num(PIDstr); % PID
        
        
        Consumedindex = Consumedindex3(PID); % ID of the consumed product
        ConsumedindexN = 	Consumedindex1(cellfun(@(s) ~isempty(strfind(PIDstr, s)), Consumedindex2));
        
        
        n_ses_units = sum(unitind & cell2mat(rasts(:,9)) == 7 & cell2mat(rasts(:,10)) == 3); % Unique number of units in rast_mrg recorded in this session and this region
        rast_mrg_sess = rasts(unitind ,:); % Units in rast_mrg recorded in this session and this region
        
        ttl_n_cond = cell2mat(rast_mrg_sess(:,9:10)); %icludes when the event (7 or 43 or 44) and the criteria (0-20)
        
        pic_rat1 = find((ttl_n_cond(:,1) == 7) & (ttl_n_cond(:,2) == 3));     %  First block judge pics
        pic_rat2 = find((ttl_n_cond(:,1) == 7) & (ttl_n_cond(:,2) == 4));     %  Second block judge pics
        
        pic_A11 = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 11)));   % Pic A1
        pic_A12 = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 12)));   % Pic A1
        pic_B1  = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 13)));   % Pic B
        pic_B2  = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 14)));   % Pic B
        pic_A21  = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 15)));
        pic_A22  = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 16)));
        
        
        [satcue,ranking,events]=load_behavioural_data(behavDir{sess});
        [RT,endResponse] = extract_reaction_times(events);
        [Red, Green, Blue, Luminance, Contrast, Hue, Saturation, ~ ] = extract_visualFeatures(outfiles{1,2});
        [Kalorien, Protein, Kohlenhydrate, Zucker, Fett, Ballaststoffe, Salz] = naehwerte(outfiles{1,2});
        
        cd ..
        
        
        for jj=[3,6]
            %   Output values (real answers)
            
            val_rat(:,jj/3)= cell2mat(satcue{jj-1}(:,14));
            val_2AFC(:,jj/3)= cell2mat(satcue{jj}(:,14));
            
            
            %   Val by Product
            RatProd (:,jj/3) = ranking{jj-1}(:,4);
            RankProd(:,jj/3) = ranking{jj}(:,3);
            
            if jj==3
                invRatProd (:,jj/3) = ranking{6-1}(:,4);
                invRankProd(:,jj/3) = ranking{6}(:,3);
            elseif jj==6
                invRatProd (:,jj/3) = ranking{3-1}(:,4);
                invRankProd(:,jj/3) = ranking{3}(:,3);
            end
            
            
            %   Input values
            %   Rating paradigm
            
            val_rat_Rat(:,jj/3)=nan(length(satcue{jj-1}(:,14)),1);
            val_rat_Rank(:,jj/3)=nan(length(satcue{jj-1}(:,14)),1);
            
            for i=1:20
                
                ind3=find(cell2mat(satcue{jj-1}(:,2))==i);
                valRat=  ranking{jj-1}(i,4); % ranking
                valRank= ranking{jj}(i,3);
                
                
                if jj==3
                    invvalRat = ranking{5}(i,4);                        % rating value (-300 to 300)
                    invvalRank = ranking{6}(i,3);                         % ranking (0 to 20)
                    
                    indinv= find(cell2mat(satcue{5}(:,2))==i); assert(length(indinv)==3);
                    invRT = RT{5}(indinv,1);
                    
                elseif jj==6
                    invvalRat = ranking{2}(i,4);                        % rating value (-300 to 300)
                    invvalRank = ranking{3}(i,3);                         % ranking (0 to 20)
                    
                    indinv= find(cell2mat(satcue{2}(:,2))==i); assert(length(indinv)==3);
                    invRT = RT{2}(indinv,1);
                    
                end
                
                
                
                
                val_rat_Rat(ind3,jj/3)= valRat;
                val_rat_Rank(ind3,jj/3)= valRank;
                invval_rat_Rat(ind3,jj/3)= invvalRat;
                invval_rat_Rank(ind3,jj/3)= invvalRank;
                
                
                stimOrderRat(ind3,jj/3) = i;
                flavourRat(ind3, jj/3) = i<11;
                
                %                     visual features
                Red_Rat (ind3,jj/3) = Red(i,1);
                Green_Rat (ind3,jj/3) = Green(i,1);
                Blue_Rat (ind3,jj/3) = Blue(i,1);
                Luminance_Rat (ind3,jj/3) = Luminance(i,1);
                Contrast_Rat (ind3,jj/3) = Contrast(i,1);
                Hue_Rat (ind3,jj/3) = Hue(i,1);
                Saturation_Rat (ind3,jj/3) = Saturation(i,1);
                
                %                     nutritional  features
                Kalorien_Rat (ind3,jj/3) = Kalorien(i,1);
                Protein_Rat (ind3,jj/3) = Protein(i,1);
                Kohlenhydrate_Rat (ind3,jj/3) = Kohlenhydrate(i,1);
                Zucker_Rat (ind3,jj/3) = Zucker(i,1);
                Fett_Rat (ind3,jj/3) = Fett(i,1);
                Ballaststoffe_Rat (ind3,jj/3) = Ballaststoffe(i,1);
                Salz_Rat (ind3,jj/3) = Salz(i,1);
                
                
                
                %% fichaxes recentes para a plantilla
                
                RTmean_Rat(ind3, jj/3)  =  mean(RT{jj-1}(ind3,1)); % mean RT for satiation effect in first task
                invRTmean_Rat(ind3, jj/3)  =  mean(invRT); % mean RT for satiation effect in first task
                
            end
            
            
            %2AFC paradigm
            
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
            
            
            %% transform the output of the likert and convert it so that the negative is A chosen and positive B chosen.
            lr_to_ab3(:,jj/3) = cell2mat(satcue{jj}(:,8));  % conversion factor LR to AB for whole experiment
            val_diffAB(:,jj/3)=   val_2AFC(:,jj/3);
            binary= ((val_diffAB(:,jj/3)>0) ==lr_to_ab3(:,jj/3));    % compute when Right chosen (positive answer). Thus, ig cnversion factor 1 and oout pos, then R chosen, and R=B, so B chosen. Whenever binary==1, B chosen.
            val_diffAB(binary,jj/3) = abs(val_diffAB(binary,jj/3));  % convert B is positive.                        if negative (L chosen) then 0 in binary, and factor 0, then L=B, thus B chosen. Thus whenever Conv factor == to find(val_2AFC>0), then B chosen.
            val_diffAB(~binary,jj/3)= -abs(val_diffAB(~binary,jj/3)); %  A in negative.                              The opposite as the previous, if pos out, R chosen, and factor 0, then L=B, A chosen, and negative as 1 ≃0
            %                 finally if output negative, L chosen, and factor 1, then L==A, so A chosen again, negative here as 1 ≃0.
            
            
            
            
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
                
                if jj==3
                    invvalRat = ranking{5}(i,4);                        % rating value (-300 to 300)
                    invvalRank = ranking{6}(i,3);                         % ranking (0 to 20)
                elseif jj==6
                    invvalRat = ranking{2}(i,4);                        % rating value (-300 to 300)
                    invvalRank = ranking{3}(i,3);                         % ranking (0 to 20)
                end
                
                
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
                
                
                
                %                     inverted order
                
                invval_A_Rat(indA, jj/3) = invvalRat;
                invval_B_Rat(indB, jj/3) = invvalRat;
                invval_L_Rat(indL, jj/3)=  invvalRat;
                invval_R_Rat(indR, jj/3)=  invvalRat;
                invval_A_Rank(indA, jj/3) = invvalRank;
                invval_B_Rank(indB, jj/3) = invvalRank;
                invval_L_Rank(indL, jj/3) = invvalRank;
                invval_R_Rank(indR, jj/3) = invvalRank;
                
                invval_C_Rat(indC, jj/3) = invvalRat;
                invval_U_Rat(indU, jj/3) = invvalRat;
                invval_C_Rank(indC, jj/3) = invvalRank;
                invval_U_Rank(indU, jj/3) = invvalRank;
                
                
                
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
                
                
                %                     visual features
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
                
                
                Red_L (indL,jj/3) = Red(i,1);
                Green_L (indL,jj/3) = Green(i,1);
                Blue_L (indL,jj/3) = Blue(i,1);
                Luminance_L (indL,jj/3) = Luminance(i,1);
                Contrast_L (indL,jj/3) = Contrast(i,1);
                Hue_L (indL,jj/3) = Hue(i,1);
                Saturation_L (indL,jj/3) = Saturation(i,1);
                
                
                
                Red_R (indR,jj/3) = Red(i,1);
                Green_R (indR,jj/3) = Green(i,1);
                Blue_R (indR,jj/3) = Blue(i,1);
                Luminance_R (indR,jj/3) = Luminance(i,1);
                Contrast_R (indR,jj/3) = Contrast(i,1);
                Hue_R (indR,jj/3) = Hue(i,1);
                Saturation_R (indR,jj/3) = Saturation(i,1);
                
                
                
                Red_C (indC,jj/3) = Red(i,1);
                Green_C (indC,jj/3) = Green(i,1);
                Blue_C (indC,jj/3) = Blue(i,1);
                Luminance_C (indC,jj/3) = Luminance(i,1);
                Contrast_C (indC,jj/3) = Contrast(i,1);
                Hue_C (indC,jj/3) = Hue(i,1);
                Saturation_C (indC,jj/3) = Saturation(i,1);
                
                Red_U (indU,jj/3) = Red(i,1);
                Green_U (indU,jj/3) = Green(i,1);
                Blue_U (indU,jj/3) = Blue(i,1);
                Luminance_U (indU,jj/3) = Luminance(i,1);
                Contrast_U (indU,jj/3) = Contrast(i,1);
                Hue_U (indU,jj/3) = Hue(i,1);
                Saturation_U (indU,jj/3) = Saturation(i,1);
                
                
                %                     nutritional features
                
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
                
                
                
                Kalorien_L (indL,jj/3) = Kalorien(i,1);
                Protein_L (indL,jj/3) = Protein(i,1);
                Kohlenhydrate_L (indL,jj/3) = Kohlenhydrate(i,1);
                Zucker_L (indL,jj/3) = Zucker(i,1);
                Fett_L (indL,jj/3) = Fett(i,1);
                Ballaststoffe_L (indL,jj/3) = Ballaststoffe(i,1);
                Salz_L (indA,jj/3) = Salz(i,1);
                
                Kalorien_R (indR,jj/3) = Kalorien(i,1);
                Protein_R (indR,jj/3) = Protein(i,1);
                Kohlenhydrate_R (indR,jj/3) = Kohlenhydrate(i,1);
                Zucker_R (indR,jj/3) = Zucker(i,1);
                Fett_R (indR,jj/3) = Fett(i,1);
                Ballaststoffe_R (indR,jj/3) = Ballaststoffe(i,1);
                Salz_R (indR,jj/3) = Salz(i,1);
                
                
                Kalorien_C (indC,jj/3) = Kalorien(i,1);
                Protein_C (indC,jj/3) = Protein(i,1);
                Kohlenhydrate_C (indC,jj/3) = Kohlenhydrate(i,1);
                Zucker_C (indC,jj/3) = Zucker(i,1);
                Fett_C (indC,jj/3) = Fett(i,1);
                Ballaststoffe_C (indC,jj/3) = Ballaststoffe(i,1);
                Salz_C (indC,jj/3) = Salz(i,1);
                
                
                Kalorien_U (indU,jj/3) = Kalorien(i,1);
                Protein_U (indU,jj/3) = Protein(i,1);
                Kohlenhydrate_U (indU,jj/3) = Kohlenhydrate(i,1);
                Zucker_U (indU,jj/3) = Zucker(i,1);
                Fett_U (indU,jj/3) = Fett(i,1);
                Ballaststoffe_U (indU,jj/3) = Ballaststoffe(i,1);
                Salz_U (indU,jj/3) = Salz(i,1);
                
            end
            
            %                 Extract RT
            
            longRT_Rat(:,jj/3) = endResponse{jj-1};
            longRT_AFC(:,jj/3) = endResponse{jj};
            
            shortRT_Rat(:,jj/3) = RT{jj-1};
            shortRT_AFC(:,jj/3) = RT{jj};
            
            
            otherRT_Rat(:,jj/3)= cell2mat(satcue{jj-1}(:,13))*1000;
            otherRT_AFC(:,jj/3)= cell2mat(satcue{jj}(:,13))*1000;
            
            diffRT_Rat(:,jj/3)= otherRT_Rat(:,jj/3) - shortRT_Rat(:,jj/3);
            diffRT_AFC(:,jj/3)= otherRT_AFC(:,jj/3) - shortRT_AFC(:,jj/3);
            
            
            
        end
        
        
        %                 inverted likert scale output
        for i= 1:20
            for hh=1:20
                if i ~=hh
                    try
                        indPair1 = find(cell2mat(satcue{3}(:,2))==i & cell2mat(satcue{3}(:,3))==hh);  indPair2 =   find(cell2mat(satcue{6}(:,2))==i & cell2mat(satcue{6}(:,3))==hh);
                        
                    end
                    try
                        indPair1 = find(cell2mat(satcue{3}(:,3))==i & cell2mat(satcue{3}(:,2))==hh);  indPair2 =   find(cell2mat(satcue{6}(:,3))==i & cell2mat(satcue{6}(:,2))==hh);
                    end
                    
                    invval_2AFC(indPair1, 1) = val_2AFC(indPair2,2);
                    invval_diffAB (indPair1, 1) = val_diffAB (indPair2,2);
                    invshortRT_AFC (indPair1,1) = shortRT_AFC (indPair2,2);
                    
                    
                    invval_2AFC (indPair2, 2) = val_2AFC(indPair1,1);
                    invval_diffAB (indPair2, 2) = val_diffAB (indPair1,1);
                    invshortRT_AFC (indPair2,2) = shortRT_AFC (indPair1,1);
                    
                    
                    invval2_C_Rat(indPair1,1) =  val_C_Rat(indPair2,2);  % the other inverse chosen and unchosen are based in the stimulus vbalue before and after,now it is based
                    invval2_C_Rat(indPair2,2) =  val_C_Rat(indPair1,1);  %inthe value of the chosen value in the trial type, so if the stimulus chosen changes, it is the value of 									          %the other stimulus when chosen.
                    invval2_C_Rank(indPair1,1) =  val_C_Rank(indPair2,2);
                    invval2_C_Rank(indPair2,2) =  val_C_Rank(indPair1,1);
                    
                    invval2_U_Rat(indPair1,1) =  val_U_Rat(indPair2,2);
                    invval2_U_Rat(indPair2,2) =  val_U_Rat(indPair1,1);
                    invval2_U_Rank(indPair1,1) =  val_U_Rank(indPair2,2);
                    invval2_U_Rank(indPair2,2) =  val_U_Rank(indPair1,1);
                    
                    
                end
            end
        end
        
        
        
        if autoregresive_corr==1
            
            val_rat=     [ val_rat(1:end-1, :)];
            val_2AFC =   [ val_2AFC(1:end-1, :)];
            val_diffAB = [ val_diffAB(1:end-1, :)];
            
            val_C_Rat=  [ val_C_Rat(1:end-1, :)];
            val_U_Rat=  [ val_U_Rat(1:end-1, :)];
            val_A_Rat=  [ val_A_Rat(1:end-1, :)];
            val_B_Rat=  [ val_B_Rat(1:end-1, :)];
            val_L_Rat=  [ val_L_Rat(1:end-1, :)];
            val_R_Rat=  [ val_R_Rat(1:end-1, :)];
            
            
            flavourC=  [ flavourC(1:end-1, :)];
            flavourU=  [ flavourU(1:end-1, :)];
            flavourA=  [ flavourA(1:end-1, :)];
            flavourB=  [ flavourB(1:end-1, :)];
            flavourL=  [ flavourL(1:end-1, :)];
            flavourR=  [ flavourR(1:end-1, :)];
            flavourRat=[ flavourRat(1:end-1, :)];
            
            
            shortRT_AFC= [ shortRT_AFC(1:end-1, :)];
            shortRT_Rat= [ shortRT_Rat(1:end-1, :)];
            diffRT_AFC= [ diffRT_AFC(1:end-1, :)];
            diffRT_Rat= [ diffRT_Rat(1:end-1, :)];
            
            
            val_C_Rank=  [ val_C_Rank(1:end-1, :)];
            val_U_Rank=  [ val_U_Rank(1:end-1, :)];
            val_A_Rank=  [ val_A_Rank(1:end-1, :)];
            val_B_Rank=  [ val_B_Rank(1:end-1, :)];
            val_L_Rank=  [ val_L_Rank(1:end-1, :)];
            val_R_Rank=  [ val_R_Rank(1:end-1, :)];
            
            
            val_rat_Rat= [ val_rat_Rat(1:end-1, :)];
            val_rat_Rank= [ val_rat_Rank(1:end-1, :)];
            invval_rat_Rat=invval_rat_Rat(1:end-1, :);
            invval_rat_Rank=invval_rat_Rank(1:end-1, :);
            
            
            
            
            
            invval_2AFC =  invval_2AFC(1:end-1, :);
            invval_diffAB =  invval_diffAB(1:end-1, :);
            invshortRT_AFC =  invshortRT_AFC(1:end-1, :);
            
            
            invval_A_Rat =  invval_A_Rat(1:end-1, :);
            invval_B_Rat =  invval_B_Rat(1:end-1, :);
            invval_L_Rat =  invval_L_Rat(1:end-1, :);
            invval_R_Rat =  invval_R_Rat(1:end-1, :);
            invval_C_Rat =  invval_C_Rat(1:end-1, :);
            invval_U_Rat =  invval_U_Rat(1:end-1, :);
            
            invval_A_Rank =  invval_A_Rank(1:end-1, :);
            invval_B_Rank =  invval_B_Rank(1:end-1, :);
            invval_L_Rank =  invval_L_Rank(1:end-1, :);
            invval_R_Rank =  invval_R_Rank(1:end-1, :);
            invval_C_Rank =  invval_C_Rank(1:end-1, :);
            invval_U_Rank =  invval_U_Rank(1:end-1, :);
            
            
            
            %                     a new type of inverse  for chosen and unchosen, based on trial type and not on stimulus identity
            invval2_C_Rank =  invval2_C_Rank(1:end-1, :);
            invval2_U_Rank =  invval2_U_Rank(1:end-1, :);
            invval2_C_Rat =  invval2_C_Rat(1:end-1, :);
            invval2_U_Rat =  invval2_U_Rat(1:end-1, :);
            
            
            
            invRTmean_Rat = invRTmean_Rat(1:end-1, :);
            
            
            RTmean_Rat = RTmean_Rat(1:end-1, :);
            
            
            Kalorien_Rat = Kalorien_Rat(1:end-1, :);
            Protein_Rat = Protein_Rat(1:end-1, :);
            Kohlenhydrate_Rat = Kohlenhydrate_Rat(1:end-1, :);
            Zucker_Rat = Zucker_Rat(1:end-1, :);
            Fett_Rat = Fett_Rat(1:end-1, :);
            Ballaststoffe_Rat = Ballaststoffe_Rat(1:end-1, :);
            Salz_Rat = Salz_Rat(1:end-1, :);
            
            Kalorien_B = Kalorien_B(1:end-1, :);
            Protein_B = Protein_B(1:end-1, :);
            Kohlenhydrate_B = Kohlenhydrate_B(1:end-1, :);
            Zucker_B = Zucker_B(1:end-1, :);
            Fett_B = Fett_B(1:end-1, :);
            Ballaststoffe_B = Ballaststoffe_B(1:end-1, :);
            Salz_B = Salz_B(1:end-1, :);
            
            Kalorien_A = Kalorien_A(1:end-1, :);
            Protein_A = Protein_A(1:end-1, :);
            Kohlenhydrate_A = Kohlenhydrate_A(1:end-1, :);
            Zucker_A = Zucker_A(1:end-1, :);
            Fett_A = Fett_A(1:end-1, :);
            Ballaststoffe_A = Ballaststoffe_A(1:end-1, :);
            Salz_A = Salz_A(1:end-1, :);
            
            
            
            stimOrderRat = stimOrderRat(1:end-1, :);
            Red_Rat = Red_Rat(1:end-1, :);
            Green_Rat = Green_Rat(1:end-1, :);
            Blue_Rat = Blue_Rat(1:end-1, :);
            Contrast_Rat = Contrast_Rat(1:end-1, :);
            Luminance_Rat = Luminance_Rat(1:end-1, :);
            Hue_Rat = Hue_Rat(1:end-1, :);
            Saturation_Rat = Saturation_Rat(1:end-1, :);
            
            stimOrderB = stimOrderB(1:end-1, :);
            Red_B = Red_B(1:end-1, :);
            Green_B = Green_B(1:end-1, :);
            Blue_B = Blue_B(1:end-1, :);
            Contrast_B = Contrast_B(1:end-1, :);
            Luminance_B = Luminance_B(1:end-1, :);
            Hue_B = Hue_B(1:end-1, :);
            Saturation_B = Saturation_B(1:end-1, :);
            
            stimOrderA = stimOrderA(1:end-1, :);
            Red_A = Red_A(1:end-1, :);
            Green_A = Green_A(1:end-1, :);
            Blue_A = Blue_A(1:end-1, :);
            Contrast_A = Contrast_A(1:end-1, :);
            Luminance_A = Luminance_A(1:end-1, :);
            Hue_A = Hue_A(1:end-1, :);
            Saturation_A= Saturation_A(1:end-1, :);
            
            
            
            stimOrderR = stimOrderR(1:end-1, :);
            Red_R = Red_R(1:end-1, :);
            Green_R = Green_R(1:end-1, :);
            Blue_R = Blue_R(1:end-1, :);
            Contrast_R = Contrast_R(1:end-1, :);
            Luminance_R = Luminance_R(1:end-1, :);
            Hue_R = Hue_R(1:end-1, :);
            Saturation_R = Saturation_R(1:end-1, :);
            
            stimOrderL = stimOrderL(1:end-1, :);
            Red_L = Red_L(1:end-1, :);
            Green_L = Green_L(1:end-1, :);
            Blue_L = Blue_L(1:end-1, :);
            Contrast_L = Contrast_L(1:end-1, :);
            Luminance_L = Luminance_L(1:end-1, :);
            Hue_L = Hue_L(1:end-1, :);
            Saturation_L= Saturation_L(1:end-1, :);
            
            
            stimOrderU = stimOrderU(1:end-1, :);
            Red_U = Red_U(1:end-1, :);
            Green_U = Green_U(1:end-1, :);
            Blue_U = Blue_U(1:end-1, :);
            Contrast_U = Contrast_U(1:end-1, :);
            Luminance_U = Luminance_U(1:end-1, :);
            Hue_U = Hue_U(1:end-1,:);
            Saturation_U = Saturation_U(1:end-1, :);
            
            stimOrderC = stimOrderC(1:end-1, :);
            Red_C = Red_C(1:end-1, :);
            Green_C = Green_C(1:end-1, :);
            Blue_C = Blue_C(1:end-1, :);
            Contrast_C = Contrast_C(1:end-1, :);
            Luminance_C = Luminance_C(1:end-1, :);
            Hue_C = Hue_C(1:end-1, :);
            Saturation_C= Saturation_C(1:end-1, :);
            
            
            
            
            
            Kalorien_U = Kalorien_U(1:end-1, :);
            Protein_U = Protein_U(1:end-1, :);
            Kohlenhydrate_U = Kohlenhydrate_U(1:end-1, :);
            Zucker_U = Zucker_U(1:end-1, :);
            Fett_U = Fett_U(1:end-1, :);
            Ballaststoffe_U = Ballaststoffe_U(1:end-1, :);
            Salz_U = Salz_U(1:end-1, :);
            
            Kalorien_L = Kalorien_L(1:end-1, :);
            Protein_L = Protein_L(1:end-1, :);
            Kohlenhydrate_L = Kohlenhydrate_L(1:end-1, :);
            Zucker_L = Zucker_L(1:end-1, :);
            Fett_L = Fett_L(1:end-1, :);
            Ballaststoffe_L = Ballaststoffe_L(1:end-1, :);
            Salz_L = Salz_L(1:end-1, :);
            
            
            
            Kalorien_R = Kalorien_R(1:end-1, :);
            Protein_R = Protein_R(1:end-1, :);
            Kohlenhydrate_R = Kohlenhydrate_R(1:end-1, :);
            Zucker_R = Zucker_R(1:end-1, :);
            Fett_R = Fett_R(1:end-1, :);
            Ballaststoffe_R = Ballaststoffe_R(1:end-1, :);
            Salz_R = Salz_R(1:end-1, :);
            
            Kalorien_C = Kalorien_C(1:end-1, :);
            Protein_C = Protein_C(1:end-1, :);
            Kohlenhydrate_C = Kohlenhydrate_C(1:end-1, :);
            Zucker_C = Zucker_C(1:end-1, :);
            Fett_C = Fett_C(1:end-1, :);
            Ballaststoffe_C = Ballaststoffe_C(1:end-1, :);
            Salz_C = Salz_C(1:end-1, :);
            
        end
        
        
        
        
        for unit = 1:n_ses_units
            
            disp(['unit number ' num2str(unitnum{sess}(unit,1)) ' out of ' num2str(size(rasts,1)/12) ' in region ' num2str(sess) ' out of ' num2str(size(outfiles,1)), ' unit ' num2str(unit) ' out of ' num2str(n_ses_units) ' in this region ' inputRegion ])
            
            
            %             A21  = rast_mrg_sess{pic_A21(unit),12};
            %             A22  = rast_mrg_sess{pic_A22(unit),12};
            A11  = rast_mrg_sess{pic_A11(unit),12};
            A12  = rast_mrg_sess{pic_A12(unit),12};
            %             B1  = rast_mrg_sess{pic_B1(unit),12};
            %             B2  = rast_mrg_sess{pic_B2(unit),12};
            %             likAFCR1  = rast_mrg_sess{lik_AFCR1(unit),12};
            %             likAFCR2  = rast_mrg_sess{lik_AFCR2(unit),12};
            ratPic1= rast_mrg_sess{pic_rat1(unit),12};
            ratPic2= rast_mrg_sess{pic_rat2(unit),12};
            %             ratLik1= rast_mrg_sess{lik_RatR1(unit),12};
            %             ratLik2= rast_mrg_sess{lik_RatR2(unit),12};
            %             bsl1 = rast_mrg_sess{baseline1(unit),12};
            %             bsl2 = rast_mrg_sess{baseline2(unit),12};
            
            
            if windowtype == 20
                
                mfrA11 = [];
                mfrA12 = [];
                mfrA21 = [];
                mfrA22 = [];
                mfrB1 = [];
                mfrB2 = [];
                mfrResp1 = [];
                mfrResp2 = [];
                
                mfrRat1 = [];
                mfrRat2 = [];
                mfrRatResp1 = [];
                mfrRatResp2 = [];
                
                for stim=1:60
                    
                    for gg=1:size(newwinwosregRating,1)
                        mfrRat1(stim,gg) = length(ratPic1{stim}(1,ratPic1{stim} > newwinwosregRating(gg,1) & ratPic1{stim} < newwinwosregRating(gg,2))) * (1000/abs(newwinwosregRating(gg,2) - newwinwosregRating(gg,1))); % Hz of this unit at each stimulus at this time window
                        mfrRat2(stim,gg) = length(ratPic2{stim}(1,ratPic2{stim} > newwinwosregRating(gg,1) & ratPic2{stim} < newwinwosregRating(gg,2))) * (1000/abs(newwinwosregRating(gg,2) - newwinwosregRating(gg,1))); % Hz of this unit at each stimulus at this time window
                        
                        mfrRatResp1(stim,gg) = length(ratLik1{stim}(1,ratLik1{stim} > newwinwosregRating(gg,1) & ratLik1{stim} < newwinwosregRating(gg,2))) * (1000/abs(newwinwosregRating(gg,2) - newwinwosregRating(gg,1))); % Hz of this unit at each stimulus at this time window
                        mfrRatResp2(stim,gg) = length(ratLik2{stim}(1,ratLik2{stim} > newwinwosregRating(gg,1) & ratLik2{stim} < newwinwosregRating(gg,2))) * (1000/abs(newwinwosregRating(gg,2) - newwinwosregRating(gg,1))); % Hz of this unit at each stimulus at this time window
                    end
                    if bsl_on
                        ratbsl1= length(ratPic1{stim}(1,ratPic1{stim}> bslwind(1) & ratPic1{stim}<bslwind(2)))* (1000/abs(bslwind(2) - bslwind(1)));
                        ratbsl2= length(ratPic2{stim}(1,ratPic2{stim}> bslwind(1) & ratPic2{stim}<bslwind(2)))* (1000/abs(bslwind(2) - bslwind(1)));
                        
                        mfrRat1(stim,:)= mfrRat1(stim,:) - ratbsl1;
                        mfrRat2(stim,:)= mfrRat2(stim,:) - ratbsl2;
                        mfrRatResp1(stim,:)= mfrRatResp1(stim,:) - ratbsl1;
                        mfrRatResp2(stim,:)= mfrRatResp2(stim,:) - ratbsl2;
                    end
                    
                    
                    if windowtype == 5
                        
                        mfrRat1(stim,gg+1) = length(ratLik1{stim}(1,ratLik1{stim} > rsp_t_likertRats(stim,1) & ratLik1{stim} < rsp_t_likertRats(stim,2))) * (1000/abs(rsp_t_likertRats(stim,2) - rsp_t_likertRats(stim,1))); % Hz of this unit at each stimulus at this time window
                        mfrRat2(stim,gg+1) = length(ratLik2{stim}(1,ratLik2{stim} > rsp_t_likertRats(stim,3) & ratLik2{stim} < rsp_t_likertRats(stim,4))) * (1000/abs(rsp_t_likertRats(stim,4) - rsp_t_likertRats(stim,3))); % Hz of this unit at each stimulus at this time window
                        mfrRat1(stim,gg+2) = length(ratLik1{stim}(1,ratLik1{stim} > rsp_t_likertRatl(stim,1) & ratLik1{stim} < rsp_t_likertRatl(stim,2))) * (1000/abs(rsp_t_likertRatl(stim,2) - rsp_t_likertRatl(stim,1))); % Hz of this unit at each stimulus at this time window
                        mfrRat2(stim,gg+2) = length(ratLik2{stim}(1,ratLik2{stim} > rsp_t_likertRatl(stim,3) & ratLik2{stim} < rsp_t_likertRatl(stim,4))) * (1000/abs(rsp_t_likertRatl(stim,4) - rsp_t_likertRatl(stim,3))); % Hz of this unit at each stimulus at this time window
                    end
                    
                end
                
                
                for stim=1:190
                    
                    
                    for gg=1:size(newwinwosreg,1)
                        
                        mfrB1(stim,gg) = length(B1{stim}(1,B1{stim} > newwinwosreg(gg,1) & B1{stim} < newwinwosreg(gg,2))) * (1000/abs(newwinwosreg(gg,2) - newwinwosreg(gg,1))); % Hz of this unit at each stimulus at this time window
                        mfrB2(stim,gg) = length(B2{stim}(1,B2{stim} > newwinwosreg(gg,1) & B2{stim} < newwinwosreg(gg,2))) * (1000/abs(newwinwosreg(gg,2) - newwinwosreg(gg,1))); % Hz of this unit at each stimulus at this time window
                        mfrA21(stim,gg) = length(A21{stim}(1,A21{stim} > newwinwosreg(gg,1) & A21{stim} < newwinwosreg(gg,2))) * (1000/abs(newwinwosreg(gg,2) - newwinwosreg(gg,1))); % Hz of this unit at each stimulus at this time window
                        mfrA22(stim,gg) = length(A22{stim}(1,A22{stim} > newwinwosreg(gg,1) & A22{stim} < newwinwosreg(gg,2))) * (1000/abs(newwinwosreg(gg,2) - newwinwosreg(gg,1))); % Hz of this unit at each stimulus at this time window
                        mfrA11(stim,gg) = length(A11{stim}(1,A11{stim} > newwinwosreg(gg,1) & A11{stim} < newwinwosreg(gg,2))) * (1000/abs(newwinwosreg(gg,2) - newwinwosreg(gg,1))); % Hz of this unit at each stimulus at this time window
                        mfrA12(stim,gg) = length(A12{stim}(1,A12{stim} > newwinwosreg(gg,1) & A12{stim} < newwinwosreg(gg,2))) * (1000/abs(newwinwosreg(gg,2) - newwinwosreg(gg,1))); % Hz of this unit at each stimulus at this time window
                        mfrResp1(stim,gg) = length(likAFCR1{stim}(1,likAFCR1{stim} > newwinwosreg(gg,1) & likAFCR1{stim} < newwinwosreg(gg,2))) * (1000/abs(newwinwosreg(gg,2) - newwinwosreg(gg,1))); % Hz of this unit at each stimulus at this time window
                        mfrResp2(stim,gg) = length(likAFCR2{stim}(1,likAFCR2{stim} > newwinwosreg(gg,1) & likAFCR2{stim} < newwinwosreg(gg,2))) * (1000/abs(newwinwosreg(gg,2) - newwinwosreg(gg,1))); % Hz of this unit at each stimulus at this time window
                        
                    end
                    if bsl_on
                        
                        afcbsl1= length(bsl1{stim}(1,bsl1{stim}> bslwind(1) & bsl1{stim}<bslwind(2)))* (1000/abs(bslwind(2) - bslwind(1)));
                        afcbsl2= length(bsl2{stim}(1,bsl2{stim}> bslwind(1) & bsl2{stim}<bslwind(2)))* (1000/abs(bslwind(2) - bslwind(1)));
                        mfrB1(stim,:) = mfrB1(stim,:) - afcbsl1;
                        mfrA11(stim,:) = mfrA11(stim,:) - afcbsl1;
                        mfrA21(stim,:) = mfrA21(stim,:) - afcbsl1;
                        mfrResp1(stim,:) = mfrResp1(stim,:)- afcbsl1;
                        mfrB2(stim,:) = mfrB2(stim,:) - afcbsl2;
                        mfrA12(stim,:) = mfrA12(stim,:) - afcbsl2;
                        mfrA22(stim,:) = mfrA22(stim,:) - afcbsl2;
                        mfrResp2(stim,:) = mfrResp2(stim,:)- afcbsl2;
                    end
                    
                    
                    if windowtype ==5
                        mfrA11(stim,gg+1) = length(likAFCR1{stim}(1,likAFCR1{stim} > rsp_t_likertAFCs(stim,1) & likAFCR1{stim} < rsp_t_likertAFCs(stim,2))) * (1000/abs(rsp_t_likertAFCs(stim,2) - rsp_t_likertAFCs(stim,1))); % Hz of this unit at each stimulus at this time window
                        mfrA12(stim,gg+1) = length(likAFCR2{stim}(1,likAFCR2{stim} > rsp_t_likertAFCs(stim,3) & likAFCR2{stim} < rsp_t_likertAFCs(stim,4))) * (1000/abs(rsp_t_likertAFCs(stim,4) - rsp_t_likertAFCs(stim,3))); % Hz of this unit at each stimulus at this time window
                        mfrA11(stim,gg+2) = length(likAFCR1{stim}(1,likAFCR1{stim} > rsp_t_likertAFCl(stim,1) & likAFCR1{stim} < rsp_t_likertAFCl(stim,2))) * (1000/abs(rsp_t_likertAFCl(stim,2) - rsp_t_likertAFCl(stim,1))); % Hz of this unit at each stimulus at this time window
                        mfrA12(stim,gg+2) = length(likAFCR2{stim}(1,likAFCR2{stim} > rsp_t_likertAFCl(stim,3) & likAFCR2{stim} < rsp_t_likertAFCl(stim,4))) * (1000/abs(rsp_t_likertAFCl(stim,4) - rsp_t_likertAFCl(stim,3))); % Hz of this unit at each stimulus at this time window
                    end
                    
                end
                
                
                for i= 1:20
                    for hh=1:20
                        if i ~=hh
                            try
                                indPair1 = find(cell2mat(satcue{3}(:,2))==i & cell2mat(satcue{3}(:,3))==hh); indPair2 =  find(cell2mat(satcue{6}(:,2))==i & cell2mat(satcue{6}(:,3))==hh);
                                
                            end
                            try
                                indPair1=find(cell2mat(satcue{3}(:,3))==i & cell2mat(satcue{3}(:,2))==hh);  indPair2 =   find(cell2mat(satcue{6}(:,3))==i & cell2mat(satcue{6}(:,2))==hh);
                            end
                            
                            
                            invmfrA11(indPair1, :) = mfrA12(indPair2,:);
                            invmfrA12 (indPair2, :) = mfrA11 (indPair1,:);
                            
                        end
                    end
                end
                
                
                mfrA11a = mfrA11(2:end,:);
                mfrA12a = mfrA12(2:end,:);
                
                invmfrA11a = invmfrA11(2:end,:);
                invmfrA12a = invmfrA12(2:end,:);
                mfrRat1a=mfrRat1(2:end,:);
                mfrRat2a=mfrRat2(2:end,:);
                
                for i=1:20
                    ind3 = find(cell2mat(satcue{jj-1}(:,2))==i);
                    mfrRat1mean(ind3,:) = mean(mfrRat1(ind3,:));     %CREATES A MEAN MFR FOR EACH STIMULUS ID OVER TIME SORTED BY ID
                    invmfrRat1mean(ind3,:) = mean(mfrRat2(ind3,:));  %CREATES A MEAN MFR FOR EACH STIMULUS ID OVER TIME FOR THE SECOND EXPERIMENTAL SESSION, THIS MAKES AN INVERSE FOR THE COMPARISON 9SINCE SORTED EQUALLY)
                    
                    mfrRatResp1mean(ind3,:) = mean(mfrRatResp1(ind3,:));     %CREATES A MEAN MFR FOR EACH STIMULUS ID OVER TIME SORTED BY ID
                    invmfrRatResp1mean(ind3,:) = mean(mfrRatResp2(ind3,:));  %CREATES A MEAN MFR FOR EACH STIMULUS ID OVER TIME FOR THE SECOND EXPERIMENTAL SESSION, THIS MAKES AN INVERSE FOR THE COMPARISON 9SINCE SORTED EQUALLY)
                    mfrOKRatResp1mean(ind3,:) = mean(mfrOKRatResp1(ind3,:));     %CREATES A MEAN MFR FOR EACH STIMULUS ID OVER TIME SORTED BY ID
                    invmfrOKRatResp1mean(ind3,:) = mean(mfrOKRatResp2(ind3,:));  %CREATES A MEAN MFR FOR EACH STIMULUS ID OVER TIME FOR THE SECOND EXPERIMENTAL SESSION, THIS MAKES AN INVERSE FOR THE COMPARISON 9SINCE SORTED EQUALLY)
                    
                end
                
                if autoregresive_corr==1
                    
                    mfrA11=mfrA11(2:end,:);
                    mfrA12=mfrA12(2:end,:);
                    invmfrA11=invmfrA11(2:end,:);
                    invmfrA12=invmfrA12(2:end,:);
                    mfrResp1=mfrResp1(2:end,:);
                    mfrResp2=mfrResp2(2:end,:);
                    mfrRat1=mfrRat1(2:end,:);
                    mfrRat2=mfrRat2(2:end,:);
                    mfrRatResp1=mfrRatResp1 (2:end,:);
                    mfrRatResp2=mfrRatResp2(2:end,:);
                    mfrOKResp1=mfrOKResp1(2:end,:);
                    mfrOKResp2=mfrOKResp2(2:end,:);
                    mfrOKRatResp1=mfrOKRatResp1(2:end,:);
                    mfrOKRatResp2=mfrOKRatResp2(2:end,:);
                    mfrRat1mean = mfrRat1mean(2:end,:);
                    invmfrRat1mean = invmfrRat1mean(2:end,:);
                    mfrRatResp1mean = mfrRatResp1mean(2:end,:);
                    invmfrRatResp1mean = invmfrRatResp1mean(2:end,:);
                    mfrOKRatResp1mean = mfrOKRatResp1mean(2:end,:);
                    invmfrOKRatResp1mean = invmfrOKRatResp1mean(2:end,:);
                end
                
                
                
                
            elseif windowtype == 10 || windowtype == 30  % gaussian kernel
                
                [mfrRat1] = kernelgauss(ratPic1, windowsize, resolution, stepsize, newwinwosregRating, stdkernel);
                [mfrRat2] = kernelgauss(ratPic2, windowsize, resolution, stepsize, newwinwosregRating, stdkernel);
                [mfrA11] = kernelgauss(A11, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                [mfrA12] = kernelgauss(A12, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                
                %inverted MFR scale output
                for i= 1:20
                    for hh=1:20
                        if i ~=hh
                            try
                                indPair1 = find(cell2mat(satcue{3}(:,2))==i & cell2mat(satcue{3}(:,3))==hh); indPair2 =  find(cell2mat(satcue{6}(:,2))==i & cell2mat(satcue{6}(:,3))==hh);
                                
                            end
                            try
                                indPair1=find(cell2mat(satcue{3}(:,3))==i & cell2mat(satcue{3}(:,2))==hh);  indPair2 =   find(cell2mat(satcue{6}(:,3))==i & cell2mat(satcue{6}(:,2))==hh);
                            end
                            
                            invmfrA11(indPair1, :) = mfrA12(indPair2,:);
                            invmfrA12 (indPair2, :) = mfrA11 (indPair1,:);
                        end
                    end
                end
                
                
                for i=1:20
                    ind3 = find(cell2mat(satcue{jj-1}(:,2))==i);
                    mfrRat1mean(ind3,:) = mean(mfrRat1(ind3,:));     %CREATES A MEAN MFR FOR EACH STIMULUS ID OVER TIME SORTED BY ID
                    invmfrRat1mean(ind3,:) = mean(mfrRat2(ind3,:));  %CREATES A MEAN MFR FOR EACH STIMULUS ID OVER TIME FOR THE SECOND EXPERIMENTAL SESSION, THIS MAKES AN INVERSE FOR THE COMPARISON 9SINCE SORTED EQUALLY)
                    
                    mfrRatResp1mean(ind3,:) = mean(mfrRatResp1(ind3,:));     %CREATES A MEAN MFR FOR EACH STIMULUS ID OVER TIME SORTED BY ID
                    invmfrRatResp1mean(ind3,:) = mean(mfrRatResp2(ind3,:));  %CREATES A MEAN MFR FOR EACH STIMULUS ID OVER TIME FOR THE SECOND EXPERIMENTAL SESSION, THIS MAKES AN INVERSE FOR THE COMPARISON 9SINCE SORTED EQUALLY)
                    mfrOKRatResp1mean(ind3,:) = mean(mfrOKRatResp1(ind3,:));     %CREATES A MEAN MFR FOR EACH STIMULUS ID OVER TIME SORTED BY ID
                    invmfrOKRatResp1mean(ind3,:) = mean(mfrOKRatResp2(ind3,:));  %CREATES A MEAN MFR FOR EACH STIMULUS ID OVER TIME FOR THE SECOND EXPERIMENTAL SESSION, THIS MAKES AN INVERSE FOR THE COMPARISON 9SINCE SORTED EQUALLY)
                    
                end
                
                if autoregresive_corr==1
                    mfrA11=mfrA11(2:end,:);
                    mfrA12=mfrA12(2:end,:);
                    invmfrA11=invmfrA11(2:end,:);
                    invmfrA12=invmfrA12(2:end,:);
                    mfrResp1=mfrResp1(2:end,:);
                    mfrResp2=mfrResp2(2:end,:);
                    mfrRat1=mfrRat1(2:end,:);
                    mfrRat2=mfrRat2(2:end,:);
                    mfrRatResp1=mfrRatResp1 (2:end,:);
                    mfrRatResp2=mfrRatResp2(2:end,:);
                    mfrOKResp1=mfrOKResp1(2:end,:);
                    mfrOKResp2=mfrOKResp2(2:end,:);
                    mfrOKRatResp1=mfrOKRatResp1(2:end,:);
                    mfrOKRatResp2=mfrOKRatResp2(2:end,:);
                    mfrRat1mean = mfrRat1mean(2:end,:);
                    invmfrRat1mean = invmfrRat1mean(2:end,:);
                    mfrRatResp1mean = mfrRatResp1mean(2:end,:);
                    invmfrRatResp1mean = invmfrRatResp1mean(2:end,:);
                    mfrOKRatResp1mean = mfrOKRatResp1mean(2:end,:);
                    invmfrOKRatResp1mean = invmfrOKRatResp1mean(2:end,:);
                end
                
                
                
            end
            
            
            oooo = tic;
            if modeltype == 1 || modeltype == 3 || modeltype == 4
                posnum = [unitnum{sess}(unit,1), unitnum{sess}(unit,1)+ size(cell2mat(unitnum'),1)];
            elseif modeltype == 2
                posnum = [unitnum{sess}(unit,1), unitnum{sess}(unit,1)+ size(cell2mat(unitnum'),1)];
            end
            
            
            if alignedtoDecision == 1
                mfrA11 =  mfrResp1 ;
                mfrA12 =  mfrResp2 ;
                mfrRat1 =  mfrRatResp1 ;
                mfrRat2 =  mfrRatResp2 ;
                mfrRat1mean =  mfrRatResp1mean;
                invmfrRat1mean = invmfrRatResp1mean;
                
            elseif alignedtoDecision == 2
                mfrA11 =  mfrOKResp1 ;
                mfrA12 =  mfrOKResp2 ;
                mfrRat1 =  mfrOKRatResp1;
                mfrRat2 =  mfrOKRatResp2;
                mfrRat1mean =  mfrOKRatResp1mean;
                invmfrRat1mean = invmfrOKRatResp1mean;
            end
            
            
            
            %% Now start all the analysis
            
            
            % % Analysis of all binary regression models
            [pvalRTreg_A1(posnum, :), corrRTreginv_A1(posnum, :), pvalAbsprefreg_A1(posnum, :), corrAbsprefreg_A1(posnum, :)] = bilinearregression ( mfrA11,  mfrA12,   shortRT_AFC,  abs(val_2AFC));
            [pvalRTreginv_A1(posnum, :), corrRTreginv_A1(posnum, :), pvalAbsprefreginv_A1(posnum, :), corrAbsprefreginv_A1(posnum, :)] = bilinearregression ( mfrA11,  mfrA12,   invshortRT_AFC,  abs(invval_2AFC));
            [pvalRTreg_A1(posnum, :), corrRTreg_A1(posnum, :), pvalRTautoreg_A1(posnum, :), corrRTautoreg_A1(posnum, :)] = bilinearregression ( mfrA11(2:end,:),  mfrA12(2:end,:),   shortRT_AFC(2:end,:),  shortRT_AFC(1:end-1,:));
            [pvalRTreginv_A1(posnum, :), corrRTreginv_A1(posnum, :), pvalRTautoreginv_A1(posnum, :), corrRTautoreginv_A1(posnum, :)] = bilinearregression ( mfrA11(2:end,:),  mfrA12(2:end,:),   invshortRT_AFC(2:end,:),  invshortRT_AFC(1:end-1,:));
            [pvalAPreg_A1(posnum, :), corrAPreg_A1(posnum, :), pvalAPautoreg_A1(posnum, :), corrAPautoreg_A1(posnum, :)] = bilinearregression ( mfrA11(2:end,:),  mfrA12(2:end,:),   abs(val_2AFC(2:end,:)),  abs(val_2AFC(1:end-1,:)));
            [pvalAPreginv_A1(posnum, :), corrAPreginv_A1(posnum, :), pvalAPautoreginv_A1(posnum, :), corrAPautoreginv_A1(posnum, :)] = bilinearregression ( mfrA11(2:end,:),  mfrA12(2:end,:),   abs(invval_2AFC(2:end,:)),  abs(invval_2AFC(1:end-1,:)));
            
            %             %% regression AP and Abs difference for reviewer number 5 and chosen value as well
            
            [pvalabsdiffreg_A1(posnum, :), corrabsdiffreg_A1(posnum, :), pvalAbsprefabsdiffreg_A1(posnum, :), corrAbsprefabsdiffreg_A1(posnum, :)] = bilinearregression ( mfrA11,  mfrA12,   abs(val_A_Rat - val_B_Rat),  abs(val_2AFC));
            [pvalabsdiffreginv_A1(posnum, :), corrabsdiffreginv_A1(posnum, :), pvalAbsprefabsdiffreginv_A1(posnum, :), corrAbsprefabsdiffreginv_A1(posnum, :)] = bilinearregression ( mfrA11,  mfrA12,    abs(invval_A_Rat - invval_B_Rat),  abs(invval_2AFC));
            
            [pvalRankabsdiffreg_A1(posnum, :), corrRankabsdiffreg_A1(posnum, :), pvalAbsprefRankabsdiffreg_A1(posnum, :), corrAbsprefRankabsdiffreg_A1(posnum, :)] = bilinearregression ( mfrA11,  mfrA12,   abs(val_A_Rank - val_B_Rank),   abs(val_2AFC));
            [pvalRankabsdiffreginv_A1(posnum, :), corrRankabsdiffreginv_A1(posnum, :), pvalAbsprefRankabsdiffreginv_A1(posnum, :), corrAbsprefRankabsdiffreginv_A1(posnum, :)] = bilinearregression ( mfrA11,  mfrA12,    abs(invval_A_Rank - invval_B_Rank),  abs(invval_2AFC));
            [pvalCvalreg_A1(posnum, :), corrCvalreg_A1(posnum, :), pvalAbsprefCvalreg_A1(posnum, :), corrAbsprefCvalreg_A1(posnum, :)] = bilinearregression ( mfrA11,  mfrA12,   val_C_Rat,  abs(val_2AFC));
            [pvalCvalreginv_A1(posnum, :), corrCvalreginv_A1(posnum, :), pvalAbsprefCvalreginv_A1(posnum, :), corrAbsprefCvalreginv_A1(posnum, :)] = bilinearregression ( mfrA11,  mfrA12,    invval_C_Rat,  abs(invval_2AFC));
            [pvalCval2reginv_A1(posnum, :), corrCval2reginv_A1(posnum, :), pvalAbsprefCval2reginv_A1(posnum, :), corrAbsprefCval2reginv_A1(posnum, :)] = bilinearregression ( mfrA11,  mfrA12,    invval2_C_Rat,  abs(invval_2AFC));
            [pvalCrankreg_A1(posnum, :), corrCrankreg_A1(posnum, :), pvalAbsprefCrankreg_A1(posnum, :), corrAbsprefCrankreg_A1(posnum, :)] = bilinearregression ( mfrA11,  mfrA12,   val_C_Rank,   abs(val_2AFC));
            [pvalCrankreginv_A1(posnum, :), corrCrankreginv_A1(posnum, :), pvalAbsprefCrankreginv_A1(posnum, :), corrAbsprefCrankreginv_A1(posnum, :)] = bilinearregression ( mfrA11,  mfrA12,    invval_C_Rank,  abs(invval_2AFC));
            [pvalCrank2reginv_A1(posnum, :), corrCrank2reginv_A1(posnum, :), pvalAbsprefCrank2reginv_A1(posnum, :), corrAbsprefCrank2reginv_A1(posnum, :)] = bilinearregression ( mfrA11,  mfrA12,    invval2_C_Rank,  abs(invval_2AFC));
            
            
            
            %% Now the Spermann correlation
            %% new variables, valuation/belief confidence/RT, new inverses for chosen and unchosen
            
            [pvalCval2_A1(posnum, :), pvalCrank2_A1(posnum, :), pvalUval2_A1(posnum, :), pvalUrank2_A1(posnum, :), ...
                pvalUCratio2_A1(posnum, :), pvalRankCUratio2_A1(posnum, :), pvalUCratioTot2_A1(posnum, :), ...
                pvalRTmean_R(posnum, :), ...
                corrCval2_A1(posnum, :), corrCrank2_A1(posnum, :), corrUval2_A1(posnum, :), corrUrank2_A1(posnum, :), ...
                corrUCratio2_A1(posnum, :), corrRankCUratio2_A1(posnum, :),  ...
                corrRTmean_R(posnum, :)]= slidereg_funNew(mfrA11, mfrA12,  mfrRat1, mfrRat2, val_C_Rat, val_U_Rat,  val_C_Rank, val_U_Rank,  RTmean_Rat);
            
            
            [pvalCval2inv_A1(posnum, :), pvalCrank2inv_A1(posnum, :), pvalUval2inv_A1(posnum, :), pvalUrank2inv_A1(posnum, :), ...
                pvalUCratio2inv_A1(posnum, :), pvalRankCUratio2inv_A1(posnum, :),  ...
                pvalRTmeaninv_R(posnum, :),...
                corrCval2inv_A1(posnum, :), corrCrank2inv_A1(posnum, :), corrUval2inv_A1(posnum, :), corrUrank2inv_A1(posnum, :),  ...
                corrUCratio2inv_A1(posnum, :), corrRankCUratio2inv_A1(posnum, :),  ...
                corrRTmeaninv_R(posnum, :)]= slidereg_funNew(mfrA11, mfrA12,  mfrRat1, mfrRat2, invval2_C_Rat, invval2_U_Rat,  invval2_C_Rank, invval2_U_Rank, invRTmean_Rat);
            
            
            
            %% change in value  correlations and productwise
            
            [pvalAval_A1(posnum, :), pvalBval_A1(posnum, :),  pvalAsal_A1(posnum, :), pvalBsal_A1(posnum, :), pvalCval_A1(posnum, :), pvalUval_A1(posnum, :), pvalabsdiff_A1(posnum, :),...
                pvaldiffAB_A1(posnum, :), pvalAbspref_A1(posnum, :), pvalprefAB_A1(posnum, :),  pvalSum_A1(posnum, :),pvalUCratio_A1(posnum, :),  ...
                pvalArank_A1(posnum, :),pvalBrank_A1(posnum, :), pvalCrank_A1(posnum, :), ...
                pvalUrank_A1(posnum, :), pvalRankabsdiff_A1(posnum, :), pvalRanksum_A1(posnum, :), pvalRankABdiff_A1(posnum, :),pvalRankCUratio_A1(posnum, :), ...
                corrAval_A1(posnum, :), corrBval_A1(posnum, :),  corrAsal_A1(posnum, :), corrBsal_A1(posnum, :),corrCval_A1(posnum, :),  corrUval_A1(posnum, :),   corrabsdiff_A1(posnum, :), ...
                corrdiffAB_A1(posnum, :), corrAbspref_A1(posnum, :),corrprefAB_A1(posnum, :), corrSum_A1(posnum, :),  corrUCratio_A1(posnum, :), ...
                corrArank_A1(posnum, :),corrBrank_A1(posnum, :), corrCrank_A1(posnum, :), ...
                corrUrank_A1(posnum, :), corrRankabsdiff_A1(posnum, :), corrRanksum_A1(posnum, :), corrRankABdiff_A1(posnum, :),corrRankCUratio_A1(posnum, :), ...
                pvalRT_A1(posnum, :), corrRT_A1(posnum, :), ...
                pvalflavourA_A1(posnum, :), pvalflavourB_A1(posnum, :), pvalflavourC_A1(posnum, :), pvalflavourU_A1(posnum, :) ,  pvalABbin_A1(posnum, :), ...
                pvalSignA_A1(posnum, :), pvalSignB_A1(posnum, :), pvalSignC_A1(posnum, :), pvalSignU_A1(posnum, :), ...
                corrflavourA_A1(posnum, :), corrflavourB_A1(posnum, :),  corrflavourC_A1(posnum, :), corrflavourU_A1(posnum, :),  corrABbin_A1(posnum, :), corrtrial_A1(posnum, :),  pvaltrial_A1(posnum, :)]=  slidereg_fun(mfrA11, mfrA12, ...
                val_A_Rat,  val_B_Rat,   val_C_Rat,    val_U_Rat, val_A_Rank, val_B_Rank, val_C_Rank, val_U_Rank, val_diffAB, shortRT_AFC, flavourA, flavourB,  flavourC, flavourU);
            
            
            
            
            [pvalLval_A1(posnum, :), pvalRval_A1(posnum, :),  pvalLsal_A1(posnum, :), pvalRsal_A1(posnum, :), ~,~, ~,...
                pvaldiffLR_A1(posnum, :), ~, pvalprefLR_A1(posnum, :),  ~,~,  ...
                pvalLrank_A1(posnum, :),pvalRrank_A1(posnum, :), ~, ...
                ~,~, ~, pvalRankLRdiff_A1(posnum, :), ~, ...
                corrLval_A1(posnum, :), corrRval_A1(posnum, :), corrLsal_A1(posnum, :), corrRsal_A1(posnum, :), ~,  ~, ~, ...
                corrdiffLR_A1(posnum, :), ~, corrprefLR_A1(posnum, :), ~, ~, ...
                corrLrank_A1(posnum, :), corrRrank_A1(posnum, :), ~, ...
                ~,~,~, corrRankLRdiff_A1(posnum, :), ~, ...
                ~, ~, ...
                pvalflavourL_A1(posnum, :), pvalflavourR_A1(posnum, :), ~,  ~,  pvalLRbin_A1(posnum, :), ...
                pvalSignL_A1(posnum, :), pvalSignR_A1(posnum, :), ~,  ~, ...
                corrflavourL_A1(posnum, :), corrflavourR_A1(posnum, :),  ~,  ~, corrLRbin_A1(posnum, :), ~,~]=  slidereg_fun(mfrA11, mfrA12, ...
                val_L_Rat,  val_R_Rat,   val_C_Rat,    val_U_Rat, val_L_Rank, val_R_Rank, val_C_Rank, val_U_Rank, val_2AFC, shortRT_AFC, flavourL, flavourR,  flavourC, flavourU);
            
            
            [pvalVal_R(posnum, :), corrVal_R(posnum, :),  pvalSal_R(posnum, :),  corrSal_R(posnum, :), pvalValmean_R(posnum, :),corrValmean_R(posnum, :),pvalSalmean_R(posnum, :),corrSalmean_R(posnum, :),...
                pvalflavourRat_R(posnum, :), corrflavourRat_R(posnum, :), ...
                pvalSign_R(posnum, :), corrSign_R(posnum, :), ...
                pvalRank_R(posnum, :), corrRank_R(posnum, :), ...
                pvalRT_R(posnum, :), corrRT_R(posnum, :)] = regRating_fun( mfrRat1, mfrRat2, val_rat,  val_rat_Rat, flavourRat, val_rat_Rank,  shortRT_Rat, stimOrderRat);
            
            
            
            
            % ineverses for analysis of the satiation break
            [pvalAvalinv_A1(posnum, :), pvalBvalinv_A1(posnum, :),  pvalAsalinv_A1(posnum, :), pvalBsalinv_A1(posnum, :), pvalCvalinv_A1(posnum, :), pvalUvalinv_A1(posnum, :), pvalabsdiffinv_A1(posnum, :),...
                pvaldiffABinv_A1(posnum, :), pvalAbsprefinv_A1(posnum, :), pvalprefABinv_A1(posnum, :),  pvalSuminv_A1(posnum, :),pvalUCratioinv_A1(posnum, :),  ...
                pvalArankinv_A1(posnum, :),pvalBrankinv_A1(posnum, :), pvalCrankinv_A1(posnum, :), ...
                pvalUrankinv_A1(posnum, :), pvalRankabsdiffinv_A1(posnum, :), pvalRanksuminv_A1(posnum, :), pvalRankABdiffinv_A1(posnum, :),pvalRankCUratioinv_A1(posnum, :), ...
                corrAvalinv_A1(posnum, :), corrBvalinv_A1(posnum, :),  corrAsalinv_A1(posnum, :), corrBsalinv_A1(posnum, :),corrCvalinv_A1(posnum, :),  corrUvalinv_A1(posnum, :),   corrabsdiffinv_A1(posnum, :), ...
                corrdiffABinv_A1(posnum, :), corrAbsprefinv_A1(posnum, :),corrprefABinv_A1(posnum, :), corrSuminv_A1(posnum, :),  corrUCratioinv_A1(posnum, :), ...
                corrArankinv_A1(posnum, :),corrBrankinv_A1(posnum, :), corrCrankinv_A1(posnum, :), ...
                corrUrankinv_A1(posnum, :), corrRankabsdiffinv_A1(posnum, :), corrRanksuminv_A1(posnum, :), corrRankABdiffinv_A1(posnum, :),corrRankCUratioinv_A1(posnum, :), ...
                pvalRTinv_A1(posnum, :), corrRTinv_A1(posnum, :),...
                pvalflavourAinv_A1(posnum, :), pvalflavourBinv_A1(posnum, :), pvalflavourCinv_A1(posnum, :), pvalflavourUinv_A1(posnum, :) ,  pvalABbininv_A1(posnum, :), corrtrialinv_A1(posnum, :), pvaltrialinv_A1(posnum, :)]=  slidereg_fun(mfrA11, mfrA12, ...
                invval_A_Rat,  invval_B_Rat,   invval_C_Rat,    invval_U_Rat, invval_A_Rank, invval_B_Rank, invval_C_Rank, invval_U_Rank, invval_diffAB, invshortRT_AFC, flavourA, flavourB,  flavourC, flavourU);
            
            
            
            [pvalLvalinv_A1(posnum, :), pvalRvalinv_A1(posnum, :),  pvalLsalinv_A1(posnum, :), pvalRsalinv_A1(posnum, :), ~,~, ~,...
                pvaldiffLRinv_A1(posnum, :), ~, pvalprefLRinv_A1(posnum, :),  ~,~,  ...
                pvalLrankinv_A1(posnum, :),pvalRrankinv_A1(posnum, :), ~, ...
                ~,~, ~, pvalRankLRdiffinv_A1(posnum, :), ~, ...
                corrLvalinv_A1(posnum, :), corrRvalinv_A1(posnum, :), corrLsalinv_A1(posnum, :), corrRsalinv_A1(posnum, :), ~,  ~, ~, ...
                corrdiffLRinv_A1(posnum, :), ~, corrprefLRinv_A1(posnum, :), ~, ~, ...
                corrLrankinv_A1(posnum, :), corrRrankinv_A1(posnum, :), ~, ...
                ~,~,~, corrRankLRdiffinv_A1(posnum, :), ~, ...
                ~, ~, ...
                pvalflavourLinv_A1(posnum, :), pvalflavourRinv_A1(posnum, :), ~,  ~,  pvalLRbininv_A1(posnum, :), ...
                pvalSignLinv_A1(posnum, :), pvalSignRinv_A1(posnum, :), ~,  ~, ...
                corrflavourLinv_A1(posnum, :), corrflavourRinv_A1(posnum, :),  ~,  ~, corrLRbininv_A1(posnum, :), ~,~]=  slidereg_fun(mfrA11, mfrA12, ...
                invval_L_Rat,  invval_R_Rat,   invval_C_Rat,    invval_U_Rat, invval_L_Rank, invval_R_Rank, invval_C_Rank, invval_U_Rank, invval_2AFC, invshortRT_AFC, flavourL, flavourR,  flavourC, flavourU);
            
            [pvalValinv_R(posnum, :), corrValinv_R(posnum, :),  pvalSalinv_R(posnum, :),  corrSalinv_R(posnum, :), pvalValmeaninv_R(posnum, :),corrValmeaninv_R(posnum, :),pvalSalmeaninv_R(posnum, :),corrSalmeaninv_R(posnum, :),...
                ~, ~, ...
                pvalSigninv_R(posnum, :), corrSigninv_R(posnum, :), ...
                pvalRankinv_R(posnum, :), corrRankinv_R(posnum, :), ...
                ~,~] = regRating_fun( mfrRat1, mfrRat2, invval_rat_Rat,  invval_rat_Rat, flavourRat, invval_rat_Rank,  shortRT_Rat, stimOrderRat);
            
            
            
            
            %% Analyses visual features
            
            [pvalARed_A1(posnum, :),pvalAGreen_A1(posnum, :), pvalABlue_A1(posnum, :), pvalAContrast_A1(posnum, :), pvalALuminance_A1(posnum, :), pvalAHue_A1(posnum, :), pvalASaturation_A1(posnum, :),...
                corrARed_A1(posnum, :),corrAGreen_A1(posnum, :), corrABlue_A1(posnum, :), corrAContrast_A1(posnum, :), corrALuminance_A1(posnum, :), corrAHue_A1(posnum, :), corrASaturation_A1(posnum, :),...
                pvalARed_Productwise_A1(posnum, :),pvalAGreen_Productwise_A1(posnum, :), pvalABlue_Productwise_A1(posnum, :), pvalAContrast_Productwise_A1(posnum, :), pvalALuminance_Productwise_A1(posnum, :), pvalAHue_Productwise_A1(posnum, :), pvalASaturation_Productwise_A1(posnum, :),...
                corrARed_Productwise_A1(posnum, :),corrAGreen_Productwise_A1(posnum, :), corrABlue_Productwise_A1(posnum, :), corrAContrast_Productwise_A1(posnum, :), corrALuminance_Productwise_A1(posnum, :), corrAHue_Productwise_A1(posnum, :), corrASaturation_Productwise_A1(posnum, :),...
                ] = corrVisual(mfrA11, mfrA12,  Red, Green, Blue, Contrast, Luminance, Hue, Saturation, stimOrderA, Red_A, Green_A, Blue_A, Contrast_A, Luminance_A, Hue_A, Saturation_A);
            
            
            [pvalBRed_A1(posnum, :),pvalBGreen_A1(posnum, :), pvalBBlue_A1(posnum, :), pvalBContrast_A1(posnum, :), pvalBLuminance_A1(posnum, :), pvalBHue_A1(posnum, :), pvalBSaturation_A1(posnum, :),...
                corrBRed_A1(posnum, :),corrBGreen_A1(posnum, :), corrBBlue_A1(posnum, :), corrBContrast_A1(posnum, :), corrBLuminance_A1(posnum, :), corrBHue_A1(posnum, :), corrBSaturation_A1(posnum, :),...
                pvalBRed_Productwise_A1(posnum, :),pvalBGreen_Productwise_A1(posnum, :), pvalBBlue_Productwise_A1(posnum, :), pvalBContrast_Productwise_A1(posnum, :), pvalBLuminance_Productwise_A1(posnum, :), pvalBHue_Productwise_A1(posnum, :), pvalBSaturation_Productwise_A1(posnum, :),...
                corrBRed_Productwise_A1(posnum, :),corrBGreen_Productwise_A1(posnum, :), corrBBlue_Productwise_A1(posnum, :), corrBContrast_Productwise_A1(posnum, :), corrBLuminance_Productwise_A1(posnum, :), corrBHue_Productwise_A1(posnum, :), corrBSaturation_Productwise_A1(posnum, :),...
                ] = corrVisual(mfrA11, mfrA12,  Red, Green, Blue, Contrast, Luminance, Hue, Saturation, stimOrderB, Red_B, Green_B, Blue_B, Contrast_B, Luminance_B, Hue_B, Saturation_B);
            
            
            
            [pvalRed_R(posnum, :),pvalGreen_R(posnum, :), pvalBlue_R(posnum, :), pvalContrast_R(posnum, :), pvalLuminance_R(posnum, :), pvalHue_R(posnum, :), pvalSaturation_R(posnum, :),...
                corrRed_R(posnum, :),corrGreen_R(posnum, :), corrBlue_R(posnum, :), corrContrast_R(posnum, :), corrLuminance_R(posnum, :), corrHue_R(posnum, :), corrSaturation_R(posnum, :),...
                pvalRed_Productwise_R(posnum, :),pvalGreen_Productwise_R(posnum, :), pvalBlue_Productwise_R(posnum, :), pvalContrast_Productwise_R(posnum, :), pvalLuminance_Productwise_R(posnum, :), pvalHue_Productwise_R(posnum, :), pvalSaturation_Productwise_R(posnum, :),...
                corrRed_Productwise_R(posnum, :),corrGreen_Productwise_R(posnum, :), corrBlue_Productwise_R(posnum, :), corrContrast_Productwise_R(posnum, :), corrLuminance_Productwise_R(posnum, :), corrHue_Productwise_R(posnum, :), corrSaturation_Productwise_R(posnum, :),...
                ] = corrVisual(mfrRat1, mfrRat2,  Red, Green, Blue, Contrast, Luminance, Hue, Saturation, stimOrderRat, Red_Rat, Green_Rat, Blue_Rat, Contrast_Rat, Luminance_Rat, Hue_Rat, Saturation_Rat);
            
            
            
            
            %% analyses nutritional features
            
            
            [pvalAKalorien_A1(posnum, :),pvalAProtein_A1(posnum, :), pvalAKohlenhydrate_A1(posnum, :), pvalAZucker_A1(posnum, :), pvalAFett_A1(posnum, :), pvalABallaststoffe_A1(posnum, :), pvalASalz_A1(posnum, :),...
                corrAKalorien_A1(posnum, :),corrAProtein_A1(posnum, :), corrAKohlenhydrate_A1(posnum, :), corrAZucker_A1(posnum, :), corrAFett_A1(posnum, :), corrABallaststoffe_A1(posnum, :), corrASalz_A1(posnum, :),...
                pvalAKalorien_Productwise_A1(posnum, :),pvalAProtein_Productwise_A1(posnum, :), pvalAKohlenhydrate_Productwise_A1(posnum, :), pvalAZucker_Productwise_A1(posnum, :), pvalAFett_Productwise_A1(posnum, :), pvalABallaststoffe_Productwise_A1(posnum, :), pvalASalz_Productwise_A1(posnum, :),...
                corrAKalorien_Productwise_A1(posnum, :),corrAProtein_Productwise_A1(posnum, :), corrAKohlenhydrate_Productwise_A1(posnum, :), corrAZucker_Productwise_A1(posnum, :), corrAFett_Productwise_A1(posnum, :), corrABallaststoffe_Productwise_A1(posnum, :), corrASalz_Productwise_A1(posnum, :),...
                ] = corrVisual(mfrA11, mfrA12,  Kalorien, Protein, Kohlenhydrate, Zucker, Fett, Ballaststoffe, Salz, stimOrderA, Kalorien_A, Protein_A, Kohlenhydrate_A, Zucker_A, Fett_A, Ballaststoffe_A, Salz_A);
            
            
            [pvalBKalorien_A1(posnum, :),pvalBProtein_A1(posnum, :), pvalBKohlenhydrate_A1(posnum, :), pvalBZucker_A1(posnum, :), pvalBFett_A1(posnum, :), pvalBBallaststoffe_A1(posnum, :), pvalBSalz_A1(posnum, :),...
                corrBKalorien_A1(posnum, :),corrBProtein_A1(posnum, :), corrBKohlenhydrate_A1(posnum, :), corrBZucker_A1(posnum, :), corrBFett_A1(posnum, :), corrBBallaststoffe_A1(posnum, :), corrBSalz_A1(posnum, :),...
                pvalBKalorien_Productwise_A1(posnum, :),pvalBProtein_Productwise_A1(posnum, :), pvalBKohlenhydrate_Productwise_A1(posnum, :), pvalBZucker_Productwise_A1(posnum, :), pvalBFett_Productwise_A1(posnum, :), pvalBBallaststoffe_Productwise_A1(posnum, :), pvalBSalz_Productwise_A1(posnum, :),...
                corrBKalorien_Productwise_A1(posnum, :),corrBProtein_Productwise_A1(posnum, :), corrBKohlenhydrate_Productwise_A1(posnum, :), corrBZucker_Productwise_A1(posnum, :), corrBFett_Productwise_A1(posnum, :), corrBBallaststoffe_Productwise_A1(posnum, :), corrBSalz_Productwise_A1(posnum, :),...
                ] = corrVisual(mfrA11, mfrA12,  Kalorien, Protein, Kohlenhydrate, Zucker, Fett, Ballaststoffe, Salz, stimOrderB, Kalorien_B, Protein_B, Kohlenhydrate_B, Zucker_B, Fett_B, Ballaststoffe_B, Salz_B);
            
            
            
            [pvalKalorien_R(posnum, :),pvalProtein_R(posnum, :), pvalKohlenhydrate_R(posnum, :), pvalZucker_R(posnum, :), pvalFett_R(posnum, :), pvalBallaststoffe_R(posnum, :), pvalSalz_R(posnum, :),...
                corrKalorien_R(posnum, :),corrProtein_R(posnum, :), corrKohlenhydrate_R(posnum, :), corrZucker_R(posnum, :), corrFett_R(posnum, :), corrBallaststoffe_R(posnum, :), corrSalz_R(posnum, :),...
                pvalKalorien_Productwise_R(posnum, :),pvalProtein_Productwise_R(posnum, :), pvalKohlenhydrate_Productwise_R(posnum, :), pvalZucker_Productwise_R(posnum, :), pvalFett_Productwise_R(posnum, :), pvalBallaststoffe_Productwise_R(posnum, :), pvalSalz_Productwise_R(posnum, :),...
                corrKalorien_Productwise_R(posnum, :),corrProtein_Productwise_R(posnum, :), corrKohlenhydrate_Productwise_R(posnum, :), corrZucker_Productwise_R(posnum, :), corrFett_Productwise_R(posnum, :), corrBallaststoffe_Productwise_R(posnum, :), corrSalz_Productwise_R(posnum, :),...
                ] = corrVisual(mfrRat1, mfrRat2,  Kalorien, Protein, Kohlenhydrate, Zucker, Fett, Ballaststoffe, Salz, stimOrderRat, Kalorien_Rat, Protein_Rat, Kohlenhydrate_Rat, Zucker_Rat, Fett_Rat, Ballaststoffe_Rat, Salz_Rat);
            
            
            
            
            
            %% just LR and CU for visual and nutritional
            
            [pvalLRed_A1(posnum, :),pvalLGreen_A1(posnum, :), pvalLBlue_A1(posnum, :), pvalLContrast_A1(posnum, :), pvalLLuminance_A1(posnum, :), pvalLHue_A1(posnum, :), pvalLSaturation_A1(posnum, :),...
                corrLRed_A1(posnum, :),corrLGreen_A1(posnum, :), corrLBlue_A1(posnum, :), corrLContrast_A1(posnum, :), corrLLuminance_A1(posnum, :), corrLHue_A1(posnum, :), corrLSaturation_A1(posnum, :),...
                pvalLRed_Productwise_A1(posnum, :),pvalLGreen_Productwise_A1(posnum, :), pvalLBlue_Productwise_A1(posnum, :), pvalLContrast_Productwise_A1(posnum, :), pvalLLuminance_Productwise_A1(posnum, :), pvalLHue_Productwise_A1(posnum, :), pvalLSaturation_Productwise_A1(posnum, :),...
                corrLRed_Productwise_A1(posnum, :),corrLGreen_Productwise_A1(posnum, :), corrLBlue_Productwise_A1(posnum, :), corrLContrast_Productwise_A1(posnum, :), corrLLuminance_Productwise_A1(posnum, :), corrLHue_Productwise_A1(posnum, :), corrLSaturation_Productwise_A1(posnum, :),...
                ] = corrVisual(mfrA11, mfrA12,  Red, Green, Blue, Contrast, Luminance, Hue, Saturation, stimOrderL, Red_L, Green_L, Blue_L, Contrast_L, Luminance_L, Hue_L, Saturation_L);
            
            
            [pvalRRed_A1(posnum, :),pvalRGreen_A1(posnum, :), pvalRBlue_A1(posnum, :), pvalRContrast_A1(posnum, :), pvalRLuminance_A1(posnum, :), pvalRHue_A1(posnum, :), pvalRSaturation_A1(posnum, :),...
                corrRRed_A1(posnum, :),corrRGreen_A1(posnum, :), corrRBlue_A1(posnum, :), corrRContrast_A1(posnum, :), corrRLuminance_A1(posnum, :), corrRHue_A1(posnum, :), corrRSaturation_A1(posnum, :),...
                pvalRRed_Productwise_A1(posnum, :),pvalRGreen_Productwise_A1(posnum, :), pvalRBlue_Productwise_A1(posnum, :), pvalRContrast_Productwise_A1(posnum, :), pvalRLuminance_Productwise_A1(posnum, :), pvalRHue_Productwise_A1(posnum, :), pvalRSaturation_Productwise_A1(posnum, :),...
                corrRRed_Productwise_A1(posnum, :),corrRGreen_Productwise_A1(posnum, :), corrRBlue_Productwise_A1(posnum, :), corrRContrast_Productwise_A1(posnum, :), corrRLuminance_Productwise_A1(posnum, :), corrRHue_Productwise_A1(posnum, :), corrRSaturation_Productwise_A1(posnum, :),...
                ] = corrVisual(mfrA11, mfrA12,  Red, Green, Blue, Contrast, Luminance, Hue, Saturation, stimOrderR, Red_R, Green_R, Blue_R, Contrast_R, Luminance_R, Hue_R, Saturation_R);
            
            [pvalLKalorien_A1(posnum, :),pvalLProtein_A1(posnum, :), pvalLKohlenhydrate_A1(posnum, :), pvalLZucker_A1(posnum, :), pvalLFett_A1(posnum, :), pvalLBallaststoffe_A1(posnum, :), pvalLSalz_A1(posnum, :),...
                corrLKalorien_A1(posnum, :),corrLProtein_A1(posnum, :), corrLKohlenhydrate_A1(posnum, :), corrLZucker_A1(posnum, :), corrLFett_A1(posnum, :), corrLBallaststoffe_A1(posnum, :), corrLSalz_A1(posnum, :),...
                pvalLKalorien_Productwise_A1(posnum, :),pvalLProtein_Productwise_A1(posnum, :), pvalLKohlenhydrate_Productwise_A1(posnum, :), pvalLZucker_Productwise_A1(posnum, :), pvalLFett_Productwise_A1(posnum, :), pvalLBallaststoffe_Productwise_A1(posnum, :), pvalLSalz_Productwise_A1(posnum, :),...
                corrLKalorien_Productwise_A1(posnum, :),corrLProtein_Productwise_A1(posnum, :), corrLKohlenhydrate_Productwise_A1(posnum, :), corrLZucker_Productwise_A1(posnum, :), corrLFett_Productwise_A1(posnum, :), corrLBallaststoffe_Productwise_A1(posnum, :), corrLSalz_Productwise_A1(posnum, :),...
                ] = corrVisual(mfrA11, mfrA12,  Kalorien, Protein, Kohlenhydrate, Zucker, Fett, Ballaststoffe, Salz, stimOrderL, Kalorien_L, Protein_L, Kohlenhydrate_L, Zucker_L, Fett_L, Ballaststoffe_L, Salz_L);
            
            
            [pvalRKalorien_A1(posnum, :),pvalRProtein_A1(posnum, :), pvalRKohlenhydrate_A1(posnum, :), pvalRZucker_A1(posnum, :), pvalRFett_A1(posnum, :), pvalRBallaststoffe_A1(posnum, :), pvalRSalz_A1(posnum, :),...
                corrRKalorien_A1(posnum, :),corrRProtein_A1(posnum, :), corrRKohlenhydrate_A1(posnum, :), corrRZucker_A1(posnum, :), corrRFett_A1(posnum, :), corrRBallaststoffe_A1(posnum, :), corrRSalz_A1(posnum, :),...
                pvalRKalorien_Productwise_A1(posnum, :),pvalRProtein_Productwise_A1(posnum, :), pvalRKohlenhydrate_Productwise_A1(posnum, :), pvalRZucker_Productwise_A1(posnum, :), pvalRFett_Productwise_A1(posnum, :), pvalRBallaststoffe_Productwise_A1(posnum, :), pvalRSalz_Productwise_A1(posnum, :),...
                corrRKalorien_Productwise_A1(posnum, :),corrRProtein_Productwise_A1(posnum, :), corrRKohlenhydrate_Productwise_A1(posnum, :), corrRZucker_Productwise_A1(posnum, :), corrRFett_Productwise_A1(posnum, :), corrRBallaststoffe_Productwise_A1(posnum, :), corrRSalz_Productwise_A1(posnum, :),...
                ] = corrVisual(mfrA11, mfrA12,  Kalorien, Protein, Kohlenhydrate, Zucker, Fett, Ballaststoffe, Salz, stimOrderR, Kalorien_R, Protein_R, Kohlenhydrate_R, Zucker_R, Fett_R, Ballaststoffe_R, Salz_R);
            
            
            
            [pvalCRed_A1(posnum, :),pvalCGreen_A1(posnum, :), pvalCBlue_A1(posnum, :), pvalCContrast_A1(posnum, :), pvalCLuminance_A1(posnum, :), pvalCHue_A1(posnum, :), pvalCSaturation_A1(posnum, :),...
                corrCRed_A1(posnum, :),corrCGreen_A1(posnum, :), corrCBlue_A1(posnum, :), corrCContrast_A1(posnum, :), corrCLuminance_A1(posnum, :), corrCHue_A1(posnum, :), corrCSaturation_A1(posnum, :),...
                pvalCRed_Productwise_A1(posnum, :),pvalCGreen_Productwise_A1(posnum, :), pvalCBlue_Productwise_A1(posnum, :), pvalCContrast_Productwise_A1(posnum, :), pvalCLuminance_Productwise_A1(posnum, :), pvalCHue_Productwise_A1(posnum, :), pvalCSaturation_Productwise_A1(posnum, :),...
                corrCRed_Productwise_A1(posnum, :),corrCGreen_Productwise_A1(posnum, :), corrCBlue_Productwise_A1(posnum, :), corrCContrast_Productwise_A1(posnum, :), corrCLuminance_Productwise_A1(posnum, :), corrCHue_Productwise_A1(posnum, :), corrCSaturation_Productwise_A1(posnum, :),...
                ] = corrVisual(mfrA11, mfrA12,  Red, Green, Blue, Contrast, Luminance, Hue, Saturation, stimOrderC, Red_C, Green_C, Blue_C, Contrast_C, Luminance_C, Hue_C, Saturation_C);
            
            
            [pvalURed_A1(posnum, :),pvalUGreen_A1(posnum, :), pvalUBlue_A1(posnum, :), pvalUContrast_A1(posnum, :), pvalULuminance_A1(posnum, :), pvalUHue_A1(posnum, :), pvalUSaturation_A1(posnum, :),...
                corrURed_A1(posnum, :),corrUGreen_A1(posnum, :), corrUBlue_A1(posnum, :), corrUContrast_A1(posnum, :), corrULuminance_A1(posnum, :), corrUHue_A1(posnum, :), corrUSaturation_A1(posnum, :),...
                pvalURed_Productwise_A1(posnum, :),pvalUGreen_Productwise_A1(posnum, :), pvalUBlue_Productwise_A1(posnum, :), pvalUContrast_Productwise_A1(posnum, :), pvalULuminance_Productwise_A1(posnum, :), pvalUHue_Productwise_A1(posnum, :), pvalUSaturation_Productwise_A1(posnum, :),...
                corrURed_Productwise_A1(posnum, :),corrUGreen_Productwise_A1(posnum, :), corrUBlue_Productwise_A1(posnum, :), corrUContrast_Productwise_A1(posnum, :), corrULuminance_Productwise_A1(posnum, :), corrUHue_Productwise_A1(posnum, :), corrUSaturation_Productwise_A1(posnum, :),...
                ] = corrVisual(mfrA11, mfrA12,  Red, Green, Blue, Contrast, Luminance, Hue, Saturation, stimOrderU, Red_U, Green_U, Blue_U, Contrast_U, Luminance_U, Hue_U, Saturation_U);
            
            [pvalCKalorien_A1(posnum, :),pvalCProtein_A1(posnum, :), pvalCKohlenhydrate_A1(posnum, :), pvalCZucker_A1(posnum, :), pvalCFett_A1(posnum, :), pvalCBallaststoffe_A1(posnum, :), pvalCSalz_A1(posnum, :),...
                corrCKalorien_A1(posnum, :),corrCProtein_A1(posnum, :), corrCKohlenhydrate_A1(posnum, :), corrCZucker_A1(posnum, :), corrCFett_A1(posnum, :), corrCBallaststoffe_A1(posnum, :), corrCSalz_A1(posnum, :),...
                pvalCKalorien_Productwise_A1(posnum, :),pvalCProtein_Productwise_A1(posnum, :), pvalCKohlenhydrate_Productwise_A1(posnum, :), pvalCZucker_Productwise_A1(posnum, :), pvalCFett_Productwise_A1(posnum, :), pvalCBallaststoffe_Productwise_A1(posnum, :), pvalCSalz_Productwise_A1(posnum, :),...
                corrCKalorien_Productwise_A1(posnum, :),corrCProtein_Productwise_A1(posnum, :), corrCKohlenhydrate_Productwise_A1(posnum, :), corrCZucker_Productwise_A1(posnum, :), corrCFett_Productwise_A1(posnum, :), corrCBallaststoffe_Productwise_A1(posnum, :), corrCSalz_Productwise_A1(posnum, :),...
                ] = corrVisual(mfrA11, mfrA12,  Kalorien, Protein, Kohlenhydrate, Zucker, Fett, Ballaststoffe, Salz, stimOrderC, Kalorien_C, Protein_C, Kohlenhydrate_C, Zucker_C, Fett_C, Ballaststoffe_C, Salz_C);
            
            
            [pvalUKalorien_A1(posnum, :),pvalUProtein_A1(posnum, :), pvalUKohlenhydrate_A1(posnum, :), pvalUZucker_A1(posnum, :), pvalUFett_A1(posnum, :), pvalUBallaststoffe_A1(posnum, :), pvalUSalz_A1(posnum, :),...
                corrUKalorien_A1(posnum, :),corrUProtein_A1(posnum, :), corrUKohlenhydrate_A1(posnum, :), corrUZucker_A1(posnum, :), corrUFett_A1(posnum, :), corrUBallaststoffe_A1(posnum, :), corrUSalz_A1(posnum, :),...
                pvalUKalorien_Productwise_A1(posnum, :),pvalUProtein_Productwise_A1(posnum, :), pvalUKohlenhydrate_Productwise_A1(posnum, :), pvalUZucker_Productwise_A1(posnum, :), pvalUFett_Productwise_A1(posnum, :), pvalUBallaststoffe_Productwise_A1(posnum, :), pvalUSalz_Productwise_A1(posnum, :),...
                corrUKalorien_Productwise_A1(posnum, :),corrUProtein_Productwise_A1(posnum, :), corrUKohlenhydrate_Productwise_A1(posnum, :), corrUZucker_Productwise_A1(posnum, :), corrUFett_Productwise_A1(posnum, :), corrUBallaststoffe_Productwise_A1(posnum, :), corrUSalz_Productwise_A1(posnum, :),...
                ] = corrVisual(mfrA11, mfrA12,  Kalorien, Protein, Kohlenhydrate, Zucker, Fett, Ballaststoffe, Salz, stimOrderU, Kalorien_U, Protein_U, Kohlenhydrate_U, Zucker_U, Fett_U, Ballaststoffe_U, Salz_U);
            
            
            
            %% direct change in mfr: correlation of the change of mfr with    the cahnge in the variable
            [pvalAbspref_RawChange_A1(posnum(1), :), pvalprefLR_RawChange_A1(posnum(1), :), pvalprefAB_RawChange_A1(posnum(1), :), pvalRT_RawChange_A1(posnum(1), :),  pvalABbin_RawChange_A1(posnum(1), :), ...
                pvalAval_RawChange_A1(posnum(1), :), pvalBval_RawChange_A1(posnum(1), :), pvalCval_RawChange_A1(posnum(1), :), pvalUval_RawChange_A1(posnum(1), :),  ...
                pvalUCratio_RawChange_A1(posnum(1), :), pvaldiffAB_RawChange_A1(posnum(1), :), pvalabsdiff_RawChange_A1(posnum(1), :), ...
                pvalSum_RawChange_A1(posnum(1), :), pvalAsal_RawChange_A1(posnum(1), :), pvalBsal_RawChange_A1(posnum(1), :), ...
                pvalArank_RawChange_A1(posnum(1), :), pvalBrank_RawChange_A1(posnum(1), :), pvalCrank_RawChange_A1(posnum(1), :), pvalUrank_RawChange_A1(posnum(1), :), ...
                pvalRankCUratio_RawChange_A1(posnum(1), :), pvalRankABdiff_RawChange_A1(posnum(1), :),  pvalRankabsdiff_RawChange_A1(posnum(1), :), ...
                pvalRanksum_RawChange_A1(posnum(1), :), ...
                corrAbspref_RawChange_A1(posnum(1), :), corrprefLR_RawChange_A1(posnum(1), :), corrprefAB_RawChange_A1(posnum(1), :), corrRT_RawChange_A1(posnum(1), :),  corrABbin_RawChange_A1(posnum(1), :), ...
                corrAval_RawChange_A1(posnum(1), :), corrBval_RawChange_A1(posnum(1), :), corrCval_RawChange_A1(posnum(1), :), corrUval_RawChange_A1(posnum(1), :), ...
                corrUCratio_RawChange_A1(posnum(1), :), corrdiffAB_RawChange_A1(posnum(1), :), corrabsdiff_RawChange_A1(posnum(1), :), ...
                corrSum_RawChange_A1(posnum(1), :), corrAsal_RawChange_A1(posnum(1), :), corrBsal_RawChange_A1(posnum(1), :), ...
                corrArank_RawChange_A1(posnum(1), :), corrBrank_RawChange_A1(posnum(1), :), corrCrank_RawChange_A1(posnum(1), :), corrUrank_RawChange_A1(posnum(1), :), ...
                corrRankCUratio_RawChange_A1(posnum(1), :), corrRankABdiff_RawChange_A1(posnum(1), :),  corrRankabsdiff_RawChange_A1(posnum(1), :), ...
                corrRanksum_RawChange_A1(posnum(1), :)] = slidereg_funChange(mfrA11, invmfrA11, val_A_Rat,  val_B_Rat,  val_C_Rat, ...
                val_U_Rat,  val_A_Rank, val_B_Rank, val_C_Rank, val_U_Rank, val_2AFC, shortRT_AFC, val_diffAB, invval_A_Rat,  invval_B_Rat,  invval2_C_Rat, ...
                invval2_U_Rat,  invval_A_Rank, invval_B_Rank, invval2_C_Rank, invval2_U_Rank, invval_2AFC, invshortRT_AFC, invval_diffAB);
            
            
            
            [pvalVal_RawChange_R, corrVal_RawChange_R, pvalSal_RawChange_R, corrSal_RawChange_R, ...
                pvalSign_RawChange_R, corrSign_RawChange_R, pvalRank_RawChange_R, corrRank_RawChange_R, ...
                pvalRT_RawChange_R, corrRT_RawChange_R] = regRating_funChange(mfrRat1mean, invmfrRat1mean, val_rat_Rat, val_rat_Rank,  RTmean_Rat,...
                invval_rat_Rat, invval_rat_Rank,  invRTmean_Rat);
            
            toc(oooo)
            
            
        end
        
    end
    
    
    
    if exist([folder_to_save '/regresspval_' inputRegion ,'.mat'])
        delete([folder_to_save '/regresspval_' inputRegion ,'.mat']);
    end
    
    h=0; g=0;
    while h<1
        try
            load([folder_to_save '/regresspval_' inputRegion ])
            h=1;
        catch
            save ([folder_to_save '/regresspval_' inputRegion ],  'pval*', '-v7.3')
            if g>0
                warning(['error while saving the pvalues in region ', inputRegion])
            end
            g=g+1;
        end
    end
    
    
    
    
    if exist([folder_to_save '/regressreg_' inputRegion ,'.mat'])
        delete([folder_to_save '/regressreg_' inputRegion ,'.mat' ]);
    end
    h=0; g=0;
    while h<1
        
        try
            load([folder_to_save '/regressreg_' inputRegion ])
            h=1;
        catch
            save ([folder_to_save '/regressreg_' inputRegion ],  'corr*', '-v7.3')
            if g>0
                warning(['error while saving the corr in region ', inputRegion])
            end
            g=g+1;
        end
    end
    
    
    
    
    clear matCorrPerm matPvalPerm
    clear pval* corr* change* MFR_m* CPD* change*
    clear -regexp ^pval ^corr ^change  ^MFR_m ^CPD ^change
    
end
cd (folder_to_save)

% pool_all_regionsChange('ALL')
pool_all_regions('ALL')














function [pvalA, pvalB, pvalAsal, pvalBsal, pvalC,  pvalU, pvalabsdiff, pvaldiffAB,  pvalAbspref, pvalpref,  pvalSum, pvalUCratio, ...
    pvalArank,pvalBrank, pvalCrank,pvalUrank, pvalRankabsdiff, pvalRanksum, pvalRankABdiff, pvalRankCUratio, ...
    corrA, corrB, corrAsal, corrBsal, corrC,  corrU,  corrabsdiff, corrdiffAB, corrAbspref, corrpref, corrSum, corrUCratio, ...
    corrArank,corrBrank, corrCrank, corrUrank, corrRankabsdiff, corrRanksum, corrRankABdiff, corrRankCUratio, ...
    pvalRT, corrRT,...
    pvalflavourA, pvalflavourB, pvalflavourC, pvalflavourU, pvalABbin, pvalSignA, pvalSignB, pvalSignC, pvalSignU,...
    corrflavourA, corrflavourB, corrflavourC, corrflavourU, corrABbin, corrSignA, corrSignB, corrSignC, corrSignU] = slidereg_fun(mfrA11, mfrA12, val_A_Rat,  val_B_Rat,  val_C_Rat, val_U_Rat,  val_A_Rank, val_B_Rank, val_C_Rank, val_U_Rank, val_2AFC, RT, flavourA, flavourB,  flavourC, flavourU)

global typeCorr modeltype
dbstop if error
var2use= val_A_Rat;

[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalA]= c;
[corrA]= a;

var2use= val_B_Rat;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalB]= c;
[corrB]= a;


% salience
var2use=abs(val_A_Rat);
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalAsal]= c;
[corrAsal]= a;

var2use=abs(val_B_Rat);
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalBsal]= c;
[corrBsal]= a;

var2use=abs(val_2AFC);
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalAbspref]= c;
[corrAbspref]= a;

var2use=val_2AFC;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalpref]= c;
[corrpref]= a;

var2use=val_A_Rat + val_B_Rat;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalSum]= c;
[corrSum]= a;

var2use=val_C_Rat;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalC]= c;
[corrC]= a;

var2use=val_U_Rat;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalU]= c;
[corrU]= a;


% binary cells firing more for left chosen or A chosen
var2use= val_2AFC>0;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   ['bin',  typeCorr]);  % not exactly pearson, remember that in hte function corr_mixed it converts the mfr to tiedrank so that it is ranked only the dependetn, and not the dummy
assert(size(c,2)==size(mfrA11,2));
[pvalABbin]= c;
[corrABbin]= a;


% differences based on mean rating
var2use= val_B_Rat-val_A_Rat;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvaldiffAB]= c;
[corrdiffAB]= a;

var2use= abs(val_B_Rat-val_A_Rat);
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalabsdiff]= c;
[corrabsdiff]= a;

var2use= RT;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalRT]= c;
[corrRT]= a;


try
    var2use= val_U_Rat./val_C_Rat;
    [a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
    assert(size(c,2)==size(mfrA11,2));
    [pvalUCratio]= c;
    [corrUCratio]= a;
catch
    [pvalUCratio]= nan(size(c));
    [corrUCratio]= nan(size(a));
end

if nargin>12
    var2use=val_A_Rat>0;
    [a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   ['bin',  typeCorr]);
    assert(size(c,2)==size(mfrA11,2));
    [pvalSignA]= c;
    [corrSignA]= a;
    
    
    
    var2use=val_B_Rat>0;
    [a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   ['bin',  typeCorr]);
    assert(size(c,2)==size(mfrA11,2));
    [pvalSignB]= c;
    [corrSignB]= a;
    
    
    var2use=val_C_Rat>0;
    [a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   ['bin',  typeCorr]);
    assert(size(c,2)==size(mfrA11,2));
    [pvalSignC]= c;
    [corrSignC]= a;
    
    var2use= val_U_Rat>0;
    [a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   ['bin',  typeCorr]);
    assert(size(c,2)==size(mfrA11,2));
    [pvalSignU]= c;
    [corrSignU]= a;
    
    
    var2use=flavourA;
    [a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   ['bin',  typeCorr]);  % not exactly pearson, remember that in hte function corr_mixed it converts the mfr to tiedrank so that it is ranked only the dependetn, and not the dummy
    assert(size(c,2)==size(mfrA11,2));
    [pvalflavourA]= c;
    [corrflavourA]= a;
    
    var2use=flavourB;
    [a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   ['bin',  typeCorr]);
    assert(size(c,2)==size(mfrA11,2));
    [pvalflavourB]= c;
    [corrflavourB]= a;
    
    var2use=flavourC;
    [a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   ['bin',  typeCorr]);
    assert(size(c,2)==size(mfrA11,2));
    [pvalflavourC]= c;
    [corrflavourC]= a;
    
    var2use=flavourU;
    [a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   ['bin',  typeCorr]);
    assert(size(c,2)==size(mfrA11,2));
    [pvalflavourU]= c;
    [corrflavourU]= a;
end

% % % % ranking

var2use=val_A_Rank;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalArank]= c;
[corrArank]= a;

var2use=val_B_Rank;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalBrank]= c;
[corrBrank]= a;

var2use=val_C_Rank;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalCrank]= c;
[corrCrank]= a;

var2use=val_U_Rank;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalUrank]= c;
[corrUrank]= a;

try
    var2use= val_U_Rank ./ val_C_Rank;
    [a,c]=corr_mixed(modeltype, mfrA11, mfrA12,  var2use(:,1), var2use(:,2), typeCorr);
    assert(size(c,2)==size(mfrA11,2));
    [pvalRankCUratio]= c;
    [corrRankCUratio]= a;
catch
    [pvalRankCUratio]= nan(size(c));
    [corrRankCUratio]= nan(size(a));
end

var2use=val_A_Rank - val_B_Rank;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalRankABdiff]= c;
[corrRankABdiff]= a;

var2use=val_A_Rank + val_B_Rank;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalRanksum]= c;
[corrRanksum]= a;

var2use=abs(val_A_Rank - val_B_Rank);
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalRankabsdiff]= c;
[corrRankabsdiff]= a;


































function [mfr1, mfr2] = delete_low_frequencies(mfr1, mfr2)

mfr1= nan(length(mfr1),1);
mfr2= nan(length(mfr2),1);































function [histoconv2] = kernelgauss (data, kern_width, resolution, stepsize, analyseswindow, stdkernel)

sizeshift = 500; % in miliseconds

analyseswindow(1) =analyseswindow(1) - sizeshift ;
analyseswindow(2) =analyseswindow(2) + sizeshift ;

t = analyseswindow(2) - analyseswindow(1) ; % add sizeshift ms for estimation and then delete them from estimates

clear histo
histo = nan(size(data,1),t/resolution);

for a = 1:size(data,1)
    %     histo(a,:) = histc(data{a,1},(-((t/2)-(resolution/2)):resolution:((t/2)-(resolution/2))));
    histo(a,:) = histc(data{a,1},((analyseswindow(1))-(resolution/2):resolution:(analyseswindow(2)-1-resolution/2)));  %% from beginnging to end and delete last ms.
end
mfr=1;


kern_length = kern_width * stdkernel;
kernel = normpdf(-kern_length:resolution:kern_length, 0, kern_width);
y= histo * ((1000/resolution)/mfr);
%
clear histoconv
for kkk=1:size(histo, 1)
    histoconv(kkk,:)= conv(y(kkk,:), kernel, 'same');
end

histoconv2= histoconv(:, sizeshift:stepsize:end-sizeshift); % add 500 ms fro estimation and then delete them from estimates





























function [pvalCregVal_A1, corrCregVal_A1, pvalUregVal_A1, corrUregVal_A1, pvalCval_A1, corrCval_A1, pvalUval_A1, corrUval_A1, pvalmodelValCU, CPD_CregVal_A1, CPD_UregVal_A1,   CPDvar_CregVal_A1 , CPDvar_UregVal_A1] = bilinearregression (mfrA11,  mfrA12,   var2use, var2use2)

for pp = 1:size(mfrA11,2)
    
    y1=zscore(tiedrank(mfrA11(:,pp)));
    stats = regstats(y1, [zscore(tiedrank(var2use(:,1))), zscore(tiedrank(var2use2(:,1)))],'linear',{'tstat', 'rsquare', 'mse', 'fstat'});
    
    pvalCregVal_A1(1,pp)= stats.tstat.pval(2);
    corrCregVal_A1(1,pp)= stats.tstat.beta(2);
    pvalUregVal_A1(1,pp)= stats.tstat.pval(3);
    corrUregVal_A1(1,pp)= stats.tstat.beta(3);
    SST(pp,1) = sum((y1- mean(y1)).^2);
    SSE(pp,1) = stats.fstat.sse;
    Rsquared(pp,1)= stats.rsquare;
    pvalmodelValCU(1,pp)= stats.fstat.pval;
    
    
    
    stats = regstats(y1, [zscore(tiedrank(var2use(:,1)))],'linear',{'tstat', 'rsquare', 'mse', 'fstat'});
    pvalCval_A1(1,pp)= stats.tstat.pval(2);
    corrCval_A1(1,pp)= stats.tstat.beta(2);
    SSEC(pp,1) = stats.fstat.sse;
    
    
    stats = regstats(y1, [zscore(tiedrank(var2use2(:,1)))],'linear',{'tstat', 'rsquare', 'mse', 'fstat'});
    pvalUval_A1(1,pp)= stats.tstat.pval(2);
    corrUval_A1(1,pp)= stats.tstat.beta(2);
    SSEU(pp,1) = stats.fstat.sse;
    
    %%
    
    
    y2=zscore(tiedrank(mfrA12(:,pp)));
    stats = regstats(y2, [zscore(tiedrank(var2use(:,2))), zscore(tiedrank(var2use2(:,2)))],'linear',{'tstat', 'rsquare', 'mse', 'fstat'});
    pvalCregVal_A1(2,pp)= stats.tstat.pval(2);
    corrCregVal_A1(2,pp)= stats.tstat.beta(2);
    pvalUregVal_A1(2,pp)= stats.tstat.pval(3);
    corrUregVal_A1(2,pp)= stats.tstat.beta(3);
    SST(pp,2) = sum((y2- mean(y2)).^2);
    SSE(pp,2) = stats.fstat.sse;
    Rsquared(pp,2)= stats.rsquare;
    pvalmodelValCU(2,pp)= stats.fstat.pval;
    
    
    
    stats = regstats(y2, [zscore(tiedrank(var2use(:,2)))],'linear',{'tstat', 'rsquare', 'mse', 'fstat'});
    pvalCval_A1(2,pp)= stats.tstat.pval(2);
    corrCval_A1(2,pp)= stats.tstat.beta(2);
    SSEC(pp,2) = stats.fstat.sse;
    
    
    stats = regstats(y2, [zscore(tiedrank(var2use2(:,2)))],'linear',{'tstat', 'rsquare', 'mse', 'fstat'});
    pvalUval_A1(2,pp)= stats.tstat.pval(2);
    corrUval_A1(2,pp)= stats.tstat.beta(2);
    SSEU(pp,2) = stats.fstat.sse;
    
end



CPD_UregVal_A1([1:2] ,:)= ((SSEC -SSE) ./ SSEC)';
CPD_CregVal_A1([1:2],:)= ((SSEU -SSE) ./ SSEU)';


CPDvar_UregVal_A1([1:2] ,:)= ((SSEC -SSE)  ./ SST)' ;
CPDvar_CregVal_A1([1:2],:)= ((SSEU -SSE)  ./ SST )';































function [pvalCRat, pvalCRank, pvalURat, pvalURank, ...
    pvalratioCURat, pvalratioCURank,  ...
    pvalRTmean_R,...
    corrCRat, corrCRank,  corrURat, corrURank,  ...
    corrratioCURat, corrratioCURank,  ...
    corrRTmean_R] = slidereg_funNew(mfrA11, mfrA12,  mfrRat1, mfrRat2, val_C_Rat, val_U_Rat,  val_C_Rank, val_U_Rank, RTmean_Rat)

global typeCorr modeltype

dbstop if error

var2use= val_C_Rat;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalCRat]= c;
[corrCRat]= a;


var2use= val_C_Rank;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalCRank]= c;
[corrCRank]= a;



var2use= val_U_Rat;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalURat]= c;
[corrURat]= a;


var2use= val_U_Rank;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalURank]= c;
[corrURank]= a;






try
    var2use= val_U_Rat./val_C_Rat;
    [a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
    assert(size(c,2)==size(mfrA11,2));
    [pvalratioCURat]= c;
    [corrratioCURat]= a;
catch
    [pvalratioCURat]= nan(size(c));
    [corrratioCURat]= nan(size(a));
end



try
    var2use= val_U_Rank./val_C_Rank;
    [a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
    assert(size(c,2)==size(mfrA11,2));
    [pvalratioCURank]= c;
    [corrratioCURank]= a;
catch
    [pvalratioCURank]= nan(size(c));
    [corrratioCURank]= nan(size(a));
end

try
    var2use= val_U_Tot./val_C_Tot;
    [a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
    assert(size(c,2)==size(mfrA11,2));
    [pvalratioCUTot]= c;
    [corrratioCUTot]= a;
catch
    [pvalratioCUTot]= nan(size(c));
    [corrratioCUTot]= nan(size(a));
end


var2use= RTmean_Rat;
[a,c]=corr_mixed(modeltype, mfrRat1, mfrRat2, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrRat1,2));
[pvalRTmean_R]= c;
[corrRTmean_R]= a;


































function  [pvalRed_A1, pvalGreen_A1,  pvalBlue_A1,  pvalContrast_A1,  pvalLuminance_A1,  pvalHue_A1,  pvalSaturation_A1, ...
    corrRed_A1, corrGreen_A1,  corrBlue_A1,  corrContrast_A1,  corrLuminance_A1,  corrHue_A1,  corrSaturation_A1, ...
    pvalRed_Productwise_A1, pvalGreen_Productwise_A1,  pvalBlue_Productwise_A1,  pvalContrast_Productwise_A1,  pvalLuminance_Productwise_A1,  pvalHue_Productwise_A1,  pvalSaturation_Productwise_A1, ...
    corrRed_Productwise_A1, corrGreen_Productwise_A1,  corrBlue_Productwise_A1,  corrContrast_Productwise_A1,  corrLuminance_Productwise_A1,  corrHue_Productwise_A1,  corrSaturation_Productwise_A1] = corrVisual(mfrA11, mfrA12,  ...
    Red, Green, Blue, Contrast, Luminance, Hue, Saturation, stimOrderA, Red_A, Green_A, Blue_A, Contrast_A, Luminance_A, Hue_A, Saturation_A)


global modeltype typeCorr

var2use= Red_A;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalRed_A1]= c;
[corrRed_A1]= a;

var2use= Green_A;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalGreen_A1]= c;
[corrGreen_A1]= a;


var2use= Blue_A;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalBlue_A1]= c;
[corrBlue_A1]= a;



var2use= Contrast_A;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalContrast_A1]= c;
[corrContrast_A1]= a;


var2use= Luminance_A;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalLuminance_A1]= c;
[corrLuminance_A1]= a;


var2use= Hue_A;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalHue_A1]= c;
[corrHue_A1]= a;


var2use= Saturation_A;
[a,c]=corr_mixed(modeltype, mfrA11, mfrA12, var2use(:,1), var2use(:,2),   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalSaturation_A1]= c;
[corrSaturation_A1]= a;






MFR1s = nan(20,size(mfrA11,2));
MFR2s = nan(20,size(mfrA11,2));


for index =  1:20
    
    for jj= 1: size(mfrA11,2)
        MFR1s (index, jj) = nanmean(mfrA11(stimOrderA(:,1) == index, jj));
        MFR2s (index,jj) = nanmean(mfrA12(stimOrderA(:,2) == index, jj));
    end
end


var2use= Red;
[a,c]=corr_mixed(modeltype, MFR1s, MFR2s, var2use, var2use,   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalRed_Productwise_A1]= c;
[corrRed_Productwise_A1]= a;

var2use= Green;
[a,c]=corr_mixed(modeltype, MFR1s, MFR2s, var2use, var2use,   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalGreen_Productwise_A1]= c;
[corrGreen_Productwise_A1]= a;


var2use= Blue;
[a,c]=corr_mixed(modeltype, MFR1s, MFR2s, var2use, var2use,   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalBlue_Productwise_A1]= c;
[corrBlue_Productwise_A1]= a;



var2use= Contrast;
[a,c]=corr_mixed(modeltype, MFR1s, MFR2s, var2use, var2use,   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalContrast_Productwise_A1]= c;
[corrContrast_Productwise_A1]= a;


var2use= Luminance;
[a,c]=corr_mixed(modeltype, MFR1s, MFR2s, var2use, var2use,   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalLuminance_Productwise_A1]= c;
[corrLuminance_Productwise_A1]= a;


var2use= Hue;
[a,c]=corr_mixed(modeltype, MFR1s, MFR2s, var2use, var2use,   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalHue_Productwise_A1]= c;
[corrHue_Productwise_A1]= a;


var2use= Saturation;
[a,c]=corr_mixed(modeltype, MFR1s, MFR2s, var2use, var2use,   typeCorr);
assert(size(c,2)==size(mfrA11,2));
[pvalSaturation_Productwise_A1]= c;
[corrSaturation_Productwise_A1]= a;































function [pvalVal, corrVal, pvalSal, corrSal, pvalflavourRat, corrflavourRat, ...
    pvalSign, corrSign, pvalRank, corrRank, ...
    pvalRT, corrRT] = regRating_fun( mfrRat1, mfrRat2, val_rat1, val_rat2,  flavourRat1, flavourRat2, val_rank1, val_rank2,  RT)



dbstop if error
global typeCorr modeltype
% if ==1 then corr, if ==2 then mixed-model



[a,c]=corr_mixed(modeltype, mfrRat1, mfrRat2, val_rat1, val_rat2,   typeCorr);
assert(size(c,2)==size(mfrRat1,2)); % just in case
[pvalVal]= c;
[corrVal]= a;


[a,c]=corr_mixed(modeltype, mfrRat1, mfrRat2, abs(val_rat1), abs(val_rat2),  typeCorr);
assert(size(c,2)==size(mfrRat1,2)); % just in case
[pvalSal]= c;
[corrSal]= a;



[a,c]=corr_mixed(modeltype, mfrRat1, mfrRat2, flavourRat1, flavourRat2,   ['bin',  typeCorr]);
assert(size(c,2)==size(mfrRat1,2)); % just in case
[pvalflavourRat]=  c;
[corrflavourRat]=  a;


[a,c]=corr_mixed(modeltype, mfrRat1, mfrRat2, val_rat1>0, val_rat2>0,  ['bin',  typeCorr]);
assert(size(c,2)==size(mfrRat1,2)); % just in case
[pvalSign]=  c;
[corrSign]=  a;


[a,c]=corr_mixed(modeltype, mfrRat1, mfrRat2, val_rank1, val_rank2,   typeCorr);
assert(size(c,2)==size(mfrRat1,2)); % just in case
[pvalRank]= c;
[corrRank]= a;


[a,c]=corr_mixed(modeltype, mfrRat1, mfrRat2, RT(:,1), RT(:,2),  typeCorr);
assert(size(c,2)==size(mfrRat1,2)); % just in case
[pvalRT]= c;
[corrRT]= a;
















function [pvalAbspref,pvalLRpref, pvalABpref, pvalRT,  pvalABbin, ...
    pvalAval, pvalBval, pvalCval, pvalUval,  pvalUCratio, pvaldiffAB, pvalabsdiff,...
    pvalSum, pvalAsal, pvalBsal, pvalArank, pvalBrank, pvalCrank, pvalUrank, ...
    pvalRankCUratio, pvalRankABdiff,  pvalRankabsdiff, pvalRanksum, ...
    corrAbspref,corrLRpref, corrABpref, corrRT,  corrABbin, ...
    corrAval, corrBval, corrCval, corrUval,  corrUCratio, corrdiffAB, corrabsdiff,...
    corrSum, corrAsal, corrBsal, corrArank, corrBrank, corrCrank, corrUrank,  ...
    corrRankCUratio, corrRankABdiff,  corrRankabsdiff, corrRanksum] = slidereg_funChange(mfr, mfrinv, val_A_Rat,  val_B_Rat,  val_C_Rat, ...
    val_U_Rat,  val_A_Rank, val_B_Rank, val_C_Rank, val_U_Rank, val_2AFC, RT, val_diffAB, invval_A_Rat,  invval_B_Rat,  invval_C_Rat, ...
    invval_U_Rat,  invval_A_Rank, invval_B_Rank, invval_C_Rank, invval_U_Rank, invval_2AFC, invRT, invval_diffAB)



global typeCorr modeltype

dbstop if error
% values

mfr = mfr - mfrinv;

var2use= val_A_Rat - invval_A_Rat;
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvalAval]= c;
[corrAval]= a;

var2use= val_B_Rat - invval_B_Rat;
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvalBval]= c;
[corrBval]= a;


% salience
var2use=abs(val_A_Rat) - abs(invval_A_Rat);
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvalAsal]= c;
[corrAsal]= a;

var2use=abs(val_B_Rat) -abs(invval_B_Rat);
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvalBsal]= c;
[corrBsal]= a;

var2use=abs(val_2AFC) -abs(invval_2AFC);
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvalAbspref]= c;
[corrAbspref]= a;

var2use=val_2AFC - invval_2AFC;
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvalLRpref]= c;
[corrLRpref]= a;


var2use=val_diffAB - invval_diffAB;
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvalABpref]= c;
[corrABpref]= a;


var2use=val_C_Rat - invval_C_Rat;
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvalCval]= c;
[corrCval]= a;

var2use=val_U_Rat - invval_U_Rat;
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvalUval]= c;
[corrUval]= a;



% binary cells firing more for left chosen or A chosen
var2use= (val_2AFC>0) - (invval_2AFC>0);
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c'); % not exactly pearson, remember that in hte function corr_mixed it converts the mfr to tiedrank so that it is ranked only the dependetn, and not the dummy
assert(size(c,1)==size(mfr,2));
[pvalABbin]= c;
[corrABbin]= a;


% differences based on mean rating
var2use=(val_A_Rat + val_B_Rat) - (invval_A_Rat + invval_B_Rat);
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvalSum]= c;
[corrSum]= a;



var2use= (val_B_Rat-val_A_Rat) - (invval_B_Rat-invval_A_Rat);
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvaldiffAB]= c;
[corrdiffAB]= a;

var2use= abs(val_B_Rat-val_A_Rat) - abs(invval_B_Rat-invval_A_Rat);
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvalabsdiff]= c;
[corrabsdiff]= a;

var2use= RT - invRT;
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvalRT]= c;
[corrRT]= a;



try
    var2use= (val_U_Rat./val_C_Rat) - (invval_U_Rat./invval_C_Rat);
    [a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
    assert(size(c,1)==size(mfr,2));
    [pvalUCratio]= c;
    [corrUCratio]= a;
catch
    [pvalUCratio]= nan(size(c));
    [corrUCratio]= nan(size(a));
end



% % % % ranking

var2use=val_A_Rank - invval_A_Rank;
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvalArank]= c;
[corrArank]= a;

var2use=val_B_Rank - invval_B_Rank;
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvalBrank]= c;
[corrBrank]= a;

var2use=val_C_Rank - invval_C_Rank;
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvalCrank]= c;
[corrCrank]= a;

var2use=val_U_Rank - invval_U_Rank;
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvalUrank]= c;
[corrUrank]= a;


try
    var2use= (val_U_Rank ./ val_C_Rank) -(invval_U_Rank ./ invval_C_Rank);
    [a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
    assert(size(c,1)==size(mfr,2));
    [pvalRankCUratio]= c;
    [corrRankCUratio]= a;
catch
    [pvalRankCUratio]= nan(size(c));
    [corrRankCUratio]= nan(size(a));
end

var2use=(val_A_Rank - val_B_Rank) - (invval_A_Rank - invval_B_Rank);
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvalRankABdiff]= c;
[corrRankABdiff]= a;

var2use=(val_A_Rank + val_B_Rank) - (invval_A_Rank + invval_B_Rank);
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvalRanksum]= c;
[corrRanksum]= a;

var2use=abs(val_A_Rank - val_B_Rank) - abs(invval_A_Rank - invval_B_Rank);
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,1)==size(mfr,2));
[pvalRankabsdiff]= c;
[corrRankabsdiff]= a;

























function [pvalVal, corrVal, pvalSal, corrSal, ...
    pvalSign, corrSign, pvalRank, corrRank, ...
    pvalRT, corrRT] = regRating_funChange(mfr, mfrinv, val_rat, val_rank,  RT,...
    invval_rat, invval_rank,  invRT);


global typeCorr modeltype
dbstop if error

mfr = mfr - mfrinv;



var2use= val_rat - invval_rat;
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,2)==size(mfrRat1,2)); % just in case
[pvalVal]= c;
[corrVal]= a;

var2use= abs(val_rat) - abs(invval_rat);
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,2)==size(mfrRat1,2)); % just in case
[pvalSal]= c;
[corrSal]= a;


var2use= val_rat>0 - invval_rat>0;
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,2)==size(mfrRat1,2)); % just in case
[pvalSign]=  c;
[corrSign]=  a;


var2use= val_rank - invval_rank;
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,2)==size(mfrRat1,2)); % just in case
[pvalRank]= c;
[corrRank]= a;


var2use= RT - invRT;
[a,c]= corr(mfr, var2use(:,1), 'type',  typeCorr, 'rows', 'c');
assert(size(c,2)==size(mfrRat1,2)); % just in case
[pvalRT]= c;
[corrRT]= a;






















function [a,c, mm, nn] = corr_mixed(type, mfrA11, mfrA12, val1, val2, typeCorr, convtaste, flav11, flav12)

trialnum1 = zscore(tiedrank(linspace(1,size(val1,1),size(val1,1))))';
trialnum2 = zscore(tiedrank(linspace(1,size(val2,1),size(val2,1))))';
if type == 2
    if strcmp(typeCorr,'Spearman') || strcmp(typeCorr,'Pearson') || strcmp(typeCorr,'Kendall')
        if nargin <= 7 || convtaste==3
            if size(val1,2)==1
                val1=repmat(val1, 1, size(mfrA11,2));
                val2=repmat(val2, 1, size(mfrA11,2));
                flav11=repmat(flav11, 1, size(mfrA11,2));
                flav12=repmat(flav12, 1, size(mfrA11,2));
            end
            
            s2u= size(mfrA11,1);
            
            for gg=1:size(mfrA11,2)
                if strcmp(typeCorr, 'Spearman')
                    mfr = [zscore(tiedrank(mfrA11(:,gg))); zscore(tiedrank(mfrA12(:,gg)))];
                    val = [zscore(tiedrank(val1(:,gg)));   zscore(tiedrank(val2(:,gg)))];
                    
                    flav = [zscore(tiedrank(flav11(:,gg)));   zscore(tiedrank(flav12(:,gg)))];
                elseif  strcmp(typeCorr, 'Pearson')
                    mfr = [zscore(mfrA11(:,gg)); zscore(mfrA12(:,gg))];
                    val = [zscore(val1(:,gg));   zscore(val2(:,gg))];
                    flav = [zscore(tiedrank(flav11(:,gg)));   zscore(tiedrank(flav12(:,gg)))];
                end
                
                try
                    [b1, ~, p1]=correlationfunctions.commoncorrSpearman (mfr(1:s2u),[val(1:s2u),flav(1:s2u)]);
                    [b2, ~, p2]=correlationfunctions.commoncorrSpearman (mfr(1+s2u:end),[val(1+s2u:end),flav(1+s2u:end)]);
                    
                    a([1,3,5,2,4,6],gg) = [b1; b2];
                    c([1,3,5,2,4,6],gg) = [p1; p2];
                    
                    
                catch  ME
                    if strcmp(ME.identifier, 'MATLAB:catenate:dimensionMismatch') || strcmp(ME.identifier, 'stats:regstats:InputSizeMismatch')
                        rethrow(ME);
                        
                    else
                        a(1,gg) = nan;
                        c(1,gg) = nan;
                        a(2,gg) = nan;
                        c(2,gg) = nan;
                    end
                end
            end
        end
        
    elseif strcmp(typeCorr,'binSpearman')  || strcmp(typeCorr,'binPearson')  || strcmp(typeCorr,'binKendall')
        
        
        if size(val1,2)==1
            val1=repmat(val1, 1, size(mfrA11,2));
            val2=repmat(val2, 1, size(mfrA11,2));
        end
        
        
        for gg=1:size(mfrA11,2)
            
            if strcmp(typeCorr, 'binSpearman')
                mfr = [zscore(tiedrank(mfrA11(:,gg))); zscore(tiedrank(mfrA12(:,gg)))];
                val = [zscore(tiedrank(val1(:,gg)));   zscore(tiedrank(val2(:,gg)))];
            elseif  strcmp(typeCorr, 'binPearson')
                mfr = [zscore(mfrA11(:,gg)); zscore(mfrA12(:,gg))];
                val = [zscore(val1(:,gg));   zscore(val2(:,gg))];
            end
            sess= [ones(size(val1,1),1) * -0.5 ; ones(size(val2,1),1) * 0.5];
            
            try
                stats = regstats(mfr,[val, sess],'interaction',{'tstat'});
                a(1,gg) = stats.tstat.beta(2);
                c(1,gg) = stats.tstat.pval(2);
                a(2,gg) = stats.tstat.beta(4);
                c(2,gg) = stats.tstat.pval(4);
            catch
                if strcmp(ME.identifier, 'MATLAB:catenate:dimensionMismatch') || strcmp(ME.identifier, 'stats:regstats:InputSizeMismatch')
                    rethrow(ME)
                else
                    a(1,gg) = nan;
                    c(1,gg) = nan;
                    a(2,gg) = nan;
                    c(2,gg) = nan;
                end
            end
        end
        
    end
    
    
elseif type == 1
    
    
    if nargin >7 && convtaste==1
        
        
        for k=1: size(val1,2)
            ind11=find(flav11(:,k)==1);
            ind12=find(flav12(:,k)==1);
            
            val11(:,k)= val1(ind11,k);
            val21(:,k)= val2(ind12,k);
            mfrA111(:,k)= mfrA11(ind11,k);
            mfrA121(:,k)= mfrA12(ind12,k);
        end
        
        val1   = val11;
        val2   = val21;
        mfrA11 = mfrA111;
        mfrA12 = mfrA121;
        
    end
    
    
    if nargin <= 7 || convtaste == 1 || convtaste == 3
        if size(val1,2)==1
            if strcmp(typeCorr,'Spearman') || strcmp(typeCorr,'Pearson') || strcmp(typeCorr,'Kendall')
                
                [a(1,:),c(1,:)]=corr([mfrA11], [val1], 'type', typeCorr, 'rows', 'c');
                [a(2,:),c(2,:)]=corr([mfrA12], [val2], 'type', typeCorr, 'rows', 'c');
                
                
                
            elseif strcmp(typeCorr,'binSpearman')  || strcmp(typeCorr,'binPearson')  || strcmp(typeCorr,'binKendall')
                
                [a(1,:),c(1,:)]=corr([mfrA11], [val1], 'type', typeCorr(4:end), 'rows', 'c');
                [a(2,:),c(2,:)]=corr([mfrA12], [val2], 'type', typeCorr(4:end), 'rows', 'c');
                
            end
            
        else
            
            if strcmp(typeCorr,'Spearman') || strcmp(typeCorr,'Pearson') || strcmp(typeCorr,'Kendall')
                
                [diag1c, diag1p] = corr([mfrA11],  [val1], 'type', typeCorr, 'rows', 'c');
                [diag2c, diag2p] = corr([mfrA12],  [val2], 'type', typeCorr, 'rows', 'c');
                a(1,:) = diag(diag1c); c(1,:) = diag(diag1p);
                a(2,:) = diag(diag2c); c(2,:) = diag(diag2p);
            elseif strcmp(typeCorr,'binSpearman')  || strcmp(typeCorr,'binPearson')  || strcmp(typeCorr,'binKendall')
                
                [diag1c, diag1p] = corr([mfrA11], [val1], 'type', typeCorr(4:end), 'rows', 'c');
                [diag2c, diag2p] = corr([mfrA12], [val2], 'type', typeCorr(4:end), 'rows', 'c');
                
                a(1,:) = diag(diag1c); c(1,:) = diag(diag1p);
                a(2,:) = diag(diag2c); c(2,:) = diag(diag2p);
            end
        end
        
    end
    
    
    if nargin >7 && convtaste==2
        
        
        
        for k=1: size(val1,2)
            ind11=find(flav11(:,k)==1);
            ind12=find(flav12(:,k)==1);
            
            val11{:,k}= val1(ind11,k);
            val21{:,k}= val2(ind12,k);
            mfrA111{:,k}= mfrA11(ind11,k);
            mfrA121{:,k}= mfrA12(ind12,k);
        end
        
        val1   = val11;
        val2   = val21;
        mfrA11 = mfrA111;
        mfrA12 = mfrA121;
        
        for gg=1:size(mfrA11,2)
            [diag1c, diag1p] = corr(mfrA11{gg},  val1{gg}, 'type', typeCorr, 'rows', 'c');
            [diag2c, diag2p] = corr(mfrA12{gg},  val2{gg}, 'type', typeCorr, 'rows', 'c');
            a(1,gg) = diag1c; c(1,gg) = diag1p;
            a(2,gg) = diag2c; c(2,gg) = diag2p;
        end
        
    end
    
elseif type == 3
    
    
    if strcmp(typeCorr,'Spearman') || strcmp(typeCorr,'Pearson') || strcmp(typeCorr,'Kendall')
        
        if nargin >7 && convtaste==1
            
            for k=1: size(val1,2)
                ind11=find(flav11(:,k)==1);
                ind12=find(flav12(:,k)==1);
                
                val11(:,k)= val1(ind11,k);
                val21(:,k)= val2(ind12,k);
                mfrA111(:,k)= mfrA11(ind11,k);
                mfrA121(:,k)= mfrA12(ind12,k);
            end
            
            val1   = val11;
            val2   = val21;
            mfrA11 = mfrA111;
            mfrA12 = mfrA121;
            
        end
        
        if nargin <= 7 || convtaste==1
            if size(val1,2)==1
                val1=repmat(val1, 1, size(mfrA11,2));
                val2=repmat(val2, 1, size(mfrA11,2));
            end
            
            
            for gg=1:size(mfrA11,2)
                if strcmp(typeCorr, 'Spearman')
                    mfr = [zscore(tiedrank(mfrA11(:,gg))); zscore(tiedrank(mfrA12(:,gg)))];
                    val = [zscore(tiedrank(val1(:,gg)));   zscore(tiedrank(val2(:,gg)))];
                elseif  strcmp(typeCorr, 'Pearson')
                    mfr = [zscore(mfrA11(:,gg)); zscore(mfrA12(:,gg))];
                    val = [zscore(val1(:,gg));   zscore(val2(:,gg))];
                end
                sess= [ones(size(val1,1),1) * -0.5 ; ones(size(val2,1),1) * 0.5];
                try
                    stats = regstats(mfr,[val, sess],'interaction',{'tstat'});
                    a(1,gg) = stats.tstat.beta(2);
                    c(1,gg) = stats.tstat.pval(2);
                    a(2,gg) = stats.tstat.beta(4);
                    c(2,gg) = stats.tstat.pval(4);
                    
                catch  ME
                    if strcmp(ME.identifier, 'MATLAB:catenate:dimensionMismatch') || strcmp(ME.identifier, 'stats:regstats:InputSizeMismatch')
                        rethrow(ME);
                        
                    else
                        a(1,gg) = nan;
                        c(1,gg) = nan;
                        a(2,gg) = nan;
                        c(2,gg) = nan;
                    end
                end
            end
        end
        
        
        if nargin >7 && convtaste==2
            
            for k=1: size(val1,2)
                ind11=find(flav11(:,k)==1);
                ind12=find(flav12(:,k)==1);
                
                val11{:,k}= val1(ind11,k);
                val21{:,k}= val2(ind12,k);
                mfrA111{:,k}= mfrA11(ind11,k);
                mfrA121{:,k}= mfrA12(ind12,k);
            end
            
            val1   = val11;
            val2   = val21;
            mfrA11 = mfrA111;
            mfrA12 = mfrA121;
            
            
            
            if size(val1,2)==1
                val1=repmat(val1, 1, size(mfrA11,2));
                val2=repmat(val2, 1, size(mfrA11,2));
            end
            
            for gg=1:size(mfrA11,2)
                if strcmp(typeCorr, 'Spearman')
                    mfr = [zscore(tiedrank(mfrA11{:,gg})); zscore(tiedrank(mfrA12{:,gg}))];
                    val = [zscore(tiedrank(val1{:,gg}));   zscore(tiedrank(val2{:,gg}))];
                elseif  strcmp(typeCorr, 'Pearson')
                    mfr = [zscore(mfrA11{:,gg}); zscore(mfrA12{:,gg})];
                    val = [zscore(val1{:,gg});   zscore(val2{:,gg})];
                end
                sess= [ones(size(val1{gg},1),1) * -0.5 ; ones(size(val2{gg},1),1) * 0.5];
                try
                    stats = regstats(mfr,[val, sess],'interaction',{'tstat'});
                    a(1,gg) = stats.tstat.beta(2);
                    c(1,gg) = stats.tstat.pval(2);
                    a(2,gg) = stats.tstat.beta(4);
                    c(2,gg) = stats.tstat.pval(4);
                catch  ME
                    if strcmp(ME.identifier, 'MATLAB:catenate:dimensionMismatch') || strcmp(ME.identifier, 'stats:regstats:InputSizeMismatch')
                        rethrow(ME)
                    else
                        a(1,gg) = nan;
                        c(1,gg) = nan;
                        a(2,gg) = nan;
                        c(2,gg) = nan;
                    end
                end
            end
        end
        
        
        
    elseif strcmp(typeCorr,'binSpearman')  || strcmp(typeCorr,'binPearson')  || strcmp(typeCorr,'binKendall')
        
        
        if size(val1,2)==1
            val1=repmat(val1, 1, size(mfrA11,2));
            val2=repmat(val2, 1, size(mfrA11,2));
        end
        
        
        for gg=1:size(mfrA11,2)
            
            if strcmp(typeCorr, 'binSpearman')
                mfr = [zscore(tiedrank(mfrA11(:,gg))); zscore(tiedrank(mfrA12(:,gg)))];
                val = [zscore(tiedrank(val1(:,gg)));   zscore(tiedrank(val2(:,gg)))];
            elseif  strcmp(typeCorr, 'binPearson')
                mfr = [zscore(mfrA11(:,gg)); zscore(mfrA12(:,gg))];
                val = [zscore(val1(:,gg));   zscore(val2(:,gg))];
            end
            sess= [ones(size(val1,1),1) * -0.5 ; ones(size(val2,1),1) * 0.5];
            
            try
                stats = regstats(mfr,[val, sess],'interaction',{'tstat'});
                a(1,gg) = stats.tstat.beta(2);
                c(1,gg) = stats.tstat.pval(2);
                a(2,gg) = stats.tstat.beta(4);
                c(2,gg) = stats.tstat.pval(4);
            catch
                if strcmp(ME.identifier, 'MATLAB:catenate:dimensionMismatch') || strcmp(ME.identifier, 'stats:regstats:InputSizeMismatch')
                    rethrow(ME)
                else
                    a(1,gg) = nan;
                    c(1,gg) = nan;
                    a(2,gg) = nan;
                    c(2,gg) = nan;
                end
            end
        end
        
    end
    
    
elseif type == 4
    
    global KNNparameter
    
    if nargin >7 && convtaste==1
        
        
        
        for k=1: size(val1,2)
            ind11=find(flav11(:,k)==1);
            ind12=find(flav12(:,k)==1);
            
            val11(:,k)= val1(ind11,k);
            val21(:,k)= val2(ind12,k);
            mfrA111(:,k)= mfrA11(ind11,k);
            mfrA121(:,k)= mfrA12(ind12,k);
        end
        
        val1   = val11;
        val2   = val21;
        mfrA11 = mfrA111;
        mfrA12 = mfrA121;
        
    end
    
    
    if nargin <= 7 || convtaste==1
        
        if size(val2,2)==1
            for hh = 1:size(mfrA11,2)
                %                 [a(1,hh)] = mutInfo([mfrA11(:,hh)], [val1]);
                %                 [a(2,hh)] = mutInfo([mfrA12(:,hh)], [val2]);
                
                
                % first imput discrete, second continuous (mfr)
                [a(1,hh)] = discrete_continuous_info_fast([val1], [mfrA11(:,hh)], KNNparameter, 2);
                [a(2,hh)] = discrete_continuous_info_fast([val2], [mfrA12(:,hh)], KNNparameter, 2);
            end
            
        else
            
            for hh = 1:size(mfrA11,2)
                %                 [a(1,hh)]= mutInfo([mfrA11(:,hh)], [val1(:,hh)]);
                %                 [a(2,hh)]= mutInfo([mfrA12(:,hh)], [val2(:,hh)]);
                
                [a(1,hh)] = discrete_continuous_info_fast([val1(:,hh)], [mfrA11(:,hh)], KNNparameter, 2);
                [a(2,hh)] = discrete_continuous_info_fast([val2(:,hh)], [mfrA12(:,hh)], KNNparameter, 2);
            end
            
        end
        
        c = nan(size(a));
        
    end
    
    
    if nargin >7 && convtaste==2
        
        
        
        for k=1: size(val1,2)
            ind11=find(flav11(:,k)==1);
            ind12=find(flav12(:,k)==1);
            
            val11{:,k}= val1(ind11,k);
            val21{:,k}= val2(ind12,k);
            mfrA111{:,k}= mfrA11(ind11,k);
            mfrA121{:,k}= mfrA12(ind12,k);
        end
        
        val1   = val11;
        val2   = val21;
        mfrA11 = mfrA111;
        mfrA12 = mfrA121;
        
        for gg = 1:size(mfrA11,2)
            [a(1,gg)] = discrete_continuous_info_fast([val1{gg}], [mfrA11{gg}],  KNNparameter, 2);
            [a(2,gg)] = discrete_continuous_info_fast([val2{gg}], [mfrA12{gg}],  KNNparameter, 2);
        end
        
        c = nan(size(a));
        
    end
end
