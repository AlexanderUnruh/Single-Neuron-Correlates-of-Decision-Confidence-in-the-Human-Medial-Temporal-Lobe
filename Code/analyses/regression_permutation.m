function regression_permutation(windowtype, modeltypeIn, typeCorrIn, extraname, windowsize, stepsize, inputRegion2, autoregresive_corr, permutations, altwindow)


% regression_permutation(10, 1, 'Spearman', 'ConfidenceRT_correctedWindow', 200, 10, {'ALL2'}, 0, 1000, 1)
% regression_permutation(10, 1, 'Spearman', 'ConfidenceRT_correctedWindow', 200, 10, {'ALL2'}, 0, 1000, 2)


intial = tic;

dbstop if error

% select the varaibles to analyze and crete some cells to resstructure data
% in the following
variableName =  {'Aval', 'Asal', 'Arank', 'Bval', 'Bsal', 'Brank', 'Abspref', 'RT', ...
    'absdiff', 'Rankabsdiff', 'flavourA', 'flavourB', 'flavourRat',  'Val', 'Sal', 'Rank',  'RT_Rat'};
var2useName  = {'val_A_Rat', 'val_A_Sal', 'val_A_Rank', 'val_B_Rat', 'val_B_Sal', 'val_B_Rank', 'val_Abspref', 'RT_AFC',...
    'val_absdiff', 'val_Rankabsdiff', 'flavourA', 'flavourB', 'flavourRat',       'val_rat', 'val_sal', 'val_rat_Rank',   'RT_Rat'};
varType = {'A1', 'R'};
listRatingVarInv = {'Val', 'Sal', 'RT_Rat'};
listRatingVar = {'Val', 'Sal', 'Rank', 'flavourRat', 'RT_Rat'};
noSatVar = {'flavourRat', 'flavourA', 'flavourB'};




%% settings

parallel_on = 1;
perm_on= 1;



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
    newwinwosreg = [-2350,1200+2350];
    newwinwosregRating = [-2350,2000+2350];
    
    newwindRespAFC = [-2350,2350];
    newwindRespRat = [-2350,2350];
    
    stdkernel= 5;
    
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


if altwindow == 0
    ana_title= ['clusterwise_permutation_',num2str(modeltypeIn), '_',typeCorrIn,'_window',num2str(windowtype),'_kernel', num2str(windowsize), 'ms_step' , num2str(stepsize), 'ms_',extraname];
elseif altwindow == 1
    ana_title= ['clusterwise_permutation_',num2str(modeltypeIn), '_',typeCorrIn,'_window',num2str(windowtype),'_kernel', num2str(windowsize), 'ms_step' , num2str(stepsize), 'ms_',extraname, '_decisionLocked'];
    newwinwosreg = newwindRespAFC;
    newwinwosregRating = newwindRespRat;
elseif altwindow == 2
    ana_title= ['clusterwise_permutation_',num2str(modeltypeIn), '_',typeCorrIn,'_window',num2str(windowtype),'_kernel', num2str(windowsize), 'ms_step' , num2str(stepsize), 'ms_',extraname, '_OKLocked'];
    newwinwosreg = newwindRespAFC;
    newwinwosregRating = newwindRespRat;
elseif altwindow == 3
    ana_title= ['clusterwise_permutation_',num2str(modeltypeIn), '_',typeCorrIn,'_window',num2str(windowtype),'_kernel', num2str(windowsize), 'ms_step' , num2str(stepsize), 'ms_',extraname, '_autocorrNew'];
    newwindRespAFC = [-2350,3600+2350];
    newwinwosreg = newwindRespAFC;
    newwinwosregRating = newwindRespRat;
end


if autoregresive_corr == 1
    ana_title= [ana_title, '_autocorr'];
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





if parallel_on
    delete(gcp('nocreate'))
    parpool(20);
end


if perm_on==1 % create permutation matrices to use the same permutations in each variable
    rng(0); % for replicability, in case more varialbe needed to be added
    
    sizepermRat = 60;
    sizepermAFC = 190;
    
    if autoregresive_corr ==1
        sizepermRat = 59;
        sizepermAFC = 189;
    end
    
    indeces_permuted60 = (1:sizepermRat)';
    num_indeces_permuted = 1;
    
    
    while (num_indeces_permuted < permutations)
        % Permutation of the indixes:
        index = randperm(sizepermRat)';
        % We check if the permutation is repetead:
        if ~any(all(indeces_permuted60 == repmat(index,1,size(indeces_permuted60,2))))
            % If it is not repeated, we put it in a new column:
            indeces_permuted60 = [indeces_permuted60,index];
            num_indeces_permuted = num_indeces_permuted + 1;
        end
    end
    
    
    indeces_permuted190 = (1:sizepermAFC)';
    num_indeces_permuted = 1;
    
    while (num_indeces_permuted < permutations)
        % Permutation of the indixes:
        index = randperm(sizepermAFC)';
        % We check if the permutation is repetead:
        if ~any(all(indeces_permuted190 == repmat(index,1,size(indeces_permuted190,2))))
            % If it is not repeated, we put it in a new column:
            indeces_permuted190 = [indeces_permuted190,index];
            num_indeces_permuted = num_indeces_permuted + 1;
        end
    end  
end





%% start analysis
for hhh=1:length(inputRegion2)
        
    clear unitnum
    inputRegion=inputRegion2{hhh};    
    inputRegion=char(inputRegion);    
    behavDir=[];
    
    %     load ([var2useNamesatiationfolder,'/Rasts matfiles/',char(inputRegion),'_rasts6.mat'], 'rasts', 'outfiles')
    load ([satiationfolder,'/', namefolderrasts, '/',char(inputRegion),'_rasts6.mat'], 'rasts', 'outfiles')
    for o=1:size(outfiles,1)
        behavDir{end+1}= outfiles{o,2};
    end
    
   
    
    %% regional analysis
    
    currentdirectory = pwd;
    % some changes for parallel processing
    lastnumber=1;
    for sess=1:size(outfiles,1)
        
        %       unitind=ismember(cell2mat(rasts(:,2:6)),cell2mat(outfiles(sess,2:6)),'rows'); % check if the unit is part of the channel, compares the name of the patient, etc and region
        unitind= all(ismember(rasts(:,2:6), outfiles(sess,2:6)),2);
        %         unitind= all(ismember(rasts(:,2:6), outfiles(sess,2:6),'rows'),2)
        n_ses_units = sum(unitind & cell2mat(rasts(:,9)) == 7 & cell2mat(rasts(:,10)) == 3);
        
        
        for a=1:n_ses_units
            unitnum{sess}(a,1)= lastnumber;
            lastnumber=lastnumber+1;
        end
    end
    
    lastnumber = 2*(lastnumber-1);
    %% preallocate matrices
    
    timelength = [length([newwinwosreg(1):stepsize: newwinwosreg(2)+2400]), length([newwinwosregRating(1):stepsize: newwinwosregRating(2)])] ; %length(newwinwosreg(1):stepsize:newwinwosreg(02))
    %     timelengthRat = length([newwinwosregRating(1):stepsize: newwinwosregRating(2)-1]);
    
    if altwindow == 1 || altwindow == 2 
        timelength = [length([newwinwosreg(1):stepsize: newwinwosreg(2)]), length([newwinwosregRating(1):stepsize: newwinwosregRating(2)])] ;
    end
    
    SizePartition = 20;
    indCreateMat = 1;
    tempposnum = nan(SizePartition, 1);
    nnss = tic;
    
    if exist([folder_to_save '/regressreg_' inputRegion ,'.mat'])
        delete([folder_to_save '/regressreg_' inputRegion ,'.mat' ]);
        delete([folder_to_save '/regresspval_' inputRegion ,'.mat' ]);
        CMat = matfile([folder_to_save '/regressreg_' inputRegion ,'.mat'], 'Writable',true);
        PMat = matfile([folder_to_save '/regresspval_' inputRegion ,'.mat'], 'Writable',true);
    else
        CMat = matfile([folder_to_save '/regressreg_' inputRegion ,'.mat'], 'Writable',true);
        PMat = matfile([folder_to_save '/regresspval_' inputRegion ,'.mat'], 'Writable',true);   
    end
    for loop = 1: length(variableName)
        if ismember(variableName{loop}, listRatingVar)
            vT = 2;
        else
            vT =1;
        end
        
        eval (['CMat.corr', variableName{loop} ,'_',varType{vT},' = nan(lastnumber, timelength(vT), permutations);'] );
        eval (['PMat.pval', variableName{loop} ,'_',varType{vT},' = nan(lastnumber, timelength(vT), permutations);'] );
        eval ([' tempcorr', variableName{loop} ,'_',varType{vT},' = nan(SizePartition, timelength(vT), permutations); '] );
        eval ([' temppval', variableName{loop} ,'_',varType{vT},' = nan(SizePartition, timelength(vT), permutations); '] );
        
        
        if ~ismember(variableName{loop}, noSatVar) % only for Dynamic variable variables
            eval (['CMat.corr', variableName{loop} ,'inv_',varType{vT},' = nan(lastnumber, timelength(vT), permutations);'] );
            eval (['PMat.pval', variableName{loop} ,'inv_',varType{vT},' = nan(lastnumber, timelength(vT), permutations);'] );
            eval (['CMat.corr', variableName{loop} ,'_DirectChange_',varType{vT},' = nan(lastnumber/2, timelength(vT), permutations);'] );
            eval (['PMat.pval', variableName{loop} ,'_DirectChange_',varType{vT},' = nan(lastnumber/2, timelength(vT), permutations);'] );
            
            eval ([' tempcorr', variableName{loop} ,'inv_',varType{vT},' = nan(SizePartition, timelength(vT), permutations); '] );
            eval ([' temppval', variableName{loop} ,'inv_',varType{vT},' = nan(SizePartition, timelength(vT), permutations); '] );
            eval ([' tempcorr', variableName{loop} ,'_DirectChange_',varType{vT},' = nan(SizePartition/2, timelength(vT), permutations); '] );
            eval ([' temppval', variableName{loop} ,'_DirectChange_',varType{vT},' = nan(SizePartition/2, timelength(vT), permutations); '] );
            
            if ismember(variableName{loop}, listRatingVar)
                if ismember(variableName{loop}, listRatingVarInv)  % only for Dynamic variable and trial by trial estimated for the satiation with averaged estimates
                    eval (['CMat.corr', variableName{loop} ,'mean_',varType{vT},' = nan(lastnumber, timelength(vT), permutations);'] );
                    eval (['PMat.pval', variableName{loop} ,'mean_',varType{vT},' = nan(lastnumber, timelength(vT), permutations);'] );
                    eval (['tempcorr', variableName{loop} ,'mean_',varType{vT},' = nan(SizePartition, timelength(vT), permutations); '] );
                    eval (['temppval', variableName{loop} ,'mean_',varType{vT},' = nan(SizePartition, timelength(vT), permutations); '] );
                end
            end
        end
    end
    
    toc(nnss)
    
    
    for sess = 1:size(outfiles,1)
        
        cd(currentdirectory)
        unitind= all(ismember(rasts(:,2:6),outfiles(sess,2:6)),2);
        actualregion=char(outfiles(sess,6));
        
        PIDstr = (outfiles{sess,2});
        PID = str2num(PIDstr); % PID
        
        
    
        
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
        
        
        lik_Rat1 = find((ttl_n_cond(:,1) == 44) & ((ttl_n_cond(:,2) == 19)));
        lik_Rat2 = find((ttl_n_cond(:,1) == 44) & ((ttl_n_cond(:,2) == 20)));
        lik_AFC1 = find((ttl_n_cond(:,1) == 43) & ((ttl_n_cond(:,2) == 17)));
        lik_AFC2 = find((ttl_n_cond(:,1) == 43) & ((ttl_n_cond(:,2) == 18)));
        
        
        lik_RatR1 = find((ttl_n_cond(:,1) == 50) & ((ttl_n_cond(:,2) == 21)));
        lik_RatR2 = find((ttl_n_cond(:,1) == 50) & ((ttl_n_cond(:,2) == 22)));
        lik_AFCR1 = find((ttl_n_cond(:,1) == 51) & ((ttl_n_cond(:,2) == 23)));
        lik_AFCR2 = find((ttl_n_cond(:,1) == 51) & ((ttl_n_cond(:,2) == 24)));
        lik_RatOK1 = find((ttl_n_cond(:,1) == 31) & ((ttl_n_cond(:,2) == 25)));
        lik_RatOK2 = find((ttl_n_cond(:,1) == 31) & ((ttl_n_cond(:,2) == 26)));
        lik_AFCOK1 = find((ttl_n_cond(:,1) == 32) & ((ttl_n_cond(:,2) == 27)));
        lik_AFCOK2 = find((ttl_n_cond(:,1) == 32) & ((ttl_n_cond(:,2) == 28)));
        
        baseline1= find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 11)));  % baseline extraction for analysis of change of firing rate, instead of mean firing rate
        baseline2 = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 12)));
        
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
            
            disp(['unit number ' num2str(unitnum{sess}(unit,1)) ' out of ' num2str(size(rasts,1)/20) ' in region ' num2str(sess) ' out of ' num2str(size(outfiles,1)), ' unit ' num2str(unit) ' out of ' num2str(n_ses_units) ' in this region ' inputRegion ])
            
            
            A21  = rast_mrg_sess{pic_A21(unit),12};
            A22  = rast_mrg_sess{pic_A22(unit),12};
            A11  = rast_mrg_sess{pic_A11(unit),12};
            A12  = rast_mrg_sess{pic_A12(unit),12};
            B1  = rast_mrg_sess{pic_B1(unit),12};
            B2  = rast_mrg_sess{pic_B2(unit),12};
            likAFCR1  = rast_mrg_sess{lik_AFCR1(unit),12};
            likAFCR2  = rast_mrg_sess{lik_AFCR2(unit),12};
            ratPic1= rast_mrg_sess{pic_rat1(unit),12};
            ratPic2= rast_mrg_sess{pic_rat2(unit),12};
            ratLik1= rast_mrg_sess{lik_RatR1(unit),12};
            ratLik2= rast_mrg_sess{lik_RatR2(unit),12};
            %             bsl1 = rast_mrg_sess{baseline1(unit),12};
            %             bsl2 = rast_mrg_sess{baseline2(unit),12};
            ratOKLik1= rast_mrg_sess{lik_RatOK1(unit),12};
            ratOKLik2= rast_mrg_sess{lik_RatOK2(unit),12};
            likAFCOK1= rast_mrg_sess{lik_AFCOK1(unit),12};
            likAFCOK2= rast_mrg_sess{lik_AFCOK2(unit),12};
            
            if windowtype == 10  % gaussian kernel
                
                if altwindow == 0
                    [mfrRat1] = kernelgauss(ratPic1, windowsize, resolution, stepsize, newwinwosregRating, stdkernel);
                    [mfrRat2] = kernelgauss(ratPic2, windowsize, resolution, stepsize, newwinwosregRating, stdkernel);
                    [mfrA11] = kernelgauss(A11, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                    [mfrA12] = kernelgauss(A12, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                    
                    [mfrB1] = kernelgauss(B1, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                    [mfrB2] = kernelgauss(B2, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                    [mfrA21] = kernelgauss(A21, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                    [mfrA22] = kernelgauss(A22, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                    
       
                    
                    
                    steps= [newwinwosreg(1):stepsize: newwinwosreg(2)];
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
                    
                elseif altwindow == 1
                    [mfrResp1] = kernelgauss(likAFCR1, windowsize, resolution, stepsize, newwindRespAFC, stdkernel);
                    [mfrResp2] = kernelgauss(likAFCR2, windowsize, resolution, stepsize, newwindRespAFC, stdkernel);
                    mfrA11 = mfrResp1;
                    mfrA12 = mfrResp2;
                    [mfrRatResp1] = kernelgauss(ratLik1, windowsize, resolution, stepsize, newwindRespRat, stdkernel);
                    [mfrRatResp2] = kernelgauss(ratLik2, windowsize, resolution, stepsize, newwindRespRat, stdkernel);
                    mfrRat1 = mfrRatResp1;
                    mfrRat2 = mfrRatResp2;
                elseif altwindow == 2
                    [mfrOKResp1] = kernelgauss(likAFCOK1, windowsize, resolution, stepsize, newwindRespAFC, stdkernel);
                    [mfrOKResp2] = kernelgauss(likAFCOK2, windowsize, resolution, stepsize, newwindRespAFC, stdkernel);
                    mfrA11 = mfrOKResp1;
                    mfrA12 = mfrOKResp2;
%                     [mfrOKRatResp1] = kernelgauss(ratOKLik1, windowsize, resolution, stepsize, newwindRespRat, stdkernel);
%                     [mfrOKRatResp2] = kernelgauss(ratOKLik2, windowsize, resolution, stepsize, newwindRespRat, stdkernel);
%                     mfrRat1 = mfrOKRatResp1;
%                     mfrRat2 = mfrOKRatResp2;
                    [mfrRatResp1] = kernelgauss(ratLik1, windowsize, resolution, stepsize, newwindRespRat, stdkernel);
                    [mfrRatResp2] = kernelgauss(ratLik2, windowsize, resolution, stepsize, newwindRespRat, stdkernel);
                    mfrRat1 = mfrRatResp1;
                    mfrRat2 = mfrRatResp2;
                end
                
                % inverted MFR scale output
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
                
                
                for jj = [3,6]
                    for i = 1:20
                        ind3 = find(cell2mat(satcue{jj-1}(:,2))==i);
                        for gg = 1:3
                            mfrRat1mean(ind3(gg), :) = mean(mfrRat1(ind3,:));
                            mfrRat2mean(ind3(gg), :) = mean(mfrRat2(ind3,:));
                        end
                    end
                    
                end
                
                
                for jj = [3,6]
                    for i = 1:20
                        ind3 = find(cell2mat(satcue{jj-1}(:,2))==i);
                        if jj==3
                            indinv = find(cell2mat(satcue{5}(:,2))==i); assert(length(indinv)==3);
                            for gg = 1:3
                                invmfrRat1mean(ind3(gg), :) = mean(mfrRat2(indinv,:));
                            end
                        elseif jj==6
                            indinv = find(cell2mat(satcue{2}(:,2))==i); assert(length(indinv)==3);
                            for gg = 1:3
                                invmfrRat2mean(ind3(gg), jj/3) = mean(mfrRat1(indinv,jj/3));
                            end
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
            
            %% %% permutation test
            
            inds = [indCreateMat, indCreateMat + SizePartition/2 ];
            tempposnum(inds, 1) =  posnum;
            

            mfrChange =  mfrA11 - invmfrA11;
            mfrChangeR = mfrRat1mean - invmfrRat1mean;
            
            
            val_A_Sal = abs(val_A_Rat);
            val_B_Sal = abs(val_B_Rat);
            invval_A_Sal = abs(invval_A_Rat);
            invval_B_Sal = abs(invval_B_Rat);
            val_Abspref = abs(val_2AFC);
            invval_Abspref = abs(invval_2AFC);
            
            
            RT_AFC = shortRT_AFC;
            invRT_AFC = invshortRT_AFC;
            RT_Rat= shortRT_Rat;
            
            
            val_ratmean = val_rat_Rat;
            invval_rat = invval_rat_Rat;
            
            val_sal = abs(val_rat);
            val_salmean = abs(val_rat_Rat);
            invval_sal = abs(invval_rat_Rat);
            
            val_absdiff = abs(val_A_Rat - val_B_Rat);
            invval_absdiff = abs(invval_A_Rat - invval_B_Rat);
            val_Rankabsdiff = abs(val_A_Rank - val_B_Rank);
            invval_Rankabsdiff = abs(invval_A_Rank - invval_B_Rank);
            
            
            
            
            for loop = 1 :length(var2useName)
                if ismember(variableName{loop}, listRatingVar)
                    vT = 2;
                else
                    vT =1;
                end
                if ismember(variableName{loop}, listRatingVar)
                    if ~ismember(variableName{loop}, noSatVar)
                        if ismember(variableName{loop}, listRatingVarInv) % If trial by trial varaible use the mean as an estimate, Dynamic variables
                            eval (['var2use= ', var2useName{loop} ,'; var2use2= ', var2useName{loop} ,'mean;  invvar2use= inv', var2useName{loop}, '; dcvar2use= var2use2(:,1) - invvar2use(:,1);']);  % if trial by trial variable use mean variable for difference
                        else % if already averaged variable like Value A, Dynamic variables
                            eval (['var2use= ', var2useName{loop} ,'; var2use2= var2use;  invvar2use= inv', var2useName{loop}, '; dcvar2use= var2use2(:,1) - invvar2use(:,1);']);
                        end
                        eval (['[tempcorr', variableName{loop} ,'_',varType{vT},'(inds, :, :), temppval', variableName{loop} ,'_',varType{vT},'(inds, :, :),    tempcorr', variableName{loop} ,'mean_',varType{vT},'(inds, :, :), temppval', variableName{loop} ,'mean_',varType{vT},'(inds, :, :),         tempcorr', variableName{loop} ,'inv_',varType{vT},'(inds, :,:), temppval', variableName{loop} ,'inv_',varType{vT},'(inds, :,:), tempcorr', variableName{loop} ,'_DirectChange_',varType{vT},'(inds(1), :,:), temppval', variableName{loop} ,'_DirectChange_',varType{vT},'(inds(1), :,:)] = do_perm_Rat (indeces_permuted60,  permutations, mfrRat1, mfrRat2, var2use, var2use2, invvar2use, mfrChangeR, dcvar2use);'])
                    else % if no satiatition for Static variables
                        eval (['var2use= ', var2useName{loop} ,';']);
                        eval (['[tempcorr', variableName{loop} ,'_',varType{vT},'(inds, :, :), temppval', variableName{loop} ,'_',varType{vT},'(inds, :, :)] = do_perm_Static (indeces_permuted60,  permutations, mfrRat1, mfrRat2, var2use);'])
                    end
                else
                    if ~ismember(variableName{loop}, noSatVar) % if Dynamic variables
                        eval (['var2use= ', var2useName{loop} ,'; invvar2use= inv', var2useName{loop}, '; dcvar2use= var2use(:,1) - invvar2use(:,1);']);
                        eval (['[tempcorr', variableName{loop} ,'_',varType{vT},'(inds, :, :), temppval', variableName{loop} ,'_',varType{vT},'(inds, :, :), tempcorr', variableName{loop} ,'inv_',varType{vT},'(inds, :,:), temppval', variableName{loop} ,'inv_',varType{vT},'(inds, :,:), tempcorr', variableName{loop} ,'_DirectChange_',varType{vT},'(inds(1), :,:), temppval', variableName{loop} ,'_DirectChange_',varType{vT},'(inds(1), :,:)] = do_perm (indeces_permuted190,  permutations, mfrA11, mfrA12, var2use, invvar2use, mfrChange, dcvar2use); '])
                    else % if no satiatition for Static variables
                        eval (['var2use= ', var2useName{loop} ,';']);
                        eval (['[tempcorr', variableName{loop} ,'_',varType{vT},'(inds, :, :), temppval', variableName{loop} ,'_',varType{vT},'(inds, :, :)] = do_perm_Static (indeces_permuted190,  permutations, mfrA11, mfrA12, var2use);'])
                    end
                end
            end           
            
            if indCreateMat == SizePartition/2
                
                disp ([num2str(indCreateMat)])
                
                ind1=  [1:SizePartition/2];
                ind2 = [1+ SizePartition/2: SizePartition];
                
                for loop = 1 :length(var2useName)
                    
                    if ismember(variableName{loop}, listRatingVar)
                        vT = 2;
                    else
                        vT =1;
                    end
                    
                    
                    
                    eval ([' CMat.corr', variableName{loop} ,'_',varType{vT},'(tempposnum(ind1), :, :) = tempcorr', variableName{loop} ,'_',varType{vT},'(ind1, :,:);                            CMat.corr', variableName{loop} ,'_',varType{vT},'(tempposnum(ind2), :, :) = tempcorr', variableName{loop} ,'_',varType{vT},'(ind2, :,:);']);
                    eval ([' PMat.pval', variableName{loop} ,'_',varType{vT},'(tempposnum(ind1), :, :) = temppval', variableName{loop} ,'_',varType{vT},'(ind1, :,:);                            PMat.pval', variableName{loop} ,'_',varType{vT},'(tempposnum(ind2), :, :) = temppval', variableName{loop} ,'_',varType{vT},'(ind2, :,:);']);
                    
                    if ~ismember(variableName{loop}, noSatVar)
                        eval ([' CMat.corr', variableName{loop} ,'inv_',varType{vT},'(tempposnum(ind1), :, :) = tempcorr', variableName{loop} ,'inv_',varType{vT},'(ind1, :,:);                             CMat.corr', variableName{loop} ,'inv_',varType{vT},'(tempposnum(ind2), :, :) = tempcorr', variableName{loop} ,'inv_',varType{vT},'(ind2, :,:);']);
                        eval ([' PMat.pval', variableName{loop} ,'inv_',varType{vT},'(tempposnum(ind1), :, :) = temppval', variableName{loop} ,'inv_',varType{vT},'(ind1, :,:);                             PMat.pval', variableName{loop} ,'inv_',varType{vT},'(tempposnum(ind2), :, :) = temppval', variableName{loop} ,'inv_',varType{vT},'(ind2, :,:);']);
                        eval ([' CMat.corr', variableName{loop} ,'_DirectChange_',varType{vT},'(tempposnum(ind1), :, :) = tempcorr', variableName{loop} ,'_DirectChange_',varType{vT},'(ind1, :,:);']);
                        eval ([' PMat.pval', variableName{loop} ,'_DirectChange_',varType{vT},'(tempposnum(ind1), :, :) = temppval', variableName{loop} ,'_DirectChange_',varType{vT},'(ind1, :,:);']);
                        if ismember(variableName{loop}, listRatingVar)
                            if ismember(variableName{loop}, listRatingVarInv)
                                eval (['CMat.corr', variableName{loop} ,'mean_',varType{vT},'(tempposnum(ind1), :, :) = tempcorr', variableName{loop} ,'mean_',varType{vT},'(ind1, :,:);                        CMat.corr', variableName{loop} ,'mean_',varType{vT},'(tempposnum(ind2), :, :) = tempcorr', variableName{loop} ,'mean_',varType{vT},'(ind2, :,:);']);
                                eval (['PMat.pval', variableName{loop} ,'mean_',varType{vT},'(tempposnum(ind1), :, :) = temppval', variableName{loop} ,'mean_',varType{vT},'(ind1, :,:);                        PMat.pval', variableName{loop} ,'mean_',varType{vT},'(tempposnum(ind2), :, :) = temppval', variableName{loop} ,'mean_',varType{vT},'(ind2, :,:);']);
                                eval (['tempcorr',  variableName{loop} ,'mean_',varType{vT},'(:,:,:) = nan;']);
                                eval (['temppval',  variableName{loop} ,'mean_',varType{vT},'(:,:,:) = nan;']);
                            end
                        end
                    end
                    
                    eval (['tempcorr', variableName{loop} ,'_',varType{vT},'(:,:,:) = nan;']);
                    eval (['temppval', variableName{loop} ,'_',varType{vT},'(:,:,:) = nan;']);
                    
                    if ~ismember(variableName{loop}, noSatVar)
                        eval (['tempcorr', variableName{loop} ,'inv_',varType{vT},'(:,:,:) = nan;']);
                        eval (['temppval', variableName{loop} ,'inv_',varType{vT},'(:,:,:) = nan;']);
                        eval (['tempcorr', variableName{loop} ,'_DirectChange_',varType{vT},'(:,:,:) = nan;']);
                        eval (['temppval', variableName{loop} ,'_DirectChange_',varType{vT},'(:,:,:) = nan;']);
                    end
                end
                
                tempposnum(:,:) =  nan;
                indCreateMat = 0;
                
            end
            indCreateMat = indCreateMat+1;
            
            toc(oooo)
            
            
        end
    end
        cd (folder_to_save)
        clusterwise_permutation_analysis(variableName, var2useName, inputRegion, lastnumber)

    clear matCorrPerm matPvalPerm '*_signedRank*' '*_cluster*'
    clear pval* corr* change* MFR_m* CPD* change*
    clear -regexp ^pval ^corr ^change  ^MFR_m ^CPD ^change
    
end
cd (folder_to_save)
disp (['end of script, Elapsed time ', num2str(toc(intial)), ' seconds'])






function [a,c, mm, nn] = corr_mixed(type, mfrA11, mfrA12, val1, val2, typeCorr, convtaste, flav11, flav12)

% Define function to use

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
        
        %         [a(1,:),c(1,:)]=corr([zscore(tiedrank(mfrA11));zscore(tiedrank(mfrA12))], [zscore(tiedrank(val1));zscore(tiedrank(val2))], 'type', 'Pearson', 'rows', 'c');
        
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
            %             [a(1,gg)] = mutInfo([mfrA11{gg}], [val1{gg}]);
            %             [a(2,gg)] = mutInfo([mfrA12{gg}], [val2{gg}]);
            %
            [a(1,gg)] = discrete_continuous_info_fast([val1{gg}], [mfrA11{gg}],  KNNparameter, 2);
            [a(2,gg)] = discrete_continuous_info_fast([val2{gg}], [mfrA12{gg}],  KNNparameter, 2);
        end
        
        c = nan(size(a));
        
    end
end




function [mfr1, mfr2] = delete_low_frequencies(mfr1, mfr2)

mfr1= nan(length(mfr1),1);
mfr2= nan(length(mfr2),1);







function [histoconv2] = kernelgauss (data, kern_width, resolution, stepsize, analyseswindow, stdkernel)

sizeshift = 500; % in miliseconds

analyseswindow(1) =analyseswindow(1) - sizeshift ;
analyseswindow(2) =analyseswindow(2) + sizeshift ;

t = analyseswindow(2) - analyseswindow(1) ; % add sizeshift ms for estimation and then delete them from estimates

histo = nan(size(data,1),t/resolution);
for a = 1:size(data,1)
    %     histo(a,:) = histc(data{a,1},(-((t/2)-(resolution/2)):resolution:((t/2)-(resolution/2))));
    histo(a,:) = histc(data{a,1},((analyseswindow(1))-(resolution/2):resolution:(analyseswindow(2)-1-resolution/2)));
end
mfr=1;


kern_length = kern_width * stdkernel;
kernel = normpdf(-kern_length:resolution:kern_length, 0, kern_width);
y= histo * ((1000/resolution)/mfr);

for kkk=1:size(histo, 1)
    histoconv(kkk,:)= conv(y(kkk,:), kernel, 'same');
end

histoconv2= histoconv(:, sizeshift:stepsize:end-sizeshift); % add 500 ms fro estimation and then delete them from estimates







function [pvalCregVal_A1, corrCregVal_A1, pvalUregVal_A1, corrUregVal_A1, pvalCval_A1, corrCval_A1, pvalUval_A1, ...
    corrUval_A1, pvalmodelValCU, CPD_CregVal_A1, CPD_UregVal_A1,   CPDvar_CregVal_A1 , CPDvar_UregVal_A1] = bilinearregression (mfrA11,  mfrA12,   var2use, var2use2)

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















function [tempCorr, tempPval, invtempCorr,invtempPval, dctempCorr, dctempPval] = do_perm (indeces_permuted190, permutations, mfrA11, mfrA12, var2use, invvar2use, mfrChange, dcvar2use)



global typeCorr modeltype
modeltypelocal = modeltype;
typeCorrlocal = typeCorr;
tempCorr = nan(2, size(mfrA11,2), permutations); % pre-allocate
tempPval = nan(2, size(mfrA11,2), permutations);
invtempCorr = nan(2, size(mfrA11,2), permutations); % pre-allocate
invtempPval = nan(2, size(mfrA11,2), permutations);
dctempCorr = nan(1, size(mfrChange,2), permutations); % pre-allocate
dctempPval = nan(1, size(mfrChange,2), permutations);


parfor perm2 = 1:permutations
    
    [a,c]=corr_mixed(modeltypelocal, mfrA11(indeces_permuted190(:,perm2),:), mfrA12(indeces_permuted190(:,perm2),:), var2use(:,1), var2use(:,2),   typeCorrlocal);
    assert(size(c,2)==size(mfrA11,2));
    [tempCorr(:, :, perm2)]= a;
    [tempPval(:, :, perm2)]= c;
    
    [a,c]=corr_mixed(modeltypelocal, mfrA11(indeces_permuted190(:,perm2),:), mfrA12(indeces_permuted190(:,perm2),:), invvar2use(:,1), invvar2use(:,2),   typeCorrlocal);
    assert(size(c,2)==size(mfrA11,2));
    [invtempCorr(:, :, perm2)]= a;
    [invtempPval(:, :, perm2)]= c;
    
    [a,c]= corr(mfrChange(indeces_permuted190(:,perm2),:), dcvar2use(:,1), 'type',  typeCorrlocal, 'rows', 'c');
    assert(size(c,1)==size(mfrChange,2));
    [dctempCorr(:, :, perm2)]= a;
    [dctempPval(:, :, perm2)]= c;
    
end














function    [tempCorr, tempPval, tempCorrmean, tempPvalmean, invtempCorr,invtempPval, dctempCorr, dctempPval] = do_perm_Rat (indeces_permuted60, permutations, mfrA11, mfrA12, var2use,var2use2, invvar2use, mfrChange, dcvar2use)

global typeCorr modeltype
modeltypelocal = modeltype;
typeCorrlocal = typeCorr;

tempCorr = nan(2, size(mfrA11,2), permutations); % pre-allocate
tempPval = nan(2, size(mfrA11,2), permutations);
tempCorrmean = nan(2, size(mfrA11,2), permutations); % pre-allocate
tempPvalmean = nan(2, size(mfrA11,2), permutations);
invtempCorr = nan(2, size(mfrA11,2), permutations); % pre-allocate
invtempPval = nan(2, size(mfrA11,2), permutations);
dctempCorr = nan(1, size(mfrChange,2), permutations); % pre-allocate
dctempPval = nan(1, size(mfrChange,2), permutations);


parfor perm2 = 1:permutations
    
    [a,c]=corr_mixed(modeltypelocal, mfrA11(indeces_permuted60(:,perm2),:), mfrA12(indeces_permuted60(:,perm2),:), var2use(:,1), var2use(:,2),   typeCorrlocal);
    assert(size(c,2)==size(mfrA11,2));
    [tempCorr(:, :, perm2)]= a;
    [tempPval(:, :, perm2)]= c;
    
    
    [a,c]=corr_mixed(modeltypelocal, mfrA11(indeces_permuted60(:,perm2),:), mfrA12(indeces_permuted60(:,perm2),:), var2use2(:,1), var2use(:,2),   typeCorrlocal);
    assert(size(c,2)==size(mfrA11,2));
    [tempCorrmean(:, :, perm2)]= a;
    [tempPvalmean(:, :, perm2)]= c;
    
    
    [a,c]=corr_mixed(modeltypelocal, mfrA11(indeces_permuted60(:,perm2),:), mfrA12(indeces_permuted60(:,perm2),:), invvar2use(:,1), invvar2use(:,2),   typeCorrlocal);
    assert(size(c,2)==size(mfrA11,2));
    [invtempCorr(:, :, perm2)]= a;
    [invtempPval(:, :, perm2)]= c;
    
    [a,c]= corr(mfrChange(indeces_permuted60(:,perm2),:), dcvar2use(:,1), 'type',  typeCorr, 'rows', 'c');
    assert(size(c,1)==size(mfrChange,2));
    [dctempCorr(:, :, perm2)]= a;
    [dctempPval(:, :, perm2)]= c;
    
end
















function [clusterSumEffect, clustersize] = do_pop_stats_Tuning (var, cvar, siVar, alpha ,pos)
%% correlation of each varaible


fraction_Abspref = mean((var(1:siVar,:,:)<sqrt(2*alpha)) & (var(siVar+1:end,:,:)<sqrt(2*alpha)) & (sign(cvar(1:siVar,:,:))==sign(cvar(1+siVar:end,:,:))));

timesize = size(fraction_Abspref,2);
assert (size (fraction_Abspref,1) == 1);


clear cluster* signedRank*
for perm = 1: size (fraction_Abspref,3)
    
    cluster = ([zeros(1,1,1), fraction_Abspref(:,:,perm)>= pos, zeros(1,1,1)] );
    findStart = (find(diff(cluster,1)==1));
    findEnd = (find(diff(cluster,1)==-1)-1);
    
    
    assert (all((findEnd - findStart +1)>0))
    if ~isempty(findEnd - findStart +1)
        [clustersize(perm)] = max(findEnd - findStart +1);  % find the biggest cluster of consecutive signifcance
        clusSumEffect =[];
        for nn = 1: length(findStart)
            clusSumEffect (nn) = sum(fraction_Abspref(:,findStart(nn):findEnd(nn),perm));
        end
        [clusterSumEffect(perm)] = max(clusSumEffect);  % find the biggest cluster of consecutive significance times effect size
    else
        clustersize(perm) = 0;
        clusterSumEffect(perm) = 0;
    end
    
end















function [clusterSumEffect, clustersize] = do_pop_stats_DirectChange (dchvar,  siVar, alpha ,pos)
%% statiation direct change



fraction_Abspref_DirectChange = mean(dchvar<alpha,1);

timesize = size(fraction_Abspref_DirectChange,2);
assert (size (fraction_Abspref_DirectChange,1) == 1);

clear cluster* signedRank*
for perm = 1: size(fraction_Abspref_DirectChange,3)
    
    cluster = ([zeros(1,1,1), fraction_Abspref_DirectChange(:,:,perm)>= pos, zeros(1,1,1)] );
    findStart = find(diff(cluster,1)==1);
    findEnd = find(diff(cluster,1)==-1)-1;
    
    assert (all((findEnd - findStart +1)>0))
    if ~isempty(findEnd - findStart +1)
        [clustersize(perm)] = max(findEnd - findStart +1);  % find the biggest cluster of consecutive signifcance
        clusSumEffect= [];
        for nn = 1: length(findStart)
            clusSumEffect (nn) = sum(fraction_Abspref_DirectChange(:,findStart(nn):findEnd(nn),perm));
        end
        [clusterSumEffect(perm), ind1] = max(clusSumEffect);  % find the biggest cluster of consecutive significance times effect size
    else
        clustersize(perm) = 0;
        clusterSumEffect(perm) = 0;
    end
    
end

























function [signedRankclustersize, signedRankclusterSumEffect, signedRankSelclustersize, signedRankSelclusterSumEffect] = do_pop_stats_Satiation (var, cvar, invvar, cinvvar, siVar, alpha)

both = ((var(1:siVar,:,:)<sqrt(2*alpha)) & (var(siVar+1:end,:,:)<sqrt(2*alpha)) & (sign(cvar(1:siVar,:,:))==sign(cvar(1+siVar:end,:,:))));
bothinv = ((invvar(1:siVar,:,:)<sqrt(2*alpha)) & (invvar(siVar+1:end,:,:)<sqrt(2*alpha)) & (sign(cinvvar(1:siVar,:,:))==sign(cinvvar(1+siVar:end,:,:))));
selectVec =   [both] | [ bothinv];

fraction_Abspref =  mean(both,1);
timesize = size(fraction_Abspref,2);
assert (size(fraction_Abspref,1) == 1);



vv= tic;
clear cluster* signedRank*
parfor perm = 1: size(fraction_Abspref,3)
    
    %% satiation permutation signed rank
    
    var2plot1 = nan(siVar, timesize);
    var2plot2 = nan(siVar, timesize);
    invvar2plot1 = nan(siVar, timesize);
    invvar2plot2 = nan(siVar, timesize);
    
    
    for ll = 1: size(var,2)
        var2plot1 (selectVec(1:siVar,ll,perm),ll) = cvar(selectVec(1:siVar,ll,perm),ll);
        var2plot2 (selectVec(1:siVar,ll,perm),ll) = cvar(selectVec(1:siVar,ll, perm),ll);
        invvar2plot1 (selectVec(1:siVar,ll,perm),ll) = cinvvar(selectVec(1:siVar,ll, perm),ll);
        invvar2plot2 (selectVec(1:siVar,ll,perm),ll) = cinvvar(selectVec(1:siVar,ll, perm),ll);
    end
    
    cvarSel = [var2plot1; var2plot2];
    cinvvarSel = [invvar2plot1; invvar2plot2];
    
    meancvar = abs(cvar(1:siVar,:, perm) + cvar(1+siVar:end ,:, perm))/2;
    meancinvvar = abs(cinvvar(1:siVar,:, perm) + cinvvar(1+siVar:end ,:, perm))/2;
    meancvarSel = abs(cvarSel(1:siVar,:) + cvarSel(1+siVar:end ,:))/2;
    meancinvvarSel = abs(cinvvarSel(1:siVar,:) + cinvvarSel(1+siVar:end ,:))/2;
    
    
    
    pvalSign= nan(size(cvar,2),1); pvalSignSel= nan(size(cvar,2),1); zval_SR_sat= nan(size(cvar,2),1); zval_SR_satSel= nan(size(cvar,2),1);
    for kk= 1:size(cvar,2)
        try
            [pvalSign(kk,1), ~, stats] = signrank(meancvar(:,kk), meancinvvar(:,kk), 'tail', 'right', 'method', 'approximate');
            zval_SR_sat(kk,1) = abs(stats.zval);
        catch
            pvalSign(kk,1)= nan;
            zval_SR_sat(kk,1)= nan;
        end
    end
    
    
    for kk= 1:size(cvar,2)
        try
            [pvalSignSel(kk,1), ~, stats] = signrank(meancvarSel(:,kk), meancinvvarSel(:,kk), 'tail', 'right', 'method', 'approximate');
            zval_SR_satSel(kk,1) = abs(stats.zval);
        catch
            pvalSignSel(kk,1)= nan;
            zval_SR_satSel(kk,1)= nan;
        end
    end
    
    
    signedRankcluster05 = [zeros(1,1); pvalSign< alpha; zeros(1,1)];
    findStart = find(diff(signedRankcluster05,1)==1);
    findEnd = find(diff(signedRankcluster05,1)==-1)-1;
    
    
    assert (all((findEnd - findStart +1)>0))
    if ~isempty(findEnd - findStart +1)
        [signedRankclustersize(perm)] = max(findEnd - findStart +1);
        
        clusterSumEffect= [];
        for nn = 1: length(findStart)
            clusterSumEffect (nn) = sum(zval_SR_sat(findStart(nn):findEnd(nn)));
        end
        [signedRankclusterSumEffect(perm)] = max(clusterSumEffect);  % find the biggest cluster of consecutive significance times effect size
    else
        signedRankclustersize(perm) = 0;
        signedRankclusterSumEffect(perm) = 0;
    end
    
    
    signedRankSelcluster05 = [zeros(1,1); pvalSignSel< alpha; zeros(1,1)];
    findStart = find(diff(signedRankSelcluster05,1)==1);
    findEnd =   find(diff(signedRankSelcluster05,1)==-1)-1;
    
    assert (all((findEnd - findStart +1)>0))
    if ~isempty(findEnd - findStart +1)
        [signedRankSelclustersize(perm)] = max(findEnd - findStart +1);
        
        clusterSumEffect= [];
        for nn = 1: length(findStart)
            clusterSumEffect (nn) = sum(zval_SR_satSel(findStart(nn):findEnd(nn)));
        end
        
        [signedRankSelclusterSumEffect(perm)] = max(clusterSumEffect);  % find the biggest cluster of consecutive significance times effect size
    else
        signedRankSelclustersize(perm) = 0;
        signedRankSelclusterSumEffect(perm) = 0;
    end
end
toc(vv)




















function [tempCorr, tempPval] = do_perm_Static (indeces_permuted190, permutations, mfrA11, mfrA12, var2use)

global typeCorr modeltype
modeltypelocal = modeltype;
typeCorrlocal = typeCorr;
tempCorr = nan(2, size(mfrA11,2), permutations); % pre-allocate
tempPval = nan(2, size(mfrA11,2), permutations);

parfor perm2 = 1:permutations
    
    [a,c]=corr_mixed(modeltypelocal, mfrA11(indeces_permuted190(:,perm2),:), mfrA12(indeces_permuted190(:,perm2),:), var2use(:,1), var2use(:,2),   typeCorrlocal);
    assert(size(c,2)==size(mfrA11,2));
    [tempCorr(:, :, perm2)]= a;
    [tempPval(:, :, perm2)]= c;
end
