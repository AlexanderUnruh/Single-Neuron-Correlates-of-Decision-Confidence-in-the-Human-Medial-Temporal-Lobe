function plot_singlecell_Corr

tic
dbstop if error
auto =0;
set(groot,'defaultAxesLinewidth', 1)
windowsize=200; resolution = 1; stepsize=10;  stdkernel= 5;
newwinwosreg = [-5000,5000];  newwinwosregRating = [-5000,5000];
rangePlot = [-5000:stepsize:7400-1];

%% folder and subfolders to save results and load data

satiationfolder = '/media/Projects/Alex/Reclustered analysis All/';
cd(satiationfolder);      %/Volumes/MTL/CS4 output')   %orginal form michael '/Volumes/Barnacle/Data/Satiation_data/Satiation_session_data'



x_ax_lim = {[-400,600],[-500,1500],[-500,2500], [-4600,7000]};


cell2save = 1;

if auto == 1
    autoregresive_corr = 1;
    name= 'auto';
elseif auto == 0
    autoregresive_corr = 0;
    name= '';
else
    error('decide what you want, man')
end


ana_title= ['plots_new_RT_AP_TicksCorrs_paper3', name];

mkdir([satiationfolder, ana_title]);
folder_to_save = [satiationfolder, ana_title];


%% define colors for plot

% blue to red
% possible_colors= [0   0  0.85;...
%     0.3 0.45 1;...
%     1 0.5,0.4;...
%     0.99   0   0 ];



% green in 75
possible_colors = [[0 0.4470 0.7410];...
    [0.6 0.8670 0.1010];... % [0.3250 0.8000  0.2980  ];...
    [0.9290 0.5940 0.1250];...
    [0.9890 0.3250 0.0980]];

% green in 75, matlab style
possible_colors = [[0 0.4470 0.7410];...
    [ 0.4660 0.6740 0.1880];... % [0.3250 0.8000  0.2980  ];...
    [0.9290 0.5940 0.1250];...
    [0.9890 0.3250 0.0980]];

% lila in 75
possible_colors = [[0 0.4470 0.7410];...
    [0.4940 0.1840 0.5560];... % [0.3250 0.8000  0.2980  ];...
    [0.9290 0.5940 0.1250];...
    [0.9290 0.3250 0.0980]];


% % % matlab profile
% possible_colors = [[0 0.4470 0.7410];...
%     [0.8500 0.3250 0.0980];...
%     [0.9290 0.5940 0.1250];...
%     [0.4940 0.1840 0.5560]];

% lila in 75
possible_colors = [[0 0.4470 0.7410];...
    [0.9290 0.5940 0.1250];...  % [0.3250 0.8000  0.2980  ];...
    [0.9290 0.3250 0.0980];...
    [0.4940 0.1840 0.5560]];


memberEC = [24,69,112, 113, 119];
memberH = [2,4,9,  18,31, 78, 101, 121, 225, 256, 296];
memberA = [188, 194, 23, 135, 88,134];
memberPHC = [];




inputRegion2 = {'A'; 'H'; 'EC'; 'PHC'};


% inputRegion2 = {'H'};
% memberH = [121];
%
inputRegion2 = {'EC'};
memberEC = [69];

for l=1:length(inputRegion2)
    mkdir([folder_to_save, '/' char(inputRegion2{l})])
    folder_to_save2{l}=[folder_to_save, '/' char(inputRegion2{l})];
end



%% start analysis
for region=1:length(inputRegion2)
    
    inputRegion = inputRegion2{region};
    
    inputRegion=char(inputRegion);
    
    behavDir=[];
    
    load ([satiationfolder,'/Rasts matfiles/',char(inputRegion),'_rasts6.mat'], 'rasts', 'outfiles')
    
    for o=1:size(outfiles,1)
        behavDir{end+1}= outfiles{o,2};
    end
    
    
    currentdirectory = pwd;
    
    
    % some changes for parallel processing
    lastnumber=1;
    for sess=1:size(outfiles,1)
        
        unitind=ismember(cell2mat(rasts(:,2:6)),cell2mat(outfiles(sess,2:6)),'rows'); % check if the unit is part of the channel, compares the name of the patient, etc and region
        n_ses_units = sum(unitind & cell2mat(rasts(:,9)) == 7 & cell2mat(rasts(:,10)) == 3);
        
        for a=1:n_ses_units
            unitnum{sess}(a,1)= lastnumber;
            lastnumber=lastnumber+1;
        end
    end
    
    
    for sess=1:size(outfiles,1)
        
        cd(currentdirectory)
        
        
        unitind=ismember(cell2mat(rasts(:,2:6)),cell2mat(outfiles(sess,2:6)),'rows'); % check if the unit is part of the channel, compares the name of the patient, etc and region
        actualregion=char(outfiles(sess,6));
        
        PID = str2num(outfiles{1,2}); % PID
        
        
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
        
        %             lik_RatOK1 = find((ttl_n_cond(:,1) == 31) & ((ttl_n_cond(:,2) == 25)));
        %             lik_RatOK2 = find((ttl_n_cond(:,1) == 31) & ((ttl_n_cond(:,2) == 26)));
        lik_AFCOK1 = find((ttl_n_cond(:,1) == 32) & ((ttl_n_cond(:,2) == 27)));
        lik_AFCOK2 = find((ttl_n_cond(:,1) == 32) & ((ttl_n_cond(:,2) == 28)));
        
        
        baseline1= find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 11)));  % baseline extraction for analysis of change of firing rate, instead of mean firing rate
        baseline2 = find((ttl_n_cond(:,1) == 7) & ((ttl_n_cond(:,2) == 12)));
        
        [satcue,ranking,events]=load_behavioural_data(behavDir{sess});
        [RT,endResponse] = extract_reaction_times(events);
        [Red, Green, Blue, Luminance, Contrast, Hue, Saturation, ~ ] = extract_visualFeatures(outfiles{1,2});
        [Kalorien, Protein, Kohlenhydrate, Zucker, Fett, Ballaststoffe, Salz] = naehwerte(outfiles{1,2});
        
        cd ..
        
        
        if autoregresive_corr==1
            clear  val_diffAB val_2AFC val_rat
            clear val_R_Rat val_L_Rat val_B_Rat val_A_Rat val_U_Rat val_C_Rat
            
            clear   flavourRat flavourR flavourL flavourB flavourA flavourU flavourC
            
            clear  shortRT_AFC shortRT_Rat diffRT_Rat   diffRT_AFC invshortRT_AFC
            
            clear val_R_Rank val_L_Rank val_B_Rank val_A_Rank val_U_Rank val_C_Rank
            clear invval_rat_Rank invval_rat_Rat val_rat_Rank val_rat_Rat
            
            clear invval_rat_Tot  val_rat_Tot
            
            
            clear Kalorien_Rat Protein_Rat Kohlenhydrate_Rat Zucker_Rat  Fett_Rat Ballaststoffe_Rat Salz_Rat
            clear Kalorien_B Protein_B Kohlenhydrate_B Zucker_B  Fett_B Ballaststoffe_B Salz_B
            clear Kalorien_A Protein_A Kohlenhydrate_A Zucker_A  Fett_A Ballaststoffe_A Salz_A
            
            
            clear stimOrderRat Red_Rat Green_Rat Blue_Rat Contrast_Rat Luminance_Rat Hue_A Saturation_Rat
            clear stimOrderB Red_B Green_B Blue_B Contrast_B Luminance_B Hue_A Saturation_B
            clear stimOrderA Red_A Green_A Blue_A Contrast_A Luminance_A Hue_A Saturation_A
        end
        
        
        
        
        
        for jj=[3,6]
            
            
            % Output values (real answers)
            
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
                valTot = ranking{jj}(i,5);
                
                if jj==3
                    invvalRat = ranking{5}(i,4);                        % rating value (-300 to 300)
                    invvalRank = ranking{6}(i,3);                         % ranking (0 to 20)
                    invvalTot = ranking{6}(i,5);                         % ranking (0 to 20)
                elseif jj==6
                    invvalRat = ranking{2}(i,4);                        % rating value (-300 to 300)
                    invvalRank = ranking{3}(i,3);                         % ranking (0 to 20)
                    invvalTot = ranking{3}(i,5);                         % ranking (0 to 20)
                end
                
                
                val_rat_Rat(ind3,jj/3)= valRat;
                val_rat_Rank(ind3,jj/3)= valRank;
                val_rat_Tot(ind3,jj/3)= valTot;
                
                invval_rat_Rat(ind3,jj/3)= invvalRat;
                invval_rat_Rank(ind3,jj/3)= invvalRank;
                invval_rat_Tot(ind3,jj/3)= invvalTot;
                
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
            
            %             stimOrderRat2(:,jj/3)= cell2mat(satcue{jj-1}(:,2));
            %
            %             if ~isequal(stimOrderRat,stimOrderRat2)
            %                 error('whaaaatt')
            %             end
            
            
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
            
            
            
            % transform the output of the likert and convert it so that
            % the negative is A chosen and positive B chosen.
            
            lr_to_ab3(:,jj/3) = cell2mat(satcue{jj}(:,8));  % conversion factor LR to AB for whole experiment
            
            val_diffAB(:,jj/3)=   val_2AFC(:,jj/3);
            binary= ((val_diffAB(:,jj/3)>0) ==lr_to_ab3(:,jj/3));    % compute when Right chosen (positive answer). Thus, ig cnversion factor 1 and oout pos, then R chosen, and R=B, so B chosen. Whenever binary==1, B chosen.
            val_diffAB(binary,jj/3) = abs(val_diffAB(binary,jj/3));  % convert B is positive.                        if negative (L chosen) then 0 in binary, and factor 0, then L=B, thus B chosen. Thus whenever Conv factor == to find(val_2AFC>0), then B chosen.
            val_diffAB(~binary,jj/3)= -abs(val_diffAB(~binary,jj/3)); %  A in negative.                              The opposite as the previous, if pos out, R chosen, and factor 0, then L=B, A chosen, and negative as 1 ≃0
            % finally if output negative, L chosen, and factor 1, then L==A, so A chosen again, negative here as 1 ≃0.
            
            
            
            
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
                valTot = ranking{jj}(i,5);                         % ranking (0 to 20)
                
                if jj==3
                    invvalRat = ranking{5}(i,4);                        % rating value (-300 to 300)
                    invvalRank = ranking{6}(i,3);                         % ranking (0 to 20)
                    invvalTot = ranking{6}(i,5);                         % ranking (0 to 20)
                elseif jj==6
                    invvalRat = ranking{2}(i,4);                        % rating value (-300 to 300)
                    invvalRank = ranking{3}(i,3);                         % ranking (0 to 20)
                    invvalTot = ranking{3}(i,5);                         % ranking (0 to 20)
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
                
                
                val_A_Tot(indA, jj/3) = valTot;
                val_B_Tot(indB, jj/3) = valTot;
                val_L_Tot(indL, jj/3) = valTot;
                val_R_Tot(indR, jj/3) = valTot;
                val_C_Tot(indC, jj/3) = valTot;
                val_U_Tot(indU, jj/3) = valTot;
                
                
                % inverted order
                
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
                
                invval_A_Tot(indA, jj/3) = invvalTot;
                invval_B_Tot(indB, jj/3) = invvalTot;
                invval_L_Tot(indL, jj/3) = invvalTot;
                invval_R_Tot(indR, jj/3) = invvalTot;
                invval_C_Tot(indC, jj/3) = invvalTot;
                invval_U_Tot(indU, jj/3) = invvalTot;
                
                
                
                
                
                indLast= find(cell2mat(satcue{jj-1}(:,2))==i);
                
                val_A_LastRat(indA,jj/3) = satcue{jj-1}{indLast(length(indLast)),14};
                val_B_LastRat(indB,jj/3) = satcue{jj-1}{indLast(length(indLast)),14};
                val_L_LastRat(indL,jj/3) = satcue{jj-1}{indLast(length(indLast)),14};
                val_R_LastRat(indR,jj/3) = satcue{jj-1}{indLast(length(indLast)),14};
                val_C_LastRat(indC, jj/3) = satcue{jj-1}{indLast(length(indLast)),14};
                val_U_LastRat(indU, jj/3) = satcue{jj-1}{indLast(length(indLast)),14};
                
                
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
            
            longRT_Rat(:,jj/3) = endResponse{jj-1};
            longRT_AFC(:,jj/3) = endResponse{jj};
            
            shortRT_Rat(:,jj/3) = RT{jj-1};
            shortRT_AFC(:,jj/3) = RT{jj};
            
            
            otherRT_Rat(:,jj/3)= cell2mat(satcue{jj-1}(:,13))*1000;
            otherRT_AFC(:,jj/3)= cell2mat(satcue{jj}(:,13))*1000;
            
            diffRT_Rat(:,jj/3)= otherRT_Rat(:,jj/3) - shortRT_Rat(:,jj/3);
            diffRT_AFC(:,jj/3)= otherRT_AFC(:,jj/3) - shortRT_AFC(:,jj/3);
            
            
            
        end
        
        
        % inverted likert scale output
        for i= 1:20
            for hh=1:20
                if i ~=hh
                    try
                        indPair1 = find(cell2mat(satcue{3}(:,2))==i & cell2mat(satcue{3}(:,3))==hh); indPair2 =  find(cell2mat(satcue{6}(:,2))==i & cell2mat(satcue{6}(:,3))==hh);
                        
                    end
                    try
                        indPair1=find(cell2mat(satcue{3}(:,3))==i & cell2mat(satcue{3}(:,2))==hh);  indPair2 =   find(cell2mat(satcue{6}(:,3))==i & cell2mat(satcue{6}(:,2))==hh);
                    end
                    
                    invval_2AFC(indPair1, 1) = val_2AFC(indPair2,2);
                    invval_diffAB (indPair1, 1) = val_diffAB (indPair2,2);
                    invshortRT_AFC (indPair1,1) = shortRT_AFC (indPair2,2);
                    
                    
                    invval_2AFC (indPair2, 2) = val_2AFC(indPair1,1);
                    invval_diffAB (indPair2, 2) = val_diffAB (indPair1,1);
                    invshortRT_AFC (indPair2,2) = shortRT_AFC (indPair1,1);
                    
                end
            end
        end
        
        
        
        if autoregresive_corr==1
            
            val_rat=     [val_rat(1:end-1, :)];
            val_2AFC =   [val_2AFC(1:end-1, :)];
            val_diffAB = [val_diffAB(1:end-1, :)];
            
            val_C_Rat=  [val_C_Rat(1:end-1, :)];
            val_U_Rat=  [val_U_Rat(1:end-1, :)];
            val_A_Rat=  [val_A_Rat(1:end-1, :)];
            val_B_Rat=  [val_B_Rat(1:end-1, :)];
            val_L_Rat=  [val_L_Rat(1:end-1, :)];
            val_R_Rat=  [val_R_Rat(1:end-1, :)];
            
            
            flavourC=  [flavourC(1:end-1, :)];
            flavourU=  [flavourU(1:end-1, :)];
            flavourA=  [flavourA(1:end-1, :)];
            flavourB=  [flavourB(1:end-1, :)];
            flavourL=  [flavourL(1:end-1, :)];
            flavourR=  [flavourR(1:end-1, :)];
            flavourRat=[flavourRat(1:end-1, :)];
            
            
            shortRT_AFC= [shortRT_AFC(1:end-1, :)];
            shortRT_Rat= [shortRT_Rat(1:end-1, :)];
            diffRT_AFC= [diffRT_AFC(1:end-1, :)];
            diffRT_Rat= [diffRT_Rat(1:end-1, :)];
            
            
            val_C_Rank=  [val_C_Rank(1:end-1, :)];
            val_U_Rank=  [val_U_Rank(1:end-1, :)];
            val_A_Rank=  [val_A_Rank(1:end-1, :)];
            val_B_Rank=  [val_B_Rank(1:end-1, :)];
            val_L_Rank=  [val_L_Rank(1:end-1, :)];
            val_R_Rank=  [val_R_Rank(1:end-1, :)];
            
            
            val_rat_Rat= [val_rat_Rat(1:end-1, :)];
            val_rat_Rank= [val_rat_Rank(1:end-1, :)];
            invval_rat_Rat=invval_rat_Rat(1:end-1, :);
            invval_rat_Rank=invval_rat_Rank(1:end-1, :);
            
            
            val_A_Tot = val_A_Tot(1:end-1, :);
            val_B_Tot = val_B_Tot(1:end-1, :);
            val_L_Tot = val_L_Tot(1:end-1, :);
            val_R_Tot = val_R_Tot(1:end-1, :);
            val_C_Tot = val_C_Tot(1:end-1, :);
            val_U_Tot = val_U_Tot(1:end-1, :);
            
            
            
            
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
            
            invval_A_Tot =  invval_A_Tot(1:end-1, :);
            invval_B_Tot =  invval_B_Tot(1:end-1, :);
            invval_L_Tot =  invval_L_Tot(1:end-1, :);
            invval_R_Tot =  invval_R_Tot(1:end-1, :);
            invval_C_Tot =  invval_C_Tot(1:end-1, :);
            invval_U_Tot =  invval_U_Tot(1:end-1, :);
            
            
            
            invval_rat_Tot = invval_rat_Tot(1:end-1, :);
            val_rat_Tot = val_rat_Tot(1:end-1, :);
            
            
            
            %             Kalorien = Kalorien(1:end-1, :);
            %             Protein = Protein(1:end-1, :);
            %             Kohlenhydrate = Kohlenhydrate(1:end-1, :);
            %             Zucker = Zucker(1:end-1, :);
            %             Fett = Fett(1:end-1, :);
            %             Ballaststoffe = Ballaststoffe(1:end-1, :);
            %             Salz = Salz(1:end-1, :);
            
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
            
            
            %             Red = Red(1:end-1, :);
            %             Green = Green(1:end-1, :);
            %             Blue = Blue(1:end-1, :);
            %             Contrast = Contrast(1:end-1, :);
            %             Luminance = Luminance(1:end-1, :);
            %             Hue = Hue(1:end-1, :);
            %             Saturation = Saturation(1:end-1, :);
            
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
            
        end
        
        %% define reaction-times-based windows
        
        shortRT_AFC=shortRT_AFC+200;
        shortRT_Rat=shortRT_Rat+200;
        
        
        rsp_t_likertAFCs = [ -(shortRT_AFC(:,1)),          zeros(length(shortRT_AFC),1) , -shortRT_AFC(:,2),             zeros(length(shortRT_AFC),1)];
        rsp_t_likertAFCl = [  zeros(length(diffRT_AFC),1), diffRT_AFC(:,1),                zeros(length(diffRT_AFC),1),  diffRT_AFC(:,2)];
        rsp_t_likertAFCt = [ -(shortRT_AFC(:,1)),          diffRT_AFC(:,1) ,              -shortRT_AFC(:,2),             diffRT_AFC(:,2)];
        
        rsp_t_likertRats = [ -shortRT_Rat(:,1),            zeros(length(shortRT_Rat),1),  -shortRT_Rat(:,2),             zeros(length(shortRT_Rat),1)];
        rsp_t_likertRatl = [  zeros(length(diffRT_Rat),1), diffRT_Rat(:,1),                zeros(length(diffRT_Rat),1),  diffRT_Rat(:,2)];
        rsp_t_likertRatt = [ -shortRT_Rat(:,1),            diffRT_Rat(:,1),               -shortRT_Rat(:,2),             diffRT_Rat(:,2)];
        
        
        
        
        
        
        
        % pic
        do_histogram={0,1}; do_spsth={0,1}; do_gaussianKern={0,1};
        coordin= {[7,3,1],[6,12,3]}; picdim={[4,2], [4,2.8], [5,2]};
        
        
        
        binwidth={50, 100};
        
        denseaxesposition=[0.8700    0.3300    0.1000    0.1000];
        nombre={['Value Rating'], ['Value A1' ],['Value B' ],['Value A2' ],['Absolute Preference' ],['Absolute Preference A2' ],['Absolute Preference B' ],['Absolute Preference A1' ],['Chosen Value Response' ],['Chosen Value A2' ],['Chosen Value B' ], ['Chosen Value A1' ], ['Value Rating Response']};
        nombre_region=[' ', actualregion];
        figdimension=[0,0,0.5,1];
        xlabelplot= {{'Time relative to stimulus onset (s)','Time relative to response onset (s)'},{'Time relative to stimulus onset (s)','Time relative to stimulus onset (s)'} ,{'Time relative to response onset (s)','Time relative to response onset (s)'}};
        
        
        
        
        %% for A2 plots
        
        
        response_windowEnd={{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}};
        endline2= {{2, 2},{0, 0},{0, 0}};
        
        xlabelplot = {{'Time relative to stimulus A1 onset (s)','Time relative to response onset (s)'},{'Time relative to stimulus onset (s)','Time relative to stimulus onset (s)'} ,{'Time relative to response onset (s)','Time relative to response onset (s)'}};
        
        
        
        
        
        clear ratResp1 ratResp2 A11 A12 B B2 A21 A22 bsl1 bsl2 likRatP1 likRatP2 likRatR1 likRatR2 likAFCP1 likAFCP2 likAFCR1 likAFCR2 likAFCOK1 likAFCOK2 likRatOK1 likRatOK2 spikeshape
        
        for unit=1:n_ses_units
            
            
            disp(['unit number ' num2str(unitnum{sess}(unit,1)) ' out of ' num2str(size(rasts,1)/20) ' in region ' num2str(sess) ' out of ' num2str(size(outfiles,1)), ' unit ' num2str(unit) ' out of ' num2str(n_ses_units) ' in this region ' inputRegion ])
            
            
            if  (strcmp(char(inputRegion), 'EC') && ismember((unitnum{sess}(unit,1)),memberEC)) ||  (strcmp(char(inputRegion), 'H') && ismember((unitnum{sess}(unit,1)),memberH)) ||  (strcmp(char(inputRegion), 'A') && ismember((unitnum{sess}(unit,1)),memberA))
                %   (strcmp(char(inputRegion), 'PHC') && ismember((unitnum{sess}(unit,1)),memberPHC))
                
                
                
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
                bsl1 = rast_mrg_sess{baseline1(unit),12};
                bsl2 = rast_mrg_sess{baseline2(unit),12};
                
                
                
                
                
                
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
                
                
                
                
                [mfrAll1] = kernelgauss(All1, windowsize, resolution, stepsize, [-5000,+7400], stdkernel);
                [mfrAll2] = kernelgauss(All2, windowsize, resolution, stepsize, [-5000,+7400], stdkernel);
                
                
                
                [mfrRat1] = kernelgauss(ratPic1, windowsize, resolution, stepsize, newwinwosregRating, stdkernel);
                [mfrRat2] = kernelgauss(ratPic2, windowsize, resolution, stepsize, newwinwosregRating, stdkernel);
                [mfrA11] = kernelgauss(A11, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                [mfrA12] = kernelgauss(A12, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                [mfrB1] = kernelgauss(B1, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                [mfrB2] = kernelgauss(B2, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                [mfrA21] = kernelgauss(A21, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                [mfrA22] = kernelgauss(A22, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                
                [mfrRatResp1] = kernelgauss(ratLik1, windowsize, resolution, stepsize, newwinwosregRating, stdkernel);
                [mfrRatResp2] = kernelgauss(ratLik2, windowsize, resolution, stepsize, newwinwosregRating, stdkernel);
                [mfrResp1] = kernelgauss(likAFCR1, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                [mfrResp2] = kernelgauss(likAFCR2, windowsize, resolution, stepsize, newwinwosreg, stdkernel);
                
                
                steps= [newwinwosregRating(1) :stepsize: newwinwosregRating(2)-1];
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
                
                
                
                %                 %                 if autoregresive_corr==1
                %                 %
                %                 %                     mfrA11=mfrA11(2:end,:);
                %                 %                     mfrA12=mfrA12(2:end,:);
                %                 %
                %                 %                     mfrResp1=mfrResp1(2:end,:);
                %                 %                     mfrResp2=mfrResp2(2:end,:);
                %                 %
                %                 %                     mfrRat1=mfrRat1(2:end,:);
                %                 %                     mfrRat2=mfrRat2(2:end,:);
                %                 %                     mfrRatResp1=mfrRatResp1 (2:end,:);
                %                 %                     mfrRatResp2=mfrRatResp2(2:end,:);
                %                 %                 end
                
                mfrAll1 = mfrA11;
                mfrAll2 = mfrA12;
                RTq1 = quantile(shortRT_AFC(:,1), [.25, 0.5, 0.75]);
                RTq2 = quantile(shortRT_AFC(:,2), [.25, 0.5, 0.75]);
                
                RTq = [(sum(shortRT_AFC(:,1) <RTq1,2)+1)*25, (sum(shortRT_AFC(:,2) <RTq2,2)+1)*25];
                
                RTqLR = RTq .*sign(val_2AFC);
                RTqAB = RTq .*sign(val_diffAB);
                
                
                spikeshape{unit}= rast_mrg_sess{pic_rat1(unit),14};
                
                unit2use = unitnum{sess}(unit,1);
                
                
                %       [toplot1(cell2save,:), toplot2(cell2save,:), signCorr(cell2save,1), timeCorr(cell2save,1), ...
                %       mfrWholeWind1(cell2save,:), mfrWholeWind2(cell2save,:), corrConf(cell2save,1), pvalConf(cell2save,1)] = extract_Vshape_pop( mfrAll1, mfrAll2, val_2AFC, val_diffAB, All1, All2);
                
                
                %%[toplot1(cell2save,:), toplot2(cell2save,:), signCorr(cell2save,1), timeCorr(cell2save,1), corrConf(cell2save,1), pvalConf(cell2save,1)] = extract_Vshape_pop( mfrAll1, mfrAll2, val_2AFC, val_diffAB);
                
                %[toplot1RT(cell2save,:), toplot2RT(cell2save,:), signCorrRT(cell2save,1), timeCorrRT(cell2save,1), corrRT(cell2save,1), pvalRT(cell2save,1)] = extract_Vshape_pop( mfrAll1, mfrAll2, RTqLR, RTqAB);
                
                
                
                
                % if  (strcmp(char(inputRegion), 'EC') && ismember((unitnum{sess}(unit,1)),memberEC)) ||  (strcmp(char(inputRegion), 'H') && ismember((unitnum{sess}(unit,1)),memberH)) ||  (strcmp(char(inputRegion), 'A') && ismember((unitnum{sess}(unit,1)),memberA))   || (strcmp(char(inputRegion), 'PHC') && ismember((unitnum{sess}(unit,1)),memberPHC))
                
                
                clockperunit=tic;
                disp(['plotting unit number ' num2str(unitnum{sess}(unit,1)) ' out of ' num2str(size(rasts,1)/20) ' in region ' num2str(sess) ' out of ' num2str(size(outfiles,1)), ' unit ' num2str(unit) ' out of ' num2str(n_ses_units) ' in this region ' char(inputRegion)])
                
                
               %{ 
                plotVshapes('Confidence ', unit, sess, mfrAll1, mfrAll2, val_2AFC, val_diffAB, nombre_region, x_ax_lim, response_windowEnd{3}, endline2{1},  spikeshape, folder_to_save2{region}, unitnum, xlabelplot{1}, {'AP', 'AP'}, 'Confidence', possible_colors, 1)
                plotall('Confidence ', unit, sess, mfrAll1, mfrAll2, abs(val_2AFC), abs(val_2AFC), nombre_region, x_ax_lim, response_windowEnd{3}, endline2{1},  spikeshape, folder_to_save2{region}, unitnum, xlabelplot{1}, {'AP', 'AP'}, rangePlot, 1, possible_colors)
                plotRasters('Confidence ', unit, sess, mfrAll1, mfrAll2, All1, All2, abs(val_2AFC), nombre_region, x_ax_lim, response_windowEnd{3}, endline2{1},  spikeshape, folder_to_save2{region}, unitnum, xlabelplot{1}, {'AP', 'AP'}, rangePlot, 1, possible_colors)
                %}
                
                %%%plotall('Reaction Time ', unit, sess, mfrAll1, mfrAll2,  shortRT_AFC, nombre_region, x_ax_lim, response_windowEnd{3}, endline2{1}, spikeshape, folder_to_save2{region}, unitnum, xlabelplot{1}, {'RT', 'RT'}, rangePlot, 0, possible_colors)
                %%%plotall2('Mean Rating ', unit, sess, mfrAll1, mfrAll2, mfrRat1, mfrRat2, val_A_Rat,  val_B_Rat, val_rat, nombre_region, x_ax_lim, response_windowEnd{3}, endline2{1}, spikeshape, folder_to_save2{region}, unitnum, xlabelplot{1}, {'mean Rating', 'mean Rating'}, rangePlot)
                %%%plotall2('Ranking ', unit, sess, mfrAll1, mfrAll2, mfrRat1, mfrRat2, val_A_Rank,  val_B_Rank, val_rat_Rank, nombre_region, x_ax_lim, response_windowEnd{3}, endline2{1},  spikeshape, folder_to_save2{region}, unitnum, xlabelplot{1}, {'Ranking', 'Ranking'}, rangePlot)
                %%%plotall2('Mean Salience ', unit, sess, mfrAll1, mfrAll2, mfrRat1, mfrRat2,  abs(val_A_Rat),  abs(val_B_Rat), abs(val_rat), nombre_region, x_ax_lim, response_windowEnd{3}, endline2{1}, spikeshape, folder_to_save2{region}, unitnum, xlabelplot{1}, {'mean Salience', 'mean Salience'}, rangePlot)
                
                
                plotWhole('Confidence ',    unit, sess, mfrAll1, mfrAll2, All1, All2, abs(val_2AFC), abs(val_2AFC), val_2AFC, val_diffAB, nombre_region, x_ax_lim, response_windowEnd{3}, endline2{1},  spikeshape, folder_to_save2{region}, unitnum, xlabelplot{1}, {'AP', 'AP'}, rangePlot, 1, possible_colors, 1, 'Confidence')
                plotWhole('Reaction Time ', unit, sess, mfrAll1, mfrAll2, All1, All2, shortRT_AFC, RTq,  RTqLR, RTqAB, nombre_region, x_ax_lim, response_windowEnd{3}, endline2{1},  spikeshape, folder_to_save2{region}, unitnum, xlabelplot{1}, {'RT', 'RT'}, rangePlot, 1, possible_colors, 2, 'Reaction Times')
               %{
                plotRasters('Reaction Time ', unit, sess, mfrAll1, mfrAll2, All1, All2, RTq, nombre_region, x_ax_lim, response_windowEnd{3}, endline2{1},  spikeshape, folder_to_save2{region}, unitnum, xlabelplot{1}, {'RT', 'RT'}, rangePlot, 1, possible_colors)
                plotall('Reaction Time ', unit, sess, mfrAll1, mfrAll2,  shortRT_AFC, RTq, nombre_region, x_ax_lim, response_windowEnd{3}, endline2{1}, spikeshape, folder_to_save2{region}, unitnum, xlabelplot{1}, {'RT', 'RT'}, rangePlot, 2  , possible_colors)
                plotVshapes('Reaction Time ', unit, sess, mfrAll1, mfrAll2, RTqLR, RTqAB, nombre_region, x_ax_lim, response_windowEnd{3}, endline2{1},  spikeshape, folder_to_save2{region}, unitnum, xlabelplot{1}, {'RT', 'RT'}, 'Reaction Times', possible_colors, 2)
                %}
                
                disp('time for plotting this unit was: ');
                toc(clockperunit)
                close all;
            end
        end
        
        
    end
    
end
%



%plot_Vshape_pop(toplot1, toplot2, signCorr, timeCorr, folder_to_save,  mfrWholeWind1, mfrWholeWind2, corrConf, pvalConf)%

% 
% plot_Vshape_pop(toplot1, toplot2, signCorr, timeCorr, folder_to_save, possible_colors, corrConf, pvalConf)
% plot_Vshape_pop(toplot1RT, toplot2RT, signCorrRT, timeCorrRT, folder_to_save, possible_colors, corrRT, pvalRT)




function plotall(nombre, unit, sess, Stim1, Stim2,  valueplot1,  valueplot2, nombre_region, x_ax_lim, response_windowEnd, endline2,   spikeshape, folder_to_save2, unitnum, xlabelplot, letra, rangePlot, legendtype, possible_colors)

dbstop if error

asfasfdasdfasdfasdf

leftPos1 = 0.075; width= 0.40; leftPos2= leftPos1 + width + 0.07;
heightPos1 = 0.2;  height1 = 0.83- heightPos1;
heightPos2 = 0.10; height2 =  0.06;
heightPos3 = 0.02;



pos11 =  [leftPos1, heightPos1, width, height1];
pos12 =  [leftPos1, heightPos2, width, height2];
pos13 =  [leftPos1, heightPos3, width, height2];

pos21 =  [leftPos2, heightPos1, width, height1];
pos22 =  [leftPos2, heightPos2, width, height2];
pos23 =  [leftPos2, heightPos3, width, height2];

poscbh = [0.953,   heightPos3,    0.0080,   heightPos2  + height2 -  heightPos3];
pospval = [poscbh(1)-0.0065,  heightPos3 + heightPos2  + height2 + 0.02, 0.02, 0.02];

LegPos = [0.87    0.7999    0.1179    0.1078];


%denseaxesposition= [0.8800 0.6495 0.1200 0.1800];
denseaxesposition = [0.83700    0.6495    0.1200    0.1800];





try
    close all
    fig1 =  figure('color', 'w', 'units', 'normalized', 'position', [0,0, 0.6, 0.4], 'visible', 'off');
    plot_corr_only2(Stim1, Stim2, valueplot1(:,1), valueplot1(:,2), x_ax_lim{4}, response_windowEnd{2}, endline2{2},     xlabelplot{1},  letra{1},   fig1, pos11, rangePlot, 1);
    plot_corr_only2(Stim1(2:end,:), Stim2(2:end,:), valueplot1(1:end-1,1), valueplot1(1:end-1,2), x_ax_lim{4}, response_windowEnd{2}, endline2{2},  xlabelplot{1},  letra{1},   fig1, pos21, rangePlot, 0);
    
    
    denseaxes=axes('Position', denseaxesposition);
    
    density_plot(spikeshape{unit});
    set(denseaxes, 'YAxisLocation', 'right')
    
    denseaxes.FontSize = 9;
    denseaxes.FontWeight = 'bold';
    denseaxes.YAxis.FontWeight = 'bold';
    %     denseaxes.YLabel.Position= [-1.000   9.0001   -1.000];
    
    saveas([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Mean Rho '], 'fig')
    %print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Mean Rho '], '-dpng', '-r200')
    %print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Mean Rho '], '-depsc', '-r500', '-opengl')
    
   set(fig1,'Units','centimeters');
    pos = get(fig1,'Position');
    
    
  %set(fig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    %print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Mean Rho '], '-dpdf', '-r500', '-opengl')
    
    
    % reso = '-r500';
    % export_fig([folder_to_save2 '/todel', reso],   '-pdf', reso,  '-nocrop')
    % export_fig([folder_to_save2 '/todelGL', reso],     '-pdf', '-eps', '-png', reso, '-OpenGL', '-nocrop')
    % export_fig([folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Mean Rho '], '-native',  '-transparent',  '-pdf', '-eps', '-png')
catch
    warning([unitnum{sess}(unit,1), ' did not plot ' nombre ' Mean Rho ' [nombre_region] ])
    
end




close all

if legendtype > 0
    
    
    try
        
        close all
        fig1 =  figure('color', 'w', 'units', 'normalized', 'position', [0,0, 0.6, 0.4], 'visible', 'off');
        
        sb1 = plot_sPSTH_only(Stim1, Stim2, valueplot2(:,1), valueplot2(:,2), x_ax_lim{4}, response_windowEnd{2}, endline2{2},     xlabelplot{1},  letra{1},   fig1, pos11, rangePlot, 1, 0, possible_colors, legendtype);
        sb2 = plot_sPSTH_only(Stim1(2:end,:), Stim2(2:end,:), valueplot2(1:end-1,1), valueplot2(1:end-1,2), x_ax_lim{4}, response_windowEnd{2}, endline2{2},  xlabelplot{1},  letra{1},   fig1, pos21, rangePlot, 0, 0, possible_colors, legendtype);
        
        [ylimtop] = max([sb1.YLim(2) ,sb2.YLim(2)]);
        
        sb1.YLim = [0, ylimtop];
        sb2.YLim = [0, ylimtop];
        sb1.Children(9).YData = [ylimtop, ylimtop];
        sb1.Children(10).YData = [ylimtop, ylimtop];
        sb1.Children(11).YData = [ylimtop, ylimtop];
        sb2.Children(9).YData = [ylimtop, ylimtop];
        sb2.Children(10).YData = [ylimtop, ylimtop];
        sb2.Children(11).YData = [ylimtop, ylimtop];
        
        
        saveas([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Firing Rate '], 'fig')
        %    print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Firing Rate '], '-dpng', '-r200')
        %   print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Firing Rate '], '-depsc', '-r500', '-opengl')
        
        
        
        %  set(fig1,'Units','centimeters');
        %  pos = get(fig1,'Position');
        %  set(fig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
        %  print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Firing Rate '], '-dpdf', '-r500', '-opengl')
        
    catch
        warning([unitnum{sess}(unit,1), ' did not plot ' nombre ' Firing Rate plot ' [nombre_region] ])
        
    end
    

    
    


    try
        close all
        fig1 =  figure('color', 'w', 'units', 'normalized', 'position', [0,0, 0.6, 0.4], 'visible', 'off');
        
        sb1 = plot_sPSTH_only(Stim1, Stim2, valueplot2(:,1), valueplot2(:,2), x_ax_lim{4}, response_windowEnd{2}, endline2{2},     xlabelplot{1},  letra{1},   fig1, pos11, rangePlot, 1,2, possible_colors, legendtype);
        sb2 = plot_sPSTH_only(Stim1(2:end,:), Stim2(2:end,:), valueplot2(1:end-1,1), valueplot2(1:end-1,2), x_ax_lim{4}, response_windowEnd{2}, endline2{2},  xlabelplot{1},  letra{1},   fig1, pos21, rangePlot, 0,2, possible_colors, legendtype);
        
        [ylimtop] = max([sb1.YLim(2) ,sb2.YLim(2)]);
        
        sb1.YLim = [0, ylimtop];
        sb2.YLim = [0, ylimtop];
        sb1.Children(9).YData = [ylimtop, ylimtop];
        sb1.Children(10).YData = [ylimtop, ylimtop];
        sb1.Children(11).YData = [ylimtop, ylimtop];
        sb2.Children(9).YData = [ylimtop, ylimtop];
        sb2.Children(10).YData = [ylimtop, ylimtop];
        sb2.Children(11).YData = [ylimtop, ylimtop];
        
        
        denseaxes=axes('Position', denseaxesposition);
        
        density_plot(spikeshape{unit});
        set(denseaxes, 'YAxisLocation', 'right')
        
        denseaxes.FontSize = 9;
        denseaxes.FontWeight = 'bold';
        
        %
        %         export_fig([folder_to_save2 '/todelGL', reso],     '-pdf', '-eps', '-png', reso, '-OpenGL')
        %
        
        saveas([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Firing and Ticks '], 'fig')
        %print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Firing and Ticks '], '-dpng', '-r200')
        %print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Firing and Ticks '], '-dpdf', '-r500', '-opengl')
        
        %set(fig1,'Units','centimeters');
        %pos = get(fig1,'Position');
        %set(fig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
        
    catch
        warning([unitnum{sess}(unit,1), ' did not plot ' nombre ' Firing Rate plus Ticks ' [nombre_region] ])
        
    end
    
    
    
    %     try
    %         denseaxesposition=[    0.83700    0.6495+0.02    0.1200    0.1800];
    %         close all
    %         fig1 =  figure('color', 'w', 'units', 'normalized', 'position', [0,0, 0.6, 0.4], 'visible', 'off');
    %
    %         sb1 = plot_sPSTH_only(Stim1, Stim2, valueplot2(:,1), valueplot2(:,2), x_ax_lim{4}, response_windowEnd{2}, endline2{2},     xlabelplot{1},  letra{1},   fig1, pos11, rangePlot, 1,1, possible_colors, legendtype);
    %         sb2 = plot_sPSTH_only(Stim1(2:end,:), Stim2(2:end,:), valueplot2(1:end-1,1), valueplot2(1:end-1,2), x_ax_lim{4}, response_windowEnd{2}, endline2{2},  xlabelplot{1},  letra{1},   fig1, pos21, rangePlot, 0,1, possible_colors, legendtype);
    %
    %         [ylimtop] = max([sb1.YLim(2) ,sb2.YLim(2)]);
    %
    %         sb1.YLim = [0, ylimtop];
    %         sb2.YLim = [0, ylimtop];
    %         sb1.Children(9).YData = [ylimtop, ylimtop];
    %         sb1.Children(10).YData = [ylimtop, ylimtop];
    %         sb1.Children(11).YData = [ylimtop, ylimtop];
    %         sb2.Children(9).YData = [ylimtop, ylimtop];
    %         sb2.Children(10).YData = [ylimtop, ylimtop];
    %         sb2.Children(11).YData = [ylimtop, ylimtop];
    %
    %
    %         denseaxes=axes('Position', denseaxesposition);
    %
    %         density_plot(spikeshape{unit});
    %         set(denseaxes, 'YAxisLocation', 'right')
    %
    %         denseaxes.FontSize = 9;
    %         denseaxes.FontWeight = 'bold';
    %
    %
    %         print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Firing and Mean Rho '], '-dpng', '-r200')
    %         print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Firing and Mean Rho '], '-depsc', '-r500', '-opengl')
    %         set(fig1,'Units','centimeters');
    %         pos = get(fig1,'Position');
    %         set(fig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    %         print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Firing and Mean Rho '], '-dpdf', '-r500', '-opengl')
    %
    %     catch
    %         warning([unitnum{sess}(unit,1), ' did not plot ' nombre ' Firing Rate plus Mean Rho ' [nombre_region] ])
    %
    %     end
    
    
end





function plotVshapes(nombre, unit, sess, Stim1, Stim2,  valueplot1, valueplot2, nombre_region, x_ax_lim, response_windowEnd, endline2,   spikeshape, folder_to_save2, unitnum, xlabelplot, letra, xlabelin, possible_colors, legendtype)

dbstop if error
asdfasdfasdfasdfsaf

leftPos1 = 0.075; width= 0.40; leftPos2= leftPos1 + width + 0.07;
heightPos1 = 0.2 - 0.03;  height1 = 0.86- heightPos1 ;
heightPos2 = 0.10; height2 =  0.06;
heightPos3 = 0.02;



pos11 =  [leftPos1, heightPos1, width, height1];
pos12 =  [leftPos1, heightPos2, width, height2];
pos13 =  [leftPos1, heightPos3, width, height2];

pos21 =  [leftPos2, heightPos1, width, height1];
pos22 =  [leftPos2, heightPos2, width, height2];
pos23 =  [leftPos2, heightPos3, width, height2];

poscbh = [0.953,   heightPos3,    0.0080,   heightPos2  + height2 -  heightPos3];
pospval = [poscbh(1)-0.0065,  heightPos3 + heightPos2  + height2 + 0.02, 0.02, 0.02];

LegPos = [0.87    0.7999    0.1179    0.1078];

denseaxesposition=[0.83700    0.6495    0.1200    0.1800];

try
    
    close all
    fig1 =  figure('color', 'w', 'units', 'normalized', 'position', [0,0, 0.6, 0.3], 'visible', 'off');
    
    sb1 = plot_Confidence_Vshape(Stim1, Stim2, valueplot1(:,1),  valueplot1(:,2),  xlabelin,  fig1, pos11,  1, possible_colors,  legendtype);
    sb2 = plot_Confidence_Vshape(Stim1, Stim2, valueplot2(:,1),  valueplot2(:,2),  xlabelin,  fig1, pos21, 0, possible_colors, legendtype);
    
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
    
    
    asdasd1= annotation('textbox',[0.12 0.858+0.03 0.15 0.05],'string','Left chosen','Units','normalized','FontName','Helvetica','FontSize',12,'FontWeight','bold','EdgeColor','none');
    asdasd2= annotation('textbox',[0.32 0.858+0.03 0.15 0.05],'string','Right chosen','Units','normalized','FontName','Helvetica','FontSize',12,'FontWeight','bold','EdgeColor','none');
    asdasd3= annotation('textbox',[0.6 0.858+0.03 0.15 0.05],'string','A chosen','Units','normalized','FontName','Helvetica','FontSize',12,'FontWeight','bold','EdgeColor','none');
    asdasd4= annotation('textbox',[0.8 0.858+0.03 0.15 0.05],'string','B chosen','Units','normalized','FontName','Helvetica','FontSize',12,'FontWeight','bold','EdgeColor','none');
    %     figure(gcf)
    
    saveas([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', V-Shape '], 'fig')
    %print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', V-Shape '], '-dpng', '-r200')
    
    %print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', V-Shape '], '-depsc', '-r500', '-opengl')
    
    %set(fig1,'Units','centimeters');
    %pos = get(fig1,'Position');
    %set(fig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    %print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', V-Shape '], '-dpdf', '-r500', '-opengl')
    
catch
    warning([unitnum{sess}(unit,1), ' did not plot ' nombre ' V shape ' [nombre_region] ])
end








function plotall2(nombre, unit, sess, Stim1, Stim2, Response1, Response2,  valueplot1, valueplot2,  valueplot3, nombre_region, x_ax_lim, response_windowEnd, endline2, spikeshape, folder_to_save2, unitnum, xlabelplot, letra, rangePlot)

dbstop if error


leftPos1 = 0.075; width= 0.40; leftPos2= leftPos1 + width + 0.07;
heightPos1 = 0.15;  height1 = 0.83- heightPos1;
heightPos2 = 0.10; height2 =  0.06;
heightPos3 = 0.02;



pos11 =  [leftPos1, heightPos1, width, height1];
pos12 =  [leftPos1, heightPos2, width, height2];
pos13 =  [leftPos1, heightPos3, width, height2];

pos21 =  [leftPos2, heightPos1, width, height1];
pos22 =  [leftPos2, heightPos2, width, height2];
pos23 =  [leftPos2, heightPos3, width, height2];

poscbh = [0.953,   heightPos3,    0.0080,   heightPos2  + height2 -  heightPos3];
pospval = [poscbh(1)-0.0065,  heightPos3 + heightPos2  + height2 + 0.02, 0.02, 0.02];

LegPos = [0.87    0.7999    0.1179    0.1078];


denseaxesposition= [0.8800 0.6495 0.1200 0.1800];

try
    close all
    fig1 =  figure('color', 'w', 'units', 'normalized', 'position', [0,0, 0.6, 0.4], 'visible', 'off');
    
    
    plot_corr_only2(Response1, Response2, valueplot3(:,1), valueplot3(:,2), x_ax_lim{4}, response_windowEnd{2}, endline2{2},  xlabelplot{1},  letra{1},  fig1, pos11, rangePlot, 1);
    plot_corr_only(Stim1, Stim2, valueplot1(:,1), valueplot1(:,2),  valueplot2(:,1), valueplot2(:,2),  x_ax_lim{4}, response_windowEnd{2}, endline2{2},  xlabelplot{1},  letra{1},  fig1, pos21, rangePlot, 0);
    
    
    
    denseaxes=axes('Position', denseaxesposition);
    
    density_plot(spikeshape{unit});
    
    
    
    denseaxes.FontSize = 9;
    denseaxes.FontWeight = 'bold';
    denseaxes.YLabel.Position= [-1.000   9.0001   -1.000];
    
    saveas([fig1], [folder_to_save2 '/' nombre ' ' [nombre_region] ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit)], 'fig');
    %     print([fig1], [folder_to_save2 '/' nombre ' ' [nombre_region] ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit)], '-dpng', '-r200')
    %     print([fig1], [folder_to_save2 '/' nombre ' ' [nombre_region] ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit)], '-depsc', '-r500', '-opengl')
    %
    %     set(fig1,'Units','centimeters');
    %     pos = get(fig1,'Position');
    %     set(fig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    %     print([fig1], [folder_to_save2 '/' nombre ' ' [nombre_region] ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit)], '-dpdf', '-r500', '-opengl')
    
catch
    warning([unitnum{sess}(unit,1), ' did not plot ' nombre ' PreSat ' [nombre_region] ])
    disp(['warning in plotting ' unitnum{sess}(unit,1), ' did not plot ' nombre ' Pics ' [nombre_region]])
end


close all









function [toplot1, toplot2, signcorr, timeCorr, corrConf, pvalConf] = extract_Vshape_pop(Stim1, Stim2,  valueplot1, valueplot2)
%function [toplot1, toplot2, signcorr, timeCorr, mfrWholeWind1, mfrWholeWind2, corrConf, pvalConf] = extract_Vshape_pop(Stim1, Stim2,  valueplot1, valueplot2,  All1, All2)

dbstop if error
asdfasdfasdfasdfasf

alpha = sqrt(0.05*2);

% define the window of interest, after B and before next trial

arange3 = [621:1:1201];
arange4 = [39:1:500];

vecrange3 =  ([ zeros(621-1,1); ones([1202-621],1); zeros(1241-1202,1)]);
vecrange4 =  ([ zeros(39-1,1); ones([501-39],1); zeros(1241-501,1)]);





% define correlation coeff an pvalues of normal and autocorr, for both
% variables

% current trial

[c11, p11] = corr(Stim1, abs(valueplot1(:,1)), 'type', 'Spearman');
[c12, p12] = corr(Stim2, abs(valueplot1(:,2)), 'type', 'Spearman');

c1= ((c11+c12)/2);
p1 =  (p11<alpha & p12<alpha & (sign(c11) == sign(c12)));




% next trial
[c21, p21] = corr(Stim1(2:end,:), abs(valueplot1(1:end-1, 1)), 'type', 'Spearman');
[c22, p22] = corr(Stim2(2:end,:), abs(valueplot1(1:end-1, 2)), 'type', 'Spearman');

c2= ((c21+c22)/2);
p2 =  (p21<alpha & p22<alpha & (sign(c21) == sign(c22)));




% define the window of correlation, only if after B and befor enext trial
% and if significant in both times

pos1 = find(vecrange3 & p1); % current trial var1

pos2 = find(vecrange4 & p2);  % next trial var1




timeCorr = (sum(vecrange3 & p1) + sum(vecrange4 & p2)) * 10;  % total time in window of interest that the cells was significant, in ms

% compute the mean MFR

% creater one vector both sessions
mfr1 = [Stim1(1:end-1,:);Stim2(1:end-1,:)]; % current trial
mfr2 = [Stim1(2:end,:); Stim2(2:end,:)]; % next trial

val1 = [valueplot1(1:end-1,:)]; val1=val1(:); % variable 1
val2 = [valueplot2(1:end-1,:)]; val2=val2(:); % variable 2


% compute MFR for befre after and both, only during the significant ticks

mean11 = nanmean(mfr1(:, pos1),2); % current trial
mean12 = nanmean(mfr2(:, pos2),2); % next trial

meanAll1 = nanmean([mfr1(:, pos1) , mfr2(:, pos2)],2);

mean21 = nanmean(mfr1(:, pos1),2); % current trial
mean22 = nanmean(mfr2(:, pos2),2); % next trial

meanAll2 = nanmean([mfr1(:, pos1) , mfr2(:, pos2)],2);



meanWhlewindowAll = nanmean([mfr1(:, 841:1190) , mfr2(:, 151:500)],2);


gg=1;
for b = [100:-25:25, -25:-25:-100]  % Print
    
    
    toplot1(gg,:)  =  nanmean(meanAll1(val1 == b,:),1);
    
    toplotSE1(gg,:) = nanstd(meanAll1(val1 == b,:),1)/sqrt(sum(val1 == b));
    
    
    toplot2(gg,:)  =  nanmean(meanAll2(val2 == b,:),1);
    toplotSE2(gg,:) = nanstd(meanAll2(val2 == b,:),1)/sqrt(sum(val2 == b));
    
    
    mfrWholeWind1(gg,:)  =  nanmean(meanWhlewindowAll(val1 == b,:),1);
    mfrWholeWind2(gg,:)  =  nanmean(meanWhlewindowAll(val2 == b,:),1);
    
    gg=gg+1;
    
end



% plot_corr_only2(Response1, Response2, valueplot3(:,1), valueplot3(:,2), x_ax_lim{4}, response_windowEnd{2}, endline2{2},  xlabelplot{1},  letra{1},  fig1, pos11, rangePlot, 1);
% plot_corr_only(Stim1, Stim2, valueplot1(:,1), valueplot1(:,2),  valueplot2(:,1), valueplot2(:,2),  x_ax_lim{4}, response_windowEnd{2}, endline2{2},  xlabelplot{1},  letra{1},  fig1, pos21, rangePlot, 0);


% define the mean correlation in the window of answer


signcorr = nanmean([c1(pos1); c2(pos2)],1);




% next trial
[corrConf, pvalConf] = corr(meanWhlewindowAll, abs(valueplot1(1:end-1, 1)), 'type', 'Spearman');

%{

mfr3 = [All1;All2];
val3 = [val1;val2];

clear plotMFR*
gg=1;
for b = [100:-25:25]  % Print
    
    plotMFR1 {gg,:}= [mfr1(val1 == b,:)];
    plotMFR2 {gg,:}= [mfr2(val2 == b,:)];
    plotMFR3 {gg,:}= [mfr3(val3 == b,:)];
    
    gg=gg+1;
end
%}




function plot_Vshape_pop(toplot1, toplot2, signcorr, timecorr, folder_to_save2, possible_colors,  corrConf, pvalConf)

dbstop if error
asdfasdf

leftPos1 = 0.075; width= 0.40; leftPos2= leftPos1 + width + 0.07;
heightPos1 = 0.15;  height1 = 0.83- heightPos1;
heightPos2 = 0.10; height2 =  0.06;
heightPos3 = 0.02;



pos11 =  [leftPos1, heightPos1, width, height1];
pos12 =  [leftPos1, heightPos2, width, height2];
pos13 =  [leftPos1, heightPos3, width, height2];

pos21 =  [leftPos2, heightPos1, width, height1];
pos22 =  [leftPos2, heightPos2, width, height2];
pos23 =  [leftPos2, heightPos3, width, height2];

poscbh = [0.953,   heightPos3,    0.0080,   heightPos2  + height2 -  heightPos3];
pospval = [poscbh(1)-0.0065,  heightPos3 + heightPos2  + height2 + 0.02, 0.02, 0.02];

LegPos = [0.87    0.7999    0.1179    0.1078];



axfontsize = 12;

set(0, 'currentFigure',fig1)



set(gcf, 'PaperPositionMode', 'auto');


folder2 = split(folder_to_save2{4}, '/PHC');
folder_to_save2 = folder2{1};

% possible_colors= [0   0  0.85;...
%     0.3 0.45 1;...
%     1 0.5,0.4;...
%     0.99   0   0;...
%     0.99   0   0;...
%     1 0.5,0.4;...
%       0.3 0.45 1;...
%     0   0  0.85];

possible_colors = [possible_colors; flipud(possible_colors)];
pos = signcorr>0;
neg = signcorr<0;



% do the fig
close all
fig1 =  figure('color', 'w', 'units', 'normalized', 'position', [0,0, 0.6, 0.4], 'visible', 'on');

subplothandle= subplot('Position', pos11);

hold on
set(gca,'FontSize',axfontsize,'FontWeight', 'bold');

for gg=1:8
    eb = errorbar([gg],nanmean(toplot1(pos,gg)), nanstd(toplot1(pos,gg))/sum(pos), nanstd(toplot1(pos,gg))/sum(pos)  ,'.', 'linewidth', 1.5, 'color', possible_colors(gg,:));
    
    plot([gg], nanmean(toplot1(pos,gg)), '.', 'markersize',20, 'color', possible_colors(gg,:));
end



ylabel('Firing Rate \pm SEM (spikes/s)')

axis tight
xlim([0.5, 8.5])
ylim([2.5,3.5])
line([4.5, 4.5], [subplothandle.YLim(1), subplothandle.YLim(2)], 'linewidth', 1.5, 'linestyle', '--', 'color', [0.6,0.6,0.6])


xlabel('Confidence')
set(gca, 'Xtick', [1:8], 'Xticklabel', {'100','75','50','25',  '25', '50', '75', '100'})
set(gca, 'Ytick', [2.5:0.5:3.5], 'Yticklabel', {'2.5', '3', '3.5'})



subplothandle2= subplot('Position', pos21);

hold on
set(gca,'FontSize',axfontsize,'FontWeight', 'bold');

for gg=1:8
    eb = errorbar([gg],nanmean(toplot2(pos,gg)), nanstd(toplot2(pos,gg))/sum(pos), nanstd(toplot2(pos,gg))/sum(pos)  ,'.', 'linewidth', 1.5, 'color', possible_colors(gg,:));
    
    plot([gg], nanmean(toplot2(pos,gg)), '.', 'markersize',20, 'color', possible_colors(gg,:));
end


axis tight
xlim([0.5, 8.5])
ylim([2.5,3.5])
line([4.5, 4.5], [subplothandle2.YLim(1), subplothandle2.YLim(2)], 'linewidth', 1.5, 'linestyle', '--', 'color', [0.6,0.6,0.6])



xlabel('Confidence')
set(gca, 'Xtick', [1:8], 'Xticklabel', {'100','75','50','25',  '25', '50', '75', '100'})
set(gca, 'Ytick', [2.5:0.5:3.5], 'Yticklabel', {'2.5', '3', '3.5'})


saveas([fig1], [folder_to_save '/ V-Shape Population Positive, ', num2str(sum(pos)), ' cells'], 'fig')
print([fig1], [folder_to_save '/ V-Shape Population Positive, ', num2str(sum(pos)), ' cells'], '-dpng', '-r200')
set(fig1,'Units','normalized');
pos2 = get(fig1,'Position');

print([fig1], [folder_to_save '/ V-Shape Population Positive, ', num2str(sum(pos)), ' cells'], '-depsc', '-r500', '-opengl')
set(fig1,'Units','centimeters');
pos = get(fig1,'Position');
set(fig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print([fig1], [folder_to_save '/ V-Shape Population Positive, ', num2str(sum(pos)), ' cells'], '-dpdf', '-r500', '-opengl')




asdasd1= annotation('textbox',[0.12 0.858 0.15 0.05],'string','Left chosen','Units','normalized','FontName','Helvetica','FontSize',12,'FontWeight','bold','EdgeColor','none');
asdasd2= annotation('textbox',[0.32 0.858 0.15 0.05],'string','Right chosen','Units','normalized','FontName','Helvetica','FontSize',12,'FontWeight','bold','EdgeColor','none');
asdasd3= annotation('textbox',[0.6 0.858 0.15 0.05],'string','A chosen','Units','normalized','FontName','Helvetica','FontSize',12,'FontWeight','bold','EdgeColor','none');
asdasd4= annotation('textbox',[0.8 0.858 0.15 0.05],'string','B chosen','Units','normalized','FontName','Helvetica','FontSize',12,'FontWeight','bold','EdgeColor','none');

set(fig1,'PaperPositionMode','Auto','PaperUnits','normalized','PaperSize', [pos2(3), pos2(4)])
saveas([fig1], [folder_to_save '/ V-Shape Population Positive with titles, ', num2str(sum(pos)), ' cells'], 'fig')
print([fig1], [folder_to_save '/ V-Shape Population Positive with titles, ', num2str(sum(pos)), ' cells'], '-dpng', '-r200')
print([fig1], [folder_to_save '/ V-Shape Population Positive with titles, ', num2str(sum(pos)), ' cells'], '-depsc', '-r500', '-opengl')

set(fig1,'Units','centimeters');
pos = get(fig1,'Position');
set(fig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

print([fig1], [folder_to_save '/ V-Shape Population Positive with titles, ', num2str(sum(pos)), ' cells'], '-dpdf', '-r500', '-opengl')




fig2 =  figure('color', 'w', 'units', 'normalized', 'position', [0,0, 0.6, 0.4], 'visible', 'on');

subplothandle= subplot('Position', pos11);

hold on
set(gca,'FontSize',axfontsize,'FontWeight', 'bold');

for gg=1:8
    eb = errorbar([gg],nanmean(toplot1(neg,gg)), nanstd(toplot1(neg,gg))/sum(neg), nanstd(toplot1(neg,gg))/sum(neg)  ,'.', 'linewidth', 1.5, 'color', possible_colors(gg,:));
    
    plot([gg], nanmean(toplot1(neg,gg)), '.', 'markersize',20, 'color', possible_colors(gg,:));
end



ylabel('Firing Rate \pm SEM (spikes/s)')

axis tight
xlim([0.5, 8.5])
ylim([2,3])
line([4.5, 4.5], [subplothandle.YLim(1), subplothandle.YLim(2)], 'linewidth', 1.5, 'linestyle', '--', 'color', [0.6,0.6,0.6])
set(gca, 'Xtick', [2:0.5:3], 'Xticklabel', {'2', '2.5', '3'})

xlabel('Confidence')
set(gca, 'Xtick', [1:8], 'Xticklabel', {'100','75','50','25',  '25', '50', '75', '100'})
set(gca, 'Ytick', [2:0.5:3], 'Yticklabel', {'2', '2.5', '3'})



subplothandle2= subplot('Position', pos21);

hold on
set(gca,'FontSize',axfontsize,'FontWeight', 'bold');

for gg=1:8
    eb = errorbar([gg],nanmean(toplot2(neg,gg)), nanstd(toplot2(neg,gg))/sum(neg), nanstd(toplot2(neg,gg))/sum(neg)  ,'.', 'linewidth', 1.5, 'color', possible_colors(gg,:));
    
    plot([gg], nanmean(toplot2(neg,gg)), '.', 'markersize',20, 'color', possible_colors(gg,:));
end


axis tight
xlim([0.5, 8.5])
ylim([2,3])

line([4.5, 4.5], [subplothandle2.YLim(1), subplothandle2.YLim(2)], 'linewidth', 1.5, 'linestyle', '--', 'color', [0.6,0.6,0.6])



xlabel('Confidence')
set(gca, 'Xtick', [1:8], 'Xticklabel', {'100','75','50','25',  '25', '50', '75', '100'})
set(gca, 'Ytick', [2:0.5:3], 'Yticklabel', {'2', '2.5', '3'})


saveas([fig2], [folder_to_save '/ V-Shape Population Negative, ', num2str(sum(neg)), ' cells'], 'fig')
print([fig2], [folder_to_save '/ V-Shape Population Negative, ', num2str(sum(neg)), ' cells'], '-dpng', '-r200')
set(fig2,'Units','normalized');
pos2 = get(fig2,'Position');
print([fig2], [folder_to_save '/ V-Shape Population Negative, ', num2str(sum(neg)), ' cells'], '-depsc', '-r500', '-opengl')
set(fig2,'Units','centimeters');
pos = get(fig2,'Position');
set(fig2,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print([fig2], [folder_to_save '/ V-Shape Population Negative, ', num2str(sum(neg)), ' cells'], '-dpdf', '-r500', '-opengl')

asdasd1= annotation('textbox',[0.12 0.858 0.15 0.05],'string','Left chosen','Units','normalized','FontName','Helvetica','FontSize',12,'FontWeight','bold','EdgeColor','none');
asdasd2= annotation('textbox',[0.32 0.858 0.15 0.05],'string','Right chosen','Units','normalized','FontName','Helvetica','FontSize',12,'FontWeight','bold','EdgeColor','none');
asdasd3= annotation('textbox',[0.6 0.858 0.15 0.05],'string','A chosen','Units','normalized','FontName','Helvetica','FontSize',12,'FontWeight','bold','EdgeColor','none');
asdasd4= annotation('textbox',[0.8 0.858 0.15 0.05],'string','B chosen','Units','normalized','FontName','Helvetica','FontSize',12,'FontWeight','bold','EdgeColor','none');

set(fig2,'PaperPositionMode','Auto','PaperUnits','normalized','PaperSize',[pos2(3), pos2(4)])
saveas([fig2], [folder_to_save '/ V-Shape Population Negative with titles, ', num2str(sum(neg)), ' cells'], 'fig')
print([fig2], [folder_to_save '/ V-Shape Population Negative with titles, ', num2str(sum(neg)), ' cells'], '-dpng', '-r200')
print([fig2], [folder_to_save '/ V-Shape Population Negative with titles, ', num2str(sum(neg)), ' cells'], '-depsc', '-r500', '-opengl')
set(fig2,'Units','centimeters');
pos = get(fig2,'Position');
set(fig2,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print([fig2], [folder_to_save '/ V-Shape Population Negative with titles, ', num2str(sum(neg)), ' cells'], '-dpdf', '-r500', '-opengl')








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










function plotRasters(nombre, unit, sess, Stim1, Stim2, Rast1, Rast2,  valueplot1,  nombre_region, x_ax_lim, response_windowEnd, endline2,   spikeshape, folder_to_save2, unitnum, xlabelplot, letra, rangePlot, doMorePlots, possible_colors)

dbstop if error
asfasdffasdf


leftPos1 = 0.075; width= 0.40; leftPos2= leftPos1 + width + 0.07;
heightPos1 = 0.2 - 0.17;  height1 = 0.99- heightPos1 ;
heightPos2 = 0.10; height2 =  0.06;
heightPos3 = 0.02;



pos11 =  [leftPos1, heightPos1, width, height1];
pos12 =  [leftPos1, heightPos2, width, height2];
pos13 =  [leftPos1, heightPos3, width, height2];

pos21 =  [leftPos2, heightPos1, width, height1];
pos22 =  [leftPos2, heightPos2, width, height2];
pos23 =  [leftPos2, heightPos3, width, height2];

poscbh = [0.953,   heightPos3,    0.0080,   heightPos2  + height2 -  heightPos3];
pospval = [poscbh(1)-0.0065,  heightPos3 + heightPos2  + height2 + 0.02, 0.02, 0.02];

LegPos = [0.87    0.7999    0.1179    0.1078];


%denseaxesposition= [1.88 0.6495 0.1200 0.1800];
denseaxesposition = [0.83700    0.6495    0.1200    0.1800];




    




try
    
    close all
    fig1 =  figure('color', 'w', 'units', 'normalized', 'position', [0,0, 0.6, 0.18], 'visible', 'off');
    markersize =6;
    
    sb1 = plot_Raster_only(Rast1, Rast2, valueplot1(:,1), valueplot1(:,2), x_ax_lim{4}, response_windowEnd{2}, endline2{2},     xlabelplot{1},  letra{1},   fig1, pos11, rangePlot, 1, 0, possible_colors, markersize);
    sb2 = plot_Raster_only(Rast1(2:end), Rast2(2:end), valueplot1(1:end-1,1), valueplot1(1:end-1,2), x_ax_lim{4}, response_windowEnd{2}, endline2{2},  xlabelplot{1},  letra{1},   fig1, pos21, rangePlot, 0, 0, possible_colors, markersize);
    
    [ylimtop] = max([sb1.YLim(2) ,sb2.YLim(2)]);
    
    sb1.YLim = [0, ylimtop];
    sb2.YLim = [0, ylimtop];
    sb1.Children(end-2).YData = [ylimtop, ylimtop];
    sb1.Children(end-2).YData = [ylimtop, ylimtop];
    sb1.Children(end-1).YData = [ylimtop, ylimtop];
    sb2.Children(end-1).YData = [ylimtop, ylimtop];
    sb2.Children(end).YData = [ylimtop, ylimtop];
    sb2.Children(end).YData = [ylimtop, ylimtop];
    
    
    saveas([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Raster only super small'], 'fig')
    
    %print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Raster only super small'], '-dpng', '-r200')
    %     print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Raster only '], '-depsc', '-r500', '-opengl')
    
    %set(fig1,'Units','centimeters');
    %pos = get(fig1,'Position');
    %set(fig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    
    %     export_fig([folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Raster only .pdf'], '-native', '-opengl', '-transparent',  '-pdf', '-eps', '-png')
    %    print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Raster only '], '-dpdf', '-r500', '-opengl')
    %print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ', Raster only super small'], '-dpdf', '-r500', '-opengl')
    
catch
    warning([unitnum{sess}(unit,1), ' did not plot ' nombre ' Raster only ' [nombre_region] ])
    
end






























%%
%%
%%






function plotWhole(nombre, unit, sess, Stim1, Stim2,  Rast1, Rast2,  valueplot1, valueplot2, valueplot3, valueplot4,  nombre_region, x_ax_lim, response_windowEnd, endline2,   spikeshape, folder_to_save2, unitnum, xlabelplot, letra, rangePlot, doMorePlots, possible_colors, legendtype ,  xlabelin)
% plotVshapes(nombre, unit, sess, Stim1, Stim2,  valueplot1, valueplot2, nombre_region, x_ax_lim, response_windowEnd, endline2,   spikeshape, folder_to_save2, unitnum, xlabelplot, letra, xlabelin, possible_colors, legendtype)
% plotall(nombre, unit, sess, Stim1, Stim2,  valueplot1,  valueplot2, nombre_region, x_ax_lim, response_windowEnd, endline2,   spikeshape, folder_to_save2, unitnum, xlabelplot, letra, rangePlot, legendtype, possible_colors)

dbstop if error



%% relativize in cm

figpos = [0,0, 30.4800, 27.4300];

pos2x = 16.6116; pos1x = 2.2860;
heightTicks = 0.5715; heightmfr = 7.2009; heightRaster = 4.9276;   heightVshape= 5.9150;
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
anot1pos = [3.6576 heightanot 4.5720 0.4286];
anot2pos = [9.7536 heightanot 4.5720 0.4286] ;
anot3pos = [18.2880 heightanot 4.5720 0.4286];
anot4pos =[24.3840 heightanot 4.5720 0.4286];




sizeletra = 26;
sizetrialtype = 15;

shiftletrax = 1.3;

heightanotA = posticksy + heightTicks + 0.8;
anot5pos = [pos1x-shiftletrax, heightanotA, 4.5720 0.4286];
heightanotB = posvshapey + heightVshape +1.7;
anot6pos = [pos1x-shiftletrax, heightanotB,4.5720 0.4286];

heightanotC = posticksy + heightTicks + 0.3;
anot7pos = [pos1x+width/2-1, heightanotC, 4.5720 0.4286];
heightanotN = posticksy + heightTicks + 0.3;
anot8pos = [pos2x+width/2-1, heightanotN,4.5720 0.4286];

 



close all
fig1 =  figure('color', 'w', 'units', 'centimeters', 'position', figpos, 'visible', 'off');


try

    
    
    markersize =6;
    
    sb1 = plot_Raster_only(Rast1, Rast2, valueplot1(:,1), valueplot1(:,2), x_ax_lim{4}, response_windowEnd{2}, endline2{2},     xlabelplot{1},  letra{1},   fig1, posRaster1, rangePlot, 1, 0, possible_colors, markersize);
    sb2 = plot_Raster_only(Rast1(2:end), Rast2(2:end), valueplot1(1:end-1,1), valueplot1(1:end-1,2), x_ax_lim{4}, response_windowEnd{2}, endline2{2},  xlabelplot{1},  letra{1},   fig1, posRaster2, rangePlot, 0, 0, possible_colors, markersize);
    
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
    

    
    
    
    
    
    

    
    
    
    [sb1, tickspos, leg1] = plot_sPSTH_only(Stim1, Stim2, valueplot2(:,1), valueplot2(:,2), x_ax_lim{4}, response_windowEnd{2}, endline2{2},     xlabelplot{1},  letra{1},   fig1, pos11, rangePlot, 1,2, possible_colors, legendtype);
    [sb2, tickspos2] = plot_sPSTH_only(Stim1(2:end,:), Stim2(2:end,:), valueplot2(1:end-1,1), valueplot2(1:end-1,2), x_ax_lim{4}, response_windowEnd{2}, endline2{2},  xlabelplot{1},  letra{1},   fig1, pos21, rangePlot, 0,2, possible_colors, legendtype);
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
    
    
    asdasd1= annotation('textbox',[0.12 0.858+0.03 0.15 0.05],'string','Left chosen','Units','normalized','FontName','Helvetica','FontSize',sizetrialtype,'FontWeight','bold','EdgeColor','none');
    asdasd2= annotation('textbox',[0.32 0.858+0.03 0.15 0.05],'string','Right chosen','Units','normalized','FontName','Helvetica','FontSize',sizetrialtype,'FontWeight','bold','EdgeColor','none');
    asdasd3= annotation('textbox',[0.6 0.858+0.03 0.15 0.05],'string','A chosen','Units','normalized','FontName','Helvetica','FontSize',sizetrialtype,'FontWeight','bold','EdgeColor','none');
    asdasd4= annotation('textbox',[0.8 0.858+0.03 0.15 0.05],'string','B chosen','Units','normalized','FontName','Helvetica','FontSize',sizetrialtype,'FontWeight','bold','EdgeColor','none');
    
    
    
    set(asdasd1, 'units', 'centimeters', 'position', anot1pos)
    set(asdasd2, 'units', 'centimeters', 'position', anot2pos)
    set(asdasd3, 'units', 'centimeters', 'position', anot3pos)
    set(asdasd4, 'units', 'centimeters', 'position', anot4pos)
    
    asdasd5= annotation('textbox',[0.12 0.858+0.03 0.15 0.05],'string','A','Units','normalized','FontName','Helvetica','FontSize',sizeletra,'FontWeight','bold','EdgeColor','none');
    asdasd6= annotation('textbox',[0.32 0.858+0.03 0.15 0.05],'string','B','Units','normalized','FontName','Helvetica','FontSize',sizeletra,'FontWeight','bold','EdgeColor','none');
    asdasd7= annotation('textbox',[0.12 0.858+0.03 0.15 0.05],'string','Current Trial','Units','normalized','FontName','Helvetica','FontSize',sizetrialtype,'FontWeight','bold','EdgeColor','none');
    asdasd8= annotation('textbox',[0.32 0.858+0.03 0.15 0.05],'string','Next Trial','Units','normalized','FontName','Helvetica','FontSize',sizetrialtype,'FontWeight','bold','EdgeColor','none');


    set(asdasd5, 'units', 'centimeters', 'position', anot5pos)
    set(asdasd6, 'units', 'centimeters', 'position', anot6pos)
    set(asdasd7, 'units', 'centimeters', 'position', anot7pos)
    set(asdasd8, 'units', 'centimeters', 'position', anot8pos)
    
    
    print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ',  Whole plot '], '-dpng', '-r500', '-opengl')
    print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ',  Whole plot '], '-depsc', '-r500', '-opengl')
    
    print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ',  Whole plot2 '], '-depsc2', '-r500', '-opengl')
    print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ',  Whole plot2 300'], '-depsc', '-r300', '-opengl')
    set(fig1,'Units','centimeters');
    pos = get(fig1,'Position');
    set(fig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    
    print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ',  Whole plot 300'], '-dpdf', '-r300', '-opengl')
    
    print([fig1], [folder_to_save2 '/' nombre ', unit ' num2str(unitnum{sess}(unit,1)) ' session ' num2str(sess), ' unit ' num2str(unit), ' ' [nombre_region], ',  Whole plot '], '-dpdf', '-r500', '-opengl')
    
catch
    warning([unitnum{sess}(unit,1), ' did not plot ' nombre ' V shape ' [nombre_region] ])
end

