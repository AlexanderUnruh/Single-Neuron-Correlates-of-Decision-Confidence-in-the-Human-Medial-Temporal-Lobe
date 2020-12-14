function behavioural_analysis

% % Analysis of behaviour
%

dbstop if error
folder2go='/media/Projects/Alex/Reclustered analysis All';
set(groot,'defaultAxesLinewidth', 1)
f={};
cd(folder2go);


DirName= dir;
for a = 1:length(DirName)
    if ~isempty(regexp( DirName(a).name, '^[0-9]{3}$', 'once')) % && (size(DirName(a).name,2) ==3)
        f{end+1,1} = [DirName(a).name];
    end
end


assert(length(f) == 12)
possible_colors = [[0 0.4470 0.7410];...
    [0.9290 0.5940 0.1250];...  % [0.3250 0.8000  0.2980  ];...
    [0.9290 0.3250 0.0980];...
    [0.4940 0.1840 0.5560]];


%% Consumed products

Consumedindex1 = [    9,    11,     20,     2,    14,     2,    16,   nan,     16,     4,     7,   17];
Consumedindex2 = {'001', '002',  '003', '004', '005', '006', '007', '008',  '009', '010', '011','012'};
Consumedindex3 = [nan(1,27), 9,nan(1,2)  11, 20,  2, 14, 2, nan(1,1), 16, nan(1,1), nan, 16,  nan(1,1), 4, nan(1,1), 7,nan(1,1),  17];
Consumedindex4 = {{'001',9}, {'002', 11},  {'003', 20}, {'004', 2}, {'005', 14}, {'006', 2}, {'007', 16}, ...
    {'008', nan}, {'009', 16}, {'010', 4}, {'011', 7}, {'012', 17}};

%% Extract the behavioural information

k=1;
for h=1:length(f)
        
    [satcue, ranking, events] = load_behavioural_data(f{h});
    [RT, endResponse] = extract_reaction_times(events);
    satcuecell{k}=satcue;
    RTcell{k}= RT;
    
    [Red, Green, Blue, Luminance, Contrast, Hue, Saturation, ~ ] = extract_visualFeatures(f{h});
    [Kalorien, Protein, Kohlenhydrate, Zucker, Fett, Ballaststoffe, Salz] = naehwerte(f{h});
    [change_of_mind{k}, change_of_mind2{k}] = find_changes_of_mind(events);
    
    
    
    for jj=[3,6]          
        % Output values (real answers)
        
        val_rat{k}(:,jj/3)= cell2mat(satcue{jj-1}(:,14));
        val_2AFC{k}(:,jj/3)= cell2mat(satcue{jj}(:,14));
                
        % Val Product
        RatProd{k} (:,jj/3) = ranking{jj-1}(:,4);
        RankProd{k}(:,jj/3) = ranking{jj}(:,3);
        
        if jj==3
            invRatProd{k} (:,jj/3) = ranking{6-1}(:,4);
            invRankProd{k}(:,jj/3) = ranking{6}(:,3);            
        elseif jj==6
            invRatProd{k} (:,jj/3) = ranking{3-1}(:,4);
            invRankProd{k}(:,jj/3) = ranking{3}(:,3);            
        end
        
        
        % Input values        
        % Rating paradigm        
        val_rat_Rat{k}(:,jj/3)=nan(length(satcue{jj-1}(:,14)),1);
        val_rat_Rank{k}(:,jj/3)=nan(length(satcue{jj-1}(:,14)),1);        
        
        for i=1:20            
            ind3=find(cell2mat(satcue{jj-1}(:,2))==i);
            valRat=  ranking{jj-1}(i,4); % ranking            
            valRank= ranking{jj}(i,3);
            
            if jj==3
                invvalRat = ranking{5}(i,4);                        % rating value (-300 to 300)
                invvalRank = ranking{6}(i,3);                         % ranking (0 to 20)
                indinv= find(cell2mat(satcue{5}(:,2))==i); assert(length(indinv)==3); % find the inverse presentations of the stimulus
                invRT = RT{5}(indinv,1); % RT in the other session
                
            elseif jj==6
                invvalRat = ranking{2}(i,4);                        % rating value (-300 to 300)
                invvalRank = ranking{3}(i,3);                         % ranking (0 to 20)                
                indinv= find(cell2mat(satcue{2}(:,2))==i); assert(length(indinv)==3);
                invRT = RT{2}(indinv,1);
            end
            
            
            val_rat_Rat{k}(ind3,jj/3)= valRat;
            val_rat_Rank{k}(ind3,jj/3)= valRank;
            
            invval_rat_Rat{k}(ind3,jj/3)= invvalRat;
            invval_rat_Rank{k}(ind3,jj/3)= invvalRank;
            
            stimOrderRat{k}(ind3,jj/3) = i;
            flavourRat{k}(ind3, jj/3) = i<11;
            
            % visual features
            Red_Rat{k} (ind3,jj/3) = Red(i,1);
            Green_Rat{k} (ind3,jj/3) = Green(i,1);
            Blue_Rat{k} (ind3,jj/3) = Blue(i,1);
            Luminance_Rat{k} (ind3,jj/3) = Luminance(i,1);
            Contrast_Rat{k} (ind3,jj/3) = Contrast(i,1);
            Hue_Rat{k} (ind3,jj/3) = Hue(i,1);
            Saturation_Rat{k} (ind3,jj/3) = Saturation(i,1);
            
            % nutritional  features
            
            Kalorien_Rat{k} (ind3,jj/3) = Kalorien(i,1);
            Protein_Rat{k} (ind3,jj/3) = Protein(i,1);
            Kohlenhydrate_Rat{k} (ind3,jj/3) = Kohlenhydrate(i,1);
            Zucker_Rat{k} (ind3,jj/3) = Zucker(i,1);
            Fett_Rat{k} (ind3,jj/3) = Fett(i,1);
            Ballaststoffe_Rat{k} (ind3,jj/3) = Ballaststoffe(i,1);
            Salz_Rat{k} (ind3,jj/3) = Salz(i,1);
            
            
            
            %% fichaxes recentes para a plantilla
            
            RTmean_Rat{k}(ind3, jj/3)  =  mean(RT{jj-1}(ind3,1)); % mean RT for satiation effect in first task
            invRTmean_Rat{k}(ind3, jj/3)  =  mean(invRT); % mean RT for for satiation effect in first task
            beliefConfidence{k}(i,jj/3) = std(val_rat{k}(ind3,jj/3)); % from der Martios paper, defined as std form value
            
        end
        
        
        % 2AFC paradigm
        
        val_A_Rat{k}(:,jj/3) = nan(length(satcue{jj}(:,14)),1);             % value of image on A given in previous rating trials (-300 to 300)
        val_B_Rat{k}(:,jj/3) = nan(length(satcue{jj}(:,14)),1);             % value of image on B given in previous rating trials (-300 to 300)
        val_A_Rank{k}(:,jj/3) = nan(length(satcue{jj}(:,14)),1);         % rank of image on A, calculated after all the choice trials (0 to 20)
        val_B_Rank{k}(:,jj/3) = nan(length(satcue{jj}(:,14)),1);         % rank of image on B, calculated after all the choice trials (0 to 20)
        
        val_L_Rat{k}(:,jj/3) = nan(length(satcue{jj}(:,14)),1);             % value of image on A given in previous rating trials (-300 to 300)
        val_R_Rat{k}(:,jj/3) = nan(length(satcue{jj}(:,14)),1);             % value of image on B given in previous rating trials (-300 to 300)
        val_L_Rank{k}(:,jj/3) = nan(length(satcue{jj}(:,14)),1);         % rank of image on A, calculated after all the choice trials (0 to 20)
        val_R_Rank{k}(:,jj/3) = nan(length(satcue{jj}(:,14)),1);         % rank of image on B, calculated after all the choice trials (0 to 20)
        
        val_C_Rat{k}(:,jj/3) = nan(length(satcue{jj}(:,14)),1);             % value of image on A given in previous rating trials (-300 to 300)
        val_U_Rat{k}(:,jj/3) = nan(length(satcue{jj}(:,14)),1);             % value of image on B given in previous rating trials (-300 to 300)
        val_C_Rank{k}(:,jj/3) = nan(length(satcue{jj}(:,14)),1);         % rank of image on A, calculated after all the choice trials (0 to 20)
        val_U_Rank{k}(:,jj/3) = nan(length(satcue{jj}(:,14)),1);         % rank of image on B, calculated after all the choice trials (0 to 20)
        
        
        
        
        % transform the output of the likert and convert it so that
        % the negative is A chosen and positive B chosen.
        
        lr_to_ab3(:,jj/3) = cell2mat(satcue{jj}(:,8));  % conversion factor LR to AB for whole experiment        
        val_SinedConfidenceAB{k}(:,jj/3)=   val_2AFC{k}(:,jj/3);
        binary= ((val_SinedConfidenceAB{k}(:,jj/3)>0) ==lr_to_ab3(:,jj/3));    % compute when Right chosen (positive answer). Thus, ig cnversion factor 1 and oout pos, then R chosen, and R=B, so B chosen. Whenever binary==1, B chosen.
        val_SinedConfidenceAB{k}(binary,jj/3) = abs(val_SinedConfidenceAB{k}(binary,jj/3));  % convert B is positive.                        if negative (L chosen) then 0 in binary, and factor 0, then L=B, thus B chosen. Thus whenever Conv factor == to find(val_2AFC>0), then B chosen.
        val_SinedConfidenceAB{k}(~binary,jj/3)= -abs(val_SinedConfidenceAB{k}(~binary,jj/3)); %  A in negative.                              The opposite as the previous, if pos out, R chosen, and factor 0, then L=B, A chosen, and negative as 1 ≃0
        
        
        
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
            
            val_A_Rat{k}(indA, jj/3) = valRat;
            val_B_Rat{k}(indB, jj/3) = valRat;
            val_L_Rat{k}(indL, jj/3)= valRat;
            val_R_Rat{k}(indR, jj/3)= valRat;
            val_A_Rank{k}(indA, jj/3) = valRank;
            val_B_Rank{k}(indB, jj/3) = valRank;
            val_L_Rank{k}(indL, jj/3) = valRank;
            val_R_Rank{k}(indR, jj/3) = valRank;
            
            val_C_Rat{k}(indC, jj/3) = valRat;
            val_U_Rat{k}(indU, jj/3) = valRat;
            val_C_Rank{k}(indC, jj/3) = valRank;
            val_U_Rank{k}(indU, jj/3) = valRank;
            
            
            % inverted order            
            invval_A_Rat{k}(indA, jj/3) = invvalRat;
            invval_B_Rat{k}(indB, jj/3) = invvalRat;
            invval_L_Rat{k}(indL, jj/3)=  invvalRat;
            invval_R_Rat{k}(indR, jj/3)=  invvalRat;
            invval_A_Rank{k}(indA, jj/3) = invvalRank;
            invval_B_Rank{k}(indB, jj/3) = invvalRank;
            invval_L_Rank{k}(indL, jj/3) = invvalRank;
            invval_R_Rank{k}(indR, jj/3) = invvalRank;
            
            invval_C_Rat{k}(indC, jj/3) = invvalRat;
            invval_U_Rat{k}(indU, jj/3) = invvalRat;
            invval_C_Rank{k}(indC, jj/3) = invvalRank;
            invval_U_Rank{k}(indU, jj/3) = invvalRank;
            
                       
            indLast= find(cell2mat(satcue{jj-1}(:,2))==i);

            flavourA{k}(indA,jj/3) = i<11;
            flavourB{k}(indB,jj/3) = i<11;
            flavourC{k}(indC,jj/3) = i<11;
            flavourU{k}(indU,jj/3) = i<11;
            flavourL{k}(indL,jj/3) = i<11;
            flavourR{k}(indR,jj/3) = i<11;
            
            stimOrderA{k}(indA,jj/3) = i;
            stimOrderB{k}(indB,jj/3) = i;
            stimOrderC{k}(indC,jj/3) = i;
            stimOrderU{k}(indU,jj/3) = i;
            stimOrderL{k}(indL,jj/3) = i;
            stimOrderR{k}(indR,jj/3) = i;
            
            
            % visual features
            Red_A{k} (indA,jj/3) = Red(i,1);
            Green_A{k} (indA,jj/3) = Green(i,1);
            Blue_A{k} (indA,jj/3) = Blue(i,1);
            Luminance_A{k} (indA,jj/3) = Luminance(i,1);
            Contrast_A{k} (indA,jj/3) = Contrast(i,1);
            Hue_A{k} (indA,jj/3) = Hue(i,1);
            Saturation_A{k} (indA,jj/3) = Saturation(i,1);
            
            Red_B{k} (indB,jj/3) = Red(i,1);
            Green_B{k} (indB,jj/3) = Green(i,1);
            Blue_B{k} (indB,jj/3) = Blue(i,1);
            Luminance_B{k} (indB,jj/3) = Luminance(i,1);
            Contrast_B{k} (indB,jj/3) = Contrast(i,1);
            Hue_B{k} (indB,jj/3) = Hue(i,1);
            Saturation_B{k} (indB,jj/3) = Saturation(i,1);
                        
            % nutritional features            
            Kalorien_A{k} (indA,jj/3) = Kalorien(i,1);
            Protein_A{k} (indA,jj/3) = Protein(i,1);
            Kohlenhydrate_A{k} (indA,jj/3) = Kohlenhydrate(i,1);
            Zucker_A{k} (indA,jj/3) = Zucker(i,1);
            Fett_A{k} (indA,jj/3) = Fett(i,1);
            Ballaststoffe_A{k} (indA,jj/3) = Ballaststoffe(i,1);
            Salz_A{k} (indA,jj/3) = Salz(i,1);
            
            
            Kalorien_B{k} (indB,jj/3) = Kalorien(i,1);
            Protein_B{k} (indB,jj/3) = Protein(i,1);
            Kohlenhydrate_B{k} (indB,jj/3) = Kohlenhydrate(i,1);
            Zucker_B{k} (indB,jj/3) = Zucker(i,1);
            Fett_B{k} (indB,jj/3) = Fett(i,1);
            Ballaststoffe_B{k} (indB,jj/3) = Ballaststoffe(i,1);
            Salz_B{k} (indB,jj/3) = Salz(i,1);
            
            beliefConfidenceChosen{k} (indC,jj/3) = beliefConfidence{k}(i,jj/3);
            beliefConfidenceUnchosen{k} (indU,jj/3) = beliefConfidence{k}(i,jj/3);
            beliefConfidenceA{k} (indA,jj/3) = beliefConfidence{k}(i,jj/3);
            beliefConfidenceB{k} (indB,jj/3) = beliefConfidence{k}(i,jj/3);
            beliefConfidenceR{k} (indL,jj/3) = beliefConfidence{k}(i,jj/3);
            beliefConfidenceL{k} (indR,jj/3) = beliefConfidence{k}(i,jj/3);
        end
        
        
     
        
        shortRT_Rat{k}(:,jj/3) = RT{jj-1};
        shortRT_AFC{k}(:,jj/3) = RT{jj};
        
        left_chosen{k}(:,jj/3) = double(val_2AFC{k}(:,jj/3) < 0);
        A_chosen{k}(:,jj/3)= double(val_SinedConfidenceAB{k}(:,jj/3) < 0);
        
        
        mean_Aconf{k}(:,jj/3) = mean(abs(val_SinedConfidenceAB{k}(logical(A_chosen{k}(:,jj/3)), jj/3)));
        mean_Bconf{k}(:,jj/3) = mean(abs(val_SinedConfidenceAB{k}(~logical(A_chosen{k}(:,jj/3)), jj/3)));        
        mean_Lconf{k}(:,jj/3) = mean(abs(val_2AFC{k}(logical(left_chosen{k}(:,jj/3)), jj/3)));
        mean_Rconf{k}(:,jj/3) = mean(abs(val_2AFC{k}(~logical(left_chosen{k}(:,jj/3)), jj/3)));        
        mean_A_RT{k}(:,jj/3) = mean(abs(shortRT_AFC{k}(logical(A_chosen{k}(:,jj/3)), jj/3)));
        mean_B_RT{k}(:,jj/3) = mean(abs(shortRT_AFC{k}(~logical(A_chosen{k}(:,jj/3)), jj/3)));        
        mean_L_RT{k}(:,jj/3) = mean(abs(shortRT_AFC{k}(logical(left_chosen{k}(:,jj/3)), jj/3)));
        mean_R_RT{k}(:,jj/3) = mean(abs(shortRT_AFC{k}(~logical(left_chosen{k}(:,jj/3)), jj/3)));
        
        
    end

    %% inverted likert scale output
    for i= 1:20
        for hh=1:20
            if i ~=hh
                try
                    indPair1 = find(cell2mat(satcue{3}(:,2))==i & cell2mat(satcue{3}(:,3))==hh); indPair2 =  find(cell2mat(satcue{6}(:,2))==i & cell2mat(satcue{6}(:,3))==hh);
                end
                try
                    indPair1=find(cell2mat(satcue{3}(:,3))==i & cell2mat(satcue{3}(:,2))==hh);  indPair2 =   find(cell2mat(satcue{6}(:,3))==i & cell2mat(satcue{6}(:,2))==hh);
                end
                
                invval_2AFC{k}(indPair1, 1) = val_2AFC{k}(indPair2,2);
                invval_SinedConfidenceAB{k} (indPair1, 1) = val_SinedConfidenceAB{k} (indPair2,2);
                invshortRT_AFC{k} (indPair1,1) = shortRT_AFC{k} (indPair2,2);
                
                
                invval_2AFC{k} (indPair2, 2) = val_2AFC{k}(indPair1,1);
                invval_SinedConfidenceAB{k} (indPair2, 2) = val_SinedConfidenceAB{k} (indPair1,1);
                invshortRT_AFC{k} (indPair2,2) = shortRT_AFC{k} (indPair1,1);
                
            end
        end
    end
    
    
    k=k+1;
end



% Valuation consistency of valuation reviewer 2, mean of std of valaution of each item, to check consistency of valaution during tasks')
val_consistency = cell2mat(beliefConfidence);
meanconsistency = mean(mean(val_consistency));
stdconsistency = std(mean(val_consistency));
meanconsistency_before = mean(mean(val_consistency(:,1:2:end)));
meanconsistency_after  =  mean(mean(val_consistency(:,2:2:end)));
stdconsistency_before = std(mean(val_consistency(:,1:2:end)));
stdconsistency_after  =  std(mean(val_consistency(:,2:2:end)));

mean_Aconf = cell2mat(mean_Aconf);
mean_Bconf = cell2mat(mean_Bconf);
mean_Rconf = cell2mat(mean_Rconf);
mean_Lconf = cell2mat(mean_Lconf);
mean_A_RT = cell2mat(mean_A_RT);
mean_B_RT = cell2mat(mean_B_RT);
mean_R_RT = cell2mat(mean_R_RT);
mean_L_RT = cell2mat(mean_L_RT);


mean_Aconf = cell2mat(mean_Aconf);
mean_Bconf = cell2mat(mean_Bconf);


%% analysis for reviewer with analysis of general devaluation of all itms, now supp Results
clear meanRatProd*
kk =1;
for ggg = 1:12
    meanRatProd(ggg) = mean(diff(RatProd{ggg}./3,1,2)); % first for each item in each product we compute hte difference, then we compute the mean of differences by item
    
    meanRatProdSug(ggg) = mean(diff(RatProd{ggg}(1:10,:)./3,1,2));
    meanRatProdSalt(ggg) = mean(diff(RatProd{ggg}(11:20,:)./3,1,2));
    
    
    if  isnan(Consumedindex(ggg))
        
    elseif Consumedindex(ggg) <= 10
        meanRatProdConsumed(kk) = mean(diff(RatProd{ggg}(1:10,:)./3,1,2));
        meanRatProdUnconsumed(kk) = mean(diff(RatProd{ggg}(11:20,:)./3,1,2));
        kk= kk+1;
    else
        meanRatProdUnconsumed(kk) = mean(diff(RatProd{ggg}(1:10,:)./3,1,2));
        meanRatProdConsumed(kk) = mean(diff(RatProd{ggg}(11:20,:)./3,1,2));
        kk= kk+1;
    end
end

meanValue_Change = mean(meanRatProd), stdValue_Change = std(meanRatProd), signrank(meanRatProd),
meanValue_Change_ConsumedTaste = mean(meanRatProdConsumed), stdValue_Change_ConsumedTaste=  std(meanRatProdConsumed), signrank(meanRatProdConsumed)
meanValue_Change_NonConsumedTaste = mean(meanRatProdUnconsumed), stdValue_Change_NonConsumedTaste = std(meanRatProdUnconsumed), signrank(meanRatProdUnconsumed)



% now for ranking not possible! See rebuttal letter
clear meanRatProd*
kk =1;
for ggg = 1:12
    meanRatProd(ggg) = mean(diff(RankProd{ggg},1,2)); % first for each item in each product we compute hte difference, then we compute the mean of differences by item
    
    meanRatProdSug(ggg) = mean(diff(RankProd{ggg}(1:10,:),1,2));
    meanRatProdSalt(ggg) = mean(diff(RankProd{ggg}(11:20,:),1,2));
    
    
    if  isnan(Consumedindex(ggg))
        
    elseif Consumedindex(ggg) <= 10
        meanRatProdConsumed(kk) = mean(diff(RankProd{ggg}(1:10,:),1,2));
        meanRatProdUnconsumed(kk) = mean(diff(RankProd{ggg}(11:20,:),1,2));
        kk= kk+1;
    else
        meanRatProdUnconsumed(kk) = mean(diff(RankProd{ggg}(1:10,:),1,2));
        meanRatProdConsumed(kk) = mean(diff(RankProd{ggg}(11:20,:),1,2));
        kk= kk+1;
    end
end

meanValue_Change = mean(meanRatProd), stdValue_Change = std(meanRatProd), signrank(meanRatProd),
meanValue_Change_ConsumedTaste = mean(meanRatProdConsumed), stdValue_Change_ConsumedTaste=  std(meanRatProdConsumed), signrank(meanRatProdConsumed)
meanValue_Change_NonConsumedTaste = mean(meanRatProdUnconsumed), stdValue_Change_NonConsumedTaste = std(meanRatProdUnconsumed), signrank(meanRatProdUnconsumed)




%% single item, consumed, plot and rank sum

clear BeforeRt BeforeRk AfterRt AfterRk
hh=1;
for kk = [1:7,9:12]
    BeforeRt(hh,1) = RatProd{kk}(Consumedindex(kk),1)./3;
    AfterRt (hh,1) = RatProd{kk}(Consumedindex(kk),2)./3;
    BeforeRk(hh,1) = RankProd{kk}(Consumedindex(kk),1);
    AfterRk(hh,1) = RankProd{kk}(Consumedindex(kk),2);
    hh=hh+1;
end

meanValue_Change = mean(AfterRt-BeforeRt), stdValue_Change = std(AfterRt-BeforeRt), signrank(AfterRt-BeforeRt),
meanValue_Change = mean(AfterRk-BeforeRk), stdValue_Change = std(AfterRk-BeforeRk), signrank(AfterRk-BeforeRk)





%% this secitons are for the reviewers
%% influence of choice of A or B in confidence

signrank(sum(cell2mat(A_chosen)), sum(~cell2mat(A_chosen)))
signrank(sum(cell2mat(left_chosen)), sum(~cell2mat(left_chosen)))

k=1;
for j =1:12
    Aconf1 = (abs(val_2AFC{1,j}(logical(A_chosen{1,j}(:,1)),1)));
    Aconf2 = (abs(val_2AFC{1,j}(logical(A_chosen{1,j}(:,2)),2)));
    k = k+1;
end


k=1;

for sess = 1:2
    
    Aconf = (abs(val_2AFC{1,j}(logical(A_chosen{1,j}(:,sess)),sess)));
    Bconf = (abs(val_2AFC{1,j}(~logical(A_chosen{1,j}(:,sess)),sess)));
    pvalsAB_conf(k,1) = ranksum(Aconf,Bconf, 'tail', 'both');
    mean_Aconf(k) = mean (Aconf);
    mean_Bconf(k) = mean (Bconf);
    k = k+1;
    
    
end
signrank(mean_Aconf,mean_Bconf, 'tail', 'both')


