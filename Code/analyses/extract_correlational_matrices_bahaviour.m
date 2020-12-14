function extract_correlational_matrices_bahaviour



%% Extract_correlational_matrices_bahaviour
% This function creates correlational matrices for the behavioural
% variables


folder2go='/media/Projects/Alex/Reclustered analysis All';
f={};
cd(folder2go);


DirName= dir;
for a = 1:length(DirName)
    if regexp(DirName(a).name, '[0-9][0-9][0-9]')
        if (length(DirName(a).name)==3)
            f{end+1,1} = [DirName(a).name];
        end
    end
end
autoregresive_corr=0;



k=1;
for h=1:length(f)
    
    
    [satcue{k}, ranking, events] = load_behavioural_data(f{h});
    [RT, endResponse] = extract_reaction_times(events);
    
    [Red, Green, Blue, Luminance, Contrast, Hue, Saturation, ~ ] = extract_visualFeatures(f{h});
    [Kalorien, Protein, Kohlenhydrate, Zucker, Fett, Ballaststoffe, Salz] = naehwerte(f{h});
    
    
    
    for jj=[3,6]
        
        stimOrderRat2{k}(:,jj/3)= cell2mat(satcue{k}{jj-1}(:,2));
        val_rat{k}(:,jj/3)= cell2mat(satcue{k}{jj-1}(:,14));
        val_2AFC{k}(:,jj/3)= cell2mat(satcue{k}{jj}(:,14));
        

        % Input values
        
        % Rating paradigm
        
        val_rat_Rat{k}(:,jj/3)=nan(length(satcue{k}{jj-1}(:,14)),1);
        
        val_rat_Rank{k}(:,jj/3)=nan(length(satcue{k}{jj-1}(:,14)),1);
        
        
        
        for i=1:20
            
            ind3=find(cell2mat(satcue{k}{jj-1}(:,2))==i);
            valRat=  ranking{jj-1}(i,4); % ranking
            
            valRank= ranking{jj}(i,3);
            
            if jj==3
                invvalRat = ranking{5}(i,4);                        % rating value (-300 to 300)
                invvalRank = ranking{6}(i,3);                         % ranking (0 to 20)
            elseif jj==6
                invvalRat = ranking{2}(i,4);                        % rating value (-300 to 300)
                invvalRank = ranking{3}(i,3);                         % ranking (0 to 20)
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
            Luminance_Rat {k}(ind3,jj/3) = Luminance(i,1);
            Contrast_Rat{k} (ind3,jj/3) = Contrast(i,1);
            Hue_Rat {k}(ind3,jj/3) = Hue(i,1);
            Saturation_Rat{k} (ind3,jj/3) = Saturation(i,1);
            
            % nutritional  features
            
            Kalorien_Rat {k}(ind3,jj/3) = Kalorien(i,1);
            Protein_Rat {k}(ind3,jj/3) = Protein(i,1);
            Kohlenhydrate_Rat {k}(ind3,jj/3) = Kohlenhydrate(i,1);
            Zucker_Rat {k}(ind3,jj/3) = Zucker(i,1);
            Fett_Rat {k}(ind3,jj/3) = Fett(i,1);
            Ballaststoffe_Rat {k}(ind3,jj/3) = Ballaststoffe(i,1);
            Salz_Rat {k}(ind3,jj/3) = Salz(i,1);
            
        end
        
    end
    

    if ~isequal(stimOrderRat{k},stimOrderRat2{k})
        error('whaaaatt')
    end
    
    
    
    % 2AFC paradigm
    
    
    for jj=[3,6]
        val_A_Rat{k}(:,jj/3) = nan(length(satcue{k}{jj}(:,14)),1);             % value of image on A given in previous rating trials (-300 to 300)
        val_B_Rat{k}(:,jj/3) = nan(length(satcue{k}{jj}(:,14)),1);             % value of image on B given in previous rating trials (-300 to 300)
        val_A_Rank{k}(:,jj/3) = nan(length(satcue{k}{jj}(:,14)),1);         % rank of image on A, calculated after all the choice trials (0 to 20)
        val_B_Rank{k}(:,jj/3) = nan(length(satcue{k}{jj}(:,14)),1);         % rank of image on B, calculated after all the choice trials (0 to 20)
        
        val_L_Rat{k}(:,jj/3) = nan(length(satcue{k}{jj}(:,14)),1);             % value of image on A given in previous rating trials (-300 to 300)
        val_R_Rat{k}(:,jj/3) = nan(length(satcue{k}{jj}(:,14)),1);             % value of image on B given in previous rating trials (-300 to 300)
        val_L_Rank{k}(:,jj/3) = nan(length(satcue{k}{jj}(:,14)),1);         % rank of image on A, calculated after all the choice trials (0 to 20)
        val_R_Rank{k}(:,jj/3) = nan(length(satcue{k}{jj}(:,14)),1);         % rank of image on B, calculated after all the choice trials (0 to 20)
        
        val_C_Rat{k}(:,jj/3) = nan(length(satcue{k}{jj}(:,14)),1);             % value of image on A given in previous rating trials (-300 to 300)
        val_U_Rat{k}(:,jj/3) = nan(length(satcue{k}{jj}(:,14)),1);             % value of image on B given in previous rating trials (-300 to 300)
        val_C_Rank{k}(:,jj/3) = nan(length(satcue{k}{jj}(:,14)),1);         % rank of image on A, calculated after all the choice trials (0 to 20)
        val_U_Rank{k}(:,jj/3) = nan(length(satcue{k}{jj}(:,14)),1);         % rank of image on B, calculated after all the choice trials (0 to 20)
        
        
        
        % transform the output of the likert and convert it so that
        % the negative is A chosen and positive B chosen.
        
        lr_to_ab3(:,jj/3) = cell2mat(satcue{k}{jj}(:,8));  % conversion factor LR to AB for whole experiment
        
        val_diffAB{k}(:,jj/3)=   val_2AFC{k}(:,jj/3);
        binary= ((val_diffAB{k}(:,jj/3)>0) ==lr_to_ab3(:,jj/3));    % compute when Right chosen (positive answer). Thus, ig cnversion factor 1 and oout pos, then R chosen, and R=B, so B chosen. Whenever binary==1, B chosen.
        val_diffAB{k}(binary,jj/3) = abs(val_diffAB{k}(binary,jj/3));  % convert B is positive.                        if negative (L chosen) then 0 in binary, and factor 0, then L=B, thus B chosen. Thus whenever Conv factor == to find(val_2AFC{k}>0), then B chosen.
        val_diffAB{k}(~binary,jj/3)= -abs(val_diffAB{k}(~binary,jj/3)); %  A in negative.                              The opposite as the previous, if pos out, R chosen, and factor 0, then L=B, A chosen, and negative as 1 ≃0
        % finally if output negative, L chosen, and factor 1, then L==A, so A chosen again, negative here as 1 ≃0.
        
        
        
        
        for i=1:20
            
            indL = find(cell2mat(satcue{k}{jj}(:,2))==i);            % when this stimulus was on left
            indR = find(cell2mat(satcue{k}{jj}(:,3))==i);           % when this stimulus was on right
            lr_to_ab = cell2mat(satcue{k}{jj}(indL,8));              % conversion factor from left_rigth to A_B, this factor saved on the column 8
            lr_to_ab2 = cell2mat(satcue{k}{jj}(indR,8));            % conversion factor from left_rigth to A_B
            
            lr_to_cu = cell2mat(satcue{k}{jj}(:,14))>0;              % conversion factor from left_rigth to A_B, just if positive or negative the output, 1 if right chosen, zero if left chosen
            
            
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
            
    
            
            indLast= find(cell2mat(satcue{k}{jj-1}(:,2))==i);

            flavourA{k}(indA,jj/3) = i<11;
            flavourB{k}(indB,jj/3) = i<11;
            flavourC{k}(indC,jj/3) = i<11;
            flavourU{k}(indU,jj/3) = i<11;
            flavourL{k}(indL,jj/3) = i<11;
            flavourR{k}(indR,jj/3) = i<11;
            
            stimOrderA(indA,jj/3) = i;
            stimOrderB(indB,jj/3) = i;
            stimOrderC(indC,jj/3) = i;
            stimOrderU(indU,jj/3) = i;
            stimOrderL(indL,jj/3) = i;
            stimOrderR(indR,jj/3) = i;
            
            
            % Extract RT{k}
            
            shortRT_Rat{k}(:,jj/3) = RT{jj-1};
            shortRT_AFC{k}(:,jj/3) = RT{jj};
            
       
            % visual features.
            
            
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
            Protein_B{k} (indB, jj/3) = Protein(i,1);
            Kohlenhydrate_B{k} (indB,jj/3) = Kohlenhydrate(i,1);
            Zucker_B{k} (indB,jj/3) = Zucker(i,1);
            Fett_B{k} (indB,jj/3) = Fett(i,1);
            Ballaststoffe_B{k} (indB,jj/3) = Ballaststoffe(i,1);
            Salz_B{k} (indB,jj/3) = Salz(i,1);
            
            
        end
    end
    
    % inverted likert scale output
    for i= 1:20
        for hh=1:20
            if i ~=hh
                try
                    indPair1 = find(cell2mat(satcue{k}{3}(:,2))==i & cell2mat(satcue{k}{3}(:,3))==hh); indPair2 =  find(cell2mat(satcue{k}{6}(:,2))==i & cell2mat(satcue{k}{6}(:,3))==hh);
                    
                end
                try
                    indPair1=find(cell2mat(satcue{k}{3}(:,3))==i & cell2mat(satcue{k}{3}(:,2))==hh);  indPair2 =   find(cell2mat(satcue{k}{6}(:,3))==i & cell2mat(satcue{k}{6}(:,2))==hh);
                end
                
                invval_2AFC{k}(indPair1, 1) = val_2AFC{k}(indPair2,2);
                invval_diffAB{k} (indPair1, 1) = val_diffAB{k} (indPair2,2);
                invshortRT_AFC{k} (indPair1,1) = shortRT_AFC{k} (indPair2,2);
                
                
                invval_2AFC{k} (indPair2, 2) = val_2AFC{k}(indPair1,1);
                invval_diffAB{k} (indPair2, 2) = val_diffAB{k} (indPair1,1);
                invshortRT_AFC{k} (indPair2,1) = shortRT_AFC{k} (indPair1,2);
                
            end
        end
    end
    
    
    
    if autoregresive_corr==1
        
        val_rat{k}=     [ val_rat{k}(1:end-1, :)];
        val_2AFC{k} =   [ val_2AFC{k}(1:end-1, :)];
        val_diffAB{k} = [ val_diffAB{k}(1:end-1, :)];
        
        val_C_Rat{k}=  [ val_C_Rat{k}(1:end-1, :)];
        val_U_Rat{k}=  [ val_U_Rat{k}(1:end-1, :)];
        val_A_Rat{k}=  [ val_A_Rat{k}(1:end-1, :)];
        val_B_Rat{k}=  [ val_B_Rat{k}(1:end-1, :)];
        val_L_Rat{k}=  [ val_L_Rat{k}(1:end-1, :)];
        val_R_Rat{k}=  [ val_R_Rat{k}(1:end-1, :)];
        
        
        flavourC{k}=  [ flavourC{k}(1:end-1, :)];
        flavourU{k}=  [ flavourU{k}(1:end-1, :)];
        flavourA{k}=  [ flavourA{k}(1:end-1, :)];
        flavourB{k}=  [ flavourB{k}(1:end-1, :)];
        flavourL{k}=  [ flavourL{k}(1:end-1, :)];
        flavourR{k}=  [ flavourR{k}(1:end-1, :)];
        flavourRat{k}=[ flavourRat{k}(1:end-1, :)];
        
        
        shortRT_AFC{k}= [ shortRT_AFC{k}(1:end-1, :)];
        shortRT_Rat{k}= [ shortRT_Rat{k}(1:end-1, :)];
        diffRT_AFC{k}= [ diffRT_AFC{k}(1:end-1, :)];
        diffRT_Rat{k}= [ diffRT_Rat{k}(1:end-1, :)];
        
        
        val_C_Rank{k}=  [ val_C_Rank{k}(1:end-1, :)];
        val_U_Rank{k}=  [ val_U_Rank{k}(1:end-1, :)];
        val_A_Rank{k}=  [ val_A_Rank{k}(1:end-1, :)];
        val_B_Rank{k}=  [ val_B_Rank{k}(1:end-1, :)];
        val_L_Rank{k}=  [ val_L_Rank{k}(1:end-1, :)];
        val_R_Rank{k}=  [ val_R_Rank{k}(1:end-1, :)];
        
        
        val_rat_Rat{k}= [ val_rat_Rat{k}(1:end-1, :)];
        val_rat_Rank{k}= [ val_rat_Rank{k}(1:end-1, :)];
        invval_rat_Rat{k}=invval_rat_Rat{k}(1:end-1, :);
        invval_rat_Rank{k}=invval_rat_Rank{k}(1:end-1, :);
        
        
    end
    k=k+1;
end








for kk=1:numel(val_2AFC)
    val_AbsPref{kk} = abs(val_2AFC{kk}) ;
    val_AbsDiff{kk}  = abs(val_A_Rat{kk} -val_B_Rat{kk}) ;
    val_AbsDiffRank{kk}  = abs(val_A_Rank{kk} - val_B_Rank{kk}) ;
    val_DiffABRat{kk}  = val_A_Rat{kk} -val_B_Rat{kk} ;
    val_DiffABRank{kk}  = val_A_Rank{kk} - val_B_Rank{kk} ;
    
    val_DiffCU{kk}  = val_C_Rat{kk} -val_U_Rat{kk} ;
    val_DiffCURank{kk}  = val_C_Rank{kk} - val_U_Rank{kk} ;
    val_DiffLR{kk}  = val_L_Rat{kk} -val_R_Rat{kk} ;
    val_DiffLRRank{kk}  = val_L_Rank{kk} - val_R_Rank{kk} ;
    val_Sum{kk}  = val_A_Rat{kk} +val_B_Rat{kk} ;
    val_SumRank{kk}  = val_A_Rank{kk} + val_B_Rank{kk} ;
    val_ABbin{kk}  = val_diffAB{kk}>0;
    val_LRbin{kk}  = val_2AFC{kk}>0;
    
    sal_A_Rat{kk}  = abs(val_A_Rat{kk});
    sal_B_Rat{kk}  = abs(val_B_Rat{kk});
    sal_rat_Rat{kk}  = abs(val_rat_Rat{kk});
    sal_rat{kk}  = abs(val_rat{kk});
    
end


variablesAFCsess = {'val_A_Rat', 'val_B_Rat', 'sal_A_Rat', 'sal_B_Rat','val_C_Rat', 'val_U_Rat', 'val_A_Rank', 'val_B_Rank', 'val_C_Rank', 'val_U_Rank', 'val_L_Rat', 'val_R_Rat',  'val_L_Rank', 'val_R_Rank',...
    'flavourA', 'flavourB', 'flavourC', 'flavourU', 'flavourL', 'flavourR', 'flavourC', 'flavourU', ...
    'val_AbsDiff', 'val_AbsDiffRank', 'val_AbsPref', 'shortRT_AFC','val_DiffABRat', 'val_DiffABRank','val_DiffLR', 'val_DiffLRRank', 'val_Sum', 'val_SumRank', ...
    'val_2AFC','val_diffAB', 'val_LRbin','val_ABbin', 'Red_A', 'Green_A' , 'Blue_A' , 'Luminance_A',  'Contrast_A',  'Hue_A', 'Saturation_A',   'Red_B', 'Green_B' , 'Blue_B' , 'Luminance_B',  'Contrast_B',  'Hue_B', 'Saturation_B',...
    'Kalorien_A', 'Protein_A',  'Kohlenhydrate_A', 'Zucker_A', 'Fett_A', 'Ballaststoffe_A',  'Salz_A', 'Kalorien_B', 'Protein_B',  'Kohlenhydrate_B', 'Zucker_B', 'Fett_B', 'Ballaststoffe_B',  'Salz_B'};



% nutritional  features
variablesRatsess = { 'val_rat', 'sal_rat', 'val_rat_Rat', 'val_rat_Rank', 'sal_rat_Rat', 'flavourRat', 'shortRT_Rat', 'Red_Rat', 'Green_Rat', 'Blue_Rat', 'Luminance_Rat', 'Contrast_Rat', 'Hue_Rat', 'Saturation_Rat',...
    'Kalorien_Rat', 'Protein_Rat', 'Kohlenhydrate_Rat' , 'Zucker_Rat', 'Fett_Rat', 'Ballaststoffe_Rat', 'Salz_Rat'};


variablesAFCsess= variablesAFCsess(:)';
variablesRatsess= variablesRatsess(:)';

for ll = 1:size(variablesAFCsess,2)
    for mm = 1: size(variablesAFCsess,2)
        name1 = variablesAFCsess{ll};
        name2 = variablesAFCsess{mm};

        for k= 1:size(eval(name1),2)
            Pcoefpre(ll,mm,k) = corr(eval([name1,'{k}(:,1)']), eval([name2,'{k}(:,1)']), 'type', 'P');
            Scoefpre(ll,mm,k) = corr(eval([name1,'{k}(:,1)']), eval([name2,'{k}(:,1)']), 'type', 'S');
            
            Pcoefpost(ll,mm,k) = corr(eval([name1,'{k}(:,2)']), eval([name2,'{k}(:,2)']), 'type', 'P');
            Scoefpost(ll,mm,k) = corr(eval([name1,'{k}(:,2)']), eval([name2,'{k}(:,2)']), 'type', 'S');
        end
        %        end
        
    end
end


close all
figure('color', 'w', 'units', 'normalized', 'position', [0,0,1, 1], 'visible', 'off');
imagesc((mean(Pcoefpre,3) + mean(Pcoefpost,3))/2, [-1,1]);
colormap jet
colorbar
gg=gca;
gg.XAxis.TickValues=[1:size(variablesAFCsess,2)];
gg.XAxis.TickLabels = variablesAFCsess;
gg.XAxis.TickDirection= 'out';
gg.XAxis.TickLabelRotation= 45;

gg.YAxis.TickDirection= 'out';
gg.YAxis.TickValues=[1:size(variablesAFCsess,2)];
gg.YAxis.TickLabels = variablesAFCsess;

print(['./ SimmilarityPearson2AFC'],'-dpng', '-r200');




close all
figure('color', 'w', 'units', 'normalized', 'position', [0,0,1, 1], 'visible', 'off');
imagesc((mean(abs(Pcoefpre),3) + mean(abs(Pcoefpost),3))/2, [0,1]);
colormap jet
colorbar
gg=gca;
gg.XAxis.TickValues=[1:size(variablesAFCsess,2)];
gg.XAxis.TickLabels = variablesAFCsess;
gg.XAxis.TickDirection= 'out';
gg.XAxis.TickLabelRotation= 45;

gg.YAxis.TickDirection= 'out';
gg.YAxis.TickValues=[1:size(variablesAFCsess,2)];
gg.YAxis.TickLabels = variablesAFCsess;

print(['./ SimmilarityABSPearson2AFC'],'-dpng', '-r200');



%% 


close all
figure('color', 'w', 'units', 'normalized', 'position', [0,0,1, 1], 'visible', 'off');
imagesc((mean(Scoefpre,3) + mean(Scoefpost,3))/2, [-1,1]);
colormap jet
colorbar
gg=gca;
gg.XAxis.TickValues=[1:size(variablesAFCsess,2)];
gg.XAxis.TickLabels = variablesAFCsess;
gg.XAxis.TickDirection= 'out';
gg.XAxis.TickLabelRotation= 45;

gg.YAxis.TickDirection= 'out';
gg.YAxis.TickValues=[1:size(variablesAFCsess,2)];
gg.YAxis.TickLabels = variablesAFCsess;

print(['./SimmilaritySpearman2AFC'],'-dpng', '-r200');




close all
figure('color', 'w', 'units', 'normalized', 'position', [0,0,1, 1], 'visible', 'off');
imagesc((mean(abs(Scoefpre),3) + mean(abs(Scoefpost),3))/2, [0,1]);
colormap jet
colorbar
gg=gca;
gg.XAxis.TickValues=[1:size(variablesAFCsess,2)];
gg.XAxis.TickLabels = variablesAFCsess;
gg.XAxis.TickDirection= 'out';
gg.XAxis.TickLabelRotation= 45;

gg.YAxis.TickDirection= 'out';
gg.YAxis.TickValues=[1:size(variablesAFCsess,2)];
gg.YAxis.TickLabels = variablesAFCsess;

print(['./SimmilarityABSSpearman2AFC'],'-dpng', '-r200');









variablesRatsess= variablesRatsess(:)';

for ll = 1:size(variablesRatsess,2)
    for mm = 1: size(variablesRatsess,2)
        name1 = variablesRatsess{ll};
        name2 = variablesRatsess{mm};

        for k= 1:size(eval(name1),2)
            Pcoefpre(ll,mm,k) = corr(eval([name1,'{k}(:,1)']), eval([name2,'{k}(:,1)']), 'type', 'P');
            Scoefpre(ll,mm,k) = corr(eval([name1,'{k}(:,1)']), eval([name2,'{k}(:,1)']), 'type', 'S');
            
            Pcoefpost(ll,mm,k) = corr(eval([name1,'{k}(:,2)']), eval([name2,'{k}(:,2)']), 'type', 'P');
            Scoefpost(ll,mm,k) = corr(eval([name1,'{k}(:,2)']), eval([name2,'{k}(:,2)']), 'type', 'S');
        end
        %        end
        
    end
end


close all
figure('color', 'w', 'units', 'normalized', 'position', [0,0,1, 1], 'visible', 'off');
imagesc((mean(Pcoefpre,3) + mean(Pcoefpost,3))/2, [-1,1]);
colormap jet
colorbar
gg=gca;
gg.XAxis.TickValues=[1:size(variablesRatsess,2)];
gg.XAxis.TickLabels = variablesRatsess;
gg.XAxis.TickDirection= 'out';
gg.XAxis.TickLabelRotation= 45;

gg.YAxis.TickDirection= 'out';
gg.YAxis.TickValues=[1:size(variablesRatsess,2)];
gg.YAxis.TickLabels = variablesRatsess;

print(['./ SimmilarityPearsonVAL'],'-dpng', '-r200');




close all
figure('color', 'w', 'units', 'normalized', 'position', [0,0,1, 1], 'visible', 'off');
imagesc((mean(abs(Pcoefpre),3) + mean(abs(Pcoefpost),3))/2, [0,1]);
colormap jet
colorbar
gg=gca;
gg.XAxis.TickValues=[1:size(variablesAFCsess,2)];
gg.XAxis.TickLabels = variablesAFCsess;
gg.XAxis.TickDirection= 'out';
gg.XAxis.TickLabelRotation= 45;

gg.YAxis.TickDirection= 'out';
gg.YAxis.TickValues=[1:size(variablesAFCsess,2)];
gg.YAxis.TickLabels = variablesAFCsess;

print(['./ SimmilarityABSPearsonVAL'],'-dpng', '-r200');



%% 


close all
figure('color', 'w', 'units', 'normalized', 'position', [0,0,1, 1], 'visible', 'off');
imagesc((mean(Scoefpre,3) + mean(Scoefpost,3))/2, [-1,1]);
colormap jet
colorbar
gg=gca;
gg.XAxis.TickValues=[1:size(variablesAFCsess,2)];
gg.XAxis.TickLabels = variablesAFCsess;
gg.XAxis.TickDirection= 'out';
gg.XAxis.TickLabelRotation= 45;

gg.YAxis.TickDirection= 'out';
gg.YAxis.TickValues=[1:size(variablesAFCsess,2)];
gg.YAxis.TickLabels = variablesAFCsess;

print(['./SimmilaritySpearmanVAL'],'-dpng', '-r200');




close all
figure('color', 'w', 'units', 'normalized', 'position', [0,0,1, 1], 'visible', 'off');
imagesc((mean(abs(Scoefpre),3) + mean(abs(Scoefpost),3))/2, [0,1]);
colormap jet
colorbar
gg=gca;
gg.XAxis.TickValues=[1:size(variablesAFCsess,2)];
gg.XAxis.TickLabels = variablesAFCsess;
gg.XAxis.TickDirection= 'out';
gg.XAxis.TickLabelRotation= 45;

gg.YAxis.TickDirection= 'out';
gg.YAxis.TickValues=[1:size(variablesAFCsess,2)];
gg.YAxis.TickLabels = variablesAFCsess;

print(['./SimmilarityABSSpearmanVAL'],'-dpng', '-r200');


