function generateplot2


%% This script generates the Figures 2 and S1
%
%

% %  Plot all behavioural plots for the paper

folder2save =  '/media/Projects/Alex/Reclustered analysis All/PaperPlots'; mkdir(folder2save)
sizeletra = 30;
folder2go='/media/Projects/Alex/Reclustered analysis All';
set(groot,'defaultAxesLinewidth', 1)

possible_colors = [[0 0.4470 0.7410];...
    [0.9290 0.5940 0.1250];...  % [0.3250 0.8000  0.2980  ];...
    [0.9290 0.3250 0.0980];...
    [0.4940 0.1840 0.5560]];


   
f={};
cd(folder2go);      %/Volumes/MTL/CS4 output')   %orginal form michael '/Volumes/Barnacle/Data/Satiation_data/Satiation_session_data'
DirName= dir;
for a = 1:length(DirName)
    if ~isempty(regexp( DirName(a).name, '^[0-9]{3}$', 'once')) % && (size(DirName(a).name,2) ==3)
        f{end+1,1} = [DirName(a).name];
    end
end
assert(length(f) == 12)



% Consumed products by subjects
Consumedindex= [9,  11, 20,  2, 14, 2, 16,  nan, 16,  4,  7,  17];
Consumedindex2 = {{'001',9}, {'002', 11},  {'003', 20}, {'004', 2}, {'005', 14}, {'006', 2}, {'007', 16}, ...
    {'008', nan}, {'009', 16}, {'010', 4}, {'011', 7}, {'012', 17}};
Consumedindex3 = [nan(1,27), 9,nan(1,2)  11, 20,  2, 14, 2, nan(1,1), 16, nan(1,1), nan, 16,  nan(1,1), 4, nan(1,1), 7,nan(1,1),  17];


%% Extract the behavioural data
k=1;
for h=1:length(f)
   
    
    [satcue, ranking, events] = load_behavioural_data(f{h});
    [RT, endResponse] = extract_reaction_times(events);
    satcuecell{k}=satcue;
    RTcell{k}= RT;
    
    [Red, Green, Blue, Luminance, Contrast, Hue, Saturation, ~ ] = extract_visualFeatures(f{h});
    [Kalorien, Protein, Kohlenhydrate, Zucker, Fett, Ballaststoffe, Salz] = naehwerte(f{h});
    

    %% Copy paste from regression_final
    
    for jj=[3,6]
        % Output values (real answers)
        val_rat{k}(:,jj/3)= cell2mat(satcue{jj-1}(:,14));
        val_2AFC{k}(:,jj/3)= cell2mat(satcue{jj}(:,14));
        
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
                invvalTot = ranking{3}(i,5);                         % ranking (0 to 20)
                
                
                indinv= find(cell2mat(satcue{2}(:,2))==i); assert(length(indinv)==3);
                invRT = RT{2}(indinv,1);
                invVar = cell2mat(satcue{2}(indinv,14));
            end
            
            
            val_rat_Rat{k}(ind3,jj/3)= valRat;
            val_rat_Rank{k}(ind3,jj/3)= valRank;
            val_rat_Tot{k}(ind3,jj/3)= valTot;
            
            invval_rat_Rat{k}(ind3,jj/3)= invvalRat;
            invval_rat_Rank{k}(ind3,jj/3)= invvalRank;
            invval_rat_Tot{k}(ind3,jj/3)= invvalTot;
            
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
            
            
            RTmean_Rat{k}(ind3, jj/3)  =  mean(RT{jj-1}(ind3,1)); % mean RT for satiation effect in first task          
           
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
        
        val_SignedConfAB{k}(:,jj/3)=   val_2AFC{k}(:,jj/3);
        binary= ((val_SignedConfAB{k}(:,jj/3)>0) ==lr_to_ab3(:,jj/3));    % compute when Right chosen (positive answer). Thus, ig cnversion factor 1 and oout pos, then R chosen, and R=B, so B chosen. Whenever binary==1, B chosen.
        val_SignedConfAB{k}(binary,jj/3) = abs(val_SignedConfAB{k}(binary,jj/3));  % convert B is positive.                        if negative (L chosen) then 0 in binary, and factor 0, then L=B, thus B chosen. Thus whenever Conv factor == to find(val_2AFC>0), then B chosen.
        val_SignedConfAB{k}(~binary,jj/3)= -abs(val_SignedConfAB{k}(~binary,jj/3)); %  A in negative.                              The opposite as the previous, if pos out, R chosen, and factor 0, then L=B, A chosen, and negative as 1 ≃0
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
            
            
            val_A_Tot{k}(indA, jj/3) = valTot;
            val_B_Tot{k}(indB, jj/3) = valTot;
            val_L_Tot{k}(indL, jj/3) = valTot;
            val_R_Tot{k}(indR, jj/3) = valTot;
            val_C_Tot{k}(indC, jj/3) = valTot;
            val_U_Tot{k}(indU, jj/3) = valTot;
            
            
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
            
            invval_A_Tot{k}(indA, jj/3) = invvalTot;
            invval_B_Tot{k}(indB, jj/3) = invvalTot;
            invval_L_Tot{k}(indL, jj/3) = invvalTot;
            invval_R_Tot{k}(indR, jj/3) = invvalTot;
            invval_C_Tot{k}(indC, jj/3) = invvalTot;
            invval_U_Tot{k}(indU, jj/3) = invvalTot;
            
            
            
            
            
            indLast= find(cell2mat(satcue{jj-1}(:,2))==i);
            
            val_A_LastRat{k}(indA,jj/3) = satcue{jj-1}{indLast(length(indLast)),14};
            val_B_LastRat{k}(indB,jj/3) = satcue{jj-1}{indLast(length(indLast)),14};
            val_L_LastRat{k}(indL,jj/3) = satcue{jj-1}{indLast(length(indLast)),14};
            val_R_LastRat{k}(indR,jj/3) = satcue{jj-1}{indLast(length(indLast)),14};
            val_C_LastRat{k}(indC, jj/3) = satcue{jj-1}{indLast(length(indLast)),14};
            val_U_LastRat{k}(indU, jj/3) = satcue{jj-1}{indLast(length(indLast)),14};
            
            
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
        
        
        
        % Extract RT
        
%         longRT_Rat(:,jj/3) = endResponse{jj-1};
%         longRT_AFC(:,jj/3) = endResponse{jj};
        
        shortRT_Rat{k}(:,jj/3) = RT{jj-1};
        shortRT_AFC{k}(:,jj/3) = RT{jj};
        
        
        otherRT_Rat(:,jj/3)= cell2mat(satcue{jj-1}(:,13))*1000;
        otherRT_AFC(:,jj/3)= cell2mat(satcue{jj}(:,13))*1000;
        
        diffRT_Rat(:,jj/3)= otherRT_Rat(:,jj/3) - shortRT_Rat{k}(:,jj/3);
        diffRT_AFC(:,jj/3)= otherRT_AFC(:,jj/3) - shortRT_AFC{k}(:,jj/3);
        
        
        
        
        dweeltime{k}(:,jj/3)= cell2mat(satcue{jj}(:,8))+1;
        left_chosen{k}(:,jj/3) = double(val_2AFC{k}(:,jj/3) < 0);
        A_chosen{k}(:,jj/3)= double(val_SignedConfAB{k}(:,jj/3) < 0);
        
        
        
        % test for interaction left right A,B
        FTMat= [sum(A_chosen{k}(:,jj/3) & left_chosen{k}(:,jj/3)), sum(~A_chosen{k}(:,jj/3) & left_chosen{k}(:,jj/3));...
            sum(A_chosen{k}(:,jj/3) & ~left_chosen{k}(:,jj/3)), sum(~A_chosen{k}(:,jj/3) & ~left_chosen{k}(:,jj/3))];
        
        [~, PvalFisher(k, jj/3)] = fishertest(FTMat, 'tail', 'both');
        
        pvalBinoAB(k, jj/3) =  myBinomTest(sum(A_chosen{k}(:,jj/3)), 190, 0.50, 'two');
        pvalBinoLR(k, jj/3) =  myBinomTest(sum(left_chosen{k}(:,jj/3)), 190, 0.50, 'two');
        
    end
    
    
    FTMat= [sum(sum(A_chosen{k} & left_chosen{k})), sum(sum(~A_chosen{k} & left_chosen{k}));...
        sum(sum(A_chosen{k} & ~left_chosen{k})), sum(sum(~A_chosen{k} & ~left_chosen{k}))];
    [~, PvalFisher2(k, 1)] = fishertest(FTMat, 'tail', 'both');
    
    pvalBinoAB2(k, 1) = myBinomTest(sum(sum(~A_chosen{k})), 380, 0.50, 'two');
    pvalBinoLR2(k, 1) = myBinomTest(sum(sum(left_chosen{k})), 380, 0.50, 'two');
    
    
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
                
                
                invval_2AFC{k}(indPair1, 1) = val_2AFC{k}(indPair2,2);
                invval_SignedConfAB{k} (indPair1, 1) = val_SignedConfAB{k} (indPair2,2);
                invshortRT_AFC{k} (indPair1,1) = shortRT_AFC{k} (indPair2,2);
                
                
                invval_2AFC{k} (indPair2, 2) = val_2AFC{k}(indPair1,1);
                invval_SignedConfAB{k} (indPair2, 2) = val_SignedConfAB{k} (indPair1,1);
                invshortRT_AFC{k} (indPair2,2) = shortRT_AFC{k} (indPair1,1);
                
            end
        end
    end
        
    k=k+1;
end





%% new analyses with new variables

% concatenate all subjects in one single vector for further analyses


zvalChosen = [];
zvalUnchosen = [];
zvalChosenConf = [];
zvalUnchosenConf = [];
accuracyABRat = [];

sessnum=[]; PID=[];
L_chosenpar=[];
A_chosenpar=[];



zPrefLR=[]; zPrefAB=[]; zConfidence=[];
PrefLR=[]; PrefAB=[]; Confidence=[];
zdiffABRat=[]; zdiffABRank=[]; 
zdiffLRRat=[]; zdiffLRRank=[]; 


zabsdiffABRat=[]; zabsdiffABRank=[]; zabsdiffABTot=[];
zsumRat=[]; zsumRank=[]; zsumTot=[];
zdweeltime2=[];indeces2avoid=[]; signsum=[];medianConf=[]; highconfidence=[];
signA=[]; signB=[]; avoidIndA=[]; avoidIndB=[];
zA_chosenpar=[];
zRT=[];
RT2=[];
dweeltime2=[];

stepSize = length(val_2AFC{1}(:,1));
diffABRat=[];


zRT_val=[]; zval= []; zmag=[];
zeffect_trial=[]; dummy_trial=[];
effect_trial=[]; PIDTrial=[];
sessTrial=[];

lengthexp = 60;


for i = 1:size(f,1)
    
    for g= 1:2
        
        sessnum(end+1:end+ stepSize,1)= g;
        PID(end+1:end+ stepSize,1)= i;
     
        L_chosenpar(end+1:end+stepSize,1) = left_chosen{i}(:, g);
        A_chosenpar(end+1:end+ stepSize,1)   = A_chosen{i}(:, g);
        zA_chosenpar(end+1:end+ stepSize,1)   = zscore(A_chosen{i}(:, g));
        
        
        PrefLR(end+1:end+ stepSize,1)= (val_2AFC{i}(:, g));
        PrefAB(end+1:end+ stepSize,1)= (val_SignedConfAB{i}(:, g));
        Confidence(end+1:end+ stepSize,1)= (abs(val_SignedConfAB{i}(:, g)));
        
        
        
        zPrefLR(end+1:end+ stepSize,1)= zscore(val_2AFC{i}(:, g));
        zPrefAB(end+1:end+ stepSize,1)= zscore(val_SignedConfAB{i}(:, g));
        zConfidence(end+1:end+ stepSize,1)= zscore(abs(val_SignedConfAB{i}(:, g)));
        
        
        % chosen and unchosen value and beleif confidence (cf. Annika Boldt 1,2, Charles Blundell 3 and Benedetto De Martino, in Bioarxiv, Confidence modulates exploration and exploitation in value-based learning)
        
        zvalChosen(end+1:end+ stepSize,1)  = zscore(val_C_Rat{i}(:,g));
        zvalUnchosen (end+1:end+ stepSize,1) = zscore(val_U_Rat{i}(:,g));
        zvalChosenConf(end+1:end+ stepSize,1)  = zscore(beliefConfidenceChosen{i}(:,g));
        zvalUnchosenConf(end+1:end+ stepSize,1)  = zscore(beliefConfidenceUnchosen{i}(:,g));
        accuracyABRat(end+1:end+ stepSize,1) = (val_A_Rat{i}(:,g) - val_B_Rat{i}(:,g))>0;
        
        
        
        % diff AB
        zvalRat  = zscore(val_A_Rat{i}(:,g) - val_B_Rat{i}(:,g));
        zvalRank = zscore(val_A_Rank{i}(:,g) - val_B_Rank{i}(:,g));
          
        zdiffABRat(end+1:end+ stepSize,1)= zvalRat;
        zdiffABRank(end+1:end+ stepSize,1)= zvalRank;

        
        %%
        valRat12  = (val_A_Rat{i}(:,g) - val_B_Rat{i}(:,g));
        diffABRat(end+1:end+ stepSize,1)= valRat12;
        %%
        
        % diff LR
        
        zvalRat  = zscore(val_L_Rat{i}(:,g) - val_R_Rat{i}(:,g));
        zvalRank = zscore(val_L_Rank{i}(:,g) - val_R_Rank{i}(:,g));
    
        
        zdiffLRRat(end+1:end+ stepSize,1)= zvalRat;
        zdiffLRRank(end+1:end+ stepSize,1)= zvalRank;

  
        zvalRat  = zscore(val_C_Rat{i}(:,g) - val_U_Rat{i}(:,g));
        zvalRank = zscore(val_C_Rank{i}(:,g) - val_U_Rank{i}(:,g));
       
        % abs diff AB
        zvalRat  = zscore(abs(val_A_Rat{i}(:,g) - val_B_Rat{i}(:,g)));
        zvalRank = zscore(abs(val_A_Rank{i}(:,g) - val_B_Rank{i}(:,g)));
     
        zabsdiffABRat(end+1:end+ stepSize,1)= zvalRat;
        zabsdiffABRank(end+1:end+ stepSize,1)= zvalRank;
 
        % sum
        
        zvalRat  = zscore(val_A_Rat{i}(:,g) + val_B_Rat{i}(:,g));
        zvalRank = zscore(val_A_Rank{i}(:,g) + val_B_Rank{i}(:,g));
   
        zsumRat(end+1:end+ stepSize,1)= zvalRat;
        zsumRank(end+1:end+ stepSize,1)= zvalRank;
    
        % other
        RT2(end+1:end+ stepSize,1)=(shortRT_AFC{i}(:, g));
        zRT(end+1:end+ stepSize,1)=zscore(shortRT_AFC{i}(:, g));
        dweeltime2(end+1:end+ stepSize,1)= (dweeltime{i}(:, g));
        zdweeltime2(end+1:end+ stepSize,1)= zscore(dweeltime{i}(:, g));
        indeces2avoid(end+1:end+ stepSize,1) = [val_L_Rat{i}(:,g) + val_R_Rat{i}(:,g)] == 0;
        signsum      (end+1:end+ stepSize,1) = [val_L_Rat{i}(:,g) + val_R_Rat{i}(:,g)]  < 0;
        
        medianConf(end+1:end+ stepSize,1)=median(abs(val_2AFC{i}(:, g)));
        highconfidence(end+1:end+ stepSize,1) = abs(val_2AFC{i}(:, g))>50;
        
        signA(end+1:end+ stepSize,1)= val_A_Rat{i}(:,g)>0;
        avoidIndA(end+1:end+ stepSize,1)= (val_A_Rat{i}(:,g))==0;
        signB(end+1:end+ stepSize,1)= val_B_Rat{i}(:,g)>0;
        avoidIndB(end+1:end+ stepSize,1)= (val_B_Rat{i}(:, g))==0;
        
        
        %      CorrDiff(i,g)=corr((val_2AFC{i}(:,g)),  (val_diff{i}(:,g)), 'type', 'Spearman');
        
        
        % Valuation task
        
            
        zRT_val(end+1:end+lengthexp,1)=zscore(shortRT_Rat{i}(:,g));
        zmag(end+1:end+lengthexp,1)=zscore(abs(val_rat {i}(:,g)));
        zval(end+1:end+lengthexp,1)=zscore((val_rat {i}(:,g)));
        
        
        
        
        dummy_trial(end+1:end+(lengthexp),1)=   val_rat_Rat {i}(:,g) >0;
        effect_trial(end+1:end+(lengthexp),1)= (val_rat_Rat {i}(:,g) >0) -0.5;
        zeffect_trial(end+1:end+(lengthexp),1)= zscore(val_rat_Rat {i}(:,g) >0);
        
        
        PIDTrial(end+1:end+(lengthexp),1)= ones(lengthexp,1) * i;
        sessTrial(end+1:end+(lengthexp),1)= ones(lengthexp,1) * g;
        
        
    end
end



%% Reaction times Models in Valuation task
% since patietnes are forced to answer diffreten than zero, no problem in
% valence analyses


tbl_lme2= table(zRT_val , zmag, zval,   dummy_trial, effect_trial, zeffect_trial , PIDTrial, sessTrial, 'VariableNames', {'RT', 'Magnitude', 'Value', 'Valence_dummy', 'Valence_effect', 'zValence_effect' , 'PID', 'sess'} );
tbl_lme2.PID=categorical(tbl_lme2.PID);
tbl_lme2.sess=categorical(tbl_lme2.sess);

lme_RT_valuation1 = fitlme(tbl_lme2, 'RT~ 1+Magnitude + (Magnitude|PID) + (Magnitude|PID:sess)');
lme_RT_valuation2 = fitlme(tbl_lme2, 'RT~ 1+Magnitude*Valence_effect  + (Magnitude|PID) + (Magnitude|PID:sess)');
lme_RT_valuation3 = fitlme(tbl_lme2, 'RT~ 1+Magnitude*zValence_effect + (Magnitude|PID) + (Magnitude|PID:sess)');








%% Reaction times Models in 2AFC task
tabl = table(zRT, zabsdiffABRat,  zConfidence, zsumRat, PID, sessnum,   A_chosenpar-0.5, zA_chosenpar, 'VariableNames', {'RT', 'absDV',  'Confidence' ,  'Sum',  'PID', 'sess',  'Achosen', 'AchosenZ'});

RT_APmdl1 =  fitlme(tabl, 'RT  ~  Confidence + (Confidence  |PID) + (Confidence|PID:sess)');
RT_APmdl2 =  fitlme(tabl, 'RT  ~  Confidence*Achosen + (Confidence*Achosen |PID)  + (Confidence*Achosen |PID:sess)');
RT_APmdl3 =  fitlme(tabl, 'RT  ~  Confidence*AchosenZ + (Confidence*AchosenZ |PID) + (Confidence*AchosenZ |PID:sess)');


%% confidence analysis

mdl1 = fitlme(tabl, 'Confidence ~ RT  + (RT |PID)  + (RT |PID:sess)');
mdl2= fitlme(tabl, 'Confidence ~ Sum  + (Sum |PID) + (Sum |PID:sess)');
mdl3= fitlme(tabl, 'Confidence ~ absDV  + (absDV |PID) + (absDV|PID:sess)');
mdl4= fitlme(tabl, 'Confidence ~ RT   + Sum + (RT   + Sum|PID)+ (RT + Sum |PID:sess)');
mdl5= fitlme(tabl, 'Confidence ~ absDV + RT  + (absDV + RT |PID)+ (absDV + RT |PID:sess)');
mdl6= fitlme(tabl, 'Confidence ~ absDV +  Sum + (absDV + Sum|PID)+ (absDV + Sum |PID:sess)');
mdl7= fitlme(tabl, 'Confidence ~ absDV + RT   + Sum + (absDV + RT + Sum|PID)+ (absDV + RT + Sum |PID:sess)');




%% Logistic regression analyses of binary choice

tabl = table(zRT, zdiffABRat, zConfidence, zsumRat , PID, sessnum, A_chosenpar, zdweeltime2-1,  zPrefAB,  'VariableNames', {'RT', 'DV', 'Confidence' ,  'Sum',  'PID', 'sess', 'Achosen', ...
    'dwelltime',  'PrefAB'} );

mdl10 = fitglme(tabl, 'Achosen ~ DV+ (DV|PID)+ (DV|PID:sess)', 'Distribution','binomial');
mdl13 = fitglme(tabl, 'Achosen ~ DV + Confidence+ (DV + Confidence|PID)+ (DV + Confidence|PID:sess)', 'Distribution','binomial');
mdl21 = fitglme(tabl, 'Achosen ~ DV+Confidence+Sum+ (DV+Confidence+Sum|PID) + (DV+Confidence+Sum|PID:sess)', 'Distribution','binomial');
mdl27 = fitglme(tabl, 'Achosen ~ DV+Confidence+ Sum+ DV:Confidence  +(DV+Confidence+ Sum+ DV:Confidence|PID)  +(DV+Confidence+ Sum+ DV:Confidence|PID:sess)', 'Distribution','binomial');
mdl28 = fitglme(tabl, 'Achosen ~ DV+Confidence+ Sum + DV:Confidence+ DV:Sum +(DV+Confidence+ Sum+ DV:Confidence+ DV:Sum |PID)  + (DV+Confidence+ Sum+ DV:Confidence+ DV:Sum |PID:sess)', 'Distribution','binomial');




%% logistic fits figure (Figure )

for lala = 25:25:100
    ind = find(Confidence==lala);
    tabl = table(diffABRat(ind)/3, PID(ind), sessnum(ind), A_chosenpar(ind), 'VariableNames', { 'DV',  'PID',  'sess', 'Achosen'} );
    mdl30{lala/25}= fitglme(tabl, 'Achosen ~ DV+ (DV|PID)+ (DV|PID:sess)', 'Distribution','binomial');
end



%% Figure 2



%% version 1
close all
yshifletra = +0.2 +1.1;
xshifletra = -2.25;
fontSizepl = 14;
fontylabel = 14;


stardist=0.05; marksize=9;
figBehav1 = figure('units', 'centimeters', 'position', [0,0, 30.6800   30], 'color', 'w', 'visible', 'off');

sub1 = subplot('position', [0.1, 0.76, 0.87, 0.21]);
set(sub1, 'units', 'centimeters', 'position', [2.5400    22.5654   27.5980    5.4062])  ;
mdl2use = mdl7;

Confidenceregplot= plot(mdl2use.Coefficients.Estimate([3,2,4]), 'ok',    'MarkerEdgeColor','k',...
    'MarkerFaceColor', 'k', 'MarkerSize', marksize );
hold on;
lengthplot = length(mdl2use.Coefficients.Estimate(2:end));
set(gca, 'Xlim', [0.5,lengthplot+0.5], 'XTick', [1:lengthplot], 'XTicklabel',{ '|DV|', mdl2use.Coefficients.Name{2}, mdl2use.Coefficients.Name{4}})
plot([0.5,lengthplot+0.5], [0,0], '--k')
errorbar([1:lengthplot], mdl2use.Coefficients.Estimate([3,2,4])', 1.96 *mdl2use.Coefficients.SE([3,2,4])', 1.96 *mdl2use.Coefficients.SE([3,2,4])'  ,'.k', 'linewidth', 1.5 );




plot([0.94, 1, 1.06], [mdl2use.Coefficients.Upper(3) + stardist, mdl2use.Coefficients.Upper(3) + stardist,  mdl2use.Coefficients.Upper(3) + stardist], '*k', 'MarkerSize',marksize, 'Marker', 'pentagram', 'MarkerFaceColor', 'k' )
plot([1.94, 2, 2.06], [mdl2use.Coefficients.Lower(2) - stardist, mdl2use.Coefficients.Lower(2) - stardist,  mdl2use.Coefficients.Lower(2) - stardist] , '*k', 'MarkerSize',marksize, 'Marker', 'pentagram', 'MarkerFaceColor', 'k'  )
plot([2.94, 3, 3.06], [mdl2use.Coefficients.Upper(4) + stardist, mdl2use.Coefficients.Upper(4) + stardist,  mdl2use.Coefficients.Upper(4) + stardist], '*k' , 'MarkerSize',marksize, 'Marker', 'pentagram', 'MarkerFaceColor', 'k' )

box off
set(gca,'FontSize', fontSizepl, 'FontWeight', 'bold')
yy1=ylabel ('Fixed-Effect Coefficients'); yy1.FontSize=fontylabel;


yy1.Units = 'centimeters';
yy1.Position = [-1.4313     2.85         0];

ylim([-0.3, 0.7])
title ('Confidence', 'fontsize', 20)

anot1 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',['A'],'Units','normalized','FontName','Helvetica','FontSize', sizeletra,'EdgeColor','none', 'fontweight', 'bold');
set(anot1, 'units', 'centimeters', 'position', [2.5400+xshifletra    21.5654+5.4062+yshifletra   2 2] )
        






stardist=0.5;

sub2 = subplot('position', [0.1, 0.45, 0.87, 0.21]);
set(sub2, 'units', 'centimeters', 'position', [2.5400    13.5848   27.5980    5.4062])   

model2use= mdl28;
Confidenceregplot= plot(model2use.Coefficients.Estimate(2:end), 'ok',    'MarkerEdgeColor','k',...
    'MarkerFaceColor', 'k', 'MarkerSize',marksize);

hold on;
lengthplot = length(model2use.Coefficients.Estimate(2:end));
set(gca, 'Xlim', [0.5,lengthplot+0.5], 'XTick', [1:lengthplot], 'XTicklabel',{model2use.Coefficients.Name{2} model2use.Coefficients.Name{3}, ...
    model2use.Coefficients.Name{4}, [model2use.Coefficients.Name{2}, ' x ', model2use.Coefficients.Name{3}], [model2use.Coefficients.Name{2}, ' x ', model2use.Coefficients.Name{4}]})
plot([0.5,lengthplot+0.5], [0,0], '--k')
errorbar([1:lengthplot], model2use.Coefficients.Estimate(2:end)', 1.96 *model2use.Coefficients.SE(2:end)', 1.96 *model2use.Coefficients.SE(2:end)'  ,'.k' , 'linewidth', 1.5)


plot([0.9, 1, 1.1], [model2use.Coefficients.Upper(2) + stardist, model2use.Coefficients.Upper(2) + stardist,  model2use.Coefficients.Upper(2) + stardist], '*k', 'MarkerSize',marksize, 'Marker', 'pentagram', 'MarkerFaceColor', 'k'  )
plot([3.9, 4, 4.1], [model2use.Coefficients.Upper(5) + stardist, model2use.Coefficients.Upper(5) + stardist,  model2use.Coefficients.Upper(5) + stardist], '*k' , 'MarkerSize',marksize, 'Marker', 'pentagram', 'MarkerFaceColor', 'k' )
plot([5.0],  [model2use.Coefficients.Lower(6) - stardist], '*k' , 'MarkerSize',marksize, 'Marker', 'pentagram', 'MarkerFaceColor', 'k' )
box off

set(gca,'FontSize', fontSizepl, 'FontWeight', 'bold')
yy2=ylabel ('Fixed-Effect Coefficients'); yy2.FontSize=fontylabel;


yy2.Units = 'centimeters';
yy2.Position = [-1.4313     2.85         0];


ylim([-3, 8])
title ('Choice', 'fontsize', 20)





anot2 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',['B'],'Units','normalized','FontName','Helvetica','FontSize', sizeletra,'EdgeColor','none', 'fontweight', 'bold');
set(anot2, 'units', 'centimeters', 'position', [2.5400+xshifletra    12.5848+5.4062+yshifletra   2 2] )
        





sub3=subplot('position', [0.1, 0.0605, 0.6, 0.32]);
set(sub3, 'units', 'centimeters', 'position', [2.5400    1.5575   19.7400    8.2381])


hold on


fplot(@(x) 1 / (1 + exp(-(x* mdl30{4}.Coefficients.Estimate(2)) )), 'Linewidth',3.2, 'color', possible_colors(1,:));
fplot(@(x) 1 / (1 + exp(-(x* mdl30{3}.Coefficients.Estimate(2)) )), 'Linewidth',3.2, 'color', possible_colors(2,:));

fplot(@(x) 1 / (1 + exp(-(x* mdl30{2}.Coefficients.Estimate(2)) )), 'Linewidth',3.2, 'color', possible_colors(3,:));

fplot(@(x) 1 / (1 + exp(-(x* mdl30{1}.Coefficients.Estimate(2)) )), 'Linewidth',3.2, 'color', possible_colors(4,:));



leg=legend('100','75','50','25');
% title('Absolute Preference')
xlabel('A Mean Rating - B Mean Rating')
set(gca, 'Fontsize', fontSizepl, 'FontWeight', 'bold')


yy3=ylabel('P(A chosen)'); yy3.FontSize=fontylabel;

yy3.Units = 'centimeters';
yy3.Position = [-1.4313    4.1143         0];

box off
% leg.Position=[0.75,0.265,0.152, 0.122];
leg.Units = 'centimeters';
leg.Position = [24.5500    6.6495    3.8608    3.4862];
leg.FontSize=fontylabel;


set(gca, 'xlim', [-200,200],'ylim', [0,1])



anot3 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',['C'],'Units','normalized','FontName','Helvetica','FontSize', sizeletra,'EdgeColor','none', 'fontweight', 'bold');
set(anot3, 'units', 'centimeters', 'position', [2.5400+xshifletra    0.5575+8.2381+yshifletra   2 2] )





print([folder2save, '/Confidence-choice plot 300'], '-dpng', '-r300')
print([folder2save, '/Confidence-choice plot'], '-dpng', '-r500')

print([folder2save, '/Confidence-choice plot'], '-depsc', '-painters')
set(figBehav1,'Units','centimeters');
pos = get(figBehav1,'Position');
set(figBehav1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print([folder2save, '/Confidence-choice plot'], '-dpdf', '-r500', '-opengl')
print([folder2save, '/Confidence-choice plot Vec'], '-dpdf', '-painters')
print([folder2save, '/Confidence-choice plot 300'], '-dpdf', '-r300', '-opengl');







%% version 2
close all
yshifletra = +0.2 +1.1;
xshifletra = -2.25;
fontSizepl = 14;
fontylabel = 14;


stardist=0.05; marksize=9;
figBehav1 = figure('units', 'centimeters', 'position', [0,0, 30.6800   30], 'color', 'w', 'visible', 'off');

sub1 = subplot('position', [0.1, 0.76, 0.87, 0.21]);
set(sub1, 'units', 'centimeters', 'position', [2.5400    22.5654   27.5980    5.4062])  ;
mdl2use = mdl7;

Confidenceregplot= plot(mdl2use.Coefficients.Estimate([3,2,4]), 'ok',    'MarkerEdgeColor','k',...
    'MarkerFaceColor', 'k', 'MarkerSize', marksize );
hold on;
lengthplot = length(mdl2use.Coefficients.Estimate(2:end));
set(gca, 'Xlim', [0.5,lengthplot+0.5], 'XTick', [1:lengthplot], 'XTicklabel',{ '|DV|', mdl2use.Coefficients.Name{2}, mdl2use.Coefficients.Name{4}})
plot([0.5,lengthplot+0.5], [0,0], '--k')
errorbar([1:lengthplot], mdl2use.Coefficients.Estimate([3,2,4])', 1.96 *mdl2use.Coefficients.SE([3,2,4])', 1.96 *mdl2use.Coefficients.SE([3,2,4])'  ,'.k', 'linewidth', 1.5 );




plot([0.94, 1, 1.06], [mdl2use.Coefficients.Upper(3) + stardist, mdl2use.Coefficients.Upper(3) + stardist,  mdl2use.Coefficients.Upper(3) + stardist], '*k', 'MarkerSize',marksize, 'Marker', 'pentagram', 'MarkerFaceColor', 'k' )
plot([1.94, 2, 2.06], [mdl2use.Coefficients.Lower(2) - stardist, mdl2use.Coefficients.Lower(2) - stardist,  mdl2use.Coefficients.Lower(2) - stardist] , '*k', 'MarkerSize',marksize, 'Marker', 'pentagram', 'MarkerFaceColor', 'k'  )
plot([2.94, 3, 3.06], [mdl2use.Coefficients.Upper(4) + stardist, mdl2use.Coefficients.Upper(4) + stardist,  mdl2use.Coefficients.Upper(4) + stardist], '*k' , 'MarkerSize',marksize, 'Marker', 'pentagram', 'MarkerFaceColor', 'k' )

box off
set(gca,'FontSize', fontSizepl, 'FontWeight', 'bold')
yy1=ylabel ('Fixed-Effect Coefficients'); yy1.FontSize=fontylabel;


yy1.Units = 'centimeters';
yy1.Position = [-1.4313     2.85         0];

ylim([-0.3, 0.7])
title ('Confidence', 'fontsize', 20)

anot1 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',['A'],'Units','normalized','FontName','Helvetica','FontSize', sizeletra,'EdgeColor','none', 'fontweight', 'bold');
set(anot1, 'units', 'centimeters', 'position', [2.5400+xshifletra    21.5654+5.4062+yshifletra   2 2] )
        






stardist=0.5;

sub2 = subplot('position', [0.1, 0.45, 0.87, 0.21]);
set(sub2, 'units', 'centimeters', 'position', [2.5400    13.5848   27.5980    5.4062])   

model2use= mdl28;
Confidenceregplot= plot(model2use.Coefficients.Estimate(2:end), 'ok',    'MarkerEdgeColor','k',...
    'MarkerFaceColor', 'k', 'MarkerSize',marksize);

hold on;
lengthplot = length(model2use.Coefficients.Estimate(2:end));
set(gca, 'Xlim', [0.5,lengthplot+0.5], 'XTick', [1:lengthplot], 'XTicklabel',{model2use.Coefficients.Name{2} model2use.Coefficients.Name{3}, ...
    model2use.Coefficients.Name{4}, [model2use.Coefficients.Name{2}, ' x ', model2use.Coefficients.Name{3}], [model2use.Coefficients.Name{2}, ' x ', model2use.Coefficients.Name{4}]})
plot([0.5,lengthplot+0.5], [0,0], '--k')
errorbar([1:lengthplot], model2use.Coefficients.Estimate(2:end)', 1.96 *model2use.Coefficients.SE(2:end)', 1.96 *model2use.Coefficients.SE(2:end)'  ,'.k' , 'linewidth', 1.5)


plot([0.9, 1, 1.1], [model2use.Coefficients.Upper(2) + stardist, model2use.Coefficients.Upper(2) + stardist,  model2use.Coefficients.Upper(2) + stardist], '*k', 'MarkerSize',marksize, 'Marker', 'pentagram', 'MarkerFaceColor', 'k'  )
plot([3.9, 4, 4.1], [model2use.Coefficients.Upper(5) + stardist, model2use.Coefficients.Upper(5) + stardist,  model2use.Coefficients.Upper(5) + stardist], '*k' , 'MarkerSize',marksize, 'Marker', 'pentagram', 'MarkerFaceColor', 'k' )
plot([5.0],  [model2use.Coefficients.Lower(6) - stardist], '*k' , 'MarkerSize',marksize, 'Marker', 'pentagram', 'MarkerFaceColor', 'k' )
box off

set(gca,'FontSize', fontSizepl, 'FontWeight', 'bold')
yy2=ylabel ('Fixed-Effect Coefficients'); yy2.FontSize=fontylabel;


yy2.Units = 'centimeters';
yy2.Position = [-1.4313     2.85         0];



ylim([-3, 8])
title ('Choice', 'fontsize', 20)



anot2 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',['B'],'Units','normalized','FontName','Helvetica','FontSize', sizeletra,'EdgeColor','none', 'fontweight', 'bold');
set(anot2, 'units', 'centimeters', 'position', [2.5400+xshifletra    12.5848+5.4062+yshifletra   2 2] )
        


sub3=subplot('position', [0.1, 0.0605, 0.6, 0.32]);
set(sub3, 'units', 'centimeters', 'position', [2.5400    1.5575   19.7400    8.2381])


hold on


fplot(@(x) 1 / (1 + exp(-((x* mdl30{4}.Coefficients.Estimate(2)) + (mdl30{4}.Coefficients.Estimate(1))))), 'Linewidth',3.2, 'color', possible_colors(1,:));
fplot(@(x) 1 / (1 + exp(-((x* mdl30{3}.Coefficients.Estimate(2)) + (mdl30{3}.Coefficients.Estimate(1))))), 'Linewidth',3.2, 'color', possible_colors(2,:));
fplot(@(x) 1 / (1 + exp(-((x* mdl30{2}.Coefficients.Estimate(2)) + (mdl30{2}.Coefficients.Estimate(1))))), 'Linewidth',3.2, 'color', possible_colors(3,:));
fplot(@(x) 1 / (1 + exp(-((x* mdl30{1}.Coefficients.Estimate(2)) + (mdl30{1}.Coefficients.Estimate(1))))), 'Linewidth',3.2, 'color', possible_colors(4,:));



leg=legend('100','75','50','25');
% title('Absolute Preference')
xlabel('A Mean Rating - B Mean Rating')
set(gca, 'Fontsize', fontSizepl, 'FontWeight', 'bold')


yy3=ylabel('P(A chosen)'); yy3.FontSize=fontylabel;

yy3.Units = 'centimeters';
yy3.Position = [-1.4313    4.1143         0];

box off
% leg.Position=[0.75,0.265,0.152, 0.122];
leg.Units = 'centimeters';
leg.Position = [24.5500    6.6495    3.8608    3.4862];
leg.FontSize=fontylabel;


set(gca, 'xlim', [-200,200],'ylim', [0,1])



anot3 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',['C'],'Units','normalized','FontName','Helvetica','FontSize', sizeletra,'EdgeColor','none', 'fontweight', 'bold');
set(anot3, 'units', 'centimeters', 'position', [2.5400+xshifletra    0.5575+8.2381+yshifletra   2 2] )

print([folder2save, '/Confidence-choice plot 2 300'], '-dpng', '-r300')
print([folder2save, '/Confidence-choice plot 2'], '-dpng', '-r500')

print([folder2save, '/Confidence-choice plot 2'], '-depsc', '-painters')

set(figBehav1,'Units','centimeters');
pos = get(figBehav1,'Position');
set(figBehav1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print([folder2save, '/Confidence-choice plot 2'], '-dpdf', '-r500', '-opengl')
print([folder2save, '/Confidence-choice plot 2 300'], '-dpdf', '-r300', '-opengl');

print([folder2save, '/Confidence-choice plot 2 Vec'], '-dpdf', '-painters')




















%% BIC for supplementary Figure 1 

toplotBIC1= [mdl3.ModelCriterion.BIC, mdl5.ModelCriterion.BIC, mdl6.ModelCriterion.BIC, mdl7.ModelCriterion.BIC ];
toplotBIC2= [mdl10.ModelCriterion.BIC, mdl13.ModelCriterion.BIC, mdl21.ModelCriterion.BIC, mdl27.ModelCriterion.BIC,  mdl28.ModelCriterion.BIC];
toplotBIC3= [lme_RT_valuation1.ModelCriterion.BIC, lme_RT_valuation3.ModelCriterion.BIC];

clearvars -except  toplotBIC1 toplotBIC2 toplotBIC3  folder2save sizeletra sizetrialtype possible_colors




%% Now analysis of the BIC and the satiation


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
        % Output values (real answers)
        
        val_rat{k}(:,jj/3)= cell2mat(satcue{k}{jj-1}(:,14));
        val_2AFC{k}(:,jj/3)= cell2mat(satcue{k}{jj}(:,14));
        
        
        
        % Input values
        
        % Rating paradigm
        
        val_rat_Rat{k}(:,jj/3)=nan(length(satcue{k}{jj-1}(:,14)),1);
        
        val_rat_Rank{k}(:,jj/3)=nan(length(satcue{k}{jj-1}(:,14)),1);
        
        
        % Extract RT{k}
       
        shortRT_Rat{k}(:,jj/3) = RT{jj-1};
        shortRT_AFC{k}(:,jj/3) = RT{jj};
        
    
        
        
        for i=1:20
            
            ind3=find(cell2mat(satcue{k}{jj-1}(:,2))==i); % extract postion from stimulus order
            valRat=  ranking{jj-1}(i,4); % ranking
            
            valRank= ranking{jj}(i,3);
      
            if jj==3
                invvalRat = ranking{5}(i,4);                        % rating value (-300 to 300)
                invvalRank = ranking{6}(i,3);                         % ranking (0 to 20)          
                indinv= find(cell2mat(satcue{k}{5}(:,2))==i); assert(length(indinv)==3); % find the inverse presentations of the stimulus
                invRT = RT{5}(indinv,1); % RT in the other session
                
            elseif jj==6
                
                invvalRat = ranking{2}(i,4);                        % rating value (-300 to 300)
                invvalRank = ranking{3}(i,3);                         % ranking (0 to 20)

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
            
            
                       
            %% fichaxes recentes para a plantilla
            
            RTmean_Rat{k}(ind3, jj/3)  =  mean(RT{jj-1}(ind3,1)); % mean RT for satiation effect in first task
            invRTmean_Rat{k}(ind3, jj/3)  =  mean(invRT); % mean RT for for satiation effect in first task
      

            
            %% productwise for satiation effect nehavioural analsis
            

            
            % visual features
            Red_Rat_Prod{k}(i,jj/3) = Red(i,1);
            Green_Rat_Prod{k}(i,jj/3) = Green(i,1);
            Blue_Rat_Prod{k}(i,jj/3) = Blue(i,1);
            Luminance_Rat_Prod{k}(i,jj/3) = Luminance(i,1);
            Contrast_Rat_Prod{k}(i,jj/3) = Contrast(i,1);
            Hue_Rat_Prod{k}(i,jj/3) = Hue(i,1);
            Saturation_Rat_Prod{k}(i,jj/3) = Saturation(i,1);
            
            % nutritional  features
            
            Kalorien_Rat_Prod{k}(i,jj/3) = Kalorien(i,1);
            Protein_Rat_Prod{k}(i,jj/3) = Protein(i,1);
            Kohlenhydrate_Rat_Prod{k}(i,jj/3) = Kohlenhydrate(i,1);
            Zucker_Rat_Prod{k}(i,jj/3) = Zucker(i,1);
            Fett_Rat_Prod{k}(i,jj/3) = Fett(i,1);
            Ballaststoffe_Rat_Prod{k}(i,jj/3) = Ballaststoffe(i,1);
            Salz_Rat_Prod{k}(i,jj/3) = Salz(i,1);
            
            % RT 
             
            RTmean_Rat_Prod{k}(i, jj/3)  =  mean(RT{jj-1}(ind3,1)); % mean RT for satiation effect in first task

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
    
        val_SignedConfAB{k}(:,jj/3)=   val_2AFC{k}(:,jj/3);
        binary= ((val_SignedConfAB{k}(:,jj/3)>0) ==lr_to_ab3(:,jj/3));    % compute when Right chosen (positive answer). Thus, ig cnversion factor 1 and oout pos, then R chosen, and R=B, so B chosen. Whenever binary==1, B chosen.
        val_SignedConfAB{k}(binary,jj/3) = abs(val_SignedConfAB{k}(binary,jj/3));  % convert B is positive.                        if negative (L chosen) then 0 in binary, and factor 0, then L=B, thus B chosen. Thus whenever Conv factor == to find(val_2AFC{k}>0), then B chosen.
        val_SignedConfAB{k}(~binary,jj/3)= -abs(val_SignedConfAB{k}(~binary,jj/3)); %  A in negative.                              The opposite as the previous, if pos out, R chosen, and factor 0, then L=B, A chosen, and negative as 1 ≃0
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
            elseif jj==6
                invvalRat = ranking{2}(i,4);                        % rating value (-300 to 300)
                invvalRank = ranking{3}(i,3);                         % ranking (0 to 20)
            end
            
            val_A_Rat{k}(indA, jj/3) = valRat;
            val_B_Rat{k}(indB, jj/3) = valRat;
            val_L_Rat{k}(indL, jj/3) = valRat;
            val_R_Rat{k}(indR, jj/3) = valRat;
            val_C_Rat{k}(indC, jj/3) = valRat;
            val_U_Rat{k}(indU, jj/3) = valRat;
            
            val_A_Rank{k}(indA, jj/3) = valRank;
            val_B_Rank{k}(indB, jj/3) = valRank;
            val_L_Rank{k}(indL, jj/3) = valRank;
            val_R_Rank{k}(indR, jj/3) = valRank;
            val_C_Rank{k}(indC, jj/3) = valRank;
            val_U_Rank{k}(indU, jj/3) = valRank;
            
            % inverted order
            
            invval_A_Rat{k}(indA, jj/3) = invvalRat;
            invval_B_Rat{k}(indB, jj/3) = invvalRat;
            invval_L_Rat{k}(indL, jj/3)=  invvalRat;
            invval_R_Rat{k}(indR, jj/3)=  invvalRat;
            invval_C_Rat{k}(indC, jj/3) = invvalRat;
            invval_U_Rat{k}(indU, jj/3) = invvalRat;
            
            invval_A_Rank{k}(indA, jj/3) = invvalRank;
            invval_B_Rank{k}(indB, jj/3) = invvalRank;
            invval_L_Rank{k}(indL, jj/3) = invvalRank;
            invval_R_Rank{k}(indR, jj/3) = invvalRank;
            invval_C_Rank{k}(indC, jj/3) = invvalRank;
            invval_U_Rank{k}(indU, jj/3) = invvalRank;
            

            
            
            
            
            indLast= find(cell2mat(satcue{k}{jj-1}(:,2))==i);
            
            val_A_LastRat{k}(indA,jj/3) = satcue{k}{jj-1}{indLast(length(indLast)),14};
            val_B_LastRat{k}(indB,jj/3) = satcue{k}{jj-1}{indLast(length(indLast)),14};
            val_L_LastRat{k}(indL,jj/3) = satcue{k}{jj-1}{indLast(length(indLast)),14};
            val_R_LastRat{k}(indR,jj/3) = satcue{k}{jj-1}{indLast(length(indLast)),14};
            val_C_LastRat{k}(indC, jj/3) = satcue{k}{jj-1}{indLast(length(indLast)),14};
            val_U_LastRat{k}(indU, jj/3) = satcue{k}{jj-1}{indLast(length(indLast)),14};
            
            
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
                invval_SignedConfAB{k} (indPair1, 1) = val_SignedConfAB{k} (indPair2,2);
                invshortRT_AFC{k} (indPair1,1) = shortRT_AFC{k} (indPair2,2);
                
                
                invval_2AFC{k} (indPair2, 2) = val_2AFC{k}(indPair1,1);
                invval_SignedConfAB{k} (indPair2, 2) = val_SignedConfAB{k} (indPair1,1);
                invshortRT_AFC{k} (indPair2,2) = shortRT_AFC{k} (indPair1,1);
                
                
                
                invval2_C_Rat{k} (indPair1,1) =  val_C_Rat{k}(indPair2,2);  % the other inverse chosen and unchosen are based in the stimulus vbalue before and after,now it is based
                invval2_C_Rat{k} (indPair2,2) =  val_C_Rat{k}(indPair1,1);  %inthe value of the chosen value in the trial type, so if the stimulus chosen changes, it is the value of 									          %the other stimulus when chosen.
                invval2_C_Rank{k} (indPair1,1) =  val_C_Rank{k}(indPair2,2);
                invval2_C_Rank{k} (indPair2,2) =  val_C_Rank{k}(indPair1,1);
                invval2_C_Tot{k} (indPair1,1) =  val_C_Tot{k}(indPair2,2);
                invval2_C_Tot{k} (indPair2,2) =  val_C_Tot{k}(indPair1,1);
                
                invval2_U_Rat{k} (indPair1,1) =  val_U_Rat{k}(indPair2,2);
                invval2_U_Rat{k} (indPair2,2) =  val_U_Rat{k}(indPair1,1);
                invval2_U_Rank{k} (indPair1,1) =  val_U_Rank{k}(indPair2,2);
                invval2_U_Rank{k} (indPair2,2) =  val_U_Rank{k}(indPair1,1);
                invval2_U_Tot{k} (indPair1,1) =  val_U_Tot{k}(indPair2,2);
                invval2_U_Tot{k} (indPair2,2) =  val_U_Tot{k}(indPair1,1);
                
            end
        end
    end
    
    
    
    if autoregresive_corr==1
        
        val_rat{k}=     [ val_rat{k}(1:end-1, :)];
        val_2AFC{k} =   [ val_2AFC{k}(1:end-1, :)];
        val_SignedConfAB{k} = [ val_SignedConfAB{k}(1:end-1, :)];
        
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
        asfd
        
        val_rat_Rat{k}= [ val_rat_Rat{k}(1:end-1, :)];
        val_rat_Rank{k}= [ val_rat_Rank{k}(1:end-1, :)];
        invval_rat_Rat{k}=invval_rat_Rat{k}(1:end-1, :);
        invval_rat_Rank{k}=invval_rat_Rank{k}(1:end-1, :);
        
        
        RTmean_Rat{k} = RTmean_Rat{k}(1:end-1, :);
        invRTmean_Rat{k} = invRTmean_Rat{k}(1:end-1, :);
    end
    
    
    k=k+1;
end




RT_AFC = shortRT_AFC;
RT_Rat = shortRT_Rat;
invRT_AFC = invshortRT_AFC;



for kk=1:numel(val_2AFC)
    sal_rat{kk} = abs(val_rat{kk});
    val_Confidence{kk} = abs(val_2AFC{kk}) ;
    
    val_AbsDiff{kk}  = abs(val_A_Rat{kk} -val_B_Rat{kk}) ;
    val_AbsDiffRank{kk}  = abs(val_A_Rank{kk} - val_B_Rank{kk});
   
    val_DiffABRat{kk}  = val_A_Rat{kk} -val_B_Rat{kk} ;
    val_DiffABRank{kk}  = val_A_Rank{kk} - val_B_Rank{kk} ;
      

    val_DiffLR{kk}  = val_L_Rat{kk} -val_R_Rat{kk} ;
    val_DiffLRRank{kk}  = val_L_Rank{kk} - val_R_Rank{kk} ;
    
    val_Sum{kk}  = val_A_Rat{kk} +val_B_Rat{kk} ;
    val_SumRank{kk}  = val_A_Rank{kk} + val_B_Rank{kk} ;
    
    val_ABbin{kk}  = double(val_SignedConfAB{kk}>0);
    val_LRbin{kk}  = double(val_2AFC{kk}>0);
    
    sal_A_Rat{kk}  = abs(val_A_Rat{kk});
    sal_B_Rat{kk}  = abs(val_B_Rat{kk});
    
    sal_rat_Rat{kk}  = abs(val_rat_Rat{kk});
    sal_rat{kk}  = abs(val_rat{kk});

    
    
    %% define inverses
    
    invval_Confidence{kk} = abs(invval_2AFC{kk}) ;
    
    invval_AbsDiff{kk}  = abs(invval_A_Rat{kk} -invval_B_Rat{kk}) ;
    invval_AbsDiffRank{kk}  = abs(invval_A_Rank{kk} - invval_B_Rank{kk});
    
    invval_DiffABRat{kk}  = invval_A_Rat{kk} -invval_B_Rat{kk} ;
    invval_DiffABRank{kk}  = invval_A_Rank{kk} - invval_B_Rank{kk} ;
        
    invval_DiffLR{kk}  = invval_L_Rat{kk} -invval_R_Rat{kk} ;
    invval_DiffLRRank{kk}  = invval_L_Rank{kk} - invval_R_Rank{kk} ;
 
    invval_Sum{kk}  = invval_A_Rat{kk} +invval_B_Rat{kk} ;
    invval_SumRank{kk}  = invval_A_Rank{kk} + invval_B_Rank{kk} ;
   
    invval_ABbin{kk}  = double(invval_SignedConfAB{kk}>0);
    invval_LRbin{kk}  = double(invval_2AFC{kk}>0);
       
    invsal_A_Rat{kk}  = abs(invval_A_Rat{kk});
    invsal_B_Rat{kk}  = abs(invval_B_Rat{kk});

    invsal_rat_Rat{kk}  = abs(invval_rat_Rat{kk});

    % new chosen inverses
    
    val2_C_Rat{kk} =  val_C_Rat{kk};
    val2_U_Rat{kk} =  val_U_Rat{kk} ;
    
    val2_C_Rank{kk} =  val_C_Rank{kk};
    val2_U_Rank{kk} =  val_U_Rank{kk} ;

       
end


%% supplementary plot  1, BIC and change in value


variablesAFCsess = {'val_A_Rat',  'sal_A_Rat', 'RTmean_Rat' ...
     'val_A_Rank', 'val_Confidence',  'RT_AFC'};

variablesAFCsessnames= {'Mean Rating',  'Unsigned Rating', 'RT Valuation', ...
     'Ranking', 'Confidence', 'RT 2AFC'};

 

clear Pcoef* Scoef*  pvalPrePostSinged  pvalPrePostRank


for ll = 1:size(variablesAFCsess,2)
    
    name1 = eval(variablesAFCsess{ll});
    name2 = eval(['inv', variablesAFCsess{ll}]);

    
    for k= 1:size(name1,2)
        name3 = name1{k}(:,1) + name1{k}(:,2);
        name4 = name2{k}(:,1) + name2{k}(:,2);
    end
    
    for k= 1:size(name1,2)
        
        
        [Scoefpre(ll,k), Scoefprepval(ll,k)]  = corr(name1{k}(:,1), name2{k}(:,1), 'type', 'S', 'rows', 'c');
        
    end
end





sizetrialtype =17;
distul = 2.0009 +0.4;
distlat= 0.2945 +0.2;
distconfchoice = 5.3213 ;

width =  8.3729 +0.4;
width12 = 10.2017 +0.4;
ypos=  18.6448 + 0.4;
height2 = 5.1452;
%%  plot BIC

heigth1= 6.8854;
xpos= 2.50 +.4; diffx1x2 = distconfchoice + width +.6;
xpos2= xpos + diffx1x2 ;


ypos2low= 1.6589;
ypos2  = ypos2low + height2 + distul; 
width2= 6.5611;


diffx1x2x3 =  distlat  + width2 +  2.0118 + .3;




figPos = [ 0   0    30.6800   28];
confidencemodel = [xpos  ypos    width  heigth1];
choicemodel = [xpos2     ypos   width12 heigth1];


yshifletra = 0.22; 
xshifletra = -2.1;


close all

figBic = figure('units', 'centimeters', 'Position' , figPos, 'color', 'w', 'visible', 'off');

% subplot('Position', [0.1300 0.1625 0.2747 0.7625])
sub1= subplot('Position', [0.1300 0.6625 0.2747  0.3012]);
set(sub1,  'units', 'centimeters', 'position', confidencemodel);

bicbar1 = bar(toplotBIC1);

ylimites = [10000, 12000];
set(gca, 'Ylim', ylimites);
set(gca, 'Xlim', [0.3,4.7]);

set(gca, 'YTick', linspace(ylimites(1), ylimites(2),5), 'YTicklabel', linspace(ylimites(1), ylimites(2),5))
bicbar1.FaceColor=[0.7,0.7,0.7];
ylabel('BIC')
set(gca, 'XTick', [1:7], 'XTicklabel',{ 'Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5', 'Model 6', 'Model 7'},'XTickLabelRotation', 45)

set(gca,'FontSize',12, 'FontWeight', 'bold')

title('Confidence', 'FontSize',sizetrialtype, 'FontWeight', 'bold')
box off





%subplot('Position', [0.5703 0.1625 0.3347 0.7625])
sub2 = subplot('Position', [0.5703 0.6625 0.3347 0.3012]);
set(sub2,  'units', 'centimeters', 'position', choicemodel);

bicbar=bar(toplotBIC2);
limites= [25500, 39500];

set(gca, 'Xlim', [0.3,5.7]);

set(gca, 'Ylim', limites);
set(gca, 'YTick', linspace(limites(1), limites(2),5), 'YTicklabel',linspace(limites(1), limites(2),5))
bicbar.FaceColor=[0.7,0.7,0.7];
ylabel('BIC')
set(gca, 'XTick', [1:size(toplotBIC2,2)], 'XTicklabel',{ 'Model 5', 'Model 6', 'Model 7', 'Model 8',  'Model 9', 'Model 10', 'Model 11'},'XTickLabelRotation', 45)




set(gca,'FontSize',12, 'FontWeight', 'bold')

title ('Choice', 'FontSize',sizetrialtype, 'FontWeight', 'bold')
box off




anot1 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',['A'],'Units','normalized','FontName','Helvetica','FontSize', sizeletra,'EdgeColor','none', 'fontweight', 'bold');
set(anot1, 'units', 'centimeters', 'position', [xpos+xshifletra,   ypos+heigth1+yshifletra,   2, 2] )



pos = [xpos,            ypos2      width2    height2;...
 xpos+ diffx1x2x3,       ypos2      width2    height2;...
xpos+(diffx1x2x3*2),    ypos2      width2    height2;...
xpos,                   ypos2low   width2    height2;...
xpos+(diffx1x2x3),      ypos2low   width2    height2;...
 xpos+(diffx1x2x3*2),   ypos2low   width2    height2];


for ll=1:6

    axes(  'units', 'centimeters', 'position', pos(ll,:));
    hold on 

    area([0.5, 12.5], [mean(Scoefpre(ll,:))+ std(Scoefpre(ll,:)), mean(Scoefpre(ll,:))+ std(Scoefpre(ll,:))],mean(Scoefpre(ll,:))- std(Scoefpre(ll,:)), 'Facecolor', [0.8,0.8,0.8],'facealpha',0.5, 'edgealpha',0.5, 'edgecolor', [0.8,0.8,0.8], 'ShowBaseLine', 'off')


    line([0.5, 12.5], [mean(Scoefpre(ll,:)),mean(Scoefpre(ll,:))], 'linestyle', '--', 'linewidth', 2, 'color', 'k')   

    plot(Scoefpre(ll,:), '.', 'markersize', 18)
    plot(8, Scoefpre(ll,8), '.', 'markersize', 18, 'color', 'r')
   
    xlim([0.45, 12.5])
       
    ylim([-0.65, 1.05])
    set(gca, 'Fontsize',12, 'fontweight', 'bold')
    
    if ll>=4
        xlabel('Patient ID')
    end
    if ll ==1 || ll==4
        ylabel('Rho')
    end
 title(variablesAFCsessnames{ll}, 'FontSize',sizetrialtype, 'FontWeight', 'bold')
end

anot2 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',['B'],'Units','normalized','FontName','Helvetica','FontSize', sizeletra,'EdgeColor','none', 'fontweight', 'bold');
set(anot2, 'units', 'centimeters', 'position', [xpos+xshifletra   ypos2+height2+yshifletra   2 2] )











toSave = [folder2save, '/Satiation by patient Spearman'];

print(toSave, '-dpng', '-r500');
print([toSave, '300'], '-dpng', '-r300');
print([toSave, '200'], '-dpng', '-r200');
print([toSave, 'Vec'], '-depsc', '-painters');
print([toSave, 'Vec'], '-depsc', '-opengl', '-r500');
set(figBic,'Units','centimeters');
pos = get(figBic,'Position');
set(figBic,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print([toSave, 'Vec'], '-dpdf', '-painters');
print(toSave, '-dpdf', '-r500', '-opengl');
print([toSave, ' 300'], '-dpdf', '-r300', '-opengl');




















%% plot BIC 2

sizetrialtype =18;

distul = 2.0009 +0.4;
distlat= 0.2945 +0.2;
distconfchoice = 5.3213 - 2.5 ;

width =  7.3729 +0.4;
width12 = (width/4)*5;
width13 = width/2;
ypos=  18.6448 + 0.4;
height2 = 5.1452 ;
%%  plot BIC

heigth1= 6.8854;
xpos= 2.50; 
diffx1x2 = distconfchoice + width;
diffx2x3 = distconfchoice + width12;
xpos2= xpos + diffx1x2;
xpos3= xpos2 + diffx2x3;

ypos2low= 1.6589;
ypos2  = ypos2low + height2 + distul; 
width2= 7.3;


diffx1x2x3 =  distlat  + width2 +  2.0118 + 0.0527;

% xpos+(diffx1x2x3*2) + width2 - (xpos + diffx1x2 + width12)








figPos = [ 0   0   30.6800   28];
confidencemodel = [xpos  ypos    width  heigth1];
choicemodel = [xpos2     ypos   width12 heigth1];
rtmodel = [xpos3     ypos   width13 heigth1];

yshifletra = 0.22; 
xshifletra = -2.1;


close all

figBic = figure('units', 'centimeters', 'Position' , figPos, 'color', 'w', 'visible', 'off');

% subplot('Position', [0.1300 0.1625 0.2747 0.7625])
sub1= subplot('Position', [0.1300 0.6625 0.2747  0.3012]);
set(sub1,  'units', 'centimeters', 'position', confidencemodel);

bicbar1 = bar(toplotBIC1);

ylimites = [10000, 12000];
set(gca, 'Ylim', ylimites);
set(gca, 'Xlim', [0.3,4.7]);

set(gca, 'YTick', linspace(ylimites(1), ylimites(2),5), 'YTicklabel', linspace(ylimites(1), ylimites(2),5))
bicbar1.FaceColor=[0.7,0.7,0.7];
ylabel('BIC')
set(gca, 'XTick', [1:7], 'XTicklabel',{ 'Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5', 'Model 6', 'Model 7'},'XTickLabelRotation', 45)

set(gca,'FontSize',12, 'FontWeight', 'bold')

title('Confidence', 'FontSize',sizetrialtype, 'FontWeight', 'bold')
box off





%subplot('Position', [0.5703 0.1625 0.3347 0.7625])
sub2 = subplot('Position', [0.5703 0.6625 0.3347 0.3012]);
set(sub2,  'units', 'centimeters', 'position', choicemodel);

bicbar=bar(toplotBIC2);
limites= [25500, 39500];

set(gca, 'Xlim', [0.3,5.7]);

set(gca, 'Ylim', limites);
set(gca, 'YTick', linspace(limites(1), limites(2),5), 'YTicklabel',linspace(limites(1), limites(2),5))
bicbar.FaceColor=[0.7,0.7,0.7];
%xlabel('BIC')
set(gca, 'XTick', [1:size(toplotBIC2,2)], 'XTicklabel',{ 'Model 5', 'Model 6', 'Model 7', 'Model 8',  'Model 9', 'Model 10', 'Model 11'},'XTickLabelRotation', 45)




set(gca,'FontSize',12, 'FontWeight', 'bold')

title ('Choice', 'FontSize',sizetrialtype, 'FontWeight', 'bold')
box off





%rt model
sub3 = axes('units', 'centimeters', 'position', rtmodel);
bicbar=bar(toplotBIC3);
limites= [4000, 4100];


set(gca, 'Xlim', [0.3,2.7]);
set(gca, 'Ylim', limites);
set(gca, 'YTick', linspace(limites(1), limites(2),5), 'YTicklabel',linspace(limites(1), limites(2),5))
bicbar.FaceColor=[0.7,0.7,0.7];
%xlabel('BIC')
set(gca, 'XTick', [1:size(toplotBIC3,2)], 'XTicklabel',{ 'Model 10', 'Model 11'},'XTickLabelRotation', 45)




set(gca,'FontSize',12, 'FontWeight', 'bold')
title ('RT', 'FontSize',sizetrialtype, 'FontWeight', 'bold')
box off



%(xpos+(diffx1x2x3*2) +width2) - (rtmodel(1) + rtmodel(3))




anot1 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',['A'],'Units','normalized','FontName','Helvetica','FontSize', sizeletra,'EdgeColor','none', 'fontweight', 'bold');
set(anot1, 'units', 'centimeters', 'position', [xpos+xshifletra,   ypos+heigth1+yshifletra,   2, 2] )



pos = [xpos,            ypos2      width2    height2;...
 xpos+ diffx1x2x3,       ypos2      width2    height2;...
xpos+(diffx1x2x3*2),    ypos2      width2    height2;...
xpos,                   ypos2low   width2    height2;...
xpos+(diffx1x2x3),      ypos2low   width2    height2;...
 xpos+(diffx1x2x3*2),   ypos2low   width2    height2];


for ll=1:6

    axes(  'units', 'centimeters', 'position', pos(ll,:));
    hold on 

    area([0.5, 12.5], [mean(Scoefpre(ll,:))+ std(Scoefpre(ll,:)), mean(Scoefpre(ll,:))+ std(Scoefpre(ll,:))],mean(Scoefpre(ll,:))- std(Scoefpre(ll,:)), 'Facecolor', [0.8,0.8,0.8],'facealpha',0.5, 'edgealpha',0.5, 'edgecolor', [0.8,0.8,0.8], 'ShowBaseLine', 'off')


    line([0.5, 12.5], [mean(Scoefpre(ll,:)),mean(Scoefpre(ll,:))], 'linestyle', '--', 'linewidth', 2, 'color', 'k')   

    plot(Scoefpre(ll,:), '.', 'markersize', 18)
    plot(8, Scoefpre(ll,8), '.', 'markersize', 18, 'color', 'r')
   
    xlim([0.45, 12.5])
       
    ylim([-0.65, 1.05])
    set(gca, 'Fontsize',12, 'fontweight', 'bold')
    
    if ll>=4
        xlabel('Patient ID')
    end
    if ll ==1 || ll==4
        ylabel('Rho')
    end
 title(variablesAFCsessnames{ll}, 'FontSize',sizetrialtype, 'FontWeight', 'bold')
end

anot2 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',['B'],'Units','normalized','FontName','Helvetica','FontSize', sizeletra,'EdgeColor','none', 'fontweight', 'bold');
set(anot2, 'units', 'centimeters', 'position', [xpos+xshifletra   ypos2+height2+yshifletra   2 2] )







toSave = [folder2save, '/Satiation by patient Spearman 3'];

print(toSave, '-dpng', '-r500');




print([toSave, '300'], '-dpng', '-r300');
print([toSave, '200'], '-dpng', '-r200');
print([toSave, 'Vec'], '-depsc', '-painters');
print([toSave, 'Vec'], '-depsc', '-opengl', '-r500');
set(figBic,'Units','centimeters');
pos = get(figBic,'Position');
set(figBic,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print([toSave, 'Vec'], '-dpdf', '-painters');
print(toSave, '-dpdf', '-r500', '-opengl');
print([toSave, ' 300'], '-dpdf', '-r300', '-opengl');





