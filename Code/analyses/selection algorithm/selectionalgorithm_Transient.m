function  selectionalgorithm_new_Transient2 (alpha,  regionIn,  direction, pooled, sorttype)


%  selectionalgorithm_new_Transient2(alpha,  regionIn, meanWindow, shift)
%  produce a selection algorithm for transient stimlus, including ony hte
%  stimulus presentations from onset to offset. do infernce with a post-hoc
%  binomial test. Results are saved as a table figure and a text document.
%  Input:
%           alpha : is the alpha value for the correlation coefficient, typically 0.05

%
%  sorttype: is the writting type by orther:
%	- 'none': no sorting, just the order they were analyzed
%	- 'prop': sorting by the highers mean proportion of cells
% selectionalgorithm_new_Transient2 (0.05,  {'ALL'; 'PHC'; 'EC'; 'H'; 'A'}, 1, 0, 'vertical', 1, 'prop')



dbstop if error
alpha2=0.05;
alphalimit = (alpha*2)^(1/2); % multiply the input alpha times 2 to correct for the sign correction made, and then sqrt for the pre post satiation condition

if ~(exist('pooled')==1)
    pooled = 1;
end

if ~(exist('sorttype')==1)
    sorttype = 'none';
end

cd('/media/Projects/Alex/Reclustered analysis All/mixed_model1_Spearman_window10_kernel200ms_step10ms__InversesVisualNutritional')
data_folder =  '/media/Projects/Alex/Reclustered analysis All/mixed_model1_Spearman_window10_kernel200ms_step10ms__InversesVisualNutritional/';
folder_to_Save = '/media/Projects/Alex/Reclustered analysis All/ Algorithm results and posthoc/'; mkdir(folder_to_Save);


if ~isequal(class(regionIn),'cell')
    error('input in regions is either not a cell, or the components are not a column vector');
end
regionIn= regionIn(:); % convert to column
if mean(ismember(regionIn,{ 'A';  'H';  'PHC'; 'EC';'ALL'; 'LA';  'RA'; 'LH'; 'LPHC'; 'LEC'; 'RH'; 'RPHC';'REC'; 'ConA';'ConH'; 'ConPHC'; 'ConEC'; 'FocA';'FocH'; 'FocPHC' ;'FocEC'; 'L'; 'R'; 'Con'; 'Foc'}))<1
    error('invalid region, check capitalization or the name')
end

if ~ismember(direction, ['vertical', 'horizontal'])
    error('invalid direction, check vertical or horizontal')
end



windowstep = 10;
windowsize = 0;
windowsizeRat = 0;
starttime=-5000; endtime=7400;
areastep=1200;
sizeAFCstim= 1000;
sizeRatstim= 2000;


clear neww*
newwinwosreg(:,1)=[starttime:windowstep:endtime-windowsize]; newwinwosreg(:,2)=newwinwosreg(:,1)+windowsize;
newwinwosregRating(:,1)=[starttime:windowstep:endtime-windowsizeRat]; newwinwosregRating(:,2)=newwinwosregRating(:,1)+windowsizeRat;
trialonsetRat=  find(newwinwosregRating(:,1)==0) - (windowsizeRat/windowstep)/2; %6;%
trialonset=     find(newwinwosreg(:,1)==0)-((windowsize/windowstep)/2); % 6; %

A1area = [trialonset,trialonset + sizeAFCstim/windowstep];
Barea = A1area + areastep/windowstep ;
A2area= Barea + areastep/windowstep;
Ratarea =  [trialonsetRat,trialonsetRat + sizeRatstim/windowstep];



namefile = [folder_to_Save,'/NEW Results TRANSIENT2 Mean Proportion of correlated cells, all windows and splitted windows, ',direction];

if pooled == 2
    namefile = [namefile, ', PreSat'];
elseif pooled == 3
    namefile = [namefile, ', PostSat'];
end
fileID = fopen(namefile, 'w');




changeINValue = '/media/Projects/Alex/Reclustered analysis All/mixed_model1_Spearman_window10_kernel200ms_step10ms__changeINvalue/';
%% multi-window  algorithm


arange1 =  501: 700;   % only stimulus time for Rat
arange2 =  501: 600;   % A1
arange3 =  621: 720;   % only for B
arange4 =  741: 840;   % only for A2



if exist('shift') ~= 1
    shift=0;
else
    shift = shift/windowstep;
end

if (shift*windowstep) == 4
    arange1 =  501: 720;   % only stimulus time for Rat
    arange2 =  501: 620;   % A1
    arange3 =  621: 740;   % only for B
    arange4 =  741: 860;   % only for A2
else
    arange1 =  arange1 + shift;   % only stimulus time for Rat
    arange2 =  arange2 + shift;   % A1
    arange3 =  arange3 + shift;   % only for B
    arange4 =  arange4 + shift;   % only for A2
end



toLoadA =  {'*Aval_A1',         '*Arank_A1',               '*Asal_A1',  ...
    '*flavourA_A1',              '*ASalz_A1',               '*AZucker_A1', ...
    '*ABallaststoffe_A1',        '*AProtein_A1',            '*AKohlenhydrate_A1',     ...
    '*AFett_A1',                 '*AKalorien_A1',           '*AHue_A1', ...
    '*ABlue_A1',                 '*ARed_A1',                '*AGreen_A1', ...
    '*AContrast_A1',             '*ALuminance_A1',          '*ASaturation_A1'};


toLoadB =  {'*Bval_A1',          '*Brank_A1',               '*Bsal_A1',  ...
    '*flavourB_A1',                 '*BSalz_A1',               '*BZucker_A1', ...
    '*BBallaststoffe_A1',           '*BProtein_A1',            '*BKohlenhydrate_A1',     ...
    '*BFett_A1',                    '*BKalorien_A1',           '*BHue_A1', ...
    '*BBlue_A1',                    '*BRed_A1',                '*BGreen_A1', ...
    '*BContrast_A1',                '*BLuminance_A1',          '*BSaturation_A1'};


toLoadR =  {'*Val_R',      '*Rank_R',            '*Sal_R', ...
    '*flavourRat_R',       '*Salz_R',            '*Zucker_R',...
    '*Ballaststoffe_R',    '*Protein_R',         '*Kohlenhydrate_R',   ...
    '*Fett_R',             '*Kalorien_R',        '*Hue_R', ...
    '*Blue_R',             '*Red_R',             '*Green_R',...
    '*Contrast_R',         '*Luminance_R',       '*Saturation_R'};







for region2=1:size(regionIn)
    
    clear pval* corr*
    region=char(regionIn{region2});
    
    
    
    
    %% Rating load process and save
    
    clear pval* corr*
    load([data_folder,'regresspval_',region ,'.mat'], toLoadR{:})
    load([data_folder,'regressreg_', region ,'.mat'], toLoadR{:})
    
    
    
    %     load change_invalue
    load([changeINValue, 'regresspval_',region,'.mat'], ['pvalRawchangeRat_', toLoadR{1}(2:end-2)], ['pvalRawchangeRat_', toLoadR{2}(2:end-2)], ['pvalRawchangeRat_', toLoadR{3}(2:end-2)]);
    load([changeINValue, 'regressreg_',region,'.mat'],  ['corrRawchangeRat_',toLoadR{1}(2:end-2)],  ['corrRawchangeRat_',toLoadR{2}(2:end-2)], ['corrRawchangeRat_',toLoadR{3}(2:end-2)]);
    changeNumR1(1,:) = nansum(pvalRawchangeRat_Val< alpha2, 1);
    changeNumR1(2,:) = nansum(pvalRawchangeRat_Rank< alpha2, 1);
    changeNumR1(3,:) = nansum(pvalRawchangeRat_Sal< alpha2, 1);
    
    
    
    
    for  yy= 1:numel(toLoadR)
        listpvalR{yy,1} = ['pval', toLoadR{yy}(2:end)];
    end
    
    
    % select cells for each variable
    for ff=1: size(listpvalR,1)
        
        clear newVal corrVal *pool  *pool *comp
        if pooled == 1
            
            clear newVal corrVal *pool  *pool *comp
            eval(['newVal=',listpvalR{ff},'< alphalimit;'])     % copy pvalues for workflow
            lopsize = size(newVal,1)/2;  % divided by 2 as pooled
            eval(['newCorr= corr',listpvalR{ff}(5:end),';']) % not asolute value % copy rhos for workflow; NOT abs value
            
            
            for ll=1:lopsize                                                               % significant in both and same sign of correlation
                pvalpool(ll,:) = newVal(ll,:) & newVal(ll+lopsize,:);                % cells significant BOTH pre and post satiation for this variable
                signpool(ll,:) = sign(newCorr(ll,:)) == sign(newCorr(ll+lopsize,:));  % find those with same sign for the correlation pre and post satiation
                corrpool(ll,:) = abs(newCorr(ll,:) + newCorr(ll+lopsize,:))/2;            % do the average correlation, its absolute value for comparison of correlation strength among variables
                corrpool(ll, isnan(corrpool(ll,:)))=0;                              % delete Nan; for cells with low Hz
            end
            
            pvalpool = pvalpool & signpool;  % both significant and same sign locations
            eval([listpvalR{ff},' =  pvalpool;'])  % force selected cells (significant both pre and post) to have the same sign
            eval(['corr',listpvalR{ff}(5:end),' =  abs(pvalpool .* corrpool);'])  % delete average-correlations of those with different sign for comparison and only for those significant BOTH pre and post
            clear newVal corrVal *pool
            
            
            
            eval([listpvalR{ff},' = ',listpvalR{ff},'(vectorPat,:);'])
            eval(['corr',listpvalR{ff}(5:end),' = corr', listpvalR{ff}(5:end),'(vectorPat,:);'])
            
            
            
            
        elseif pooled == 2
            
            clear newVal corrVal *pool  *pool *comp
            eval(['newVal=',listpvalR{ff},'< alpha;'])     % copy pvalues for workflow
            lopsize = size(newVal,1)/2;  % divided by 2 as pooled
            eval(['newCorr= corr',listpvalR{ff}(5:end),';']) % not asolute value   % copy rhos for workflow; NOT abs value
            
            eval([listpvalR{ff},' =  newVal(1:lopsize,:);'])  % force selected cells (significant both pre and post) to have the same sign
            eval(['corr',listpvalR{ff}(5:end),' =  abs(newCorr (1:lopsize,:));'])
            
        elseif pooled == 3
            
            clear newVal corrVal *pool  *pool *comp
            eval(['newVal=',listpvalR{ff},'< alpha;'])     % copy pvalues for workflow
            lopsize = size(newVal,1)/2;  % divided by 2 as pooled
            eval(['newCorr= corr',listpvalR{ff}(5:end),';']) % not asolute value   % copy rhos for workflow; NOT abs value
            
            eval([listpvalR{ff},' =  newVal(1+lopsize:end,:);'])  % force selected cells (significant both pre and post) to have the same sign
            eval(['corr',listpvalR{ff}(5:end),' =  abs(newCorr (1+lopsize:end,:));'])
            
        end
        
    end
    
    clear newVal corrVal *pool *comp pureBestRBinom  explainedfitsNumR bestfitsNumR explainedfitsNumA bestfitsNumA explainedfitsNumB bestfitsNumB
    
    
    
    
    
    % from here on is the algorithm itself
    
    explainedfitsR = zeros(size(listpvalR,1), size(eval(listpvalR{1}),2));
    bestfitsR = zeros(size(listpvalR,1), size(eval(listpvalR{1}),2));
    explainedNumR = zeros(size(listpvalR,1), size(eval(listpvalR{1}),2));
    bestNumR = zeros(size(listpvalR,1), size(eval(listpvalR{1}),2));
    pureBestR = zeros(size(listpvalR,1), size(eval(listpvalR{1}),2));
    pureNumR = zeros(size(listpvalR,1), size(eval(listpvalR{1}),2));
    
    
    
    
    
    for gg= 1: size(listpvalR,1)  % loop over all varaibles to do binary comparisons
        
        pureBestRBinom{gg} = zeros(size(listpvalR,1),  size(eval(listpvalR{gg}),2));
        pureBestRBinom2{gg} = zeros(size(listpvalR,1),  size(eval(listpvalR{gg}),2));
        pureNoOverlapRBinom{gg} = zeros(size(listpvalR,1),  size(eval(listpvalR{gg}),2));
        pureNoOverlapRBinom2{gg} = zeros(size(listpvalR,1),  size(eval(listpvalR{gg}),2));
        
        
        p1 = eval(listpvalR{gg}); %rename for easier reading
        c1 = eval(['corr',listpvalR{gg}(5:end)]);
        
        
        for bb= 1: size(listpvalR,1)
            if gg ~= bb %only if different variables
                
                
                p2 = eval(listpvalR{bb});
                c2 = eval(['corr',listpvalR{bb}(5:end)]);
                
                
                
                explainedfitsNumR{gg}(bb,:)= sum(p1,1) - sum(p2,1);
                bestfitsNumR{gg}(bb,:)    =  sum(p1 & (c1 > c2),1) - sum(p2 & (c2 >c1),1);
                pureBestRBinom{gg}(bb,:)    =  sum(p1 & (c1 > c2),1);
                pureBestRBinom2{gg}(bb,:)    =  sum(p2 & (c2 > c1),1);
                
                pureNoOverlapRBinom{gg}(bb,:)    =  sum(p1 & (~p2),1);
                pureNoOverlapRBinom2{gg}(bb,:)    =  sum(p2 & (~p1),1);
                
            end
            
        end
        
        explainedfitsR(gg,:)  = sum(explainedfitsNumR{gg}>0,1);
        bestfitsR(gg,:)       = sum(bestfitsNumR{gg}>0,1);
        explainedNumR(gg,:)   = sum(explainedfitsNumR{gg},1);
        bestNumR(gg,:)        = sum(bestfitsNumR{gg},1);
        pureBestR(gg,:)       = sum(pureBestRBinom{gg}, 1);
        pureNumR(gg,:)        = sum(p1 ,1);
        
    end
    
    
    
    
    
    
    
    
    
    
    %% A load process and save
    
    clear pval* corr*
    load([data_folder,'regresspval_',region ,'.mat'], toLoadA{:})
    load([data_folder,'regressreg_', region ,'.mat'], toLoadA{:})
    
    %     % load change_invalue
    %     load([changeINValue, 'regresspval_',region,'.mat'], ['pvalRawchangeA_', toLoadR{1}(2:end-2)], ['pvalRawchangeA_', toLoadR{2}(2:end-2)], ['pvalRawchangeA_', toLoadR{3}(2:end-2)]);
    %     load([changeINValue, 'regressreg_',region,'.mat'],  ['corrRawchangeA_',toLoadR{1}(2:end-2)],  ['corrRawchangeA_',toLoadR{2}(2:end-2)], ['corrRawchangeA_',toLoadR{3}(2:end-2)]);
    %     changeNumA1(1,:) = nansum(pvalRawchangeA_Val< alpha2, 1);
    %     changeNumA1(2,:) = nansum(pvalRawchangeA_Rank< alpha2, 1);
    %     changeNumA1(3,:) = nansum(pvalRawchangeA_Sal< alpha2, 1);
    %
    %
    
    
    for  yy= 1:numel(toLoadA)
        listpvalA{yy,1} = ['pval', toLoadA{yy}(2:end)];
    end
    
    
    % select cells for each variable
    for ff=1: size(listpvalA,1)
        
        clear newVal corrVal *pool  *pool *comp
        
        if pooled == 1
            clear newVal corrVal *pool  *pool *comp
            eval(['newVal=',listpvalA{ff},'< alphalimit;'])     % copy pvalues for workflow
            lopsize = size(newVal,1)/2;  % divided by 2 as pooled
            eval(['newCorr= corr',listpvalA{ff}(5:end),';']) % not asolute value % copy rhos for workflow; NOT abs value
            
            
            for ll=1:lopsize                                                               % significant in both and same sign of correlation
                pvalpool(ll,:) = newVal(ll,:) & newVal(ll+lopsize,:);                % cells significant BOTH pre and post satiation for this variable
                signpool(ll,:) = sign(newCorr(ll,:)) == sign(newCorr(ll+lopsize,:));  % find those with same sign for the correlation pre and post satiation
                corrpool(ll,:) = abs(newCorr(ll,:) + newCorr(ll+lopsize,:))/2;            % do the average correlation, its absolute value for comparison of correlation strength among variables
                corrpool(ll, isnan(corrpool(ll,:)))=0;                              % delete Nan; for cells with low Hz
            end
            
            pvalpool = pvalpool & signpool;  % both significant and same sign locations
            eval([listpvalA{ff},' =  pvalpool;'])  % force selected cells (significant both pre and post) to have the same sign
            eval(['corr',listpvalA{ff}(5:end),' =  abs(pvalpool .* corrpool);'])  % delete average-correlations of those with different sign for comparison and only for those significant BOTH pre and post
            clear newVal corrVal *pool
            
            
            eval([listpvalA{ff},' = ',listpvalA{ff},'(vectorPat,:);'])
            eval(['corr',listpvalA{ff}(5:end),' = corr', listpvalA{ff}(5:end),'(vectorPat,:);'])
            
            
            
            
        elseif pooled == 2
            
            clear newVal corrVal *pool  *pool *comp
            eval(['newVal=',listpvalA{ff},'< alpha;'])     % copy pvalues for workflow
            lopsize = size(newVal,1)/2;  % divided by 2 as pooled
            eval(['newCorr= corr',listpvalA{ff}(5:end),';']) % not asolute value   % copy rhos for workflow; NOT abs value
            
            eval([listpvalA{ff},' =  newVal(1:lopsize,:);'])  % force selected cells (significant both pre and post) to have the same sign
            eval(['corr',listpvalA{ff}(5:end),' =  abs(newCorr (1:lopsize,:));'])
            
        elseif pooled == 3
            
            clear newVal corrVal *pool  *pool *comp
            eval(['newVal=',listpvalA{ff},'< alpha;'])     % copy pvalues for workflow
            lopsize = size(newVal,1)/2;  % divided by 2 as pooled
            eval(['newCorr= corr',listpvalA{ff}(5:end),';']) % not asolute value   % copy rhos for workflow; NOT abs value
            
            eval([listpvalA{ff},' =  newVal(1+lopsize:end,:);'])  % force selected cells (significant both pre and post) to have the same sign
            eval(['corr',listpvalA{ff}(5:end),' =  abs(newCorr (1+lopsize:end,:));'])
            
        end
    end
    
    clear newVal corrVal *pool pureBestABinom explainedfitsNumA bestfitsNumA
    
    
    
    
    
    % from here on is the algorithm itself
    
    
    explainedfitsA = zeros(size(listpvalA,1), size(eval(listpvalA{1}),2));
    bestfitsA = zeros(size(listpvalA,1), size(eval(listpvalA{1}),2));
    explainedNumA = zeros(size(listpvalA,1), size(eval(listpvalA{1}),2));
    bestNumA = zeros(size(listpvalA,1), size(eval(listpvalA{1}),2));
    pureBestA = zeros(size(listpvalA,1), size(eval(listpvalA{1}),2));
    pureNumA = zeros(size(listpvalA,1), size(eval(listpvalA{1}),2));
    
    
    
    for gg= 1: size(listpvalA,1)  % loop over all varaibles to do binary comparisons
        
        pureBestABinom{gg} = zeros(size(listpvalA,1),  size(eval(listpvalA{gg}),2));
        pureBestABinom2{gg} = zeros(size(listpvalA,1),  size(eval(listpvalA{gg}),2));
        pureNoOverlapABinom{gg} = zeros(size(listpvalA,1),  size(eval(listpvalA{gg}),2));
        pureNoOverlapABinom2{gg} = zeros(size(listpvalA,1),  size(eval(listpvalA{gg}),2));
        
        
        clear p1 c1
        p1 = eval(listpvalA{gg}); %rename for easier reading
        c1 = eval(['corr',listpvalA{gg}(5:end)]);
        
        
        for bb= 1: size(listpvalA,1)
            if gg ~= bb %only if different variables
                
                clear p2 c2
                p2 = eval(listpvalA{bb});
                c2 = eval(['corr',listpvalA{bb}(5:end)]);
                
                
                
                explainedfitsNumA{gg}(bb,:) =  sum(p1,1) - sum(p2,1);
                bestfitsNumA{gg}(bb,:)      =  sum(p1 & (c1 > c2),1) - sum(p2 & (c2 >c1),1);
                pureBestABinom{gg}(bb,:)    =  sum(p1 & (c1 > c2),1);
                pureBestABinom2{gg}(bb,:)    =  sum(p2 & (c2 > c1),1);
                
                pureNoOverlapABinom{gg}(bb,:)    =  sum(p1 & (~p2),1);
                pureNoOverlapABinom2{gg}(bb,:)    =  sum(p2 & (~p1),1);
            end
            
        end
        
        explainedfitsA(gg,:)  = sum(explainedfitsNumA{gg}>0,1);
        bestfitsA(gg,:)       = sum(bestfitsNumA{gg}>0,1);
        explainedNumA(gg,:)   = sum(explainedfitsNumA{gg},1);
        bestNumA(gg,:)        = sum(bestfitsNumA{gg},1);
        pureBestA(gg,:)       = sum(pureBestABinom{gg}, 1);
        pureNumA(gg,:)        = sum(p1 ,1);
        
    end
    
    
    
    
    
    
    
    
    %% B load process and save
    
    clear pval* corr*
    load([data_folder,'regresspval_',region ,'.mat'], toLoadB{:})
    load([data_folder,'regressreg_', region ,'.mat'], toLoadB{:})
    
    
    
    
    
    %     % load change_invalue
    %     load([changeINValue, 'regresspval_',region,'.mat'], ['pvalRawchangeB_', toLoadR{1}(2:end-2)], ['pvalRawchangeB_', toLoadR{2}(2:end-2)], ['pvalRawchangeB_', toLoadR{3}(2:end-2)]);
    %     load([changeINValue, 'regressreg_',region,'.mat'],  ['corrRawchangeB_',toLoadR{1}(2:end-2)],  ['corrRawchangeB_',toLoadR{2}(2:end-2)],  ['corrRawchangeB_',toLoadR{3}(2:end-2)]);
    %     changeNumB1(1,:) = nansum(pvalRawchangeB_Val< alpha2, 1);
    %     changeNumB1(2,:) = nansum(pvalRawchangeB_Rank< alpha2, 1);
    %     changeNumB1(3,:) = nansum(pvalRawchangeB_Sal< alpha2, 1);
    
    
    for  yy= 1:numel(toLoadB)
        listpvalB{yy,1} = ['pval', toLoadB{yy}(2:end)];
    end
    
    
    % select cells for each variable
    for ff=1: size(listpvalB,1)
        
        clear newVal corrVal *pool  *pool *comp
        
        if pooled == 1
            clear newVal corrVal *pool  *pool *comp
            eval(['newVal=',listpvalB{ff},'< alphalimit;'])     % copy pvalues for workflow
            lopsize = size(newVal,1)/2;  % divided by 2 as pooled
            eval(['newCorr= corr',listpvalB{ff}(5:end),';']) % not asolute value % copy rhos for workflow; NOT abs value
            
            
            for ll=1:lopsize                                                               % significant in both and same sign of correlation
                pvalpool(ll,:) = newVal(ll,:) & newVal(ll+lopsize,:);                % cells significant BOTH pre and post satiation for this variable
                signpool(ll,:) = sign(newCorr(ll,:)) == sign(newCorr(ll+lopsize,:));  % find those with same sign for the correlation pre and post satiation
                corrpool(ll,:) = abs(newCorr(ll,:) + newCorr(ll+lopsize,:))/2;            % do the average correlation, its absolute value for comparison of correlation strength among variables
                corrpool(ll, isnan(corrpool(ll,:)))=0;                              % delete Nan; for cells with low Hz
            end
            
            pvalpool = pvalpool & signpool;  % both significant and same sign locations
            eval([listpvalB{ff},' =  pvalpool;'])  % force selected cells (significant both pre and post) to have the same sign
            eval(['corr',listpvalB{ff}(5:end),' =  abs(pvalpool .* corrpool);'])  % delete average-correlations of those with different sign for comparison and only for those significant BOTH pre and post
            clear newVal corrVal *pool
            
            
            if restricted > 0
                eval([listpvalB{ff},' = ',listpvalB{ff},'(vectorPat,:);'])
                eval(['corr',listpvalB{ff}(5:end),' = corr', listpvalB{ff}(5:end),'(vectorPat,:);'])
            end
            
            
            
            
        elseif pooled == 2
            
            clear newVal corrVal *pool  *pool *comp
            eval(['newVal=',listpvalB{ff},'< alpha;'])     % copy pvalues for workflow
            lopsize = size(newVal,1)/2;  % divided by 2 as pooled
            eval(['newCorr= corr',listpvalB{ff}(5:end),';']) % not asolute value   % copy rhos for workflow; NOT abs value
            
            eval([listpvalB{ff},' =  newVal(1:lopsize,:);'])  % force selected cells (significant both pre and post) to have the same sign
            eval(['corr',listpvalB{ff}(5:end),' =  abs(newCorr (1:lopsize,:));'])
            
        elseif pooled == 3
            
            clear newVal corrVal *pool  *pool *comp
            eval(['newVal=',listpvalB{ff},'< alpha;'])     % copy pvalues for workflow
            lopsize = size(newVal,1)/2;  % divided by 2 as pooled
            eval(['newCorr= corr',listpvalB{ff}(5:end),';']) % not asolute value   % copy rhos for workflow; NOT abs value
            
            eval([listpvalB{ff},' =  newVal(1+lopsize:end,:);'])  % force selected cells (significant both pre and post) to have the same sign
            eval(['corr',listpvalB{ff}(5:end),' =  abs(newCorr (1+lopsize:end,:));'])
            
        end
        
    end
    
    clear newVal corrVal *pool pureBestBBinom explainedfitsNumB bestfitsNumB
    
    
    
    
    
    % from here on is the algorithm itself
    
    explainedfitsB = zeros(size(listpvalB,1), size(eval(listpvalB{1}),2));
    bestfitsB = zeros(size(listpvalB,1), size(eval(listpvalB{1}),2));
    explainedNumB = zeros(size(listpvalB,1), size(eval(listpvalB{1}),2));
    bestNumB = zeros(size(listpvalB,1), size(eval(listpvalB{1}),2));
    pureBestB = zeros(size(listpvalB,1), size(eval(listpvalB{1}),2));
    pureNumB = zeros(size(listpvalB,1), size(eval(listpvalB{1}),2));
    
    
    
    
    for gg= 1: size(listpvalB,1)  % loop over all varaibles to do binary comparisons
        
        pureBestBBinom{gg} = zeros(size(listpvalB,1),  size(eval(listpvalB{gg}),2));
        pureBestBBinom2{gg} = zeros(size(listpvalB,1),  size(eval(listpvalB{gg}),2));
        pureNoOverlapBBinom{gg} = zeros(size(listpvalB,1),  size(eval(listpvalB{gg}),2));
        pureNoOverlapBBinom2{gg} = zeros(size(listpvalB,1),  size(eval(listpvalB{gg}),2));
        
        
        
        p1 = eval(listpvalB{gg}); %rename for easier reading
        c1 = eval(['corr',listpvalB{gg}(5:end)]);
        
        
        for bb= 1: size(listpvalB,1)
            if gg ~= bb %only if different variables
                
                
                p2 = eval(listpvalB{bb});
                c2 = eval(['corr',listpvalB{bb}(5:end)]);
                
                
                explainedfitsNumB{gg}(bb,:) = sum(p1,1) - sum(p2,1);
                bestfitsNumB{gg}(bb,:)      =  sum(p1 & (c1 > c2),1) - sum(p2 & (c2 >c1),1);
                pureBestBBinom{gg}(bb,:)    =  sum(p1 & (c1 > c2),1);
                pureBestBBinom2{gg}(bb,:)    =  sum(p2 & (c2 > c1),1);
                
                pureNoOverlapBBinom{gg}(bb,:)    =  sum(p1 & (~p2),1);
                pureNoOverlapBBinom2{gg}(bb,:)    =  sum(p2 & (~p1),1);
                
            end
        end
        
        explainedfitsB(gg,:)  = sum(explainedfitsNumB{gg}>0,1);
        bestfitsB(gg,:)       = sum(bestfitsNumB{gg}>0,1);
        explainedNumB(gg,:)   = sum(explainedfitsNumB{gg},1);
        bestNumB(gg,:)        = sum(bestfitsNumB{gg},1);
        pureBestB(gg,:)       = sum(pureBestBBinom{gg}, 1);
        pureNumB(gg,:)        = sum(p1 ,1);
        
    end
    
    
    
    
    %% select windows to analyze
    
    
    explainedNumR1      = nansum(explainedNumR(:,arange1), 2);
    explainedNumA1      = nansum(explainedNumA(:,arange2), 2);
    explainedNumB1      = nansum(explainedNumB(:,arange3), 2);
    explainedNumA2      = nansum(explainedNumA(:,arange4), 2);
    
    
    
    bestNumR1      = nansum(bestNumR(:,arange1), 2);
    bestNumA1      = nansum(bestNumA(:,arange2), 2);
    bestNumB1      = nansum(bestNumB(:,arange3), 2);
    bestNumA2      = nansum(bestNumA(:,arange4), 2);
    
    
    
    
    
    %% change in value sum cells
    
    pureNumR1      =  nanmean(pureNumR(:,arange1), 2);
    pureNumA1      =  nanmean(pureNumA(:,arange2), 2);
    pureNumB1      =  nanmean(pureNumB(:,arange3), 2);
    pureNumA2      =  nanmean(pureNumA(:,arange4), 2);
    
    pureBestR1      =  nanmean(pureBestR(:,arange1), 2);
    pureBestA1      =  nanmean(pureBestA(:,arange2), 2);
    pureBestB1      =  nanmean(pureBestB(:,arange3), 2);
    pureBestA2      =  nanmean(pureBestA(:,arange4), 2);
    
    
    
    
    
    %% add up the different windows
    
    sumExp1 =  explainedNumR1 + explainedNumA1 + explainedNumB1 + explainedNumA2;
    sumBestFits1 =  bestNumR1 + bestNumA1 + bestNumB1 + bestNumA2;
    sumpureNum1 = sumpureNum1/4;
    sumpureBest1 = sumpureBest1/4;
    
    if isequal(sorttype, 'prop')
        clear sortedfitsPureNum*
        numberfitsPureNum1{region2} = sumpureNum1/ lopsize;
        [~, sortedfitsPureNum1{region2}] = sort(sumpureNum1,'descend');
        vec2for = sortedfitsPureNum1{region2}';
    else
        
        vec2for = 1:numel(toLoadR);
    end
    
    
    for bb= vec2for
        fprintf(fileID, ['region,   ', region, '\n']);
        
        fprintf(fileID, ['%3s  %14s, TOTAL    percentage   %8.2f,   p-value binomial test, no correction   %11.5G, \n'], region, toLoadR{bb}(2:end-2), round(100*(floor(sumpureNum1(bb))/lopsize),2), 1-binocdf(sumpureNum1(bb)-1, lopsize, 0.05));
        fprintf(fileID, ['%3s  %14s, RATING   percentage   %8.2f,   p-value binomial test, no correction   %11.5G, \n'], region, toLoadR{bb}(2:end-2), round(100*(floor(pureNumR1(bb))/lopsize),2), 1-binocdf(pureNumR1(bb)-1, lopsize, 0.05));
        fprintf(fileID, ['%3s  %14s, A1       percentage   %8.2f,   p-value binomial test, no correction   %11.5G, \n'], region, toLoadR{bb}(2:end-2), round(100*(floor(pureNumA1(bb))/lopsize),2), 1-binocdf(pureNumA1(bb)-1, lopsize, 0.05));
        fprintf(fileID, ['%3s  %14s, B        percentage   %8.2f,   p-value binomial test, no correction   %11.5G, \n'], region, toLoadR{bb}(2:end-2), round(100*(floor(pureNumB1(bb))/lopsize),2), 1-binocdf(pureNumB1(bb)-1, lopsize, 0.05));
        fprintf(fileID, ['%3s  %14s, A2       percentage   %8.2f,   p-value binomial test, no correction   %11.5G, \n'], region, toLoadR{bb}(2:end-2), round(100*(floor(pureNumA2(bb))/lopsize),2), 1-binocdf(pureNumA2(bb)-1, lopsize, 0.05));
      
        fprintf(fileID, ['\n\n']);
    end
    
    
    
    % sort the variables
    clear sortedfitsPureNum*
    numberfitsPureNum1{region2} = sumpureNum1;
    [~, sortedfitsPureNum1{region2}] = sort(sumpureNum1,'descend');
    
    
    
    
    %% post-hoc test, do binomial tests and resort with algorithm

    
    binoBest1  = nan(size(listpvalB,1));
    binoNoOverlap = nan(size(listpvalB,1));
    

    for gg= 1: size(listpvalB,1)
        for bb= 1: size(listpvalB,1)
            if gg ~= bb
                
                p1 = floor((sum(pureBestRBinom{gg}(bb,arange1))  /numel(arange1)   + sum(pureBestABinom{gg}(bb,arange2)) /numel(arange2) +...
                    sum(pureBestBBinom{gg}(bb,arange3)) /numel(arange3)   + sum(pureBestABinom{gg}(bb,arange4)) /numel(arange4)) /4);
                p2 = floor((sum(pureBestRBinom2{gg}(bb,arange1)) /numel(arange1) + sum(pureBestABinom2{gg}(bb,arange2))  /numel(arange2) +...
                    sum(pureBestBBinom2{gg}(bb,arange3)) /numel(arange3)  + sum(pureBestABinom2{gg}(bb,arange4)) /numel(arange4)) /4);
                
                binoBest1(gg,bb) = 1- binocdf(p1-1 , p1+p2, 0.5);
                
                
                p1 = floor((sum(pureNoOverlapRBinom{gg}(bb,arange1))  /numel(arange1)  + sum(pureNoOverlapABinom{gg}(bb,arange2)) /numel(arange2) ...
                    + sum(pureNoOverlapBBinom{gg}(bb,arange3)) /numel(arange3)  + sum(pureNoOverlapABinom{gg}(bb,arange4)) /numel(arange4)) /4);
                p2 = floor((sum(pureNoOverlapRBinom2{gg}(bb,arange1)) /numel(arange1) + sum(pureNoOverlapABinom2{gg}(bb,arange2))  /numel(arange2) ...
                    + sum(pureNoOverlapBBinom2{gg}(bb,arange3)) /numel(arange3) + sum(pureNoOverlapABinom2{gg}(bb,arange4))  /numel(arange4)) /4);
                
                binoNoOverlap(gg,bb) = 1- binocdf(p1-1 , p1+p2, 0.5);
                
                
                
            end
        end
    end
    
    %% floor all
    sumpureNum1 = floor(sumpureNum1);
    
    BestNum1{region2} = nan(size(listpvalB,1));
    BestBest1{region2} = nan(size(listpvalB,1));
        
    kk = 1;
    for gg = 1: size(listpvalB,1)    % loop over all varaibles to do binary comparisons
        
        mm = 1;
        for bb = 1: size(listpvalB,1)
            
            if gg ~= bb    % only if different variables
                
                if kk == mm
                    mm = mm + 1;
                end                
                BestNum1{region2}(kk,mm) =  1 - binocdf(sumpureNum1(gg)-1, floor(sumpureNum1(gg) + sumpureNum1(bb)), 0.5);
                BestBest1{region2}(kk,mm) = binoBest1(gg,bb);                
                mm = mm +1;                
            end
            
        end
        kk = kk +1;
    end
    
   
   
    
    
    %%  save in a textfile the results
    
    searchind = 1;
    
    insert = [];
    inserts=[];
    
    for ff= 1:6
        insert = [insert, '%16.5f,          '];
        inserts = [inserts, '%16s,          '];
    end
    
    meansorted = [];
    for ff= 1:size(toLoadR,2)
        meansorted= [meansorted, ',   ',  [toLoadR{1,sortedfitsPureNum1{region2}(ff)}(2:end)] ];
    end
    
    
    fprintf(fileID,['\n\n\n sorted by higher mean \n']);
    fprintf(fileID, meansorted);
    
        fprintf(fileID,['\n\n Post-Hoc Binomal Test \n']);
    
    if strcmp(direction, 'horizontal')
        fprintf(fileID,['\n sorted in horizontal, so the one in row bigger than the one in column \n']);
        BestBest1{region2} = BestBest1{region2}';
    elseif strcmp(direction, 'vertical')
        fprintf(fileID,['\n sorted in vertical, so the one in row bigger than the one in column \n']);
    end
    
    while searchind< size(toLoadR,2)
             
        fprintf(fileID,['%21s', inserts, '\n' ], '',toLoadR{searchind:searchind+5});
        for ff= 1:size(BestBest1{region2})
            fprintf(fileID,[ '%16s%5s', insert, '\n'],toLoadR{ff},'' ,  BestBest1{region2}(searchind:searchind+5,ff));
        end
        
        fprintf(fileID, ['\n\n\n\n\n']);
        searchind = searchind+6;
    end
    
    
    
    fprintf(fileID, ['\n\n\n\n\n\n\n\n\n\n\n\n']);
    
end


fclose(fileID);
