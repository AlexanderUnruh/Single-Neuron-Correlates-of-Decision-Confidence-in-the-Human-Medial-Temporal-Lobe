function selectionalgorithm_new_Persistent_auto20ms(alpha,  regionIn,  direction)


%% This code performs the selection algorithm to compare binary variables for the persistent neural activity
%

dbstop if error
cd '/media/Projects/Alex/Reclustered analysis All 2/20s_1_Spearman_window10_kernel200ms_step10ms__ConfidenceRT/'
fold2save2 = '/media/Projects/Alex/Reclustered analysis All/ Algorithm results and posthoc/';
mkdir(fold2save2)

for pooled = 1:3
    
    clearvars -except fold2save2 alpha  regionIn  direction pooled
    namefile =  [fold2save2, 'NEW Results PERSISTENT Mean Proportion of correlated cells, all windows and splitted windows, ',direction, 'after likert +2.35 (classical) seconds, aplha = ', num2str(alpha),'.txt'];
    if ~(exist('pooled')==1)
        pooled = 1;
    end
    
    if ~isequal(class(regionIn),'cell')
        error('input in regions is either not a cell, or the components are not a column vector');
    end
    regionIn= regionIn(:);
    if mean(ismember(regionIn,{ 'A';  'H';  'PHC'; 'EC';'ALL'; 'LA';  'RA'; 'LH'; 'LPHC'; 'LEC'; 'RH'; 'RPHC';'REC'; 'ConA';'ConH'; 'ConPHC'; 'ConEC'; 'FocA';'FocH'; 'FocPHC' ;'FocEC'; 'L'; 'R'; 'Con'; 'Foc'}))<1
        error('invalid region, check capitalization or the name')
    end
    
    if ~ismember(direction, ['vertical', 'horizontal'])
        error('invalid direction, check vertical or horizontal')
    end
    
    
    
    cd('/media/Projects/Alex/Reclustered analysis All 2/20s_1_Spearman_window10_kernel200ms_step10ms__ConfidenceRT/')
    
    
    alphalimit = (alpha*2)^(1/2); % mul;tiply the input alpha times 2 to correct for the sign correction made, and then sqrt for the pre post satiation condition
    
    if pooled==2
        namefile = [namefile, ', PreSat'];
    elseif pooled ==3
        namefile = [namefile, ', PostSat'];
    end
    
    
    fileID = fopen(namefile, 'w');
    
    
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
    
    
    data_folder = '/media/Projects/Alex/Reclustered analysis All 2/20s_1_Spearman_window10_kernel200ms_step10ms__ConfidenceRT/';
    toLoad =  {'*RT_A1',   '*Abspref_A1',  ...
        '*absdiff_A1', '*Rankabsdiff_A1', ...
        '*Cval2_A1',  '*Crank2_A1', '*Uval2_A1', '*Urank2_A1', ...
        '*Sum_A1', '*Ranksum_A1', ...
        '*prefAB_A1',  '*diffAB_A1', '*RankABdiff_A1', ...
        '*ABbin_A1', '*LRbin_A1'...
        '*prefLR_A1',  '*diffLR_A1', '*RankLRdiff_A1', ...
        '*UCratio2_A1', '*RankCUratio2_A1', ...
        '*Lval_A1',  '*Rval_A1',  '*Lrank_A1', '*Rrank_A1', '*Lsal_A1', '*Rsal_A1', ...
        '*Aval_A1',  '*Bval_A1',  '*Arank_A1', '*Brank_A1', '*Asal_A1', '*Bsal_A1',   ...
        '*flavourA_A1', '*flavourB_A1',...
        '*flavourC_A1', '*flavourU_A1', ...
        '*flavourL_A1', '*flavourR_A1',...
        '*ASalz_A1', '*BSalz_A1', '*AZucker_A1', '*BZucker_A1', '*AProtein_A1', '*BProtein_A1', '*AKohlenhydrate_A1', '*BKohlenhydrate_A1', '*AFett_A1', '*BFett_A1', '*AKalorien_A1', '*BKalorien_A1', '*ABallaststoffe_A1', '*BBallaststoffe_A1', ...
        '*AHue_A1', '*BHue_A1', '*ABlue_A1', '*BBlue_A1', '*ARed_A1', '*BRed_A1', '*AGreen_A1', '*BGreen_A1', '*AContrast_A1', '*BContrast_A1', '*ALuminance_A1', '*BLuminance_A1',  '*ASaturation_A1', '*BSaturation_A1',...
        '*CSalz_A1', '*USalz_A1', '*CZucker_A1', '*UZucker_A1', '*CProtein_A1', '*UProtein_A1', '*CKohlenhydrate_A1', '*UKohlenhydrate_A1', '*CFett_A1', '*UFett_A1', '*CKalorien_A1', '*UKalorien_A1', '*CBallaststoffe_A1', '*UBallaststoffe_A1', ...
        '*CHue_A1', '*UHue_A1', '*CBlue_A1', '*UBlue_A1', '*CRed_A1', '*URed_A1', '*CGreen_A1', '*UGreen_A1', '*CContrast_A1', '*UContrast_A1', '*CLuminance_A1', '*ULuminance_A1',  '*CSaturation_A1', '*USaturation_A1',...
        '*LSalz_A1', '*RSalz_A1', '*LZucker_A1', '*RZucker_A1', '*LProtein_A1', '*RProtein_A1', '*LKohlenhydrate_A1', '*RKohlenhydrate_A1', '*LFett_A1', '*RFett_A1', '*LKalorien_A1', '*RKalorien_A1', '*LBallaststoffe_A1', '*RBallaststoffe_A1', ...
        '*LHue_A1', '*RHue_A1', '*LBlue_A1', '*RBlue_A1', '*LRed_A1', '*RRed_A1', '*LGreen_A1', '*RGreen_A1', '*LContrast_A1', '*RContrast_A1', '*LLuminance_A1', '*RLuminance_A1',  '*LSaturation_A1', '*RSaturation_A1',...
        };
    
    
    
    %% multi-window  algorithm
    window_of_analysis = [2361: 2361+235];
    
    for region2 = 1:size(regionIn)
        
        clear pval* corr*
        region = char(regionIn{region2});
        
        load([data_folder,'regresspval_',region ,'.mat'], toLoad{:})
        load([data_folder,'regressreg_', region ,'.mat'], toLoad{:})
        
        for  yy= 1:numel(toLoad)
            listpval{yy,1} = ['pval', toLoad{yy}(2:end)];
        end
        
        for ff=1: size(listpval,1)
            
            if pooled == 1
                eval(['newVal=',listpval{ff},'< alphalimit;'])     % copy pvalues for workflow
                lopsize = size(newVal,1)/2;  % divided by 2 as pooled
                eval(['newCorr= corr',listpval{ff}(5:end),';']) % not asolute value   % copy rhos for workflow; NOT abs value
                
                
                % process d ata
                clear *pool *comp
                
                % significant in both and same sign of correlation
                for ll=1:lopsize
                    pvalpool(ll,:) = (newVal(ll,:) & newVal(ll+lopsize,:));                % cells significant BOTH pre and post satiation for this variable
                    signpool(ll,:) = (sign(newCorr(ll,:)) == sign(newCorr(ll+lopsize,:)));  % find those with same sign for the correlation pre and post satiation
                    corrpool(ll,:) = abs(newCorr(ll,:) + newCorr(ll+lopsize,:))/2;            % do the average correlation, its absolute value for comparison of correlation strength among variables
                    corrpool(ll, isnan(corrpool(ll,:)))=0;                              % delete Nan; for cells with low Hz
                end
                pvalpool = pvalpool & signpool;
                eval([listpval{ff},' =  pvalpool;'])  % force selected cells (significant both pre and post) to have the same sign
                eval(['corr',listpval{ff}(5:end),' =  abs(pvalpool .* corrpool);'])  % delete average-correlations of those with different sign for comparison and only for those significant BOTH pre and post
                clear newVal corrVal *pool
                
            elseif pooled == 2
                
                eval(['newVal=',listpval{ff},'< alpha;'])     % copy pvalues for workflow
                lopsize = size(newVal,1)/2;  % divided by 2 as pooled
                eval(['newCorr= corr',listpval{ff}(5:end),';']) % not asolute value   % copy rhos for workflow; NOT abs value
                
                eval([listpval{ff},' =  newVal(1:lopsize,:);'])  % force selected cells (significant both pre and post) to have the same sign
                eval(['corr',listpval{ff}(5:end),' =  abs(newCorr (1:lopsize,:));'])
                
            elseif pooled == 3
                
                eval(['newVal=',listpval{ff},'< alpha;'])     % copy pvalues for workflow
                lopsize = size(newVal,1)/2;  % divided by 2 as pooled
                eval(['newCorr= corr',listpval{ff}(5:end),';']) % not asolute value   % copy rhos for workflow; NOT abs value
                
                eval([listpval{ff},' =  newVal(1+lopsize:end,:);'])  % force selected cells (significant both pre and post) to have the same sign
                eval(['corr',listpval{ff}(5:end),' =  abs(newCorr (1+lopsize:end,:));'])
                
            end
            
            
            eval([listpval{ff},' = ',listpval{ff},'(vectorPat,:);'])
            eval(['corr',listpval{ff}(5:end),' = corr', listpval{ff}(5:end),'(vectorPat,:);'])
            
            
            clear newVal corrVal *pool
            
        end
        
        
        
        
        % from here on is the algorithm itself
        
        explainedfits = zeros(size(listpval,1), size(eval(listpval{1}),2));
        bestfits = zeros(size(listpval,1), size(eval(listpval{1}),2));
        explainedNum = zeros(size(listpval,1), size(eval(listpval{1}),2));
        bestNum = zeros(size(listpval,1), size(eval(listpval{1}),2));
        pureBest = zeros(size(listpval,1), size(eval(listpval{1}),2));
        pureNum = zeros(size(listpval,1), size(eval(listpval{1}),2));
        
        
        
        for gg= 1: size(listpval,1)  % loop over all varaibles to do binary comparisons
            
            pureBestBinom{gg} = zeros(size(listpval,1),  size(eval(listpval{gg}),2));
            pureBestBinom2{gg} = zeros(size(listpval,1),  size(eval(listpval{gg}),2));
            pureNoOverlapBinom{gg} = zeros(size(listpval,1),  size(eval(listpval{gg}),2));
            pureNoOverlapBinom2{gg} = zeros(size(listpval,1),  size(eval(listpval{gg}),2));
            
            
            p1= eval(listpval{gg}); %rename for easier reading
            c1= eval(['corr',listpval{gg}(5:end)]);
            
            for bb= 1: size(listpval,1)
                if gg ~= bb %only if different variables
                    
                    
                    p2= eval(listpval{bb});
                    c2= eval(['corr',listpval{bb}(5:end)]);
                    
                    explainedfitsNum{gg}(bb,:)= sum(p1,1) - sum(p2,1);
                    bestfitsNum{gg}(bb,:)    =  sum(p1 & (c1 > c2),1) - sum(p2 & (c2 >c1),1);
                    bestNumBinom{gg}(bb,:)    =  sum(p1 & (c1 > c2),1);
                    
                    pureBestBinom{gg}(bb,:)    =  sum(p1 & (c1 > c2),1);
                    pureBestBinom2{gg}(bb,:)    =  sum(p2 & (c2 > c1),1);
                    pureNoOverlapBinom{gg}(bb,:)    =  sum(p1 & (~p2),1);
                    pureNoOverlapBinom2{gg}(bb,:)    =  sum(p2 & (~p1),1);
                    
                    
                else
                    %                   eval(['explainedfits', num2str(gg),'(bb,:)= nan(size(p1,2),1);'])
                    %                   eval(['bestfits', num2str(gg),'(bb,:)= nan(size(p1,2),1);'])
                end
                
            end
            
            
            explainedfits(gg,:)  =sum(explainedfitsNum{gg}>0,1);
            bestfits(gg,:)       =sum(bestfitsNum{gg}>0,1);
            explainedNum(gg,:)   =sum(explainedfitsNum{gg},1);
            bestNum(gg,:)       =sum(bestfitsNum{gg},1);
            
            pureBest(gg,:)       = sum(pureBestBinom{gg}, 1);
            pureNum(gg,:)        = sum(p1 ,1);
            pureNoOverlap(gg,:)  = sum(pureNoOverlapBinom{gg}, 1);
            
        end
        
        clear pval* corr*
        
        
        %% sum the windows
        
        explainedNum1 = nansum(explainedNum(:,window_of_analysis), 2);
        bestNum1      = nansum(bestNum(:,window_of_analysis), 2);
        pureNum1      =  nanmean(pureNum(:,window_of_analysis), 2);
        pureBest1     =  nanmean(pureBest(:,window_of_analysis), 2);
        
        %% add up the different windows
        
        for bb= 1:numel(toLoad)
            fprintf(fileID, ['region,   ', region, '\n']);
            fprintf(fileID, ['%5s  %18s, Persistent window, algorithm    percentage  %8.4f,   p-value binomial test, no correction  %11.5G, \n'], region, toLoad{bb}, round(100*(floor(pureNum1(bb))/lopsize),2), 1-binocdf(pureNum1(bb)-1, lopsize, 0.05));
            fprintf(fileID, ['\n\n']);
        end
        
        
        clear sortedfitsPureNum*
        numberfitsPureNum1{region2} = pureNum1;
        [~, sortedfitsPureNum1{region2}] = sort(pureNum1,'descend');
        
        
        
        %% Binomial test
        
        %% post-hoc test, do binomial tests and resort with algorithm
        
        
        binoBest1  = nan(size(listpval,1));
        binoNoOverlap1 = nan(size(listpval,1));
        
        clear p11 p111 p22 p222
        %% floor all
        
        pureNum1 = floor(pureNum1);
        
        for gg = 1: size(listpval,1)
            for bb = 1:size(listpval,1)
                
                if gg ~= bb
                    
                    
                    p1 = floor(mean(pureBestBinom{gg}(bb,window_of_analysis))) ;
                    p2 = floor(mean(pureBestBinom2{gg}(bb,window_of_analysis))) ;
                    p11 (gg,bb) = p1;
                    p22 (gg,bb) = p2;
                    binoBest1(gg,bb) = 1- binocdf(floor(p1)-1 , floor(p1)+floor(p2), 0.5);
                    
                    
                    p1 = floor(mean(pureNoOverlapBinom{gg}(bb,window_of_analysis)));
                    p2 = floor(mean(pureNoOverlapBinom2{gg}(bb,window_of_analysis)));
                    p111 (gg,bb) = p1; % number of no overlap
                    p222 (gg,bb) = p2; % number of no overlap
                    binoNoOverlap1(gg,bb) = 1- binocdf(floor(p1)-1 , floor(p1)+floor(p2), 0.5);
                    
                end
            end
        end
        
        
        %% display the direct comparisons between RT and confidence against the best alteranive models
        disp(['p11, big variable vs small variable'])
        p11(1:6,1:6)
        disp(['p22, small variable vs big variable'])
        p22(1:6,1:6)
        disp(['p111, noOverlap big variable vs small variable'])
        p111(1:6,1:6)
        disp(['p222, noOverlap small variable vs big variable'])
        p222(1:6,1:6)
        
        
        %% do binomial test without algorithm
        
        BestNum1{region2} = nan(size(listpval,1));
        BestBest1{region2} = nan(size(listpval,1));
        
        kk = 1;
        for gg = 1: size(listpval,1)    % loop over all varaibles to do binary comparisons
            
            mm = 1;
            for bb = 1: size(listpval,1)
                
                if gg ~= bb    % only if different variables
                    
                    if kk == mm
                        mm = mm + 1;
                    end
                    
                    BestNum1{region2}(kk,mm) =  1 - binocdf(pureNum1(gg)-1, floor(pureNum1(gg) + pureNum1(bb)), 0.5);
                    BestBest1{region2}(kk,mm) = binoBest1(gg,bb);
                    mm = mm +1;
                    
                end
            end
            kk = kk +1;
        end
        
        
        
        % explained algorithm
        
        kk = 1;
        for gg = (sortedfitsPureNum1{region2})'    % loop over all varaibles to do binary comparisons
            
            mm = 1;
            for bb = (sortedfitsPureNum1{region2})'
                
                if gg ~= bb    % only if different variables
                    
                    if kk == mm
                        mm = mm + 1;
                    end
                    
                    postPNNum1{region2}(kk,mm) = 1 - binocdf(pureNum1(gg)-1, floor(pureNum1(gg) + pureNum1(bb)), 0.5);
                    %postPNBest1{region2}(kk,mm) = 1 - binocdf(sumpureBest1(gg)-1, sumpureBest1(gg) + sumpureBest1(bb), 0.5);
                    postPNBest1{region2}(kk,mm) =  binoBest1(gg,bb);
                    
                    mm = mm +1;
                end
            end
            kk = kk +1;
        end
        
        
        
        
        
        
        %% transform to string to pllot
        %%  save in a textfile the results
        
        searchind = 1;
        insert = [];
        inserts=[];
        
        for ff= 1:6
            insert = [insert, '%18.5G,          '];
            inserts = [inserts, '%18s,          '];
        end
        
        meansorted = [];
        for ff= 1:size(toLoad,2)
            meansorted= [meansorted, ',   ',  [toLoad{1,sortedfitsPureNum1{region2}(ff)}(2:end)] ];
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
        
        
        while searchind< size(toLoad,2)
            
            if searchind+5 <= size(toLoad,2)
                fprintf(fileID,['%23s', inserts, '\n' ], '',toLoad{searchind:searchind+5});
                for ff= 1:size(BestBest1{region2})
                    fprintf(fileID,[ '%18s%5s', insert, '\n'],toLoad{ff},'' ,  BestBest1{region2}(searchind:searchind+5,ff));
                end
                fprintf(fileID, ['\n\n\n\n\n']);
                searchind = searchind+6;
            else
                fprintf(fileID,['%23s', inserts, '\n' ], '',toLoad{searchind:end});
                for ff= 1:size(BestBest1{region2})
                    fprintf(fileID,[ '%18s%5s', insert, '\n'],toLoad{ff},'' ,  BestBest1{region2}(searchind:end,ff));
                end
                fprintf(fileID, ['\n\n\n\n\n']);
                searchind = searchind+6;
            end
        end
        
        
        fprintf(fileID, ['\n\n\n\n\n\n\n\n\n\n\n\n']);
        
    end
    
    fclose(fileID)     
end
