function clusterwise_permutation_analysis(variableName, var2useName, inputRegion, lastnumber)


%% if error use this piece function to compute the cluster size in clusterwise permutation test

varType = {'A1', 'R'};
listRatingVarInv = {'Val', 'Sal', 'RT_Rat'};
listRatingVar = {'Val', 'Sal', 'Rank', 'flavourRat', 'RT_Rat'};
noSatVar = {'flavourRat', 'flavourA', 'flavourB'};


folder_to_save = pwd;
CMat = matfile([folder_to_save '/regressreg_' inputRegion ,'.mat']);
PMat = matfile([folder_to_save '/regresspval_' inputRegion ,'.mat']);
vT = 1;
varType{vT} = 'A1';

if exist([folder_to_save '/permdistribution_' inputRegion ,'.mat'])
    delete([folder_to_save '/permdistribution_' inputRegion ,'.mat']);
end

% find the binomial limit for significance
alpha =0.05;
siVar = lastnumber/2;
clear kk
for n=1:siVar
    kk(n)=1-binocdf(n-1,siVar,alpha);
end

pos05 = find(kk< 0.05);     pos05 = pos05(1)/siVar;
pos01 = find(kk< 0.01);     pos01 = pos01(1)/siVar;
pos001= find(kk< 0.001);    pos001= pos001(1)/siVar;



for loop = 1: length(variableName)
    if ismember(variableName{loop}, listRatingVar)
        vT = 2;
    else
        vT =1;
    end
    eval (['pvar = PMat.pval', variableName{loop}, '_',varType{vT},'; cvar = CMat.corr', variableName{loop} ,'_',varType{vT},';'] );
    eval(['[', variableName{loop}, '_clusterSumEffect05, ', variableName{loop}, '_clustersize05] = do_pop_stats_Tuning(pvar, cvar,  siVar, alpha, pos05);']);
    eval(['[', variableName{loop}, '_clusterSumEffect01, ', variableName{loop}, '_clustersize01] = do_pop_stats_Tuning(pvar, cvar,  siVar, alpha, pos01);']);
    eval(['[', variableName{loop}, '_clusterSumEffect001, ', variableName{loop}, '_clustersize001] = do_pop_stats_Tuning(pvar, cvar,  siVar, alpha, pos001);']);
    
    if ~ismember(variableName{loop}, noSatVar) % if Dynamic Variables
        eval (['dchvar = PMat.pval', variableName{loop}, '_DirectChange_',varType{vT},'; invvar = PMat.pval', variableName{loop} ,'inv_',varType{vT},'; cinvvar = CMat.corr', variableName{loop} ,'inv_',varType{vT},';']);
        eval(['[', variableName{loop}, '_clusterSumEffect_DirectChange05, ', variableName{loop}, '_clustersize_DirectChange05] = do_pop_stats_DirectChange (dchvar,  siVar, alpha, pos05);']);
        eval(['[', variableName{loop}, '_clusterSumEffect_DirectChange01, ', variableName{loop}, '_clustersize_DirectChange01] = do_pop_stats_DirectChange (dchvar,  siVar, alpha, pos01);']);
        eval(['[', variableName{loop}, '_clusterSumEffect_DirectChange001, ', variableName{loop}, '_clustersize_DirectChange001] = do_pop_stats_DirectChange (dchvar, siVar, alpha, pos001);']);
        if ismember(variableName{loop}, listRatingVar)
            if ismember(variableName{loop}, listRatingVarInv)
                eval(['pvar = PMat.pval', variableName{loop}, 'mean_',varType{vT},'; cvar = CMat.corr', variableName{loop} ,'mean_',varType{vT},';']) % for those variables trial-by-trial in rating trials change the variable sfor satiation
            end
        end
        eval(['[', variableName{loop}, '_signedRankclustersize, ', variableName{loop}, '_signedRankclusterSumEffect, ', variableName{loop}, '_signedRankSelclustersize, ', variableName{loop}, '_signedRankSelclusterSumEffect] = do_pop_stats_Satiation (pvar, cvar, invvar, cinvvar, siVar, alpha);']);
    end
    
end





h=0; g=0;
while h<1
    try
        load([folder_to_save '/permdistribution_' inputRegion ])
        h=1;
        
    catch
        save ([folder_to_save '/permdistribution_' inputRegion ],  '*_signedRank*', '*_cluster*', '-v7.3')
        if g>0
            warning(['error while saving the pvalues in region ', inputRegion])
        end
        g=g+1;
    end
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
for perm = 1: size(fraction_Abspref,3)
    
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

for perm2 = 1:permutations
    
    [a,c]=corr_mixed(modeltypelocal, mfrA11(indeces_permuted190(:,perm2),:), mfrA12(indeces_permuted190(:,perm2),:), var2use(:,1), var2use(:,2),   typeCorrlocal);
    assert(size(c,2)==size(mfrA11,2));
    [tempCorr(:, :, perm2)]= a;
    [tempPval(:, :, perm2)]= c;
end