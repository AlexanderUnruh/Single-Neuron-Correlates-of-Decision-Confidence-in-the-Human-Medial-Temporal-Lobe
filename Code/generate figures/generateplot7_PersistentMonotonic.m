function generateplot7_PersistentMonotonic

%% this function plots histogram Figure S5
%

dbstop if error

folder2save =  '/media/Projects/Alex/Reclustered analysis All/PaperPlots'; mkdir(folder2save)
autocorr_folder =  '/media/Projects/Alex/Reclustered analysis All/mixed_model1_Spearman_window10_kernel200ms_step10ms__InversesVisualNutritional_autocorr/';
fold =  '/media/Projects/Alex/Reclustered analysis All/mixed_model1_Spearman_window10_kernel200ms_step10ms__InversesVisualNutritional/';
cd(fold)
autocorr_folderreg =  '/media/Projects/Alex/Reclustered analysis All/mixed_model1_Spearman_window10_kernel200ms_step10ms__regression_AP_RT_subpop_autocorr/';
foldreg =  '/media/Projects/Alex/Reclustered analysis All/mixed_model1_Spearman_window10_kernel200ms_step10ms__regression_AP_RT_subpop/';


global sizeletra sizetrialtype
sizeletra = 30;
sizetrialtype = 20;
shiftletra = 1.5;
shiftletrax = 2.3;
shifttrialtype = 0.5 + 0.5;
shiftpvalx= +0.2961-0.23 -0.15;
shiftpvaly= 1.3017;
shiftlegx= 2.80;
shiftlegy=5.5891;
sizelegx= 3.5;

shiftAllx =  1.6711;
shiftAlly =  0.1839;
shiftSelx =   2.2851;
shiftSely =  shiftAlly;
shfitcbhx= 0.2438;
shiftanot7x = 1.9;
shiftanot8x = 1.9;

possible_colors = [[0 0.4470 0.7410];...
    [0.9290 0.5940 0.1250];...  % [0.3250 0.8000  0.2980  ];...
    [0.9290 0.3250 0.0980];...
    [0.4940 0.1840 0.5560]];



heightTicks = 0.6858; heightmain = 7.2009; width =  12.1920;

accumulator=0;


% xaxis
pos2x = 16.6116; pos1x = 2.2860; diffx1x2 = pos2x -pos1x;

shiftall = 0.2;
pos1x = 2.2860 + shiftall;
pos2x = pos1x + diffx1x2;


sub2yadd = 2.2860 - 0.4;
diff32 = 0.9144;
sub3yadd = sub2yadd + diff32;

totalheight = heightmain  + heightTicks + heightTicks  + diff32 + sub2yadd + 1 +.2;

pos2xor = pos2x;
pos1xor = pos1x;





%% PERSISTENCE PLOT
figpos = [0,0, 30, 12];

regtype =1;
sorted = 0; % if sorted == 1 then perssisten cells sorted by lenght of persistent activity. if sorted == 1 by ere cell ID sorted
for  figures = 1:2

    for  plottype = 1:2
        close all
        fig1 =  figure('color', 'w', 'units', 'centimeters', 'position', figpos, 'visible', 'off');
        if figures ==2
            plot_persistent(fig1, {'ALL'}, autocorr_folder, fold, {'RT'}, {'RT'}, plottype, regtype, sorted)
        elseif figures ==1
            plot_persistent(fig1, {'ALL'}, autocorr_folder, fold,  {'Abspref'}, {'Confidence'}, plottype, regtype, sorted)
        end
    end
end
sadasdasdasd


%% MONOTONY PLOT

figpos = [0, 0, 15, 15];

regtype = 1;
for figures = 1:2
    for plottype = 2:-1:1
        close all
        fig1 =  figure('color', 'w', 'units', 'centimeters', 'position', figpos, 'visible', 'off');
        if figures ==2
            plot_monotonic(fig1, {'ALL'}, autocorr_folder, fold, {'RT'}, {'RT'}, plottype, regtype)
        elseif figures ==1
            plot_monotonic(fig1, {'ALL'}, autocorr_folder, fold,  {'Abspref'}, {'Confidence'}, plottype, regtype)
        end
    end
end





function plot_persistent(fig1, region2, folder_auto, fold, var2use, var2name, plottype, regtype, Sorted)


set(0, 'currentFigure',fig1)


dbstop if error
%plot_paper_inverses_autocorr_reduced(0.05)



set(groot,'defaultAxesLinewidth', 1)
alpha = 0.05;


AreaRat = [501, 701];
AreaA1  = [501, 601];
AreaB   = [621, 721];
AreaA2  = [741, 841];


windowstep = 10;
trialonset = 501;
tickpossize= 1000/windowstep;
leftReduction = (3600+2350)/windowstep-1;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr
rightReduction =  trialonset-(2350)/windowstep;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr

alphaSigned = 0.15; % top value for the satiation effect plot
heightplot = 0.2;

ggg =  1;

listRatingCombine = { 'Val', 'Sal', 'RT_R'};
listRating = {'Val', 'Sal', 'Valmean', 'Salmean', 'Rank', 'RTmean', 'RT_R', 'Sign'};



for g=1%:size(region2,1)
    
    region = region2{g}; load(['../Rasts matfiles/patientVector_', region]);
    vectorPat = eval(['patientVector_', region]); % extract patient ID of each cell
    vectorPat = ismember([vectorPat; vectorPat], [2,3,4,5,6,1,7,8,9,10,11,12]);
    

    
    for cc= 1%: size(var2use,2)
        
        variable = var2use{cc};
        disp([variable, '  ',  region])
        
        if regtype == 1
            if ismember(variable, listRating)
                if strcmp(variable, 'RT_R')
                    toLoad= {['*RTmeaninv_R']; ['*RT_R'] ; ['*RTmean_R']};  % compare trial by trial RT with inverse of mean RT, like in Val, Sal or Rank in rating trials.
                else
                    toLoad= {['*', variable, 'inv_R']; ['*', variable,'_R']}; % all the rest
                    if ismember(variable, listRatingCombine)
                        toLoad= {['*', variable, 'inv_R']; ['*', variable,'_R']; ['*', variable,'mean_R']}; % all the rest
                    end
                end
            else
                toLoad= {['*', variable, 'inv_A1']; ['*', variable, '_A1']};
            end
            load([fold,'/regresspval_',region,'.mat'],toLoad{:});
            load([fold,'/regressreg_', region,'.mat'],toLoad{:});
            
            invvar =  eval(['pval', toLoad{1}(2:end)]);
            var=      eval(['pval', toLoad{2}(2:end)]);
            cinvvar=  eval(['corr', toLoad{1}(2:end)]);
            cvar=     eval(['corr', toLoad{2}(2:end)]);
            
        elseif regtype == 2
            
            toLoad= {['*', variable, 'reginv_A1']; ['*', variable, 'reg_A1']; ['*', variable, 'Autoreginv_A1']; ['*', variable, 'Autoreg_A1']};
            
            load([fold,'/regresspval_',region,'.mat'],toLoad{:});
            load([fold,'/regressreg_', region,'.mat'],toLoad{:});
            invvar =  eval(['pval', toLoad{1}(2:end)]);
            var=      eval(['pval', toLoad{2}(2:end)]);
            cinvvar=  eval(['corr', toLoad{1}(2:end)]);
            cvar=     eval(['corr', toLoad{2}(2:end)]);
                        
            invvara =  eval(['pval', toLoad{3}(2:end)]);
            vara    =  eval(['pval', toLoad{4}(2:end)]);
            cinvvara=  eval(['corr', toLoad{3}(2:end)]);
            cvara=     eval(['corr', toLoad{4}(2:end)]);
            
        end
   
        
        if ismember(variable, listRatingCombine)
            cvar =  eval(['corr', toLoad{3}(2:end)]);
        end
        
        
        %% define the limits of the plot
        x_ax_lim = [rightReduction, trialonset + leftReduction];
        
        
        
        invvar = invvar(vectorPat, :);
        var = var (vectorPat, :);
        cinvvar= cinvvar(vectorPat, :);
        cvar = cvar(vectorPat, :);
        
        
        siVar= size(var,1)/2;
        siVar2 = siVar;
        
        
        
        
        cvar1 = cvar(1:siVar, :);
        cvar2 = cvar(1+siVar:end, :);
        cinvvar1 = cinvvar(1:siVar, :);
        cinvvar2 = cinvvar(1+siVar:end, :);
        
        
        both = ((var(1:siVar,:)<sqrt(2*alpha)) & (var(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvar(1:siVar,:))==sign(cvar(1+siVar:end,:))));
        bothinv = ((invvar(1:siVar,:)<sqrt(2*alpha)) & (invvar(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cinvvar(1:siVar,:))==sign(cinvvar(1+siVar:end,:))));
        
        
        selectVec =   [both] | [ bothinv];
        
        var2plot1 = nan(size(cvar1));
        var2plot2 = nan(size(cvar1));
        invvar2plot1 = nan(size(cvar1));
        invvar2plot2 = nan(size(cvar1));
        
        
        for ll = 1: size(var,2)
            
            var2plot1 (selectVec(:,ll),ll) = cvar1(selectVec(:,ll),ll);
            var2plot2 (selectVec(:,ll),ll) = cvar2(selectVec(:,ll),ll);
            invvar2plot1 (selectVec(:,ll),ll) = cinvvar1(selectVec(:,ll),ll);
            invvar2plot2 (selectVec(:,ll),ll) = cinvvar2(selectVec(:,ll),ll);
            
        end
        
        cvarSel = [var2plot1; var2plot2];
        cinvvarSel = [invvar2plot1; invvar2plot2];
        
        siVar3 = (sum(selectVec));
        
        
        
        smeancvar = (cvar(1:siVar,:) + cvar(1+siVar:end ,:))/2;
        smeancinvvar = (cinvvar(1:siVar,:) + cinvvar(1+siVar:end ,:))/2;
        
        meancvar = abs(cvar(1:siVar,:) + cvar(1+siVar:end ,:))/2;
        meancinvvar = abs(cinvvar(1:siVar,:) + cinvvar(1+siVar:end ,:))/2;
        
        meancvarSel = abs(cvarSel(1:siVar,:) + cvarSel(1+siVar:end ,:))/2;
        meancinvvarSel = abs(cinvvarSel(1:siVar,:) + cinvvarSel(1+siVar:end ,:))/2;
        
        
        toplot2= ((var(1:siVar,:)<sqrt(2*alpha)) & (var(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvar(1:siVar,:))==sign(cvar(1+siVar:end,:))));
        
        

        
        
        
        %% ticks
        
        Xtickspos=[[trialonset:-tickpossize: 1], [trialonset:tickpossize:size(cvar,2)]];  Xtickspos = sort(Xtickspos(2:end)); % from the trial aligment, go forward and backwards determining te position and sort
        Xtickslabelsnames = [[0:-1:-(length([trialonset:-tickpossize: 1])-1)], [0:1:length(trialonset:tickpossize:size(cvar,2))-1]]; Xtickslabelsnames = sort(Xtickslabelsnames(2:end));
        
        for hh=1:length(Xtickslabelsnames)
            aa{hh}=num2str(Xtickslabelsnames(hh));
        end
        Xtickslabelsnames=[];Xtickslabelsnames=aa;
        
        
        
        clim = [-.25, .25];
        if Sorted == 1
            [~, ind] = sort(nansum(toplot2(:,trialonset:x_ax_lim(2)),2), 'ascend');
        else
            ind = 1: size(toplot2,1);
        end
        subplot(1,2,1)
        hold on
        if plottype == 1
            colormap(flipud(gray))
            imagesc(toplot2(ind,:))
        elseif plottype == 2
            colormap(flipud(parula))
            imagesc(meancvar, clim)
            colorbar
        end
        ylabel('Unit ID')
        
                ytop = size(toplot2,1);
        if strcmp(toLoad{1}(end-1:end), '_R')
            area(AreaRat, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            
            xlabel('Time from stimulus onset (s)')
            
            xlim(x_ax_lim)
        else
            
            %             area(AreaA1, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            %             area(AreaB, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.3,0.3,0.3], 'Edgecolor', 'none')
            %             area(AreaA2, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            line([AreaA1(1),AreaA1(1)], [0,ytop], 'linestyle', '-','color', [0.1   0.1    0.1], 'linewidth', 1.3)
            line([AreaB(1),AreaB(1)], [0,ytop], 'linestyle', '-','color', [0.1   0.1    0.1], 'linewidth', 1.3)
            line([AreaA2(1),AreaA2(1)], [0,ytop], 'linestyle', '-','color',  [0.1   0.1    0.1], 'linewidth', 1.3)
            xlabel('Time relative to stimulus A1 onset (s)')
            
            xlim(x_ax_lim)
        end

        ylim([0, ytop])
        
        set(gca, 'Xtick', Xtickspos, 'Xticklabels',Xtickslabelsnames);
        set(gca, 'FontSize',12, 'fontweight', 'bold')
        
        
        
        
        
        %% autocorr plot
        
        if regtype == 1
            clear corr* pval* varchange toplot2
            
            
            
            load([folder_auto, '/regresspval_',region,'.mat'],toLoad{:});
            load([folder_auto, '/regressreg_',region,'.mat'],toLoad{:});
            
            
            invvar =  eval(['pval', toLoad{1}(2:end)]);
            var =      eval(['pval', toLoad{2}(2:end)]);
            cinvvar =  eval(['corr', toLoad{1}(2:end)]);
            cvar =     eval(['corr', toLoad{2}(2:end)]);
            if ismember(variable, listRatingCombine)
                cvar =  eval(['corr', toLoad{3}(2:end)]);
            end
        else
            
            invvar = invvara;
            var =     vara;
            cinvvar =   cinvvara;
            cvar = cvara;
        end
        
        
        
        invvar = invvar(vectorPat, :);
        var = var (vectorPat, :);
        cinvvar= cinvvar(vectorPat, :);
        cvar = cvar(vectorPat, :);
        
        
        siVar= size(var,1)/2;
        siVar2 = siVar;
        
        
        
        
        cvar1 = cvar(1:siVar, :);
        cvar2 = cvar(1+siVar:end, :);
        cinvvar1 = cinvvar(1:siVar, :);
        cinvvar2 = cinvvar(1+siVar:end, :);
        
        
        both = ((var(1:siVar,:)<sqrt(2*alpha)) & (var(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvar(1:siVar,:))==sign(cvar(1+siVar:end,:))));
        bothinv = ((invvar(1:siVar,:)<sqrt(2*alpha)) & (invvar(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cinvvar(1:siVar,:))==sign(cinvvar(1+siVar:end,:))));
        
        
        selectVec =   [both] | [ bothinv];
        
        var2plot1 = nan(size(cvar1));
        var2plot2 = nan(size(cvar1));
        invvar2plot1 = nan(size(cvar1));
        invvar2plot2 = nan(size(cvar1));
        
        for ll = 1: size(var,2)
            
            var2plot1 (selectVec(:,ll),ll) = cvar1(selectVec(:,ll),ll);
            var2plot2 (selectVec(:,ll),ll) = cvar2(selectVec(:,ll),ll);
            invvar2plot1 (selectVec(:,ll),ll) = cinvvar1(selectVec(:,ll),ll);
            invvar2plot2 (selectVec(:,ll),ll) = cinvvar2(selectVec(:,ll),ll);
        end
        
        cvarSel = [var2plot1; var2plot2];
        cinvvarSel = [invvar2plot1; invvar2plot2];
        
        siVar3 = (sum(selectVec));
        
        
        
        
        smeancvar = (cvar(1:siVar,:) + cvar(1+siVar:end ,:))/2;
        smeancinvvar = (cinvvar(1:siVar,:) + cinvvar(1+siVar:end ,:))/2;
        
        meancvar = abs(cvar(1:siVar,:) + cvar(1+siVar:end ,:))/2;
        meancinvvar = abs(cinvvar(1:siVar,:) + cinvvar(1+siVar:end ,:))/2;
        
        meancvarSel = abs(cvarSel(1:siVar,:) + cvarSel(1+siVar:end ,:))/2;
        meancinvvarSel = abs(cinvvarSel(1:siVar,:) + cinvvarSel(1+siVar:end ,:))/2;
        
        
        toplot2= ((var(1:siVar,:)<sqrt(2*alpha)) & (var(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvar(1:siVar,:))==sign(cvar(1+siVar:end,:))));
        
        
 
      
        subplot(1,2,2)
        hold on
        
        if plottype == 1
            colormap(flipud(gray))
            imagesc(toplot2(ind,:))
        elseif plottype == 2
            colormap(flipud(parula))
            imagesc(meancvar, clim)
            colorbar
        end
        
                ytop = size(toplot2,1);
        if strcmp(toLoad{1}(end-1:end), '_R')
            area(AreaRat, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            
            xlabel('Time from stimulus onset (s)')
            
            xlim(x_ax_lim)
        else
            
            %             area(AreaA1, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            %             area(AreaB, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.3,0.3,0.3], 'Edgecolor', 'none')
            %             area(AreaA2, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            line([AreaA1(1),AreaA1(1)], [0,ytop], 'linestyle', '-','color', [0.1   0.1    0.1], 'linewidth', 1.3)
            line([AreaB(1),AreaB(1)], [0,ytop], 'linestyle', '-','color', [0.1   0.1    0.1], 'linewidth', 1.3)
            line([AreaA2(1),AreaA2(1)], [0,ytop], 'linestyle', '-','color',  [0.1   0.1    0.1], 'linewidth', 1.3)
            xlabel('Time relative to stimulus A1 onset (s)')
            
            xlim(x_ax_lim)
        end
        
        ylim([0, ytop])
        
        set(gca, 'Xtick',Xtickspos, 'Xticklabels',Xtickslabelsnames);
        set(gca, 'FontSize',12, 'fontweight', 'bold')
        
        
        
        
        
  
        global sizeletra sizetrialtype
        
        anot7 = annotation('textbox',[0.215  0.97   0.3    0.04],   'string',['Current Trial'],'Units','normalized','FontName','Helvetica','FontSize',sizetrialtype,'EdgeColor','none', 'fontweight', 'bold');
        anot8 = annotation('textbox',[0.665  0.97   0.3    0.04],   'string',['Next Trial'],'Units','normalized','FontName','Helvetica','FontSize',sizetrialtype,'EdgeColor','none', 'fontweight', 'bold');
        
        
        
        folder2save =  '/media/Projects/Alex/Reclustered analysis All/PaperPlots'; mkdir(folder2save)
        if Sorted == 1 
        toSave = [folder2save, '/PersistentplotSorted_', var2name{:}, num2str(plottype)];
        else 
            toSave = [folder2save, '/Persistentplot_', var2name{:}, num2str(plottype)];
        end
        
        print([toSave, '200'], '-dpng', '-r200');
        print([toSave, '300'], '-dpng', '-r300');
        
        print([toSave, 'Vec'], '-depsc', '-painters');
        print([toSave, 'Vec'], '-depsc', '-opengl', '-r500');
        set(fig1,'Units','centimeters');
        pos = get(fig1,'Position');
        set(fig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
        
        print([toSave, 'Vec'], '-dpdf', '-painters');
        print(toSave, '-dpdf', '-r500', '-opengl');
        print([toSave, ' 300'], '-dpdf', '-r300', '-opengl');
        
        
    end
end








function plot_monotonic(fig1, region2, folder_auto, fold, var2use, var2name, plottype, regtype)


set(0, 'currentFigure',fig1)
dbstop if error

set(groot,'defaultAxesLinewidth', 1)
alpha = 0.05;

winOfIntereston = 501+120*3;
winOfInterestoff = winOfIntereston+235;


listRatingCombine = {'Val', 'Sal', 'RT_R'};
listRating = {'ValTot', 'SalTot', 'Val', 'Sal', 'Valmean', 'Salmean', 'Rank', 'RTmean', 'RTstd', 'Valstd', 'RT_R', 'Sign', 'SignTot',   'Negsplit', 'Possplit', 'NegsplitTot', 'PossplitTot',  'ValSweetsplit', 'SalSweetsplit', 'ValSweetsplitTot', 'SalSweetsplitTot', 'RankSaltysplit', ...
    'ValSweetsplit', 'SalSweetsplit', 'ValSweetsplitTot', 'SalSweetsplitTot', 'RankSweetsplit'};



for g = 1
    
    region = region2{g}; load(['../Rasts matfiles/patientVector_', region]);
    vectorPat = eval(['patientVector_', region]); % extract patient ID of each cell
    vectorPat = ismember([vectorPat; vectorPat], [1, 2,3,4,5,6,7,8,9,10,11,12]);
    

    
    for cc= 1%: size(var2use,2)
        
        variable = var2use{cc};
        disp([variable, '  ',  region])
        
        if regtype == 1
            if ismember(variable, listRating)
                if strcmp(variable, 'RT_R')
                    toLoad= {['*RTmeaninv_R']; ['*RT_R'] ; ['*RTmean_R']};  % compare trial by trial RT with inverse of mean RT, like in Val, Sal or Rank in rating trials.
                else
                    toLoad= {['*', variable, 'inv_R']; ['*', variable,'_R']}; % all the rest
                    if ismember(variable, listRatingCombine)
                        toLoad= {['*', variable, 'inv_R']; ['*', variable,'_R']; ['*', variable,'mean_R']}; % all the rest
                    end
                end
            else
                toLoad= {['*', variable, 'inv_A1']; ['*', variable, '_A1']};
            end
            load([fold,'/regresspval_',region,'.mat'],toLoad{:});
            load([fold,'/regressreg_', region,'.mat'],toLoad{:});
            
            var=      eval(['pval', toLoad{2}(2:end)]);
            cvar=     eval(['corr', toLoad{2}(2:end)]);
            
        elseif regtype == 2
            toLoad= {['*', variable, 'reginv_A1']; ['*', variable, 'reg_A1']; ['*', variable, 'Autoreginv_A1']; ['*', variable, 'Autoreg_A1']};
            
            load([fold,'/regresspval_',region,'.mat'],toLoad{:});
            load([fold,'/regressreg_', region,'.mat'],toLoad{:});
            var=      eval(['pval', toLoad{2}(2:end)]);
            cvar=     eval(['corr', toLoad{2}(2:end)]);
            
            vara=      eval(['pval', toLoad{4}(2:end)]);
            cvara=     eval(['corr', toLoad{4}(2:end)]);   
        end
        
        
        
        if ismember(variable, listRatingCombine)
            cvar =  eval(['corr', toLoad{3}(2:end)]);
        end
        
        
        %% define the limits of the plot

        var = var (vectorPat, :);
        cvar = cvar(vectorPat, :);
        
        
        siVar= size(var,1)/2;
        siVar2 = siVar;
        
        
        
        
        cvar1 = cvar(1:siVar, :);
        cvar2 = cvar(1+siVar:end, :);
        
        both = ((var(1:siVar,:)<sqrt(2*alpha)) & (var(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvar(1:siVar,:))==sign(cvar(1+siVar:end,:))));
 
        selectVec =   [both]; % only if correlated for the normal one
        
        var2plot1 = nan(size(cvar1));
        var2plot2 = nan(size(cvar1));
     
        for ll = 1: size(var,2)
            
            var2plot1 (selectVec(:,ll),ll) = cvar1(selectVec(:,ll),ll);
            var2plot2 (selectVec(:,ll),ll) = cvar2(selectVec(:,ll),ll);

        end
        
        
        cvarSel = [var2plot1; var2plot2];
        siVar3 = (sum(selectVec));       

        smeancvar = (cvar(1:siVar,:) + cvar(1+siVar:end ,:))/2;
        meancvar = abs(cvar(1:siVar,:) + cvar(1+siVar:end ,:))/2;  
        smeancvarSel = (cvarSel(1:siVar,:) + cvarSel(1+siVar:end ,:))/2;
        meancvarSel = abs(cvarSel(1:siVar,:) + cvarSel(1+siVar:end ,:))/2;


        if plottype==1
            numcells = nanmean(smeancvar(:,winOfIntereston:winOfInterestoff),2);
            toadd = '_Allbins';
            titlename = [var2name{1}];
        elseif plottype ==2
            numcells = nanmean(smeancvarSel(:,winOfIntereston:winOfInterestoff),2);
            toadd = '_Selbins';
            titlename = [var2name{1}];
        end
        
       
        histogram(numcells, 'BinWidth',0.01, 'Normalization','probability')
        xlim([-0.5,0.5]); 
        set(gca, 'XTick', [-.5, -.25, 0, 0.25, .5], 'YTick', [0:0.02:0.06])
        ylabel('Fraction of Units')
        xlabel('Mean Spearman Rho')
        
        set(gca, 'FontSize',12, 'fontweight', 'bold')
        
        
        folder2save =  '/media/Projects/Alex/Reclustered analysis All/PaperPlots'; mkdir(folder2save)
       
        toSave = [folder2save, '/Histogramplot_', var2name{:}, toadd];
        
        
        title (titlename) 
        
        print([toSave, '200'], '-dpng', '-r200');
        print([toSave, '300'], '-dpng', '-r300');
        
        print([toSave, 'Vec'], '-depsc', '-painters');
        print([toSave, 'Vec'], '-depsc', '-opengl', '-r500');
        set(fig1,'Units','centimeters');
        pos = get(fig1,'Position');
        set(fig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
        
        print([toSave, 'Vec'], '-dpdf', '-painters');
        print(toSave, '-dpdf', '-r500', '-opengl');
        print([toSave, ' 300'], '-dpdf', '-r300', '-opengl');
        
        
    end
end


