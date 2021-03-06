function generateplot3_clusterwisePermutation


%% This script generates Figures of the Transient neural activity, value, salience ... with clusterwise correction
dbstop if error

folder2save =  '/media/Projects/Alex/Reclustered analysis All/PaperPlots'; mkdir(folder2save)
% changefolder = '/media/raid/Alex/Reclustered analysis All/mixed_model1_Spearman_window10_kernel200ms_step10ms__changeINvalue/';
fold =  '/media/Projects/Alex/Reclustered analysis All/mixed_model1_Spearman_window10_kernel200ms_step10ms__InversesVisualNutritional/';
cd(fold)

StimIDfold = '/media/raid/Alex/Reclustered analysis All/mixed_model1_Spearman_window10_kernel200ms_step10ms_StimulusIdentity_Anova';

clustPermFold = '/media/raid/Alex/Reclustered analysis All/clusterwise_permutation_1_Spearman_window10_kernel200ms_step10ms_Val_Sal_Rank_2AFC_correctedWindow/';
DirectChangeFold = '/media/Projects/Alex/Reclustered analysis All/mixed_model1_Spearman_window10_kernel200ms_step10ms__changeINvalue/';
clustPermFold2 = '/media/raid/Alex/Reclustered analysis All/clusterwise_permutation_1_Spearman_window10_kernel200ms_step10ms_TasteValSalRT_correctedWindow';
DirectChangeFold2 = '/media/raid/Alex/Reclustered analysis All/mixed_model1_Spearman_window10_kernel200ms_step10ms_ChangeInVariables/';


global sizeletra sizetrialtype type_clusterSumEffect typeclusterSumEffect_DirectChange
global alpha alphaBinom alphaPerm alphaSigned

% Alpha values from tests
alpha = 0.05;       % Alpha from correlation over trials
alphaSigned = 0.05; % Alpha from signed rank test
alphaBinom = 0.001;  % Alpha from binomial test across neurons
alphaPerm = 0.01;   % Alpha from clusterwise permutation across time


stralphaBinom= num2str(alphaBinom);
type_clusterSumEffect = ['clusterSumEffect', stralphaBinom(3:end)];
typeclusterSumEffect_DirectChange  =  ['clusterSumEffect_DirectChange', stralphaBinom(3:end)];



sizeletra = 30;
sizetrialtype = 20;


shiftletray = 1.75;
shiftletrax = 2.3;

shifttrialtype = 0.6  +0.5;
sizelegx= 3.5 -.2;
shiftanot6x = 3.8;
shiftanot8x = 2.1;
shiftanot7x = shiftanot8x;



shiftAllx =  1.6711;
shiftAlly =  0.1839;
shiftSelx =   2.2851;
shiftSely =  shiftAlly;




cbhpos = [33.8531     0.2286    0.2845    1.6002];
pospval = [33.6220      2.2860 0.7112 0.2286];
legpos = [   30.7758    9.0163    4.5508    1.9447];



shiftcbhx = 33.8531- (23.6474 +  9.9568);
shiftpvalx =   33.6220- (23.6474 +  9.9568) -.2;
shiftlegx =   30.7758 - (23.6474 +  9.9568) + 0.6;


shiftpvaly =  2.2860 -   1.1430 +0.2;
shiftlegy =  9.0163 -  3.4290 - 0.2 ;



possible_colors = [[0 0.4470 0.7410];...
    [0.9290 0.5940 0.1250];...  % [0.3250 0.8000  0.2980  ];...
    [0.9290 0.3250 0.0980];...
    [0.4940 0.1840 0.5560]];



possible_colors =  [[0    0.4470    0.7410];...
    [0.9290    0.5940    0.1250];...
    [0.9290    0.3250    0.0980];...
    [0.4940    0.1840    0.5560];...
    [0.4660    0.6740    0.1880];...
    [0.3010    0.7450    0.9330];...
    [0.6350    0.0780    0.1840]];






heightTicks = 0.6858; heightmain = 7.2009; width2= 7.7; width1 = (width2/1162) * 922;

% xaxis
pos3x= 23.6474; pos2x =  11.9126; pos1x =  2.6670; diffx1x2 = pos2x -pos1x -1; diffx2x3 = pos3x -pos2x-1;

shiftall = 0.2;
%pos1x =  2.6670 +shiftall; %shift a bit
pos1x =  2.06 +shiftall; %shift a bit


interdist= 2.5;
diffx1x2 =  width1 + interdist;
diffx2x3 =  width2 + interdist;
pos2x = pos1x + diffx1x2;
pos3x = pos2x + diffx2x3;
sub2yadd = 2.2860 - 0.4;
diff32 = 0.9144;
sub3yadd = sub2yadd + diff32;

pos1xor=    pos1x;
pos2xor =   pos2x;
totalheight = heightmain  + heightTicks + heightTicks  + diff32 + sub2yadd + 1 +.2;








%%inlludeing ranking

% plot only the change in value


figpos = [0,0, 31.4800+shiftall, 35.8];
close all
fig1 =  figure('color', 'w', 'units', 'centimeters', 'position', figpos, 'visible', 'off');
accumulator = 0;

for  figures = 1:3
    % mainplots  and satiation
    mainy = 1.4290 + accumulator;
    sub2y=  mainy - sub2yadd;
    sub3y=  mainy - sub3yadd;
    
    
    pos31    =  [pos3x    mainy  width2    heightmain];
    pos32    =  [pos3x    sub2y  width2    heightTicks];
    pos33    =  [pos3x    sub3y  width2    heightTicks];
    pos21    =  [pos2x    mainy  width2    heightmain];
    pos22    =  [pos2x    sub2y  width2    heightTicks];
    pos23    =  [pos2x    sub3y  width2    heightTicks];
    pos11    =  [pos1x    mainy  width1    heightmain];
    pos12    =  [pos1x    sub2y  width1    heightTicks];
    pos13    =  [pos1x    sub3y  width1    heightTicks];
    
    
    
    
    
    % letter for figure
    heightanotA = mainy + heightmain + shiftletray;
    anot5pos = [pos1x-shiftletrax, heightanotA, 4.5720 0.4286];
    
    % current next trial, size 12
    heightanotC = mainy + heightmain + shifttrialtype;
    
    anot6pos = [pos1x+width1/2-shiftanot6x, heightanotC,  8.5720 0.4286];
    anot7pos = [pos2x+width2/2-shiftanot7x heightanotC,  78.5720 0.4286];
    anot8pos = [pos3x+width2/2-shiftanot8x, heightanotC, 78.5720 0.4286];
    
    
    % anotations satitation
    anot2x = pos1x - shiftAllx - shiftall; anot2y = sub2y + shiftAlly;
    anot3x = pos1x - shiftSelx - shiftall; anot3y = sub3y + shiftSely;
    cbhx   =  pos3x+width2 + shiftcbhx;  cbhy   =  sub3y ;
    anot4x =  pos3x+width2 +shiftpvalx;  anot4y =  sub2y + shiftpvaly;
    
    cbhpos = [cbhx, cbhy,  0.2438    1.6002];
    anot2pos = [anot2x, anot2y,  9.2234    0.4889]  ; %'All'
    anot3pos = [anot3x, anot3y,  9.2234    0.4889]  ; %'Seected';
    anot4pos = [anot4x, anot4y,  0.2, 0.2]  ; %'pval'
    
    % legend
    legx =  pos3x+width2 + shiftlegx ; legy = mainy + shiftlegy;
    legpos =[legx, legy, sizelegx   1.4817];
    

    if figures == 1
% do_the_plot(fig1, {'ALL'}, fold, clustPermFold, clustPermFold2, DirectChangeFold,  DirectChangeFold2, {'val'}, {'Value'}, pos11, pos12,pos13, pos21, pos22, pos23, pos31, pos32, pos33, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot6pos , anot7pos, anot8pos , legpos, possible_colors([1,2], :), 'A', 0, 1 ,1, 0, 0)
        do_the_plot(fig1, {'ALL'}, fold, clustPermFold, clustPermFold2, DirectChangeFold,  DirectChangeFold2, {'sal'}, {'Unsigned Rating'}, pos11, pos12,pos13, pos21, pos22, pos23, pos31, pos32, pos33, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot6pos , anot7pos, anot8pos , legpos, possible_colors([4], :), 'C', 1, 1, 0, 0, 0)
    elseif figures == 2
        do_the_plot(fig1, {'ALL'}, fold, clustPermFold, clustPermFold2, DirectChangeFold,  DirectChangeFold2, {'rank'}, {'Ranking'}, pos11, pos12,pos13, pos21, pos22, pos23, pos31, pos32, pos33, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot6pos , anot7pos, anot8pos , legpos, possible_colors([2,2], :), 'B', 0, 1, 0, 0, 0)
    elseif figures == 3
        do_the_plot(fig1, {'ALL'}, fold, clustPermFold, clustPermFold2, DirectChangeFold,  DirectChangeFold2, {'val'}, {'Rating'}, pos11, pos12,pos13, pos21, pos22, pos23, pos31, pos32, pos33, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot6pos , anot7pos, anot8pos , legpos, possible_colors([1,2], :), 'A', 0, 1 ,1, 0, 0)
    end
    accumulator = accumulator + totalheight;
end

toSave = [folder2save, '/Clusterwise Permutation test, value_satiation ranking population'];
print([toSave, '200'], '-dpng', '-r200');


% % %     print(toSave, '-dpng', '-r500');
print([toSave, '300'], '-dpng', '-r300');



% plot only the change in value


figpos = [0,0, 31.4800+shiftall, 35.8];
close all
fig1 =  figure('color', 'w', 'units', 'centimeters', 'position', figpos, 'visible', 'off');
accumulator = 0;

for  figures = 1:3
    % mainplots  and satiation
    mainy = 1.4290 + accumulator;
    sub2y=  mainy - sub2yadd;
    sub3y=  mainy - sub3yadd;
    
    
    pos31    =  [pos3x    mainy  width2    heightmain];
    pos32    =  [pos3x    sub2y  width2    heightTicks];
    pos33    =  [pos3x    sub3y  width2    heightTicks];
    pos21    =  [pos2x    mainy  width2    heightmain];
    pos22    =  [pos2x    sub2y  width2    heightTicks];
    pos23    =  [pos2x    sub3y  width2    heightTicks];
    pos11    =  [pos1x    mainy  width1    heightmain];
    pos12    =  [pos1x    sub2y  width1    heightTicks];
    pos13    =  [pos1x    sub3y  width1    heightTicks];
    
    
 
    % letter for figure
    heightanotA = mainy + heightmain + shiftletray;
    anot5pos = [pos1x-shiftletrax, heightanotA, 4.5720 0.4286];
    
    % current next trial, size 12
    heightanotC = mainy + heightmain + shifttrialtype;
    
    anot6pos = [pos1x+width1/2-shiftanot6x, heightanotC,  8.5720 0.4286];
    anot7pos = [pos2x+width2/2-shiftanot7x heightanotC,  78.5720 0.4286];
    anot8pos = [pos3x+width2/2-shiftanot8x, heightanotC, 78.5720 0.4286];
    
    
    % anotations satitation
    anot2x = pos1x - shiftAllx - shiftall; anot2y = sub2y + shiftAlly;
    anot3x = pos1x - shiftSelx - shiftall; anot3y = sub3y + shiftSely;
    cbhx   =  pos3x+width2 + shiftcbhx;  cbhy   =  sub3y ;
    anot4x =  pos3x+width2 +shiftpvalx;  anot4y =  sub2y + shiftpvaly;
    
    cbhpos = [cbhx, cbhy,  0.2438    1.6002];
    anot2pos = [anot2x, anot2y,  9.2234    0.4889]  ; %'All'
    anot3pos = [anot3x, anot3y,  9.2234    0.4889]  ; %'Seected';
    anot4pos = [anot4x, anot4y,  0.2, 0.2]  ; %'pval'
    
    % legend
    legx =  pos3x+width2 + shiftlegx ; legy = mainy + shiftlegy;
    legpos =[legx, legy, sizelegx   1.4817];
    

    if figures == 1
% do_the_plot(fig1, {'ALL'}, fold, clustPermFold, clustPermFold2, DirectChangeFold,  DirectChangeFold2, {'val'}, {'Value'}, pos11, pos12,pos13, pos21, pos22, pos23, pos31, pos32, pos33, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot6pos , anot7pos, anot8pos , legpos, possible_colors([1,2], :), 'A', 0, 1 ,1, 0, 0)
        do_the_plot(fig1, {'ALL'}, fold, clustPermFold, clustPermFold2, DirectChangeFold,  DirectChangeFold2, {'flavour'}, {'Taste'}, pos11, pos12,pos13, pos21, pos22, pos23, pos31, pos32, pos33, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot6pos , anot7pos, anot8pos , legpos, possible_colors([4], :), 'C', 1, 1, 0, 0, 0)
    elseif figures == 2
        do_the_plot(fig1, {'ALL'}, fold, clustPermFold, clustPermFold2, DirectChangeFold,  DirectChangeFold2, {'sal'}, {'Unsigned Value'}, pos11, pos12,pos13, pos21, pos22, pos23, pos31, pos32, pos33, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot6pos , anot7pos, anot8pos , legpos, possible_colors([2,2], :), 'B', 0, 1, 0, 0, 0)
    elseif figures == 3
        do_the_plot(fig1, {'ALL'}, fold, clustPermFold, clustPermFold2, DirectChangeFold,  DirectChangeFold2, {'val'}, {'Value'}, pos11, pos12,pos13, pos21, pos22, pos23, pos31, pos32, pos33, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot6pos , anot7pos, anot8pos , legpos, possible_colors([1,2], :), 'A', 0, 1 ,1, 0, 0)
    end
    accumulator = accumulator + totalheight;
end

toSave = [folder2save, '/Clusterwise Permutation test, value_satiation population'];
print([toSave, '200'], '-dpng', '-r200');





    print(toSave, '-dpng', '-r500');
print([toSave, '300'], '-dpng', '-r300');

    print([toSave, 'Vec'], '-depsc', '-painters');
    print([toSave, 'Vec'], '-depsc', '-opengl', '-r500');
    set(fig1,'Units','centimeters');
    pos = get(fig1,'Position');
    set(fig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

    print([toSave, 'Vec'], '-dpdf', '-painters');
    print(toSave, '-dpdf', '-r500', '-opengl');
    print([toSave, ' 300'], '-dpdf', '-r300', '-opengl');







%% only direct change


        
        % Alpha values from tests
        alpha = 0.05;       % Alpha from correlation over trials
        alphaSigned = 0.05; % Alpha from signed rank test
        alphaBinom = 0.05;  % Alpha from binomial test across neurons
        alphaPerm = 0.01;   % Alpha from clusterwise permutation across time
        
        stralphaBinom= num2str(alphaBinom);
        type_clusterSumEffect = ['clusterSumEffect', stralphaBinom(3:end)];
        typeclusterSumEffect_DirectChange  =  ['clusterSumEffect_DirectChange', stralphaBinom(3:end)];
        % difference in value only
        
        
        totalheight2 = heightmain + sub2yadd  + 1 +.2;
        pos1x = pos1xor;
        pos2x = pos1x + diffx1x2;
        pos3x = pos2x + diffx2x3;
        accumulator = 0 ;
        figpos = [0,0, 31.4800+shiftall, 31.5];
        close all
        fig1 =  figure('color', 'w', 'units', 'centimeters', 'position', figpos, 'visible', 'off');
        accumulator = 0;
        
        for  figures = 1:3
            % mainplots  and satiation
            mainy = 1.4290 + accumulator;
            sub2y=  mainy - sub2yadd;
            sub3y=  mainy - sub3yadd;
            
            
            pos31    =  [pos3x    mainy  width2    heightmain];
            pos32    =  [pos3x    sub2y  width2    heightTicks];
            pos33    =  [pos3x    sub3y  width2    heightTicks];
            pos21    =  [pos2x    mainy  width2    heightmain];
            pos22    =  [pos2x    sub2y  width2    heightTicks];
            pos23    =  [pos2x    sub3y  width2    heightTicks];
            pos11    =  [pos1x    mainy  width1    heightmain];
            pos12    =  [pos1x    sub2y  width1    heightTicks];
            pos13    =  [pos1x    sub3y  width1    heightTicks];
            
            
            
            
            
            % letter for figure
            heightanotA = mainy + heightmain + shiftletray;
            anot5pos = [pos1x-shiftletrax, heightanotA, 4.5720 0.4286];
            
            % current next trial, size 12
            heightanotC = mainy + heightmain + shifttrialtype;
            
            anot6pos = [pos1x+width1/2-shiftanot6x, heightanotC,  8.5720 0.4286];
            anot7pos = [pos2x+width2/2-shiftanot7x heightanotC,  78.5720 0.4286];
            anot8pos = [pos3x+width2/2-shiftanot8x, heightanotC, 78.5720 0.4286];
            
            
            % anotations satitation
            anot2x = pos1x - shiftAllx - shiftall; anot2y = sub2y + shiftAlly;
            anot3x = pos1x - shiftSelx - shiftall; anot3y = sub3y + shiftSely;
            cbhx   =  pos3x+width2 + shiftcbhx;  cbhy   =  sub3y ;
            anot4x =  pos3x+width2 +shiftpvalx;  anot4y =  sub2y + shiftpvaly;
            
            cbhpos = [cbhx, cbhy,  0.2438    1.6002];
            anot2pos = [anot2x, anot2y,  9.2234    0.4889]  ; %'All'
            anot3pos = [anot3x, anot3y,  9.2234    0.4889]  ; %'Seected';
            anot4pos = [anot4x, anot4y,  0.2, 0.2]  ; %'pval'
            
            % legend
            legx =  pos3x+width2 + shiftlegx ; legy = mainy + shiftlegy;
            legpos =[legx, legy, sizelegx   1.4817];
            
            
            if figures == 1
                % do_the_plot(fig1, {'ALL'}, fold, clustPermFold, clustPermFold2, DirectChangeFold,  DirectChangeFold2, {'val'}, {'Value'}, pos11, pos12,pos13, pos21, pos22, pos23, pos31, pos32, pos33, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot6pos , anot7pos, anot8pos , legpos, possible_colors([1,2], :), 'A', 0, 1 ,1, 0, 0)
                do_the_plot(fig1, {'ALL'}, fold, clustPermFold, clustPermFold2, DirectChangeFold,  DirectChangeFold2, {'sal'}, {'Unsigned Rating'}, pos11, pos12,pos13, pos21, pos22, pos23, pos31, pos32, pos33, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot6pos , anot7pos, anot8pos , legpos, possible_colors([4,4], :), 'C', 1, 1, 0, 0, 1)
            elseif figures == 2
                do_the_plot(fig1, {'ALL'}, fold, clustPermFold, clustPermFold2, DirectChangeFold,  DirectChangeFold2, {'rank'}, {'Ranking'}, pos11, pos12,pos13, pos21, pos22, pos23, pos31, pos32, pos33, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot6pos , anot7pos, anot8pos , legpos, possible_colors([2,2], :), 'B', 0, 1, 0, 0, 1)
            elseif figures == 3
                do_the_plot(fig1, {'ALL'}, fold, clustPermFold, clustPermFold2, DirectChangeFold,  DirectChangeFold2, {'val'}, {'Rating'}, pos11, pos12,pos13, pos21, pos22, pos23, pos31, pos32, pos33, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot6pos , anot7pos, anot8pos , legpos, possible_colors([1,1], :), 'A', 0, 1 ,1, 0, 1)
            end
            accumulator = accumulator + totalheight2;
        end
        
        toSave = [folder2save, '/Clusterwise Permutation test, direct change Value and Salience population,  alphaBinom' num2str(alphaBinom), ', alphaPerm', num2str(alphaPerm)];
        print([toSave, '200'], '-dpng', '-r200');
            
        




function do_the_plot (fig1, region2,  fold, clustPermFold, clustPermFold2, DirectChangeFold, DirectChangeFold2, var2use, var2name,...
    pos11, pos12,pos13, pos21, pos22, pos23, pos31, pos32, pos33, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot6pos , anot7pos,  anot8pos , ...
    legpos, color2use, letra, xlabelon, legon, plottrialtype, plotChangeProportion, onlyPlotChange)



dbstop if error
global alpha alphaBinom alphaPerm alphaSigned type_clusterSumEffect typeclusterSumEffect_DirectChange
set(0, 'currentFigure',fig1)
set(groot,'defaultAxesLinewidth', 1)
transparency= 0.2;
AreaRat = [501, 701];
AreaA1  = [501, 601];
AreaB   = [621, 721];
AreaA2  = [741, 841];



windowstep = 10;
trialonset = 501;
tickpossize= 1000/windowstep;
leftReduction = (3600+2350)/windowstep-1;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr
rightReduction =  trialonset-(2350)/windowstep;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr

heightSigned = 0.15; % top value for the satiation effect plot
heightplot = 0.2;


var2use1 = {'rank'}; %'val', 'rank', 'valTot', 'sal', 'salTot', 'valmean', 'salmean'
var2use2= {'Salz', 'Blue',  'Green',  'Red',  'Luminance',  'Contrast', 'Hue', 'Saturation', 'Fett',  'Zucker',  'Kalorien',  'Kohlenhydrate',  'Protein', 'Ballaststoffe'};
var2use3 = {'flavour'};
var2use4 = {'val', 'sal'};



if ismember(var2use, var2use1)
    ggg=1;
elseif ismember(var2use, var2use4)
    ggg=4;
elseif ismember(var2use, var2use2)
    ggg=2;
    plotChangeProportion = 0; onlyPlotChange = 0;
elseif ismember(var2use, var2use3)
    ggg=3;
    plotChangeProportion = 0; onlyPlotChange = 0;
end



for do_single = 0%  0:2; %do_single==2 is is a single unit, do_single==1 is is a multi unit only, do_single==0 is is a  all SU and MU
    
    for restricted = 0%0:2
        
        for g=1:size(region2,1)
            
            region = region2{g};
            load(['../Rasts matfiles/patientVector_', region]);
            
            vectorPat = eval(['patientVector_', region]);

            vectorPat1 = ismember([vectorPat], [2,3,4,5,6,1,7,8,9,10,11,12]);
            vectorPat = ismember([vectorPat; vectorPat], [2,3,4,5,6,1,7,8,9,10,11,12]);
            
            
            for cc= 1: size(var2use,2)
                
                variable2 = var2use {cc};
                variable = var2use {cc};
                disp([var2use {cc}, ' region,', region, ', restricted ,',  restricted])
                
                variable2(1) = upper(variable2(1));
                
                %% define the limits of the plot
                leftReduction = (2000+2350)/windowstep-1;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr
                rightReduction =  trialonset-(2350)/windowstep;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr
                x_ax_lim = [rightReduction, trialonset + leftReduction];
                
                
                if ggg == 3
                    variable2 = [variable, 'Rat'];
                end
                toLoad = {['*', variable2,'_R']};
                if ggg == 1
                    toLoad = {['*', variable2,'_R']; ['*', variable2, 'inv_R']};
                elseif ggg == 4
                    toLoad = {['*', variable2,'_R']; ['*', variable2, 'inv_R']; ['*', variable2,'mean_R']};  % only if val or sal, take cvar from mean values for satiation analyses.
                end
                
                
                load([fold, '/regresspval_',region,'.mat'],toLoad{:});
                load([fold, '/regressreg_', region,'.mat'],toLoad{:});
                
                var=      eval(['pval', toLoad{1}(2:end)]);
                cvar=     eval(['corr', toLoad{1}(2:end)]);
                
                
                if ggg ==3
                    toLoad2 = {[variable, 'Rat_', type_clusterSumEffect]};
                else
                    toLoad2 = {[variable2, '_', type_clusterSumEffect]; [variable2, '_', typeclusterSumEffect_DirectChange]; [variable2, '_signedRankclusterSumEffect']; [variable2, '_signedRankSelclusterSumEffect']};
                end
                
                load([clustPermFold2,'/permdistribution_',region,'2.mat'],toLoad2{1:end});
                clustVar        =  eval([toLoad2{1}]);
                
                if ggg == 4 || ggg == 1
                    clustVarDC      =  eval([toLoad2{2}]);
                    clustVarSR_All  =  eval([toLoad2{3}]);
                    clustVarSR_Sel  =  eval([toLoad2{4}]);
                end
                
                
                
                if ggg == 1 || ggg == 4
                    invvar =  eval(['pval', toLoad{2}(2:end)]);
                    cinvvar=  eval(['corr', toLoad{2}(2:end)]);
                    if plotChangeProportion ||  onlyPlotChange
                        load([DirectChangeFold, 'regresspval_',region,'.mat'], ['pvalRawchangeRat_', variable2]);
                        load([DirectChangeFold, 'regressreg_',region,'.mat'],  ['corrRawchangeRat_', variable2]);
                        varchange=      eval(['pvalRawchangeRat_', variable2]);
                        cvarchange=     eval(['corrRawchangeRat_', variable2]);
                    end
                    
                    if ggg == 4
                        varmean  =  eval(['pval', toLoad{3}(2:end)]);
                        cvarmean =  eval(['corr',  toLoad{3}(2:end)]);
                    end
                end
                
                
                
                var = var (vectorPat, :);
                cvar = cvar(vectorPat, :);
                if ggg  == 1 || ggg == 4
                    invvar = invvar(vectorPat, :);
                    cinvvar= cinvvar(vectorPat, :);
                    if plotChangeProportion ||  onlyPlotChange
                        varchange = varchange(vectorPat1, :);
                        cvarchange = cvarchange(vectorPat1, :);
                    end
                    if ggg == 4
                        varmean = varmean (vectorPat, :);
                        cvarmean = cvarmean(vectorPat, :);
                    end
                end
                
                if ggg ~= 4
                    varmean = var;
                    cvarmean = cvar;
                end
                
                
                siVar= size(var,1)/2;
                siVar2 = siVar;
                
                cvar1 = cvarmean(1:siVar, :);
                cvar2 = cvarmean(1+siVar:end, :);
                
                both = ((varmean(1:siVar,:)<sqrt(2*alpha)) & (varmean(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvarmean(1:siVar,:))==sign(cvarmean(1+siVar:end,:))));
                
                if ggg  == 1 || ggg == 4 % compute mean values for satiation analyses
                    
                    cinvvar1 = cinvvar(1:siVar, :);
                    cinvvar2 = cinvvar(1+siVar:end, :);
                    
                    bothinv = ((invvar(1:siVar,:)<sqrt(2*alpha)) & (invvar(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cinvvar(1:siVar,:))==sign(cinvvar(1+siVar:end,:))));
                    selectVec =   [both] | [ bothinv]; % for selected cells
                    
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
                    
                    % varchange = [pchange];
                    % cvarchange = [cchange];
                    
                    
                    siVar3 = (sum(selectVec));
                    
                    
                    meancvar = abs(cvarmean(1:siVar,:) + cvarmean(1+siVar:end ,:))/2;
                    meancinvvar = abs(cinvvar(1:siVar,:) + cinvvar(1+siVar:end ,:))/2;
                    
                    meancinvvarSel = abs(cinvvarSel(1:siVar,:) + cinvvarSel(1+siVar:end ,:))/2;
                    meancvarSel = abs(cvarSel(1:siVar,:) + cvarSel(1+siVar:end ,:))/2;
                end
                
                
                
                
                %% ticks
                Xtickspos=[[trialonset:-tickpossize: 1], [trialonset:tickpossize:size(cvar,2)]];  Xtickspos = sort(Xtickspos(2:end)); % from the trial aligment, go forward and backwards determining te position and sort
                Xtickslabelsnames = [[0:-1:-(length([trialonset:-tickpossize: 1])-1)], [0:1:length(trialonset:tickpossize:size(cvar,2))-1]]; Xtickslabelsnames = sort(Xtickslabelsnames(2:end));
                
                for hh=1:length(Xtickslabelsnames)
                    aa{hh}=num2str(Xtickslabelsnames(hh));
                end
                Xtickslabelsnames=[]; Xtickslabelsnames=aa;
                
                
                
                
                %% plot Proportion of cells
                ax11 = subplot('position', pos11);
                set(ax11, 'units', 'centimeters', 'position', pos11)
                hold on
                
                %title('Current trial')
                
                clear kk
                for n=1:siVar
                    kk(n)=1-binocdf(n-1,siVar,alpha);
                end
                
                pos= find(kk< alphaBinom); pos= pos(1)/siVar;
                ylabel('Proportion of significant cells')
                
                
                if xlabelon ==1
                    xlabel('Time from stimulus onset (s)')
                end
                
                area(AreaRat, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
                
                
                if onlyPlotChange ||  plotChangeProportion
                    % direct change plot
                    toplotDC= varchange < alpha;
                    [toplotDC2] = create_permVector (toplotDC, clustVarDC, pos, x_ax_lim, alphaPerm);
                    plot(mean(toplotDC), 'color', [color2use(2,:), transparency], 'linewidth',1.5)
                    plot(toplotDC2, 'color', [color2use(2,:)], 'linewidth', 1.5)
                    
                end
                if onlyPlotChange ~= 1
                    %  proportion plot
                    toplot2= ((var(1:siVar,:)<sqrt(2*alpha)) & (var(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvar(1:siVar,:))==sign(cvar(1+siVar:end,:))));
                    [toplot] = create_permVector (toplot2, clustVar, pos, x_ax_lim, alphaPerm);
                    plot(sum(toplot2)./siVar, 'color', [color2use(1,:), transparency], 'linewidth',1.5)
                    plot(toplot, 'color', [color2use(1,:)], 'linewidth', 1.5)
                end
                %binomial line
                line(x_ax_lim, [alpha,alpha], 'linestyle', '--','color','k', 'linewidth', 1.3)
                line(x_ax_lim, [pos,pos], 'linestyle', '--', 'color','r', 'linewidth', 1.3)
                xlim(x_ax_lim)
                
                % settings
                ylim([0, heightplot])
                set(gca, 'Xtick',Xtickspos, 'Xticklabels',Xtickslabelsnames);
                set(gca, 'FontSize',12, 'fontweight', 'bold')
                
                
                
                %% pval plot Satiation
                
                if (ggg == 1 || ggg == 4) && (onlyPlotChange ~= 1)
                    
                    xaxis= 1:size(var,2);
                    clear pvalSign zvalues
                    for kk= 1:size(cvar,2)
                        try
                            [pvalSign(1, kk), ~, stats ] = signrank(meancvar(:,kk), meancinvvar(:,kk), 'tail', 'right', 'method','approximate');
                            zvalues (1,kk) = stats.zval;
                        catch
                            pvalSign(1,kk)= nan;
                            zvalues (1,kk) = nan;
                        end
                    end
                    
                    [toplot, CreateVector] = create_permVectorSat(pvalSign, clustVarSR_All, alphaSigned, x_ax_lim, alphaPerm, zvalues);
                    
                    ax12 = subplot('position', pos12);
                    set(ax12, 'units', 'centimeters', 'position', pos12)
                    colormap(flipud(parula));
                    assert(isequal(size(xaxis) , size(pvalSign)))
                    imagesc(xaxis,1:5,pvalSign,[0 heightSigned]);
                    xlim(x_ax_lim)
                    axis off
                    hold off
                    
                    %%  clusterwise correction
                    
                    pos14 = pos12 + [0, pos12(4)+0.005, 0, -0.5];
                    ax14 = subplot('position', pos14);
                    set(ax14, 'units', 'centimeters', 'position', pos14)
                    plot(CreateVector, 'k', 'linewidth', 2)
                    assert(isequal(size(xaxis) , size(toplot)))
                    xlim(x_ax_lim)
                    axis off
                    hold off
                    
                    
                    clear pvalSignSel zvalues stats
                    for kk= 1:size(cvar,2)
                        try
                            [pvalSignSel(1, kk), ~, stats] = signrank(meancvarSel(:,kk), meancinvvarSel(:,kk), 'tail', 'right','method','approximate');
                            zvaluesSel (1, kk) = stats.zval;
                            
                        catch
                            pvalSignSel(1, kk)= nan;
                            zvaluesSel (1, kk) = nan;
                        end
                    end
                    
                    [toplot, CreateVector] = create_permVectorSat(pvalSignSel, clustVarSR_Sel, alphaSigned, x_ax_lim, alphaPerm, zvaluesSel);
                    
                    ax13 =  subplot('position',  pos13);
                    set(ax13, 'units', 'centimeters', 'position', pos13)
                    colormap(flipud(parula));
                    assert(isequal(size(xaxis) , size(pvalSignSel)))
                    imagesc(xaxis,1:5,pvalSignSel,[0 heightSigned]);
                    xlim(x_ax_lim)
                    axis off
                    hold off
                    
                    %%  clusterwise correction
                    pos15 = pos13 + [0, -pos14(4)-0.005, 0, -0.5];
                    ax15 = subplot('position', pos15);
                    set(ax15, 'units', 'centimeters', 'position', pos15)
                    plot(CreateVector, 'k', 'linewidth', 2)
                    assert(isequal(size(xaxis) , size(toplot)))
                    xlim(x_ax_lim)
                    axis off
                    hold off
                    
                    
                    %%  temporal correlation vetween signals
                    xlimcorr= x_ax_lim(1):x_ax_lim(end);
                    sim1 = corr( sum(toplot2(:,xlimcorr))' , pvalSign(:, xlimcorr)', 'type', 'Spearman', 'rows', 'c');
                    simSel1 = corr( sum(toplot2(:,xlimcorr))' , pvalSignSel(:, xlimcorr)', 'type', 'Spearman', 'rows', 'c');
                    %         sim1DC1 = corr( sum(toplot2(:,xlimcorr))' , sum(toplotDC(:, xlimcorr))', 'type', 'Spearman', 'rows', 'c');
                    
                end
                
                
                
                
                
                
                clear  var cvar varmean  cvarmean varinv cvarinv
                %% plot A
                
                leftReduction = (3600+2350)/windowstep-1;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr
                rightReduction =  trialonset-(2350)/windowstep;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr
                
                if ~isempty(regexp (variable, 'mean'))
                    ind = regexp (variable, 'mean',  'start'): regexp (variable, 'mean',  'end');
                    variable (ind)=[];
                end
                
                
                clear corr* pval*
                
                if ggg == 4 ||  ggg== 1
                    toLoad = {['*A', variable,'_A1'];['*A', variable, 'inv_A1']};
                else
                    toLoad = {['*', variable,'A_A1']};
                end
                
                
                load([fold, '/regresspval_',region,'.mat'],toLoad{:});
                load([fold, '/regressreg_', region,'.mat'],toLoad{:});
                var=      eval(['pval', toLoad{1}(2:end)]);
                cvar=     eval(['corr', toLoad{1}(2:end)]);
                
                
                
                if ggg ==3
                    toLoad2 = {[variable, 'A_', type_clusterSumEffect]};
                    load([clustPermFold2,'/permdistribution_',region,'2.mat'],toLoad2{1:end});
                else
                    toLoad2 = {['A',variable, '_', type_clusterSumEffect]; ['A',variable, '_', typeclusterSumEffect_DirectChange]; ['A',variable, '_signedRankclusterSumEffect']; ['A',variable, '_signedRankSelclusterSumEffect']};
                    load([clustPermFold,'/permdistribution_',region,'2.mat'],toLoad2{1:end});
                end
                clustVar  =  eval([toLoad2{1}]);
                
                if ggg == 4 || ggg == 1
                    clustVarDC      =  eval([toLoad2{2}]);
                    clustVarSR_All  =  eval([toLoad2{3}]);
                    clustVarSR_Sel  =  eval([toLoad2{4}]);
                    
                    invvar =  eval(['pval', toLoad{2}(2:end)]);
                    cinvvar=  eval(['corr', toLoad{2}(2:end)]);
                    
                    if plotChangeProportion ||  onlyPlotChange
                        toLoad3 = {['*A', variable,'_RawChange_A1']};
                        load([DirectChangeFold2, 'regresspval_',region,'.mat'],toLoad3{1});
                        load([DirectChangeFold2, 'regressreg_',region,'.mat'],toLoad3{1});
                        varchange  =   eval(['pval', toLoad3{1}(2:end)]);
                        cvarchange =   eval(['corr', toLoad3{1}(2:end)]);
                    end
                end
                
                
                %% define the limits of the plot
                x_ax_lim = [0 + rightReduction, trialonset + leftReduction];
                
                
                
                
                var = var (vectorPat, :);
                cvar = cvar(vectorPat, :);
                
                if ggg == 1 || ggg == 4
                    invvar = invvar(vectorPat, :);
                    cinvvar= cinvvar(vectorPat, :);
                    if plotChangeProportion ||  onlyPlotChange
                        varchange = varchange(vectorPat1, :);
                        cvarchange = cvarchange(vectorPat1, :);
                    end
                end
                
                
                
                siVar= size(var,1)/2;
                siVar2 = siVar;
                cvar1 = cvar(1:siVar, :);
                cvar2 = cvar(1+siVar:end, :);
                both = ((var(1:siVar,:)<sqrt(2*alpha)) & (var(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvar(1:siVar,:))==sign(cvar(1+siVar:end,:))));
                
                if ggg == 1 || ggg == 4
                    cinvvar1 = cinvvar(1:siVar, :);
                    cinvvar2 = cinvvar(1+siVar:end, :);
                    
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
                    
                    meancinvvar = abs(cinvvar(1:siVar,:) + cinvvar(1+siVar:end ,:))/2;
                    meancinvvarSel = abs(cinvvarSel(1:siVar,:) + cinvvarSel(1+siVar:end ,:))/2;
                    meancvarSel = abs(cvarSel(1:siVar,:) + cvarSel(1+siVar:end ,:))/2;
                end
                
                meancvar = abs(cvar(1:siVar,:) + cvar(1+siVar:end ,:))/2;
                
                
                
                
                Xtickspos=[[trialonset:-tickpossize: 1], [trialonset:tickpossize:size(cvar,2)]];  Xtickspos = sort(Xtickspos(2:end)); % from the trial aligment, go forward and backwards determining te position and sort
                Xtickslabelsnames = [[0:-1:-(length([trialonset:-tickpossize: 1])-1)], [0:1:length(trialonset:tickpossize:size(cvar,2))-1]]; Xtickslabelsnames = sort(Xtickslabelsnames(2:end));
                
                for hh=1:length(Xtickslabelsnames)
                    aa{hh}=num2str(Xtickslabelsnames(hh));
                end
                Xtickslabelsnames=[];Xtickslabelsnames=aa;
                
                
                
                
                %% plot Proportion of cells
                
                ax21 = subplot('position', pos21);
                set(ax21, 'units', 'centimeters', 'position', pos21)
                hold on
                
                %title('Current trial')
                
                clear kk
                for n=1:siVar
                    kk(n)=1-binocdf(n-1,siVar,alpha);
                end
                
                pos= find(kk< alphaBinom); pos= pos(1)/siVar;
                %                 ylabel('Proportion of significant cells')
                
                
                
                area(AreaA1, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
                area(AreaB, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.3,0.3,0.3], 'Edgecolor', 'none')
                area(AreaA2, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
                if xlabelon ==1
                    xlabel('Time relative to stimulus A1 onset (s)')
                end
                xlim(x_ax_lim)
                
                
                %
                
                
                %  proportion plot
                
                
                if onlyPlotChange ||  plotChangeProportion
                    %    % direct change plot
                    toplotDC= varchange < alpha;
                    [toplotDC2] = create_permVector (toplotDC, clustVarDC, pos, x_ax_lim, alphaPerm);
                    plot(mean(toplotDC), 'color', [color2use(2,:), transparency], 'linewidth',1.5)
                    plot(toplotDC2, 'color', [color2use(2,:)], 'linewidth', 1.5)
                    
                end
                if onlyPlotChange ~= 1
                    toplot2= ((var(1:siVar,:)<sqrt(2*alpha)) & (var(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvar(1:siVar,:))==sign(cvar(1+siVar:end,:))));
                    [toplot] = create_permVector (toplot2, clustVar, pos, x_ax_lim, alphaPerm);
                    plot(sum(toplot2)./siVar, 'color', [color2use(1,:), transparency], 'linewidth',1.5)
                    plot(toplot, 'color', [color2use(1,:)], 'linewidth', 1.5)
                end
                %binomial line
                line(x_ax_lim, [alpha,alpha], 'linestyle', '--','color','k', 'linewidth', 1.3)
                line(x_ax_lim, [pos,pos], 'linestyle', '--', 'color','r', 'linewidth', 1.3)
                xlim(x_ax_lim)
                % settings
                ylim([0, heightplot])
                set(gca, 'Xtick',Xtickspos, 'Xticklabels',Xtickslabelsnames);
                set(gca, 'FontSize',12, 'fontweight', 'bold')
                
                
                
                
                
                
                %% pval plot
                
                
                if (ggg == 1 || ggg == 4) && (onlyPlotChange ~= 1)
                    
                    
                    %% pval plot
                    
                    xaxis= 1:size(var,2);
                    
                    
                    clear pvalSign zvalues stats
                    for kk= 1:size(cvar,2)
                        
                        try
                            [pvalSign(1, kk), ~, stats] = signrank(meancvarSel(:,kk), meancinvvarSel(:,kk), 'tail', 'right', 'method','approximate');
                            zvalues (1,kk) = stats.zval;
                        catch
                            pvalSign(1,kk)= nan;
                            zvalues (1,kk) = nan;
                        end
                        
                    end
                    
                    
                    [toplot, CreateVector] = create_permVectorSat(pvalSign, clustVarSR_All, alphaSigned, x_ax_lim, alphaPerm, zvalues);
                    
                    ax22 = subplot('position', pos22);
                    set(ax22, 'units', 'centimeters', 'position', pos22)
                    colormap(flipud(parula));
                    assert(isequal(size(xaxis) , size(pvalSign)))
                    imagesc(xaxis,1:5,pvalSign,[0 heightSigned]);
                    xlim(x_ax_lim)
                    axis off
                    hold off
                    
                    %  clusterwise correction
                    
                    pos24 = pos22 + [0, pos22(4)+0.005, 0, -0.5];
                    ax24 = subplot('position', pos24);
                    set(ax24, 'units', 'centimeters', 'position', pos24)
                    plot(CreateVector, 'k', 'linewidth', 2)
                    assert(isequal(size(xaxis) , size(toplot)))
                    xlim(x_ax_lim)
                    axis off
                    hold off
                    
                    
                    clear pvalSignSel zvalues stats
                    for kk= 1:size(cvar,2)
                        try
                            [pvalSignSel(1, kk), ~, stats] = signrank(meancvarSel(:,kk), meancinvvarSel(:,kk), 'tail', 'right', 'method','approximate');
                            zvaluesSel (1, kk) = stats.zval;
                            
                        catch
                            pvalSignSel(1, kk)= nan;
                            zvaluesSel (1, kk) = nan;
                        end
                    end
                    
                    [toplot, CreateVector] = create_permVectorSat(pvalSignSel, clustVarSR_Sel, alphaSigned, x_ax_lim, alphaPerm, zvaluesSel);
                    
                    ax23 =  subplot('position',  pos23);
                    set(ax23, 'units', 'centimeters', 'position', pos23)
                    colormap(flipud(parula));
                    assert(isequal(size(xaxis) , size(pvalSignSel)))
                    imagesc(xaxis,1:5,pvalSignSel,[0 heightSigned]);
                    xlim(x_ax_lim)
                    axis off
                    hold off
                    
                    %%  clusterwise correction
                    pos25 = pos23 + [0, -pos24(4)-0.005, 0, -0.5];
                    ax25 = subplot('position', pos25);
                    set(ax25, 'units', 'centimeters', 'position', pos25)
                    plot(CreateVector, 'k', 'linewidth', 2)
                    assert(isequal(size(xaxis) , size(toplot)))
                    xlim(x_ax_lim)
                    axis off
                    hold off
                    
                    
                    %%  temporal correlation vetween signals
                    xlimcorr= x_ax_lim(1):x_ax_lim(end);
                    sim2 = corr( sum(toplot2(:,xlimcorr))' , pvalSign(:, xlimcorr)', 'type', 'Spearman', 'rows', 'c');
                    simSel2 = corr( sum(toplot2(:,xlimcorr))' , pvalSignSel(:, xlimcorr)', 'type', 'Spearman', 'rows', 'c');
                    %         sim1DC2 = corr( sum(toplot2(:,xlimcorr))' , sum(toplotDC(:, xlimcorr))', 'type', 'Spearman', 'rows', 'c');
                end
                
                
                
                
                
                
                
                %% plot B
                
                
                clear  var cvar varmean  cvarmean varinv cvarinv
                
                
                
                
                leftReduction = (3600+2350)/windowstep-1;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr
                rightReduction =  trialonset-(2350)/windowstep;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr
                x_ax_lim = [0 + rightReduction, trialonset + leftReduction];
                
                if ~isempty(regexp (variable, 'mean'))
                    ind = regexp (variable, 'mean',  'start'): regexp (variable, 'mean',  'end');
                    variable (ind)=[];
                end
                
                
                clear corr* pval*
                
                if ggg== 4 ||  ggg== 1
                    toLoad = {['*B', variable,'_A1'];['*B', variable, 'inv_A1']};
                else
                    toLoad = {['*', variable,'B_A1']};
                end
                
                
                load([fold, '/regresspval_',region,'.mat'],toLoad{:});
                load([fold, '/regressreg_', region,'.mat'],toLoad{:});
                var=      eval(['pval', toLoad{1}(2:end)]);
                cvar=     eval(['corr', toLoad{1}(2:end)]);
                
                
                
                %% define the limits of the plot
                
                if ggg ==3
                    toLoad2 = {[variable, 'B_', type_clusterSumEffect]};
                    load([clustPermFold2,'/permdistribution_',region,'2.mat'],toLoad2{1:end});
                else
                    toLoad2 = {['B',variable, '_', type_clusterSumEffect]; ['B',variable, '_', typeclusterSumEffect_DirectChange]; ['B',variable, '_signedRankclusterSumEffect']; ['B',variable, '_signedRankSelclusterSumEffect']};
                    load([clustPermFold,'/permdistribution_',region,'2.mat'],toLoad2{1:end});
                end
                
                clustVar        =  eval([toLoad2{1}]);
                
                
                if ggg == 4 || ggg == 1
                    clustVarDC      =  eval([toLoad2{2}]);
                    clustVarSR_All  =  eval([toLoad2{3}]);
                    clustVarSR_Sel  =  eval([toLoad2{4}]);
                    
                    invvar =  eval(['pval', toLoad{2}(2:end)]);
                    cinvvar=  eval(['corr', toLoad{2}(2:end)]);
                    
                    if plotChangeProportion ||  onlyPlotChange
                        toLoad3 = {['*B', variable,'_RawChange_A1']};
                        load([DirectChangeFold2, 'regresspval_',region,'.mat'],toLoad3{1});
                        load([DirectChangeFold2, 'regressreg_',region,'.mat'],toLoad3{1});
                        varchange  =   eval(['pval', toLoad3{1}(2:end)]);
                        cvarchange =   eval(['corr', toLoad3{1}(2:end)]);
                    end
                    
                end
                
                
                var = var (vectorPat, :);
                cvar = cvar(vectorPat, :);
                
                if ggg == 1 || ggg == 4
                    invvar = invvar(vectorPat, :);
                    cinvvar= cinvvar(vectorPat, :);
                    if plotChangeProportion ||  onlyPlotChange
                        varchange = varchange(vectorPat1, :);
                        cvarchange = cvarchange(vectorPat1, :);
                    end
                end
                
                
                
                siVar= size(var,1)/2;
                siVar2 = siVar;
                cvar1 = cvar(1:siVar, :);
                cvar2 = cvar(1+siVar:end, :);
                both = ((var(1:siVar,:)<sqrt(2*alpha)) & (var(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvar(1:siVar,:))==sign(cvar(1+siVar:end,:))));
                
                
                if ggg == 1 || ggg == 4
                    cinvvar1 = cinvvar(1:siVar, :);
                    cinvvar2 = cinvvar(1+siVar:end, :);
                    
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
                    
                    meancinvvar = abs(cinvvar(1:siVar,:) + cinvvar(1+siVar:end ,:))/2;
                    meancinvvarSel = abs(cinvvarSel(1:siVar,:) + cinvvarSel(1+siVar:end ,:))/2;
                    meancvarSel = abs(cvarSel(1:siVar,:) + cvarSel(1+siVar:end ,:))/2;
                end
                
                meancvar = abs(cvar(1:siVar,:) + cvar(1+siVar:end ,:))/2;
                
                
                
                
                
                
                Xtickspos=[[trialonset:-tickpossize: 1], [trialonset:tickpossize:size(cvar,2)]];  Xtickspos = sort(Xtickspos(2:end)); % from the trial aligment, go forward and backwards determining te position and sort
                Xtickslabelsnames = [[0:-1:-(length([trialonset:-tickpossize: 1])-1)], [0:1:length(trialonset:tickpossize:size(cvar,2))-1]]; Xtickslabelsnames = sort(Xtickslabelsnames(2:end));
                
                for hh=1:length(Xtickslabelsnames)
                    aa{hh}=num2str(Xtickslabelsnames(hh));
                end
                Xtickslabelsnames=[];Xtickslabelsnames=aa;
                
                
                
                
                %% plot Proportion of cells
                
                ax31 = subplot('position', pos31);
                set(ax31, 'units', 'centimeters', 'position', pos31)
                hold on
                
                %title('Current trial')
                
                clear kk
                for n=1:siVar
                    kk(n)=1-binocdf(n-1,siVar,alpha);
                end
                
                pos= find(kk< alphaBinom); pos= pos(1)/siVar;
                %                 ylabel('Proportion of significant cells')
                line(x_ax_lim, [alpha,alpha], 'linestyle', '--','color','k', 'linewidth', 1.3)
                line(x_ax_lim, [pos,pos], 'linestyle', '--', 'color','r', 'linewidth', 1.3)
                xlim(x_ax_lim)
                
                
                if onlyPlotChange ||  plotChangeProportion
                    toplotDC= varchange < alpha;
                    [toplotDC2] = create_permVector (toplotDC, clustVarDC, pos, x_ax_lim, alphaPerm);
                    plot(toplotDC2, 'color', [color2use(2,:)], 'linewidth', 1.5)
                end
                if onlyPlotChange ~= 1
                    toplot2= ((var(1:siVar,:)<sqrt(2*alpha)) & (var(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvar(1:siVar,:))==sign(cvar(1+siVar:end,:))));
                    [toplot] = create_permVector (toplot2, clustVar, pos, x_ax_lim, alphaPerm);
                    plot(toplot, 'color', [color2use(1,:)], 'linewidth', 1.5)
                    
                end
                
                
                area(AreaA1, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
                area(AreaB, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.3,0.3,0.3], 'Edgecolor', 'none')
                area(AreaA2, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
                if xlabelon ==1
                    xlabel('Time relative to stimulus A1 onset (s)')
                end
                xlim(x_ax_lim)
                
                
                
                if onlyPlotChange ||  plotChangeProportion
                    %    % direct change plot
                    toplotDC= varchange < alpha;
                    [toplotDC2] = create_permVector (toplotDC, clustVarDC, pos, x_ax_lim, alphaPerm);
                    plot(mean(toplotDC), 'color', [color2use(2,:), transparency], 'linewidth',1.5)
                    plot(toplotDC2, 'color', [color2use(2,:)], 'linewidth', 1.5)
                    
                end
                if onlyPlotChange ~= 1
                    %  proportion plot
                    
                    toplot2= ((var(1:siVar,:)<sqrt(2*alpha)) & (var(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvar(1:siVar,:))==sign(cvar(1+siVar:end,:))));
                    plot(sum(toplot2)./siVar, 'color', [color2use(1,:), transparency], 'linewidth',1.5)
                    [toplot] = create_permVector (toplot2, clustVar, pos, x_ax_lim, alphaPerm);
                    plot(toplot, 'color', [color2use(1,:)], 'linewidth', 1.5)
                end
                
                
                
                %binomial line
                line(x_ax_lim, [alpha,alpha], 'linestyle', '--','color','k', 'linewidth', 1.3)
                line(x_ax_lim, [pos,pos], 'linestyle', '--', 'color','r', 'linewidth', 1.3)
                xlim(x_ax_lim)
                
                % settings
                ylim([0, heightplot])
                set(gca, 'Xtick',Xtickspos, 'Xticklabels',Xtickslabelsnames);
                set(gca, 'FontSize',12, 'fontweight', 'bold')
                
                
                
                
                if legon ==1
                    if (ggg == 1 || ggg == 4) && (plotChangeProportion == 1) && onlyPlotChange == 0
                        leg =legend('Chance', 'Significance', 'Change', var2name{cc}); % [var2name{cc},'Change']
                    else
                        leg =legend('Chance', 'Significance',  var2name{cc});
                    end
                    set(leg, 'units','centimeters', 'Position', legpos ,'FontSize',12);
                end
                
                
                %% pval plot
                
                if (ggg == 1 || ggg == 4) && (onlyPlotChange ~= 1)
                    
                    %% pval plot
                    xaxis= 1:size(var,2);
                    clear pvalSign zvalues
                    for kk= 1:size(cvar,2)
                        try
                            [pvalSign(1, kk), ~, stats ] = signrank(meancvar(:,kk), meancinvvar(:,kk), 'tail', 'right', 'method','approximate');
                            zvalues (1,kk) = stats.zval;
                        catch
                            pvalSign(1,kk)= nan;
                            zvalues (1,kk) = nan;
                        end
                    end
                    
                    [toplot, CreateVector] = create_permVectorSat(pvalSign, clustVarSR_All, alphaSigned, x_ax_lim, alphaPerm, zvalues);
                    
                    ax32 = subplot('position', pos32);
                    set(ax32, 'units', 'centimeters', 'position', pos32)
                    colormap(flipud(parula));
                    assert(isequal(size(xaxis) , size(pvalSign)))
                    imagesc(xaxis,1:5,pvalSign,[0 heightSigned]);
                    xlim(x_ax_lim)
                    axis off
                    hold off
                    
                    %%  clusterwise correction
                    
                    pos34 = pos32 + [0, pos32(4)+0.005, 0, -0.5];
                    ax34 = subplot('position', pos34);
                    set(ax34, 'units', 'centimeters', 'position', pos34)
                    plot(CreateVector, 'k', 'linewidth', 2)
                    assert(isequal(size(xaxis) , size(toplot)))
                    xlim(x_ax_lim)
                    axis off
                    hold off
                    
                    
                    clear pvalSignSel zvalues stats
                    for kk= 1:size(cvar,2)
                        try
                            [pvalSignSel(1, kk), ~, stats] = signrank(meancvarSel(:,kk), meancinvvarSel(:,kk), 'tail', 'right', 'method','approximate');
                            zvaluesSel (1, kk) = stats.zval;
                            
                        catch
                            pvalSignSel(1, kk)= nan;
                            zvaluesSel (1, kk) = nan;
                        end
                    end
                    [toplot, CreateVector] = create_permVectorSat(pvalSignSel, clustVarSR_Sel, alphaSigned, x_ax_lim, alphaPerm, zvaluesSel);
                    
                    
                    ax33 =  subplot('position',  pos33);
                    set(ax33, 'units', 'centimeters', 'position', pos33)
                    colormap(flipud(parula));
                    assert(isequal(size(xaxis) , size(pvalSignSel)))
                    imagesc(xaxis,1:5,pvalSignSel,[0 heightSigned]);
                    xlim(x_ax_lim)
                    axis off
                    hold off
                    
                    %%  clusterwise correction
                    pos35 = pos33 + [0, -pos34(4)-0.005, 0, -0.5];
                    ax35 = subplot('position', pos35);
                    set(ax35, 'units', 'centimeters', 'position', pos35)
                    plot(CreateVector, 'k', 'linewidth', 2)
                    assert(isequal(size(xaxis) , size(toplot)))
                    xlim(x_ax_lim)
                    axis off
                    hold off
                    
                    
                    %%  temporal correlation vetween signals
                    xlimcorr= x_ax_lim(1):x_ax_lim(end);
                    sim3 = corr( sum(toplot2(:,xlimcorr))' , pvalSign(:, xlimcorr)', 'type', 'Spearman', 'rows', 'c');
                    simSel3 = corr( sum(toplot2(:,xlimcorr))' , pvalSignSel(:, xlimcorr)', 'type', 'Spearman', 'rows', 'c');
                    %         sim1DC2 = corr( sum(toplot2(:,xlimcorr))' , sum(toplotDC(:, xlimcorr))', 'type', 'Spearman', 'rows', 'c');
                    
                    
                    
                    
                    disp(['mean correlation over time between fraction of significant neurons and Satiation according to All SignedRank =  ' , num2str((sim1 + sim2 + sim3)/3)])
                    disp(['mean correlation over time between fraction of significant neurons and Satiation according to Selected SignedRank =  ' , num2str((simSel1 + simSel2 + simSel3)/3)])
                    %          disp(['mean correlation over time between fraction of significant neurons and Satiation according to Direct Change =  ' , num2str((sim1DC1 + sim1DC2)/2)])
                    
                end
                
                
                
                
                
                
                
                global sizeletra sizetrialtype
                anot5 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',[letra],'Units','normalized','FontName','Helvetica','FontSize',sizeletra,'EdgeColor','none', 'fontweight', 'bold');
                
                
                
                
                if (ggg == 1 || ggg == 4) && (onlyPlotChange ~= 1)
                    
                    anot2 = annotation('textbox',[0.02      0 + 0.02   0.3    0.04],   'string',['All'],'Units','normalized','FontName','Helvetica','FontSize',12,'EdgeColor','none', 'fontweight', 'bold');
                    anot3 = annotation('textbox',[0.0005    0 + 0.02   0.3    0.04],   'string',['Selected'],'Units','normalized','FontName','Helvetica','FontSize',12,'EdgeColor','none', 'fontweight', 'bold');
                    
                    set(anot2, 'units', 'centimeters', 'position', anot2pos)
                    set(anot3, 'units', 'centimeters', 'position', anot3pos)
                    
                    cbh = colorbar;
                    cbh.FontSize= 11;
                    cbh.Units = 'centimeters';
                    cbh.Position= cbhpos;
                    cbh.FontWeight= 'bold';
                    anot4 = annotation('textbox', [0,0,0.1,0.1],   'string',['pval'],'Units','normalized','FontName','Helvetica','FontSize',11,'EdgeColor','none','fontweight', 'bold');
                    
                    if xlabelon ~=1
                        set(anot4, 'units', 'centimeters', 'position', anot4pos)
                    end
                end
                
                
                if plottrialtype ==1
                    anot6 = annotation('textbox',[0.00003    0+0.02   0.7    0.04],   'string',['Stimulus Valuation'],'Units','normalized','FontName','Helvetica','FontSize',sizetrialtype,'EdgeColor','none', 'fontweight', 'bold');
                    anot7 = annotation('textbox',[0.00003    0+0.02   0.3    0.04],   'string',['Stimulus A'],'Units','normalized','FontName','Helvetica','FontSize',sizetrialtype,'EdgeColor','none', 'fontweight', 'bold');
                    anot8 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',['Stimulus B'],'Units','normalized','FontName','Helvetica','FontSize',sizetrialtype,'EdgeColor','none', 'fontweight', 'bold');
                    set(anot6, 'units', 'centimeters', 'position', anot6pos)
                    set(anot7, 'units', 'centimeters', 'position', anot7pos)
                    set(anot8, 'units', 'centimeters', 'position', anot8pos)
                end
                
                
                
                set(anot5, 'units', 'centimeters', 'position', anot5pos)
                
                
                
                %disp(' corrAll'), (sim1+sim2+sim3)/3,
                %disp('corrSel'), (simSel1+simSel2+simSel3)/3,
                
            end
        end
    end
end







function do_the_plot2 (fig1, region2, DirectChangeFold, fold, alphaBinom, var2use, var2name, ...
    pos11, pos12,pos13, pos21, pos22, pos23, pos31, pos32, pos33, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot6pos , anot7pos,  anot8pos , legpos, color2use, letra, xlabelon, legon, plottrialtype, plotChangeProportion, onlyPlotChange)


alpha = 0.05;
set(0, 'currentFigure',fig1)
dbstop if error

set(groot,'defaultAxesLinewidth', 1)

AreaRat = [501, 701];
AreaA1  = [501, 601];
AreaB   = [621, 721];
AreaA2  = [741, 841];

windowstep = 10;
trialonset = 501;
tickpossize= 1000/windowstep;


heightSigned = 0.15; % top value for the satiation effect plot
heightplot = 0.2;

var2use1 = {'rank'}; %'val', 'rank', 'valTot', 'sal', 'salTot', 'valmean', 'salmean'
var2use2= {'Salz', 'Blue',  'Green',  'Red',  'Luminance',  'Contrast', 'Hue', 'Saturation', 'Fett',  'Zucker',  'Kalorien',  'Kohlenhydrate',  'Protein', 'Ballaststoffe'};
var2use3 = {'flavour'};
var2use4 = {'val', 'sal'};


if ismember(var2use, var2use1)
    ggg=1;
elseif ismember(var2use, var2use4)
    ggg=4;
elseif ismember(var2use, var2use2)
    ggg=2;
    plotChangeProportion = 0; onlyPlotChange = 0;
elseif ismember(var2use, var2use3)
    ggg=3;
    plotChangeProportion = 0; onlyPlotChange = 0;
end



for do_single = 0%  0:2; %do_single==2 is is a single unit, do_single==1 is is a multi unit only, do_single==0 is is a  all SU and MU
    
    for restricted = 0%0:2
        
        for g=1:size(region2,1)
            
            region = region2{g};
            load(['../Rasts matfiles/patientVector_', region]);
            
            vectorPat = eval(['patientVector_', region]);
            vectorPat1 = ismember([vectorPat], 2,3,4,5,6,1,7,8,9,10,11,12);
            vectorPat = ismember([vectorPat; vectorPat], 2,3,4,5,6,1,7,8,9,10,11,12);
            
            for cc= 1: size(var2use,2)
                
                
                variable2 = var2use {cc};
                variable = var2use {cc};
                disp([var2use {cc}, ' region,', region, ', restricted ,',  restricted])
                
                toLoad = {['*', variable, '_Rat_anova']};
                
                
                load([fold, '/regresspval_',region,'.mat'],toLoad{:});
                load([fold, '/regressreg_', region,'.mat'],toLoad{:});
                var=      eval(['pval', toLoad{1}(2:end)]);
                cvar=     eval(['corr', toLoad{1}(2:end)]);
                
                leftReduction = (2000+2350)/windowstep-1;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr
                rightReduction =  trialonset-(2350)/windowstep;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr
                
                
                %% define the limits of the plot
                x_ax_lim = [rightReduction, trialonset + leftReduction];
                
                siVar= size(var,1)/2;
                siVar2 = siVar;
                toplot2 = (var(1:siVar,:)<sqrt(alpha))   &   (var(siVar+1:end,:)<sqrt(alpha));
                
                
                %% ticks
                
                Xtickspos=[[trialonset:-tickpossize: 1], [trialonset:tickpossize:size(cvar,2)]];  Xtickspos = sort(Xtickspos(2:end)); % from the trial aligment, go forward and backwards determining te position and sort
                Xtickslabelsnames = [[0:-1:-(length([trialonset:-tickpossize: 1])-1)], [0:1:length(trialonset:tickpossize:size(cvar,2))-1]]; Xtickslabelsnames = sort(Xtickslabelsnames(2:end));
                
                for hh=1:length(Xtickslabelsnames)
                    aa{hh}=num2str(Xtickslabelsnames(hh));
                end
                Xtickslabelsnames=[]; Xtickslabelsnames=aa;
                
                
                
                
                
                
                %% plot Proportion of cells
                
                ax11 = subplot('position', pos11);
                set(ax11, 'units', 'centimeters', 'position', pos11)
                hold on
                
                
                clear kk
                for n=1:siVar
                    kk(n)=1-binocdf(n-1,siVar,alpha);
                end
                pos= find(kk< alphaBinom); pos= pos(1)/siVar;
                
                
                ylabel('Proportion of significant cells')
                if xlabelon ==1
                    xlabel('Time from stimulus onset (s)')
                end
                
                
                area(AreaRat, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
                
                if onlyPlotChange ~= 1
                    plot(sum(toplot2)./siVar, 'color', color2use(1,:), 'linewidth', 1.5)
                end
                line(x_ax_lim, [alpha,alpha], 'linestyle', '--','color','k', 'linewidth', 1.3)
                line(x_ax_lim, [pos,pos], 'linestyle', '--', 'color','r', 'linewidth', 1.3)
                xlim(x_ax_lim)
                
                ylim([0, heightplot])
                set(gca, 'Xtick',Xtickspos, 'Xticklabels',Xtickslabelsnames);
                set(gca, 'FontSize',12, 'fontweight', 'bold')
                
                
                
                
                
                leftReduction = (3600+2350)/windowstep-1;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr
                rightReduction =  trialonset-(2350)/windowstep;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr
                
                
                
                
                
                %% plot A
                
                clear corr* pval* var cvar
                toLoad = {['*', variable, '_A_anova']};
                load([fold, '/regresspval_',region,'.mat'],toLoad{:});
                load([fold, '/regressreg_', region,'.mat'],toLoad{:});
                
                var =     eval(['pval', toLoad{1}(2:end)]);
                cvar=     eval(['corr', toLoad{1}(2:end)]);
                
                
                
                
                %% define the limits of the plot
                x_ax_lim = [rightReduction, trialonset + leftReduction];
                
                siVar= size(var,1)/2;
                siVar2 = siVar;
                cvar1 = cvar(1:siVar, :);
                cvar2 = cvar(1+siVar:end, :);
                toplot2 = (var(1:siVar,:)<sqrt(alpha))   &   (var(siVar+1:end,:)<sqrt(alpha));
                
                
                Xtickspos=[[trialonset:-tickpossize: 1], [trialonset:tickpossize:size(cvar,2)]];  Xtickspos = sort(Xtickspos(2:end)); % from the trial aligment, go forward and backwards determining te position and sort
                Xtickslabelsnames = [[0:-1:-(length([trialonset:-tickpossize: 1])-1)], [0:1:length(trialonset:tickpossize:size(cvar,2))-1]]; Xtickslabelsnames = sort(Xtickslabelsnames(2:end));
                
                
                for hh=1:length(Xtickslabelsnames)
                    aa{hh}=num2str(Xtickslabelsnames(hh));
                end
                Xtickslabelsnames=[];Xtickslabelsnames=aa;
                
                
                
                
                
                ax21 = subplot('position', pos21);
                set(ax21, 'units', 'centimeters', 'position', pos21)
                hold on
                
                
                clear kk
                for n=1:siVar
                    kk(n)=1-binocdf(n-1,siVar,alpha);
                end
                
                pos= find(kk< alphaBinom); pos= pos(1)/siVar;
                
                if xlabelon == 1
                    xlabel('Time relative to stimulus A1 onset (s)')
                end
                
                
                
                area(AreaA1, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
                area(AreaB, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.3,0.3,0.3], 'Edgecolor', 'none')
                area(AreaA2, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
                
                
                if onlyPlotChange ~= 1
                    plot(sum(toplot2)./siVar, 'color', color2use(1,:), 'linewidth', 1.5)
                end
                line(x_ax_lim, [alpha,alpha], 'linestyle', '--','color','k', 'linewidth', 1.3)
                line(x_ax_lim, [pos,pos], 'linestyle', '--', 'color','r', 'linewidth', 1.3)
                
                
                xlim(x_ax_lim)
                ylim([0, heightplot])
                set(gca, 'Xtick',Xtickspos, 'Xticklabels',Xtickslabelsnames);
                set(gca, 'FontSize',12, 'fontweight', 'bold')
                
                
                
                
                
                %% plot B
                
                clear corr* pval* var cvar
                toLoad = {['*', variable, '_B_anova']};
                load([fold, '/regresspval_',region,'.mat'],toLoad{:});
                load([fold, '/regressreg_', region,'.mat'],toLoad{:});
                
                var=      eval(['pval', toLoad{1}(2:end)]);
                cvar=     eval(['corr', toLoad{1}(2:end)]);
                
                
                
                %% define the limits of the plot
                x_ax_lim = [rightReduction, trialonset + leftReduction];
                
                
                
                siVar= size(var,1)/2;
                siVar2 = siVar;
                cvar1 = cvar(1:siVar, :);
                cvar2 = cvar(1+siVar:end, :);
                toplot2 = (var(1:siVar,:)<sqrt(alpha))   &   (var(siVar+1:end,:)<sqrt(alpha));
                
                
                
                
                Xtickspos=[[trialonset:-tickpossize: 1], [trialonset:tickpossize:size(cvar,2)]];  Xtickspos = sort(Xtickspos(2:end)); % from the trial aligment, go forward and backwards determining te position and sort
                Xtickslabelsnames = [[0:-1:-(length([trialonset:-tickpossize: 1])-1)], [0:1:length(trialonset:tickpossize:size(cvar,2))-1]]; Xtickslabelsnames = sort(Xtickslabelsnames(2:end));
                
                
                for hh=1:length(Xtickslabelsnames)
                    aa{hh}=num2str(Xtickslabelsnames(hh));
                end
                Xtickslabelsnames=[];Xtickslabelsnames=aa;
                
                
                
                
                
                ax31 = subplot('position', pos31);
                set(ax31, 'units', 'centimeters', 'position', pos31)
                hold on
                
                
                clear kk
                for n=1:siVar
                    kk(n)=1-binocdf(n-1,siVar,alpha);
                end
                
                pos= find(kk< alphaBinom); pos= pos(1)/siVar;
                
                if xlabelon == 1
                    xlabel('Time relative to stimulus A1 onset (s)')
                end
                
                
                
                area(AreaA1, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
                area(AreaB, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.3,0.3,0.3], 'Edgecolor', 'none')
                area(AreaA2, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
                
                
                if onlyPlotChange ~= 1
                    plot(sum(toplot2)./siVar, 'color', color2use(1,:), 'linewidth', 1.5)
                end
                line(x_ax_lim, [alpha,alpha], 'linestyle', '--','color','k', 'linewidth', 1.3)
                line(x_ax_lim, [pos,pos], 'linestyle', '--', 'color','r', 'linewidth', 1.3)
                
                
                xlim(x_ax_lim)
                ylim([0, heightplot])
                set(gca, 'Xtick',Xtickspos, 'Xticklabels',Xtickslabelsnames);
                set(gca, 'FontSize',12, 'fontweight', 'bold')
                
                
                
                
                
                global sizeletra sizetrialtype
                anot5 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',[letra],'Units','normalized','FontName','Helvetica','FontSize',sizeletra,'EdgeColor','none', 'fontweight', 'bold');
                set(anot5, 'units', 'centimeters', 'position', anot5pos)
                
                
                if plottrialtype ==1
                    anot6 = annotation('textbox',[0.00003    0+0.02   0.7    0.04],   'string',['Stimulus Valuation'],'Units','normalized','FontName','Helvetica','FontSize',sizetrialtype,'EdgeColor','none', 'fontweight', 'bold');
                    anot7 = annotation('textbox',[0.00003    0+0.02   0.3    0.04],   'string',['Stimulus A'],'Units','normalized','FontName','Helvetica','FontSize',sizetrialtype,'EdgeColor','none', 'fontweight', 'bold');
                    anot8 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',['Stimulus B'],'Units','normalized','FontName','Helvetica','FontSize',sizetrialtype,'EdgeColor','none', 'fontweight', 'bold');
                    set(anot6, 'units', 'centimeters', 'position', anot6pos)
                    set(anot7, 'units', 'centimeters', 'position', anot7pos)
                    set(anot8, 'units', 'centimeters', 'position', anot8pos)
                end
                
                
            end
        end
    end
end













function [toplot] = create_permVector (toplot2, clustVar, pos, x_ax_lim, alphaPerm)

toplot = mean(toplot2(:,x_ax_lim(1):x_ax_lim(2)),1);
cluster = ([zeros(1,1,1), toplot>= pos, zeros(1,1,1)] );
findStart = (find(diff(cluster,1)==1));
findEnd = (find(diff(cluster,1)==-1)-1);


assert (all((findEnd - findStart +1)>0))
if ~isempty(findEnd - findStart +1)
    [clustersize, indsize] = sort(findEnd - findStart +1, 'descend');  % find the biggest cluster of consecutive signifcance
    clusSumEffect =[];
    for nn = 1: length(findStart)
        clusSumEffect (nn) = sum(toplot(:,findStart(nn):findEnd(nn)));
    end
    [clusterSumEffect, indsizeSum] = sort(clusSumEffect, 'descend');  % find the biggest cluster of consecutive significance times effect size
    
    SortedClusPerm= sort(clustVar, 'descend');
    indSigClus = find(clusterSumEffect > SortedClusPerm(alphaPerm*1000));
    indClus2 =  indsizeSum(indSigClus);
    
else
    clustersize = 0;
    clusterSumEffect = 0;
end

CrateVector = zeros(size(toplot));
if ~isempty(findEnd - findStart +1)
    for nn = 1: length(indClus2)
        CrateVector(findStart(indClus2(nn)): findEnd(indClus2(nn))) = 1;
    end
end

toplot(~CrateVector) = nan;
toplot = [nan(1, x_ax_lim(1)-1), toplot, nan(1,size(toplot2,2) - x_ax_lim(2))];






function [toplot, CreateVector] = create_permVectorSat (pvaluesSR, clustVar, pos, x_ax_lim, alphaPerm, zvaluesSR)

pvaluesSR= pvaluesSR(:)';
toplot = pvaluesSR(:,x_ax_lim(1):x_ax_lim(2));

toplot2 = zvaluesSR(:,x_ax_lim(1):x_ax_lim(2));

cluster = ([zeros(1,1,1), toplot<pos, zeros(1,1,1)] );
findStart = (find(diff(cluster,1)==1));
findEnd = (find(diff(cluster,1)==-1)-1);


assert (all((findEnd - findStart +1)>0))
if ~isempty(findEnd - findStart +1)
    [clustersize, indsize] = sort(findEnd - findStart +1, 'descend');  % find the biggest cluster of consecutive signifcance
    clusSumEffect =[];
    for nn = 1: length(findStart)
        clusSumEffect (nn) = sum(toplot2(:,findStart(nn):findEnd(nn)));
    end
    [clusterSumEffect, indsizeSum] = sort(clusSumEffect, 'descend');  % find the biggest cluster of consecutive significance times effect size
    
    SortedClusPerm= sort(clustVar, 'descend');
    indSigClus = find(clusterSumEffect > SortedClusPerm(alphaPerm*1000));
    indClus2 =  indsizeSum(indSigClus);
    
else
    clustersize = 0;
    clusterSumEffect = 0;
end

CreateVector = zeros(size(toplot));
if ~isempty(findEnd - findStart +1)
    for nn = 1: length(indClus2)
        CreateVector(findStart(indClus2(nn)): findEnd(indClus2(nn))) = 1;
    end
end

toplot(~CreateVector) = 1;
toplot = [ones(1, x_ax_lim(1)-1), toplot, ones(1,size(pvaluesSR,2) - x_ax_lim(2))];
CreateVector= [zeros(1, x_ax_lim(1)-1), CreateVector, zeros(1,size(pvaluesSR,2) - x_ax_lim(2))];
CreateVector (CreateVector==0)=nan;


















