function generateplot4_clusterwisePermutation


%% this function plots RT and confidence plots aligned to trial onset
% Figure 4, Figure S3

dbstop if error


global sizeletra sizetrialtype type_clusterSumEffect typeclusterSumEffect_DirectChange
global alpha alphaBinom alphaPerm alphaSigned

folder2save =  '/media/Projects/Alex/Reclustered analysis All/PaperPlots'; mkdir(folder2save)
autocorr_fold =  '/media/Projects/Alex/Reclustered analysis All/mixed_model1_001Spearman_window10_kernel200ms_step10ms__InversesVisualNutritional_autocorr/';
fold =  '/media/Projects/Alex/Reclustered analysis All/mixed_model1_Spearman_window10_kernel200ms_step10ms__InversesVisualNutritional/';
cd(fold)


clustPermFold = '/media/raid/Alex/Reclustered analysis All/clusterwise_permutation_1_Spearman_window10_kernel200ms_step10ms_ConfidenceRTAbsDiff_correctedWindow';
clustPermFold_auto = '/media/raid/Alex/Reclustered analysis All/clusterwise_permutation_1_Spearman_window10_kernel200ms_step10ms_ConfidenceRTAbsDiff_correctedWindow_autocorr';
DirectChangeFold_auto = '/media/raid/Alex/Reclustered analysis All/mixed_model1_Spearman_window10_kernel200ms_step10ms_ChangeInVariables_autocorr/';
DirectChangeFold = '/media/raid/Alex/Reclustered analysis All/mixed_model1_Spearman_window10_kernel200ms_step10ms_ChangeInVariables/';
clustPermFold2 = '/media/raid/Alex/Reclustered analysis All/clusterwise_permutation_1_Spearman_window10_kernel200ms_step10ms_TasteValSalRT_correctedWindow';
clustPermFold2_auto = '/media/raid/Alex/Reclustered analysis All/clusterwise_permutation_1_Spearman_window10_kernel200ms_step10ms_TasteValSalRT_correctedWindow_autocorr';


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





%% Figure for reviewer, not in paper

accumulator = 0;
pos1x = pos1xor;
pos2x = pos1x + diffx1x2;

figpos = [0,0, 30.4800+shiftall, 25.4];
close all
fig1 =  figure('color', 'w', 'units', 'centimeters', 'position', figpos, 'visible', 'off');

mainystart = 3.5;
for  figures = 1:2
    
    mainy = mainystart + accumulator;
    sub2y=  mainy - sub2yadd;
    sub3y=  mainy - sub3yadd;
    
    if  figures == 1
        sub2y=  mainy - sub2yadd - 0.2;
        sub3y=  mainy - sub3yadd - 0.2;
    end
    
    pos21    =  [pos2x    mainy  width    heightmain];
    pos22    =  [pos2x    sub2y  width    heightTicks];
    pos23    =  [pos2x    sub3y  width    heightTicks];
    pos11    =  [pos1x    mainy  width    heightmain];
    pos12    =  [pos1x    sub2y  width    heightTicks];
    pos13    =  [pos1x    sub3y  width    heightTicks];
    
    
    
    % letter for figure
    heightanotA = mainy + heightmain + shiftletra;
    anot5pos = [pos1x-shiftletrax, heightanotA, 4.5720 0.4286];
    
    
    % current next trial, size 12
    heightanotC = mainy + heightmain + shifttrialtype;
    anot7pos = [pos1x+width/2-shiftanot7x, heightanotC, 6.5720 0.4286];
    heightanotN = mainy + heightmain + shifttrialtype ;
    anot8pos = [pos2x+width/2-shiftanot8x, heightanotN,6.5720 0.4286];
    
    % anotations satitation
    anot2x = pos1x - shiftAllx - shiftall; anot2y = sub2y + shiftAlly;
    anot3x = pos1x - shiftSelx - shiftall; anot3y = sub3y + shiftSely;
    cbhx   =  pos2x+width + shfitcbhx;  cbhy   =  sub3y ;
    anot4x =  pos2x+width +shiftpvalx;  anot4y =  sub2y + shiftpvaly;
    
    cbhpos = [cbhx, cbhy,  0.2438    1.6002];
    anot2pos = [anot2x, anot2y,  9.2234    0.4889]  ; %'All'
    anot3pos = [anot3x, anot3y,  9.2234    0.4889]  ; %'Seected';
    anot4pos = [anot4x, anot4y,  0.2, 0.2]  ; %'pval'
    
    % legend
    legx =  pos2x+width - shiftlegx ; legy = mainy + shiftlegy;
    legpos =[legx, legy, sizelegx    1.4817];
    
    
    if figures == 1
        do_the_plot(fig1, {'ALL'}, autocorr_fold, fold, clustPermFold,  clustPermFold_auto, DirectChangeFold, DirectChangeFold_auto, {'Rankabsdiff'}, {'Abs Diff (Rk)'}, pos11, pos12,pos13, pos21, pos22, pos23, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot7pos, anot8pos , legpos, possible_colors([4,3], :), 'B', 1, 1, 0, 1, 0, 0)
    elseif figures == 2
        do_the_plot(fig1, {'ALL'}, autocorr_fold, fold, clustPermFold, clustPermFold_auto, DirectChangeFold, DirectChangeFold_auto, {'absdiff'}, {'Abs Diff (Rt)'}, pos11, pos12,pos13, pos21, pos22, pos23, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot7pos, anot8pos , legpos, possible_colors([3,3], :), 'A', 0, 1, 1, 1, 0, 0)
    end
    
    accumulator = accumulator + totalheight;
end



toSave = [folder2save, '/Clusterwise Permutation test, Abs. Difference'];
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






%% Figure 4 in paper
accumulator = 0;
pos1x = pos1xor;
pos2x = pos1x + diffx1x2;


figpos = [0,0, 30.4800+shiftall, 25.4];
close all
fig1 =  figure('color', 'w', 'units', 'centimeters', 'position', figpos, 'visible', 'off');

mainystart = 3.5;
for  figures = 1:2
    
    mainy = mainystart + accumulator;
    sub2y=  mainy - sub2yadd;
    sub3y=  mainy - sub3yadd;
    
    if  figures == 1
        sub2y=  mainy - sub2yadd - 0.2;
        sub3y=  mainy - sub3yadd - 0.2;
    end
    
    pos21    =  [pos2x    mainy  width    heightmain];
    pos22    =  [pos2x    sub2y  width    heightTicks];
    pos23    =  [pos2x    sub3y  width    heightTicks];
    pos11    =  [pos1x    mainy  width    heightmain];
    pos12    =  [pos1x    sub2y  width    heightTicks];
    pos13    =  [pos1x    sub3y  width    heightTicks];
    
    
    
    % letter for figure
    heightanotA = mainy + heightmain + shiftletra;
    anot5pos = [pos1x-shiftletrax, heightanotA, 4.5720 0.4286];
    
    
    % current next trial, size 12
    heightanotC = mainy + heightmain + shifttrialtype;
    anot7pos = [pos1x+width/2-shiftanot7x, heightanotC, 6.5720 0.4286];
    heightanotN = mainy + heightmain + shifttrialtype ;
    anot8pos = [pos2x+width/2-shiftanot8x, heightanotN,6.5720 0.4286];
    
    % anotations satitation
    anot2x = pos1x - shiftAllx - shiftall; anot2y = sub2y + shiftAlly;
    anot3x = pos1x - shiftSelx - shiftall; anot3y = sub3y + shiftSely;
    cbhx   =  pos2x+width + shfitcbhx;  cbhy   =  sub3y ;
    anot4x =  pos2x+width +shiftpvalx;  anot4y =  sub2y + shiftpvaly;
    
    cbhpos = [cbhx, cbhy,  0.2438    1.6002];
    anot2pos = [anot2x, anot2y,  9.2234    0.4889]  ; %'All'
    anot3pos = [anot3x, anot3y,  9.2234    0.4889]  ; %'Seected';
    anot4pos = [anot4x, anot4y,  0.2, 0.2]  ; %'pval'
    
    % legend
    legx =  pos2x+width - shiftlegx ; legy = mainy + shiftlegy;
    legpos =[legx, legy, sizelegx    1.4817];
    
    
    if figures == 1
        do_the_plot(fig1, {'ALL'}, autocorr_fold, fold, clustPermFold,  clustPermFold_auto, DirectChangeFold, DirectChangeFold_auto, {'RT'}, {'RT'}, pos11, pos12,pos13, pos21, pos22, pos23, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot7pos, anot8pos , legpos, possible_colors([2,3], :), 'B', 1, 1, 0, 1, 0, 0)
    elseif figures == 2
        do_the_plot(fig1, {'ALL'}, autocorr_fold, fold, clustPermFold, clustPermFold_auto, DirectChangeFold, DirectChangeFold_auto, {'Abspref'}, {'Confidence'}, pos11, pos12,pos13, pos21, pos22, pos23, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot7pos, anot8pos , legpos, possible_colors([1,3], :), 'A', 0, 1, 1, 1, 0, 0)
    end
    
    accumulator = accumulator + totalheight;
end



toSave = [folder2save, '/Clusterwise Permutation test, confidence and RT'];
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







%% Supp. Figure S3, RT and salience in the Valuation Task


accumulator = 0;accumulator2 = 0 ;
pos1x = pos1xor;
pos2x = pos1x + diffx1x2;
figpos = [0,0, 30.4800+shiftall, 25];
close all
fig1 =  figure('color', 'w', 'units', 'centimeters', 'position', figpos, 'visible', 'off');



mainystart = 3.5;

for  figures = 1:2
    
    % mainplots  and satiation
    
    mainy = mainystart + accumulator;
    sub2y=  mainy - sub2yadd;
    sub3y=  mainy - sub3yadd;
    if  figures == 1
        sub2y=  mainy - sub2yadd - 0.2;
        sub3y=  mainy - sub3yadd - 0.2;
    end
    
    
    pos21    =  [pos2x    mainy  width    heightmain];
    pos22    =  [pos2x    sub2y  width    heightTicks];
    pos23    =  [pos2x    sub3y  width    heightTicks];
    pos11    =  [pos1x    mainy  width    heightmain];
    pos12    =  [pos1x    sub2y  width    heightTicks];
    pos13    =  [pos1x    sub3y  width    heightTicks];
    

    
    % letter for figure
    heightanotA = mainy + heightmain + shiftletra;
    anot5pos = [pos1x-shiftletrax, heightanotA, 4.5720 0.4286];
    
    %heightanotB = posvshapey + heightVshape +1.4;
    %anot6pos = [pos1x-1, heightanotB,4.5720 0.4286];
    
    
    % current next trial, size 12
    heightanotC = mainy + heightmain + shifttrialtype;
    anot7pos = [pos1x+width/2-shiftanot7x, heightanotC, 6.5720 0.4286];
    heightanotN = mainy + heightmain + shifttrialtype;
    anot8pos = [pos2x+width/2-shiftanot8x, heightanotN,6.5720 0.4286];
    
    % anotations satitation
    anot2x = pos1x - shiftAllx - shiftall; anot2y = sub2y + shiftAlly;
    anot3x = pos1x - shiftSelx - shiftall; anot3y = sub3y + shiftSely;
    cbhx   =  pos2x+width + shfitcbhx;  cbhy   =  sub3y ;
    anot4x =  pos2x+width +shiftpvalx;  anot4y =  sub2y + shiftpvaly;
    
    cbhpos = [cbhx, cbhy,  0.2438    1.6002];
    anot2pos = [anot2x, anot2y,  9.2234    0.4889]  ; %'All'
    anot3pos = [anot3x, anot3y,  9.2234    0.4889]  ; %'Seected';
    anot4pos = [anot4x, anot4y,  0.2, 0.2]  ; %'pval'
    
    % legend
    legx =  pos2x+width - shiftlegx ; legy = mainy + shiftlegy;
    legpos =[legx, legy, sizelegx    1.4817];
    
   
    if figures == 1
        do_the_plot2(fig1, {'ALL'},autocorr_fold, fold,  clustPermFold2, clustPermFold2_auto, DirectChangeFold, DirectChangeFold_auto, {'RT_R'}, {'RT'}, pos11, pos12,pos13, pos21, pos22, pos23, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot7pos, anot8pos , legpos, possible_colors(2, :), 'B', 1, 1, 0, 1)
        accumulator2 = 0 ;
    elseif figures == 2
        do_the_plot2(fig1, {'ALL'}, autocorr_fold, fold, clustPermFold2, clustPermFold2_auto, DirectChangeFold, DirectChangeFold_auto,  {'Sal'}, {'Confidence'}, pos11, pos12,pos13, pos21, pos22, pos23, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot7pos, anot8pos , legpos, possible_colors(1, :), 'A', 0, 1, 1, 1)
        accumulator2 = 0 ;
        accumulator2 = - (heightTicks + heightTicks  + diff32 );
    end
    
    
    accumulator = accumulator + totalheight + accumulator2;
end





toSave = [folder2save, '/Clusterwise Permutation test, Supplementary Sal RT Valuation Task'];
print([toSave, '200'], '-dpng', '-r200');

print([toSave, '300'], '-dpng', '-r300');
print(toSave, '-dpng', '-r500');

print([toSave, 'Vec'], '-depsc', '-painters');
print([toSave, 'Vec'], '-depsc', '-opengl', '-r500');
set(fig1,'Units','centimeters');
pos = get(fig1,'Position');
set(fig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

print([toSave, 'Vec'], '-dpdf', '-painters');
print(toSave, '-dpdf', '-r500', '-opengl');
print([toSave, ' 300'], '-dpdf', '-r300', '-opengl');







%% Plots direct change of firinge rate for reviwers, not in paper

% Alpha values from tests
alpha = 0.05;       % Alpha from correlation over trials
alphaSigned = 0.05; % Alpha from signed rank test
alphaBinom = 0.05;  % Alpha from binomial test across neurons
alphaPerm = 0.01;   % Alpha from clusterwise permutation across time

stralphaBinom= num2str(alphaBinom);
type_clusterSumEffect = ['clusterSumEffect', stralphaBinom(3:end)];
typeclusterSumEffect_DirectChange  =  ['clusterSumEffect_DirectChange', stralphaBinom(3:end)];

accumulator = 0;
pos1x = pos1xor;
pos2x = pos1x + diffx1x2;

totalheight2 = heightmain +sub2yadd + 1 +.2;

figpos = [0,0, 30.4800+shiftall, 42.5];
close all
fig1 =  figure('color', 'w', 'units', 'centimeters', 'position', figpos, 'visible', 'off');

mainystart =  1.4290 ;
for  figures = 1:4
    % mainplots  and satiation
    mainy = mainystart + accumulator;
    sub2y=  mainy - sub2yadd;
    sub3y=  mainy - sub3yadd;
    
    if  figures == 1
        sub2y=  mainy - sub2yadd - 0.2;
        sub3y=  mainy - sub3yadd - 0.2;
    end
    
    pos21    =  [pos2x    mainy  width    heightmain];
    pos22    =  [pos2x    sub2y  width    heightTicks];
    pos23    =  [pos2x    sub3y  width    heightTicks];
    pos11    =  [pos1x    mainy  width    heightmain];
    pos12    =  [pos1x    sub2y  width    heightTicks];
    pos13    =  [pos1x    sub3y  width    heightTicks];
    
   
    
    % letter for figure
    heightanotA = mainy + heightmain + shiftletra;
    anot5pos = [pos1x-shiftletrax, heightanotA, 4.5720 0.4286];
    
    %heightanotB = posvshapey + heightVshape +1.4;
    %anot6pos = [pos1x-1, heightanotB,4.5720 0.4286];
    
    
    % current next trial, size 12
    heightanotC = mainy + heightmain + shifttrialtype;
    anot7pos = [pos1x+width/2-shiftanot7x, heightanotC, 6.5720 0.4286];
    heightanotN = mainy + heightmain + shifttrialtype ;
    anot8pos = [pos2x+width/2-shiftanot8x, heightanotN,6.5720 0.4286];
    
    % anotations satitation
    anot2x = pos1x - shiftAllx - shiftall; anot2y = sub2y + shiftAlly;
    anot3x = pos1x - shiftSelx - shiftall; anot3y = sub3y + shiftSely;
    cbhx   =  pos2x+width + shfitcbhx;  cbhy   =  sub3y ;
    anot4x =  pos2x+width +shiftpvalx;  anot4y =  sub2y + shiftpvaly;
    
    cbhpos = [cbhx, cbhy,  0.2438    1.6002];
    anot2pos = [anot2x, anot2y,  9.2234    0.4889]  ; %'All'
    anot3pos = [anot3x, anot3y,  9.2234    0.4889]  ; %'Seected';
    anot4pos = [anot4x, anot4y,  0.2, 0.2]  ; %'pval'
    
    % legend
    legx =  pos2x+width - shiftlegx ; legy = mainy + shiftlegy;
    legpos =[legx, legy, sizelegx    1.4817];
    
    
    if figures == 1
        do_the_plot(fig1, {'ALL'}, autocorr_fold, fold, clustPermFold,  clustPermFold_auto, DirectChangeFold, DirectChangeFold_auto, {'Rankabsdiff'}, {'Abs Diff (Rk)'}, pos11, pos12,pos13, pos21, pos22, pos23, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot7pos, anot8pos , legpos, possible_colors([4,4], :), 'D', 1, 1, 0, 1, 0, 1)
    elseif figures == 2
        do_the_plot(fig1, {'ALL'}, autocorr_fold, fold, clustPermFold, clustPermFold_auto, DirectChangeFold, DirectChangeFold_auto, {'absdiff'}, {'Abs Diff (Rt)'}, pos11, pos12,pos13, pos21, pos22, pos23, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot7pos, anot8pos , legpos, possible_colors([3,3], :), 'C', 0, 1, 0, 1, 0, 1)
    elseif figures == 3
        do_the_plot(fig1, {'ALL'}, autocorr_fold, fold, clustPermFold,  clustPermFold_auto, DirectChangeFold, DirectChangeFold_auto, {'RT'}, {'RT'}, pos11, pos12,pos13, pos21, pos22, pos23, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot7pos, anot8pos , legpos, possible_colors([2,2], :), 'B', 0, 1, 0, 1, 0, 1)
        %         do_the_plot(fig1, {'ALL'}, autocorr_fold, fold, clustPermFold,  clustPermFold_auto, DirectChangeFold, DirectChangeFold_auto,{'Abspref'}, {'Confidence'}, pos11, pos12,pos13, pos21, pos22, pos23, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot7pos, anot8pos , legpos, possible_colors([1,3], :), 'A', 1, 1, 1, 1)
    elseif figures == 4
        do_the_plot(fig1, {'ALL'}, autocorr_fold, fold, clustPermFold, clustPermFold_auto, DirectChangeFold, DirectChangeFold_auto, {'Abspref'}, {'Confidence'}, pos11, pos12,pos13, pos21, pos22, pos23, cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot7pos, anot8pos , legpos, possible_colors([1,1], :), 'A', 0, 1, 1, 1, 0, 1)
    end
    
    accumulator = accumulator + totalheight2;
end



toSave = [folder2save, '/Clusterwise Permutation test, Direct Change Confidence, RT, Abs. Difference , alphaBinom ', num2str(alphaBinom), ', alphaPerm', num2str(alphaPerm)];
print([toSave, '200'], '-dpng', '-r200');









function do_the_plot(fig1, region2, fold_auto, fold, clustPermFold, clustPermFold_auto, DirectChangeFold, DirectChangeFold_auto, var2use, var2name, pos11, pos12,pos13, pos21, pos22, pos23,...
    cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot7pos, anot8pos , legpos, color2use, letra, xlabelon, legon, plottrialtype, regtype, plotChangeProportion, onlyPlotChange)

global alpha alphaBinom alphaPerm alphaSigned type_clusterSumEffect typeclusterSumEffect_DirectChange
set(0, 'currentFigure',fig1)


dbstop if error
%plot_paper_inverses_autocorr_reduced(0.05)



set(groot,'defaultAxesLinewidth', 1);
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


listRatingCombine = { 'Val', 'Sal', 'RT_R'};

listRating = {'ValTot', 'SalTot', 'Val', 'Sal', 'Valmean', 'Salmean', 'Rank', 'RTmean', 'RTstd', 'Valstd', 'RT_R', 'Sign', 'SignTot',   'Negsplit', 'Possplit', 'NegsplitTot', 'PossplitTot',  'ValSweetsplit', 'SalSweetsplit', 'ValSweetsplitTot', 'SalSweetsplitTot', 'RankSaltysplit', ...
    'ValSweetsplit', 'SalSweetsplit', 'ValSweetsplitTot', 'SalSweetsplitTot', 'RankSweetsplit'};


% for do_single = 0 % 0:2; %do_single==2 is is a single unit, do_single==1 is is a multi unit only, do_single==0 is is a  all SU and MU
%
%     for restricted = 0%:2

for g=1%:size(region2,1)
    
    region = region2{g}; load(['../Rasts matfiles/patientVector_', region]);
    vectorPat = eval(['patientVector_', region]); % extract patient ID of each cell
    vectorPat = ismember([vectorPat; vectorPat], [2,3,4,5,6,1,7,8,9,10,11,12]);
    
    
    
    for cc= 1%: size(var2use,2)
        
        variable = var2use{cc};
        disp([variable, '  ',  region])
        
        
        % change window of interest for valauation task
        if ismember(variable, listRating)
            leftReduction = (2000+2350)/windowstep-1;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr
            rightReduction =  trialonset-(2350)/windowstep;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr
        end
        
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
                toLoad= {['*', variable, 'inv_A1']; ['*', variable, '_A1']; [variable, '_', type_clusterSumEffect]; [variable, '_', typeclusterSumEffect_DirectChange]; [variable, '_signedRankclusterSumEffect']; [variable, '_signedRankSelclusterSumEffect']};
            end
            load([fold,'/regresspval_',region,'.mat'],toLoad{1:2});
            load([fold,'/regressreg_', region,'.mat'],toLoad{1:2});
            
            invvar =  eval(['pval', toLoad{1}(2:end)]);
            var=      eval(['pval', toLoad{2}(2:end)]);
            cinvvar=  eval(['corr', toLoad{1}(2:end)]);
            cvar=     eval(['corr', toLoad{2}(2:end)]);
            
            
            load([clustPermFold,'/permdistribution_',region,'2.mat'],toLoad{3:end});
            
            clustVar        =  eval([toLoad{3}]);
            clustVarDC      =  eval([toLoad{4}]);
            clustVarSR_All  =  eval([toLoad{5}]);
            clustVarSR_Sel  =  eval([toLoad{6}]);
            
            
            toLoad2 = {['*', variable,'_RawChange_A1']};
            load([DirectChangeFold, 'regresspval_',region,'.mat'],toLoad2{1});
            load([DirectChangeFold, 'regressreg_',region,'.mat'],toLoad2{1});
            
            varchange  =   eval(['pval', toLoad2{1}(2:end)]);
            %             cvarchange =   eval(['corr', toLoad2{1}(2:end)]);
            
        elseif regtype == 2
            toLoad= {['*', variable, 'reginv_A1']; ['*', variable, 'reg_A1']; ['*', variable, 'Autoreginv_A1']; ['*', variable, 'Autoreg_A1']};
            
            load([fold,'/regresspval_',region,'.mat'],toLoad{:});
            load([fold,'/regressreg_', region,'.mat'],toLoad{:});
            invvar =  eval(['pval', toLoad{1}(2:end)]);
            var=      eval(['pval', toLoad{2}(2:end)]);
            cinvvar=  eval(['corr', toLoad{1}(2:end)]);
            cvar=     eval(['corr', toLoad{2}(2:end)]);
            
            
            invvara =  eval(['pval', toLoad{3}(2:end)]);
            vara=      eval(['pval', toLoad{4}(2:end)]);
            cinvvara=  eval(['corr', toLoad{3}(2:end)]);
            cvara=     eval(['corr', toLoad{4}(2:end)]);
            
        end
        
        
        if ismember(variable, listRatingCombine)
            cvar =  eval(['corr', toLoad{3}(2:end)]);
        end
        
        
        %% define the limits of the plot
        x_ax_lim = [0 + rightReduction, trialonset+  leftReduction];
        
        
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
        
        meancvar = abs(cvar(1:siVar,:) + cvar(1+siVar:end ,:))/2;
        meancinvvar = abs(cinvvar(1:siVar,:) + cinvvar(1+siVar:end ,:))/2;
        
        meancvarSel = abs(cvarSel(1:siVar,:) + cvarSel(1+siVar:end ,:))/2;
        meancinvvarSel = abs(cinvvarSel(1:siVar,:) + cinvvarSel(1+siVar:end ,:))/2;
        
        
        
        
        
        %% ticks
        
        Xtickspos=[[trialonset:-tickpossize: 1], [trialonset:tickpossize:size(cvar,2)]];  Xtickspos = sort(Xtickspos(2:end)); % from the trial aligment, go forward and backwards determining te position and sort
        Xtickslabelsnames = [[0:-1:-(length([trialonset:-tickpossize: 1])-1)], [0:1:length(trialonset:tickpossize:size(cvar,2))-1]]; Xtickslabelsnames = sort(Xtickslabelsnames(2:end));
        
        for hh=1:length(Xtickslabelsnames)
            aa{hh}=num2str(Xtickslabelsnames(hh));
        end
        Xtickslabelsnames=[];Xtickslabelsnames=aa;
        
        
        
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
        line(x_ax_lim, [alpha,alpha], 'linestyle', '--','color','k', 'linewidth', 1.3)
        line(x_ax_lim, [pos,pos], 'linestyle', '--', 'color','r', 'linewidth', 1.3)
        ylabel('Proportion of significant cells')
        
        
        
        
        % plot Areas
        if strcmp(toLoad{1}(end-1:end), '_R')
            area(AreaRat, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            if xlabelon ==1
                xlabel('Time from stimulus onset (s)')
            end
            xlim(x_ax_lim)
        else
            
            area(AreaA1, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            area(AreaB, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.3,0.3,0.3], 'Edgecolor', 'none')
            area(AreaA2, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            if xlabelon ==1
                xlabel('Time relative to stimulus A1 onset (s)')
            end
            xlim(x_ax_lim)
        end
        
        
        
        % direct change plot
        
        
        if onlyPlotChange ||  plotChangeProportion
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
            
            range2use = trialonset+3600/windowstep:trialonset+((3600+2350-1)/windowstep);
            num121 =  mean(mean(toplot2(:, range2use)));
            p121 =  myBinomTest(num121*siVar, siVar, alpha, 'one');
            disp (['total fraction of units for direct change for variable  ', variable, ':  ', num2str(num121), '  Pvalue binomial test: ', num2str(p121) ])
        end
        %binomial line
        line(x_ax_lim, [alpha,alpha], 'linestyle', '--','color','k', 'linewidth', 1.3)
        line(x_ax_lim, [pos,pos], 'linestyle', '--', 'color','r', 'linewidth', 1.3)
        
        % settings
        ylim([0, heightplot])
        set(gca, 'Xtick',Xtickspos, 'Xticklabels',Xtickslabelsnames);
        set(gca, 'FontSize',12, 'fontweight', 'bold')
        
        
        
        
        
        
        
        
        %% pval plot
        
        if (onlyPlotChange ~= 1)
            
            xaxis= 1:size(var,2);
            
            clear pvalSign zvalues
            for kk= 1:size(cvar,2)
                try
                    [pvalSign(1, kk), ~, stats ] = signrank(meancvar(:,kk), meancinvvar(:,kk), 'tail', 'right');
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
                    [pvalSignSel(1, kk), ~, stats] = signrank(meancvarSel(:,kk), meancinvvarSel(:,kk), 'tail', 'right');
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
        
        
        
        
        %% autocorr plot
        
        
        
        
        %%
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
                toLoad= {['*', variable, 'inv_A1']; ['*', variable, '_A1']; [variable, '_', type_clusterSumEffect]; [variable, '_', typeclusterSumEffect_DirectChange]; [variable, '_signedRankclusterSumEffect']; [variable, '_signedRankSelclusterSumEffect']};
            end
            load([fold_auto,'/regresspval_',region,'.mat'],toLoad{1:2});
            load([fold_auto,'/regressreg_', region,'.mat'],toLoad{1:2});
            
            invvar =  eval(['pval', toLoad{1}(2:end)]);
            var=      eval(['pval', toLoad{2}(2:end)]);
            cinvvar=  eval(['corr', toLoad{1}(2:end)]);
            cvar=     eval(['corr', toLoad{2}(2:end)]);
            
            
            load([clustPermFold_auto,'/permdistribution_',region,'2.mat'],toLoad{3:end});
            
            clustVar        =  eval([toLoad{3}]);
            clustVarDC      =  eval([toLoad{4}]);
            clustVarSR_All  =  eval([toLoad{5}]);
            clustVarSR_Sel  =  eval([toLoad{6}]);
            
            toLoad2 = {['*', variable,'_RawChange_A1']};
            
            load([DirectChangeFold_auto, 'regresspval_',region,'.mat'],toLoad2{1});
            load([DirectChangeFold_auto, 'regressreg_',region,'.mat'],toLoad2{1});
            
            varchange  =   eval(['pval', toLoad2{1}(2:end)]);
            %             cvarchange =   eval(['corr', toLoad2{1}(2:end)]);
        elseif regtype == 2
            toLoad= {['*', variable, 'reginv_A1']; ['*', variable, 'reg_A1']; ['*', variable, 'Autoreginv_A1']; ['*', variable, 'Autoreg_A1']};
            
            load([fold,'/regresspval_',region,'.mat'],toLoad{:});
            load([fold,'/regressreg_', region,'.mat'],toLoad{:});
            invvar =  eval(['pval', toLoad{1}(2:end)]);
            var=      eval(['pval', toLoad{2}(2:end)]);
            cinvvar=  eval(['corr', toLoad{1}(2:end)]);
            cvar=     eval(['corr', toLoad{2}(2:end)]);
            
            
            invvara =  eval(['pval', toLoad{3}(2:end)]);
            vara=      eval(['pval', toLoad{4}(2:end)]);
            cinvvara=  eval(['corr', toLoad{3}(2:end)]);
            cvara=     eval(['corr', toLoad{4}(2:end)]);
            
        end
        
        if ismember(variable, listRatingCombine)
            cvar =  eval(['corr', toLoad{3}(2:end)]);
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
        
        meancvar = abs(cvar(1:siVar,:) + cvar(1+siVar:end ,:))/2;
        meancinvvar = abs(cinvvar(1:siVar,:) + cinvvar(1+siVar:end ,:))/2;
        
        meancvarSel = abs(cvarSel(1:siVar,:) + cvarSel(1+siVar:end ,:))/2;
        meancinvvarSel = abs(cinvvarSel(1:siVar,:) + cinvvarSel(1+siVar:end ,:))/2;
        
        
        
        
        
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
        line(x_ax_lim, [alpha,alpha], 'linestyle', '--','color','k', 'linewidth', 1.3)
        line(x_ax_lim, [pos,pos], 'linestyle', '--', 'color','r', 'linewidth', 1.3)
        %ylabel('Proportion of significant cells')
        
        
        
        if onlyPlotChange ||  plotChangeProportion
            % direct change plot
            toplotDC= varchange < alpha;
            [toplotDC2] = create_permVector (toplotDC, clustVarDC, pos, x_ax_lim, alphaPerm);
            plot(toplotDC2, 'color', [color2use(2,:)], 'linewidth', 1.5)
        end
        if onlyPlotChange ~= 1
            %  proportion plot
            toplot2= ((var(1:siVar,:)<sqrt(2*alpha)) & (var(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvar(1:siVar,:))==sign(cvar(1+siVar:end,:))));
            [toplot] = create_permVector (toplot2, clustVar, pos, x_ax_lim, alphaPerm);
            plot(toplot, 'color', [color2use(1,:)], 'linewidth', 1.5)
            range2use = trialonset+3600/windowstep:trialonset+((3600+2350-1)/windowstep);
            num121 =  mean(mean(toplot2(:, range2use)));
            p121 =  myBinomTest(num121*siVar, siVar, alpha, 'one');
            disp (['total fraction of units for variable  ', variable, ':  ', num2str(num121), '  Pvalue binomial test: ', num2str(p121)])
        end
        
        % plot Areas
        if strcmp(toLoad{1}(end-1:end), '_R')
            area(AreaRat, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            if xlabelon ==1
                xlabel('Time from stimulus onset (s)')
            end
            xlim(x_ax_lim)
        else
            
            area(AreaA1, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            area(AreaB, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.3,0.3,0.3], 'Edgecolor', 'none')
            area(AreaA2, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            if xlabelon ==1
                xlabel('Time relative to stimulus A1 onset (s)')
            end
            xlim(x_ax_lim)
        end
        
        
        if onlyPlotChange ||  plotChangeProportion
            % direct change plot
            toplotDC= varchange < alpha;
            [toplotDC2] = create_permVector (toplotDC, clustVarDC, pos, x_ax_lim, alphaPerm);
            plot(mean(toplotDC), 'color', [color2use(2,:), transparency], 'linewidth',1.5)
            plot(toplotDC2, 'color', [color2use(2,:)], 'linewidth', 1.5)
        end
        if (onlyPlotChange ~= 1)
            %  proportion plot
            toplot2= ((var(1:siVar,:)<sqrt(2*alpha)) & (var(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvar(1:siVar,:))==sign(cvar(1+siVar:end,:))));
            [toplot] = create_permVector (toplot2, clustVar, pos, x_ax_lim, alphaPerm);
            plot(sum(toplot2)./siVar, 'color', [color2use(1,:), transparency], 'linewidth',1.5)
            plot(toplot, 'color', [color2use(1,:)], 'linewidth', 1.5)
        end
        
        %binomial line
        line(x_ax_lim, [alpha,alpha], 'linestyle', '--','color','k', 'linewidth', 1.3)
        line(x_ax_lim, [pos,pos], 'linestyle', '--', 'color','r', 'linewidth', 1.3)
        
        % settings
        ylim([0, heightplot])
        set(gca, 'Xtick',Xtickspos, 'Xticklabels',Xtickslabelsnames);
        set(gca, 'FontSize',12, 'fontweight', 'bold')
        
        
        if legon ==1
            if (plotChangeProportion == 1) && onlyPlotChange == 0
                leg =legend('Chance', 'Significance',  'Change', var2name{cc}); % [var2name{cc},'Change']
            else
                leg =legend('Chance', 'Significance',  var2name{cc});
            end
            set(leg, 'units','centimeters', 'Position', legpos ,'FontSize',12);
        end
        
        
        
        
        %% pval plot
        if (onlyPlotChange ~= 1)
            xaxis= 1:size(var,2);
            
            clear pvalSign zvalues
            for kk= 1:size(cvar,2)
                try
                    [pvalSign(1, kk), ~, stats ] = signrank(meancvar(:,kk), meancinvvar(:,kk), 'tail', 'right');
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
            
            %%  clusterwise correction
            
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
                    [pvalSignSel(1, kk), ~, stats] = signrank(meancvarSel(:,kk), meancinvvarSel(:,kk), 'tail', 'right');
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
            %sim1DC2 = corr( sum(toplot2(:,xlimcorr))' , sum(toplotDC(:, xlimcorr))', 'type', 'Spearman', 'rows', 'c');
            
            disp(['mean correlation over time between fraction of significant neurons and Satiation according to All SignedRank =  ' , num2str((sim1 + sim2)/2)])
            disp(['mean correlation over time between fraction of significant neurons and Satiation according to Selected SignedRank =  ' , num2str((simSel1 + simSel2)/2)])
            %          disp(['mean correlation over time between fraction of significant neurons and Satiation according to Direct Change =  ' , num2str((sim1DC1 + sim1DC2)/2)])
            
            
            cbh = colorbar;
            cbh.FontSize= 11;
            cbh.Units = 'centimeters';
            cbh.Position= cbhpos;
            cbh.FontWeight= 'bold';
            anot4 = annotation('textbox', [0,0,0.1,0.1],   'string',['pval'],'Units','normalized','FontName','Helvetica','FontSize',11,'EdgeColor','none','fontweight', 'bold');
            
            
            anot2 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',['All'],'Units','normalized','FontName','Helvetica','FontSize',12,'EdgeColor','none', 'fontweight', 'bold');
            anot3 = annotation('textbox',[0.00003    0+0.02   0.3    0.04],   'string',['Selected'],'Units','normalized','FontName','Helvetica','FontSize',12,'EdgeColor','none', 'fontweight', 'bold');
            set(anot2, 'units', 'centimeters', 'position', anot2pos)
            set(anot3, 'units', 'centimeters', 'position', anot3pos)
            set(anot4, 'units', 'centimeters', 'position', anot4pos)
        end
        
        global sizeletra sizetrialtype
        anot5 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',[letra],'Units','normalized','FontName','Helvetica','FontSize',sizeletra,'EdgeColor','none', 'fontweight', 'bold');
        
        if plottrialtype ==1
            anot7 = annotation('textbox',[0.00003    0+0.02   0.3    0.04],   'string',['Current Trial'],'Units','normalized','FontName','Helvetica','FontSize',sizetrialtype,'EdgeColor','none', 'fontweight', 'bold');
            anot8 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',['Next Trial'],'Units','normalized','FontName','Helvetica','FontSize',sizetrialtype,'EdgeColor','none', 'fontweight', 'bold');
            set(anot7, 'units', 'centimeters', 'position', anot7pos)
            set(anot8, 'units', 'centimeters', 'position', anot8pos)
        end
        
        
        set(anot5, 'units', 'centimeters', 'position', anot5pos)
        
        
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













function do_the_plot2(fig1, region2, fold_auto, fold, clustPermFold, clustPermFold_auto, DirectChangeFold, DirectChangeFold_auto, var2use, var2name, pos11, pos12,pos13, pos21, pos22, pos23,...
    cbhpos, anot2pos, anot3pos, anot4pos, anot5pos , anot7pos, anot8pos , legpos, color2use, letra, xlabelon, legon, plottrialtype, regtype)

global alpha alphaBinom alphaPerm alphaSigned type_clusterSumEffect typeclusterSumEffect_DirectChange
set(0, 'currentFigure',fig1)
dbstop if error
set(groot,'defaultAxesLinewidth', 1);
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


listRatingCombine = { 'Val', 'Sal', 'RT_R'};

listRating = {'ValTot', 'SalTot', 'Val', 'Sal', 'Valmean', 'Salmean', 'Rank', 'RTmean', 'RTstd', 'Valstd', 'RT_R', 'Sign', 'SignTot',   'Negsplit', 'Possplit', 'NegsplitTot', 'PossplitTot',  'ValSweetsplit', 'SalSweetsplit', 'ValSweetsplitTot', 'SalSweetsplitTot', 'RankSaltysplit', ...
    'ValSweetsplit', 'SalSweetsplit', 'ValSweetsplitTot', 'SalSweetsplitTot', 'RankSweetsplit'};


% for do_single = 0 % 0:2; %do_single==2 is is a single unit, do_single==1 is is a multi unit only, do_single==0 is is a  all SU and MU
%
%     for restricted = 0%:2

for g=1%:size(region2,1)
    
    region = region2{g}; load(['../Rasts matfiles/patientVector_', region]);
    vectorPat = eval(['patientVector_', region]); % extract patient ID of each cell
    vectorPat = ismember([vectorPat; vectorPat], [2,3,4,5,6,1,7,8,9,10,11,12]);
    
    
    for cc= 1%: size(var2use,2)
        
        variable = var2use{cc};
        disp([variable, '  ',  region])
        
        
        % change window of interest for valauation task
        if ismember(variable, listRating)
            leftReduction = (2000+2350)/windowstep-1;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr
            rightReduction =  trialonset-(2350)/windowstep;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr
        end
        
        
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
                toLoad= {['*', variable, 'inv_A1']; ['*', variable, '_A1']; [variable, '_', type_clusterSumEffect]; [variable, '_', typeclusterSumEffect_DirectChange]; [variable, '_signedRankclusterSumEffect']; [variable, '_signedRankSelclusterSumEffect']};
            end
            load([fold,'/regresspval_',region,'.mat'],toLoad{1:2});
            load([fold,'/regressreg_', region,'.mat'],toLoad{1:2});
            
            invvar =  eval(['pval', toLoad{1}(2:end)]);
            var=      eval(['pval', toLoad{2}(2:end)]);
            cinvvar=  eval(['corr', toLoad{1}(2:end)]);
            cvar=     eval(['corr', toLoad{2}(2:end)]);
            
            if ismember(variable, listRatingCombine)
                load([fold,'/regresspval_',region,'.mat'],toLoad{3});
                load([fold,'/regressreg_', region,'.mat'],toLoad{3});
                varmean=      eval(['pval', toLoad{3}(2:end)]);
                cvarmean=     eval(['corr', toLoad{3}(2:end)]);
            end
            
            if strcmp(variable, 'RT_R')
                toLoad2 = {['RT_Rat_', type_clusterSumEffect]; ['RT_Rat_', typeclusterSumEffect_DirectChange]; ['RT_Rat_signedRankclusterSumEffect']; ['RT_Rat_signedRankSelclusterSumEffect']};
            else
                toLoad2 = {[variable, '_', type_clusterSumEffect]; [variable, '_', typeclusterSumEffect_DirectChange]; [variable, '_signedRankclusterSumEffect']; [variable, '_signedRankSelclusterSumEffect']};
            end
            
            load([clustPermFold,'/permdistribution_',region,'2.mat'],toLoad2{1:end});
            
            clustVar        =  eval([toLoad2{1}]);
            clustVarDC      =  eval([toLoad2{2}]);
            clustVarSR_All  =  eval([toLoad2{3}]);
            clustVarSR_Sel  =  eval([toLoad2{4}]);
            
            
            
            %             if strcmp(variable, 'RT_R')
            %                     toLoad3 = {['*RT_RawChange_A1']};
            %             else
            %                 toLoad3 = {['*', variable,'_RawChange_A1']};
            %             end
            %
            %             load([DirectChangeFold, 'regresspval_',region,'.mat'],toLoad3{1});
            %             load([DirectChangeFold, 'regressreg_',region,'.mat'],toLoad3{1});
            %
            %             varchange  =   eval(['pval', toLoad3{1}(2:end)]);
            %             cvarchange =   eval(['corr', toLoad3{1}(2:end)]);
            
        elseif regtype == 2
            toLoad= {['*', variable, 'reginv_A1']; ['*', variable, 'reg_A1']; ['*', variable, 'Autoreginv_A1']; ['*', variable, 'Autoreg_A1']};
            
            load([fold,'/regresspval_',region,'.mat'],toLoad{:});
            load([fold,'/regressreg_', region,'.mat'],toLoad{:});
            invvar =  eval(['pval', toLoad{1}(2:end)]);
            var=      eval(['pval', toLoad{2}(2:end)]);
            cinvvar=  eval(['corr', toLoad{1}(2:end)]);
            cvar=     eval(['corr', toLoad{2}(2:end)]);
            
            
            invvara =  eval(['pval', toLoad{3}(2:end)]);
            vara=      eval(['pval', toLoad{4}(2:end)]);
            cinvvara=  eval(['corr', toLoad{3}(2:end)]);
            cvara=     eval(['corr', toLoad{4}(2:end)]);
            
        end
        
        
        if ~ismember(variable, listRatingCombine)
            cvarmean = cvar;
            varmean = var;
        end
        
        
        %% define the limits of the plot
        x_ax_lim = [0 + rightReduction, trialonset+  leftReduction];
        
        
        invvar = invvar(vectorPat, :);
        var = var (vectorPat, :);
        cinvvar= cinvvar(vectorPat, :);
        cvar = cvar(vectorPat, :);
        
        varmean = varmean (vectorPat, :);
        cvarmean = cvarmean (vectorPat, :);
        
        siVar= size(var,1)/2;
        siVar2 = siVar;
        
        
        cvar1 = cvarmean(1:siVar, :);
        cvar2 = cvarmean(1+siVar:end, :);
        cinvvar1 = cinvvar(1:siVar, :);
        cinvvar2 = cinvvar(1+siVar:end, :);
        
        
        both = ((varmean(1:siVar,:)<sqrt(2*alpha)) & (varmean(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvarmean(1:siVar,:))==sign(cvarmean(1+siVar:end,:))));
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
        
        meancvar = abs(cvarmean(1:siVar,:) + cvarmean(1+siVar:end ,:))/2;
        meancinvvar = abs(cinvvar(1:siVar,:) + cinvvar(1+siVar:end ,:))/2;
        
        meancvarSel = abs(cvarSel(1:siVar,:) + cvarSel(1+siVar:end ,:))/2;
        meancinvvarSel = abs(cinvvarSel(1:siVar,:) + cinvvarSel(1+siVar:end ,:))/2;
        
        
        
        
        
        %% ticks
        
        Xtickspos=[[trialonset:-tickpossize: 1], [trialonset:tickpossize:size(cvar,2)]];  Xtickspos = sort(Xtickspos(2:end)); % from the trial aligment, go forward and backwards determining te position and sort
        Xtickslabelsnames = [[0:-1:-(length([trialonset:-tickpossize: 1])-1)], [0:1:length(trialonset:tickpossize:size(cvar,2))-1]]; Xtickslabelsnames = sort(Xtickslabelsnames(2:end));
        
        for hh=1:length(Xtickslabelsnames)
            aa{hh}=num2str(Xtickslabelsnames(hh));
        end
        Xtickslabelsnames=[];Xtickslabelsnames=aa;
        
        
        
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
        
        
        
        
        if strcmp(toLoad{1}(end-1:end), '_R')
            area(AreaRat, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            if xlabelon ==1
                xlabel('Time from stimulus onset (s)')
            end
            xlim(x_ax_lim)
        else
            
            area(AreaA1, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            area(AreaB, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.3,0.3,0.3], 'Edgecolor', 'none')
            area(AreaA2, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            if xlabelon ==1
                xlabel('Time relative to stimulus A1 onset (s)')
            end
            xlim(x_ax_lim)
        end
        
        
        % direct change plot
        
        %         toplotDC= varchange < alpha;
        %         [toplotDC2] = create_permVector (toplotDC, clustVarDC, pos, x_ax_lim, alphaPerm);
        %         plot(mean(toplotDC), 'color', [color2use(2,:), transparency], 'linewidth',1.5)
        %         plot(toplotDC2, 'color', [color2use(2,:)], 'linewidth', 1.5)
        
        %  proportion plot
        toplot2= ((var(1:siVar,:)<sqrt(2*alpha)) & (var(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvar(1:siVar,:))==sign(cvar(1+siVar:end,:))));
        [toplot] = create_permVector (toplot2, clustVar, pos, x_ax_lim, alphaPerm);
        plot(sum(toplot2)./siVar, 'color', [color2use(1,:), transparency], 'linewidth',1.5)
        plot(toplot, 'color', [color2use(1,:)], 'linewidth', 1.5)
        
        %binomial line
        line(x_ax_lim, [alpha,alpha], 'linestyle', '--','color','k', 'linewidth', 1.3)
        line(x_ax_lim, [pos,pos], 'linestyle', '--', 'color','r', 'linewidth', 1.3)
        
        % settings
        ylim([0, heightplot])
        set(gca, 'Xtick',Xtickspos, 'Xticklabels',Xtickslabelsnames);
        set(gca, 'FontSize',12, 'fontweight', 'bold')
        
        
        
        
        
        
        
        
        %% pval plot
        
        xaxis= 1:size(var,2);
        
        clear pvalSign zvalues
        for kk= 1:size(cvar,2)
            try
                [pvalSign(1, kk), ~, stats ] = signrank(meancvar(:,kk), meancinvvar(:,kk), 'tail', 'right');
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
                [pvalSignSel(1, kk), ~, stats] = signrank(meancvarSel(:,kk), meancinvvarSel(:,kk), 'tail', 'right');
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
        
        
        
        
        
        clear var cvar varmean  cvarmean varinv cvarinv
        
        %% autocorr plot
        
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
                toLoad= {['*', variable, 'inv_A1']; ['*', variable, '_A1']; [variable, '_', type_clusterSumEffect]; [variable, '_', typeclusterSumEffect_DirectChange]; [variable, '_signedRankclusterSumEffect']; [variable, '_signedRankSelclusterSumEffect']};
            end
            load([fold_auto,'/regresspval_',region,'.mat'],toLoad{1:2});
            load([fold_auto,'/regressreg_', region,'.mat'],toLoad{1:2});
            
            invvar =  eval(['pval', toLoad{1}(2:end)]);
            var=      eval(['pval', toLoad{2}(2:end)]);
            cinvvar=  eval(['corr', toLoad{1}(2:end)]);
            cvar=     eval(['corr', toLoad{2}(2:end)]);
            
            if ismember(variable, listRatingCombine)
                load([fold_auto,'/regresspval_',region,'.mat'],toLoad{3});
                load([fold_auto,'/regressreg_', region,'.mat'],toLoad{3});
                varmean=      eval(['pval', toLoad{3}(2:end)]);
                cvarmean=     eval(['corr', toLoad{3}(2:end)]);
            end
            
            if strcmp(variable, 'RT_R')
                toLoad2 = {['RT_Rat_', type_clusterSumEffect]; ['RT_Rat_', typeclusterSumEffect_DirectChange]; ['RT_Rat_signedRankclusterSumEffect']; ['RT_Rat_signedRankSelclusterSumEffect']};
            else
                toLoad2 = {[variable, '_', type_clusterSumEffect]; [variable, '_', typeclusterSumEffect_DirectChange]; [variable, '_signedRankclusterSumEffect']; [variable, '_signedRankSelclusterSumEffect']};
            end
            
            load([clustPermFold_auto,'/permdistribution_',region,'2.mat'],toLoad2{1:end});
            
            clustVar        =  eval([toLoad2{1}]);
            clustVarDC      =  eval([toLoad2{2}]);
            clustVarSR_All  =  eval([toLoad2{3}]);
            clustVarSR_Sel  =  eval([toLoad2{4}]);
            
            %
            %
            %             toLoad = {['*', variable,'_RawChange_A1']};
            %             if strcmp(variable, 'RT_R')
            %                 toLoad3 = {['*RT_RawChange_A1']};
            %             end
            %             load([DirectChangeFold_auto, 'regresspval_',region,'.mat'],toLoad3{1});
            %             load([DirectChangeFold_auto, 'regressreg_',region,'.mat'],toLoad3{1});
            %
            %             varchange  =   eval(['pval', toLoad{1}(2:end)]);
            %             %             cvarchange =   eval(['corr', toLoad3{1}(2:end)]);
        end
        
        if ~ismember(variable, listRatingCombine)
            1
            cvarmean = cvar;
            varmean = var;
        end
        
        
        
        
        invvar = invvar(vectorPat, :);
        var = var (vectorPat, :);
        cinvvar= cinvvar(vectorPat, :);
        cvar = cvar(vectorPat, :);
        
        varmean = varmean (vectorPat, :);
        cvarmean = cvarmean (vectorPat, :);
        
        siVar= size(var,1)/2;
        siVar2 = siVar;
        
        
        cvar1 = cvarmean(1:siVar, :);
        cvar2 = cvarmean(1+siVar:end, :);
        cinvvar1 = cinvvar(1:siVar, :);
        cinvvar2 = cinvvar(1+siVar:end, :);
        
        
        both = ((varmean(1:siVar,:)<sqrt(2*alpha)) & (varmean(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvarmean(1:siVar,:))==sign(cvarmean(1+siVar:end,:))));
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
        
        meancvar = abs(cvarmean(1:siVar,:) + cvarmean(1+siVar:end ,:))/2;
        meancinvvar = abs(cinvvar(1:siVar,:) + cinvvar(1+siVar:end ,:))/2;
        
        meancvarSel = abs(cvarSel(1:siVar,:) + cvarSel(1+siVar:end ,:))/2;
        meancinvvarSel = abs(cinvvarSel(1:siVar,:) + cinvvarSel(1+siVar:end ,:))/2;
        
        
        
        
        
        
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
        line(x_ax_lim, [alpha,alpha], 'linestyle', '--','color','k', 'linewidth', 1.3)
        line(x_ax_lim, [pos,pos], 'linestyle', '--', 'color','r', 'linewidth', 1.3)
        %         ylabel('Proportion of significant cells')
        
        
        
        toplot2= ((var(1:siVar,:)<sqrt(2*alpha)) & (var(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvar(1:siVar,:))==sign(cvar(1+siVar:end,:))));
        [toplot] = create_permVector (toplot2, clustVar, pos, x_ax_lim, alphaPerm);
        plot(toplot, 'color', [color2use(1,:)], 'linewidth', 1.5)
        
        % plot Areas
        if strcmp(toLoad{1}(end-1:end), '_R')
            area(AreaRat, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            if xlabelon ==1
                xlabel('Time from stimulus onset (s)')
            end
            xlim(x_ax_lim)
        else
            
            area(AreaA1, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            area(AreaB, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.3,0.3,0.3], 'Edgecolor', 'none')
            area(AreaA2, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            if xlabelon ==1
                xlabel('Time relative to stimulus A1 onset (s)')
            end
            xlim(x_ax_lim)
        end
        
        
        %         % direct change plot
        %         toplotDC= varchange < alpha;
        %         [toplotDC2] = create_permVector (toplotDC, clustVarDC, pos, x_ax_lim, alphaPerm);
        %         plot(mean(toplotDC), 'color', [color2use(2,:), transparency], 'linewidth',1.5)
        %         plot(toplotDC2, 'color', [color2use(2,:)], 'linewidth', 1.5)
        
        %  proportion plot
        
        toplot2= ((var(1:siVar,:)<sqrt(2*alpha)) & (var(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvar(1:siVar,:))==sign(cvar(1+siVar:end,:))));
        [toplot] = create_permVector (toplot2, clustVar, pos, x_ax_lim, alphaPerm);
        plot(sum(toplot2)./siVar, 'color', [color2use(1,:), transparency], 'linewidth',1.5)
        plot(toplot, 'color', [color2use(1,:)], 'linewidth', 1.5)
        
        %binomial line
        line(x_ax_lim, [alpha,alpha], 'linestyle', '--','color','k', 'linewidth', 1.3)
        line(x_ax_lim, [pos,pos], 'linestyle', '--', 'color','r', 'linewidth', 1.3)
        
        % settings
        ylim([0, heightplot])
        set(gca, 'Xtick',Xtickspos, 'Xticklabels',Xtickslabelsnames);
        set(gca, 'FontSize',12, 'fontweight', 'bold')
        
        
        if legon ==1
            %             if (ggg == 1 || ggg == 4) && (plotChangeProportion == 1) && onlyPlotChange == 0
            %                 leg =legend('Chance', 'Significance',  var2name{cc}, 'Change'); % [var2name{cc},'Change']
            %             else
            leg =legend('Chance', 'Significance',  var2name{cc});
            %             end
            set(leg, 'units','centimeters', 'Position', legpos ,'FontSize',12);
        end
        
        
        %% pval plot
        
        xaxis= 1:size(var,2);
        
        clear pvalSign zvalues
        for kk= 1:size(cvar,2)
            try
                [pvalSign(1, kk), ~, stats ] = signrank(meancvar(:,kk), meancinvvar(:,kk), 'tail', 'right');
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
        
        %%  clusterwise correction
        
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
                [pvalSignSel(1, kk), ~, stats] = signrank(meancvarSel(:,kk), meancinvvarSel(:,kk), 'tail', 'right');
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
        
        
        disp(['mean correlation over time between fraction of significant neurons and Satiation according to All SignedRank =  ' , num2str((sim1 + sim2)/2)])
        disp(['mean correlation over time between fraction of significant neurons and Satiation according to Selected SignedRank =  ' , num2str((simSel1 + simSel2)/2)])
        %          disp(['mean correlation over time between fraction of significant neurons and Satiation according to Direct Change =  ' , num2str((sim1DC1 + sim1DC2)/2)])
        
        
        cbh = colorbar;
        cbh.FontSize= 11;
        cbh.Units = 'centimeters';
        cbh.Position= cbhpos;
        cbh.FontWeight= 'bold';
        anot4 = annotation('textbox', [0,0,0.1,0.1],   'string',['pval'],'Units','normalized','FontName','Helvetica','FontSize',11,'EdgeColor','none','fontweight', 'bold');
        
        
        anot2 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',['All'],'Units','normalized','FontName','Helvetica','FontSize',12,'EdgeColor','none', 'fontweight', 'bold');
        anot3 = annotation('textbox',[0.00003    0+0.02   0.3    0.04],   'string',['Selected'],'Units','normalized','FontName','Helvetica','FontSize',12,'EdgeColor','none', 'fontweight', 'bold');
        
        global sizeletra sizetrialtype
        anot5 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',[letra],'Units','normalized','FontName','Helvetica','FontSize',sizeletra,'EdgeColor','none', 'fontweight', 'bold');
        
        if plottrialtype ==1
            anot7 = annotation('textbox',[0.00003    0+0.02   0.3    0.04],   'string',['Current Trial'],'Units','normalized','FontName','Helvetica','FontSize',sizetrialtype,'EdgeColor','none', 'fontweight', 'bold');
            anot8 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',['Next Trial'],'Units','normalized','FontName','Helvetica','FontSize',sizetrialtype,'EdgeColor','none', 'fontweight', 'bold');
            set(anot7, 'units', 'centimeters', 'position', anot7pos)
            set(anot8, 'units', 'centimeters', 'position', anot8pos)
        end
        
        set(anot2, 'units', 'centimeters', 'position', anot2pos)
        set(anot3, 'units', 'centimeters', 'position', anot3pos)
        set(anot4, 'units', 'centimeters', 'position', anot4pos)
        set(anot5, 'units', 'centimeters', 'position', anot5pos)
        
        
    end
end















function do_the_plot3(fig1, region2,  fold, clustPermFold,  DirectChangeFold,  var2use, var2name, pos11, pos12,pos13, ...
    cbhpos, anot2pos, anot3pos, anot4pos, anot5pos, legpos, color2use, letra, xlabelon, legon, ylabelon, regtype, plotChangeProportion, onlyPlotChange)

global alpha alphaBinom alphaPerm alphaSigned type_clusterSumEffect typeclusterSumEffect_DirectChange
set(0, 'currentFigure',fig1)


dbstop if error
%plot_paper_inverses_autocorr_reduced(0.05)



set(groot,'defaultAxesLinewidth', 1);
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


listRatingCombine = { 'Val', 'Sal', 'RT_R'};

listRating = {'ValTot', 'SalTot', 'Val', 'Sal', 'Valmean', 'Salmean', 'Rank', 'RTmean', 'RTstd', 'Valstd', 'RT_R', 'Sign', 'SignTot',   'Negsplit', 'Possplit', 'NegsplitTot', 'PossplitTot',  'ValSweetsplit', 'SalSweetsplit', 'ValSweetsplitTot', 'SalSweetsplitTot', 'RankSaltysplit', ...
    'ValSweetsplit', 'SalSweetsplit', 'ValSweetsplitTot', 'SalSweetsplitTot', 'RankSweetsplit'};


% for do_single = 0 % 0:2; %do_single==2 is is a single unit, do_single==1 is is a multi unit only, do_single==0 is is a  all SU and MU
%
%     for restricted = 0%:2

for g=1%:size(region2,1)
    
    region = region2{g}; load(['../Rasts matfiles/patientVector_', region]);
    vectorPat = eval(['patientVector_', region]); % extract patient ID of each cell
    vectorPat = ismember([vectorPat; vectorPat], [2,3,4,5,6,1,7,8,9,10,11,12]);
    
    
    
    for cc= 1%: size(var2use,2)
        
        variable = var2use{cc};
        disp([variable, '  ',  region])
        
        
        % change window of interest for valauation task
        if ismember(variable, listRating)
            leftReduction = (2000+2350)/windowstep-1;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr
            rightReduction =  trialonset-(2350)/windowstep;  % reduce the range of the plots 400 ms to reduce bad estimation of the mfr
        end
        
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
                toLoad= {['*', variable, 'inv_A1']; ['*', variable, '_A1']; [variable, '_', type_clusterSumEffect]; [variable, '_', typeclusterSumEffect_DirectChange]; [variable, '_signedRankclusterSumEffect']; [variable, '_signedRankSelclusterSumEffect']};
            end
            load([fold,'/regresspval_',region,'.mat'],toLoad{1:2});
            load([fold,'/regressreg_', region,'.mat'],toLoad{1:2});
            
            invvar =  eval(['pval', toLoad{1}(2:end)]);
            var=      eval(['pval', toLoad{2}(2:end)]);
            cinvvar=  eval(['corr', toLoad{1}(2:end)]);
            cvar=     eval(['corr', toLoad{2}(2:end)]);
            
            
            load([clustPermFold,'/permdistribution_',region,'2.mat'],toLoad{3:end});
            
            clustVar        =  eval([toLoad{3}]);
            clustVarDC      =  eval([toLoad{4}]);
            clustVarSR_All  =  eval([toLoad{5}]);
            clustVarSR_Sel  =  eval([toLoad{6}]);
            
            
            toLoad2 = {['*', variable,'_RawChange_A1']};
            load([DirectChangeFold, 'regresspval_',region,'.mat'],toLoad2{1});
            load([DirectChangeFold, 'regressreg_',region,'.mat'],toLoad2{1});
            
            varchange  =   eval(['pval', toLoad2{1}(2:end)]);
            %             cvarchange =   eval(['corr', toLoad2{1}(2:end)]);
            
        elseif regtype == 2
            toLoad= {['*', variable, 'reginv_A1']; ['*', variable, 'reg_A1']; ['*', variable, 'Autoreginv_A1']; ['*', variable, 'Autoreg_A1']};
            
            load([fold,'/regresspval_',region,'.mat'],toLoad{:});
            load([fold,'/regressreg_', region,'.mat'],toLoad{:});
            invvar =  eval(['pval', toLoad{1}(2:end)]);
            var=      eval(['pval', toLoad{2}(2:end)]);
            cinvvar=  eval(['corr', toLoad{1}(2:end)]);
            cvar=     eval(['corr', toLoad{2}(2:end)]);
            
            
            invvara =  eval(['pval', toLoad{3}(2:end)]);
            vara=      eval(['pval', toLoad{4}(2:end)]);
            cinvvara=  eval(['corr', toLoad{3}(2:end)]);
            cvara=     eval(['corr', toLoad{4}(2:end)]);
            
        end
        
        
        if ismember(variable, listRatingCombine)
            cvar =  eval(['corr', toLoad{3}(2:end)]);
        end
        
        
        %% define the limits of the plot
        x_ax_lim = [0 + rightReduction, trialonset+  leftReduction];
        
        
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
        
        meancvar = abs(cvar(1:siVar,:) + cvar(1+siVar:end ,:))/2;
        meancinvvar = abs(cinvvar(1:siVar,:) + cinvvar(1+siVar:end ,:))/2;
        
        meancvarSel = abs(cvarSel(1:siVar,:) + cvarSel(1+siVar:end ,:))/2;
        meancinvvarSel = abs(cinvvarSel(1:siVar,:) + cinvvarSel(1+siVar:end ,:))/2;
        
        
        
        
        
        %% ticks
        
        Xtickspos=[[trialonset:-tickpossize: 1], [trialonset:tickpossize:size(cvar,2)]];  Xtickspos = sort(Xtickspos(2:end)); % from the trial aligment, go forward and backwards determining te position and sort
        Xtickslabelsnames = [[0:-1:-(length([trialonset:-tickpossize: 1])-1)], [0:1:length(trialonset:tickpossize:size(cvar,2))-1]]; Xtickslabelsnames = sort(Xtickslabelsnames(2:end));
        
        for hh=1:length(Xtickslabelsnames)
            aa{hh}=num2str(Xtickslabelsnames(hh));
        end
        Xtickslabelsnames=[];Xtickslabelsnames=aa;
        
        
        
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
        line(x_ax_lim, [alpha,alpha], 'linestyle', '--','color','k', 'linewidth', 1.3)
        line(x_ax_lim, [pos,pos], 'linestyle', '--', 'color','r', 'linewidth', 1.3)
        if ylabelon
            ylabel('Proportion of significant cells')
        end
        
        
        
        if onlyPlotChange ||  plotChangeProportion
            toplotDC= varchange < alpha;
            [toplotDC2] = create_permVector (toplotDC, clustVarDC, pos, x_ax_lim, alphaPerm);
            plot(toplotDC2, 'color', [color2use(2,:)], 'linewidth', 1.5)
            
        end
        if onlyPlotChange ~= 1
            %  proportion plot
            toplot2= ((var(1:siVar,:)<sqrt(2*alpha)) & (var(siVar+1:end,:)<sqrt(2*alpha)) & (sign(cvar(1:siVar,:))==sign(cvar(1+siVar:end,:))));
            [toplot] = create_permVector (toplot2, clustVar, pos, x_ax_lim, alphaPerm);
            plot(toplot, 'color', [color2use(1,:)], 'linewidth', 1.5)
        end
        
        
        % plot Areas
        if strcmp(toLoad{1}(end-1:end), '_R')
            area(AreaRat, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            if xlabelon ==1
                xlabel('Time from stimulus onset (s)')
            end
            xlim(x_ax_lim)
        else
            
            area(AreaA1, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            area(AreaB, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.3,0.3,0.3], 'Edgecolor', 'none')
            area(AreaA2, [1,1], 0,  'Facealpha', 0.5, 'Facecolor', [0.6,0.6,0.6], 'Edgecolor', 'none')
            if xlabelon ==1
                xlabel('Time relative to stimulus A1 onset (s)')
            end
            xlim(x_ax_lim)
        end
        
        
        
        % direct change plot
        
        
        if onlyPlotChange ||  plotChangeProportion
            toplotDC= varchange < alpha;
            [toplotDC2] = create_permVector (toplotDC, clustVarDC, pos, x_ax_lim, alphaPerm);
            plot(mean(toplotDC), 'color', [color2use(2,:), transparency], 'linewidth',1.5)
            plot(toplotDC2, 'color', [color2use(2,:)], 'linewidth', 1.5)
            
            
            range2use = trialonset+3600/windowstep:trialonset+((3600+2350-1)/windowstep);
            num121 =  mean(mean(toplotDC(:, range2use)));
            p121 =  myBinomTest(num121*siVar, siVar, alpha, 'one');
            disp (['total fraction of units for direct change for variable  ', variable, ':  ', num2str(num121), '  Pvalue binomial test: ', num2str(p121) ])
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
        
        
        
        %% pval plot
        
        if (onlyPlotChange ~= 1)
            
            xaxis= 1:size(var,2);
            
            clear pvalSign zvalues
            for kk= 1:size(cvar,2)
                try
                    [pvalSign(1, kk), ~, stats ] = signrank(meancvar(:,kk), meancinvvar(:,kk), 'tail', 'right');
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
                    [pvalSignSel(1, kk), ~, stats] = signrank(meancvarSel(:,kk), meancinvvarSel(:,kk), 'tail', 'right');
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
            
            
            
            cbh = colorbar;
            cbh.FontSize= 11;
            cbh.Units = 'centimeters';
            cbh.Position= cbhpos;
            cbh.FontWeight= 'bold';
            anot4 = annotation('textbox', [0,0,0.1,0.1],   'string',['pval'],'Units','normalized','FontName','Helvetica','FontSize',11,'EdgeColor','none','fontweight', 'bold');
            
            
            anot2 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',['All'],'Units','normalized','FontName','Helvetica','FontSize',12,'EdgeColor','none', 'fontweight', 'bold');
            anot3 = annotation('textbox',[0.00003    0+0.02   0.3    0.04],   'string',['Selected'],'Units','normalized','FontName','Helvetica','FontSize',12,'EdgeColor','none', 'fontweight', 'bold');
            set(anot2, 'units', 'centimeters', 'position', anot2pos)
            set(anot3, 'units', 'centimeters', 'position', anot3pos)
            set(anot4, 'units', 'centimeters', 'position', anot4pos)
        end
        
        
        if legon ==1
            if (plotChangeProportion == 1) && onlyPlotChange == 0
                leg =legend('Chance', 'Significance',  'Change', var2name{cc}); % [var2name{cc},'Change']
            else
                leg =legend('Chance', 'Significance',  var2name{cc});
            end
            set(leg, 'units','centimeters', 'Position', legpos ,'FontSize',12);
        end
        
        global sizeletra sizetrialtype
        anot5 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',[letra],'Units','normalized','FontName','Helvetica','FontSize',sizeletra,'EdgeColor','none', 'fontweight', 'bold');
        set(anot5, 'units', 'centimeters', 'position', anot5pos)
        
        
    end
end
