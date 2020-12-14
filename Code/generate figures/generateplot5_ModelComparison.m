function generateplot5_ModelComparison

%% This creates the Figure S4
%

dbstop if error
set(groot,'defaultAxesLinewidth', 1)
saveFolder = '/media/Projects/Alex/Reclustered analysis All/PaperPlots';
mkdir(saveFolder)

possible_colors =  [[0    0.4470    0.7410];...
    [0.9290    0.5940    0.1250];...
    [0.9290    0.3250    0.0980];...
    [0.4940    0.1840    0.5560];...
    [0.4660    0.6740    0.1880];...
    [0.3010    0.7450    0.9330];...
    [0.6350    0.0780    0.1840]];

 


%% Regional analysis of different variables

% nutritional features AND Value

close all
markerstars= 3.7;
heigthstar = 0.5/100;
shiftletra=2;

star3dist = 0.025;
star2dist = 0.012;

figpos = [0,0, 39, 23.4];



toplot2=[[9.88, 9.88, 9.88, 9.28, 9.16, 9.04,  8.55];...
    [12.45, 10.94, 10.94, 10.94, 10.94, 11.32,  9.43];...   
    [9.18, 9.51,  9.51, 8.85, 8.52, 7.87, 7.87];...
    [8.55, 9.21,  9.21, 9.21, 7.89,  8.55, 7.89];...
    [6.48, 7.41,  7.41, 6.48, 7.41, 6.48,  6.48]]./100;


close all
pos1= [2.1 13.0303 36.3920 8.7829];
pos2= [2.1 0.8318 36.3920 8.7829];
fig2 =  figure('color', 'w', 'units', 'centimeters', 'position', figpos, 'visible', 'on');

ax1 = subplot('position',  pos1);
set(ax1, 'units', 'centimeters', 'position', pos1)



sizeletra = 30;
anot5 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',['A'],'Units','normalized','FontName','Helvetica','FontSize',sizeletra,'EdgeColor','none', 'fontweight', 'bold');

anot5pos =[pos1(1)-shiftletra,pos1(2)+pos1(4)-0.2, 2, 2];
set(anot5, 'units', 'centimeters', 'position', anot5pos)
                        
                            
barplot = bar(toplot2);
line([0.5, 5.5], [0.05,0.05], 'color', 'k', 'linestyle', '--', 'linewidth', 1.5)


barplot(1).FaceColor  = possible_colors(1,:);
barplot(2).FaceColor  = possible_colors(5,:);
barplot(3).FaceColor  = possible_colors(2,:);
barplot(4).FaceColor  = possible_colors(6,:);
barplot(5).FaceColor  = possible_colors(3,:);
barplot(6).FaceColor  = possible_colors(7,:);
barplot(7).FaceColor  = possible_colors(4,:);

          

%
%
% barplot(1,2).FaceColor = [0.9,0.9,0.9];
% barplot(1,1).FaceColor = [0.45,0.45,0.45];
% barplot(1,1).FaceColor = [0,0,0];


set(gca, 'xticklabels', {'ALL', 'A', 'H','EC','PHC'})

box off



% barWhiskerBridge(toplot2, [zeros(numel(toplot2),1)] )
ylabel('proportion of neurons')
xlabel('')

ylim([0, 0.171])
xlim([0.5, 5.5])



vecshift = [-0.345,-0.23, -0.115, 0, + 0.115, +0.23,+.345];

vvecXpos= [[1,1,1,1,1,1,1] + vecshift;...
    [1,1,1,1,1,1,nan]*2+ vecshift;...
    [nan,1,1,nan,nan,1,1]*3 + vecshift;...
    [1,1,1,1,nan,1,nan]*4 + vecshift;...
    [nan,nan,nan,nan,nan, nan,nan]*5 + vecshift];

hold on
for gg= 1:35
    plot(vvecXpos(gg),toplot2(gg)+heigthstar, 'k','Marker', 'p', 'MarkerFaceColor', 'k', 'Markersize', markerstars)
end



vecshift1 = vecshift - star3dist;

vvecXpos= [[1,1,1,1,1,1,1] + vecshift1;...
    [1,1,1,1,1,1,nan]*2+ vecshift1;...
    [nan,1,1,nan,nan,nan,nan]*3 + vecshift1;...
    [nan,nan,nan,nan,nan, nan,nan]*4 + vecshift1;...
    [nan,nan,nan,nan,nan, nan,nan]*5 + vecshift1];




hold on
for gg= 1:35
    plot(vvecXpos(gg),toplot2(gg)+heigthstar, 'k','Marker', 'p', 'MarkerFaceColor', 'k', 'Markersize', markerstars)
end


vecshift1 = vecshift+ star3dist;

vvecXpos= [[1,1,1,1,1,1,1] + vecshift1;...
    [1,1,1,1,1,1,nan]*2+ vecshift1;...
    [nan, 1,1,nan,nan,nan,nan]*3 + vecshift1;...
    [nan,nan,nan,nan,nan, nan,nan]*4 + vecshift1;...
    [nan,nan,nan,nan,nan, nan,nan]*5 + vecshift1];


hold on
for gg= 1:35
    plot(vvecXpos(gg),toplot2(gg)+ heigthstar, 'k','Marker', 'p', 'MarkerFaceColor', 'k', 'Markersize', markerstars)
end






%% 2 stars
vecshift1 = vecshift - star2dist;


vvecXpos= [[ nan,nan,nan,nan,nan, nan,nan] + vecshift1;...
    [nan,nan,nan,nan,nan, nan,1]*2 + vecshift1;...
    [1,nan,nan,1,1,nan,nan]*3 + vecshift1;...
    [nan,nan,nan,nan,nan, nan,nan]*4 + vecshift1;...
    [nan,nan,nan,nan,nan, nan,nan]*5 + vecshift1];

hold on
for gg= 1:35
    plot(vvecXpos(gg),toplot2(gg)+heigthstar, 'k','Marker', 'p', 'MarkerFaceColor', 'k', 'Markersize', markerstars)
end

vecshift1 = vecshift + star2dist;

vvecXpos= [[ nan,nan,nan,nan,nan, nan,nan] + vecshift1;...
    [nan,nan,nan,nan,nan, nan,1]*2 + vecshift1;...
    [1,nan,nan,1,1,nan,nan]*3 + vecshift1;...
    [nan,nan,nan,nan,nan, nan,nan]*4 + vecshift1;...
    [nan,nan,nan,nan,nan, nan,nan]*5 + vecshift1];

hold on
for gg= 1:35
    plot(vvecXpos(gg),toplot2(gg)+heigthstar, 'k','Marker', 'p', 'MarkerFaceColor', 'k', 'Markersize', markerstars)
end



leg= legend({'Taste', 'Value Ranking', 'Value Rating',...
    'Sugar', 'Salt', 'Fiber', ...
    'Unsigned Value'}, 'NumColumns',7);
%  22.2483 20.9021 0.8464
% legpos=[anot5pos =[pos1(1)-shiftletra,pos1(2)+pos1(4),2,2];
% set(anot5, 'units', 'centimeters', 'position', anot5pos)]
legpos = [10  pos1(2)+pos1(4)    pos1(2)+pos1(4)-1    0.8464]; 
set(leg,'units', 'centimeters','Position', legpos ,'FontSize',12); % 'units','centimeters',
set(gca, 'FontSize',12, 'fontweight', 'bold')


toplot2=[13.01, 11.81,  9.76,   9.16,  9.04,  7.95, 7.95;...
            12.08, 11.32,  8.68,   9.06,  8.30,  8.30, 7.17; ...
            13.44, 11.15,  9.84,   8.85,  8.52,  6.89,  7.21;...
            13.82, 15.13, 12.50,  10.53, 11.84, 10.53, 10.53; ...
            11.11,  9.26,  7.41,   6.48,  7.41,  4.63, 6.48;...
            ]./100;

ax2 = subplot('position',  pos2);
set(ax2, 'units', 'centimeters', 'position', pos2)


sizeletra = 30;
anot5 = annotation('textbox',[0.02    0+0.02   0.3    0.04],   'string',['B'],'Units','normalized','FontName','Helvetica','FontSize',sizeletra,'EdgeColor','none', 'fontweight', 'bold');

anot5pos =[pos2(1)-shiftletra,pos2(2)+pos2(4)-0.2, 2 2];
set(anot5, 'units', 'centimeters', 'position', anot5pos)
                   
barplot = bar(toplot2);
line([0.5, 5.5], [0.05,0.05], 'color', 'k', 'linestyle', '--', 'linewidth', 1.5)



barplot(1).FaceColor  = possible_colors(2,:);
barplot(2).FaceColor  = possible_colors(1,:);
barplot(3).FaceColor  = possible_colors(3,:);
barplot(4).FaceColor  = possible_colors(4,:);
barplot(5).FaceColor  = possible_colors(5,:);
barplot(6).FaceColor  = possible_colors(6,:);
barplot(7).FaceColor  = possible_colors(7,:);


set(gca, 'xticklabels', {'ALL', 'A',  'H', 'EC','PHC'})

box off



% barWhiskerBridge(toplot2, [zeros(numel(toplot2),1)] )
ylabel('proportion of neurons')
xlabel('')

ylim([0, 0.171])
xlim([0.5, 5.5])


vvecXposid= [[3,3,3,3,3,3,3];...
    [3,3,2,2,1,1,0];...
    [3,3,3,2,2,0,0];...
    [3,3,3,2,3,2,2];...
    [2,1,0,0,0,0,0]];


vecshift = [-0.345,-0.23, -0.115, 0, + 0.115, +0.23,+.345];
vecshift1 = vecshift;
vvecXpos= [(vvecXposid(1,:)==1 |vvecXposid(1,:)==3)  + vecshift1;...
    (vvecXposid(2,:)==1 |vvecXposid(2,:)==3)*2+ vecshift1;...
    (vvecXposid(3,:)==1 |vvecXposid(3,:)==3)*3 + vecshift1;...
    (vvecXposid(4,:)==1 |vvecXposid(4,:)==3)*4 + vecshift1;...
    (vvecXposid(5,:)==1 |vvecXposid(5,:)==3)*5 + vecshift1];

hold on
for gg= 1:35
    plot(vvecXpos(gg),toplot2(gg)+heigthstar, 'k','Marker', 'p', 'MarkerFaceColor', 'k', 'Markersize', markerstars)
end



vecshift1 = vecshift - star3dist;
vvecXpos= [(vvecXposid(1,:)==3)  + vecshift1;...
    (vvecXposid(2,:)==3)*2+ vecshift1;...
    (vvecXposid(3,:)==3)*3 + vecshift1;...
    (vvecXposid(4,:)==3)*4 + vecshift1;...
    (vvecXposid(5,:)==3)*5 + vecshift1];

vvecXpos(vvecXpos==0)=nan;

hold on
for gg= 1:35
    plot(vvecXpos(gg),toplot2(gg)+heigthstar, 'k','Marker', 'p', 'MarkerFaceColor', 'k', 'Markersize', markerstars)
end

vecshift1 = vecshift+ star3dist;
vvecXpos= [(vvecXposid(1,:)==3)  + vecshift1;...
    (vvecXposid(2,:)==3)*2+ vecshift1;...
    (vvecXposid(3,:)==3)*3 + vecshift1;...
    (vvecXposid(4,:)==3)*4 + vecshift1;...
    (vvecXposid(5,:)==3)*5 + vecshift1];


hold on
for gg= 1:35
    plot(vvecXpos(gg),toplot2(gg)+ heigthstar, 'k','Marker', 'p', 'MarkerFaceColor', 'k', 'Markersize', markerstars)
end






%% 2 stars
vecshift1 = vecshift - star2dist;


vvecXpos= [(vvecXposid(1,:)==2)  + vecshift1;...
    (vvecXposid(2,:)==2)*2+ vecshift1;...
    (vvecXposid(3,:)==2)*3 + vecshift1;...
    (vvecXposid(4,:)==2)*4 + vecshift1;...
    (vvecXposid(5,:)==2)*5 + vecshift1];

hold on
for gg= 1:35
    plot(vvecXpos(gg),toplot2(gg)+heigthstar, 'k','Marker', 'p', 'MarkerFaceColor', 'k', 'Markersize', markerstars)
end

vecshift1 = vecshift + star2dist;


vvecXpos= [(vvecXposid(1,:)==2)  + vecshift1;...
    (vvecXposid(2,:)==2)*2+ vecshift1;...
    (vvecXposid(3,:)==2)*3 + vecshift1;...
    (vvecXposid(4,:)==2)*4 + vecshift1;...
    (vvecXposid(5,:)==2)*5 + vecshift1];

hold on
for gg= 1:35
    plot(vvecXpos(gg),toplot2(gg)+heigthstar, 'k','Marker', 'p', 'MarkerFaceColor', 'k', 'Markersize', markerstars)
end



leg = legend({'Reaction Time', 'Confidence', 'Unsigned Difference', ...
    'Chosen Value', 'Unchosen/Chosen Ratio', ...
    'Summed Value', 'Unchosen Value'}, 'NumColumns',7);

legpos = [16  pos2(2)+pos2(4)   pos2(2)+pos2(4)-1    0.8464]; 
set(leg, 'units', 'centimeters','Position', legpos ,'FontSize',12); % 'units','centimeters',
set(gca, 'FontSize',12, 'fontweight', 'bold')


toSave = [saveFolder, '/plotCombinedValNatandRTConf'];

print([toSave, '200'], '-dpng', '-r200');


print(toSave, '-dpng', '-r500');
print([toSave, '300'], '-dpng', '-r300');

print([toSave, 'Vec'], '-depsc', '-painters');
print([toSave, 'Vec'], '-depsc', '-opengl', '-r500');
set(fig2,'Units','centimeters');
pos = get(fig2,'Position');
set(fig2,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

print([toSave, 'Vec'], '-dpdf', '-painters');
print(toSave, '-dpdf', '-r500', '-opengl');
print([toSave, ' 300'], '-dpdf', '-r300', '-opengl');


name = [saveFolder, '/plotCombinedValNatandRTConf300']
file=imread([name, '.png']);  
imwrite(file,[name, '.tiff'],'tiff');
















%% Plot of thr 20 best variables
 
names={'RT',   'Confidence',   ...
    'Abs. Difference (Rk.)',   'Abs. Difference (Rt.)',...
    'Chosen Value (Rk.)',   'Unchosen/Chosen Ratio (Rk.)',   'Chosen Value (Rt.)',   'Chosen Sugar', ...
    'Summed Value (Rk.)',   'Unchosen Value (Rk.)',   'Chosen Taste',   'Summed Value (Rt.)',   'Unchosen Value (Rt.)',   ...
    'Stimulus B Unsigned Rating',   'Stimulus B Value (Rk.)',   'Unchosen/Chosen ratio (Rt.)',   'Stimulus B Fiber', ...
    'Chosen Fiber',   'Right Value (Rk.)',   'Left Value (Rk.)'};

pvalvsconf = [0.048197,0.038827,  0.020699,  0.024458, 0.017118, 0.0042656, 0.0021505, 0.0065388,  0.0040438,   0.0034148,0.0030594, 0.0004531,  0.0011222,  4.952E-05, 0.00063438, 0.00058628, 0.00013909, 0.00038215];
pvalvsRt = [  0.010753, 0.0085843, 0.0053396, 0.0035442, 0.0041589,0.00047215,  0.00047215,  0.00036812,   0.00018752, 0.0002433,   6.9521E-05, 3.2074E-05,  6.9521E-05, 3.2074E-05, 3.2074E-05, 4.2742E-05 ,  2.8699E-05,    1.4115E-05 ];

pvalvsconfBinary = [0,0,1     1     1     1     1     2     2     2     2     2     2     3     2     3     3     3     3     3];
pvalvsRtBinary    = [0,0,1     2     2     2     2     3     3     3     3     3     3     3     3     3     3     3     3     3];
toplot2 = [13.01, 11.81,  9.76, 9.52,  9.16, 9.04,  8.92, 8.07,   7.95, 7.95, 7.71, 7.71, 7.71, ...
    7.47,  7.47, 7.35,  7.23, 7.11, 7.11, 6.99] ./100;
toplot2pval = [3, 3,  3, 3,  3, 3,  3, 3,   3, 3, 3, 3, 3,      2,        2,       2,      2,      2, 2, 2];

markerstars= 6;%3.5;
heigthstar = 0.5;


star3dist = 0.2;
star2dist = 0.1;

figpos = [0,0, 26.6800, 20];


% toplot2=[13.0835, 11.8869,  9.8234,  9.5487,  9.1576,  9.1112,  8.9667];

close all
fig1 =  figure('color', 'w', 'units', 'centimeters', 'position', figpos, 'visible', 'on');

sb1 = subplot(1,1,1);
set(sb1, 'units', 'centimeters', 'position', [2.0243 5.4 23.9911 12.6767])

barplot = bar(toplot2, 'facecolor', possible_colors(4,:)); hold on
barplot = bar(1, toplot2(1), 'facecolor', possible_colors(2,:));
barplot = bar(2, toplot2(2), 'facecolor', possible_colors(1,:));
lineChance = line([0, 21], [0.05,0.05], 'color', 'k', 'linestyle', '--', 'linewidth', 1.5);

alpha = 0.05;
sizePlot = 830;

clear hh
for n=1:sizePlot
    hh(n)=1-binocdf(n-1,sizePlot,alpha);
end

pos= find(hh< 0.001); pos= pos(1)/sizePlot;
lineSignif = line([0, 21], [pos,pos], 'color', 'r', 'linestyle', '--', 'linewidth', 1.5);

        
sb1.XTick=[1:20];
sb1.XTickLabel=names;
set(gca,'XTickLabelRotation',45)

box off

ylabel('proportion of neurons')
xlabel('')

ylim([0, 0.15])
xlim([0, 21])



vvecXposid= toplot2pval;%[3,3,3,3,3,3,3];





vecshift = [1:20];

%% Confidence comparison


vvecXposid = pvalvsconfBinary;

colorstar = possible_colors(1,:); % 'r';
heigthstar = 0.01;

vecshift = [1:20];
vecshift1 = vecshift;
vvecXpos= [(vvecXposid(1,:)==1 |vvecXposid(1,:)==3)  .* vecshift1];

vvecXpos(vvecXpos==0)=nan;

hold on
for gg= 1:20
    pstarConfidence = plot(vvecXpos(gg),toplot2(gg)+heigthstar, 'color', colorstar,'Marker', 'p', 'MarkerFaceColor', colorstar, 'Markersize', markerstars);
end



vecshift1 = vecshift - star3dist;
vvecXpos= [(vvecXposid(1,:)==3)  .* vecshift1];

vvecXpos(vvecXpos==0)=nan;

hold on
for gg= 1:20
    plot(vvecXpos(gg),toplot2(gg)+heigthstar, 'color', colorstar,'Marker', 'p', 'MarkerFaceColor', colorstar, 'Markersize', markerstars)
end

vecshift1 = vecshift+ star3dist;
vvecXpos= [(vvecXposid(1,:)==3)  .* vecshift1];
vvecXpos(vvecXpos==0)=nan;

hold on
for gg= 1:20
    plot(vvecXpos(gg),toplot2(gg)+ heigthstar, 'color', colorstar,'Marker', 'p', 'MarkerFaceColor', colorstar, 'Markersize', markerstars)
end






%% 2 stars
vecshift1 = vecshift - star2dist;


vvecXpos= [(vvecXposid(1,:)==2)  .* vecshift1];

vvecXpos(vvecXpos==0)=nan;
hold on
for gg= 1:20
   plot(vvecXpos(gg),toplot2(gg)+heigthstar, 'color', colorstar,'Marker', 'p', 'MarkerFaceColor', colorstar, 'Markersize', markerstars)
end

vecshift1 = vecshift + star2dist;


vvecXpos= [(vvecXposid(1,:)==2)  .* vecshift1];

vvecXpos(vvecXpos==0)=nan;
hold on
for gg= 1:20
    plot(vvecXpos(gg),toplot2(gg)+heigthstar, 'color', colorstar,'Marker', 'p', 'MarkerFaceColor', colorstar, 'Markersize', markerstars)
end






%% RT


vvecXposid = pvalvsRtBinary;

colorstar = possible_colors(2,:);%'b';
heigthstar = 0.015;

vecshift = [1:20];
vecshift1 = vecshift;
vvecXpos= [(vvecXposid(1,:)==1 |vvecXposid(1,:)==3)  .* vecshift1];

vvecXpos(vvecXpos==0)=nan;

hold on
for gg= 1:20
     pstarRT =plot(vvecXpos(gg),toplot2(gg)+heigthstar, 'color', colorstar,'Marker', 'p', 'MarkerFaceColor',  colorstar, 'Markersize', markerstars);
end



vecshift1 = vecshift - star3dist;
vvecXpos= [(vvecXposid(1,:)==3)  .* vecshift1];

vvecXpos(vvecXpos==0)=nan;

hold on
for gg= 1:20
    plot(vvecXpos(gg),toplot2(gg)+heigthstar, 'color', colorstar,'Marker', 'p', 'MarkerFaceColor',  colorstar, 'Markersize', markerstars)
end

vecshift1 = vecshift+ star3dist;
vvecXpos= [(vvecXposid(1,:)==3)  .* vecshift1];
vvecXpos(vvecXpos==0)=nan;

hold on
for gg= 1:20
    plot(vvecXpos(gg),toplot2(gg)+ heigthstar, 'color', colorstar,'Marker', 'p', 'MarkerFaceColor', colorstar, 'Markersize', markerstars)
end






%% 2 stars
vecshift1 = vecshift - star2dist;


vvecXpos= [(vvecXposid(1,:)==2)  .* vecshift1];

vvecXpos(vvecXpos==0)=nan;
hold on
for gg= 1:20
    plot(vvecXpos(gg),toplot2(gg)+heigthstar, 'color', colorstar,'Marker', 'p', 'MarkerFaceColor', colorstar, 'Markersize', markerstars)
end

vecshift1 = vecshift + star2dist;


vvecXpos= [(vvecXposid(1,:)==2)  .* vecshift1];

vvecXpos(vvecXpos==0)=nan;
hold on
for gg= 1:20
    plot(vvecXpos(gg),toplot2(gg)+heigthstar, 'color', colorstar,'Marker', 'p', 'MarkerFaceColor', colorstar, 'Markersize', markerstars)
end


leg= legend( [pstarRT,  pstarConfidence, lineSignif, lineChance],{'Model comparison, RT better', 'Model comparison, Confidence better', 'Significance level, binomial test' , 'Chance level'});

set(gca, 'FontSize',12, 'fontweight', 'bold')
 set(leg,'FontSize',12);


toSave = [saveFolder, '/plotRTConfidencAlternative2'];




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
%}