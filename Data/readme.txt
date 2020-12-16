All neural and behavioural data can be found in the following link: 

https://zenodo.org/record/4327565



==========================================================================================================================================================================================================
%
% Here is some more info about the content of other relevant variables
% within the code but you'll probably not need this info if you're just
% running your analysis:
%
% Satcue columns:
% 1......Trial order (just ascending numbers from 1 through as many images as were shown)
% 2......The ordinal number of the left image on screen (i.e. if you have a total of 20 images in your 'Used Stimuli' foder then this will be a number between 1 and 20 irrespective of the respective file name).
% 3......The ordinal number of the right image on screen (i.e. if you have a total of 20 images in your 'Used Stimuli' foder then this will be a number between 1 and 20 irrespective of the respective file name).
% 4......The file name of the left image on screen
% 5......The file name of the right image on screen
% 6......
% 7......
% 8......Show left image first (1) or right image first (0) 
% 9......Actual measured time fixation cross was on screen
% 10.....Actual measured time between fixation cross and stimulus presentation
% 11.....Actual measured time stimuli were on screen
% 12.....Actual measured time between the stimuli and the presentation of the Likert scales
% 13.....Actual measured reaction time (time used to answer the Likert scale)
% 14.....The subject's response in the Likert scale (negative number = left, positive number = right as long as your srand variable is always equal to 2, which it should be here)
% 15.....
% 16.....Continuous interpretation of the Likert scales whether sweet (image number < 100) or salty (image number > 100) was chosen. Here a sweet choice is represented as positive numbers and a salty choice as negative numbers. Beware though that the choice is not necessarily mutually exclusive as it is possible for that two sweet stimuli were shown at the same time!
% 17.....Binary interpretation of the Likert scales whether sweet (image number < 100) or salty (image number > 100) was chosen. Here a sweet choice is represented as 1 and a salty choice as 0. Beware though that the choice is not necessarily mutually exclusive as it is possible for that two sweet stimuli were shown at the same time!
%
%satcue for rating trials only 13 first columns, 
%values form likert scale 9 steps, basically from -100 to +100 in 25
%unit-steps, thus 9 steps



%=======================================================================================================
% Ranking and satcue variable has that 1 and 3 empty, 2-rating1, 3-2AC1, 5-rating2, 6-2AFC2. 
%=======================================================================================================

%in ranking already everything by stimulus, not by trial as in satcue

% Ranking of 2AFC part:    %%%%% inside variables 3 and 6
% 1 ordinal
% 2 Stimulus ID
% 3 rank
% 4 rating (absolut deviation in the direction of this individual stimulus, 0 - 1900 (in steps of 25)?)
%
% Ranking of rating part:   %%%inside varaible 2 and 5
% 1 ordinal
% 2 Stimulus ID
% 3 number of selections (0 - 3, +1 for each presentation & pos rating)
% 4 absolute rating over 3 trials (-300 - 300 (in steps of 25), good to bad)
%



%% 
% TTLs of Satiation experiment:
% TTL 43 = Likert is on


%%%

%     cond_ttls = [...
%         7,0;...     % Any sat pic shown
%         7,1;...     % 7 after 3                                     All judge pics
%         7,2;...     % 7 after 4,5,6                                 All comp pics
%         7,3;...     % 7 after 3 without 4,5,6 before it             First block judge pics
%         7,4;...     % 7 after 3 with 4,5,6 before it                Second block judge pics
%         7,5;...     % 7 after 4,5,6 with 7 after 3 after it         First block comp pics
%         7,6;...     % 7 after 4,5,6 without 7 after 3 after it      Second block comp pics
%         7,7;...     % 7 after 4                                     A1 pic during the comp part
%         7,8;...     % 7 after 5                                     B pic during the comp part
%         7,9;...     % 7 after 6                                     A2 pic during the comp part
%         7,10;...    % 7 after 4 or 6                                First or third satiation pic during the comp part
%      
%                
%         7,11;...    % 7 after 4 without 7 after 3 after it          A1 pic first block
%         7,12;...    % 7 after 4 with 7 after 3 after it             A1 pic second block
%         7,13;...    % 7 after 5 without 7 after 3 after it          B pic first block
%         7,14;...    % 7 after 5 with 7 after 3 after it             B pic second block
%         7,15;...    % 7 after 6 without 7 after 3 after it          A2 pic first block
%         7,16;...    % 7 after 6 with 7 after 3 after it             A2 pic second block
%                     
%         43,0;...    % Any 43                                        Any 2AFC Likert
%         43,17;...   % 43 with 7 after 3 after it                    Rating Likert first block 
%         43,18;...   % 43 without 7 after 3 after it                 Rating Likert first block 
%         44,0;...    % Any 44                                        Any rating Likert
%         44,19;...   % 44 without 4,5,6 before it                    Rating Likert first block 
%         44,20];     % 44 with 4,5,6 before it                       Rating Likert first block 
%%%


==========================================================================================================================================================================================================
% Rasts from neural data


% Columns in rast variable:
% 1     Experiment name
% 2     Patient number (randomized for privacy)
% 3     Session number
% 4     Recording date (deleted for privacy)
% 5     Recording time (deleted for privacy)
% 6     Brain area (recording location)
% 7     Channel number
% 8     Cluster number
% 9     Event number (TTL)      % 7 means picture shown, 43 means likert
% 10    TTL inclusion criteria  %from 0 to 10, each diff pics combinations
% 11    Time points of TTLs
% 12    Raster plot
% 13    Spike times for this cluster

