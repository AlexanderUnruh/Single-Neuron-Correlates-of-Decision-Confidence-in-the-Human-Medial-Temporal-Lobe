function [rast_mrg, files, rast_bsl_mrg, positionLeft, positionRight, positionContra, positionFocal, PID] = load_rasts(f)



dbstop if error

% location of focal andcontralateral sides of the epileptic focus

focalPIDL  = [1,3,4,6,8,9,10,11,12];
focalPIDR  = [2,4,5,7];
contraPIDL = [2,5,7];
contraPIDR = [1,3,6,8,9,10,11,12];




% Concatenate all rast files (but only the relevant units, TTLs and conditions)
rast_mrg = {};
rast_bsl_mrg = {};
files = cell(0,6);

positionLeft = [];
positionRight = [];
positionContra = [];
positionFocal = [];
PID = [];

for a = 1:length(f)
    clear rast
    load(f{a})
    
    Rastfilename = regexp(f{a}, '/', 'split');
    
    
    
    if exist('rast') && ~isempty(rast)
        rast_rsp = rast(ismember(cell2mat(rast(:,9)),7) | ismember(cell2mat(rast(:,9)),43) | ismember(cell2mat(rast(:,9)),44) | ismember(cell2mat(rast(:,9)),50)| ismember(cell2mat(rast(:,9)),51) | ismember(cell2mat(rast(:,9)),31) | ismember(cell2mat(rast(:,9)),32),:); %take all rast for evetn 7 ==image || event 43==likert-2AFC || event 44= likert rating
        rast_mrg = [rast_mrg;rast_rsp]; %includes not all the information of the rasts, so patient, experiment, PSTH...
        files(end+1,1:6) = rast(1,1:6);
        
        % this if-statements part creates a vector indicating if the cell is
        % right or left, focal or contralateral
        %numCells  = size(rast_rsp,1)/20;  % normal version
        numCells  = size(rast_rsp,1)/12;   % 20 seconds version
        
        %numCells  = size(rast_rsp,1)/23;
        
        PID = [PID(:); ones(numCells,1) * str2num(Rastfilename{end}(19:21))];
        if strcmp(Rastfilename{end}(1), 'L')
            
            positionLeft =  [positionLeft(:); ones(numCells,1)];
            positionRight= [positionRight(:); zeros(numCells,1)];
            
            if ismember(str2num(Rastfilename{end}(19:21)), focalPIDL)
                
                positionFocal= [positionFocal(:); ones(numCells,1)];
                positionContra= [positionContra(:); zeros(numCells,1)];
            end
            
            if ismember(str2num(Rastfilename{end}(19:21)), contraPIDL)
                positionContra= [positionContra(:); ones(numCells,1)];
                positionFocal= [positionFocal(:); zeros(numCells,1)];
            end
            
            
            
        elseif strcmp(Rastfilename{end}(1), 'R')
            
            positionRight =  [positionRight(:); ones(numCells,1)];
            positionLeft =  [positionLeft(:); zeros(numCells,1)];
            
            if ismember(str2num(Rastfilename{end}(19:21)), focalPIDR)
                positionFocal= [positionFocal(:); ones(numCells,1)];
                positionContra= [positionContra(:); zeros(numCells,1)];
            end
            
            if ismember(str2num(Rastfilename{end}(19:21)), contraPIDR)
                positionContra= [positionContra(:); ones(numCells,1)];
                positionFocal= [positionFocal(:); zeros(numCells,1)];
                
            end
        end
        
        if nargout>2
            rast_bsl = rast(ismember(cell2mat(rast(:,9)),7) & ismember(cell2mat(rast(:,10)),7),:); % choice A1 trials as baseline for the B1 analysis only ,as B and A2 have no baseline time: condition 7,7 == A1
            rast_bsl_mrg = [rast_bsl_mrg;rast_bsl];
        end
    end
end

rast_mrg(all(cellfun(@isempty,rast_mrg),2),:) = [];       % remove empty rows

if nargout>2
    rast_bsl_mrg(all(cellfun(@isempty,rast_bsl_mrg),2),:) = [];       % remove empty rows
end

