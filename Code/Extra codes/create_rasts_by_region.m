function create_rasts_by_region(inputRegion, savefolder)


% create_rasts_by_region('A', './Rasts matfiles')

dbstop if error
assert(nargin==2, 'this function needs 2 inputs, one with the region code, and the second with the foldername to save all the data');


satiationfolder='/media/Projects/Alex/Reclustered analysis All 2';
cd (satiationfolder)

namefiles=[];
DirName=dir;
behavDir=[];


focalPIDL  = [01,03,04,06,08,09,10,11,12];
focalPIDR  = [02,04,05,07];
contraPIDL = [02,05,07];
contraPIDR = [01,03,06,08,09,10,11,12];




for a=1:length(DirName)
    
    if regexp(DirName(a).name, '[0-9][0-9][0-9]')
        
        if length(DirName(a).name)==3
            %         if isempty(regexp(DirName(a).name, '039', 'once')) % don't analyze pat 39, he did not eat.
            
            cd (DirName(a).name)
            PID =  str2num(DirName(a).name);
            
            pathfolders=dir;
            for i=1:length(pathfolders)
                
                if regexp(pathfolders(i).name, 'Rast__')
                    
                    cd (pathfolders(i).name)
                    filesfolder=dir;
                    
                    for j=1:length(filesfolder)
                        if ~isempty(regexp(filesfolder(j).name, '_Ra_'))
                            
                            if strcmp(inputRegion(1),'L')
                                if  ~isempty(regexp(filesfolder(j).name, ['L([AMP])?' ,char(inputRegion(2:end)), '___'], 'once'))
                                    behavDir{end+1}=DirName(a).name;
                                    namefiles{end+1}=[satiationfolder, '/', DirName(a).name, '/', pathfolders(i).name '/', filesfolder(j).name];
                                    
                                end
                                
                            elseif strcmp(inputRegion(1),'R')
                                if   ~isempty(regexp(filesfolder(j).name, ['R([AMP])?' ,char(inputRegion(2:end)), '___'], 'once'))
                                    behavDir{end+1}=DirName(a).name;
                                    namefiles{end+1}=[satiationfolder, '/', DirName(a).name, '/', pathfolders(i).name '/', filesfolder(j).name];
                                    
                                end
                                
                            elseif length(inputRegion)>2 && strcmp(inputRegion(1:3),'Con')
                                if  ~isempty(regexp(filesfolder(j).name, ['L([AMP])?' ,char(inputRegion(4:end)), '___'], 'once')) && ismember(PID, contraPIDL)
                                    behavDir{end+1}=DirName(a).name;
                                    namefiles{end+1}=[satiationfolder, '/', DirName(a).name, '/', pathfolders(i).name '/', filesfolder(j).name];
                                    
                                end
                                
                                if  ~isempty(regexp(filesfolder(j).name, ['R([AMP])?' ,char(inputRegion(4:end)), '___'], 'once')) && ismember(PID, contraPIDR)
                                    behavDir{end+1}=DirName(a).name;
                                    namefiles{end+1}=[satiationfolder, '/', DirName(a).name, '/', pathfolders(i).name '/', filesfolder(j).name];
                                    
                                end
                                
                            elseif length(inputRegion)>2 && strcmp(inputRegion(1:3),'Foc')
                                if  ~isempty(regexp(filesfolder(j).name, ['L([AMP])?' ,char(inputRegion(4:end)), '___'], 'once')) && ismember(PID, focalPIDL)
                                    behavDir{end+1}=DirName(a).name;
                                    namefiles{end+1}=[satiationfolder, '/', DirName(a).name, '/', pathfolders(i).name '/', filesfolder(j).name];
                                    
                                end
                                
                                if  ~isempty(regexp(filesfolder(j).name, ['R([AMP])?' ,char(inputRegion(4:end)), '___'], 'once')) && ismember(PID, focalPIDR)
                                    behavDir{end+1}=DirName(a).name;
                                    namefiles{end+1}=[satiationfolder, '/', DirName(a).name, '/', pathfolders(i).name '/', filesfolder(j).name];
                                    
                                end
                                
                                
                                
                                
                            elseif strcmp(inputRegion,'ALL')
                                
                                if	isempty(regexp(filesfolder(j).name, ['TO'], 'once'))
                                    behavDir{end+1}=DirName(a).name;
                                    namefiles{end+1}=[satiationfolder, '/', DirName(a).name, '/', pathfolders(i).name '/', filesfolder(j).name];
                                end
                                
                                
                            elseif ~(strcmp(inputRegion(1),'L') || strcmp(inputRegion(1),'R') || (length(inputRegion)>2 && strcmp(inputRegion(1:3),'Con')) || (length(inputRegion)>2 && strcmp(inputRegion(1:3),'Foc')))
                                
                                if ~strcmp(inputRegion,'H')
                                    if ~isempty(regexp(filesfolder(j).name, ['R' char(inputRegion) '___'], 'once')) || ~isempty(regexp(filesfolder(j).name, ['L' ,char(inputRegion), '___'], 'once'))
                                        behavDir{end+1}=DirName(a).name;
                                        namefiles{end+1}=[satiationfolder, '/', DirName(a).name, '/', pathfolders(i).name '/', filesfolder(j).name];
                                        
                                    end
                                elseif strcmp(inputRegion,'H')
                                    if ~isempty(regexp(filesfolder(j).name, ['R[AMP]' char(inputRegion) '___'], 'once')) || ~isempty(regexp(filesfolder(j).name, ['L[AMP]' ,char(inputRegion), '___'], 'once'))
                                        behavDir{end+1}=DirName(a).name;
                                        namefiles{end+1}=[satiationfolder, '/', DirName(a).name, '/', pathfolders(i).name '/', filesfolder(j).name];
                                        
                                    end
                                end
                                
                                
                                
                            end
                        end
                    end
                    
                end
                
            end
            cd ../..
            
            %             end
        end
    end
end
%  behavDir(:,cellfun(@isempty,behavDir)) = [];


aspdfiuhasifhaspdf

[rasts, outfiles, ~, positionLeft, positionRight, positionContra, positionFocal, posPID]=load_rasts(namefiles);


mkdir([satiationfolder,'/', savefolder, '/'])


if  ismember(inputRegion, {'A', 'H', 'EC', 'PHC'})
    save ([satiationfolder,'/', savefolder, '/positionVectors_',char(inputRegion), ], 'positionLeft', 'positionRight', 'positionContra', 'positionFocal', 'PID')
end



h=0;
while h<1
    try
        save ([satiationfolder,'/', savefolder, '/',char(inputRegion),'_rasts6.mat'], 'rasts', 'outfiles', '-v6')
        h=1;
    catch
    end
end

cd(satiationfolder)





% toDel

% % % % if isequal(inputregion(1:3), 'Foc') ||isequal(inputregion(1:3), 'Con')
% % % %
% % % %     h=0;
% % % %     while h<1
% % % %         try
% % % %             save ([satiationfolder,'/', savefolder, '/',char(inputRegion([1:3, 4:end])),'_rasts6.mat'], 'rasts', 'outfiles', '-v6')
% % % %             h=1;
% % % %         catch
% % % %         end
% % % %     end
% % % % else


