function extract_units

load cluster_info.mat

namefile=pwd;
[~,numPat] = fileparts(namefile);

numpat2=numPat(1:3);
unitVectorLA=[];
unitVectorRA=[];

ChanPIDvect={};

fornum=1;
for h=1:size(cluster_info,2)
     if ~ isempty(regexp(cluster_info{2,h}, 'RA[1-9]'))
         unitVectorRA = [unitVectorRA(:)', cluster_info{1,h}]';
         for j=1:length( cluster_info{1,h})
         ChanPIDvect{fornum}=[numpat2, cluster_info{2,h}];
         fornum=fornum+1;
         end
     end
     if ~ isempty(regexp(cluster_info{2,h}, 'LA[1-9]'))
         unitVectorLA = [unitVectorLA(:)', cluster_info{1,h}]';
         for j=1:length( cluster_info{1,h})
             ChanPIDvect{fornum}=[numpat2, cluster_info{2,h}];
             fornum=fornum+1;
         end
     end
    
end

eval(['unitVector' numpat2 ' = [unitVectorLA ; unitVectorRA];'])

eval(['ChanPIDvect' numpat2 ' = ChanPIDvect;'])


