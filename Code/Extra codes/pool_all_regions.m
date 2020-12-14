function pool_all_regions(whattotpool)


% This function runs over all matfiles containing anaylses of the variables 
% and pools them or splits them as requested (contralateral, focal, all...)
% DO NOT CONFUSE WITH create_rasts_by_region (which comibnes the rasts files directly)
%

dbstop if error

folder_to_save= pwd;

%% read all regions and load pvalues and correlaton/beta coefficients in the matfiles; extract them keeping the structure and in presatiation post satiation.


if strcmp(whattotpool,'ALL')
    inputRegion= whattotpool;
    
    jaja={'A', 'H', 'EC', 'PHC'};
elseif strcmp(whattotpool,'L')
    inputRegion= whattotpool;
    
    jaja={'LA', 'LH', 'LEC', 'LPHC'};
elseif strcmp(whattotpool,'R')
    inputRegion= whattotpool;
    
    jaja={'RA', 'RH', 'REC', 'RPHC'};
    
elseif strcmp(whattotpool,'Con')
    inputRegion= whattotpool;
    
    jaja={'ConA', 'ConH', 'ConEC', 'ConPHC'};
    
elseif strcmp(whattotpool,'Foc')
    inputRegion= whattotpool;
    
    jaja={'FocA', 'FocH', 'FocEC', 'FocPHC'};
end



% filename = ['regresspval_PHC.mat'];
% m = matfile(filename);



load (['regresspval_PHC.mat'])
listpval=whos('pval*');
clear pval*


if exist([ '/regresspval_' inputRegion ,'.mat'])  || exist([folder_to_save '/regressreg_' inputRegion ,'.mat'])
    delete([folder_to_save '/regresspval_' inputRegion ,'.mat']);
    delete([folder_to_save '/regressreg_' inputRegion ,'.mat' ]);
     delete([folder_to_save '/regressCPD_' inputRegion ,'.mat' ]);
end



matpval = matfile([folder_to_save '/regresspval_' inputRegion ,'.mat'], 'Writable',true);
matcorr = matfile([folder_to_save '/regressreg_' inputRegion ,'.mat'], 'Writable',true);
 




for gg= 1:size(listpval,1)
    if isempty(regexp(listpval(gg).name, 'model', 'once'))
    pos=1;

    val2use =  listpval(gg).name;
    corr2use = ['corr', listpval(gg).name(5:end)];
%     CPD2use = ['CPD_', listpval(gg).name(5:end)];
    
    for region = 1:size(jaja,2)
        
        region2 = jaja{region};
        load (['regresspval_',region2,'.mat'], val2use)
        load (['regressreg_',region2,'.mat'], corr2use)


        cellsize = size(eval(val2use), 1);
        
        ppre(pos:pos-1 + cellsize/2,:) = eval([val2use,'(1:cellsize/2,:)']);
        cpre(pos:pos-1 + cellsize/2,:) = eval([corr2use, '(1:cellsize/2,:)']);
        ppost(pos:pos-1 + cellsize/2,:)= eval([val2use, '(1+cellsize/2:end,:)']);
        cpost(pos:pos-1 + cellsize/2,:)= eval([corr2use,'(1+cellsize/2:end,:)']);
        
        
        pos = pos + cellsize/2;
        
    end
    
    
    val2useVal =  [ppre; ppost];
    corr2useVal = [cpre; cpost];
%     CPD2useVal = [cpre; cpost];
    
    clear ppre  cpre ppost cpost
    
    
    matpval.(val2use) = val2useVal;
    matcorr.(corr2use) = corr2useVal;
%     matCPD.(corr2use) = CPD2useVal;

   clear pval* corr*
    end
end




%% This part only for certain analyses
%{

for gg= 1: size(listpval,1)
    if ~isempty(regexp(listpval(gg).name, 'model', 'once'))
    
    pos=1;

    val2use =  listpval(gg).name;

    for region = 1:size(jaja,2)
        
        region2 = jaja{region};
        load (['regresspval_',region2,'.mat'], val2use)

        cellsize = size(eval(val2use), 1);
        
        ppre(pos:pos-1 + cellsize/2,:) = eval([val2use,'(1:cellsize/2,:)']);

        ppost(pos:pos-1 + cellsize/2,:)= eval([val2use, '(1+cellsize/2:end,:)']);

        
        pos = pos + cellsize/2;
        
    end

   val2useVal =  [ppre; ppost];
     clear ppre  cpre ppost cpost
  
    matpval.(val2use) = val2useVal;
    
    clear pval* corr*
    end
end





matCPD = matfile([folder_to_save '/regressCPD_' inputRegion ,'.mat'], 'Writable',true);

clear CPD*
load (['regressCPD_PHC.mat'])
listCPD=whos('CPD*');

clear CPD*




for gg= 1: size(listCPD,1)
    
    
    pos=1;
    
    
    CPD2use = listCPD(gg).name;
    
    for region = 1:size(jaja,2)
        
        region2 = jaja{region};
        load (['regressCPD_',region2,'.mat'], CPD2use)
        
        
        
        cellsize = size(eval(CPD2use), 1);
        
        ppre(pos:pos-1 + cellsize/2,:) = eval([CPD2use,'(1:cellsize/2,:)']);
        
        ppost(pos:pos-1 + cellsize/2,:)= eval([CPD2use, '(1+cellsize/2:end,:)']);
        
        
        pos = pos + cellsize/2;
        
        
    end
    
    CPD2useVal = [ppre; ppost];
    clear ppre  cpre ppost cpost
    
    matCPD.(CPD2use) = CPD2useVal;
    
    clear pval* corr*
end
%}
%%



%% create new varaibles, with the same names and structure but pooled

% for gg= 1: size(listpval,1)
%     eval([listpval(gg).name, ' =  [ppre', num2str(gg), '; ppost', num2str(gg), '];'])
%     eval(['corr',listpval(gg).name(5:end), ' =  [cpre', num2str(gg), '; cpost', num2str(gg), '];'])
%     
% end


% 
% 
% h=0; g=0;
% while h<1
%     try
%         load([folder_to_save '/regresspval_' inputRegion ])
%         h=1;
%     catch
%         save ([folder_to_save '/regresspval_' inputRegion ],  'pval*', '-v7.3')
%         if g>0
%             warning(['error while saving the pvalues in region ', inputRegion])
%         end
%         g=g+1;
%     end
% end
% 
% h=0; g=0;
% while h<1
%     
%     try
%         load([folder_to_save '/regressreg_' inputRegion ])
%         h=1;
%     catch
%         save ([folder_to_save '/regressreg_' inputRegion ],  'corr*', '-v7.3')
%         if g>0
%             warning(['error while saving the corr in region ', inputRegion])
%         end
%         g=g+1;
%     end
% end
% 

     
%      
%      
%     save ([folder_to_save '/regressreg_' inputRegion ],  'corr*', '-v7.3', '-append')
%     h=0; g=0;
%     while h<1
%         try
%            save ([folder_to_save '/regresspval_' inputRegion ],  'pval*', '-v7.3', '-append')
%           save ([folder_to_save '/regresspval_' inputRegion ],  'pval*', '-v7.3', '-append')
%             h=1;
%         catch
%             save ([folder_to_save '/regresspval_' inputRegion ],  'pval*', '-v7.3', '-append')
%             if g>0
%                 warning(['error while saving the pvalues in region ', inputRegion])
%             end
%             g=g+1;
%         end
%     end
%     
%     h=0; g=0;
%     while h<1
%         
%         try
%             save ([folder_to_save '/regressreg_' inputRegion ],  'corr*', '-v7.3', '-append')
%             h=1;
%         catch
%             save ([folder_to_save '/regressreg_' inputRegion ],  'corr*', '-v7.3', '-append')
%             if g>0
%                 warning(['error while saving the corr in region ', inputRegion])
%             end
%             g=g+1;
%         end
%     end
    

