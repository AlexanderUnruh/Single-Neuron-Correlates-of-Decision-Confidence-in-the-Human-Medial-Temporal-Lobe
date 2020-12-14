function pool_all_regionsChange(whattotpool)


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
    delete([folder_to_save '/regressMFR_' inputRegion ,'.mat' ]);
end



matpval = matfile([folder_to_save '/regresspval_' inputRegion ,'.mat'], 'Writable',true);
matcorr = matfile([folder_to_save '/regressreg_' inputRegion ,'.mat'],  'Writable',true);
 




for gg= 1:size(listpval,1)
    if isempty(regexp(listpval(gg).name, 'model', 'once'))
	    pos=1;


    val2use =  listpval(gg).name;
    corr2use = ['corr', listpval(gg).name(5:end)];

	  
	    for region = 1:size(jaja,2)
		
		region2 = jaja{region};
		load (['regresspval_',region2,'.mat'], val2use)
		load (['regressreg_',region2,'.mat'], corr2use)


		cellsize = size(eval(val2use), 1);
		
		ppre(pos:pos-1 + cellsize,:) = eval(val2use);
		cpre(pos:pos-1 + cellsize,:) = eval(corr2use);

		
		pos = pos + cellsize;
		
	    end
	    
    
	    val2useVal =  [ppre];
	    corr2useVal = [cpre];
	

	    
	    clear ppre  cpre ppost cpost
	    
	    
	    matpval.(val2use) = val2useVal;
	    matcorr.(corr2use) = corr2useVal;


	   clear pval* corr*
    end
end


 

