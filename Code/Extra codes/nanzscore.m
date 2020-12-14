function [xOut] =  nanzscore (xIn, type, dim)

if nargin ==1
    type =0;
    dim =1;
elseif nargin ==2
    
    dim =1;
end
   nanmean = nansum(xIn, dim)./ sum(~isnan(xIn), dim);
   xOut = (xIn-nanmean)./nanstd(xIn,type,dim);
    
    