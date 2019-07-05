function Rsq = comp_rsq(data,model)

SST = nansum((data - nanmean(data)).^2);
SSres = nansum((data - model).^2);
Rsq = 1 - SSres/SST; 

end