function corrdim = corrdim(mpimage,x,y)

a = find(mpimage(:)~=0); 

xloc = x(a); 
yloc = y(a); 

X = [xloc yloc]'; 

D = sqdist(X,X); 
D = D(D>0); 

minval = min(D(:)); 
eps = logspace(log10(minval),log10(max(D(:))),50); 

%%
for i = 1:length(eps)
    
    C(i) = sum(sum(D < eps(i)))/numel(D); 
    
end

T = polyfit(log(eps),log(C),1); 

corrdim = T(1); 



