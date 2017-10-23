function [sbins,sbinsN]=throwoutbins(mybinsN,ntol,Nt)
nanbins=isnan(mybinsN);
% get the bins that have data points in them
[nonzbins,foo,ibins]=unique(mybinsN(~nanbins));
%count the number of observations in each bin
if length(nonzbins)==1
    Ncount=length(mybinsN);
else
    [Ncount,foo]=hist(mybinsN,nonzbins);
end
if sum(Ncount)+sum(nanbins)~=length(mybinsN)
    error('all of the observations are not accounted for in the surviving bins')
end
Ncount=Ncount';
percentbins=Ncount/Nt;
ikeep=percentbins>=ntol;
sbins=nonzbins(ikeep);
sbinsN=Ncount(ikeep);
end