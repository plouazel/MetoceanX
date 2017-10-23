function myedges=getfirstedges(nP,nE)
eCount=zeros(1,nP);
for ii=1:nP   
    if isfield(nE,'count')   
        Edges(ii).loc=getNptsfromPDF(nE(ii).xs,nE(ii).ys,round(nE(ii).count));
        eCount(ii)=length(Edges(ii).loc)-1;
    else
       Edges(ii).loc=nE(ii).xs;
       eCount(ii)=length(Edges(ii).loc)-1;
    end
end
nB=max(eCount);
myedges=nan(nB+1,nP);
for ii=1:nP
    myedges(1:eCount(ii)+1,ii)=Edges(ii).loc';
end
%myedges=[e1,e2,e3,e4,e5];

if min(max(myedges,[],1))<1
    ni=find(max(myedges,[],1)<1);
    for ii=1:length(ni)
        warning(sprintf('1-D maxs do not encompass all of the data for Parameter = %d',ni(ii)))
    end
end
if sum(myedges(1,:)<0)>0
    ni=find(mins<myedges(1,:));
    for ii=1:length(ni)
        warning(sprintf('1-D mins do not encompass all of the data for Parameter = %d',ni(ii)))
    end
end
end