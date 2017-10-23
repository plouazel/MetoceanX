function [allbinsN]=findmybin(rawdata,myedges,idir)
[Nb,Np1]=size(myedges);
[Nt,Np2]=size(rawdata);
if Np1~=Np2 || Np1 ~= length(idir)
    error('number of parameters in data and edges do not match up')
end
Np=Np1;
allbinsP=nan(Nt,Np);
nbinsP=sum(~isnan(myedges),1)-1; % 1 x nP array - remove nans
nsub=ones(1,Np); nsub(Np)=0;
allbinsN=zeros(Nt,1);
for jj=1:Np
    je=myedges(~isnan(myedges(:,jj)),jj );
    Njb=length(je)-1;
    %add a little bit to je(end) so that we can use < (and not <=) on the
    %last 

    dataj=repmat(rawdata(:,jj),[1 Njb]);
    if idir(jj)
        jeL=repmat(je(1:end-1)',[Nt 1]);
        jeH=repmat(je(2:end)',[Nt 1]);  
        circ=2*pi;
        if max(abs(rawdata(:,jj)))>30
            circ=360;
        end
        logsC = zeros(size(dataj));
        for nn=-2:2
            logsL=dataj>=jeL + nn*circ;
            logsH=dataj<jeH + nn*circ;
            logsC = logsC | (logsL & logsH);
        end
    else
        je(end)=je(end)*1.01;
        jeL=repmat(je(1:end-1)',[Nt 1]);
        jeH=repmat(je(2:end)',[Nt 1]);  
        logsL=dataj>=jeL;
        logsH=dataj<jeH;
        logsC = logsL & logsH;
        if sum(logsL,2)==0
            aL=find(sum(logsL,2)==0);
            error(sprintf('Min of the bins does not encompass all of the data for Parameter %d at Observation %d',jj,aL))
        elseif sum(logsH,2)==0
            aH=find(sum(logsH,2)==0);
            error(sprintf('Max of the bins does not encompass all of the data for Parameter %d at Observation %d',jj,aL))
        end
        
    end
   
    [rlogs,clogs,unos]=find(logsC); 
    % for some reason "find" mixes up the order of observations 
    %resort
    [rowsa,ia]=sort(rlogs);
    binsN=clogs(ia);
    allbinsP(:,jj)=binsN;

    nplace=prod(nbinsP(1,jj+1:end),2);
    nHi=(allbinsP(:,jj)-nsub(jj)).*repmat(nplace,[Nt 1]);
    allbinsN=allbinsN+nHi;
end  
end