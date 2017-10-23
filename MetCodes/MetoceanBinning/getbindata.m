function [bindata,binN]=getbindata(sbins,allbins,rawdata)
Nbins=length(sbins);
[Nt,Np]=size(rawdata);
bindata=nan(Nbins,2*Np); % means, stdevs 
binN=nan(Nbins,1); % means, stdevs 
for ii=1:Nbins
    idata=allbins==sbins(ii);
    bindata(ii,:)=[ mean(rawdata(idata,:),1) std(rawdata(idata,:),[],1) ];
    binN(ii)=sum(idata); % same vector as sbinsN
end
end