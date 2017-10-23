function [rdata,binN,q]=getbintable(sbins,allbins,rawdata,norms,offsets)
%[Wdfit,Hsfit,Tpfit]=getWindSeaFits(met); %Winddir-> Chopdir, Windspeed->ChopHeight, Chopheight-> ChopPer

Nbins=length(sbins);
[Nt,Np]=size(rawdata);
bindata=nan(Nbins,2*Np); % means, 
rdata=nan(Nbins,Np); % means , un-normalized
binN=nan(Nbins,1); % means, 
jpD=nan(Nt,Np);
k=2;

for ii=1:Nbins
    idata=allbins==sbins(ii);
    nibin=sum(idata); % number of data points in bin same vector as sbinsN
    bindata(ii,:)=[ mean(rawdata(idata,:),1) std(rawdata(idata,:),[],1) ];
    binN(ii)=nibin;
    dataun=(rawdata(idata,:).*repmat(norms,[nibin 1]))+repmat(offsets,[nibin 1]);
    rdata(ii,:)=mean(dataun,1);
    jpD(logical(idata),:)=rawdata(idata,:)-repmat( mean(rawdata(idata,:),1),[nibin 1]);
    %wHs=feval(Hsfit,rdata(ii,5)); % based on wind speed
    %wTp=feval(Tpfit,rdata(ii,5)); % based on wind speed
    %wD=round(feval(Wdfit,rdata(ii,4))); % based on wind direction
    %tdata=[ii binN(ii) 10000*binN(ii)/Nt rdata(ii,1) rdata(ii,2) round(rdata(ii,3)) round(rdata(ii,4)) rdata(ii,5) wHs wTp wD];
%     if iwrite
%     	fprintf(fid,fmt,tdata);
%     end
end
q=100*(1-(1/Nt)*sum( ( (1./Np) *sum(jpD.^k,2)).^(1/k),1));
end