function Dat=whatsMyState(IPTfile,varargin)
if nargin>0
    iptname=IPTfile;
else
    iptname=getIPTname(mfilename);
end
run(iptname);

datstr=[dattype '_{' casestr '}'];

switch casestr
    case 'SSS'
        iSSS=met.all.vspdhh<=VspdRated+1 & met.all.vspdhh>=VspdRated-1;
        Ryrs = [1 10 50];
    case 'ESS'
        iSSS=ones(length(met.all.time),1);
        Ryrs = [1];
end
noutliers=0;
Dat=cell(length(dirbins)-1,3);
for ii=1:length(dirbins)-1
    th1=mod(dirbins(ii),360);
    if ii==length(dirbins)
        th1=dirbins(1)-th0; th2=dirbins(ii)-th0+.01;
        titlestr='Omni';
        if strcmp(casestr,'SSS')
            noutliers=0;
        end
    else
        th2=dirbins(ii+1);
        titlestr=num2str(mean([dirbins(ii) th2]));
    end
    if dirbins(ii)<0
        iDir = met.all.vdir>=th1 | met.all.vdir<th2;
    else
        iDir = met.all.vdir>=th1 & met.all.vdir<th2;
    end
    iUse= iSSS & iDir;
    dat=met.all.(dattype);
    if exist('iUse','var')
        datTp=met.all.Tp(iUse);
    else
        datTp=met.all.Tp;
    end
    rDat=getRyr(met.all.time,dat,datstr,iUse,noutliers,Ryrs,1,iplot);
    hf=gcf;
    set(hf,'name',[casestr '_Hs_' titlestr])
    title(['WFM ' datstr num2str(th1) '< \theta_{wave} <' num2str(th2)])
    h1=gca;
    set(h1,'xlim',[0 ceil(max(met.all.Hs)*1.4)])
    if strcmp(dattype,'Hs')
        xe=[ones(length(dat(iUse)),1) dat(iUse)];
        b=xe\datTp; % this is some linear regression
        maxT=b(1)+rDat*b(2);
        iDat={titlestr, rDat(end) ,maxT(end)};
    end
    Dat(ii,:) = iDat; %(ii-1)*ndat+1:(ii-1)*ndat+ndat
end

end