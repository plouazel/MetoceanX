function Wts=getWeights(rawdata,maxs,offsets,iget,Fname)
iplot=0;
%[x1g,x2g,x3g,x4g,x5g]=ndgrid( myedges(~isnan(myedges(:,1)),1),myedges(~isnan(myedges(:,2)),2),myedges(~isnan(myedges(:,3)),3),myedges(~isnan(myedges(:,4)),4),myedges(~isnan(myedges(:,5)),5));
[nE,nP]=size(rawdata);

if iget
    F=getWtdata(Fname); % F better be in the same order as myedges, last column of F is the weights
    [nB,nP2]=size(F);
    if nP~=nP2-1
        error(['weight data needs to have ' num2str(nP) ' dimensions, not ' num2str(nP2-1) 'dimensions'])
    end
    data=[(F(:,1:end-1)-repmat(offsets,[nB 1]))./repmat(maxs,[nB 1]) F(:,end)]; %normalize centroids 
    %[x1,x2,x3,x4,x5]=ndgrid([0 1],[0 1], [0 1],[0 1],[0 1]); %set 1-D weights?
    %WTF to do outside of the range of the centroids? how to extrapolate?
    data=[data ;ones(1,nP) 10];
    %[x1g,x2g,x3g,x4g,x5g]=ndgrid(data(:,1),data(:,2),data(:,3),data(:,4),data(:,5));
    %F=scatteredInterpolant({data(:,1),data(d:,2),data(:,3),data(:,4),data(:,5) },data(:,6));
    rbfCoeff=rbfcreate(data(:,1:end-1)',data(:,end)','RBFfunction','linear');
    Wts=rbfinterp(rawdata',rbfCoeff);
    Wts=Wts';
    Wts(Wts<0)=ones(sum(Wts<0),1);
else
    Wts=ones(nE,1);
end
%Wts=max([Wts ones(nE,1)],[],2);
if nP==2 && iplot
    vizWts(data,rawdata,maxs,offsets,Wts);
end
end

function vizWts(data,rawdata,maxs,offsets,Wts)

%re-normalize
% how the heck do we visualize this? only in 2D
[nB,nP1]=size(data);
newdata=[data(:,1:2).*repmat(maxs,[nB 1])+repmat(offsets,[nB 1]) data(:,end)];
[nB2,nP2]=size(rawdata);
alldata=rawdata.*repmat(maxs,[nB2 1])+repmat(offsets,[nB2 1]);

[xq,yq] = meshgrid([0:.01:1]*maxs(1)+offsets(1), [0:.01:1]*maxs(2)+offsets(2));
vq = griddata(alldata(:,1),alldata(:,2),Wts,xq,yq);
    figure(1)
    mesh(xq,yq,vq);
    hold on
    plot3(newdata(:,1),newdata(:,2),newdata(:,3),'ko',alldata(:,1),alldata(:,2),Wts,'r.');
    h = gca;
    hold off

end
function [F]=getWtdata(fname)
F=csvread(fname);
[nT,nP2]=size(F);
if nP2-1==2
    [foo,iF]=unique(F(:,1));
    F=F(iF,:);
end
% [nE,nP]=size(myedges);
% fid=fopen(fname,'r');
% hstr1=repmat('%s,',[1 nP+1]);
% headers=fgetl(fid,hstr1);
% dstr1=repmat('%3.3f,',[1 nP+1]);
% ii=1;
% while ~feof(fid)
%    F(:,ii)=fgetl(fid,dstr1) ;
%    ii=ii+1;
% end
% %F=textscan(fid,dstr1,'Delimiter', ',');
end