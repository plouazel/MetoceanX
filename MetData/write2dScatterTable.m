function varargout=write2dScatterTable(IPTfile,varargin)
if nargin<1
    IPTfile = uigetfile('*ipt.m','Select input file');
end
run(IPTfile);

if ~exist(tdir,'dir')
    mkdir(tdir)
end

met=readhindcast(metname,metdir,0);

%define sampling rate (to reduce data size)
nS=1;
if ~exist('WaveType','var')
    WaveType='combo';
end
if strcmp(WaveType,'combo')
    metfield = 'all';
else
    metfield = WaveType;
end
datamat=[ met.(metfield).wdir(1:nS:end),  met.(metfield).Tp(1:nS:end), met.(metfield).Hs(1:nS:end)]; % wdir is TOWARDS going +CW (oceanographic convention)

idir = [1 0 0];
[ndata,nP]=size(datamat);
rawdata=datamat;%./repmat(maxs,[ndata 1]);
%change the wave bins so it goes around the platform
wPn=(wP-wBin/2); % wP is in OrcaFlex convention
mywave=mod(-rawdata(:,1),360); % convert to OrcaFlex/WAMIT convention, TOWARDS going +CCW
%mywave(mywave>=360+wPn)=mywave(mywave>=360+wPn)-360;
mydata=[mywave rawdata(:,2:end)];
%[rawdata,maxs,offsets]=preparedata4bins(mydata); %nobservation x Np matrix of all the time series of the data
% somehow get the first bins
myedges=getmyedges(wBin,wPn,Hsbin,Tpbin); % Np x Nbins size matrix of the edges of the initial bins
[nE,foo]=size(myedges);

allbinsN=findmybin(mydata,myedges,idir); % [Nt x 1] array where each row is a scalar of the index of the corresponding bin for each parameter (in columns)

%each wave direction bin has nE x nE Hs, Tp bins 
nE=zeros(size(datamat,2),1);
for kk=1:nE
    nE(kk) = length(myedges(~isnan(myedges(:,kk)),kk))-1;
end
 nW=length(myedges(~isnan(myedges(:,1)),1))-1;
 nT=length(myedges(~isnan(myedges(:,2)),2))-1;
 nH=length(myedges(~isnan(myedges(:,3)),3))-1;
mymat=zeros(nT,nH,nW);
%mymat = zeros(nE(2:end)',nE(3));
%nM=nT*nH;
nM=prod(nE(2:end));

%% write to a table
fid=fopen(tname,'w+');
if exist('strh','var')
    for jj=1:length(strh)
        fprintf(fid,'%s \n',strh{jj});
    end
end
strd='%2.6f,';
fmt=[repmat(strd,[1,nE(3)+1]) ' \n'];
[colstr,rowstr]=getrcstrs(myedges(1:nE(2),2),myedges(1:nE(3),3));
%% Write to a matfile
WaveScatter.Project = ProjectName;
WaveScatter.MetData = metname;
WaveScatter.Prob{nW}=nan(nE(2:end)');

%% fill in the tables
for ii=1:nE(1)
   cbinI= allbinsN>=1+nM*(ii-1) & allbinsN <= nM*ii;
   [counts,places]=hist(allbinsN(cbinI),unique(allbinsN(cbinI)));
   %fill in the current matrix
   for jj=1:length(counts)
       elem=places(jj)-nM*(ii-1);%mod(places(jj),nM);%;       
       rowj=ceil(elem/nH);
       colj=elem-nH*(rowj-1);
       if counts(jj)
           mymat(rowj,colj,ii)=counts(jj);
       end
   end
   fprintf(fid,'\n');
   fprintf(fid,'%s %d to %d %s','OrcaFlex (TOWARDS, +CCW) Wave Heading: ',[round(myedges(ii,1) ) round(myedges(ii+1,1))],'deg');
   fprintf(fid,'\n');
   fprintf(fid,'%s',colstr{:});
   fprintf(fid,'\n');
   
   %matfile
   %HeadName = strrep(   sprintf('%1.2f',mean([(myedges(ii,1) ) (myedges(ii+1,1))])),'.','_');
   WaveScatter.Wdir1(ii) = myedges(ii,1) ;
   WaveScatter.Wdir2(ii) = myedges(ii+1,1) ;
   WaveScatter.Wdir(ii) = mean( [WaveScatter.Wdir1(ii) WaveScatter.Wdir2(ii)]);
   WaveScatter.Hs = [mean( [Hsbin(1:end-1); Hsbin(2:end) ],1) Hsbin(end)]; % this isn't really right.. what should the last bin be?
   WaveScatter.Tp = [mean( [Tpbin(1:end-1); Tpbin(2:end) ],1) Tpbin(end)];
   probmat = zeros(nT,nH);
   for kk=1:nT
      crow=squeeze(mymat(kk,:,ii))*100/ndata;
      tdata=[rowstr{kk} sprintf(fmt,[crow sum(crow)])];
      fprintf(fid,'%s',tdata);
      probmat(kk,:) = crow;
      %WaveScatter(ii).Prob(kk,:) = crow;
   end
   WaveScatter.Prob{ii} = probmat;
   cmat=squeeze(mymat(:,:,ii))*100/ndata;
   colsum=sum(cmat,1);
   tdata=['Sum, ' sprintf(fmt,[colsum sum(colsum)]) ];
   fprintf(fid,'%s',tdata);
%    sum(cbinI)
%    sum(counts)
%    datamat(cbinI,1)
end
%sum table
summat=squeeze(sum(mymat,3));
fprintf(fid,'\n');
fprintf(fid,'%s','All Headings: ');
fprintf(fid,'\n');
fprintf(fid,'%s',colstr{:});
fprintf(fid,'\n');
for kk=1:nT
  crow=squeeze(summat(kk,:))*100/ndata;
  tdata=[rowstr{kk} sprintf(fmt,[crow sum(crow)])];
  fprintf(fid,'%s',tdata);
end
cmat=summat*100/ndata;
colsum=sum(cmat,1);
tdata=['Sum, ' sprintf(fmt,[colsum sum(colsum)]) ];
fprintf(fid,'%s',tdata);
fclose(fid);

% matfile
save(matname,'WaveScatter')
varargout{1} = matname;
end
function [colstr,rowstr]=getrcstrs(eT,eH)
nT=length(eT);
nH=length(eH);
for jj=1:nH+2
    if jj<=nH && jj>1
        colstr{jj}=sprintf('%1.1f-%1.1f (m),',eH(jj-1),eH(jj));
    elseif jj==1
        colstr{jj}=sprintf('%s','foo,');
    elseif jj==nH+1
        colstr{jj}=sprintf('%d +,', eH(jj-1));
    elseif jj==nH+2
        colstr{jj}=sprintf('%s,', 'Sum');
    end
end
for jj=1:nT
    if jj<nT
        rowstr{jj}=sprintf('%d-%d (s),',eT(jj),eT(jj+1));
    else
        rowstr{jj}=sprintf('%d +,', eT(jj));
    end
end
end


function myedges=getmyedges(wBin,wPn,Hsbin,Tpbin)
nP=3;%1D edges
Hsedge=[Hsbin,inf];%./maxs(3);
nHs=length(Hsedge)-1;
Tpedge=[Tpbin,inf];%./maxs(2);
nTp=length(Tpedge)-1;
if wPn>180
    wPn=wPn-360;
elseif wPn<-180
    wPn = wPn+360;
end
Wdedge=[wPn:wBin:360+wPn];% OrcaFlex conventions
nWd=length(Wdedge)-1;

%Nbins=nHs*nTp*nWd*nVd*nVsp;
%pad with nans to form an array-> they will be removed lates
%turn into N-D edges
padN=max([nWd nTp nHs])*ones(1,nP)-[nWd nTp nHs];

e3=[Hsedge' ; nan(padN(3),1)];
e2=[Tpedge' ; nan(padN(2),1)];
e1=[Wdedge' ; nan(padN(1),1)];
myedges=[e1,e2,e3];
end

