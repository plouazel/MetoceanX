function varargout = write1dScatterTable(IPTfile)
run(IPTfile);

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
dw = mean(diff(wBin));
wPn=(wP-dw/2); % wP is in OrcaFlex convention

if strcmp(MetType,'Wave')
    datamat=[ mod(-met.(metfield).wdir(1:nS:end),360),  met.(metfield).Tp(1:nS:end), met.(metfield).Hs(1:nS:end)]; % wdir is TOWARDS going +CW (oceanographic convention) -> Orca
    myvars = {'wdir','Tp','Hs'};
    myunits = {'deg','s','m'};
    myedges=getmyedges(wPn,wBin,Tpbin,Hsbin); % Np x Nbins size matrix of the edges of the initial bins
elseif strcmp(MetType,'Wind')
    datamat=[ mod(180-met.all.vdir(1:nS:end),360),  met.all.vspd(1:nS:end)]; % vdir is FROM going +CW (meteorological convention) -> Orca
    myvars = {'vdir','vspd'};
    myunits = {'deg','m/s'};
    myedges=getmyedges(wPn,wBin,Vspdbin); % Np x Nbins size matrix of the edges of the initial bins
elseif strcmp(MetType,'Current')
        datamat=[ mod(180- met.all.Udir(1:nS:end),360),  met.all.Uspd(1:nS:end)]; % vdir is FROM going +CW (meteorological convention) -> Orca
        myvars = {'udir','uspd'};
        myunits = {'deg','m/s'};
        myedges=getmyedges(wPn,wBin,Uspdbin); % Np x Nbins size matrix of the edges of the initial bins
end
DirStr='(TOWARDS, +CCW)';
strd='%2.6f,';
[ndata,nD] = size(datamat);
idir = [1 zeros(1,nD-1)];
nE=zeros(nD,1);

for kk=1:nD
    nE(kk) = length(myedges(~isnan(myedges(:,kk)),kk))-1;
end

allbinsN=findmybin(datamat,myedges,idir); % [Nt x 1] array where each row is a scalar of the index of the corresponding bin for each parameter (in columns)

nM=prod(nE(2:end));
if nD>2
    mymat = zeros([nE(2:end)' nE(1)]);
    [colstr,rowstr]=getrcstrs(myunits(2:end),myedges(1:nE(2),2),myedges(1:nE(3),3));
else
     mymat = zeros(nE(2:end)',1, nE(1));
     [colstr,rowstr]=getrcstrs(myunits(2:end),myedges(1:nE(2),2));
end

Scatter.Project = ProjectName;
Scatter.MetData = metname;

%% write to a table
fmt=[repmat(strd,[1,size(mymat,2)+1]) ' \n'];

fid=fopen(tname,'w+');
if exist('strh','var')
    for jj=1:length(strh)
        fprintf(fid,'%s \n',strh{jj});
    end
end


%% Grab the data
for ii=1:nE(1)
    % get the observations in between the current bin
   cbinI= allbinsN>=1+nM*(ii-1) & allbinsN <= nM*ii;
   [counts,places]=hist(allbinsN(cbinI),unique(allbinsN(cbinI)));
   % write dir
   
   fprintf(fid,'\n');
   fprintf(fid,'%s %d to %d %s',[DirStr MetType ' Heading: '],[round(myedges(ii,1) ) round(myedges(ii+1,1))],'deg');
   fprintf(fid,'\n');
   fprintf(fid,'%s',colstr{:});
   fprintf(fid,'\n');
      %fill in the current matrix
   for jj=1:length(counts)
       elem=places(jj)-nM*(ii-1);%mod(places(jj),nM);%;  
       if nD>2
           nC = nE(3);
       else
           nC = 1;
       end
       rowj=ceil(elem/nC);
       
        colj=elem-nC*(rowj-1);
%        else
%            colj=1;
%        end

       if counts(jj)
           mymat(rowj,colj,ii)=counts(jj);
       end
       
   end
   % TABLE
   for kk=1:size(mymat,1)
      crow=squeeze(mymat(kk,:,ii))*100/ndata;
      tdata=[rowstr{kk} sprintf(fmt,[crow sum(crow)])];
      fprintf(fid,'%s',tdata);
      %WaveScatter(ii).Prob(kk,:) = crow;
   end
   % SUMTABLE
   cmat=squeeze(mymat(:,:,ii))*100/ndata;
   colsum=sum(cmat,1);
   tdata=['Sum, ' sprintf(fmt,[colsum sum(colsum)]) ];
   fprintf(fid,'%s',tdata);
   
   % Store
   Scatter.dir1(ii) = myedges(ii,1) ;
   Scatter.dir2(ii) = myedges(ii+1,1) ;
   %Scatter.wdir(ii) = mean( [Scatter.dir1(ii) Scatter.dir2(ii)]);
       %if nE(
%    for kk=1:nT
%       crow=squeeze(mymat(kk,:,ii))*100/ndata;
%       tdata=[rowstr{kk} sprintf(fmt,[crow sum(crow)])];
%       fprintf(fid,'%s',tdata);
%       probmat(kk,:) = crow;
%       %WaveScatter(ii).Prob(kk,:) = crow;
%    end    
   Scatter.Prob{ii} = squeeze(mymat(:,:,ii))*100/ndata;
   
end

summat=squeeze(sum(mymat,3));
fprintf(fid,'\n');
fprintf(fid,'%s','All Headings: ');
fprintf(fid,'\n');
fprintf(fid,'%s',colstr{:});
fprintf(fid,'\n');
for kk=1:size(mymat,1)
  crow=squeeze(summat(kk,:))*100/ndata;
  tdata=[rowstr{kk} sprintf(fmt,[crow sum(crow)])];
  fprintf(fid,'%s',tdata);
end
cmat=summat*100/ndata;
colsum=sum(cmat,1);
tdata=['Sum, ' sprintf(fmt,[colsum sum(colsum)]) ];
fprintf(fid,'%s',tdata);
fclose(fid);
%Scatter.E = nan(nD-1,max(nE(2:end)));

for pp=1:nD
 %  Scatter.E(pp,1:nE(pp)) = [mean( [myedges(1:nE(pp)-1,pp) myedges(2:nE(pp),pp) ],2); myedges(nE(pp),pp)]'; % this isn't really right.. what should the last bin be?
   if pp>1 % should use idir...
        Scatter.(myvars{pp}) =[mean( [myedges(1:nE(pp)-1,pp) myedges(2:nE(pp),pp) ],2); myedges(nE(pp),pp)]'; 
   else
       Scatter.(myvars{pp}) =[mean( [myedges(1:nE(pp),pp) myedges(2:nE(pp)+1,pp) ],2)];
   end
end
% matfile
save(matname,'Scatter')
varargout{1} = matname;

end

function myedges=getmyedges(wPn,bin1,bin2,bin3,varargin)
% bin1 must be directions!! -> should use idir?
if wPn>180
    wPn=wPn-360;
elseif wPn<-180
    wPn = wPn+360;
end

if nargin<4
    nP=2;
    nBs = [length(bin1) length(bin2)+1];
elseif nargin ==4;
    nP=3;
    nBs = [length(bin1) length(bin2)+1 length(bin3)+1];
end
nBmax = max(nBs);

myedges = nan(nBmax,nP);
for pp=1:nP
    eval(['pbins = bin' num2str(pp) ';']) % ugly!!!
    if pp>1
        myedges(1:nBs(pp)-1,pp) = pbins;
        myedges(nBs(pp),pp) = inf;
    else
        % directional! -> take into account platform heading
       pbins = bin1 + wPn;% OrcaFlex conventions
       myedges(1:nBs(pp),pp) = pbins; 
    end
end

end

function [colstr,rowstr]=getrcstrs(units,e1,e2,varargin)
n1=length(e1);


for jj=1:n1
    if jj<n1
        rowstr{jj}=sprintf('%1.1f-%1.1f(%s),',e1(jj),e1(jj+1),units{1});
    else
        rowstr{jj}=sprintf('%1.1f +,', e1(jj));
    end
end
if nargin>2
n2=length(e2);
for jj=1:n2+2
    if jj<=n2 && jj>1
        colstr{jj}=sprintf('%1.1f-%1.1f (%s),',e2(jj-1),e2(jj),units{2});
    elseif jj==1
        colstr{jj}=sprintf('%s','foo,');
    elseif jj==n2+1
        colstr{jj}=sprintf('%1.1f +,', e2(jj-1));
    elseif jj==n2+2
        colstr{jj}=sprintf('%s,', 'Sum');
    end
end
else
    colstr{1}='';
end
end