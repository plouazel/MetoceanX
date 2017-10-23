function bindata=mda2(obsdata,n,bindata,idir,varargin)

if nargin<4
    Lk=2; % Euclidian
else
    Lk=varargin{1};
end
if nargin<5
    mydir=userpath;
    matdir=[mydir(1:end-1) filesep];
else
    matdir=varargin{2};
end
[N,p]=size(obsdata);
% initialize
dmin=nan(N,1); imin=nan(N,1);
dmax=nan(N,1); imax=nan(N,1);
dRold = nan(N,1);
nold = size(bindata,1);
nremain = N;
if n>nremain
    n = nremain;
end
bindata=[bindata; nan(n,p)];
kstart = max([ find(isnan(bindata(:,1)),1) 2]) ; % don't startat1
kfinish = size(bindata,1);
n2save=1000;
initfile=[matdir sprintf('Distmax_%dd%06dN.mat',p,N)];

if exist(initfile,'file') && N > 50e3
   disp(['Loading ' initfile])
   load(initfile);
else
    for kk=kstart:kfinish
       if kk==2
           disp(sprintf('Calculating distances for %d observations',N))
           for ii=1:N
            [dmin(ii),imin(ii),dPmax,dmax(ii),imax(ii)] = getmydist(obsdata(ii,:),obsdata,idir,1,Lk);
              %[dmax(ii,:),imax(ii,:)]=getNfurthest(obsdata(ii,:),obsdata,R,n,Lk);
              is=mod(ii-1,n2save)+1;

           end
       
         [maxmax,iPstar]=nanmax(dmax(:,1)); % should be 2 that have equal distance

          bindata(1,:) = obsdata(iPstar,:);
          bindata(2,:) = obsdata(imax(iPstar,1),:);
          obsdata(iPstar,:)=nan(1,p);
          obsdata(imax(iPstar),:)=nan(1,p);
          ibin(1) = iPstar;
          ibin(2) = imax(iPstar,1);
       
      
      
       else
           %dmax(ibin(1:kk-1),:)
           disp(sprintf('Running mda on %d/%d observations',kk-kstart+1,n))
           for ii=1:N
    %            if kk==3
    %                dRold(ii) = sum((obsdata(ii,:)-bindata(kk-2,:)).^Lk,2).^(1/Lk);
    %            end
               if kk==kstart || kk==3
                   % you've entered this statement for the first time
                   dRs=sum((repmat(obsdata(ii,:),[kk-2 1])-bindata(1:kk-2,:)).^Lk,2).^(1/Lk); % distance betwen arbitrary obs and all previous clusters
                   dRold(ii) = min(dRs);
               end
                dR=sum((obsdata(ii,:)-bindata(kk-1,:)).^Lk,2).^(1/Lk); % the distance between an arbitrary obs and the last observation that was put into the set
                dmin(ii) = min([dR dRold(ii)]);
           end
           [maxmin,iPstar] = max(dmin);
           bindata(kk,:) = obsdata(iPstar,:);
           ibin(kk) = iPstar;
           dRold = dmin;
       end
    end
    if N > 50e3
        save(initfile,'bindata')
    end
end
       
end
function [dmax,imax]=getNfurthest(idata,bindata,R,n,Lk)
% get n distances, which pass the distance threshold, 
% return indices of points which are far away
dmax=nan(1,n);
imax=nan(1,n);
[Nbins,Np]=size(bindata); %means
dallP=repmat(idata,[Nbins 1])-bindata; %Nbins x Np, distance of data point to all centroids of bins
dall=sum(dallP.^Lk,2).^(1/Lk); %N-D distance of data point, Euclidian distance
itop = dall > R; %threshold check
npass=sum(itop);
passind=find(itop);
[dsort,isort]=sort(dall(itop),'descend');
if npass<n
    % pad with nans
    
    dmax(1,1:npass) = dsort; imax(1,1:npass)=passind(isort);
else
    imax(1,:)=passind(isort(1:n));
    dmax(1,:)=dsort(1:n);
    
end

end