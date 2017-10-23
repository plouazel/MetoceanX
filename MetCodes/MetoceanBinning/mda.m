function bindata=mda(obsdata,n,varargin)
% maximum dissimilarity algorithm
% Computer Aided Design of Experiments
% Author(s): R. W. Kennard and L. A. Stone
% Technometrics, Vol. 11, No. 1 (Feb., 1969), pp. 137-148
% http://www.libpls.net/publication/KS_1969.pdf
% INPUTS:
% input data that is N x p
% n is the design 
% OUTPUTS:
% bindata is a design set of n points chosen from the obsdata
if nargin<3
    Lk=2; % Euclidian
else
    Lk=varargin{1};
end
if nargin<4
    mydir=userpath;
    matdir=[mydir(1:end-1) filesep];
else
    matdir=varargin{2};
end
[N,p]=size(obsdata);
means=mean(obsdata,1);
stds=std(obsdata,[],1);
dmin=nan(N,1); imin=nan(N,1);
dmax=nan(N,1); imax=nan(N,1);
dstd=0.5;
for nstd=0:dstd:2
    iex = sum(obsdata > repmat(means,[N 1]) + repmat(nstd*stds,[N 1]) | obsdata < repmat(means,[N 1]) - repmat(nstd*stds,[N 1]),2);
    ikeep= iex >= p-1 ;
    if sum(ikeep) < .1*N
         break
    end
end

if sum(ikeep)==0 && nstd>0
    nstd=nstd-dstd;
end
% means+nstd*stds
% means-nstd*stds
% obsdata(4652,:)
% sqrt(sum((diff(obsdata([1589 13427],:),1,1)).^2))
% sqrt(sum((diff(obsdata([4652 9370],:),1,1)).^2))
% pause
% n2save=1000;
% ilast=0;
initfile=[matdir sprintf('Distmax%06d.mat',N)];
if exist(initfile,'file') && N > 50e3
    load(initfile);
else
    disp(sprintf('Calculating distances for %d observations',N))
    for ii=1:N
        irun=1;
        idata = obsdata(ii,:);
        iobsdata = obsdata;
        iobsdata(ii,:)=nan(1,p); % get rid of i'th observation
        if N>100
            % try to only look at a small sample
            if sum(idata > means + nstd*stds | idata < means - nstd*stds)<p-1
                irun=0;
            end
        end
        if irun
%             is=mod(ii-1,n2save)+1;
%             if N-ii<n2save && is==1 && ilast==0
%                 ilast=1;
%             elseif ilast>0
%                 ilast=0;
%             end
%             if is==1
%                 if ~ilast
%                     %initialize
%                 	Distdata(n2save).d=nan(N,1);
%                 else
%                     Distdata(N-ii).d=nan(N,1);
%                 end
%             end
            [dmin(ii),imin(ii),dPmax,dmax(ii),imax(ii)] = getmydist(idata,iobsdata,1,Lk);
            %disp(sprintf('Calculating distances for %d/%d observations',ii,N))
            
%            [dmin(ii),imin(ii),dPmax,dmax(ii),imax(ii),dall] = getmydist(idata,iobsdata,1,Lk);
            %Distdata(is).d=uint16(dall/sqrt(p)*2^16);
%             if is==n2save || ii==N
%                 disp(sprintf('Saving distances for %d/%d observations',ii,N))
%                 save([matdir sprintf('Distdata%04d.mat',round(ii/100))],'Distdata') 
%                 clear Distdata
%             end
        end
    end
    if N > 50e3
        save(initfile,'dmax','imax')
    end
end

bindata=nan(n,p);

%iunbin=1:N;
ibin=zeros(n,1);
iobsdata=obsdata;

for kk=2:n
   disp(sprintf('Running mda on %d/%d observations',kk,n))
   if kk==2
      [maxmax,iPstar]=nanmax(dmax); % should be 2
      
      bindata(1,:) = obsdata(iPstar,:);
      bindata(2,:) = obsdata(imax(iPstar),:);
      iobsdata(iPstar,:)=nan(1,p);
      iobsdata(imax(iPstar),:)=nan(1,p);
      ibin(1) = iPstar;
      ibin(2) = imax(iPstar);
   else
       dmin=nan(N,1); imin=nan(N,1);
       dmax=nan(N,1); imax=nan(N,1);
       %iunbin = iunbin(~ismember(iunbin,ibin(1:kk-1)));
       iobsdata(ibin(1:kk-1),:)=nan(kk-1,p);
       
       for ii=1:N
%            is=mod(ii-1,n2save)+1;
%            if is==1
%                if N-ii<n2save
%                    load([matdir sprintf('Distdata%04d.mat',round(N/100))]);
%                else
%                    sprintf('Distdata%04d.mat',ii+n2save-1)
%                    load([matdir sprintf('Distdata%04d.mat',(ii+n2save-1)/100)]);
%                end
%            end
%            ddata = Distdata(is).d*sqrt(p)/2^16;
           if ~ismember(ii,ibin)
               idata = obsdata(ii,:);
               [dmin(ii),imin(ii),dPmax,dmax(ii),imax(ii)] = getmydist(idata,bindata(1:kk-1,:),1,Lk);
               %dmin(ii) = min(ddata(ibin(1:kk-1)));
               %dmax(ii) = max(ddata(ibin(1:kk-1)));
           end
       end
       [maxmin,iPstar]=nanmax(dmin);
       bindata(kk,:) = obsdata(iPstar,:) ;
       ibin(kk) = iPstar;
   end
end
end