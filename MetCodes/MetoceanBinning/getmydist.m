function [dmin,imin,dPmax,varargout]=getmydist(idata,bindata,idir,wdata,Lk,varargin)%wdata,Lk)
% OUTPUTS:
% d = minimum distance between current bin and all other centroids
% ibin = the index number of this minimum bin (using the vector bindata)
% dPmax = the max single dimension (weighted) distance between this
% centroid and the current observation (to use for distance checking)

% INPUTS
% idata is the n-Dim centroid of the current bin
% bindata is the means and stdevs of all bins
% wdata = scalar representing the 'weight' function at the centroid of the bin in question
% dF = dimensions to not allow observed data to be re-binned in
if nargin<4
    wdata=1;
end
if nargin<5
    Lk = 2; % Euclidian
end
if ~isempty(bindata)
    [Nbins,Np]=size(bindata); %means
    %Np=np2/2;
    %centroids=bindata(:,1:Np);
    wvec=wdata*ones(1,Np); % weight all dimensions the same
    dallP=repmat(wvec,[Nbins 1]).*abs((repmat(idata,[Nbins 1])-bindata)); %Nbins x Np, distance of data point to all centroids of bins
    if sum(dallP<0)
        error('your weights are negative!')
    end
    dall=sum(dallP.^Lk,2).^(1/Lk); %N-D distance of data point, Euclidian distance
    [dmin,imin] = nanmin(dall,1); %find the index of the bin that minimizes the Euclidian distance between the current observation and all of the centroids 
    [dmax,imax] = nanmax(dall,1); %find the index of the bin that maximizes the Euclidian distance between the current observation and all of the centroids 
    dPmax=max(dallP(imin,:),[],2); % Chebyshev distance of the bin that minimizes the Euclidian distance
else
    dmin=inf;
    imin=[];
    dPmax=inf;
    dmax=inf;
    imax=[];
end
varargout{1} = dmax;
varargout{2} = imax;
varargout{3} = dall;
end