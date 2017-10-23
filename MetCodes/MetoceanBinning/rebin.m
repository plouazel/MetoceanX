function [newbins,dPs,varargout]=rebin(sbins,allbins,rawdata,bindata,dtol,Wts,niter,dF,plt,myedges,irerebin,dmax,varargin) %sbins,allbinsN,rawdata,bindata,dtol,W
% take a set of bins (sbins), find out which observations have not been
% binned
% then define some distance function and find the 'nearest' bin centroid to
% that unbinned observation
% if there is a 1D that is too large, don't rebin 
% Lk is the type of distance , Lk=2 is Euclidian distance
iEdge=1;
if nargin<10
    iEdge=0;
end
%irerebin=1; % see if a new bin is closer to an observation, than when it had previously passed a tolerance test.
if nargin<11
     % you're running the old algorithm
     irerebin=0;
end
if nargin<12
    dmax=0;
    dnewmax= 0;
end
Lk=2;
%clrs={[1 0 0],[1 1 0],[.6 .6 0],[0 1 1],[0 1 0],[0 .6 .6],[0 0 1], [1 0 1],[.6 0 .6],};
clrs={[1 0 0],[1 1 0],[0 1 1],[0 1 0],[0 0 1], [1 0 1],};
clrs1=cellfun(@(x) x.*[.65 .8 .5],clrs,'un',0);
clrs2=cellfun(@(x) x.*[.5 .65 .8],clrs,'un',0);
clrs3=cellfun(@(x) x.*[.8 .5 .65],clrs,'un',0);
colorstr={clrs{:},clrs2{:},clrs3{:},clrs1{:}};
Nbins=length(sbins);
[Nt,Np]=size(rawdata);
lrebins=~ismember(allbins,sbins); %logical vector, length of all observations, of the bins who have failed the population test
Nrebins=sum(lrebins); % number of bins that have failed population test
irebins=find(lrebins); % vector containing indices of the bins that have failed population tests
newbins=allbins; % start the new vector of bins as the ol done
dPs=zeros(Nt,1);
%% FUNNY PLOTTING INTERLUDE
iplot= isfield(plt,'tmpfolder') & isfield(plt,'tofolder');
if iplot
   if niter==1
        basename=sprintf('Unweighted Bins in 2D: Iteration %d',niter-1);
       h0=figure('name',basename);
       clf
        ix=1;
        iy=2;
       plot(rawdata(:,ix),rawdata(:,iy),'k.')
       %title(sprintf('Number of Iterations: %d',niter))
       xlabel('normalized x-distance [-]')
       ylabel('normalized y-distance [-]')
       if ~isempty(myedges)   
            set(gca,'xtick',myedges(~isnan(myedges(:,ix)),ix ),'ytick',myedges(~isnan(myedges(:,iy)),iy ));
            set(gca,'xticklabel',round(100*myedges(~isnan(myedges(:,ix)),ix ))/100,'yticklabel',round(100*myedges(~isnan(myedges(:,iy)),iy ))/100);
       end
        grid on
        SaveAllFig(plt.tofolder)
   end
   
       basename = sprintf('Unweighted Bins in 2D: Iteration %da',niter);
       %addpath('/Users/samkanner/Documents/MATLAB/handy_MATLAB_codes/')
       h1=figure('name',basename);
       clf
        ix=1;
        iy=2;
       plot(rawdata(:,ix),rawdata(:,iy),'k.')
       %title(sprintf('Number of Iterations: %d',niter))
       xlabel('normalized x-distance [-]')
       ylabel('normalized y-distance [-]')
       if ~isempty(myedges)   
            set(gca,'xtick',myedges(~isnan(myedges(:,ix)),ix ),'ytick',myedges(~isnan(myedges(:,iy)),iy ));
            set(gca,'xticklabel',round(100*myedges(~isnan(myedges(:,ix)),ix ))/100,'yticklabel',round(100*myedges(~isnan(myedges(:,iy)),iy ))/100);
       end
        grid on
       
        %exportfigure(basename,plt.tofolder,plt.tmpfolder);
        hold on
   
   
   if min(min(Wts))<1
       [xq,yq] = meshgrid([0:.01:1], [0:.01:1]);
        vq = griddata(rawdata(:,ix),rawdata(:,iy),Wts,xq,yq);
        v=[.5:1:10.5];
       [C,h]=contour(xq,yq,vq,v);
       clabel(C,h)
   end
   for ii=1:Nbins
        idata=find(allbins==sbins(ii));
        plot(rawdata(idata,ix),rawdata(idata,iy),'.','markerfacecolor',colorstr{ii},'markeredgecolor',colorstr{ii})
        hLine=scatter(bindata(ii,ix),bindata(ii,iy),200,'o','markerfacecolor',.8*colorstr{ii},'markeredgecolor','k');%,'markersize',24);
        %
        ndata=length(idata);
        ndig=floor(log10(ndata))+1;
        text(bindata(ii,ix)-ndig*.004,bindata(ii,iy),sprintf('%d',ndata),'fontsize',6,'color','w')
        % draw the boundary around the cluster
        hold on
        if niter==1 && ~isinf(dtol)
            rectangle('position',[bindata(ii,ix)-dtol,bindata(ii,iy)-dtol,dtol*2,dtol*2],'edgecolor','r','linestyle','--','facecolor','none')
        end
%         if isempty(dF)
%             text(bindata(ii,ix)-ndig*.004,bindata(ii,iy)-ndig*.004,sprintf('Cluster %d',ii),'fontsize',6,'color','k')
%         end
        
   end
    xlim([0 1])
    ylim([0 1])
    axis equal
    SaveAllFig(plt.tofolder)
    basename=sprintf('Unweighted Bins in 2D: Iteration %d',niter);
    set(0, 'CurrentFigure', h1)
    set(h1,'name',basename)
    hold on
end

%% ACTUAL CODE
if irerebin
    Nuse=Nt;
else
    Nuse=Nrebins;
end
di=zeros(Nt,1); hl=zeros(Nt,1);
for ii=1:Nuse
    if irerebin
        iuse= ii;
    else
        iuse= irebins(ii);
    end
    if ~isempty(dF)
        % go through each unbinned data point 
       iallowbins=getallowedbins( rawdata(iuse,:),bindata,myedges,dF);
       okbins=bindata(logical(iallowbins),1:Np);
    else
        okbins=bindata(:,1:Np);
    end
   [dmin,ibin,dPmax]= getmydist(rawdata(iuse,:),okbins,Wts(iuse),Lk);
   if dPmax<dtol
        %ibini=find(ismember(bindata,okbins(ibin,:),'rows'));
        ibini = FindRow(bindata(:,1:Np), okbins(ibin,:));            
        oldbinnum = newbins(iuse);  
        % store the distance
        di(iuse) =dmin;
        dmaxold=dmax;
        dmax=max([dmin dmaxold]); % keep track of the longest distance an observation has to get clustered

        try
            newbins(iuse)=sbins(ibini); %if it passes the tolerance test, replace the bin number
        catch
            okbins(ibin,:)
            ibini
            bindata(ibini,:)
            pause
        end
        % leave the distance (dPs) equal to 0
        %% 2Plot?
        if iplot
            if isnan(oldbinnum) || (oldbinnum ~= newbins(iuse)) %|| dmin==dmax
            curr_obs=rawdata(iuse,:);
            min_cen=bindata(ibini,1:Np); % just the means
            plotdata=[curr_obs;min_cen];
            
            hl(iuse) = plot(plotdata(:,ix),plotdata(:,iy),'k-'); %draw the line from the obsv data to the rebinned centroid
            %text(mean(plotdata(:,1)),mean(plotdata(:,2)),sprintf('%1.2f',dPmax))
            end
        end
   else
        newbins(iuse)=nan; %else set it to nan
        dPs(iuse)=dPmax;
   end
end %% END ACTUAL CODE
varargout{1} = dmax;

%% START FINAL PLOTTING
if iplot

%     set(hl(logical(di == dmax)),'color','r')
%     xdm=get(hl(logical(di == dmax)),'xd'); ydm=get(hl(logical(di == dmax)),'yd');
%     dbarmax='$\bar{d}_{max}$';
%     ht=text(mean(xdm),mean(ydm),sprintf('%s = %1.3f',dbarmax,dmax));
%     set(ht,'interpreter','latex');
    nprow=4;
    if isempty(dF)
    rectangle('position',[.5 0.01 .62 .05*nprow+.01],'facecolor','w','edgecolor','k')
    icol=-1;
    newloc = nan(Nbins,Np);
    dmaxold = 0;
    for ii=1:Nbins
        idata=find(newbins==sbins(ii));
        ndata=length(idata);
        newloc(ii,:) = mean(rawdata(idata,:),1);
        % find the max distance between newloc and its assoc observation
        [foo,ibar,baz,dmax,imax]= getmydist(newloc(ii,:),rawdata(idata,:),Wts(iuse),Lk);
        if dmax> dmaxold
            ilongd = idata(imax);
            ilongbin = sbins(ii);
            dmaxold=dmax;
        end
        plot(rawdata(idata,ix),rawdata(idata,iy),'.','markerfacecolor',colorstr{ii},'markeredgecolor',colorstr{ii})
        if ndata>1
            scatter(newloc(ii,ix),newloc(ii,iy),80,'o','markerfacecolor',.8*colorstr{ii},'markeredgecolor','k');%,'markersize',24);
            line([newloc(ii,ix) bindata(ii,ix)],[newloc(ii,iy) bindata(ii,iy)],'color',.8*colorstr{ii},'linestyle','--','linewidth',2)
        end
        
        hLine=scatter(bindata(ii,ix),bindata(ii,iy),80,'o','markerfacecolor','k','markeredgecolor','k');%,'markersize',24);
        text(bindata(ii,ix)-ndig*.004,bindata(ii,iy),sprintf('%d',ndata),'fontsize',6,'color','w')
        icol = icol + ~mod(ii-1,nprow);
        irow = mod(ii-1,nprow);
        scatter(.54+icol*.2, .05*nprow-.01-.05*irow,200,'o','markerfacecolor',.8*colorstr{ii},'markeredgecolor','k')
        text(.57+icol*.2, .05*nprow-.01-.05*irow,sprintf('Cluster %d',ii))
        % new loc
    end
    % redo distance calculation based on newloc
%     di=zeros(Nt,1); ibin=zeros(Nt,1);
%     for ii=1:Nbins
%         if irerebin
%             iuse= ii;
%         else
%             iuse= irebins(ii);
%         end
%         [dmin,ibin(iuse),dPmax]= getmydist(rawdata(iuse,:),newloc,Wts(iuse),Lk);
%         di(iuse) =dmin;
%         dmaxold=dmax;
%         dmax=max([dmin dmaxold]); % keep track of the longest distance an observation has to get clustered
%     end
        %ilongd = find(di == dmax); ilongbin = ibin(ilongd);
        plotdata=[rawdata(ilongd,:);newloc(ilongbin,:)];
            
        hr = plot(plotdata(:,ix),plotdata(:,iy),'r-'); 
        dbarmax='$\bar{d}_{max}$';
        ht=text(mean(plotdata(:,1)),mean(plotdata(:,2)),sprintf('%s = %1.3f',dbarmax,dmaxold));
        set(ht,'interpreter','latex');
    end
    xlim([0 1])
    ylim([0 1])
    axis equal
    SaveAllFig(plt.tofolder)
    %basename=sprintf('UnweightedBins2D_N%d',niter);
    %exportfigure(basename,plt.tofolder,plt.tmpfolder);
end    
hold off
end
function passbins=getallowedbins(idata,bindata,myedges,dF)
[Nb,Np2]=size(bindata);
Np=Np2/2;
   iedge=nan(1,length(dF));
   okbins=nan(Nb,length(dF));
    %find all bins that are within the dFedges
    for jj=1:length(dF)
        dd=dF(jj);
        %find the edges that the idata point is inbetween
        je=myedges(~isnan(myedges(:,dd)),dd );
        Njb=length(je)-1;
        
        je(end)=je(end)*1.001;
        jeL=je(1:end-1);
        jeH=je(2:end);

        iedge(jj)=find(idata(dd)>=jeL & idata(dd)<jeH);
        % find all the bins that are within that edge as well
        %nokbins=0;
        ii=0;
        %while nokbins<1;
            iL=max([ iedge(jj)-ii 1]);
            iH=min([iedge(jj)+1+ii  Njb+1]);        
            okbins(:,jj)=bindata(:,dd)>=je(iL) & bindata(:,dd)<je(iH);
            %nokbins=sum(okbins(:,jj));
         %   disp(['Cannot find allowable re-bin, expanding allowed bins for raw data point: ' num2str(idata)]) 
          %  ii=ii+1;
        %end
    end
    passbins=sum(okbins,2)>length(dF)-1; %

%    passbins=ones(Nb,1);

end

function Index = FindRow(S, s)
nCol  = size(S, 2);
match = S(:, 1) == s(1);
for col = 2:nCol
   match(match) = S(match, col) == s(col);
end
Index = find(match);
end