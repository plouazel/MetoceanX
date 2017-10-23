function rDat=getRyr(time,dat,dattype,iUse,nout,Ryrs,nS,iplot,varargin)
if nargin<4
    iUse=ones(length(time),1);
end
if nargin<5
    nout = 0;
end
if nargin<6
    Ryrs=[1 10 50];
end
if nargin<7
    nS=1;
end
if nargin<8
   iplot=1; 
end
    

% fname=[mettype '_rawdata'];
% fulldatafile=which(fname);
% [metpath,metname,foo]=fileparts(fulldatafile);
% if isempty(fulldatafile)
%     [metpath,filename,foo]=uigetfile('*.mat',['Please find ' fname '.mat']);
% end
% met=readhindcast(metpath,mettype,0);

maxdat = runWeibullExe(time,dat,iUse,nout,nS,Ryrs);
maxs=maxdat(:,1);
rDat=maxs;

%% POST-PROCESS
weibullfile=which('weibull.exe'); % if you've made it this far, then weibull.exe exists
[wdir,wfile,foo]=fileparts(weibullfile);
wdir = [wdir filesep];
load([wdir 'temp_in2.dat'],'-ascii')
datp=flipud(temp_in2(:,1)); % why flip?
prankp=flipud(temp_in2(:,2)); % why flip?
rawdat=[datp prankp];
load([wdir 'temp_out2.dat'],'-ascii') ;
%create fits
xfit  = temp_out2(:,1);
pfit  = temp_out2(:,2);

pfiti=pfit(end:-1:1);
xfiti=xfit(end:-1:1);
fitdat=[xfiti pfiti];
if iplot
       % plot algorithm
    hF=figure('name',['Extreme ' dattype]);
    cstr='k';
    set(hF,'position',2*[.5 .5 4 2.75],'units','inches');
    set(hF,'position',2*[.5 .5 4 2.75],'units','inches'); % have to repeat bc stupid windows does not get it right
    h1=gca;
    set(h1,'position',5*[.2 .1 1.3 .9],'units', 'inches');
    set(h1,'position',5*[.2 .1 1.3 .9],'units', 'inches'); % have to repeat bc stupid windows does not get it right
    h2=semilogy(rawdat(:,1),rawdat(:,2),[cstr 'o'],fitdat(:,1),fitdat(:,2),[cstr '-'],maxdat(:,1),maxdat(:,2),['r' 'x']);
    grid on
    xlabel([dattype ' (m)'],'FontSize',12,'FontWeight','demi');
    ylabel('P(Exceedence)','FontSize',12,'FontWeight','demi');
    for ii=1:length(Ryrs)
        text(maxdat(ii,1)+.3,maxdat(ii,2), sprintf('%d-yr, %2.1f %s',Ryrs(ii),maxdat(ii,1),'m'))
    end
    
end
end