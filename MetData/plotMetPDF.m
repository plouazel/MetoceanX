function plotMetPDF(varargin)
MetType = 'Wind';
BiType = 'combo'; % only used for waves
if strcmp(MetType,'Wind')
    vars2plot = {'vdir','vspd'}; %{'Wdir','Tp','Hs'}
    units2plot = {'FROM +CW, deg','m/s'};
    BiType = '';
elseif strcmp(MetType,'Current')
    vars2plot = {'udir','uspd'}; %{'Wdir','Tp','Hs'}
    units2plot = {'FROM +CCW, deg','m/s'};
    BiType ='';
elseif strcmp(MetType,'Waves')
    vars2plot = {'Wdir','Tp','Hs'};
    units2plot = {'FROM +CCW, deg','m','s'};    
end
dmatstr= '.mat';
wP = 0; % do not make things more confusing by setting wP~=0
%metsuffix= ['_' BiType '_Orient_' sprintf('%+04d',wP)  '_ScatterTable' dmatstr];
metsuffix=['_' MetType BiType '_Orient_' sprintf('%+04d',wP)  '_ScatterTable' dmatstr];

nMet=0;
if nargin>0
    % you have entered a list of input files
    if nargin==1
        if iscell(varargin{1})
            allmets=varargin{1};
            nuMet=length(allmets);
        else
            nuMet=1;
        end
    else
        nuMet=nargin; 
    end
    for jj=1:nuMet
        if iscell(varargin{1})
            %then we are inputing a cell array of strings
            fullfile=allmets{jj};                    
        else
            fullfile=varargin{jj};
        end 
        dotmat = strfind(fullfile,dmatstr);
        if ~isempty(dotmat)
            fullfile = fullfile(1:dotmat);
        end
        slashes=strfind(fullfile,filesep);
        if isempty(slashes)
            metname = fullfile;
            metfile = which([metname metsuffix]);
            slashes=strfind(metfile,filesep);
            metdir = metfile(1:slashes(end));
        else
            metname = fullfile(slashes(end)+1:end);
            metdir = fullfile(1:slashes(end));
            metfile = [metdir metname metsuffix];
        end
        if exist(metfile,'file')
            nMet=nMet+1;
            load(metfile);%readhindcast(metname,metdir,0);
            if exist('WaveScatter','var')
                Scatter= WaveScatter;
            end
            Scatters(nMet) = Scatter;
            mymet{nMet} = metfile;
        else
            warning([fullfile ' does not exist. Moving on...'])
        end  
    end
else
    % choose met files manually
end
nD = length(vars2plot);
for jj=1:nMet
    

   for kk= 1:nD
       
      nE(kk) = length(Scatters(jj).(vars2plot{kk}));
   end
   % re-order probability stuff so directional is at the end
    if nD>2
        ProbMat = zeros([nE(2:end) nE(1)]);
    else
        ProbMat = zeros(nE(2:end),1, nE(1));
    end
    Vars2Plot={vars2plot{2:end},vars2plot{1}};
    Units2Plot={units2plot{2:end},units2plot{1}};

%        if ~isempty(strfind(vars2plot{kk},'wdir')
%            towdata = Scatters(jj).(vars2plot{kk});
%            fromdata = mod(180-Scatters(jj).Wdir,360); % FROM +CW
%        end
    for kk= 1:nD
        for nn = 1:nE(1)
            ProbMat(:,:,nn) = Scatters(jj).Prob{nn};
        end
        ProbDF(jj).(Vars2Plot{kk}) = Scatters(jj).(Vars2Plot{kk});
        if ~isempty(strfind(Vars2Plot{kk},'dir'))
            ProbDF(jj).(Vars2Plot{kk}) = mod(180-Scatters(jj).(Vars2Plot{kk}),360); % FROM +CW            
        end
        ProbDF(jj).([Vars2Plot{kk} 'Unit']) = Units2Plot{kk};
        for pp=1:length(ProbDF(jj).(Vars2Plot{kk}))
            pj = find(ProbDF(jj).(Vars2Plot{kk})(pp) == Scatters(jj).(Vars2Plot{kk}));
%             if ~isempty(strfind(Vars2Plot{kk},'dir'))
%                 Scatters(jj).(Vars2Plot{kk})
%                 ProbDF(jj).(Vars2Plot{kk})(pp)
%                 pj
%                 pause
%             end
            if kk==1    
                ProbDF(jj).(([Vars2Plot{kk} 'Prob']))(:,pp) = sum(sum(ProbMat(pj,:,:)))/100; % there has to be a more elegant way to do this...
            elseif kk==2 && nD>2
                ProbDF(jj).(([Vars2Plot{kk} 'Prob']))(:,pp) = sum(sum(ProbMat(:,pj,:)))/100;
            elseif kk==3 || (kk==2 && nD==2)
                ProbDF(jj).(([Vars2Plot{kk} 'Prob']))(:,pp) = sum(sum(ProbMat(:,:,pj)))/100;
            end
        end
     end
%    for kk=1:nD
%        % sum over other directions.. 
%        if nD==2
%            ProbDF(jj).(([vars2plot{kk} 'Prob'])) = sum(ProbMat,kk);
%        elseif nD==3
%            ProbDF(jj).(([vars2plot{kk} 'Prob'])) = sum(sum(ProbMat,kk));
%        end
%    end
%    nW = length(Scatters(jj).Wdir);
%    nT = length(Scatters(jj).Tp);
%    nH = length(Scatters(jj).Hs);
%    ProbMat = zeros(nT,nH,nW);
%    for nn = 1:nW
%       ProbMat(:,:,nn) = Scatters(jj).Prob{nn};
%    end
%    %Tp
%       ProbDF(jj).Tp = Scatters(jj).Tp;
%    for kk = 1:length(ProbDF(jj).Tp)
%         ProbDF(jj).TpProb(:,kk) = sum(sum(ProbMat(kk,:,:)))/100;
%    end
%    ProbDF(jj).TpUnit = 's';
%    %Hs
%    ProbDF(jj).Hs = Scatters(jj).Hs;
%    ProbDF(jj).HsUnit = 'm';
%    for kk = 1:length(ProbDF(jj).Hs)
%         ProbDF(jj).HsProb(kk,:) = sum(sum(ProbMat(:,kk,:)))/100;
%    end
%    %Wdir 
%    ProbDF(jj).Wdir = mod(180-Scatters(jj).Wdir,360); % FROM +CW
%    ProbDF(jj).WdirUnit = 'FROM +CW, deg';
%    for kk = 1:length(ProbDF(jj).Wdir)
%        jk=find(ProbDF(jj).Wdir(kk) == Scatters(jj).Wdir);
%         ProbDF(jj).WdirProb(:,kk) = sum(sum(ProbMat(:,:,jk)))/100;
%    end   
   ProbDF(jj).Name = [Scatters(jj).Project '_' Scatters(jj).MetData];

end
plotPDF(ProbDF,vars2plot,BiType)
end

function plotPDF(ProbDF,vars2plot,BiType)
colorstr={[0 0 1], [1 0 1],[1 0 0],[1 1 0],[0 1 0],[0 1 1],[.6 .6 0],[0 .6 .6]};
markerstr = {'+','o','d','s','^','v','<','>'};
for kk=1:length(vars2plot)
    kvar = vars2plot{kk};
    kvarp = [kvar 'Prob'];
    figure('name', ['Probability Distribution of ' BiType ' ' kvar ])
    hold on
    for jj=1:length(ProbDF)
       xd = ProbDF(jj).(kvar);
       yd = ProbDF(jj).(kvarp);
       [xds,isort] = sort(xd);

       kunit = ProbDF(jj).([kvar 'Unit']);
       legstr{jj} = strrep(ProbDF(jj).Name,'_','-');
        plot(xds,yd(isort), 'color', colorstr{jj},'linestyle', '-','marker', markerstr{jj})
        
    end
    
   
    xlabel([kvar ' [' kunit ']']);
    ylabel('Probability')
    grid on
    legend(legstr)
    hold off
end
end