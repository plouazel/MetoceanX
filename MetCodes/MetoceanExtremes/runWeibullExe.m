function maxdat = runWeibullExe(time,dat1,iUse,nout,nS,Ryrs)
ndX=length(dat1);
dt=round(mode(diff(time))*3600*24); %[sec] original time step

time=time(iUse);
dat1=dat1(iUse);
time=time(nS:nS:end);
dat1=dat1(nS:nS:end);
nd= length(dat1); % before truncation
maxsize=99999-nout; 
if length(dat1)>maxsize
    dat1=dat1(1:maxsize); % truncate the data series
    warning('Time series had to be truncated due to memory constraints. Blame Exxon.')
end
dat1=sort(dat1,1,'ascend');
fakedt=dt*ndX/nd; % sample rate, if you have downsampled 
lam=1/fakedt;
prank=flipud(1-exp(-lam*(1-(nd-[1:nd]')/nd)));
% add in the outliers
dat1=[dat1; repmat(dat1(end),[nout 1])];
prank=[prank; repmat(prank(end),[nout 1])];
    myweibullexe=which('weibull.exe'); % only returns the first one it finds. Use '-all' to return multiple
    if isempty(myweibullexe)
        error('Please put weibull.exe on your matlab path')
    else
        [mydir,foo,exe]=fileparts(myweibullexe);
        temp1 = [1, nd+nout, nd, fakedt]';
        dat1file=[mydir filesep 'temp_in1.dat'];
        save(dat1file, 'temp1', '-ascii');   
        
        temp2 = [dat1 prank];
        dat2file=[mydir filesep 'temp_in2.dat'];
        save(dat2file, 'temp2', '-ascii');
        
        disp('Running weibull.exe. Thank you Exxon.')
        cdir=pwd;
        cd(mydir);
        if ispc 
            system(myweibullexe);
        else
            %if you're on a mac try installing xcode (5GB!), homebrew and then wine
            % https://www.davidbaumgold.com/tutorials/wine-mac/
            %setenv('PATH', [getenv('PATH') ':/usr/local/bin']) % matlab's path does not line up with system's path
            runstr=['!wine ' myweibullexe];
            eval(runstr)
        end
        load('temp_out1.dat') ;
        
        cd(cdir);
        alpha = temp_out1(1);
        beta  = temp_out1(2);
        theta = temp_out1(3);
        xmpm  = temp_out1(4);
        xp    = temp_out1(5);
        %create fits
%         xfit  = temp_out2(:,1);
%         pfit  = temp_out2(:,2);
%         pfiti=pfit(end:-1:1);
%         xfiti=xfit(end:-1:1);
%         fitdat=[xfiti pfiti];
        ps=1./(3600*24*365*Ryrs);
        maxs= alpha*(-log(ps)).^(1/beta)+ theta;
        maxdat=[maxs' ps'];
        
    end
end