function met=readhindcast(fn,ldir,irun)
fname=[fn '_rawdata'];
%ldir = '..\hindcast\';
%% METOCEAN STANDARD follows from Oceanweather INC:

%% WIND/CURRENT DIRECTION = From which the wind/water is blowing/flowing, clockwise from true north in degrees (meteorological convention).

%% WAVE DIRECTION = To which waves are traveling, clockwise from north in degrees (oceanographic convention).
if irun
    % Which metocean file
    if strncmpi(fn,'IHC',3)
        %% sea/swell data from other file
        DataName = ['WFA-' fn];

        Data2=importdata([ldir fn '_sea_swell.txt']);
        Data2=Data2.data;
        nd2=length(Data2(:,1));
        yyyy2=Data2(:,1);
        mm2=Data2(:,2);
        dd2=Data2(:,3);
        HH2=Data2(:,4);
        MM2=zeros(nd2,1);
        SS2=zeros(nd2,1);
        chopHs=Data2(:,5);
        chopTp=Data2(:,6);
        chopWdir=mod(Data2(:,7)-180,360); % make consistent with OceanWeather, direction TOWARDS wave propagation, clockwise 
        swellHs=Data2(:,8);
        swellTp=Data2(:,9);
        swellWdir=mod(Data2(:,10)-180,360); % make consistent with OceanWeather, direction TOWARDS wave propagation, clockwise 
        time2=datenum(yyyy2,mm2,dd2,HH2,MM2,SS2);
        %% Original Data
        Data=importdata([ldir fn(1:3) '.txt']);
        Data=Data.data;
        nd=length(Data(:,1));
        yyyy=Data(:,1);
        mm=Data(:,2);
        dd=Data(:,3);
        HH=Data(:,4);
        MM=zeros(nd,1);
        SS=zeros(nd,1);
        
        time1=datenum(yyyy,mm,dd,HH,MM,SS);
        iu=ismember(time1,time2);
        time=time1(iu);
        Hs=Data(iu,5);
        Tp=Data(iu,7);
        Wdir=mod(Data(iu,8)-180,360); % make consistent with OceanWeather
        Vspd=Data(iu,9);
        Vdir=Data(iu,10);
        
    elseif  strncmpi(fn,'Bor',3)
        DataName = 'Borselle';
        %We note that the wind and wave directions are given using the nautical (or meteorological) convention: the coming from direction, clockwise from the North (coming from North=0°N, East=90°N, South=180°N, West=270°N).
        Data=importdata([ldir fn '.csv'],',',10);
        Data=Data.data;
        %fid=fopen([ldir fn '.csv'],'r');
        %Data=textscan(fid,'%d,%d,%d,%d,%2.2f,%d,%2.2f,%2.2f,%d,%1.2f,%1.2f,%1.2f','headerLines', 9);
        yyyy=Data(:,1);
        mm=Data(:,2);
        dd=Data(:,3);
        hh=Data(:,4);
        time=datenum(yyyy,mm,dd,hh,zeros(length(hh),1),zeros(length(hh),1));
        Vspd=Data(:,5);
        Vdir=Data(:,6);
        Hs=Data(:,7);
        Tp=Data(:,8);
        Wdir=mod(Data(:,9)-180,360); % TOWARDS, make consistent with Oceanweather
        % water level in m MSL
        % flow velocity.
    elseif strncmpi(fn,'Eif',3)
        DataName = 'Eiffage';
        nHL=3;
        Data=importdata([ldir fn 'Waves' '.txt'],' ',nHL);
        WaveData=Data.data;
        Hs=WaveData(:,1); Tp=WaveData(:,2);
        Wdir=WaveData(:,3); % FROM, +CW Email from <Julien.LARGUIER@eiffage.com> on 02/22/17 % Wave Dir is FROM, +CW Email from <Julien.LARGUIER@eiffage.com> on 02/22/17
        % convert to oceanographic convention
        Wdir = mod(Wdir-180,360); % TOWARDS, make consistent with Oceanweather
        %WaveDataTimeStr=Data.textdata(nHL+1:end,:);
        mmddyyyy=Data.textdata(nHL+1:end,1);
        hhmmss=Data.textdata(nHL+1:end,2);
        waveDateStr=strcat(mmddyyyy,{' '}, hhmmss);
        waveDate=datenum(waveDateStr,'dd.mm.yyyy HH:MM:SS'); % 3-hr, starting 1990
        % CURRENT
        Data=importdata([ldir fn 'Cur0m' '.txt'],' ',nHL);
        CurData=Data.data;
        Uspd=CurData(:,3);
        Udir=CurData(:,4); % FROM, +CW Email from <Julien.LARGUIER@eiffage.com> on 02/22/17
        mmddyyyy=Data.textdata(nHL+1:end,1);
        hhmmss=Data.textdata(nHL+1:end,2);
        curDateStr=strcat(mmddyyyy,{' '}, hhmmss);
        curDate=datenum(curDateStr,'dd.mm.yyyy HH:MM:SS');% 3 hr, starting 1993
        % WIND
        nHL=4;
        Data=importdata([ldir fn 'Wind_100m.txt'],' ',nHL);
        WindData=Data.data;
        yyyymmdd=WindData(:,1);
        HHMM=WindData(:,2);
        Vspd=WindData(:,3);
        Vdir=WindData(:,4); % FROM, +CW Email from <Julien.LARGUIER@eiffage.com> on 02/22/17
        yyyy = floor(yyyymmdd/10000);
        mm = floor((yyyymmdd-yyyy*10000)/100);
        dd = rem(yyyymmdd,100);
        HH=floor(HHMM/100);
        MM=rem(HHMM,100);
        windDate=datenum(yyyy,mm,dd,HH,MM,MM); % 1 hr, starting 1990
        %% DOWNSAMPLE multiple sources
        time=waveDate;
        [iTime]=ismember(time,waveDate); % missing a few years from the wave data..
        time=time(iTime);
        %Uspd=Uspd(iTime);
        %Udir=Udir(iTime);
        iWind=ismember(windDate,time);
        Vdir=Vdir(iWind);%=interp1(windDate,Vdir,time);  
        Vspd=Vspd(iWind);%interp1(windDate,Vspd,time);
        [iWave]=ismember(waveDate,time);
        Hs=Hs(iWave);%interp1(waveDate,Hs,time);      
        Tp=Tp(iWave);%interp1(waveDate,Tp,time); 
        Wdir=Wdir(iWave);%interp1(waveDate,Wdir,time); 
    elseif strncmpi(fn,'WAM',3)
        DataName = 'KOWL';
        nHL=1;
        cname = [ldir fn  '.csv'];
        delim = ',';
        fid = fopen(cname,'r');
        fgetl(fid);
        nline=0;
        while ~feof(fid);
            nline=nline+1;
            csvline = fgetl(fid);
            idelim = strfind(csvline,delim);
            datef = csvline(1:idelim(1)-1);
            islsh = strfind(datef,'/');
            if isempty(islsh)
                time(nline,1) = str2double(datef);
            else
                isp = strfind(datef,' ');
                icol = strfind(datef,':');
                dd = str2double(datef(1:islsh(1)-1));
                mm = str2double(datef(islsh(1)+1:islsh(2)-1));
                yyyy = str2double(datef(islsh(2)+1:isp(1)-1));
                HH = str2double(datef(isp(1)+1:icol(1)-1));
                MM = str2double(datef(icol(1)+1:end));
                time(nline,1)=datenum(yyyy,mm,dd,HH,MM,0);
            end
        end
        Data=importdata([ldir fn  '.csv'],',',nHL);
        %M = csvread([ldir fn  '.csv'],1,0,'A2..A171174') ;
        Data = Data.data;
        ncol=1;
        Vspd = Data(:,2-ncol);% At 100m!!!
        Vdir = Data(:,4-ncol);% FROM which the wind is "approaches"
        Hs = Data(:,5-ncol);
        Tp = Data(:,6-ncol);
        Wdir = Data(:,8-ncol); % From??? NOT SPECIFIED HOORAY!!
        Wdir = mod(Wdir+180,360); % TOWARD +CW
        chopHs= Data(:,10-ncol);
        chopTp = Data(:,11-ncol);
        chopWdir = Data(:,12-ncol);
        chopWdir = mod(chopWdir+180,360);
        swellHs = Data(:,13-ncol);
        swellTp = Data(:,14-ncol);
        swellWdir = Data(:,15-ncol);
        swellWdir = mod(swellWdir+180,360);  % TOWARD +CW
        Uspd = Data(:,16-ncol);
        Udir = Data(:,17-ncol); % toward which the current is flowing
        Udir = mod(Udir+180,360); % FROM which the current is flowing
        fclose(fid);
    else  
        DataName ='';
        Data=importdata([ldir fn '.txt']);
        Data=Data.data;
        ikeep=Data(:,11)<361;
        yyyymm=num2str(Data(ikeep,2),'%06u');
        ddHHMM=num2str(Data(ikeep,3),'%06u');
        Vdir=Data(ikeep,4); %[deg]
        Vspd=Data(ikeep,5); %[m/s]
        Hs=4*sqrt(Data(ikeep,6)); %[m]
        Tp=Data(ikeep,7); %[s]
        Wdir=Data(ikeep,8); %[deg] 
        chopHs=4*sqrt(Data(ikeep,9)); %[m]
        chopTp=Data(ikeep,10); %[s]
        chopWdir=Data(ikeep,11); %[deg]
        swellHs=4*sqrt(Data(ikeep,12)); %[m]
        swellTp=Data(ikeep,13); %[s]
        swellWdir=Data(ikeep,14); %[deg] 
        whitespace=zeros(length(Data(ikeep,1)),1)+32;
        whitespace=char(whitespace);
        time=datenum([yyyymm(:,1:4) yyyymm(:,5:6) ddHHMM(:,1:2) whitespace ddHHMM(:,3:4) ddHHMM(:,5:6)],'yyyymmdd HHMM');
        fclose(fid);
    end

    met.all.time=time;

    met.all.vdir=Vdir;
    met.all.vspd=Vspd;

    met.all.Hs=Hs;
    met.all.Tp=Tp;
    met.all.wdir=Wdir;

    if exist('chopHs','var')
        met.chop.Hs=chopHs;
        met.chop.Tp=chopTp;
        met.chop.wdir=chopWdir;
    end
    if exist('swellHs','var')
        met.swell.Hs=swellHs;
        met.swell.Tp=swellTp;
        met.swell.wdir=swellWdir;
    end
    
    if exist('Uspd','var')
        met.all.Uspd=Uspd;
        met.all.Udir=Udir;
    end

    %met=transformdata4christian(met); % adds a met.all.vtow field
    vtowards=met.all.vdir-180;
    vtowards(vtowards<-90)=vtowards(vtowards<-90)+360;
    vtowards(met.all.vdir>180 & vtowards<90)=vtowards(met.all.vdir>180 & vtowards<90)+360;
    met.all.vtow=vtowards;
    met.name = DataName;
%     data.time.month=str2num(datestr(data.time.raw,'mm'));
%     data.time.hour=str2num(datestr(data.time.raw,'HH'));
%     data.time.day=str2num(datestr(data.time.raw,'dd'));
%     data.time.year=str2num(datestr(data.time.raw,'yyyy'));

    save([ldir fname '.mat'],'-struct','met')
else
    met=load([ldir fname '.mat']);
end

end