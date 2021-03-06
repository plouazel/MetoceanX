% Input parameters
write2dfcn=which('write2dScatterTable'); 
[allmetdir,wname,wext] = fileparts(write2dfcn);

ProjectName = 'WFKOWL';%'WFKOWL'; %'WFA'
metname='WAM10_5700N_0181W';%'WAM10_5700N_0181W';%'IHCv2';%;%;

MetType = 'Current';
if strcmp(MetType,'Wave')
    HsType = 'chop'; % use 'combo' (default),  'swell', or 'chop'
else
    HsType = '';
end
wP=0;%-20; %deg WindFloat heading (+CCW relative to due N, like OrcaFlex convention), only used for "phase origin" of wave direction


metdir = [allmetdir filesep ProjectName filesep]; 
% Table dir/filenames
tdir=[metdir  'ScatterTables' filesep];

tablename=[ ProjectName '_' metname '_' MetType HsType '_Orient_' sprintf('%+04d',wP)  '_ScatterTable'];
% output filenames
tname=[tdir tablename '.csv'];
matname = [tdir tablename '.mat'];

%% DISCRETIZATION
wBin=[0:30:360];% deg
Tpbin = 0:1:25;%0:1:15;%;% 
Hsbin =0:1:13;%0:.5:6;% ;%;%% 
Vspdbin = 0:2:28;
Uspdbin = 0:.05:0.6;


%% HEADER for .csv
%print header
if strcmp(ProjectName,'WFM')
    strh{1}='N.B: Wdir (wave heading TOWARDS)= N=0deg E=270 deg S= 180 deg W=90deg';
    strh{2}= 'Data from hindcast set 1/1/1990 00:00 to 12/31/2015 24:00, from Effiage';
    strh{3}='Location: Lat 42.86389n, Lon 3.282501w, Depth 200m, Hub-Height 100m';
    strh{4}='Table lists percentages';
elseif strcmp(ProjectName,'WFA') && strcmp(metname,'M3010613')
    strh{1} = 'N.B: Wdir (wave heading TOWARDS)= N=0deg E=270 deg S= 180 deg W=90deg ';
    strh{2} ='Data from hindcast set 1/1/1954 00:00 to 12/31/2008 21:00, created by Oceanweather Inc. on 11/06/2009 ';
    strh{3} = 'Location: M3 Grid Point 10613  Lat 41.5n  Long 9.0w  Depth 159.5327m ';
    strh{4} ='Table lists percentages';
else
    strh{1}=ProjectName;
    strh{2} = 'N.B: Wdir (wave heading TOWARDS)= N=0deg E=270 deg S= 180 deg W=90deg ';
end