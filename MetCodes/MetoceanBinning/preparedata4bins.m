function [data,maxs,offsets]=preparedata4bins(datamat,typestr)
%normalize by max so no parameter is unfairly weighted in distance function
%data=[met.all.Hs(1:nS:end), met.all.Tp(1:nS:end), met.all.wdir(1:nS:end) met.all.vdir(1:nS:end) met.all.vspd(1:nS:end)];
switch typestr
    case 'c'
        data=transformdata4christian(datamat);
        disp('Transforming all directions to TOWARDS')
    case 'a'
        data=transformdata4alan(datamat);
        disp('Transforming all directions to WAMIT convention')
    case 'v'
        data=transformdata4vestas(datamat);
        disp('Transforming all directions to Vestas convention')
    case ''
        data=datamat;
        disp('Not transforming anything')
end
%make all data go from 0 to 1
[ndata,foo]=size(data);
offsets=min(data);
newdata=data-repmat(offsets,[ndata 1]);
maxs=max(newdata);
data=newdata./repmat(maxs,[ndata 1]);
end
function data=transformdata4vestas(data)
mydir=data(:,3); %towards, clockwise from north
theirdir=mod(mydir-180,360); %from, clockwise from north
data=[data(:,1:2) theirdir  data(:,4:end)];
[~,Np]=size(data);
if Np>5
    mychopdir=data(:,8); %towards, clockwise from north
    theirchopdir=mod(mychopdir-180,360); %from, clockwise from north
    data=[data(:,1:2) theirdir  data(:,4:end-1) theirchopdir];
end
end
function data=transformdata4alan(data)
% the data I have looks like
%         180
%          ^
%          |
% 90<------  -----> 270           
%          |
%          0
% WAMITs convention
%          0
%          ^
%          |
% 90<------  -----> 270
%          |
%         180
% Theta_Wave = Theta_Wind - Theta_Wave
mywind=data(:,4);
theirwind=mod(-mywind-180,360);

mywave=data(:,3);
twave=mod(-mywave-180,360); %% IS THIS RIGHT?!?!
theirwave=theirwind-twave;
theirwave(theirwave<-180)=theirwave(theirwave<-180)+360;
data=[data(:,1:2) theirwave theirwind data(:,5)];
end

function data=transformdata4christian(data)
% the wind data I have looks like (it is direction FROM)
%         180
%          ^
%          |
% 90<------  -----> 270           
%          |
%          0
% Direction TOWARDS
%          0
%          ^
%          |
% 270<------  -----> 90
%          |
%         180
% Theta_Wave = Theta_Wind - Theta_Wave
mywind=data(:,4);
theirwind=mod(mywind-180,360);

mywavedir=data(:,3);
%twave=mod(mywave-180,360);
%theirwave=theirwind-mywavedir;
%theirwave(theirwave<-180)=theirwave(theirwave<-180)+360;
data=[data(:,1:2) mywavedir theirwind data(:,5)];
end