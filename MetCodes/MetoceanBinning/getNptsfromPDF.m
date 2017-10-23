function [mypts]=getNptsfromPDF(myfunx,myfuny,N)
% get N points on the interval [0,1] based on the PDF of how you want the pts distributed
myfunx
if min(myfunx)~=0 || (max(myfunx)<1-1e-3 || max(myfunx)>1+1e-3)
	error('myfunx must be on the closed interval [0,1]')
end

if N<2
	error('must have at least 2 points on the interval')
end
%scale function so that the integral is 1 on the interval [0,1]
myint=trapz(myfunx,myfuny);
mytol=1e-4;
dy=mytol*100;
while (myint-1)-mytol>0 || (myint-1)+mytol<0
	sy=sign(myint-1);
	myfuny=myfuny*(1-sy*dy);
	myint=trapz(myfunx,myfuny);
end
% partition the function such that there are N
targetint=1/(N-1);
%lay down a linear grid to start

newpts=linspace(0,1,N);
intvs=[newpts(1:end-1) newpts(2:end)];
newys=interp1(myfunx,myfuny,newpts);
myints=crappyint(myfunx,myfuny,newpts,10);
offs=myints-targetint;
dt=mytol;
kspr=.1;
while max(offs)-mytol>0 || min(offs)+mytol<0 
	forces=kspr*offs;
	dx=diff(forces);
	newpts(2:end-1)=newpts(2:end-1)+dx; %do not worry about the last integral
	newys=interp1(myfunx,myfuny,newpts);
	myints=crappyint(myfunx,myfuny,newpts,10);
	offs=myints-targetint;
end
mypts=newpts;
end

function myints=crappyint(myfunx,myfuny,newpts,Nref)
N=length(newpts);
myints=zeros(1,N-1);
ifast=1;
if ifast
for ii=1:N-1
	iint=linspace(newpts(ii),newpts(ii+1),2^Nref);
	yint=interp1(myfunx,myfuny,iint);
	myints(ii)= trapz(iint,yint);
end 
else
myints=diff(newpts).*(newys(2:end)+newys(1:end-1))/2; %single trap
end
end