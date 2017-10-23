fname='WAM10_5700N_0181W';
basename='/Users/skanner/GoogleDrive/ToCode/WFKOWL_metocean/';
dattype='Hs';
casestr='ESS';
fdir=[basename 'Data' filesep];
iplot=1;
met=readhindcast(fname,fdir,0);
met.shear = 0.14;
met.VspdH=10; %[m]
met.HubH=100;
met.all.vspdhh=met.all.vspd;%*(met.HubH/met.VspdH)^met.shear; %sometime this exists in the met file?

VspdRated=11;

th0=-30/2;
dirbins=th0+[0:30:360];