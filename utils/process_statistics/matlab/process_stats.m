%
clear
clc
format shortG
%#ok<*NOPTS>
%#ok<*ST2NM>
%#ok<*SAGROW>
%

%%Parameters (read from simulation input files)
input_file = importdata('./dns.in');
Reyn = cell2mat(strtrim(erase(input_file(6),"! uref, lref, rey")));
Reyn = str2num(Reyn); u_scale = Reyn(1); delta = Reyn(2); Re = Reyn(3);

visc = u_scale*delta/Re;

input_file = importdata('./heat.in');
Pr = cell2mat(strtrim(erase(input_file(3),"! Pr               (Prandtl number)")));
Pr = str2num(Pr);
alpha = visc/Pr;

%%Cell- and face-centered grid coordinates
y = load('./data/grid.out');
yc = y(2:end-1,3); yf = y(1:end-1,2);
clear y

%%Simulation time
t = load('./data/time.out');
time = t(:,3);
clear t

%%Starting index for reading volume forcing data
I = find(time>=0 & time<=1);
io = I(1);
time = time (io:end);
clear I

%%Mean pressure gradient and bulk velocity
f = load('./data/forcing.out');
dpdx = f(io:end,2); dpdy = f(io:end,3); dpdz = f(io:end,4);
u_bulk = f(io:end,5); v_bulk = f(io:end,6); w_bulk = f(io:end,7);

%%Friction Reynolds number
u_tau = sqrt(delta*mean(-dpdx));
Re_tau = u_tau*delta/visc

%Bulk Reynolds number
Re_bulk = mean(u_bulk)*delta/visc

%%Fluid fraction (IBM)
if isfile('./data/psi.out')
 f = load('./data/psi_u.out');
 u_frac = f(:,2);
 clear f
 f = load('./data/psi_v.out');
 v_frac = f(:,2);
 clear f
 f = load('./data/psi_w.out');
 w_frac = f(:,2);
 clear f
 f = load('./data/psi.out');
 p_frac = f(:,2);
 clear f
end

%%Number of data files per field variable (velocities and scalar)
num_files =  dir('./data/velstats*');
num_files2 =  dir('./data/tmpstats*');

%%Starting index for reading statistics data
clear io
io = 1;

%%Read momentum field data and store in arrays
for i = io:length(num_files)
    temp = load(strcat('./data/',num_files(i).name));
   %Regular double-averaged quantities
    ut(:,i) = temp(:,2);
    vt(:,i) = temp(:,3);
    wt(:,i) = temp(:,4);
    pt(:,i) = temp(:,5);
    u_2t(:,i) = temp(:,6);
    v_2t(:,i) = temp(:,7);
    w_2t(:,i) = temp(:,8);
    p_2t(:,i) = temp(:,9);
    uwt(:,i) = temp(:,10);
   %Fluid-phase-averaged quantities
    uft(:,i) = temp(:,11);
    vft(:,i) = temp(:,12);
    wft(:,i) = temp(:,13);
    pft(:,i) = temp(:,14);
    uf_2t(:,i) = temp(:,15);
    vf_2t(:,i) = temp(:,16);
    wf_2t(:,i) = temp(:,17);
    pf_2t(:,i) = temp(:,18);
    ufwft(:,i) = temp(:,19);
end

%Read temperature field data
if (~isempty(num_files2))
 clear temp
 for i = io:length(num_files2)
     temp = load(strcat('./data/',num_files2(i).name));
    %Regular double-averaged quantities
     tmpt(:,i)    = temp(:,2);
     tmp_2t(:,i)  = temp(:,3);
     utmpt(:,i)   = temp(:,4);
     vtmpt(:,i)   = temp(:,5);
     wtmpt(:,i)   = temp(:,6);
    %Fluid-phase-averaged quantities
     tmpft(:,i)   = temp(:,7);
     tmpf_2t(:,i) = temp(:,8);
     uftmpft(:,i) = temp(:,9);
     vftmpft(:,i) = temp(:,10);
     wftmpft(:,i) = temp(:,11);
    %Solid-phase-averaged quantities
     tmpst(:,i)   = temp(:,12);
     tmps_2t(:,i) = temp(:,13);
 end
end

%%Remove empty array entries
ut(:,1:io)   = [];
vt(:,1:io)   = [];
wt(:,1:io)   = [];
pt(:,1:io)   = [];
u_2t(:,1:io) = [];
v_2t(:,1:io) = [];
w_2t(:,1:io) = [];
p_2t(:,1:io) = [];
uwt(:,1:io)  = [];
%
uft(:,1:io)   = [];
vft(:,1:io)   = [];
wft(:,1:io)   = [];
pft(:,1:io)   = [];
uf_2t(:,1:io) = [];
vf_2t(:,1:io) = [];
wf_2t(:,1:io) = [];
pf_2t(:,1:io) = [];
ufwft(:,1:io) = [];
%
if (~isempty(num_files2))
 tmpt(:,1:io)   = [];
 tmp_2t(:,1:io) = [];
 utmpt(:,1:io)  = [];
 vtmpt(:,1:io)  = [];
 wtmpt(:,1:io)  = [];
%
 tmpft(:,1:io)   = [];
 tmpf_2t(:,1:io) = [];
 uftmpft(:,1:io) = [];
 vftmpft(:,1:io) = [];
 wftmpft(:,1:io) = [];
%
 tmpst(:,1:io)   = [];
 tmps_2t(:,1:io) = [];
end
    
%%Time-averaging
u_mean = mean(ut,2);                %Streamwise velocity
v_mean = mean(vt,2);                %Spanwise velocity
w_mean = mean(wt,2);                %Wall-normal velocity
p_mean = mean(pt,2);                %Pressure
u_rms  = sqrt(mean(u_2t,2));        %r.m.s streamwise velocity fluctuations
v_rms  = sqrt(mean(v_2t,2));        %r.m.s spanwise velocity fluctuations
w_rms  = sqrt(mean(w_2t,2));        %r.m.s wall-normal velocity fluctuations
p_rms  = sqrt(mean(p_2t,2));        %r.m.s pressure fluctuations
uw     = mean(uwt,2);               %Reynolds shear stress
%   
uf_mean = mean(uft,2);              %Streamwise velocity (fluid-phase averaged)
vf_mean = mean(vft,2);              %Spanwise velocity (fluid-phase averaged)
wf_mean = mean(wft,2);              %Wall-normal velocity (fluid-phase averaged)
pf_mean = mean(pft,2);              %Pressure (fluid-phase averaged)
uf_rms  = sqrt(mean(uf_2t,2));      %r.m.s streamwise velocity fluctuations (fluid-phase averaged)
vf_rms  = sqrt(mean(vf_2t,2));      %r.m.s spanwise velocity fluctuations (fluid-phase averaged)
wf_rms  = sqrt(mean(wf_2t,2));      %r.m.s wall-normal velocity fluctuations (fluid-phase averaged)
pf_rms  = sqrt(mean(pf_2t,2));      %r.m.s pressure fluctuations (fluid-phase averaged)
ufwf    = mean(ufwft,2);            %Reynolds shear stress (fluid-phase averaged)
%
if (~isempty(num_files2))
 tmp_mean = mean(tmpt,2);           %Temperature
 tmp_rms  = sqrt(mean(tmp_2t,2));   %r.m.s temperature fluctuations
 utmp     = mean(utmpt,2);          %Streamwise temperature flux
 vtmp     = mean(vtmpt,2);          %Spanwise temperature flux
 wtmp     = mean(wtmpt,2);          %Wall-normal temperature flux
%
 tmpf_mean = mean(tmpft,2);         %Temperature (fluid-phase averaged)
 tmpf_rms  = sqrt(mean(tmpf_2t,2)); %r.m.s temperature fluctuations (fluid-phase averaged)
 uftmpf    = mean(uftmpft,2);       %Streamwise temperature flux (fluid-phase averaged)
 vftmpf    = mean(vftmpft,2);       %Spanwise temperature flux (fluid-phase averaged)
 wftmpf    = mean(wftmpft,2);       %Wall-normal temperature flux (fluid-phase averaged)
%
 tmps_mean = mean(tmpst,2);         %Temperature (solid-phase averaged)
 tmps_rms  = sqrt(mean(tmps_2t,2)); %r.m.s temperature fluctuations (solid-phase averaged)
end

%%Time-averaged 3D fields (for calcualting dispersive quantaties)
if isfile('./data/u_mean.bin')
 precision = 'single';     % precision of the real-valued data
 r0 = [0.,0.,0.];          % domain origin
 non_uniform_grid = true;
 %
 geofile  = "./data/geometry.out";
 data = dlmread(geofile);
 ng   = data(1,:);
 l    = data(2,:);
 dl   = l./ng;
 %
 xp = linspace(r0(1)+dl(1)/2.,r0(1)+l(1),ng(1));
 yp = linspace(r0(2)+dl(2)/2.,r0(2)+l(2),ng(2));
 zp = linspace(r0(3)+dl(3)/2.,r0(3)+l(3),ng(3));
 xu = xp + dl(1)/2.;
 yv = yp + dl(2)/2.;
 zw = zp + dl(3)/2.;
 %
 if(non_uniform_grid)
     f   = fopen('./data/grid.bin');
     grid_z = fread(f,[ng(3),4],precision);
     fclose(f);
     zp = r0(3) + grid_z(:,3)';
     zw = r0(3) + grid_z(:,4)';
 end
 %
 filenamei = ('./data/u_mean.bin');
 datau      = zeros([ng(1),ng(2),ng(3)]);
 f = fopen(filenamei);
 datau(:,:,:) = reshape(fread(f,ng(1)*ng(2)*ng(3),precision),[ng(1),ng(2),ng(3)]);
 fclose(f);
 for i = 1:ng(1)
  for j = 1:ng(2)
    for k = 1:ng(3)
     datau(i,j,k) = datau(i,j,k) - u_mean(k);
   end
  end
 end
 ud_mean(:) = mean(mean(datau,1),2);         %Dispersive streamwise velocity
 ud_rms(:) = sqrt(mean(mean(datau.^2,1),2)); %r.m.s. dispersive streamwise velocity
end
%
if isfile('./data/v_mean.bin')
 filenamei = ('./data/v_mean.bin');
 datav      = zeros([ng(1),ng(2),ng(3)]);
 f = fopen(filenamei);
 datav(:,:,:) = reshape(fread(f,ng(1)*ng(2)*ng(3),precision),[ng(1),ng(2),ng(3)]);
 fclose(f);
 for i = 1:ng(1)
  for j = 1:ng(2)
    for k = 1:ng(3)
     datav(i,j,k) = datav(i,j,k) - v_mean(k);
   end
  end
 end
 vd_mean(:) = mean(mean(datav,1),2);         %Dispersive streamwise velocity
 vd_rms(:) = sqrt(mean(mean(datav.^2,1),2)); %r.m.s. dispersive streamwise velocity
 data = datau.*datav;
 udvd(:) = mean(mean(data,1),2);
end
%
if isfile('./data/w_mean.bin')
 filenamei = ('./data/w_mean.bin');
 dataw      = zeros([ng(1),ng(2),ng(3)]);
 f = fopen(filenamei);
 dataw(:,:,:) = reshape(fread(f,ng(1)*ng(2)*ng(3),precision),[ng(1),ng(2),ng(3)]);
 fclose(f);
 for i = 1:ng(1)
  for j = 1:ng(2)
    for k = 1:ng(3)
     dataw(i,j,k) = dataw(i,j,k) - w_mean(k);
   end
  end
 end
 wd_mean(:) = mean(mean(dataw,1),2);          %Dispersive wall-normal velocity
 wd_rms(:)  = sqrt(mean(mean(dataw.^2,1),2)); %r.m.s. dispersive wall-normal velocity
 data = datau.*dataw;
 udwd(:) = mean(mean(data,1),2);              %Dispersive shear-stress
end
%
if isfile('./data/s_mean.bin')
 filenamei = ('./data/s_mean.bin');
 datas      = zeros([ng(1),ng(2),ng(3)]);
 f = fopen(filenamei);
 datas(:,:,:) = reshape(fread(f,ng(1)*ng(2)*ng(3),precision),[ng(1),ng(2),ng(3)]);
 fclose(f);
 for i = 1:ng(1)
  for j = 1:ng(2)
    for k = 1:ng(3)
     datas(i,j,k) = datas(i,j,k) - tmp_mean(k);
   end
  end
 end
 tmpd_mean(:) = mean(mean(datas,1),2);          %Dispersive temperature
 tmpd_rms(:)  = sqrt(mean(mean(datas.^2,1),2)); %r.m.s. dispersive temperature
 data = dataw.*datas;
 wdtmpd(:) = mean(mean(data,1),2);              %r.m.s. dispersive wall-normal heat flux
 clear datau datav dataw datas
end

% Clear temporary arrays from memory
clear temp ut vt wt pt u_rmst v_rmst w_rmst p_rmst uwt
clear udt vdt wdt pdt u_2t v_2t w_2t p_2t udwdt
clear uf vf wf pf uf_2t vf_2t wf_2t pf_2t ufwft
if (~isempty(num_files2))
 clear tmpt tmp_2t utmpt vtmpt wtmpt
 clear tmpd tmpd_2t udtmpdt vdtmpdt wdtmpdt
 clear tmpft tmpf_2t uftmpft vftmpft wftmpft
 clear tmpst tmps_2t
end