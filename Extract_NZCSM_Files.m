clc
clear

tic 
%% read info from NZCSM - all files
ncvars =  {'density_rr', 'model_qcl', 'time0'};
projectdir = pwd;
dinfo = dir(fullfile(projectdir, '*umnsaa*.nc')); %calls only umnsaa files
num_files = length(dinfo);
filenames = fullfile(projectdir, {dinfo.name});
time0 = cell(num_files, 1);
umnsaa = dir(fullfile(projectdir, '*umnsaa_cb003*.nc'));
hybridt70 = ncread(fullfile(projectdir, umnsaa.name), 'hybridt70');
hybridt71_pivot = ncread(fullfile(projectdir, umnsaa.name), 'hybridt71');
sfc = dir(fullfile(projectdir, '*nzcsm_003*.nc'));
orog = ncread(fullfile(projectdir, sfc.name), 'orog_model');
lat = ncread(fullfile(projectdir, sfc.name), 'lat');
lon = ncread(fullfile(projectdir, sfc.name), 'lon');
hybridt71 = zeros(size(hybridt71_pivot));
R_earth = 6371229;
transmission_line_height = 10;
for i = 1:size(hybridt71,1)-1
    hybridt71(i+1) = hybridt71_pivot(i);
end
hybridt71(1) = hybridt71_pivot(end);

lwc_f = zeros(1100,1100,141);
density_f = zeros(1100,1100,141);

p = 1;
file_num = 1;
for K = 1:num_files
    tic
    file = filenames{K};
    density = ncread(file, ncvars{1});
    lwc_pivot = ncread(file, ncvars{2});
    time0{K} = ncread(file, ncvars{3});
    lwc = zeros(size(lon,1),size(lat,1),size(hybridt71,1),size(time0{K},1));
    for i = 1:size(hybridt71,1)-1
        lwc(:,:,i+1,:) = lwc_pivot(:,:,i,:);
    end
    lwc(:,:,1,:) = lwc_pivot(:,:,end,:);
    for m = 1:size(time0{K},1)
        disp(p)
        for i = 1:size(lat,1)
            for j = 1:size(lon,1)
                lwc_eq = griddedInterpolant(hybridt71,squeeze(lwc(i,j,:,m)));
                lwc_f(i,j,p) = lwc_eq(orog(i,j) + transmission_line_height);
                density_eq = griddedInterpolant(hybridt70,squeeze(density(i,j,:,m)));
                density_f(i,j,p) = density_eq(orog(i,j) + transmission_line_height);
            end
        end
       p = p + 1;
    end
    toc
end

lwc = lwc_f.*density_f*1000./((R_earth+orog).^2);

% lw = lwc_f;
% lw(lw==0) = NaN;
% figure
% h=pcolor(lon, lat, (lw(:,:,5))'); %must transpose
% set(h, 'LineStyle', 'none')

ncvars =  {'sfc_zonal_wind', 'sfc_merid_wind', 'sfc_rh', 'sfc_temp', 'sfc_air_press', 'sfc_visibility', 'time1', 'low_cloud_base', 'sfc_fog', 'rain_rate', 'snowfall_rate'};
projectdir = pwd;
dinfo = dir(fullfile(projectdir, '*nwpsfc*.nc')); %calls only sfc files, ignores stdl 
num_files = length(dinfo);
filenames = fullfile(projectdir, {dinfo.name});
time1 = cell(num_files, 1);
U = cell(num_files, 1);
rel_humid = cell(num_files, 1);
T_a = cell(num_files, 1);
P_static = cell(num_files, 1);
V_m = cell(num_files, 1);
cloud_height = cell(num_files, 1);
fog = cell(num_files, 1);
rain = cell(num_files, 1);
snow = cell(num_files, 1);
for K = 1:num_files
    disp(K)
    file = filenames{K};
    U{K} = sqrt(ncread(file, ncvars{1}).^2 + ncread(file, ncvars{2}).^2); %converts to magnitude of wind
    rel_humid{K} = ncread(file, ncvars{3});
    T_a{K} = ncread(file, ncvars{4})-273.15; %in deg C
    P_static{K} = ncread(file, ncvars{5});
    V_m{K} = ncread(file, ncvars{6});
    time1{K} = ncread(file, ncvars{7});
    cloud_height{K} = ncread(file, ncvars{8});
    fog{K} = ncread(file, ncvars{9});
    rain{K} = ncread(file, ncvars{10});
    snow{K}= ncread(file, ncvars{11});   
end


cloud_height = cat(3, cloud_height{:}) - orog + transmission_line_height; %calculates cloud base height relative to the terrain

% concatenate into 3D matrix

% if lat{1}==lat{2}
%     lat = lat{1};
% else
%     error('latitude data is not consistent')
% end
% 
% if lon{1}==lon{2}
%     lon = lon{1};
% else
%     error('longitude data is not consistent')
% end

U = cat(3, U{:});
rel_humid = cat(3, rel_humid{:});
T_a = cat(3, T_a{:});
P_static = cat(3, P_static{:})/1000;
V_m = cat(3, V_m{:});
t = cat(1, time1{:});
fog = cat(3, fog{:});
rain = cat(3, rain{:});
snow = cat(3, snow{:});

% finds and adjusts the time to its offset
t1 = ncread(filenames{1}, 'time1');
tunit1 = ncreadatt(filenames{1}, 'time1' , 'units');
tparts1 = textscan(tunit1, '%s since %s', 1);
tparts1{2} = datetime(tparts1{2}, 'inputformat', 'yyyy-MM-dd');
t2 = ncread(filenames{2}, 'time1');
tunit2 = ncreadatt(filenames{2}, 'time1' , 'units');
tparts2 = textscan(tunit2, '%s since %s', 1);
tparts2{2} = datetime(tparts2{2}, 'inputformat', 'yyyy-MM-dd');

if tparts1{2}==tparts2{2}
    if strcmpi(tparts1{1}, 'hours')==1
        t = datetime(tparts1{2}, 'Format','dd-MMM-yyyy HH:mm:ss') + hours(12) + hours(t);
    elseif strcmpi(tparts1{1}, 'minutes')==1
        t = datetime(tparts1{2}, 'Format','dd-MMM-yyyy HH:mm:ss') + minutes(12) + minutes(t);
    elseif strcmpi(tparts1{1}, 'seconds')==1
        t = datetime(tparts1{2}, 'Format','dd-MMM-yyyy HH:mm:ss') + seconds(12) + seconds(t);
    else
        error('the time is unable to be converted')
    end
else
    error('these files have different reference times')
end

if seconds(t(2)-t(1))==seconds(t(3)-t(2))
    delta_t = seconds(t(2)-t(1)); %calculates timestep in seconds
else
    error('time step is not constant')
end

save('2013190612UTC_input', 'cloud_height', 'delta_t', 'fog', 'lat', 'lon', 'lwc', 'orog', 'P_static', 'rain', 'rel_humid', 'snow', 't', 'T_a', 'U', 'V_m')

toc

