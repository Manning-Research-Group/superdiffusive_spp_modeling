function [datas,pcell,collapse_times,mins] = data_gather(xyzs_id,xyzs_id_columns,mu)
%{
This function generates a matrix of speeds, rotational diffusion
coefficients and run times for trajectories gathered through the ACTIVE
nuclei tracking software. A simple modification can be made to incorporate
other trajectories by inputting them as an xyzs_id matrix with four
columns corresponding to x position, y position, frame index and cell
identification number. 

INPUT: 
xyzs_id: output from the ACTIVE nuclei tracking software
xyzs_id_columns: additional output from the ACTIVE tracking software, if
you are creating your own matrix of positions then set this equal to the
number of columns of your matrix, where the columns should be x position, y
position, frame index and cell index, respectively. 
mu: estimated from the mean-squared displacement (MSD) where alpha is the slope
of the loglog plot of MSD and timescale, where mu = 3 - alpha. 

OUTPUT: 
datas: matrix of cell trajectory data with columns run length, speed,
rotational diffusion and orientation angle
pcell: Matlab cell structure containing displacement probability
information
collapse_times: Frames at which p_cell is calculated
min: frame indicies where the Canny edge detection algorithm identified a
tumble

Note: If you don't have experimental images that you will be generating
trajectory data with, you can create a very simple Brownian particle
simulation in Matlab. Example code: 

dt = 1;
v0 = 1;
a = 1;
ts = 100;
np = 100;
dtmax = ts-1;
posvec = cell(np,1);
for i = 1:np
    rng('shuffle')
    %disp(num2str(i))
    xpos = 0;
    ypos = 0;
    theta = rand(1)*2*pi;
    thetaRan = a * (randn(1,ts));
    posvec{i} = zeros(ts,3);
    for t=1:ts
        dtheta = thetaRan(t);
        theta = theta + dtheta*sqrt(dt);
        xpos = xpos + v0*cos(theta)*dt;
        ypos = ypos + v0*sin(theta)*dt;
        theta = mod(theta,2*pi);
        posvec{i}(t,1) = xpos;
        posvec{i}(t,2) = ypos;
        posvec{i}(t,3) = theta;
    end
end
%}
frameindx = 3;
cellindx = 4;
cell_tracks = [xyzs_id(:,1) xyzs_id(:,2) xyzs_id(:,xyzs_id_columns-1) xyzs_id(:,xyzs_id_columns)];

% look at the data to find the maximum frame and number of cells
nframes = max(cell_tracks(:,frameindx));
ts = nframes;
ncells = max(cell_tracks(:,cellindx));
np = ncells;
ltracks = zeros(ncells,1);
%construct matricies with all the x and y coordinates for each cell
xpos = zeros(ncells,nframes);
ypos = zeros(ncells,nframes);
zpos = zeros(ncells,nframes);
for i=1:ncells
    %boolean matrix identifying all cells with id i
    boolcell = (cell_tracks(:,cellindx) == (i));
    ltracks(i) = nnz(boolcell); % length of that track
    xpos(i,cell_tracks(boolcell,frameindx)) = cell_tracks(boolcell,1);
    ypos(i,cell_tracks(boolcell,frameindx)) = cell_tracks(boolcell,2);    
end

% To account for small gaps in cell tracking, we reconstruct
%trajectories if the cell has not been missing for more than F frames
%(where F=5 is default) 
F = 5;
xcell = cell(ncells,1);
ycell = cell(ncells,1);
for i = 1:ncells
    x = xpos(i,1:end);
    y = ypos(i,1:end);
    f = find(x);
    f = [0 f inf];
    d = diff(f);
    for j = 2:length(d)
        if d(j) > 1 && d(j) < F
            ind = f(j);
            for k = 1:d(j)-1
                if ind == 0
                    x(ind + k) = x(f(j+1)) - (x(f(j+2))-x(f(j+1)));
                    y(ind + k) = y(f(j+1)) - (y(f(j+2))-y(f(j+1)));
                    break;
                else
                    x(ind + k) = x(ind) + k*(x(ind + d(j)) - x(ind))/d(j);
                    y(ind + k) = y(ind) + k*(y(ind + d(j)) - y(ind))/d(j);
                end
            end
        end
    end
    xcell{i} = x;
    ycell{i} = y;
end

% Calculate pcell and collapse_times for prt_collapse.m
dtt=1.0;
collapse_times = [round(ts/5)/dtt round(2*ts/5)/dtt round(3*ts/5)/dtt round(4*ts/5)/dtt];
pcell = cell(length(collapse_times),1);
for ii = 1:length(collapse_times)
    nsampled = 0;
    dt = collapse_times(ii);
    p = zeros(nframes-dt-2,ncells);
    for i = 1:ncells
        for j = find(xcell{i}) 
            if j+dt < ts && xcell{i}(j) ~= 0 && xcell{i}(j+dt) ~= 0
                dr = sqrt((xcell{i}(j+dt) - xcell{i}(j))^2 + (ycell{i}(j+dt) - ycell{i}(j))^2);
                p(j,i) = dr;
                nsampled = nsampled + 1;
            end
        end
    end
    pcell{ii} = p;
end

% We impose a limit cutoff=20 of how few frames a cell may be present in
%consecutively to be included in our analysis in order to increase accuracy
cutoff = 20;
posx = cellfun(@(x) plateau(x.'), xcell,'UniformOutput',0);
posy = cellfun(@(x) plateau(x.'), ycell,'UniformOutput',0);
ind = cellfun(@(x)iscell(x),posx);
posx = posx(find(1-ind),1);
posy = posy(find(1-ind),1);
ind = cellfun(@(x) length(x),posx);
ind = find(ind > 10);
posx = posx(ind,1);
posy = posy(ind,1);
vx = cellfun(@diff,posx,'UniformOutput',0);
vy = cellfun(@diff,posy,'UniformOutput',0);
dr = cell(length(posx),1);
dtheta = cell(length(posx),1);
for i = 1:length(posx)
    if length(posx{i}) > cutoff
        dx = diff(posx{i}(:,1));
        dy = diff(posy{i}(:,1));
        dr{i} = sqrt(dx.^2 + dy.^2);
        dtheta{i} = atan2(dy,dx);
    end
end

% We utilize Otsu's method (from otsu.m) to calculate the Canny thresholds
ts = zeros(length(dtheta),1);
for i = 1:length(dtheta)
    try
        d = dtheta{i};
        d = diff(dtheta{i}+pi);
        d = min(2*pi-abs(d),abs(d)); 
        [thresh] = otsu(d,0);
        ts(i) = thresh/max(d);
    end
end
threshold = mean(ts(ts~=0));
runs = [];
mins = cell(length(dtheta),1);
for i = 1:length(dtheta)
    try
        d = diff(dtheta{i}+pi);
        d = min(2*pi-abs(d),abs(d));
        [ dData, minmax, stats ] = AnalyzeEdgesNP(d, 3, ts(i)/2);
        mins{i} = minmax;
        minmax = [1 minmax.' 1];
        d = diff(find(minmax));
        runs = [runs d(d~=1)];
        %runs = [runs d];
    end
end
mins = cellfun(@find,mins,'UniformOutput',0);
t0 = mean(runs)/(mu-1);
disp('done')
 
dr = dr(~cellfun('isempty',dr));
dtheta = dtheta(~cellfun('isempty',dtheta));
dtheta1 = cellfun(@(x) x+pi,dtheta,'UniformOutput',0);
dtheta1 = cellfun(@(x) diff(x),dtheta1,'UniformOutput',0);
dtheta1 = cellfun(@(x) min(2*pi-abs(x),abs(x)),dtheta1,'UniformOutput',0);
npoints = cellfun(@sum,mins);
npoints = sum(npoints) + 2*length(npoints);
datasT = [];
phis = [];
for i = 1:length(mins)
    try
        ind = mins{i};
        d = diff(ind);
        v = [];
        % calculate speed based on displacment between frames
        for j = 1:length(ind)-1
            v = [v sum(dr{i}(ind(j):ind(j+1)))/(d(j)+1)];
        end
        % estimate rotational diffusion based on variance of orientations
        a = zeros(length(ind)-1,1);
        for j = 2:length(ind)
            a(j-1) = var(dtheta1{i}(ind(j-1):ind(j)));%sqrt(var(dtheta1{i}(ind(j-1):ind(j))))/4;
        end
        
        ind2 = [1 mins{i}.' length(posx{i})].';
        phi = zeros(length(ind2)-1,1);
        test = zeros(length(ind2)-1,5);
        for j = 1:length(ind2)-1
            x = posx{i}(ind2(j):ind2(j+1),1);
            y = posy{i}(ind2(j):ind2(j+1),1);
            if length(x) > 10
                p = polyfit(x,y,1);
                y1 = p(1)*x(1) + p(2);
                y2 = p(1)*x(end) + p(2);
                phi(j) = atan2(y2-y1,x(end)-x(1));
                test(j,:) = [x(1) x(end) y1 y2 atan2(y2-y1,x(end)-x(1))];
            elseif length(x) < 10
                phi(j) = 0;
            end
        end
        phis = [phis.' mod(diff(phi),2*pi).'].';
        T = dtheta{i}(ind);
        datasT = [datasT.' [d v.' a T(2:end) ].'].';
    end
end
ind_good = datasT(:,2) ~= 0;
%datas = [datas.' datasT(ind_good,:).'].';
datas = datasT(ind_good,:);

% Calculate pcell and collapse_times for prt_collapse.m
Ts = max(xyzs_id(:,3));
collapse_times = [round(Ts/5)/dtt round(2*Ts/5)/dtt round(3*Ts/5)/dtt round(4*Ts/5)/dtt];
pcell = cell(length(collapse_times),1);
for ii = 1:length(collapse_times)
    nsampled = 0;
    dt = collapse_times(ii);
    p = zeros(nframes-dt-2,ncells);
    for i = 1:ncells
        for j = find(xcell{i}) 
            if j+dt < Ts && xcell{i}(j) ~= 0 && xcell{i}(j+dt) ~= 0
                dr = sqrt((xcell{i}(j+dt) - xcell{i}(j))^2 + (ycell{i}(j+dt) - ycell{i}(j))^2);
                p(j,i) = dr;
                nsampled = nsampled + 1;
            end
        end
    end
    pcell{ii} = p;
end



% Generate and save the file needed to run spp_gen with the appropriate
% information. 
pp = pcell{1};
pp = pp(pp~=0);
[CDF,X] = ksdensity(pp,'function','cdf','support','positive');
pH = diff(CDF)./mean(diff(X));
xH = X(1:end-1);
tPRT = collapse_times;

posx = cellfun(@(x) plateau(x.'), xcell,'UniformOutput',0);
posy = cellfun(@(x) plateau(x.'), ycell,'UniformOutput',0);
ind = cellfun(@(x)iscell(x),posx);
posx = posx(find(1-ind),1);
posy = posy(find(1-ind),1);
ind = cellfun(@(x) length(x),posx);
ind = find(ind > cutoff);
posx = posx(ind,1);
posy = posy(ind,1);
fr =1; 
vx = cellfun(@(x) [diff(x(:,1))/fr x(1:end-1,2)],posx,'UniformOutput',0);
vy = cellfun(@(x) [diff(x(:,1))/fr x(1:end-1,2)],posy,'UniformOutput',0);
v = cellfun(@(x,y) [sqrt(x(:,1).^2 + y(:,1).^2) x(:,2)],vx,vy,'UniformOutput',0);
v = cell2mat(v);
vel = zeros(max(v(:,2)),3);
for i = 1:max(v(:,2))
    temp = v(v(:,2)==i,1);
    vel(i,:) = [i mean(temp) std(temp)];
end

ts = max(cellfun(@length,vx))+1;
vx = cellfun(@(x) x(:,1),vx,'UniformOutput',0);
vy = cellfun(@(x) x(:,1),vy,'UniformOutput',0);
np = length(vx);
acorr_cell = cell(ts-1,1);
for dt = 0:ts-1
    for i = 1:np
        tm = length(vx{i})+1;
        for t = 1:tm-dt
            try
                v1 = sqrt(vx{i}(t)^2 + vy{i}(t)^2);
                v2 = sqrt(vx{i}(t+dt)^2 + vy{i}(t+dt)^2);
                acorr_cell{dt+1}(i,t) = vx{i}(t)*vx{i}(t+dt)/(v1*v2) + vy{i}(t)*vy{i}(t+dt)/(v1*v2);
            end
        end
    end
end
ntrials = 20;
acorr = zeros(ntrials,ts-1);
np = np-1;
for n = 1:ntrials
    ind = ceil(np*rand(np,1));
    m = cellfun(@(x) x(ind(ind<size(x,1)),:),acorr_cell,'UniformOutput',0);
    m = cellfun(@(x) mean(x(x~=0)),m,'UniformOutput',0);
    m = cell2mat(m).';
    acorr(n,:) = m;
end

np = length(posx);
tmin = 1;
tmax = ts;
msd_cell = cell(tmax-tmin,1);
for dt=1:(tmax-tmin)
    msd_cell{dt} = zeros(i,tmax-dt-tmin);
    dr2 = zeros(ts-dt,np);
    %disp(dt)
    for i = 1:np
        tm = length(posx{i});
        for t=tmin:(tm-dt)
            try
                x = posx{i}(t,1);
                y = posy{i}(t,1);
                
                xdt = posx{i}(t+dt,1);
                ydt = posy{i}(t+dt,1);
                if x~=0 && y~=0 && xdt~=0 && ydt~=0
                    dr2(t,i)=(xdt - x).^2 + (ydt-y).^2;
                    msd_cell{dt}(i,t) = (xdt - x).^2 + (ydt-y).^2;
                end
            end
        end
    end
end

ntrials = 20;
msd = zeros(ntrials,ts-1);
for n = 1:ntrials
    ind = ceil(np*rand(np,1));
    m = cellfun(@(x) x(ind(ind<size(x,1)),:),msd_cell,'UniformOutput',0);
    m = cellfun(@(x) mean(x(x~=0)),m,'UniformOutput',0);
    m = cell2mat(m).';
    msd(n,:) = log10(m);
end
save('traj_info.mat','msd','acorr','xH','pH','tPRT')
