function [Q1,Q2,Q3,posvec,pcell,collapse_times] = spp_gen(Type,P,plot_toggle)
%{ 
Input: 
plot_toggle = if this is 0, there will be no plots generated, if this is
any other number the code will generate plots for MSD, VACF and P(r(t))
compared to data. 

Type = 1 for a Levy walk, 2 for a heterogeneous SPP model
The Levy walk model needs four inputs, two to describe the shape of the
power law distribution (exponent and mean run time), one constant speed
parameter and finally a parameter to control the amount of positional error
added to trajectories. 
P = vector of inputs for either model. 
P for Levy walk = [power_law_exponent mean_run_time constant_speed
positional_error]

The heterogeneous SPP model has seven inputs, five related to the bivariate
distribution of speed and rotational noise (two means, two standard
deviations and a correlation parameter), the mean run time for an
exponential distribution, and a parameter to control the amount of
positional error on trajectories. 
P for heterogeneous SPP = [mu_speed sigma_speed mu_noise positional_error
sigma_noise correlation mean_run_time]
A note on positional error: while this is mentioned in the corresponding
manuscript for this code, it merits an additional comment. In the mouse
fibroblast data we worked on there was a sharp decrease in the velocity
auto-correlation function after the first timestep. We determined that this
was due to positional errors from imaging noise, and added this to the code
to replicate that effect to match the VACF. You can set this number to zero
if you do not have this issue. 

Output; 
Q1, Q2, Q3 - linear regression goodness to fit for the mean-squared
displacement, velocity auto-correlation function and displacement
probability distribution to experimental data, respectively. 
posvec - matrix of position and orientation data for the chosen model 

Additional notes: 
It is important to note that there are several parameters in this code that
the user is free, and encouraged to, change. The number of particles (np)
can impact your statistics as well as the integration timestep size (dtt). 

Acknowledgements: 
The MSD calculation portion of this code was inspired by Kehl.m, which can
be found on the MATLAB file exchange. 
Kehl.m measures the mean squared displacement (MSD) from a trajectory.
Kehl is written by Maxime Deforet, May 21 2013. MSKCC.
%}

%% global parameters
if nargin < 3
    plot_toggle = 0;
end
dtt = 1;
% IMPORTANT: this must be run in the same directory as traj_info.mat
% generated by data_gather.m or you need to add a direct path to the load
% function
load('traj_info.mat')
msd = mean(msd);
acorr = mean(acorr);
f = find(acorr < 0);
if isempty(f) == 0
    f = f(1)-1;
elseif isempty(f) == 1
    f = length(acorr);
end
acorr = log10(acorr);
T = length(msd)+1;
ts = round(T/dtt);
np = 100;
t = 1:length(msd);
t1=logspace(0,log10(t(end)),length(t));
msd = interp1(log10(t),msd,log10(t1));
vacf = interp1(log10(t),acorr,log10(t1));



%% Levy walk
% parameters v0, mu, t0, delta
if Type == 1
    rng('shuffle')
    mu = P(1);
    t0 = P(2);
    v0 = P(3);
    posvec = cell(np,1);
    Theta = cell(np,1);
    for i = 1:np
        theta = rand(1)*2*pi;
        xpos = 0;
        ypos = 0;
        posvec{i} = zeros(ts,3);
        Theta{i} = zeros(ts,1);
        cran = rand(1);
        c = round(cran);
        n = 0;
        for t=1:ts
            if n == 0
                if c == 0
                    check = 0;
                    while check == 0
                        u = rand;
                        n = t0*(u.^(-1/mu)-1);
                        if n > 1
                            n = round(n);
                            check = 1;
                        end
                    end
                    a = 0;
                    c = 1;
                elseif c == 1
                    n = 1;
                    a = 2*pi;
                    c = 0;
                end
            end
            dtheta = a*randn(1);
            theta = mod(theta + dtheta*sqrt(dtt),2*pi);
            xpos = xpos + v0*cos(theta)*dtt;
            ypos = ypos + v0*sin(theta)*dtt;
            posvec{i}(t,1) = xpos;
            posvec{i}(t,2) = ypos;
            posvec{i}(t,3) = t;
            Theta{i}(t,1) = theta;
            n = n - 1;
        end
    end
end

%% Heterogeneous SPP
if Type == 2
    rng('shuffle')
    muv = P(1);
    sv = P(2);
    mun = P(3);
    sn = P(5);
    pp = P(6);
    tau = P(7);
    posvec = cell(np,1);
    Theta = cell(np,1);
    for i = 1:np
        theta = rand(1)*2*pi;
        xpos = 0;
        ypos = 0;
        posvec{i} = zeros(ts,3);
        Theta{i} = zeros(ts,1);
        n = 0;
        v0 = 0;
        cran = rand(1);
        c = round(cran);
        for t=1:ts
            if n == 0
                if c == 0
                    check = 0;
                    while check == 0
                        para = mvnRND(muv,mun,sv,sn,pp);
                        arun = para(2);
                        v0 = para(1);
                        a = arun;
                        check = 1;
                    end
                    check = 0;
                    while check == 0
                        if tau == 0
                            n = 1;
                            check = 1;
                        end
                        if tau > 0
                            n = exprnd(tau);
                            if n > 1
                                n = round(n); 
                                check = 1;
                            end
                        end
                    end
                    c = 1;
                elseif c == 1
                    n = 1;
                    a = 2*pi;
                    c = 0;
                end
            end
            dtheta = a*randn(1);
            theta = mod(theta + dtheta*sqrt(dtt),2*pi);
            xpos = xpos + v0*cos(theta)*dtt;
            ypos = ypos + v0*sin(theta)*dtt;
            posvec{i}(t,1) = xpos;
            posvec{i}(t,2) = ypos;
            posvec{i}(t,3) = t;
            Theta{i}(t,1) = theta;
            n = n - 1;
        end
    end
end
disp('Done generating trajectories...')
%% add positional error
Dt = P(4);
for i = 1:length(posvec)
    dpos = Dt*randn(length(posvec{1}),2);
    posvec{i}(:,[1 2]) = posvec{i}(:,[1 2]) + dpos;
end

%% Generate pcell for P(r(t)) collapse in prt_collapse.m
% Feel free to change the times in collapse_times to something more
% visually appealing, this was just a generic way of including times to
% calculate P(r(t))
collapse_times = [round(ts/5)/dtt round(2*ts/5)/dtt round(3*ts/5)/dtt round(4*ts/5)/dtt];
ts = length(posvec{1});
np = length(posvec);
pcell = cell(length(collapse_times),1);
for ii = 1:length(collapse_times)
    dt = collapse_times(ii);
    p = zeros(ts-dt-2,np);
    nsampled = 0;   
    for i = 1:length(posvec)
        for j = 1:ts-dt-1 
            x = posvec{i}(j,1);
            xdt = posvec{i}(j+dt,1);
            y = posvec{i}(j,2);
            ydt = posvec{i}(j+dt,2);
            dr = sqrt((xdt - x)^2 + (ydt - y)^2);
            p(j,i) = dr;
            nsampled = nsampled + 1;
        end
    end
    pcell{ii} = p(:);
end
%% MSD and VACF calculation
Trajectory = posvec;
T =  size(Trajectory{1},1); % T is the number of point in the trajectory;
N = size(Trajectory,1); % N is the number of cells
[ I j ] = find(triu(ones(T), 1)); % list of indices of possible pairings
pos = cell2mat(Trajectory);
xpos = reshape(pos(:,1),T,N).'; % rows = particles, columns = time
ypos = reshape(pos(:,2),T,N).';
D = (xpos(:,I)-xpos(:,j)).^2 + (ypos(:,I)-ypos(:,j)).^2;


vx = diff(xpos,1,2);
vy = diff(ypos,1,2);
[ i, J ] = find(triu(ones(T-1), 1)); % list of indices of possible pairings
V = (vx(:,i).*vx(:,J) + vy(:,i).*vy(:,J))./(sqrt(vx(:,i).^2+vy(:,i).^2).*sqrt(vx(:,J).^2+vy(:,J).^2));
dt = -( Trajectory{1}(I,end) - Trajectory{1}(j,end) );
dvt = -( Trajectory{1}(i,end) - Trajectory{1}(J,end) );

% Then we sort dt in ascending order, and sort D in the same way.
[DT,idx]=sort(dt(:));
DD = D(:,idx);
[DVT,idx]=sort(dvt(:));
VV = V(:,idx);


% Now we have to compute the mean DD for each possible value of DT.
% Let's get the first and last indices of each set of DT
First_idx=find(DT-circshift(DT,1)~=0);
Last_idx=find(DT-circshift(DT,-1)~=0);
tau=DT(First_idx);
tau1=logspace(0,log10(tau(end)),length(tau));
C = cumsum([zeros(size(DD,1),1) DD],2);
div = repmat((Last_idx-First_idx+1).',size(C,1),1);
MSD = mean(log10((C(:,Last_idx+1)-C(:,First_idx))./div));
MSDmat = ((C(:,Last_idx+1)-C(:,First_idx))./div);
First_idx=find(DVT-circshift(DVT,1)~=0);
Last_idx=find(DVT-circshift(DVT,-1)~=0);
CV = cumsum([zeros(size(VV,1),1) VV],2); 
div = repmat((Last_idx-First_idx+1).',size(CV,1),1);
VACF = [0 mean(log10((CV(:,Last_idx+1)-CV(:,First_idx))./div))];
VACFmat = ((CV(:,Last_idx+1)-CV(:,First_idx))./div);
MSD = interp1(log10(tau),MSD,log10(tau1));
VACF = interp1(log10(tau),VACF,log10(tau1));
tau = log10(tau1);
disp('Done calculating statistics...')
%%
% Deprecated plotting without error, can enable but not recommended. 
% if plot_toggle == 1
% figure
% scatter(tau,MSD,'b')
% hold on
% scatter(tau,msd,'r')
% xlabel('timescale','FontSize',30)
% ylabel('MSD','FontSize',30)
% set(gca,'FontSize',24)
% figure
% scatter(tau,VACF,'b')
% hold on
% scatter(tau,vacf,'r')
% xlabel('timescale','FontSize',30)
% ylabel('VACF','FontSize',30)
% set(gca,'FontSize',24)
% end
% Q1 = ((MSD-msd).^2)./msd;
% Q2 = ((VACF(1:f)-vacf(1:f)).^2)./vacf(1:f);

ntrials = 20;
MSDbs = zeros(ntrials,size(MSDmat,2));
VACFbs = zeros(ntrials,size(VACFmat,2));
for i = 1:ntrials
    idx = ceil(rand(size(MSDmat,1),1)*size(MSDmat,1));
    MSDbs(i,:) = mean(MSDmat(idx,:));
    VACFbs(i,:) = mean(VACFmat(idx,:));
end

f = find(std(log10(VACFbs)) > 0.5);
if isempty(f) == 0
    f = f(1);
elseif isempty(f) == 1
    f = length(VACF);
end
Q1 = 1 - sum((MSD-msd).^2)/sum((msd-mean(msd)).^2);
Q2 = 1 - sum((VACF(1:f)-vacf(1:f)).^2)/sum((vacf(1:f)-mean(vacf(1:f))).^2);

t = tPRT;
ts = length(posvec{1});
np = length(posvec);
pcell = cell(length(t),1);
for ii = 1:length(t)
    dt = t(ii);
    p = zeros(ts-dt-2,np);
    nsampled = 0;
    for i = 1:length(posvec)
        for j = 1:ts-dt-1 
            x = posvec{i}(j,1);
            xdt = posvec{i}(j+dt,1);
            y = posvec{i}(j,2);
            ydt = posvec{i}(j+dt,2);
            dr = sqrt((xdt - x)^2 + (ydt - y)^2);
            p(j,i) = dr;
            nsampled = nsampled + 1;
        end
    end
    pcell{ii} = p(:);
end

pp = pcell{1};
pp = pp(pp~=0);
[CDF,X] = ksdensity(pp,'function','cdf','support','positive');
PDF = diff(CDF)./mean(diff(X));
X = X(1:end-1);

vq = interp1(X,PDF,xH);
f = 1-isnan(vq);
f = find(f);
Q3 = 1 - sum((vq(f)-pH(f)).^2)/sum((pH(f)-mean(pH(f))).^2);
if plot_toggle ~= 0
   figure(1)
   loglog(xH,pH)
   hold on
   loglog(X,PDF,'g')
   xlabel('timescale','FontSize',30)
   ylabel('P(r(t))','FontSize',30)
   set(gca,'FontSize',24)
end

if plot_toggle ~= 0
    load('traj_info.mat','msd','acorr')
    figure(2)
    clr = parula(3);
    hplot = 1;
    if hplot == 1
    hold on
    lineProps.col = {'r'};
    mseb(log10(1:size(MSDbs,2)),mean((msd)),std((msd)),lineProps,0);
    end
    hold on
    lineProps.col = {clr(2,:)};
    mseb(log10(1:size(MSDbs,2)),mean(log10(MSDbs)),std(log10(MSDbs)),lineProps,0);
    xlabel('timescale','FontSize',30)
    ylabel('MSD','FontSize',30)
    set(gca,'FontSize',24)
    figure(3)
    if hplot == 1figure(2)
    clr = parula(3);
    hplot = 1;
    if hplot == 1
    hold on
    lineProps.col = {'r'};
    mseb(log10(1:size(MSDbs,2)),mean((msd)),std((msd)),lineProps,0);
    end
    hold on
    lineProps.col = {clr(2,:)};
    mseb(log10(1:size(MSDbs,2)),mean(log10(MSDbs)),std(log10(MSDbs)),lineProps,0);
    xlabel('timescale','FontSize',30)
    ylabel('MSD','FontSize',30)
    set(gca,'FontSize',24)
    figure(3)
    if hplot == 1
    hold on
    lineProps.col = {'r'};
    mseb(log10(1:size(MSDbs,2)),mean(log10(acorr)),std(log10(acorr)),lineProps,0);
    end
    hold on
    lineProps.col = {clr(2,:)};
    mseb(log10(1:size(MSDbs,2)),[0 mean(log10(VACFbs))],[0 std(log10(VACFbs))],lineProps,0);
    xlabel('timescale','FontSize',30)
    ylabel('VACF','FontSize',30)
    set(gca,'FontSize',24)
    hold on
    lineProps.col = {'r'};
    mseb(log10(1:size(MSDbs,2)),mean(log10(acorr)),std(log10(acorr)),lineProps,0);
    end
    hold on
    lineProps.col = {clr(2,:)};
    mseb(log10(1:size(MSDbs,2)),[0 mean(log10(VACFbs))],[0 std(log10(VACFbs))],lineProps,0);
    xlabel('timescale','FontSize',30)
    ylabel('VACF','FontSize',30)
    set(gca,'FontSize',24)
end
end

