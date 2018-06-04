function prt_collapse(pcell,ntrials,collapse_times)
%{
This function finds the best-fit collapse of the displacement probability
distribution calculated at times collapse_times, using ntrials number of
iterations in a bootstrap error propagation. The displacements needed to
construct the distriubtions are in pcell and generated by either
data_gather.m for experimental data or spp_gen.m for simulations. 
%}
bootcell = cell(ntrials,1);
gammin = 0.40;
gammax = 1;
gaminc = 0.01;
gam = [0 gammin:gaminc:gammax];
t = collapse_times;
for nnn = 1:ntrials
    probcell = cell(length(gam),1);
    for g = 1:length(gam)
        disp([num2str(nnn) '   ' num2str(g)])
        str = num2str(gam(g));
        str = str(str ~= '.');
        prob = cell(length(t),1);
        for ii = 1:length(t)
            dt = t(ii);
            p = pcell{ii};
            p = p(p~=0);
            p = p(ceil(length(p)*rand(length(p),1)));
            p = p./(dt)^gam(g);
            [cdf,x] = ksdensity(p,'function','cdf','support','positive');
            pdf = diff(cdf)./mean(diff(x));
            x = x(1:end-1);
            probcell{g,1}{ii}(:,1) = x;
            probcell{g,2}{ii}(:,2) = pdf;
        end
    end
    bootcell{nnn} = probcell;
end

pdferr = cell(length(bootcell{1}),1);
for i = 1:length(bootcell)
    for  j = 1:length(bootcell{i})
        for k = 1:length(bootcell{i}{j})
            for l = 1:length(bootcell{i}{j}{k})
                pdferr{j,1}{k}(l,i) = bootcell{i}{j,1}{k}(l);
                pdferr{j,2}{k}(l,i) = bootcell{i}{j,2}{k}(l,2);
            end
        end
    end
end
pdfall = pdferr;

pdferr = cell(length(pdfall),2);
for i = 1:length(pdfall)
    pdferr{i,1} = cell(length(pdfall{i}),1);
    pdferr{i,2} = cell(length(pdfall{i}),1);
    for j = 1:length(pdfall{i})
        xq = mean(pdfall{i,1}{j}.');
        m = mean(pdfall{i,2}{j}.');
        s = std(pdfall{i,2}{j}.');
        pdferr{i,2}{j} = zeros(size(pdfall{i,2}{j},1),3);
        pdferr{i,2}{j}(:,2) = m;
        pdferr{i,2}{j}(:,3) = s;
        pdferr{i,1}{j}(:,1) = xq;
    end
end
%%
gam = [0 gammin:gaminc:gammax];
gamvec = zeros(length(pdferr),2);
chi2 = [];
for jj = 1:length(pdferr)
    logbins = pdferr{jj,1};
    prob = pdferr{jj,2};
    vec = zeros(length(prob),2);
    for i = 1:length(prob)
        vec(i,1) = logbins{i}(1,1);
        vec(i,2) = logbins{i}(end,1);
    end
    rang = 1;
    xmin = log10(vec(end,1));
    xmax = log10(vec(1,2));
    lsvec = [];
    n = 0;
    lsq = cell(length(t),2);
    for ind = 1:length(t)
        lsq{ind,1} = zeros(length(t),length(prob{ind}(:,1)));
        lsq{ind,2} = zeros(length(t),length(prob{ind}(:,1)));
        for i = 1:length(prob{ind}(:,1))
            ref = log10(prob{ind}(i,2));
            refx = log10(logbins{ind}(i,1));
            for j = 1:length(t)
                val = interp1(log10(logbins{j}(2:75,1)),log10(prob{j}(2:75,2)),refx);
                err = interp1(logbins{j}(2:75,1),prob{j}(2:75,3),10^refx)/(10^val);
                lsq{ind,1}(j,i) = (ref - val)^2;
                lsq{ind,2}(j,i) = err^2;
                lsvec = [lsvec ((ref - val)^2)];
                chi2 = [chi2 ((ref - val)^2)/err];
                n = n+1;
            end
        end
    end
    nel = numel(lsvec);
    lsvec = lsvec(lsvec~=0);
    lsvec = lsvec(isnan(lsvec) == 0);
    lsvec = lsvec(isinf(lsvec) == 0);
    gamvec(jj,:) = [sum(lsvec)/numel(lsvec) gam(jj)];
end

figure
plot(gamvec(2:end,2),gamvec(2:end,1))
xlabel('\gamma','FontSize',20)
ylabel('\chi^2','FontSize',20)

%%
plot_toggle = 1;
if plot_toggle == 1
    
    min_ind = gamvec(gamvec(:,1)==min(gamvec(:,1)),2);
    nn = 1;
    for gamm = [0 min_ind]
        str = num2str(gamm);
        str = str(str ~= '.');
        clr = parula(length(t));
        if nn == 2
            j = round((gamm-gammin)*100+1)+1;
            figure('name',['scaled_gamma_' ],'numbertitle','off')
            hold on
            title(['Scaling factor of gamma = ' num2str(gamm)])
            hold on
        elseif gamm == 0
            j = 1;
            figure('name','no_scaling','numbertitle','off')
            hold on
            title('No scaling factor')
            hold on
        end
        
        for ii = 1:length(t)
            logbins = log10(pdferr{j,1}{ii}(:,1));
            Prob = log10(pdferr{j,2}{ii}(:,2));
            err = pdferr{j,2}{ii}(:,3)./pdferr{j,2}{ii}(:,2);
            hold on
            if gamm == 0
                hold on
                errorbar(logbins(2:end),Prob(2:end),err(2:end),'Color',clr(ii,:))
            end
            if gamm > 0
                hold on
                title(['gamma = ' num2str(gamm)])
                errorbar(logbins(2:end),Prob(2:end),err(2:end),'Color',clr(ii,:))
                
            end
        end
        xlabel('\rho(t)','FontSize',20)
        ylabel('P(\rho(t))','FontSize',20)
        nn = nn + 1;
    end
end