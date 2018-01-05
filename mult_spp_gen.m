function [] = mult_spp_gen(Type,filename)
% I must stress that this is a template for running spp_gen.m multiple
% times to sweep phase space. This is extremely time consuming if run on
% local machines and I highly suggest utilizing your university cluster
% services to split the overall phase space sweep into smaller portions for
% higher efficiency. The goal here is to sweep along each parameter (4 for
% Levy walk, 7 for heterogeneous SPP) and create a phase space. We found
% that all of our velocity auto-correlation functions had a similar initial
% dropoff and thus were able to just fit parameter 4 (positional error) and
% make that a constant, saving a tremendous amount of computation time.
% Input: Type = 1 for Levy walk and 2 for heterogeneous SPP, inputs for
% spp_gen.m
% filename = save name for the matricies of phase space
% plot_toggle = set to 0 by default (advise against changing when sweeping
% phase space)
% IMPORTANT NOTE: Please edit p1, p2, p3, p4, p5, p6, and p7 to look at the
% appropriate areas of phase space. Use the information in datas generated
% by data_gather.m to guide your choice of parameter sweep. The code is
% initially set up to sweep parameters logarithmically but in some cases
% this is too general and a closer examination of phase space is necessary.
% It should be noted that this process can be very time consuming, as a
% disclamer. 
% After finding the configuration of parameters that best approximates your
% data you can then construct the appropriate P vector and run that in
% spp_gen.m with plots enabled to generate plots of MSD, VACF and P(r(t)). 
plot_toggle = 0;

if Type == 1   
p1 = [logspace(-1,0,5)]; 
p2 = [logspace(-2,0,5)];
p3 = [logspace(0,1,5)];
p4 = [0 0.2 0.4 0.6 0.8 1]; 
chi1 = zeros(length(p1),length(p2),length(p3),length(p4));
chi2 = zeros(length(p1),length(p2),length(p3),length(p4));
chi3 = zeros(length(p1),length(p2),length(p3),length(p4));
for i = 1:length(p1)
    for j = 1:length(p2)
        for k = 1:length(p3)
            for l = 1:length(p4)
                        P = [p1(i) p2(j) p3(k) p4(l)];
                        [Q1, Q2, Q3] = spp_gen(Type,P,plot_toggle);
                        chi1(i,j,k,l) = Q1;
                        chi2(i,j,k,l) = Q2;
                        chi3(i,j,k,l) = Q3;
            end
        end
    end
end 

% create plots showing change in parameters around minimum
% feel free to change how this plot is presented to your own specifications
clr = lines(4);
[~, minidx] = max(chi1(:));
[i, j, k, l] = ind2sub( size(chi1), minidx );
figure(3)
hold on
plot(p1,1-squeeze(chi1(:,j,k,l)),'Color',clr(1,:),'LineWidth',4)
hold on
plot(p2,1-squeeze(chi1(i,:,k,l)),'Color',clr(2,:),'LineWidth',4)
hold on
plot(p3,1-squeeze(chi1(i,j,:,l)),'Color',clr(3,:),'LineWidth',4)
hold on
plot(p4,1-squeeze(chi1(i,j,k,:)),'Color',clr(4,:),'LineWidth',4)
legend('\mu','t_0','v_0','\Delta')
set(gca,'FontSize',24)
xlabel('Parameter','FontSize',30)
ylabel('1-R^2','FontSize',30)
set(gca,'YScale','log')

[~, minidx] = max(chi2(:));
[i, j, k, l] = ind2sub( size(chi2), minidx );
figure(2)
hold on
plot(p1,1-squeeze(chi2(:,j,k,l)),'Color',clr(1,:),'LineWidth',4)
hold on
plot(p2,1-squeeze(chi2(i,:,k,l)),'Color',clr(2,:),'LineWidth',4)
hold on
plot(p3,1-squeeze(chi2(i,j,:,l)),'Color',clr(3,:),'LineWidth',4)
hold on
plot(p4,1-squeeze(chi2(i,j,k,:)),'Color',clr(4,:),'LineWidth',4)
legend('\mu','t_0','v_0','\Delta')
set(gca,'FontSize',24)
xlabel('Parameter','FontSize',30)
ylabel('1-R^2','FontSize',30)
set(gca,'YScale','log')

[~, minidx] = max(chi3(:));
[i, j, k, l] = ind2sub( size(chi3), minidx );
figure(3)
hold on
plot(p1,1-squeeze(chi3(:,j,k,l)),'Color',clr(1,:),'LineWidth',4)
hold on
plot(p2,1-squeeze(chi3(i,:,k,l)),'Color',clr(2,:),'LineWidth',4)
hold on
plot(p3,1-squeeze(chi3(i,j,:,l)),'Color',clr(3,:),'LineWidth',4)
hold on
plot(p4,1-squeeze(chi3(i,j,k,:)),'Color',clr(4,:),'LineWidth',4)
legend('\mu','t_0','v_0','\Delta')
set(gca,'FontSize',24)
xlabel('Parameter','FontSize',30)
ylabel('1-R^2','FontSize',30)
set(gca,'YScale','log')

end

if Type == 2
p1 = [logspace(-1,0,5)];
p2 = [logspace(-2,0,5)];
p3 = [logspace(0,1,5)];
p4 = [0 0.2 0.4 0.6 0.8 1]; 
p5 = [logspace(-2,0,5)];
p6 = [-1 -0.5 0 0.5 1];
p7 = [logspace(0,3,5)];
chi1 = zeros(length(p1),length(p2),length(p3),length(p4),length(p5),length(p6),length(p7));
chi2 = zeros(length(p1),length(p2),length(p3),length(p4),length(p5),length(p6),length(p7));
chi3 = zeros(length(p1),length(p2),length(p3),length(p4),length(p5),length(p6),length(p7));
for i = 1:length(p1)
    for j = 1:length(p2)
        for k = 1:length(p3)
            for l = 1:length(p4)
                for m = 1:length(p5)
                    for n = 1:length(p6)
                        for o = 1:length(p7)
                        P = [p1(i) p2(j) p3(k) p4(l) p5(m) p6(n) p7(o)];
                        [Q1, Q2, Q3] = spp_gen(Type,P,plot_toggle);
                        chi1(i,j,k,l,m,n,o) = Q1;
                        chi2(i,j,k,l,m,n,o) = Q2;
                        chi3(i,j,k,l,m,n,o) = Q3;
                        end
                    end
                end
            end
        end
    end
end

% create plots showing change of parameters around minimum
clr = lines(7);
figure(1)
[~, minidx] = max(chi1(:));
[i, j, k, l, m, n, o] = ind2sub( size(chi1), minidx );
figure(1)
hold on
plot(p1,1-squeeze(chi1(:,j,k,l,m,n,o)),'Color',clr(1,:),'LineWidth',4)
hold on
plot(p2,1-squeeze(chi1(i,:,k,l,m,n,o)),'Color',clr(2,:),'LineWidth',4)
hold on
plot(p3,1-squeeze(chi1(i,j,:,l,m,n,o)),'Color',clr(3,:),'LineWidth',4)
hold on
plot(p5,1-squeeze(chi1(i,j,k,l,:,n,o)),'Color',clr(5,:),'LineWidth',4)
hold on
plot(p6,1-squeeze(chi1(i,j,k,l,m,:,o)),'Color',clr(6,:),'LineWidth',4)
hold on
plot(p7,1-squeeze(chi1(i,j,k,l,m,n,:)),'Color',clr(7,:),'LineWidth',4)
hold on
plot(p4,1-squeeze(chi1(i,j,k,:,m,n,o)),'Color',clr(4,:),'LineWidth',4)
legend('\mu_v','\sigma_v','\mu_n','\sigma_n','P','\tau','\Delta')
set(gca,'FontSize',24)
xlabel('Parameter','FontSize',30)
ylabel('1-R^2','FontSize',30)
set(gca,'YScale','log')

figure(2)
[~, minidx] = max(chi2(:));
[i, j, k, l, m, n, o] = ind2sub( size(chi2), minidx );
figure(1)
hold on
plot(p1,1-squeeze(chi2(:,j,k,l,m,n,o)),'Color',clr(1,:),'LineWidth',4)
hold on
plot(p2,1-squeeze(chi2(i,:,k,l,m,n,o)),'Color',clr(2,:),'LineWidth',4)
hold on
plot(p3,1-squeeze(chi2(i,j,:,l,m,n,o)),'Color',clr(3,:),'LineWidth',4)
hold on
plot(p5,1-squeeze(chi2(i,j,k,l,:,n,o)),'Color',clr(5,:),'LineWidth',4)
hold on
plot(p6,1-squeeze(chi2(i,j,k,l,m,:,o)),'Color',clr(6,:),'LineWidth',4)
hold on
plot(p7,1-squeeze(chi2(i,j,k,l,m,n,:)),'Color',clr(7,:),'LineWidth',4)
hold on
plot(p4,1-squeeze(chi2(i,j,k,:,m,n,o)),'Color',clr(4,:),'LineWidth',4)
legend('\mu_v','\sigma_v','\mu_n','\sigma_n','P','\tau','\Delta')
set(gca,'FontSize',24)
xlabel('Parameter','FontSize',30)
ylabel('1-R^2','FontSize',30)
set(gca,'YScale','log')

figure(3)
[~, minidx] = max(chi3(:));
[i, j, k, l, m, n, o] = ind2sub( size(chi3), minidx );
figure(1)
hold on
plot(p1,1-squeeze(chi3(:,j,k,l,m,n,o)),'Color',clr(1,:),'LineWidth',4)
hold on
plot(p2,1-squeeze(chi3(i,:,k,l,m,n,o)),'Color',clr(2,:),'LineWidth',4)
hold on
plot(p3,1-squeeze(chi3(i,j,:,l,m,n,o)),'Color',clr(3,:),'LineWidth',4)
hold on
plot(p5,1-squeeze(chi3(i,j,k,l,:,n,o)),'Color',clr(5,:),'LineWidth',4)
hold on
plot(p6,1-squeeze(chi3(i,j,k,l,m,:,o)),'Color',clr(6,:),'LineWidth',4)
hold on
plot(p7,1-squeeze(chi3(i,j,k,l,m,n,:)),'Color',clr(7,:),'LineWidth',4)
hold on
plot(p4,1-squeeze(chi3(i,j,k,:,m,n,o)),'Color',clr(4,:),'LineWidth',4)
legend('\mu_v','\sigma_v','\mu_n','\sigma_n','P','\tau','\Delta')
set(gca,'FontSize',24)
xlabel('Parameter','FontSize',30)
ylabel('1-R^2','FontSize',30)
set(gca,'YScale','log')

end
% save matricies to current folder if no directory path is given in
% filename
save(filename,'chi1','chi2','chi3','p1','p2','p3','p4','p5','p6','p7')


