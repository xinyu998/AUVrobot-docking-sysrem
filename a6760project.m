%Localization of the AUV
%Particle f
%MAE 6760 Project
%Xinyu Liu
clear all;close all;
rng('default')

rng(100);
%
% Two cases to consider: select either baseline or swervy
scenario_type='path2';

[Uacc,Uomega]=get_controlinputs(scenario_type);

%state vector: x,y 2D position, velocity, heading
%simulate no noise system
[Xnonoise,n,t,dt,nt]=simulate_2Dcar(Uacc,Uomega);
nk=nt;
%
%plots
plot_birdseyeview(Xnonoise,[],[],'True trajectory View');
ii_plot=[3 4];
plot_estimator(t,Xnonoise,[],[],ii_plot,'Truth: Velocity/Heading States');
%% Create IMU measurements (accel, rate gyro) and 2D GPS

%Create accel/rate gyro measurements with bias
nw=2;
Qw=diag([0.1^2 0.04^2]);
w=sqrtm(Qw)*randn(nw,nt);
Zacc=Uacc+w(1,:);
Zrg=Uomega+w(2,:);

%Create 2D GPS like measurements 
nz=2;
R=eye(nz)*1^2;
v=sqrtm(R)*randn(nz,nt);
Z=[Xnonoise(1:2,:)] + v;
Hgps=[eye(2) zeros(2,2)]; %output matrix is linear
%% SENSOR MODEL: lane detector
%two mode Gaussian
mx=0;
sx=1;
%
pdfx = makedist('Normal',mx,sx^2);
%sensor likelihood model, assuming mx1 or mx2 are correct
Lx = makedist('Normal',0,sx^2);
%% SENSOR MODEL: Odometry
%one mode Gaussian
sy=1;
my=0;
pdfy = makedist('Normal',my,sy);
%sensor likelihood model, assuming my is correct
Ly = makedist('Normal',0,sy);
%% PF
%Initialize state using sensor model
ns=1000;
x0=random(pdfx,ns,1);
y0=random(pdfy,ns,1);

X=zeros(ns,2,nk);
X(:,:,1)=[x0 y0];
W=zeros(ns,nk);
W(:,1) = ones(ns,1)/ns;
xEst=zeros(nk,2);xEst(1,:)=mean(X(:,:,1),1);
xSig=zeros(nk,2);xSig(1,:)=std(X(:,:,1),1);

figure(10);
hold on;

hPoints = plot(X(:,1,1),X(:,2,1),'b.','markersize',10);
hEst = plot(xEst(1,1),xEst(1,2),'r*','markersize',10);
hold off;
title('Particle Distribution');


axis([-5 105 -50 5]);

Neff_thresh=ns; %will always resample (enables SIR)
%Neff_thresh=ns/2; %standard practice

for k=2:nk,
    
    Xprior=X(:,:,k-1);
    
    %PREDICTION STEP
    V=Xnonoise(3,k);
    theta=Xnonoise(4,k);
    u=ones(ns,1)*dt*[V*cos(theta) V*sin(theta)];%sample process noise
    
    Xpred=Xprior+u; %predition step for all particles: k -> k+1
    
    %UPDATE STEP
    Zcurrent=Z(:,k)'; %pull off current measurement
    Zhat = Xpred; %calculate estimated measurement: assumes H=eye(2) or z=x+v 
    
    Inn = ones(ns,1)*Zcurrent - Zhat; %calculate innovations for each particle
    
    %calculate likelihood weighting for each particle
    for ipart=1:ns,
        L(ipart,1) = pdf(Lx,Inn(ipart,1))*pdf(Ly,Inn(ipart,2)); %likelihoods of two indepedent msts multiply
    end
    
    %UPDATE WEIGHTS
    Wk_unnorm=W(:,k-1).*L; %update weight with unnormalized likelihood
    Wk=Wk_unnorm/sum(Wk_unnorm); %normalize all weights to sum to 1

    %RESAMPLING    
    Neff=1/[Wk'*Wk]; %calculate effective number of particles (varies between 1 and ns)
    if Neff<Neff_thresh, %resample if under threshold
        %
        CDF = cumsum(Wk)/sum(Wk); %create a running sum function of the weights
        CDF_plus=CDF+rand(ns,1)*1E-6; %for cases when a particle weight went to zero
        %randomly (uniform) choose likely (better) particles...
        iSelect  = rand(ns,1);
        %find particle corresponding to each y value
        iNextGeneration = interp1(CDF_plus,1:ns,iSelect,'nearest','extrap');
        %copy selected particles for next generation
        X(:,:,k) = Xpred(iNextGeneration,:);
        W(:,k) = ones(ns,1)/ns;
    else,
        %
        X(:,:,k) = Xpred;
        W(:,k)=Wk;
    end
        
    %our estimate is simply the mean of the particles
    xEst(k,:) = sum(Wk.*X(:,:,k),1);
    xSig(k,:) = sqrt([sum(Wk.*(X(:,1,k)-xEst(k,1)).^2,1) sum(Wk.*(X(:,2,k)-xEst(k,2)).^2,1)]);
    %for bootstrap/SIR (always resample), where all weights are the same
    %xEst(k,:) = mean(X(:,:,k),1);
    %xSig(k,:) = std(X(:,:,k),1);
    
%     figure(10);
%     set(hPoints,'XData',X(:,1,k));
%     set(hPoints,'YData',X(:,2,k));
%     set(hEst,'XData',xEst(k,1));
%     set(hEst,'YData',xEst(k,2));
%     drawnow;    
%     pause(0.1);
%     
%     if record_video,
%         F(k) = getframe(10); 
%         writeVideo(vidfile,F(k));
%     end
%P [2x2x350]
P = zeros(2,2,nk);
for i=1:nk
    Sx= xSig(i,1);
    Sy= xSig(i,2);
    P(:,:,i)= [Sx 0 ; 0 Sy];
end
ii_plot=[1 2];
plot_estimator_error(t,Xnonoise,xEst',P,ii_plot,'PF: North/East States (no bias estimation)');

plot_birdseyeview(Xnonoise,xEst',P,'PF: Birds Eye (NO-bias estimation)');

%% Simulation Function calls
function [Uacc,Uomega]=get_controlinputs(scenario_type);
%
%scenario_type='baseline';
%scenario_type='swervy';
%
%   generates clean acceleration and rotational rate control inputs
%
%
%define inputs: acceleration and heading rate
if strcmp(scenario_type,'baseline'),
    Uacc=[ones(1,60) zeros(1,100) -ones(1,60) zeros(1,20) ones(1,60) zeros(1,50)];
    Uomega=[zeros(1,240) -ones(1,40)*pi/2/4 zeros(1,70)];
elseif strcmp(scenario_type,'path2'),
    Uacc=[zeros(1,350)];
    Uomega=[zeros(1,350)];  
else,
    return;
end
%
end

function [Xnonoise,n,t,dt,nt]=simulate_2Dcar(Uacc,Uomega);
%
%   simulates a 2D car using acceleration and rotational rate control inputs
%   state vector: x,y 2D position, velocity, heading
%
dt=0.1;
nt=length(Uacc);
t=[0:dt:dt*(nt-1)];
n=4;
Xnonoise=zeros(n,nt); 
%random create nonoise data?
for k=1:(nt-1),
    Vk=Xnonoise(3,k); %Velocity 
    Tk=Xnonoise(4,k);   %Time 
    Xnonoise(:,k+1) = Xnonoise(:,k) +...
        dt*[Vk*cos(Tk);Vk*sin(Tk);Uacc(k);Uomega(k)];   
end
end
%% Plotting functions

function plot_estimator(t,x1,x2,P2,ii_plot,title_name);
% x1 is the true value or reference comparison
% x2,P2 is the estimator state and covariance
% ii_plot: 2x1 vector of which states to plot
%
axis_names={'North (m)','East (m)','Velocity (m/sec)','Heading (rad)','Accel Bias (m/sec^2)','RG Bias (rad/sec)'};
figure;subplot(122);
ii_x1=[];ii_x2=[];ii_P2=[]; %for legend
%
for i=1:length(ii_plot),
    ii=ii_plot(i);
    subplot(1,2,i);
    hold on;
    if ~isempty(x1),
        plot(t,x1(ii,:),'color',[0 0.5 0]);ii_x1=1;
    end  
    if ~isempty(x2),
        plot(t,x2(ii,:),'b-');ii_x2=2;
    end  
    if ~isempty(P2)
        plot(t,x2(ii,:)'-2*sqrt(squeeze(P2(ii,ii,:))),'b:');
        plot(t,x2(ii,:)'+2*sqrt(squeeze(P2(ii,ii,:))),'b:');ii_P2=3;
    end
    hold off
    xlabel('time (sec)');ylabel(axis_names(ii));grid;
    xlim([0 35]);set(gca,'xtick',[0:5:35]);
end
legend_names={'true state','estimate','2\sigma bound'};
legend(legend_names{ii_x1},legend_names{ii_x2},legend_names{ii_P2},'Location','South');
%
sgtitle(title_name);
PrepFigPresentation(gcf);
end

function plot_estimator_error(t,x1,x2,P2,ii_plot,title_name);
% x1 is the true value or reference comparison
% x2,P2 is the estimator state and covariance
% ii_plot: 2x1 vector of which states to plot
%
axis_names={'North (m)','East (m)','Velocity (m/sec)','Heading (rad)','Accel Bias (m/sec^2)','RG Bias (rad/sec)'};
figure;subplot(122);
%
for i=1:length(ii_plot),
    ii=ii_plot(i);
    subplot(1,2,i);
    err=x2(ii,:)-x1(ii,:);
    plot(t,err,'b-');
    hold on;
    if ~isempty(P2)
        plot(t,err'-2*sqrt(squeeze(P2(ii,ii,:))),'b:');
        plot(t,zeros(length(t),1),'r--');
        plot(t,err'+2*sqrt(squeeze(P2(ii,ii,:))),'b:');
    end
    hold off
    xlabel('time (sec)');ylabel(axis_names(ii));grid;
    xlim([0 35]);set(gca,'xtick',[0:5:35]);
end
legend('estimator error','2\sigma bound','zero error','Location','South');
%
sgtitle(title_name);
PrepFigPresentation(gcf);
end

function plot_birdseyeview(x1,x2,P2,title_name);
% x1 is the true value or reference comparison
% x2,P2 is the estimator state and covariance
% 
ii_x1=[];ii_x2=[];ii_P2=[]; %for legend
figure;
if ~isempty(x1),
    plot(x1(1,:),x1(2,:),'color',[0 0.5 0]);ii_x1=1;
end
hold on;
if ~isempty(x2),
    plot(x2(1,:),x2(2,:),'b-');ii_x2=2;
end
if ~isempty(P2),
    iell=[2 10:10:350];
    for i=1:length(iell),
        ii=iell(i);
        [Xe,Ye] = calculateEllipseCov(x2([1 2],ii),P2([1 2],[1 2],ii),3);
        plot(Xe,Ye,'m-');
    end
    ii_P2=3;
end
xlabel('North (m)');ylabel('East (m)');grid;
hold off;
legend_names={'true trajectory','estimated trajectory','3\sigma bound'};
legend(legend_names{ii_x1},legend_names{ii_x2},legend_names{ii_P2},'Location','South')
%
title(title_name);
PrepFigPresentation(gcf);
end

function [Xe,Ye] = calculateEllipseCov(X, P, nsig, steps) 
    %# This functions returns points to draw an ellipse 
    %# 
    %#  @param X     x,y coordinates 
    %#  @param P     covariance matrix 
    %# 
 
    error(nargchk(2, 3, nargin)); 
    if nargin<3, nsig = 1; end 
    if nargin<4, steps = 36; end 
    
    [U,S,V]=svd(P);
    s1=sqrt(S(1,1));s2=sqrt(S(2,2));angle=acos(U(1,1))*180/pi;
    x=X(1);
    y=X(2);

    %scale by nsig
    s1=nsig*s1;
    s2=nsig*s2;

    beta = angle * (pi / 180); 
    sinbeta = sin(beta); 
    cosbeta = cos(beta); 
 
    alpha = linspace(0, 360, steps)' .* (pi / 180); 
    sinalpha = sin(alpha); 
    cosalpha = cos(alpha); 
 
    Xe = x + (s1 * cosalpha * cosbeta - s2 * sinalpha * sinbeta); 
    Ye = y + (s1 * cosalpha * sinbeta + s2 * sinalpha * cosbeta); 
 
end 

function PrepFigPresentation(fignum);
%
% prepares a figure for presentations
%
% Fontsize: 14
% Fontweight: bold
% LineWidth: 2
% 

figure(fignum);
fig_children=get(fignum,'children'); %find all sub-plots

for i=1:length(fig_children),
    
    set(fig_children(i),'FontSize',16);
    set(fig_children(i),'FontWeight','bold');
    
    fig_children_children=get(fig_children(i),'Children');
    set(fig_children_children,'LineWidth',2);
end
end