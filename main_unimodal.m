function main_unimodal(option,varargin)

switch option
    case 'plot_only'
        plot_unimodal(varargin{1},[3]);
        return
    case 'analysis_only'
        post_process_unimodal(varargin{1});
        return
    otherwise
        gen_raw_data_unimodal;
end
end


function gen_raw_data_unimodal(~)

p.ntrials = 8*3; %for parallel processing,use multiples of 8
p.target = 'unimodal';
p.sigma = 1;
p.dt = 1e-3; %1e-2 is no good for calculating Riesz derivative
T = 1e3;

%AIM: 1. Sample trajectory. 2. histogram 2. ACF. 3. MSE
%For FHMC,HMC,

a = [1.2 2];
solver = {'Hamiltonian2','Langevin'}; %'Hamiltonian','Underdamped'};
%for accurate results don't need the other two... too slow

X = zeros(T/p.dt,p.ntrials,length(solver),length(a));

disp('Generating raw data...')
for i = 1:length(a) %group fractional/nonfractional together as this is less familiar to people
    aa = a(i);
    tic
    for j = 1:length(solver)
        p.methods.solver = solver{j};
        temp = zeros(T/p.dt,p.ntrials);
        parfor k = 1:p.ntrials
            [temp(:,k),t] = fHMC(T,aa,p);
        end
        clc
        fprintf('Solver: %u/%u\n',j,length(solver))
        fprintf('alpha: %u/%u\n',i,length(a))
        toc
        X(:,:,j,i) = temp;
    end
    
    
    
end

n = floor(T/p.dt); %number of samples
t = (0:n-1)*p.dt;

save main_unimodal_raw_data.mat T a solver X t p
end

function post_process_unimodal(Z)
n=50;
binedge = linspace(-5,5,n+1);
bin = binedge(1:end-1) + 0.5*(binedge(2)-binedge(1));
disp('Evaluating histogram, ACF, and MSE')
H = zeros(n,length(Z.solver),length(Z.a));

m = 10/Z.p.dt;%1e4;
tau = (0:m)*Z.p.dt;
ACF = zeros(length(tau),length(Z.solver),length(Z.a));
MSE = zeros(length(tau),length(Z.solver),length(Z.a));

t = (0:5e4)*Z.p.dt;
X = zeros(length(t),length(Z.solver),length(Z.a)); %sample X for visualization

tic
for i = 1:length(Z.a)
    for j = 1:length(Z.solver)
        x = Z.X(:,:,j,i); x = x(:);
        H(:,j,i) = histcounts(x,binedge,'normalization','pdf');
        temp_acf = 0;
        temp_mse = 0;
        for k = 1:Z.p.ntrials
            temp_acf = temp_acf + autocorr(Z.X(:,k,j,i),m);
            temp_mse = temp_mse + myMSE(Z.X(:,k,j,i),m);
        end
        ACF(:,j,i) = temp_acf/Z.p.ntrials;
        MSE(:,j,i) = temp_mse/Z.p.ntrials;
        X(:,j,i) = Z.X(1:length(t),1,j,i);
        clc
        toc
        fprintf('%u / %u\n',i,length(Z.a))
        fprintf('%u / %u\n',j,length(Z.solver))
    end
end
save main_unimodal_post_process.mat H ACF MSE tau bin binedge t X
end


function plot_unimodal(Q,index)

close all
if any(index == 0) %sample trajectory
    f0 = myfigure;
    ttl = {'Fractional Hamiltonian Monte Carlo','Fractional Langevin Monte Carlo','Hamiltonian Monte Carlo','Langevin Monte Carlo'};
    for k = 1:4
        subplot(4,1,k)
        [j,i]=ind2sub([2 2],k);
        plot(Q.t(1:10:end),Q.X(1:10:end,j,i),'.','color',mycolor(1,k),'markersize',1);
        title(ttl(k))
        if k<4
            set(gca,'xtick',[],'XColor','none');
        else
            xlabel('t')
        end
        ylabel('x_t')
        %subplotmod;
        box off
        set(gca,'TickDir','out')
        %hold on
        %plot(xx,yy,'k--','linewidth',1.5)        
        ylim([-4 4])
        xlim([Q.t(1) Q.t(end)])
    end
    pos =get(f0,'Position');
    set(gcf,'position',[pos(1) pos(2) pos(3) 600])
    export_fig(f0,'figures/fig_main_unimodal_postprocess_sample_path.pdf','-pdf','-nocrop','-transparent','-painters');
end

if any(index == 1)    %ACF
    f1 = myfigure;    
    % plot(Q.tau,Q.ACF(:,1,1),'-','color',mycolor(1,2),'linewidth',1.5)
    % hold on
    % plot(Q.tau,Q.ACF(:,2,1),'-','color',mycolor(1,1),'linewidth',1.5)
    % plot(Q.tau,Q.ACF(:,1,2),'--','color',mycolor(1,2) ,'linewidth',1.5)
    % plot(Q.tau,Q.ACF(:,2,2),'--','color',mycolor(1,1),'linewidth',1.5)
    
    %color scheme default
    plot(Q.tau,Q.ACF(:,1:2,1),Q.tau,Q.ACF(:,1:2,2),'--','linewidth',1.5)
    hold on
    plot(Q.tau([1 end]),[0 0],'color',[0.4 0.4 0.4],'HandleVisibility','off')
    
    xlim([0 8])
    xlabel('Time Lag')
    ylabel('Autocorrelation function')
    legend('FHMC','FLMC','HMC','LMC','box','off')
    subplotmod;
    offsetAxes;
    
    export_fig(f1,'figures/fig_main_unimodal_postprocess_acf.pdf','-pdf','-nocrop','-transparent','-painters');
end
if any(index == 2)
    f2 = myfigure;
    gs = makedist('normal');
    xx = linspace(-3,3,51);
    yy = pdf(gs,xx);
    
    for k = 1:4
        subplot(4,1,k)
        [j,i]=ind2sub([2 2],k);
        bar(Q.bin,Q.H(:,j,i),1,'FaceColor',mycolor(1,k),'EdgeColor',mycolor(1,k),'FaceAlpha',0.5);
        %bar(bins,hc,1,'FaceColor',mycolor(2),'EdgeColor',mycolor(0));hold on
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        hold on
        plot(xx,yy,'k--','linewidth',1.5)        
        ylim([0 0.5])
        xlim([-4 4])
        set(gcf,'position',[100 100 200 600])
    end
    %export_fig(f2,'figures/fig_main_unimodal_postprocess_histogram.pdf','-pdf','-nocrop','-transparent','-painters');
    print(gcf, '-dpdf', 'figures\fig_main_unimodal_postprocess_histogram.pdf'); 
end

if any(index == 3)
    f3 = myfigure;
    plot(Q.tau,Q.MSE(:,1:2,1)./Q.MSE(1,1:2,1),Q.tau,Q.MSE(:,1:2,2)./Q.MSE(1,1:2,2),'--','linewidth',1.5)
    hold on
    plot(Q.tau([1 end]),[0 0],'color',[0.4 0.4 0.4],'HandleVisibility','off')
    
    xlim([0 8])
    xlabel('Time Lag')
    ylabel('Mean Squared Error')
    legend('FHMC','FLMC','HMC','LMC','box','off')
    subplotmod;
    offsetAxes;
    
    export_fig(f3,'figures/fig_main_unimodal_postprocess_mse.pdf','-pdf','-nocrop','-transparent','-painters');
end

end