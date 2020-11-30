function main_bimodal(option,varargin)

switch option
    case 'plot_only'
        plot_bimodal(varargin{1},[2]);
        return
    case 'analysis_only'
        post_process_bimodal(varargin{1});
        return
    case 'sample_path'
        gen_sample_path_bimodal;
    otherwise
        exit_time;
end
end


function gen_sample_path_bimodal(~)

%p.ntrials = 8*3; %for parallel processing,use multiples of 8
p.target = 'bimodal';
p.sigma = 0.32;
p.dt = 1e-3; %1e-2 is no good for calculating Riesz derivative
T = 1e4;

%AIM: 1. Sample trajectory. 2. histogram

a = [1.2 2];
p.location = 2.5; %modal location
solver = {'Hamiltonian2','Langevin'}; %'Hamiltonian','Underdamped'};
%for accurate results don't need the other two... too slow

X = zeros(T/p.dt,length(solver),length(a)); %save 1 example is ok

disp('Generating raw data...')
for i = 1:length(a) %group fractional/nonfractional together as this is less familiar to people
    aa = a(i);
    tic
    for j = 1:length(solver)
        p.methods.solver = solver{j};
        
        disp('Simulating...')
        %parfor k = 1:p.ntrials
        [temp,t] = fHMC(T,aa,p);
        %end
        
        clc
        fprintf('Solver: %u/%u\n',j,length(solver))
        fprintf('alpha: %u/%u\n',i,length(a))
        toc
         
        X(:,j,i) = temp;
    end
   
end

n = floor(T/p.dt); %number of samples
t = (0:n-1)*p.dt;

save main_bimodal_raw_data_sample_path.mat T a solver X t p
end

function exit_time()

%x0 = linspace(0.5,2.5,21); %locations of bimodal distribution
%a = [1.2 2];
%p.T = 1e3;

x0 = linspace(0.5,1.5,21);
a = [2]; %do gaussian only
p.T = 1e4; %do more time since time is longer

solver = {'Hamiltonian2','Langevin'};
p.target = 'bimodal';
p.ntrials = 8*3;
p.dt = 1e-3;
p.sigma = 0.32;

flag = true(size(a));
tic



binedge = linspace(-5,5,101);
dx = (binedge(2)-binedge(1))*0.5;
bin = binedge(2:end) - dx;

H = zeros(length(bin), length(x0),length(solver),length(a)); %histogram
T = zeros(length(x0),length(solver),length(a),p.ntrials);

tic
disp('Beginning simulation...')
for m = 1:length(solver)
    p.methods.solver = solver{m};
    
    for i = 1:length(x0) %modal separation
        p.location = x0(i); %2.5;
        
        for k = 1:length(a)
            
            if ~flag(k)
                T(i,m,k,:) = NaN;
                continue
            end
            
            aa = a(k);
            TT = p.T;
            
            h = zeros(length(bin), p.ntrials);
            texit = zeros(1,p.ntrials);
            
            parfor j = 1:p.ntrials
                [X,t] = fHMC(TT,aa,p);
                h(:,j) = histcounts(X,binedge,'normalization','pdf');
                try
                    texit(j) = mean([TrapTime(X>0);TrapTime(X<0)]); %mean exit time
                catch ME
                    texit(j) = NaN;
                    disp('Unknown error!')
                end
            end
            H(:,i,m,k) = mean(h,2);
            T(i,m,k,:) = texit;
            
            if all(isnan(texit)) %T(k,i,j)>1e5 %if mean trap time is too long, skip subsequent simulations
                %nan => TrapTime output empty array
                flag(k) = false;
            end
            
            clc
            fprintf('Solver: %u/%u\n',m,length(solver))
            fprintf('alpha: %u/%u\n',k,length(a))
            fprintf('Mode location: %u/%u\n',i,length(x0))
            toc
            
        end
        
    end
end

%save main_bimodal_mean_exit_time.mat H bin T p x0 a solver
save main_bimodal_mean_exit_time_gaussian.mat H bin T p x0 a solver

end


function plot_bimodal(Q,flag)

close all
clc

if any(flag==0)
load main_bimodal_mean_exit_time_gaussian.mat T x0 p
%load exit_time_HMC_vs_FHMC_corrected.mat T x0
T(T==0)=nan;
T = mean(T,4)*p.dt;
x0 = x0*2;x0=x0(:);
xx = linspace(x0(1),x0(end),101);

f0 = myfigure;
for i =2:-1:1
    indx = ~isnan(T(:,i));
    linestyle = {'--','-'};
    markerstyle = {'x','o'};
    try
        f = fit(x0(indx),T(indx,i),'exp1');
        yy = f(xx);        
        plot(xx(yy<100),yy(yy<100),linestyle{i},'color',mycolor(1,i+2),'HandleVisibility','off','linewidth',1)
        hold on
    catch
    end
    plot(x0(indx),T(indx,i),markerstyle{i},'color',mycolor(1,i+2),'linewidth',1);
end

load main_bimodal_mean_exit_time.mat T x0 p
T(T==0)=nan;
T = mean(T,4)*p.dt;
x0 = x0*2;x0=x0(:);
xx = linspace(x0(1),x0(end),101);

for i =2:-1:1
    f = fit(x0(:),T(:,i,1),'poly1');
    linestyle = {'--','-'};
    markerstyle = {'x','o'};
    plot(xx,f(xx),linestyle{i},'color',mycolor(1,i),'HandleVisibility','off','linewidth',1)
    hold on
    %indx = 1:2:length(x0);
    plot(x0,T(:,i,1),markerstyle{i},'color',mycolor(1,i),'linewidth',1);    
end


subplotmod;

xlim([1 5])
ylim([0 100])
xlabel('Modal Separation')
ylabel('Mean Exit Time')

legend('LD','HD','FLD','FHD','box','off','location','ne')
offsetAxes;
export_fig(f0,'figures/fig_main_bimodal_met.pdf','-pdf','-nocrop','-transparent','-painters');
end

if any(flag==1)
    %load main_bimodal_raw_data_sample_path.mat
    f1 = myfigure;
    ttl = {'Fractional Hamiltonian Dynamics','Fractional Langevin Dynamics','Hamiltonian Dynamics','Langevin Dynamics'};
    
    Tmax = 2e2;
    tindx = (1:40:(Tmax/Q.p.dt))+2e4;
    
    for k = 1:4
        subplot(4,1,k)
        [j,i]=ind2sub([2 2],k);
        plot(Q.t(tindx)-Q.t(tindx(1)),Q.X(tindx,j,i),'.','color',mycolor(1,k),'markersize',1);
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
        xlim([0 Tmax])
        set(gca,'linewidth',1);
        
    end
    pos =get(f1,'Position');
    set(gcf,'position',[pos(1) pos(2) pos(3) 600])
    
    export_fig(f1,'figures/fig_main_bimodal_sample_path.pdf','-pdf','-nocrop','-transparent','-painters');
end

if any(flag ==2) %Histogram
    f2 = myfigure;
    gs1 = makedist('normal','mu',2.5,'sigma',0.32);
    gs2 = makedist('normal','mu',-2.5,'sigma',0.32);
    
    xx = linspace(-4,4,101);
    yy = 0.5*pdf(gs1,xx)+0.5*pdf(gs2,xx);
    
    binedge = linspace(-4,4,101);
    bin = (binedge(2)-binedge(1))*0.5 + binedge(1:end-1);
    
    for k = 1:4
        subplot(4,1,k)
        [j,i]=ind2sub([2 2],k);
        hc = histcounts(Q.X(:,j,i),binedge,'normalization','pdf');
        %bar(Q.bin,Q.H(:,j,i),1,'FaceColor',mycolor(1,k),'EdgeColor',mycolor(1,k),'FaceAlpha',0.5);
        bar(bin,hc,1,'FaceColor',mycolor(1,k),'EdgeColor',mycolor(1,k),'FaceAlpha',0.5);hold on
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        hold on
        plot(xx,yy,'k','linewidth',1)        
        if k<0;%k <2.5
            ylim([0 1.5/2])
        else
            ylim([0 1.5])
        end
        xlim([-4 4])
        set(gcf,'position',[100 100 200 600])
    end
    %export_fig(f2,'figures/fig_main_unimodal_postprocess_histogram.pdf','-pdf','-nocrop','-transparent','-painters');
    print(gcf, '-dpdf', 'figures\fig_main_bimodal_histogram.pdf'); 
end

end