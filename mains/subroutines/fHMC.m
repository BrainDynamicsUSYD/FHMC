function [X,t] = fHMC(T,a,p)
%Levy Monte Carlo
%INPUT: T = time span (s), a = Levy characteristic exponent
% p = other parameters
%OUTPUT: X = samples, t = time
%
p.depth = 1/p.sigma/p.sigma;
p.rmax = 5; %cutoff
p.xmax = 5; %cutoff
p.x0 = 0;
p.r0 = 1;

%p.methods.drift = 'approx'; %exact (numerically unstable)


%
dt = p.dt;%1e-3; %integration time step (s)
dta = dt.^(1/a); %fractional integration step

n = floor(T/dt); %number of samples
t = (0:n-1)*dt;
x = p.x0; %initial condition
r = p.r0;
X = zeros(n,1); %position
%R = zeros(n,1); %momentum

switch p.target
    case 'unimodal'
        U = @(x) p.depth*x.*x/2;
        dU = @(x) p.depth*x;
    case 'bimodal'
        s1 = p.location;%2.5;
        s2 = -s1;
        w = 0.5;
        sigma2 = p.sigma.*p.sigma;%0.32;%*(0.32/0.54);
        U = @(x) -log( 1e-8 +w*exp(-(x-s1).*(x-s1)/2/sigma2) + (1-w)*exp(-(x-s2).*(x-s2)/2/sigma2))+0.5;
        dx = 1e-3;
        dU = @(x) (U(x+dx)-U(x-dx))/dx/2;
end

if a<2 %approximation of drift term
    ba = @(x) frac_drift(U,dU,x,a); %approximate full Riesz D with partial Dx                
else
    ba = @(x) -dU(x); %approximation
end

c = gamma(a-1)*(gamma(a/2).^(-2)); %<--should be this!
% if a < 2  %use this for a=1.2 as it seems to give better result???
%     c = gamma(a+1)*(gamma(a/2+1).^(-2)); %a = a-2;
% end
%tic
switch p.methods.solver
    case 'Hamiltonian'
        for i = 1:n
            while true
                %diffusion term
                g = stblrnd(a,0,1,0); %Levy r.v. wt expo a and scale 1                
                xnew = x + c*r*dt;
                rnew = r -c*dU(x)*dt - c*r*dt + g*dta;
                %accept criterion
                if abs(rnew)<p.rmax %5
                    x = xnew;
                    r = rnew;
                    break %accept
                end
            end
            X(i) = x;
        end
    case 'Langevin'
        for i = 1:n
            while true
                %diffusion term
                g = stblrnd(a,0,1,0); %Levy r.v. wt expo a and scale 1
                %NB. when a=2, g~sqrt(2)*gaussian r.v. as it should be
                xnew = x + ba(x)*dt + g*dta;
                %accept criterion
                if abs(xnew)<p.xmax %5
                    x = xnew;
                    break %accept
                end
            end
            X(i) = x;
        end
    case 'Hamiltonian2'
        for i = 1:n
            while true
                %diffusion term
                g = stblrnd(a,0,1,0); %Levy r.v. wt expo a and scale 1                
                xnew = x + ba(x)*dt + r*dt + g*dta;
                rnew = r + ba(x)*dt;
                %accept criterion
                if abs(xnew)<p.xmax %5
                    x = xnew;
                    r = rnew;
                    break %accept
                end
            end
            X(i) = x;
        end
    case 'Underdamped'
        sas = makedist('Stable','alpha',1/a,'beta',0,'gam',1,'delta',0);
        dv = 1e-3;
        G = @(v) (pdf(sas,v+dv)-pdf(sas,v-dv))/dv/2; %rough approximation
        for i = 1:n
            while true
                %diffusion term
                g = stblrnd(a,0,1,0); %Levy r.v. wt expo a and scale 1
                
                rnew = r - dU(x)*dt - r*dt + g*dta;
                xnew = x + G(rnew)*dt;
                
                %accept criterion
                if abs(rnew)<p.rmax %5
                    x = xnew;
                    r = rnew;
                    break %accept
                end
            end
            X(i) = x;
        end
        
end
%toc
