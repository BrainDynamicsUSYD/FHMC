function T = TrapTime(t)
%INPUT: t = binary variable
%OUTPUT: Time spent in one mode
%t = x < 0;
dt = diff(t);
t_in = find(dt == 1);
t_out = find(dt == -1);

T = zeros(size(t_in));
for j =1 :length(T)
    if t_in(j)==1 %ignore influence from initial condition
        T(j)=NaN;
        continue
    end
    
    jj = find(t_out>t_in(j),1);
    if isempty(jj)  %ignore terminal condition
        T(j) = NaN;
    else
        T(j) = t_out(jj)-t_in(j);
    end
end
T(isnan(T))=[];