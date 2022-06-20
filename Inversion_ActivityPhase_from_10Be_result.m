% Rock glacier age profiles
% Written by R. Anderson and modified by B.Lehmann 23 Feb 2022
% contact: lehmann.benj@gmai.com

clear all
clc
clf

%% Import data
[num,text]=xlsread('Cosmo_data.xlsx');
a_S  = num(1:19,3);
ea_S = num(1:19,4);
x_S  = num(1:19,2);

figure(1)

%% Choice of the constant

tmax = 12100; % [year] this is time into the simulation, max 10Be age 
inh = 2159.6; % [year] maximum inheritance from residence on headwall from linear regression
%inh = 10;
xinit = 100;  % [m] maximum range for rockfall onto rock glacier surface
%xinit = 1;  % [m] maximum range for rockfall onto rock glacier surface
t3 = 3400; 
t2m = 3600;

%% Parametrization of the inversion

dt = 200;
t = dt:dt:tmax;
imax = length(t);

TT = 100;

t4 = min(a_S(9:19))*1000;
t2 = t2m:(t4-t2m)/(TT-1):t4;

inh_vec = 0:inh/(TT-1):inh;

for j=1:TT

age = inh*rand(size(t)); % this is initial age
distance = xinit*rand(size(t));

age_save = age;
distance_save = distance;

speed = zeros(size(t));

v1 = max(x_S)/(tmax-t2(j));
v2 = max(x_S(1:8))/t3; 
T1 = tmax-t2(j);
T2 = tmax-t3;
speed(t<=T1) = v1; % [m/a] early in the simulation, so longer ago
speed(t>=T2) = v2; % [m/a] later in the simulation so closer to the present

stopped = find(speed==0);

for i = 1:imax

age(1:i) = age(1:i)+dt;

% need to prevent sites that have stopped from restarting
if(i>stopped(end))
    distance(stopped(end):i) = distance(stopped(end):i) + (speed(i)*dt);
else
    distance(1:i) = distance(1:i) + (speed(i)*dt);
end
end

subplot(1,3,1)
plot((age-inh)/1000,speed,'linewidth',1,'Color',[0.5,0.5,0.5])

hold on
axis([0 tmax/1000 0 max(speed)+0.05])

subplot(1,3,2)
plot(distance(t>T2),age(t>T2)/1000,'o','linewidth',3,'MarkerEdgeColor',[0.5,0.5,0.5])
hold on
maxx = max(distance(t>T2));

plot(distance(t<T1 & distance>maxx),age(t<T1 & distance>maxx)/1000,'+','linewidth',3,'MarkerEdgeColor',[0.5,0.5,0.5])

age_interp1 = interp1(distance(t>T2),age(t>T2)/1000,x_S);
age_interp2 = interp1(distance(t<T1 & distance>maxx),age(t<T1 & distance>maxx)/1000,x_S);

chi1(j) = nansum(((a_S-age_interp1).^2)/((var(a_S))^2));
chi2(j) = nansum(((a_S-age_interp2).^2)/((var(a_S))^2));

end

%% Chi square calculation

freedome= numel(a_S)-4;
chi1_v = (chi1./freedome);
chi2_v = (chi2./freedome);
chi1_norm = chi1_v./max(chi1_v);
chi2_norm = chi2_v./max(chi2_v);


M1=1./exp(0.5*chi1_norm);
M2=1./exp(0.5*chi2_norm);

M1_norm = M1./max(M1);
M2_norm = M2./max(M2);

M_combined = M1_norm.*M2_norm;
M_norm = M_combined./max(M_combined);

%% Resampling the likelihood by comparing to a random value between 0 and 1

jt=0;
for it=1:TT
    R=rand;
    if (M2(it)>R)
        jt        = jt+1;
        s_M(jt) = M_norm(it);
        s_M2(jt) = M2(it);
        s_t2(jt)= t2(it);
    end
end

[value_BF,index_BF]=max(s_M2);
%% extract 1d PDFs and confidence intervals

nbin             = 20;

[n,xout]         =   hist(s_t2,nbin);

xwork            =   cumsum(n/sum(n));

ix               =   find(xwork>0.175,1);
t_1sd            =   xout(ix);
t_1sd_pd         =   n(ix)/sum(n);

ix               =   find(xwork>0.025,1);
t_2sd            =   xout(ix);
t_2sd_pd         =   n(ix)/sum(n);

ix               =   find(xwork>0.825,1);
t_1su            =   xout(ix);
t_su_pd          =   n(ix)/sum(n);

ix               =   find(xwork>0.925,1);
t_2su            =   xout(ix);
t_2su_pd         =   n(ix)/sum(n);

ix               =   find(xwork>0.50,1);
t2_median         =   xout(ix);
t_median_pd      =   n(ix)/sum(n);

%% Calcul result for t_median 

speed_median = zeros(size(t));
age_median = age_save;
distance_median = distance_save;

v1_median = max(x_S)/(tmax-t2_median);
v2_median = max(x_S(1:8))/t3; 
T1_median = tmax-t2_median;
T2_median = tmax-t3;
speed_median(t<=T1_median) = v1_median; % [m/a] early in the simulation, so longer ago
speed_median(t>=T2_median) = v2_median; % [m/a] later in the simulation so closer to the present

stopped = find(speed_median==0);

for i_median = 1:imax

age_median(1:i_median) = age_median(1:i_median)+dt;

if(i_median>stopped(end))
    distance_median(stopped(end):i_median) = distance_median(stopped(end):i_median) + (speed_median(i_median)*dt);
else
    distance_median(1:i_median) = distance_median(1:i_median) + (speed_median(i_median)*dt);
end
end
maxx_median = max(distance_median(t>T2_median));

%% Calcul result for t Best Fit 

speed_BF = zeros(size(t));
age_BF = age_save;
distance_BF = distance_save;

v1_BF = max(x_S)/(tmax-t2(index_BF));
v2_BF = max(x_S(1:8))/t3; 
T1_BF = tmax-t2(index_BF);
T2_BF = tmax-t3;
speed_BF(t<=T1_BF) = v1_BF; % [m/a] early in the simulation, so longer ago
speed_BF(t>=T2_BF) = v2_BF; % [m/a] later in the simulation so closer to the present

stopped = find(speed_BF==0);

for i_BF = 1:imax

age_BF(1:i_BF) = age_BF(1:i_BF)+dt;

if(i_BF>stopped(end))
    distance_BF(stopped(end):i_BF) = distance_BF(stopped(end):i_BF) + (speed_BF(i_BF)*dt);
else
    distance_BF(1:i_BF) = distance_BF(1:i_BF) + (speed_BF(i_BF)*dt);
end
end
maxx_BF = max(distance_BF(t>T2_BF));

%%
speed_1su = zeros(size(t));
age_1su = age_save;
distance_1su = distance_save;

v1_1su = max(x_S)/(tmax-t_1su);
v2_1su = max(x_S(1:8))/t3; 
T1_1su = tmax-t_1su;
T2_1su = tmax-t3;
speed_1su(t<=T1_1su) = v1_1su; % [m/a] early in the simulation, so longer ago
speed_1su(t>=T2_1su) = v2_1su; % [m/a] later in the simulation so closer to the present

stopped = find(speed_1su==0);

for i_1su= 1:imax

age_1su(1:i_1su) = age_1su(1:i_1su)+dt;

if(i_1su>stopped(end))
    distance_1su(stopped(end):i_1su) = distance_1su(stopped(end):i_1su) + (speed_1su(i_1su)*dt);
else
    distance_1su(1:i_1su) = distance_1su(1:i_1su) + (speed_1su(i_1su)*dt);
end
end
maxx_1su = max(distance_1su(t>T2_1su));

%%
speed_1sd = zeros(size(t));
age_1sd = age_save;
distance_1sd = distance_save;

v1_1sd = max(x_S)/(tmax-t_1sd);
v2_1sd = max(x_S(1:8))/t3; 
T1_1sd = tmax-t_1sd;
T2_1sd = tmax-t3;
speed_1sd(t<=T1_1sd) = v1_1sd; % [m/a] early in the simulation, so longer ago
speed_1sd(t>=T2_1sd) = v2_1sd; % [m/a] later in the simulation so closer to the present

stopped = find(speed_1sd==0);

for i_1sd= 1:imax

age_1sd(1:i_1sd) = age_1sd(1:i_1sd)+dt;

if(i_1sd>stopped(end))
    distance_1sd(stopped(end):i_1sd) = distance_1sd(stopped(end):i_1sd) + (speed_1sd(i_1sd)*dt);
else
    distance_1sd(1:i_1sd) = distance_1sd(1:i_1sd) + (speed_1sd(i_1sd)*dt);
end
end
maxx_1sd = max(distance_1sd(t>T2_1sd));

%% Plotting

subplot(1,3,1)
plot((age_median-inh)/1000,speed_median,'g','linewidth',3)

xlabel('Time [ka]','fontname','arial','fontsize',18)
ylabel('Rock glacier speed [m/a]','fontname','arial','fontsize',18)
set(gca,'fontsize',18,'fontname','arial') 
    
subplot(1,3,2)
errorbar(x_S(1:8),a_S(1:8),ea_S(1:8),ea_S(1:8),'rs','MarkerFaceColor','r')
errorbar(x_S(9:19),a_S(9:19),ea_S(9:19),ea_S(9:19),'bs','MarkerFaceColor','b')

plot(distance_median(t>T2_median),age_median(t>T2_median)/1000,'gd','MarkerFaceColor','g')
plot(distance_median(t<T1_median & distance_median>maxx_median),age_median(t<T1_median & distance_median>maxx_median)/1000,'gd','MarkerFaceColor','g') 

xlabel('Horizontal distance from the headwall [m]','fontname','arial','fontsize',18)
ylabel('Age including inheritance [ka]','fontname','arial','fontsize',18)
set(gca,'fontsize',18,'fontname','arial')
axis([0 1800 0 15])

subplot(1,3,3)

t2_median_vec = t2_median/1000*ones(100,1);
t2_1su_vec = t_1su/1000*ones(100,1);
t2_1sd_vec = t_1sd/1000*ones(100,1);
t2_BF_vec = t2(index_BF)/1000*ones(100,1);

likeH        = 0:1/(100-1):1;


plot(t2_median_vec,likeH,'g','linewidth',3)
hold on
plot(s_t2/1000,s_M,'k')
plot(t2_1su_vec,likeH,'g--','linewidth',3)
plot(t2_1sd_vec,likeH,'g--','linewidth',3)
% plot(t2_BF_vec,likeH,'k','linewidth',3)


axis([(min(t2)/1000)-0.5 (max(t2)/1000)+0.5 min(s_M2)-0.01 1])
xlabel('Explored time value [a]','fontname','arial','fontsize',18)
ylabel('Normalized likelihood','fontname','arial','fontsize',18)
legend('Median value ±1sigma','fontname','arial','fontsize',12, 'Location', 'Southeast')

%%
disp('Median value')
disp(['From ' num2str(tmax/1000,4) ' to ' num2str(t2_median/1000,3) ' ± ' num2str((abs(t2_median-t_1su))/1000,3) ...
    ' ka, surface velocity = ' num2str(v1_median,3) ' ± ' num2str(max(abs(v1_median-v1_1su),abs(v1_median-v1_1sd)),3) ' m/a'])
disp(['From ' num2str(t3/1000,4) ' ka, surface velocity = ' num2str(v2_median,2) ' m/a'])


