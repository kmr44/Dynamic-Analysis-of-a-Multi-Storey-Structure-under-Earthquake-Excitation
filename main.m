clear; clc; close all;

%Input
N=5;                                %Number of stories
Ks=1.44e8;                           %Stiffness
Ms=2.0e8;                           %Mass of each floor
floor_height = 3;                   

%Using Ground Motion Data (Acceleration(m/s2) vs Time)
data = readmatrix('Dataset.csv');   
t = data(:,1);                       
ag = data(:,2);                      

dt = mean(diff(t));   % Time step (s)
fs = 1/dt;            % Sample Frequency (Hz)

%Plot 1: Input Ground Acceleration (from CSV)
figure;
plot(t, ag);
xlabel('Time (s)'); ylabel('a_g (m/s^2)');
title('Input Ground Acceleration (from CSV)');
grid on;


M=Ms*diag(ones(N,1));                  %Mass matrix
K=MultiStory_Stiffness(Ks,N);          %Stiffness matrix

%Damping matrix
[V,D] = eig(K,M);
wn = sqrt(diag(D));                % natural frequencies (rad/s)
w1 = wn(1); w2 = wn(min(2,N));
zeta = 0.05;                       % 5% damping
A = [1/(2*w1), w1/2; 1/(2*w2), w2/2];
x = A\[zeta; zeta];
alpha = x(1); beta = x(2);
C = alpha*M + beta*K;

% Force matrix
ag = ag(:)';                      
f = zeros(N, length(ag));         
for i = 1:N
    f(i,:) = -Ms * ag;            
end

Result=MDOF_simulation(M,C,K,f,fs);    %MDOF solver

%Acceleration, Velocty, Time Response Plot
t=[0:1/fs:(length(ag)-1)/fs];
figure;
for i=1:1:N   
    subplot(N,3,3*i-2)
    plot(t,Result.Acceleration(i,:)); xlabel('Time(s)'); ylabel(strcat('Acc',num2str(i)));
    subplot(N,3,3*i-1)
    plot(t,Result.Velocity(i,:)); xlabel('Time(s)'); ylabel(strcat('Vel',num2str(i)));
    subplot(N,3,3*i)
    plot(t,Result.Displacement(i,:)); xlabel('Time(s)'); ylabel(strcat('Disp',num2str(i)));        
end
    
%Modeshapes Plot
figure;
xscale = 10.0e5;
yscale = 25;
for i=1:1:N
    plot([0 ;xscale*Result.Parameters.ModeShape(:,i)],yscale*(0:N),'*-');
    hold on
end
hold off
xlabel('Amplitude');
ylabel('Floor');
grid on
daspect([1 1 1]);
title('Modeshapes');



%Interstory Drift
u = Result.Displacement;                     
interstory_drift = diff(u,1,1);              
max_drift = max(abs(interstory_drift),[],2); 
drift_ratio = max_drift / floor_height;      

%Drift Ratios Plot 
figure;
bar(1:length(drift_ratio), drift_ratio);
xlabel('Story Level');
ylabel('Max Drift Ratio');
title('Interstory Drift Ratio per Story');
grid on;



a = Result.Acceleration;           % (N x time) floor accelerations
m = Ms * ones(N,1);                % (N x 1) mass per floor
base_shear = (m' * a);             % (1 x time) total shear at base
base_shear = base_shear(:)';       % ensure row vector for plotting
max_base_shear = max(abs(base_shear));

%Base Shear Time History Plot 
figure;
plot(t(1:length(base_shear)), base_shear, 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Base Shear (N)');
title('Base Shear Time History');
grid on;

disp(['Maximum Base Shear = ', num2str(max_base_shear), ' N']);




function Result=MDOF_simulation(M,C,K,f,fs)

n=size(f,1);
dt=1/fs; %sampling rate
[Vectors, Values]=eig(K,M);
Freq=sqrt(diag(Values))/(2*pi); % undamped natural frequency
steps=size(f,2);

Mn=diag(Vectors'*M*Vectors); % uncoupled mass
Cn=diag(Vectors'*C*Vectors); % uncoupled damping
Kn=diag(Vectors'*K*Vectors); % uncoupled stifness
wn=sqrt(diag(Values));
zeta=Cn./(2*sqrt(Mn.*Kn));   % damping ratio
wd=wn.*sqrt(1-zeta.^2);

fn=Vectors'*f; % generalized input force matrix

t=[0:dt:dt*steps-dt];

%forced vibration

for i=1:1:n
    
    h(i,:)=(1/(Mn(i)*wd(i))).*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t);                                                                                            %transfer function of displacement
    hd(i,:)=(1/(Mn(i)*wd(i))).*(-zeta(i).*wn(i).*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t)+wd(i).*exp(-zeta(i)*wn(i)*t).*cos(wd(i)*t));                              %transfer function of velocity
    hdd(i,:)=(1/(Mn(i)*wd(i))).*((zeta(i).*wn(i))^2.*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t)-zeta(i).*wn(i).*wd(i).*exp(-zeta(i)*wn(i)*t).*cos(wd(i)*t)-wd(i).*((zeta(i).*wn(i)).*exp(-zeta(i)*wn(i)*t).*cos(wd(i)*t))-wd(i)^2.*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t)); %transfer function of acceleration
    
    qq=conv(fn(i,:),h(i,:))*dt;
    qqd=conv(fn(i,:),hd(i,:))*dt;
    qqdd=conv(fn(i,:),hdd(i,:))*dt;
    
    q(i,:)=qq(1:steps); % modal displacement
    qd(i,:)=qqd(1:steps); % modal velocity
    qdd(i,:)=qqdd(1:steps); % modal acceleration
    
    
end

x=Vectors*q;   %Displacement
v=Vectors*qd;  %Vecloity
a=Vectors*qdd; %Acceleration


% Free vibration

xi=zeros(n,1);  % Displacement initial condition
vi=zeros(n,1);  % Velocity initial condition

xno=Vectors'*M*xi./Mn;
vno=Vectors'*M*vi./Mn;

for i=1:1:n

AA=(vno(i)+xno(i).*zeta(i).*wn(i))./wd(i);
BB=xno(i);

qf(i,:)=exp(-zeta(i)*wn(i)*t).*(AA.*sin(wd(i)*t)+BB.*cos(wd(i)*t));
qdf(i,:)=wd(i)*exp(-zeta(i)*wn(i)*t).*(AA.*cos(wd(i)*t)-BB.*sin(wd(i)*t))-zeta(i).*wn(i).*exp(-zeta(i)*wn(i)*t).*(AA.*sin(wd(i)*t)+BB.*cos(wd(i)*t));
qddf(i,:)=wd(i)^2*exp(-zeta(i)*wn(i)*t).*(-AA.*sin(wd(i)*t)-BB.*cos(wd(i)*t))-2*zeta(i).*wn(i).*wd(i).*exp(-zeta(i)*wn(i)*t).*(AA.*cos(wd(i)*t)-BB.*sin(wd(i)*t))+zeta(i)^2.*wn(i)^2*exp(-zeta(i)*wn(i)*t).*(-AA.*sin(wd(i)*t)-BB.*cos(wd(i)*t));


end

x=x+Vectors*qf;
v=v+Vectors*qdf;
a=a+Vectors*qddf;

Result.Displacement=x;
Result.Velocity=v;
Result.Acceleration=a;
Result.Parameters.Freq=Freq;
Result.Parameters.DampRatio=zeta*100;
Result.Parameters.ModeShape=Vectors;
end




function K=MultiStory_Stiffness(Ks,N)

KK=Ks*ones(N+1,1);
KK(N+1)=0;
K=zeros(N,N);

for i=1:1:N
for j=1:1:N              
if i==j 
K(i,j)=KK(i)+KK(i+1);
end            
if i==j+1
K(i,j)=-KK(i);
end      
if j==i+1
K(i,j)=-KK(j);
end        
end
end

end