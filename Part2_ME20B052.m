%Project Part 2 - ME20B052 Sneha M S
% g,h,k --> direction of centre of Earth, Star A and Star B from the body frame of the satellite respectively
% <variable>_wn --> data without sensor noise

%% Modelling the System
J=[480 0 0; 0 640 0; 0 0 720]; %Moment of Inertia Tensor
n=zeros(3,1000); %Moments provided by Thrusters;
g=zeros(3,1000);
k=zeros(3,1000);
h=zeros(3,1000);
Q=zeros(4,1000); %Rotation Quaternion
w=zeros(3,1000); %Angular velocity
a=zeros(3,1000); %Angular Acceleration
g_wn=zeros(3,1000);
k_wn=zeros(3,1000);
h_wn=zeros(3,1000);
w_wn=zeros(3,1000);

%disp("Initialize the altitude quaternion");
%Q(1,1)=input("Enter the value of q0: ");
%Q(2,1)=input("Enter the value of q1: ");
%Q(3,1)=input("Enter the value of q2: ");
%Q(4,1)=input("Enter the value of q3: ");
Q(:,1)=[1;0;0;0];
Q(1,1)=Q(1,1);
Q(2,1)=Q(2,1);
Q(3,1)=Q(3,1);

%disp("Initialize the angular velocity");
%w(1,1)=input("Enter the value of wx: ");
%w(2,1)=input("Enter the value of wy: ");
%w(3,1)=input("Enter the value of wz: ");
w_wn(1,1)=1;
w_wn(2,1)=10;
w_wn(3,1)=1;

%disp("Enter the direction of Star A measured in ground frame");
%h1=input("Enter the value of hx: ");
%h2=input("Enter the value of hy: ");
%h3=input("Enter the value of hz: ");
h0=[0;1;0]; %With respect to ground frame
ru=rand;
rn=normrnd(0,1,3,1);
h_wn(:,1)=quat_to_matrix((Q(:,1))')*h0;
h(:,1)=(eye(3,3)-((0.01*(ru-1)).*crossprodv(w(:,1))))*quat_to_matrix((Q(:,1))')*(0.1.*rn);

%disp("Enter the direction of Star B measured in ground frame");
%k1=input("Enter the value of kx: ");
%k2=input("Enter the value of ky: ");
%k3=input("Enter the value of kz: ");
k0=[0;-1;0]; %With respect to ground frame
ru=rand;
rn=normrnd(0,1,3,1);
k_wn(:,1)=quat_to_matrix((Q(:,1))')*k0;
k(:,1)=(eye(3,3)-((0.01*(ru-1)).*crossprodv(w(:,1))))*quat_to_matrix((Q(:,1))')*k0+(0.1.*rn);

g0=[0;0;1]; %With respect to ground frame
ru=rand;
rn=normrnd(0,1,3,1);
g_wn(:,1)=quat_to_matrix((Q(:,1))')*g0;
g(:,1)=(eye(3,3)-((0.01*(ru-1)).*crossprodv(w(:,1))))*quat_to_matrix((Q(:,1))')*g0+(0.1.*rn);

Kp=[1000 0 0; 0 0 0; 0 0 1000];
w_cd=[0 5 0]';
Kd=2000;

dt=0.01;
total_time=10;
time=0;
for i =1:(total_time/dt)-1
    time=time+dt;
    rn=normrnd(0,3,1);
    w(:,i)=w_wn(:,i)+rn;
    %n(:,i+1)=-Kp*Q(2:4,i);
    %n(:,i+1)=-Kp*Q(2:4,i)-Kd.*(w(:,i)-w_cd);
    a(:,i)=J\(n(:,i+1)-cross(w_wn(:,i),(J*w_wn(:,i))));
    w_wn(:,i+1)=w_wn(:,i)+a(:,i).*dt;
    temp=[0,(w_wn(:,i+1))'];
    q_dot=(0.5.*quatprod(temp,(Q(:,i))'));
    Q(:,i+1)=Q(:,i)+(q_dot).*dt;
    Q(:,i+1)=Norm(Q(:,i+1));
    ru=rand;
    rn=normrnd(0,1,3,1);
    C=quat_to_matrix((Q(:,i+1))');
    h(:,i+1)=(eye(3,3)-((0.01*(2*ru-1)).*crossprodv(w(:,1))))*C'*h0+(0.1.*rn);
    ru=rand;
    rn=normrnd(0,1,3,1);
    k(:,i+1)=(eye(3,3)-((0.01*(2*ru-1)).*crossprodv(w(:,1))))*C'*k0+(0.1.*rn);
    ru=rand;
    rn=normrnd(0,1,3,1);
    g(:,i+1)=(eye(3,3)-((0.01*(2*ru-1)).*crossprodv(w(:,1))))*C'*g0+(0.1.*rn);
    h_wn(:,i+1)=uni(C'*h0);
    k_wn(:,i+1)=uni(C'*k0);
    g_wn(:,i+1)=uni(C'*g0);
    if(time==0.1)
        break;
    end
end

% Plotting Graphs
axis=0:dt:total_time;
axis(:,201)=[ ];
t=tiledlayout(2,2);
title(t,"Sensor Outputs");

nexttile
plot(axis,g,axis,g_wn);
title('Home Tracker')
xlabel("Time (in s)")
legend('ĝx','ĝy','ĝz','gx','gy','gz')

nexttile
plot(axis,h,axis,h_wn);
title('Star A Tracker')
xlabel("Time (in s)")
legend('ĥx','ĥy','ĥz','hx','hy','hz')

nexttile
plot(axis,k,axis,k_wn);
title('Star B Tracker')
xlabel("Time (in s)")
legend('ǩx','ǩy','ǩz','kx','ky','kz')

nexttile
plot(axis,w);
hold on
plot(axis,w_wn);
title('Gyroscope Readings')
xlabel("Time (in s)")
ylabel("ω(t) (in rad/s)")
legend('ŵx','ŵy','ŵz','ωx','ωy','ωz')

%% Attitude Estimation from Sensor Data - Triad Algorithm
%Home Tracker Sensor Output - g_wn
%Star A Tracker Sensor output - h_wn
%Measurements made in ground frame - g0 and h0
Q_triad=zeros(4,1000); 
G=[uni(g0),uni(cross(g0,h0)),uni(cross(g0,cross(g0,h0)))];
for i=1:(total_time/dt)
    u=g_wn(:,i);
    v=h_wn(:,i);
    B=[uni(u),uni(cross(u,v)),uni(cross(u,cross(u,v)))];
    C_triad=G*transpose(B);
    Q_triad(:,i)=Norm(Matrix_to_Quat(C_triad));
end

%% Attitude Estimation from Sensor Data - Davenport's Q Method
% Vectors - h0,h_wn (Star A); g0,g_wn (Earth); k0,k_wn (Star B)
Q_qmethod=zeros(4,1000);
w1=2/3; w2=2/3; w3=2/3;
v1_N=g0;
v2_N=h0;
v3_N=k0;
for i=1:(total_time/dt)
    v1_B=g_wn(:,i);
    v2_B=h_wn(:,i);
    v3_B=k_wn(:,i);
    B=w1.*(v1_B*v1_N')+w2.*(v2_B*v2_N')+w3.*(v3_B*v3_N');
    sigma=trace(B);
    S=B+B';
    Z=[ B(2,3) - B(3,2) ;
    B(3,1) - B(1,3) ;
    B(1,2) - B(2,1) ];
    K=[ sigma Z'; Z (S-sigma*eye(3))];
    [eigvec,eigval]=eig(K);
    max=eigval(1,1);
    ind=1;
    for j=2:4
        if(eigval(j,j)>max)
            max=eigval(j,j);
            ind=j;
        end
    end
    Q_qmethod(:,i)=eigvec(:,ind);
end
%% Functions 
%To find the rotation Matrix corresponding to the Quaternion
function C = quat_to_matrix(q)
    C(1)=(q(1)^2)+(q(2)^2)-(q(3)^2)-(q(4)^2);
    C(2,2)=(q(1)^2)-(q(2)^2)+(q(3)^2)-(q(4)^2);
    C(3,3)=(q(1)^2)-(q(2)^2)-(q(3)^2)+(q(4)^2);
    C(1,2)=2*((q(2)*q(3))-(q(4)*q(1)));
    C(2)=2*((q(2)*q(3))+(q(4)*q(1)));
    C(3)=2*((q(2)*q(4))+(q(3)*q(1)));
    C(1,3)=2*((q(2)*q(4))-(q(3)*q(1)));
    C(2,3)=2*((q(3)*q(4))-(q(2)*q(1)));
    C(3,2)=2*((q(3)*q(4))+(q(2)*q(1)));
end

%To find the product of two quaternions
function prod = quatprod(p,q)
    p0=p(1)*q(1)-p(2)*q(2)-p(3)*q(3)-p(4)*q(4);
    p1=p(1)*q(2)+p(2)*q(1)-p(3)*q(4)+p(4)*q(3);
    p2=p(1)*q(3)+p(2)*q(4)+p(3)*q(1)-p(4)*q(2);
    p3=p(1)*q(4)-p(2)*q(3)+p(3)*q(2)+p(4)*q(1);
    prod=[p0;p1;p2;p3];
end

%To Find the Cross Product Vector
function V = crossprodv(u)
    V(1,1)=0;
    V(1,2)=-u(3);
    V(1,3)=u(2);
    V(2,1)=u(3);
    V(2,2)=0;
    V(2,3)=-u(1);
    V(3,1)=-u(2);
    V(3,2)=u(1);
    V(3,3)=0;
end

%To normalize a quaternion
function n=Norm(q)
    sum=sqrt(q(1)^2+q(2)^2+q(3)^2+q(4)^2);
    n=q./sum;
end

%Obtaining a Unit Vector
function a = uni(b)
norm=0;
for i=1:3
    norm=norm+(b(i)^2);
end
norm=sqrt(norm);
a=b./norm;
end

%Obtaining Magnitude of a Vector
function a = mag(b)
norm=0;
for i=1:3
    norm=norm+(b(i)^2);
end
norm=sqrt(norm);
a=norm;
end

function q=Matrix_to_Quat(C)
    q(1)=0.25*sqrt(1+C(1,1)+C(2,2)+C(3,3));
    q(2)=0.25*sqrt(1+C(1,1)-C(2,2)-C(3,3));
    q(3)=0.25*sqrt(1-C(1,1)+C(2,2)-C(3,3));
    q(4)=0.25*sqrt(1-C(1,1)-C(2,2)+C(3,3));    
end

function R = rqprod(v)
    R=[0 v'; -v -crossprodv(v)];
end

function L = lqprod(u)
    L=[0 -u'; u -crossprodv(u)];
end

function q_1 = qinverse(q)
    q_1(1,1)=q(1,1);
    q_1(2,1)=-q(2,1);
    q_1(3,1)=-q(3,1);
    q_1(4,1)=-q(4,1);
end

function A = state_transition_Matrix(w,dt)
    A(1,1)=1;
    A(1,2)=-w(1,1)*dt/2;
    A(1,3)=-w(2,1)*dt/2;
    A(1,4)=-w(3,1)*dt/2;
    A(2,1)=w(1,1)*dt/2;
    A(2,2)=1;
    A(2,3)=-w(3,1)*dt/2;
    A(2,4)=w(2,1)*dt/2;
    A(3,1)=w(2,1)*dt/2;
    A(3,2)=w(3,1)*dt/2;
    A(3,3)=1;
    A(3,4)=-w(1,1)*dt/2;
    A(4,1)=w(3,1)*dt/2;
    A(4,2)=-w(2,1)*dt/2;
    A(4,3)=w(1,1)*dt/2;
    A(4,4)=1;
end