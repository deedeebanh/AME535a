clc
clear


%% User Inputs

nx=11;
ny=11;
re=100;
xmax=50;
h=1;
itmax=50;
uin=1;
dt=0.01;

%%
dx=xmax/(nx-1);
dy=h*2/(ny-1);
c=dt*uin/dx;
while (c>1)
    display('Error! CFL Condition > 1');
    break;
end
u=zeros(ny,nx);
v=zeros(ny,nx);
p=zeros(ny,nx);

x=0:dx:xmax;
y=-1:dy:1;

%Exact Solution
[uex qex]=exactvelocity(nx,ny);

%Initiate velocity u-component
u=uex;

%%
firstiter=1;
iter=0;

while (iter<=itmax)
    % Non-linear step
     if (firstiter==1)
        [u,v,fu,fv] = euler(u,v,nx,ny,dx,dy,re,dt,uin);
         firstiter=firstiter+1;
     else
         [u,v,fu,fv] = CONVEC(u,v,nx,ny,dx,dy,re,dt,fu,fv,uin);
     end
    % Pressure
    [p u v]=PRES(u,v,p,nx,ny,dx,dy,dt,uin,iter);
    % Viscous step
     u=VISC(u,nx,ny,dx,dy,re,dt,uin,1);
     v=VISC(v,nx,ny,dx,dy,re,dt,uin,2);
    %--------------------------------------
    iter=iter+1;
end

%%
%RMS
%For u-component
sum=0;
for k=2:ny-1
    error=uex(k)-u(k,nx);
    sum=sum+(error*error);
end
nn=(ny+1)/2;
an=ny-2;
u_rms=sqrt(sum/an);

%For the flow rate
i=0;
ii=0;
q=zeros(nx);
for k=1:nx
    ii=ii+1;
    for j=1:ny
        sum=i+u(j,ii);
        if (j==1)
            i=0;
        else
            i=sum;
        end
    end
    q(ii)=i*dy;
end
sum=0;
for j=1:nx
    error=qex(j)-q(j);
    sum=sum+(error*error);
end
an=nx;
q_rms=sqrt(sum/an);

fprintf('nx=%d; ny=%d; re=%d;\n',nx,ny,re)
fprintf('xmax=%d; h=%d; itmax=%d; uin=%d; dt=%.3e\n\n',xmax,h,itmax,uin,dt)
fprintf('RMS for u-component at the outflow is: %.4e\n',u_rms);
fprintf('RMS for the flow rate: %.4e\n\n',q_rms);
%%
% Plots

% Vertical Components
figure(1)
plot(u(:,2),y,'-x',u(:,(nx+1)/2),y,'-x',u(:,nx-1),y,'-x',u(:,nx),y,'-x',uex(:,nx),y,'LineWidth',1.4)
legend('j=2','j=(nx+1)/2','j=nx-1','j=nx','exact solution');
xlabel('u'),ylabel('y');title('Vertical Profiles for u');
figure(2)
plot(v(:,2),y,'--',v(:,(nx+1)/2),y,v(:,nx-1),y,v(:,nx),y,'LineWidth',1.4)
legend('j=2','j=(nx+1)/2','j=nx-1','j=nx');
xlabel('v'),ylabel('y');title('Vertical Profiles for v');

% Along the centerline
figure(3)
plot(x,u(nn,:),xmax,uex(nn,nx),'o','LineWidth',2)
xlabel('x'),ylabel('u');title('Centerline for u');
figure(4)
plot(x,v((ny+1)/2,:),'LineWidth',2)
xlabel('x'),ylabel('v');title('Centerline for v');
figure(5)
plot(x,q,x,qex,'--','LineWidth',2)
xlabel('x'),ylabel('q');title('Centerline for q');

%% To test subroubtine matrix.m, uncomment the code below

% nx=[11 21 31];
% ny=nx;
% for i=1:3
%     rms=artificialpressure(nx(i),ny(i));
%     fprintf('Mesh: %dx%d \nrms: %.4e\n',nx(i),ny(i),rms);
% end

%%
