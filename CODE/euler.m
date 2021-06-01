function [u,v,fu,fv] = euler(u,v,nx,ny,dx,dy,re,dt,uin)
%First-step: use Euler method to initiate Fu and Fv
fu = u_f(u,v,nx,ny,dx,dy,re);
fv = v_f(u,v,nx,ny,dx,dy);
u = u + dt*fu; %update u
v = v + dt*fv; %update v
u(2:ny-1,1)=uin;
end