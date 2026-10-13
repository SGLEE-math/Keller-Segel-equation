clear; clc; close all;

nx=100; ny=nx; Lx=-2; Rx=2; Ly=Lx; Ry=Rx;
h=(Rx-Lx)/nx; x=(Lx+0.5*h:h:Rx-0.5*h)'; y=x';
dt=0.1*h;
max_it=100; ns=max_it/100;

Ix=speye(nx); Iy=speye(ny);
LapX=-2*Ix;
for ix=1:nx-1
    LapX(ix,ix+1)=1; LapX(ix+1,ix)=1;
end
LapX(1,1)=LapX(1,1)+LapX(1,2); LapX(nx,nx)=LapX(nx,nx)+LapX(nx,nx-1);
LapX=LapX/h^2;

LapY=-2*Iy;
for iy=1:ny-1
    LapY(iy,iy+1)=1; LapY(iy+1,iy)=1;
end
LapY(1,1)=LapY(1,1)+LapY(1,2); LapY(ny,ny)=LapY(ny,ny)+LapY(ny,ny-1);
LapY=LapY/h^2;

alp=30; a=0.6;
u=1+exp(-alp*((x-a).^2+(y-a).^2))+exp(-alp*((x+a).^2+(y-a).^2)) ...
    +exp(-alp*((x-a).^2+(y+a).^2))+exp(-alp*((x+a).^2+(y+a).^2));
% u=3+2*exp(-20*(sqrt(x.^2+y.^2)-0.7).^2);
v=u;
mass(1)=sum(u,'all')*h^2;
ene(1)=sum(0.5*v.^2-0.5*v.*(LapX*v+(LapY*v')')-u.*v+u.*log(u)-u,'all')*h^2;
M(1)=max(u,[],"all"); m(1)=min(u,[],"all");

surf(x,y,u'); shading interp;
grid on; set(gca,'FontSize',24); box on; axis([Lx Rx Ly Ry 0.9 2.2]);
clim([0.85 1.1]); colormap("jet");
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$u$','Interpreter','latex');
drawnow;
% exportgraphics(gca,'result2d_initial.eps');

Av=(1+dt)*speye(nx*ny)-dt*(kron(Iy,LapX)+kron(LapY,Ix));
dAv=decomposition(Av,'chol');
jx=(1:nx-1)'; jy=(1:ny-1)';
rowX=[jx;jx;jx+1;jx+1];
colX=[jx;jx+1;jx;jx+1];
rowY=[jy;jy;jy+1;jy+1];
colY=[jy;jy+1;jy;jy+1];
theta=0.99;
for it=1:max_it
    maxGS=1000; itGS=0; tol=1e-8; ou=u; tu=0*u; nu=0*u;
    ouxp=0.5*(ou([2:end,end],:)+ou); ouxn=0.5*(ou([1,1:end-1],:)+ou);
    ouyp=0.5*(ou(:,[2:end,end])+ou); ouyn=0.5*(ou(:,[1,1:end-1])+ou);
    
    Kx=cell(ny,1); Ky=cell(nx,1);
    for iy=1:ny
        face=ouxp(1:nx-1,iy)/h^2;
        Kx{iy}=sparse(rowX,colX,[-face;face;face;-face],nx,nx);
    end
    for ix=1:nx
        face=ouyp(ix,1:ny-1).'/h^2;
        Ky{ix}=sparse(rowY,colY,[-face;face;face;-face],ny,ny);
    end
    while true
        itGS=itGS+1;

        uyp=u(:,[2:end,end]); uyn=u(:,[1,1:end-1]);
        vyp=v(:,[2:end,end]); vyn=v(:,[1,1:end-1]);
        for iy=1:ny
            Bx=Kx{iy}*spdiags(1./u(:,iy),0,nx,nx);
            cof=Ix-dt*Bx;
            sor=ou(:,iy)+(dt/h^2) ...
                *(ouxp(:,iy).*(log(u([2:end,end],iy))-v([2:end,end],iy)) ...
                +ouxn(:,iy).*(log(u([1,1:end-1],iy))-v([1,1:end-1],iy)) ...
                +ouyp(:,iy).*(log(uyp(:,iy))-vyp(:,iy)) ...
                +ouyn(:,iy).*(log(uyn(:,iy))-vyn(:,iy)) ...
                -(ouxp(:,iy)+ouxn(:,iy)+ouyp(:,iy)+ouyn(:,iy)) ...
                .*(log(u(:,iy))-v(:,iy)));
            tu(:,iy)=cof\sor;
        end

        for ix=1:nx
            By=Ky{ix}*spdiags(1./u(ix,:).',0,ny,ny);
            cof=Iy-dt*By;
            nu(ix,:) = (cof\tu(ix,:).').';
        end
        nu=ks_positive_step(u,nu,theta);
        err=norm(nu(:)-u(:),'inf');
        u=nu;
        if (err<tol)||(itGS>maxGS)
            % fprintf('it=%d, err=%1.3e\n',itGS,err)
            break;
        end
    end

    sor=v+dt*u;
    v=reshape(dAv\sor(:),nx,ny);

    mass(it+1)=sum(u,'all')*h^2;
    ene(it+1)=sum(0.5*v.^2-0.5*v.*(LapX*v+(LapY*v')')-u.*v+u.*log(u)-u,'all')*h^2;
    M(it+1)=max(u,[],"all"); m(it+1)=min(u,[],"all");

    if mod(it,ns)==0
        clf;
        surf(x,y,u'); shading interp;
        grid on; set(gca,'FontSize',24); box on; axis([Lx Rx Ly Ry 0.9 2.2]);
        clim([0.85 1.1]); colormap("jet");
        xlabel('$x$','Interpreter','latex');
        ylabel('$y$','Interpreter','latex');
        zlabel('$u$','Interpreter','latex');
        drawnow;
        % if it==10
        %     exportgraphics(gca,'result2d_10dt.eps');
        % elseif it==50
        %     exportgraphics(gca,'result2d_50dt.eps');
        % elseif it==100
        %     exportgraphics(gca,'result2d_100dt.eps');
        % end
    end
end

clf;
plot(dt*(0:max_it),M, dt*(0:max_it),m, 'LineWidth',5);
grid on; set(gca,'FontSize',24); box on; axis([0 dt*max_it 0.9 2.2]);
xlabel('$t$','Interpreter','latex');
legend('${\max}(u^n)$','${\min}(u^n)$','interpreter','latex','Location','east');
% exportgraphics(gcf,'max_min_2d.eps');

clf;
plot(dt*(0:max_it),mass, dt*(0:max_it),ene, 'LineWidth',5);
grid on; set(gca,'FontSize',24); box on; axis([0 dt*max_it -25 20]);
xlabel('$t$','Interpreter','latex');
legend('$\mathcal{M}^n_h$','$$\mathcal{E}^n_h$$','interpreter','latex','Location','east');
% exportgraphics(gcf,'mass_energy_2d.eps');
