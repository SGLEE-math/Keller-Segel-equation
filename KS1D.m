clear; clc; close all;

nx=100; Lx=-1; Rx=1; h=(Rx-Lx)/nx; x=(Lx+0.5*h:h:Rx-0.5*h)'; dt=0.1*h;
max_it=10; ns=max_it/10; 

I=speye(nx);
L=-2*speye(nx);
for ix=1:nx-1
    L(ix,ix+1)=1; L(ix+1,ix)=1;
end
L(1,1)=L(1,1)+L(1,2); L(nx,nx)=L(nx,nx)+L(nx,nx-1); L=L/h^2;

u=2*(2.1+cos(pi*x)); v=u;
mass(1)=sum(u,'all')*h;
ene(1)=sum(0.5*v.^2-0.5*v.*(L*v)-u.*v+u.*log(u)-u,'all')*h;

plot(x,u, x,v);
drawnow;

theta=0.99;
Av=(1+dt)*I-dt*L;
dAv=decomposition(Av,'chol');
for it=1:max_it
    maxGS=1000; itGS=0; tol=1e-10; ou=u; 
    oup=0.5*(ou([2:end,end])+ou); oun=0.5*(ou([1,1:end-1])+ou);

    j=(1:nx-1)';
    a=oup(1:nx-1)/h^2;
    K=sparse([j; j+1; j; j+1], [j; j+1; j+1; j], [-a; -a; a; a],nx,nx);
    while true
        itGS=itGS+1;

        B=K*spdiags(1./u,0,nx,nx);
        cof=I-dt*B; sor=ou+dt*K*(log(u)-v);
        nu=cof\sor; nu=ks_positive_step(u,nu,theta);
        
        err=norm(nu-u,'inf');
        u=nu;
        if (err<tol)||(itGS>maxGS)
            % fprintf('it=%d, err=%1.3e\n',itGS,err)
            break;
        end
    end

    sor=v+dt*u;
    v=dAv\sor;

    mass(it+1)=sum(u,'all')*h;
    ene(it+1)=sum(0.5*v.^2-0.5*v.*(L*v)-u.*v+u.*log(u)-u,'all')*h;

    if mod(it,ns)==0
        clf;
        subplot(1,2,1);
        plot(x,u, x,v);
        grid on;

        subplot(1,2,2);
        plot((0:it)/max_it,mass/mass(1));
        hold on;
        plot((0:it)/max_it,ene/ene(1));
        grid on;
        drawnow;
    end
end
