function u=ks_positive_step(u,candidate,theta)
d=candidate-u;
alpha=1; negative = d<0;
if any(negative(:))
    alpha=min(1,theta*min(-u(negative)./d(negative)));
end
u=u+alpha*d;
end
