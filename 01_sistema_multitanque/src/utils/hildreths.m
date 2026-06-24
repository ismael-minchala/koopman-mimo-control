function DeltaU=hildreths(e,Nc,M,H,gamma,f_ite)
    P=M*inv(2*H)*M';
    d =gamma(:,1)-M*inv(2*H)*2*(f_ite);
    [nx,mx] = size(d);
    x_ini = zeros(nx,mx);
    lamb = x_ini;
    al = 10;
    for km =1:650
        lambda_p = lamb;
        for i=1:nx
            w = P(i,:)*lamb - P(i,i)*lamb(i,1);
                w = w + d(i,1);
                la = -w/P(i,i);
                lamb(i,1) = max(0,la);
            end
            al = (lamb - lambda_p)'*(lamb - lambda_p);
            if al < 10e-36
                break
            else
                eta=zeros(e*Nc);
            end
        end
        eta= inv(2*H)*(2*f_ite)-inv(2*H)*M'*lamb;
        DeltaU(:,1)=eta; 
end