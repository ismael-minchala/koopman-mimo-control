function [u]=calc_accion_control(e,Nc,H,Phi,BarRs,F,Xf,u,M,gamma)
    f_ite=Phi'*(BarRs-F*Xf);
    eta=inv(H)*f_ite;
    kk = 0;
%Se comprueba si la solución no viola restricciones
    for i=1:size(M,1)
        if M(i,:)*eta > gamma(i,1)     
            kk = kk + 1;         
        end
    end
     if kk==0
       DeltaU=eta;                   %Solución analítica   
     else
        DeltaU=hildreths(e,Nc,M,H,gamma,f_ite);
     end
         
    deltau = DeltaU(1:e,1);
    u = u+deltau;
    
end