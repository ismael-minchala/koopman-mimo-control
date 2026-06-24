function [u]=calc_accion_control(e,Nc,H,Phi,BarRs,F,Xf,u,M,gamma)
    f_ite=Phi'*(BarRs-F*Xf);
    eta=inv(H)*f_ite;
    kk = 0;
