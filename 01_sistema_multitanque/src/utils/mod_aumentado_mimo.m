function [Am,Bm,Cm,Dm]=mod_aumentado_mimo(Ad_cont,Bd_cont,Cd_cont,Dd_cont,n1,q)
Am=[Ad_cont zeros(q,n1)';Cd_cont*Ad_cont eye(q,q)];
Bm=[Bd_cont;Cd_cont*Bd_cont];
Cm=[zeros(q,n1) eye(q,q)];
Dm=Dd_cont;
end