function [F,Phi]=mat_f_phi(Am,Bm,Cm,Np,Nc,q,e)
h(:,:)=Cm;
F(:,:)=Cm*Am;
for i=1:Np-1
h(q*i+1:q*i+q,:)=h(q*i+1-q:q*i,:)*Am;
F(q*i+1:q*i+q,:)= F(q*i+1-q:q*i,:)*Am;
end
v=h*Bm;
Phi=zeros(q*Np,e*Nc);
Phi(:,1:e)=v; 
for i=1:Nc-1
    Phi(:,e*i+1:e*i+e)=[zeros(i*q,e);v(1:q*(Np-i),:)]; 
end
