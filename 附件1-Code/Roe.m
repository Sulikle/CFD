function F_ij=Roe(v_n_l,v_n_r,Vl,Ul,Vr,Ur,n_vector)
global gamma;

%焓
H_i=0.5*(Vl(2,1)^2+Vl(3,1)^2)+Vl(5,1)^2/(gamma-1);
H_j=0.5*(Vr(2,1)^2+Vr(3,1)^2)+Vr(5,1)^2/(gamma-1);
%calculate Roe average state
R_ij=sqrt(Vr(1,1)/Vl(1,1));
rho_ij=R_ij*Vl(1,1);
U_ij=(R_ij*Vr(2,1)+Vl(2,1))/(R_ij+1);
V_ij=(R_ij*Vr(3,1)+Vl(3,1))/(R_ij+1);
H_ij=(R_ij*H_j+H_i)/(R_ij+1);
C_ij=sqrt((gamma-1)*(H_ij-(U_ij^2+V_ij^2)/2));

%calculate lambda
lambda(1,1)=dot([U_ij,V_ij],n_vector);
lambda(1,2)=lambda(1,1);
lambda(1,3)=lambda(1,1)+C_ij;
lambda(1,4)=lambda(1,1)-C_ij;
%modify
epsilon=zeros(1,4);
epsilon(1,1)=max([0,lambda(1,1)-v_n_l,v_n_r-lambda(1,1)]);
epsilon(1,2)=max([0,lambda(1,2)-v_n_l,v_n_r-lambda(1,2)]);
epsilon(1,3)=max([0,lambda(1,3)-v_n_l,v_n_r-lambda(1,3)]);
epsilon(1,4)=max([0,lambda(1,4)-v_n_l,v_n_r-lambda(1,4)]);

abs_lambda=zeros(1,4);
for i=1:4
    if abs(lambda(1,i))>=epsilon(1,i)
        abs_lambda(1,i)=abs(lambda(1,i));
    else
        abs_lambda(1,i)=0.5*(lambda(1,i)^2/epsilon(1,i)+epsilon(1,i));
    end
end
r1=[1,U_ij,V_ij,0.5*(U_ij^2+V_ij^2)];
r2=[0,C_ij*n_vector(1,2),-C_ij*n_vector(1,1),C_ij*(U_ij*n_vector(1,2)-V_ij*n_vector(1,1))];
r3=(rho_ij/(2*C_ij))*[1,U_ij+C_ij*n_vector(1,1),V_ij+C_ij*n_vector(1,2),H_ij+C_ij*lambda(1,1)];
r4=(rho_ij/(2*C_ij))*[1,U_ij-C_ij*n_vector(1,1),V_ij-C_ij*n_vector(1,2),H_ij-C_ij*lambda(1,1)];
R=[r1',r2',r3',r4'];
A=diag(abs_lambda);
deltaW=[Vr(1,1)-Vl(1,1)-(Vr(4,1)-Vl(4,1))/C_ij^2;(Vr(2,1)-Vl(2,1))*n_vector(1,2)-(Vr(3,1)-Vl(3,1))*n_vector(1,1);(v_n_r-v_n_l)+(Vr(4,1)-Vl(4,1))/(rho_ij*C_ij);-(v_n_r-v_n_l)+(Vr(4,1)-Vl(4,1))/(rho_ij*C_ij)];
f_i=[Vl(1,1)*v_n_l;v_n_l*Vl(1,1)*Vl(2,1)+Vl(4,1)*n_vector(1,1);v_n_l*Vl(1,1)*Vl(3,1)+Vl(4,1)*n_vector(1,2);v_n_l*(Ul(4,1)+Vl(4,1))];
%f_i=[Vl(1,1)*Vl(2,1);Vl(1,1)*Vl(2,1)^2+Vl(4,1);Vl(1,1)*Vl(2,1)*Vl(3,1);Vl(2,1)*(Ul(4,1)+Vl(4,1))]*n_vector(1,1)+[Vl(1,1)*Vl(3,1);Vl(1,1)*Vl(2,1)*Vl(3,1);Vl(1,1)*Vl(3,1)^2+Vl(4,1);Vl(3,1)*(Ul(4,1)+Vl(4,1))]*n_vector(1,2);
f_j=[Vr(1,1)*v_n_r;v_n_r*Vr(1,1)*Vr(2,1)+Vr(4,1)*n_vector(1,1);v_n_r*Vr(1,1)*Vr(3,1)+Vr(4,1)*n_vector(1,2);v_n_r*(Ur(4,1)+Vr(4,1))];
%f_j=[Vr(1,1)*Vr(2,1);Vr(1,1)*Vr(2,1)^2+Vr(4,1);Vr(1,1)*Vr(2,1)*Vr(3,1);Vr(2,1)*(Ur(4,1)+Vr(4,1))]*n_vector(1,1)+[Vr(1,1)*Vr(3,1);Vr(1,1)*Vr(2,1)*Vr(3,1);Vr(1,1)*Vr(3,1)^2+Vr(4,1);Vr(3,1)*(Ur(4,1)+Vr(4,1))]*n_vector(1,2);
F_ij=0.5*(f_i+f_j)-0.5*R*A*deltaW;

end