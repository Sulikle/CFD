function Vcurrent=convert(Ucurrent)
%Vcurrent store rho u v p a
global gamma;
Vcurrent(1,:)=Ucurrent(1,:);
Vcurrent(2,:)=Ucurrent(2,:)./Ucurrent(1,:);
Vcurrent(3,:)=Ucurrent(3,:)./Ucurrent(1,:);
Vcurrent(4,:)=(gamma-1)*(Ucurrent(4,:)-0.5*Vcurrent(1,:).*(Vcurrent(2,:).^2+Vcurrent(3,:).^2));
Vcurrent(5,:)=sqrt(gamma*Vcurrent(4,:)./Vcurrent(1,:));
end