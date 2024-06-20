function OK=DGP0
%% Preproceeding
%global paramaters
global Timeadvance;
global Ma_far;
global V_far;
global gamma;
global Ndim;
global Nelement;
global Npoint;
global tol;
% Geometry
global INPOEL;
global COORD;
global BCOND;
global ESUEL;
global INTFAC;
global ESUP1;
global ESUP2;
global Striangle;


%DGP0
Ucurrent=zeros(4,Nelement);%rho rhou rhov rhoe
Vcurrent=zeros(5,Nelement);%rho u v p a 
%-TVDRK3
Uhold=zeros(4,Nelement);
Unhold=zeros(4,Nelement);
afa=[1,3/4,1/3];gamma1=[0,1/4,2/3];belta=[1,1/4,2/3];
Unext=zeros(4,Nelement);


%% Proceeding
%initial condition
%远场初值
for ie=1:Nelement
    Ucurrent(1,ie)=V_far(1,1);
    Ucurrent(2,ie)=V_far(2,1);
    Ucurrent(3,ie)=V_far(3,1);
    Ucurrent(4,ie)=0.5+1/(gamma*(gamma-1)*Ma_far^2);
end
U_0=Ucurrent;
Vcurrent=convert(Ucurrent);

deltat=Deltat(ESUEL,INPOEL,INTFAC,Vcurrent,Ucurrent,Striangle,COORD,BCOND);
%deltat=0.001;
%time iteration
nowtime=0;
endtimes=1e5;
for itimes=1:endtimes
%calculate RHS
    RHS=CRHSP0(Vcurrent,Ucurrent);
%-time advance
    if Timeadvance==1
        for ie=1:Nelement
            Unext(:,ie)=Ucurrent(:,ie)+deltat*RHS(:,ie)/Striangle(1,ie);
        end
    elseif Timeadvance==2
        Uhold=Ucurrent;
    for istage=1:3
        Unhold=afa(istage).*Ucurrent+gamma1(istage).*Uhold+belta(istage)*deltat*(RHS./Striangle);
        Uhold=Unhold;
        if istage~=3
           Vhold=convert(Uhold);
           RHS=CRHSP0(Vhold,Uhold);
        end
    end
    Unext=Uhold;
    end
    
%-break cycle
    JUDGE=zeros(Ndim+2,1);
for ie=1:Nelement
    U_tol=(Unext(:,ie)-Ucurrent(:,ie)).^2;
    JUDGE=JUDGE+(Unext(:,ie)-Ucurrent(:,ie)).^2./U_tol*Striangle(1,ie);
end
    judge=max(sqrt(JUDGE));
    if judge<tol
        break;
    end

%-next level
    Ucurrent=Unext;
    Vcurrent=convert(Ucurrent);
%-next time
    deltat=Deltat(ESUEL,INPOEL,INTFAC,Vcurrent,Ucurrent,Striangle,COORD,BCOND);
    nowtime=nowtime+deltat
itimes

if itimes/1000==round(itimes/1000)

rho_point=zeros(1,Npoint);
rhou_point=zeros(1,Npoint);
rhov_point=zeros(1,Npoint);
rhoe_point=zeros(1,Npoint);
v_point=zeros(2,Npoint);
P_point=zeros(1,Npoint);
a_point=zeros(1,Npoint);
for ip=1:Npoint
    weight=0;
    for ie2=ESUP2(ip)+1:ESUP2(ip+1)
        ie=ESUP1(ie2);
        weight=weight+Striangle(ie);
        rho_point(1,ip)=rho_point(1,ip)+Ucurrent(1,ie)*Striangle(ie);
        rhou_point(1,ip)=rhou_point(1,ip)+Ucurrent(2,ie)*Striangle(ie);
        rhov_point(1,ip)=rhov_point(1,ip)+Ucurrent(3,ie)*Striangle(ie);
        rhoe_point(1,ip)=rhoe_point(1,ip)+Ucurrent(4,ie)*Striangle(ie);

    end
    rho_point(1,ip)=rho_point(1,ip)/weight;
    rhou_point(1,ip)=rhou_point(1,ip)/weight;
    rhov_point(1,ip)=rhov_point(1,ip)/weight;
    rhoe_point(1,ip)=rhoe_point(1,ip)/weight;

end
for ip=1:Npoint
    V=convert([rho_point(1,ip);rhou_point(1,ip);rhov_point(1,ip);rhoe_point(1,ip)]);
    v_point(1:2,ip)=V(2:3,1);
    P_point(1,ip)=V(4,1);
    a_point(1,ip)=V(5,1);    

end

S_far=V_far(4,1)/(V_far(1,1)^gamma);
e=0;
for ie=1:Nelement
    Square=Striangle(1,ie);
    rho_ie=Ucurrent(1,ie);
    rhou_ie=Ucurrent(2,ie);
    rhov_ie=Ucurrent(3,ie);
    p_ie=(gamma-1)*(Ucurrent(4,ie)-0.5*(rhou_ie^2/rho_ie+rhov_ie^2/rho_ie));
    S_ie=p_ie/rho_ie^gamma;
    fi=((S_ie-S_far)/S_far)^2;
    e=e+fi*Square;
end
e=sqrt(e);



%% Postproceeding
%导出网格
% 打开文件以进行写入
DATA_mix=[COORD',v_point',rho_point',P_point',a_point'];
DATA = fopen('DGP0_data.dat', 'w');
% 使用 fprintf 将数据写入文件
fprintf(DATA,'Variables=x,y,u,v,rho,P,a\n');
fprintf(DATA,'Zone n=');fprintf(DATA,num2str(Npoint));fprintf(DATA,',e=');fprintf(DATA,num2str(Nelement));fprintf(DATA,',f=fepoint,et=triangle\n');
fprintf(DATA,'%f\t%f\t%f\t%f\t%f\t%f\t%f\n',DATA_mix');
fprintf(DATA,'%d\t%d\t%d\n',INPOEL);
% 关闭文件
fclose(DATA);
end


end

rho_point=zeros(1,Npoint);
rhou_point=zeros(1,Npoint);
rhov_point=zeros(1,Npoint);
rhoe_point=zeros(1,Npoint);
v_point=zeros(2,Npoint);
P_point=zeros(1,Npoint);
a_point=zeros(1,Npoint);
for ip=1:Npoint
    weight=0;
    for ie2=ESUP2(ip)+1:ESUP2(ip+1)
        ie=ESUP1(ie2);
        weight=weight+Striangle(ie);
        rho_point(1,ip)=rho_point(1,ip)+Ucurrent(1,ie)*Striangle(ie);
        rhou_point(1,ip)=rhou_point(1,ip)+Ucurrent(2,ie)*Striangle(ie);
        rhov_point(1,ip)=rhov_point(1,ip)+Ucurrent(3,ie)*Striangle(ie);
        rhoe_point(1,ip)=rhoe_point(1,ip)+Ucurrent(4,ie)*Striangle(ie);

    end
    rho_point(1,ip)=rho_point(1,ip)/weight;
    rhou_point(1,ip)=rhou_point(1,ip)/weight;
    rhov_point(1,ip)=rhov_point(1,ip)/weight;
    rhoe_point(1,ip)=rhoe_point(1,ip)/weight;

end
for ip=1:Npoint
    V=convert([rho_point(1,ip);rhou_point(1,ip);rhov_point(1,ip);rhoe_point(1,ip)]);
    v_point(1:2,ip)=V(2:3,1);
    P_point(1,ip)=V(4,1);
    a_point(1,ip)=V(5,1);    

end

S_far=V_far(4,1)/(V_far(1,1)^gamma);
e=0;
for ie=1:Nelement
    Square=Striangle(1,ie);
    rho_ie=Ucurrent(1,ie);
    rhou_ie=Ucurrent(2,ie);
    rhov_ie=Ucurrent(3,ie);
    p_ie=(gamma-1)*(Ucurrent(4,ie)-0.5*(rhou_ie^2/rho_ie+rhov_ie^2/rho_ie));
    S_ie=p_ie/rho_ie^gamma;
    fi=((S_ie-S_far)/S_far)^2;
    e=e+fi*Square;
end
e=sqrt(e);



%% Postproceeding
%导出网格
% 打开文件以进行写入
INPOEL=INPOEL';
DATA_mix=[COORD',v_point',rho_point',P_point',a_point'];
DATA = fopen('DGP0_data.dat', 'w');
% 使用 fprintf 将数据写入文件
fprintf(DATA,'Variables=x,y,u,v,rho,P,a\n');
fprintf(DATA,'Zone n=');fprintf(DATA,num2str(Npoint));fprintf(DATA,',e=');fprintf(DATA,num2str(Nelement));fprintf(DATA,',f=fepoint,et=triangle\n');
fprintf(DATA,'%f\t%f\t%f\t%f\t%f\t%f\t%f\n',DATA_mix');
fprintf(DATA,'%d\t%d\t%d\n',INPOEL');
% 关闭文件
fclose(DATA);
OK=[1,itimes,log10(e),log10(min(Striangle))];
end


