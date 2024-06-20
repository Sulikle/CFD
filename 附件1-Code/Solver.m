clc
clear
close all
%% Receive information
global Ma_far;
global afa_attack;
global gamma;
global CFL;
global tol;
global V_far;
global Timeadvance;
global Flux;
global rDG;
global Mesh;
% Geometry
global INPOEL;
global COORD;
global BCOND;
global ESUEL;
global INTFAC;
global ESUP1;
global ESUP2;
global Striangle;

%read solver
file = fopen('Solver.dat', 'rb');
for i = 1:2
    line=fgetl(file);
end
Mode_section = fscanf(file, '%lf',[1,5]);
Mesh=Mode_section(1,1);
Timeadvance=Mode_section(1,2);
Flux=Mode_section(1,3);
Pn=Mode_section(1,4);Pm=Mode_section(1,5);

for i = 1:4
    line=fgetl(file);
end
Mode_section = fscanf(file, '%lf',[1,3]);
tol=Mode_section(1,1);CFL=Mode_section(1,2);rDG=Mode_section(1,3);

for i = 1:4
    line=fgetl(file);
end
Mode_section = fscanf(file, '%lf',[1,3]);
Ma_far=Mode_section(1,1);afa_attack=Mode_section(1,2)*pi/180;gamma=Mode_section(1,3);

V_far=[1;cos(afa_attack);sin(afa_attack);1/(gamma*Ma_far^2);sqrt(gamma*(1/(gamma*Ma_far^2)))];
fclose(file);

%% Input data
addpath(genpath('/Users/mac/Documents/CFD-Exercise-main/研一下/Project2/pro2_mesh'))
addpath(genpath('/Users/mac/Documents/CFD-Exercise-main/研一下/Project2/DGP0'))
addpath(genpath('/Users/mac/Documents/CFD-Exercise-main/研一下/Project2/DGP1'))
addpath(genpath('/Users/mac/Documents/CFD-Exercise-main/研一下/Project2/rDGP0P1'))
global Ndim;
global Nelement;
global Npoint;
global Nnode;
global NBface;
global Nface;
%Input mesh
if Mesh==1.1
    feflo_bump = fopen('pro2_mesh/feflo.domn.cylinder.coarse', 'rb');%open the file
elseif Mesh==1.2    
    feflo_bump = fopen('pro2_mesh/feflo.domn.cylinder.medium', 'rb');%open the file
elseif Mesh==1.3
    feflo_bump = fopen('pro2_mesh/feflo.domn.cylinder.fine', 'rb');%open the file
elseif Mesh==1.4
    feflo_bump = fopen('pro2_mesh/feflo.domn.cylinder.vfine', 'rb');%open the file
elseif Mesh==2
    feflo_bump = fopen('pro2_mesh/feflo.domn.naca0012', 'rb');%open the file
elseif Mesh==3
    feflo_bump = fopen('pro2_mesh/feflo.domn.bump', 'rb');%open the file
end

skiprows=fscanf(feflo_bump, '%lf',1);% skip some rows
for i = 1:skiprows+2
    line=fgetl(feflo_bump);
end
Ndata = fscanf(feflo_bump, '%lf',[2,1]);
%read dimention&Nnode
Ndim=Ndata(1,1);Nnode=Ndata(2,1);
for i = 1:2
    line=fgetl(feflo_bump);
end
Ndata = fscanf(feflo_bump, '%lf',[3,1]);
Nelement=Ndata(1,1);Npoint=Ndata(2,1);NBface=Ndata(3,1);%Information of elements&points&boundary

for i = 1:3
    line=fgetl(feflo_bump);
end
%read data per row
C = textscan(line, '%f');
numbers = length(str2num(line));
INPOEL = fscanf(feflo_bump, '%lf',[numbers,Nelement]);
INPOEL=[C{1,1},INPOEL];
INPOEL=INPOEL(2:4,:);%每一列存储的是这个单元的第1，2，3个点

%-store the coordinates of the points————COORD
for i = 1:1% skip some rows
    line=fgetl(feflo_bump);
end
COORD = fscanf(feflo_bump, '%lf',[Ndim+1,Npoint]);
COORD=COORD(2:3,:);%每一列存储的是这个点的横纵坐标

%-store Boundary conditions
for i = 1:2+Npoint+2
    line=fgetl(feflo_bump);
end
%由于部分数据缺失，因此BOCND采取这种读取方式
BCOND=zeros(3,NBface);
for i=1:NBface
    C = textscan(line, '%f');
    BCOND(:,i)=C{1,1}(2:4,1);
    line=fgetl(feflo_bump);
end

% C = textscan(line, '%f');
% numbers = length(str2num(line));
% BCOND=fscanf(feflo_bump, '%lf',[numbers,Nface]);
% BCOND=[C{1,1},BCOND];
% BCOND=BCOND(2:4,:);%每一列存储的是该边界处的左右点以及flag，2代表固定面，4代表无穷远处
fclose(feflo_bump);
 
%-Elements surrounding points
ESUP2=zeros(1,Npoint+1);
for ie=1:Nelement
    for in=1:Nnode
        ESUP2(INPOEL(in,ie)+1)=ESUP2(INPOEL(in,ie)+1)+1;%斜对面存储该point涉及到的单元数
    end
end

for IPOIN=2:Npoint+1
    ESUP2(IPOIN)=ESUP2(IPOIN)+ESUP2(IPOIN-1);%将这些单元标号
end
MESUP=ESUP2(Npoint+1);
ESUP1=zeros(1,MESUP);
for ie=1:Nelement
    for in=1:Nnode
        IPOIN=INPOEL(in,ie);%点的标号
        ISTOR=ESUP2(IPOIN)+1;%改点所用单元的初地址
        ESUP2(IPOIN)=ISTOR;%更新地址
        ESUP1(ISTOR)=ie;%记录单元
    end
end
for IPOIN=Npoint+1:-1:2
    ESUP2(IPOIN)=ESUP2(IPOIN-1);
end
ESUP2(1)=0;

%Elements surrounding elements
ESUEL=zeros(3,Nelement);
LPOIN=zeros(1,Npoint);
%-Auxiliary array
LNOFA=zeros(2,3);
LNOFA(1,1)=2;LNOFA(2,1)=3;
LNOFA(1,2)=3;LNOFA(2,2)=1;
LNOFA(1,3)=1;LNOFA(2,3)=2;

for ie=1:Nelement
    for iface=1:Nnode
        LHELP=INPOEL(LNOFA(:,iface),ie);
        LPOIN(1,LHELP(1,1))=1;LPOIN(1,LHELP(2,1))=1;
        IPOIN=INPOEL(LNOFA(1,iface),ie);%select a point of the face
        FLAG=0;
        for ISTOR=ESUP2(IPOIN)+1:ESUP2(IPOIN+1)
            JELEM=ESUP1(ISTOR);
            if(JELEM~=ie)
                for JFAEL=1:Nnode%loop over the faces of element JELEM
                    ICOUN=0;
                    for JOFA=1:Nnode-1
                        JPOIN=INPOEL(LNOFA(JOFA,JFAEL),JELEM);
                        if(LPOIN(1,JPOIN)==1)
                            ICOUN=ICOUN+1;
                        end
                    end
                        if(ICOUN==Nnode-1)
                            ESUEL(iface,ie)=JELEM;
                            ESUEL(JFAEL,JELEM)=ie;
                            FLAG=1;
                        end
                 end
            end
        end

        if(FLAG==0)
            for IFACE=1:NBface
                if((BCOND(1,IFACE)==LHELP(1,1)&&BCOND(2,IFACE)==LHELP(2,1))||(BCOND(2,IFACE)==LHELP(1,1)&&BCOND(1,IFACE)==LHELP(2,1)))
                    ESUEL(iface,ie)=IFACE+Nelement;
                    break;
                end
            end
        end
        LPOIN=zeros(1,Npoint);
    end
end


Nface=NBface;%Numbers of all faces
for ie=1:Nelement
    for IFAEL=1:Nnode
        JELEM=ESUEL(IFAEL,ie);
        if(JELEM>ie&&JELEM<=Nelement)
            Nface=Nface+1;
        end
    end
end

NIFFA=7;%number of information of face
INTFAC=zeros(NIFFA,Nface);%interface array 
iface=1;
for ie=1:Nelement
    for IFAEL=1:Nnode
        JELEM=ESUEL(IFAEL,ie);
        if(JELEM>ie)
            INTFAC(2,iface)=JELEM;%store right element
            INTFAC(1,iface)=ie;%store left element
            INTFAC(3,iface)=INPOEL(LNOFA(1,IFAEL),ie);%store left point of face
            INTFAC(4,iface)=INPOEL(LNOFA(2,IFAEL),ie);%store right point of face
            if(JELEM>Nelement)
                INTFAC(5,iface)=BCOND(3,JELEM-Nelement);%judge whether this face is fixed or far or internal
            else
                INTFAC(5,iface)=-1;
            end
            iel=ie;ier=JELEM;
            %calculate V_n
            ip1=INPOEL(LNOFA(1,IFAEL),ie);%points in face
            ip2=INPOEL(LNOFA(2,IFAEL),ie);
             for IPL=1:Nnode
                 if(INPOEL(IPL,iel)~=ip1&&INPOEL(IPL,iel)~=ip2)
                     ip3=INPOEL(IPL,iel);%another point in iel
                 end
             end
             Judge_vector1=[COORD(1,ip2)-COORD(1,ip1),COORD(2,ip2)-COORD(2,ip1)];
             Judge_vector2=[COORD(1,ip3)-COORD(1,ip1),COORD(2,ip3)-COORD(2,ip1)];%该单元内部点与边界点构成的向量
             V_a=COORD(1,ip2)-COORD(1,ip1);V_b=COORD(2,ip2)-COORD(2,ip1);
             if V_b==0
                  V_product1=[0,1];V_product2=[0,-1];
             else
                  V_product1=[-V_b,V_a]/norm([-V_b,V_a]);%单位法线向量
                  V_product2=[V_b,-V_a]/norm([V_b,-V_a]);
             end

             if dot(Judge_vector2,V_product1)<0
                 V_n=V_product1;
             elseif dot(Judge_vector2,V_product2)<0
                 V_n=V_product2;
             end
             INTFAC(6,iface)=V_n(1,1);
             INTFAC(7,iface)=V_n(1,2);
             tau_ij=sqrt((COORD(1,ip2)-COORD(1,ip1))^2+(COORD(2,ip2)-COORD(2,ip1))^2);
             INTFAC(8,iface)=tau_ij;

            iface=iface+1;
        end
    end
end

%calculate squares of elements
Striangle=zeros(1,Nelement);
for ie=1:Nelement
    ipoin1=INPOEL(1,ie);ipoin2=INPOEL(2,ie);ipoin3=INPOEL(3,ie);
    HPcycle=[ipoin1,ipoin2,ipoin3,ipoin1,ipoin2];c_ie=zeros(1,3);
    for k=1:Nnode%先计算c
           c_ie(k)=COORD(1,HPcycle(k+1))*COORD(2,HPcycle(k+2))-COORD(1,HPcycle(k+2))*COORD(2,HPcycle(k+1));
    end
    Striangle(1,ie)=0.5*sum(c_ie);
end


%% Method

if Pn==0&&Pm==0
    OK=DGP0;
elseif Pn==1&&Pm==1
    OK=DGP1;
elseif Pn==0&&Pm==1
    OK=rDGP0P1;
end

