function test=test(ok)
tf=[-sqrt(15)/5,0,sqrt(15)/5];Wf=[5/9,8/9,5/9];%-Gauss in face
x1=0;x2=1;y1=0;y2=1;
tau=sqrt(2);
cos_theta=(x2-x1)/tau;
sin_theta=(y2-y1)/tau;
test=0;
if ok==1
    for i=1:3
        ri=0.5*tau*tf(i)+0.5*tau;
        fi=x1+ri*cos_theta+y1+ri*sin_theta;
        test=test+Wf(i)*fi;
    end
    test=test*0.5*tau;

end
end