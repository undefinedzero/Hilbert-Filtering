clear all
clc;clf;
%设置时域
dt=pi/8000;
t=-30*pi:dt:30*pi;
 
%设置采样参数
T=pi/320;
ws=2*pi/T;
tmin=-1000*T;tmax=1000*T;
wc=pi/T;
t1=tmin:T:tmax;

%频率设置
w0=300;w1=2;w2=3;
%幅度设置
p0=1;p1=2;p2=3;
%基带信号
ai=p1*cos(w1*t1);
bi=p2*sin(w2*t1);

%已调信号
x=ai.*p0.*cos(w0.*t1)+bi.*p0.*sin(w0.*t1);
%希尔伯特变换
y=hilbert(x);

%%参数估计%%
a=abs(y);%包络估计
b=unwrap(angle(y));%相位估计
%w0估计
c=diff(b)*320/pi;
sum=0;
for k=1/4*(length(t1)-1):3/4*(length(t1)-1)
    sum=sum+c(k);
end
we0=sum/((length(t1)-1)/2);%时间均值
we1=(max(c(1/4*(length(t1)-1):3/4*(length(t1)-1)))+min(c(1/4*(length(t1)-1):3/4*(length(t1)-1))))/2;%作差法
we2=(max(c(1.5/4*(length(t1)-1):2.5/4*(length(t1)-1)))+min(c(1.5/4*(length(t1)-1):2.5/4*(length(t1)-1))))/2;%作差法
%基带信号的估计
x10=a.*sin(b-we0*t1);
x20=a.*cos(b-we0*t1);
x11=a.*sin(b-we1*t1);
x21=a.*cos(b-we1*t1);
x12=a.*sin(b-we2*t1);
x22=a.*cos(b-we2*t1);
an=0;bn=0;cn=0;x12n=0;x22n=0;

%内插重建
L=length(t);
A=wc*T/pi;
for m=1:L
    y1=A*sinc(wc*(t(m)-t1)/pi);
    an(m)=a*y1';
    bn(m)=b*y1';
    cn(m)=c*y1(1:end-1)';
    x12n(m)=x12*y1';
    x22n(m)=x22*y1';
end

%w0估计
sum=0;
for k=1/4*(length(t)-1):3/4*(length(t)-1)
    sum=sum+cn(k);
end
we01=sum/((length(t)-1)/2);%时间均值
we11=(max(cn(1/4*(length(t)-1):3/4*(length(t)-1)))+min(cn(1/4*(length(t)-1):3/4*(length(t)-1))))/2;%作差法
we21=(max(cn(1.5/4*(length(t)-1):2.5/4*(length(t)-1)))+min(cn(1.5/4*(length(t)-1):2.5/4*(length(t)-1))))/2;%作差法

subplot(3,2,4);
plot(t,an);axis([-10,10,1.2*min(a),1.2*max(a)]);title('包络估计');grid on;
subplot(3,2,2);
plot(t,bn);axis([-10,10,1.2*min(b),1.2*max(b)]);title('相位估计');grid on;
subplot(3,2,3);
plot(t,cn);axis([-10,10,-1.2*w0,1.2*w0]);title('w0估计');grid on;
subplot(3,2,1);
plot(t1,x);axis([-10,10,1.2*min(x),1.2*max(x)]);title('已调信号');grid on;
subplot(3,2,5);
plot(t,x12n);axis([-10,10,1.2*min(x12),1.2*max(x12)]);title('基波信号a(t)估计');grid on;
subplot(3,2,6);
plot(t,x22n);axis([-10,10,1.2*min(x22),1.2*max(x22)]);title('基波信号b(t)估计');grid on;
