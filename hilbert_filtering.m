clear all
clc;clf;
%����ʱ��
dt=pi/8000;
t=-30*pi:dt:30*pi;
 
%���ò�������
T=pi/320;
ws=2*pi/T;
tmin=-1000*T;tmax=1000*T;
wc=pi/T;
t1=tmin:T:tmax;

%Ƶ������
w0=300;w1=2;w2=3;
%��������
p0=1;p1=2;p2=3;
%�����ź�
ai=p1*cos(w1*t1);
bi=p2*sin(w2*t1);

%�ѵ��ź�
x=ai.*p0.*cos(w0.*t1)+bi.*p0.*sin(w0.*t1);
%ϣ�����ر任
y=hilbert(x);

%%��������%%
a=abs(y);%�������
b=unwrap(angle(y));%��λ����
%w0����
c=diff(b)*320/pi;
sum=0;
for k=1/4*(length(t1)-1):3/4*(length(t1)-1)
    sum=sum+c(k);
end
we0=sum/((length(t1)-1)/2);%ʱ���ֵ
we1=(max(c(1/4*(length(t1)-1):3/4*(length(t1)-1)))+min(c(1/4*(length(t1)-1):3/4*(length(t1)-1))))/2;%���
we2=(max(c(1.5/4*(length(t1)-1):2.5/4*(length(t1)-1)))+min(c(1.5/4*(length(t1)-1):2.5/4*(length(t1)-1))))/2;%���
%�����źŵĹ���
x10=a.*sin(b-we0*t1);
x20=a.*cos(b-we0*t1);
x11=a.*sin(b-we1*t1);
x21=a.*cos(b-we1*t1);
x12=a.*sin(b-we2*t1);
x22=a.*cos(b-we2*t1);
an=0;bn=0;cn=0;x12n=0;x22n=0;

%�ڲ��ؽ�
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

%w0����
sum=0;
for k=1/4*(length(t)-1):3/4*(length(t)-1)
    sum=sum+cn(k);
end
we01=sum/((length(t)-1)/2);%ʱ���ֵ
we11=(max(cn(1/4*(length(t)-1):3/4*(length(t)-1)))+min(cn(1/4*(length(t)-1):3/4*(length(t)-1))))/2;%���
we21=(max(cn(1.5/4*(length(t)-1):2.5/4*(length(t)-1)))+min(cn(1.5/4*(length(t)-1):2.5/4*(length(t)-1))))/2;%���

subplot(3,2,4);
plot(t,an);axis([-10,10,1.2*min(a),1.2*max(a)]);title('�������');grid on;
subplot(3,2,2);
plot(t,bn);axis([-10,10,1.2*min(b),1.2*max(b)]);title('��λ����');grid on;
subplot(3,2,3);
plot(t,cn);axis([-10,10,-1.2*w0,1.2*w0]);title('w0����');grid on;
subplot(3,2,1);
plot(t1,x);axis([-10,10,1.2*min(x),1.2*max(x)]);title('�ѵ��ź�');grid on;
subplot(3,2,5);
plot(t,x12n);axis([-10,10,1.2*min(x12),1.2*max(x12)]);title('�����ź�a(t)����');grid on;
subplot(3,2,6);
plot(t,x22n);axis([-10,10,1.2*min(x22),1.2*max(x22)]);title('�����ź�b(t)����');grid on;
