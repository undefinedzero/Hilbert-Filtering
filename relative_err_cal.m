fs=1000;
T=1/fs;
t1=-10;t2=10;
t=t1:T:t2;

%Ƶ������
w0=300;w1=2;w2=3;
%��������
p0=1;p1=2;p2=3;
%�����ź�
ai=p1*cos(w1*t);
bi=p2*sin(w2*t);

%�ѵ��ź�
x=ai.*p0.*cos(w0.*t)+bi.*p0.*sin(w0.*t);
fai=300*t+unwrap(asin(ai./sqrt(ai.^2+bi.^2)))-w0*t1;
%ϣ�����ر任
y=hilbert(x);

%%��������%%
a=abs(y);%�������
b=unwrap(angle(y));%��λ����
%w0����
c=diff(b)*1000;
sum=0;
for k=1/4*(length(t)-1):3/4*(length(t)-1)
    sum=sum+c(k);
end
we=sum/(1/2*(length(t)-1));
%�����źŵĹ���
x1=a.*sin(b-we*t);
x2=a.*cos(b-we*t);

v1=max(x1)-min(x1);
v2=max(x2)-min(x2);

%��ͼ
subplot(2,3,1);
plot(t,a);axis([-10,10,1.2*min(a),1.2*max(a)]);
title('�������');grid on;
subplot(2,3,2);
plot(t,b);axis([-10,10,1.2*min(b),1.2*max(b)]);
title('��λ����');grid on;
subplot(2,3,3);
plot(t(1:20000),c);axis([-10,10,-1.2*w0,1.2*w0]);
title('w0����');grid on;
subplot(2,3,4);
plot(t,(a-sqrt(ai.^2+bi.^2))./a);axis([-10,10,-0.1,0.1]);
title('����������');grid on;
subplot(2,3,5);
plot(t,(b-fai)./b);axis([-10,10,-0.1,0.1]);
title('��λ������');grid on;
subplot(2,3,6);
plot(t(1:20000),(c-w0)./c);axis([-10,10,-0.1,0.1]);
title('w0������');grid on;
