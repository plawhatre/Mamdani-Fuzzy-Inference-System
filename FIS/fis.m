clc,clear,close all

%input
cx=input('Enter the distance');
cr=input('Enter road condition');

%Distance
l_x=zeros(1,1001);
m_x=zeros(1,1001);
h_x=zeros(1,1001);
vh_x=zeros(1,1001);
x=0:0.1:100;
c=0;
for i=0:0.1:100
    c=c+1;
    l_x(c)=trap(i,0,0,8,33);
    m_x(c)=trap(i,25,33,42,50);
    h_x(c)=trap(i,42,50,58,66);
    vh_x(c)=trap(i,58,66,100,100);
    if i==cx
        if l_x(c)~=0
            x_c1=1;
            x_v1=l_x(c);
        else
            x_c1=nan;
            x_v1=nan;
        end
        if m_x(c)~=0
            x_c2=2;
            x_v2=m_x(c);
        else
            x_c2=nan;
            x_v2=nan;
        end
        if h_x(c)~=0
            x_c3=3;
            x_v3=h_x(c);
        else
            x_c3=nan;
            x_v3=nan;
        end
        if vh_x(c)~=0
            x_c4=4;
            x_v4=vh_x(c);
        else
            x_c4=nan;
            x_v4=nan;
        end
        fx=[x_c1,x_c2,x_c3,x_c4;x_v1,x_v2,x_v3,x_v4];
    end
end
figure
subplot(3,1,1)
plot(x,l_x);
axis([0 100 0 1.1]);
hold on
plot(x,m_x);
hold on
plot(x,h_x);
hold on
plot(x,vh_x);
title('input1--->distance')
legend('low','medium','large','very large')

%Road condition
l_r=zeros(1,1001);
m_r=zeros(1,1001);
h_r=zeros(1,1001);
r=0:0.001:1;
c=0;
for i=0:0.001:1
    c=c+1;
    l_r(c)=trap(i,0,0,0.1,0.3);
    m_r(c)=trap(i,0.1,0.3,0.5,0.6);
    h_r(c)=trap(i,0.5,0.6,1,1);
    if i==cr
        if l_r(c)~=0
            r_c1=1;
            r_v1=l_r(c);
        else
            r_c1=nan;
            r_v1=nan;
        end
        if m_r(c)~=0
            r_c2=2;
            r_v2=m_r(c);
        else
            r_c2=nan;
            r_v2=nan;
        end
        if h_r(c)~=0
            r_c3=3;
            r_v3=h_r(c);
        else
            r_c3=nan;
            r_v3=nan;
        end
        fr=[r_c1,r_c2,r_c3,;r_v1,r_v2,r_v3];
    end
end
subplot(3,1,2)
plot(r,l_r);
axis([0 1 0 1.1]);
hold on
plot(r,m_r);
hold on
plot(r,h_r);
title('input2--->road condition')
legend('low','medium','high')

%Speed
l_v=zeros(1,1201);
m_v=zeros(1,1201);
h_v=zeros(1,1201);
vh_v=zeros(1,1201);
v=0:0.1:120;
c=0;
for i=0:0.1:120
    c=c+1;
    l_v(c)=trap(i,0,0,10,40);
    m_v(c)=trap(i,30,40,50,60);
    h_v(c)=trap(i,50,60,70,80);
    vh_v(c)=trap(i,70,80,120,120);
end
subplot(3,1,3)
plot(v,l_v);
axis([0 120 0 1.1]);
hold on
plot(v,m_v);
hold on
plot(v,h_v);
hold on
plot(v,vh_v);
title('output--->speed')
legend('low','medium','high','very high')


%class identification
r1=[];
r2=[];
r3=[];
rx=[];
rr=[];
for i=1:4
    if isnan(fx(1,i))==0
        r1=[r1 i];
        rx=[rx;i fx(2,i)];
    end
end
for i=1:3
    if isnan(fr(1,i))==0
        r2=[r2 i];
        rr=[rr;i fr(2,i)];
    end
end
rt=[];
for i=1:length(r1)
    for j=1:length(r2)
        rt=[rt;r1(i),r2(j)];
        tmp=rule(r1(i),r2(j));
        r3=[r3 tmp];
    end
end

%pre-Centroid
p1=rx(:,2);
p2=rr(:,2);
p3=rx(:,1);
p4=rr(:,1);
p=[];
for i=1:length(p1)
    for j=1:length(p2)
        p=[p;p1(i) p2(j)];
    end
end
p(:,3)=min(p(:,1),p(:,2));
p(:,4)=r3;
pz=[];
for i=1:length(p3)
    for j=1:length(p4)
        pz=[pz;p3(i) p4(j)];
    end
end
p=[p,pz];

%in1
fl_x=[];
fm_x=[];
fh_x=[];
fvh_x=[];
for i=1:length(p(:,5))
    if p(i,5)==1
        c=0;
        for j=0:0.1:100
            c=c+1;
            fl_x(1,c,i)=trap(j,0,0,8,33);
            fl_x(2,c,i)=j;
        end
    end
    if p(i,5)==2
        c=0;
        for j=0:0.1:100
            c=c+1;
            fm_x(1,c,i)=trap(j,25,33,42,50);
            fm_x(2,c,i)=j;
        end
    end
    if p(i,5)==3
        c=0;
        for j=0:0.1:100
            c=c+1;
            fh_x(1,c,i)=trap(j,42,50,58,66);
            fh_x(2,c,i)=j;
        end
    end
    if p(i,5)==4
        c=0;
        for j=0:0.1:100
            c=c+1;
            fvh_x(1,c,i)=trap(j,58,66,100,100);
            fvh_x(2,c,i)=j;
        end
    end
end

%in2
fl_r=[];
fm_r=[];
fh_r=[];
for i=1:length(p(:,6))
    if p(i,6)==1
        c=0;
        for j=0:0.001:1
            c=c+1;
            fl_r(1,c,i)=trap(j,0,0,0.1,0.3);
            fl_r(2,c,i)=j;
        end
    end
    if p(i,6)==2
        c=0;
        for j=0:0.001:1
            c=c+1;
            fm_r(1,c,i)=trap(i,0.1,0.3,0.5,0.6);
            fm_r(2,c,i)=j;
        end
    end
    if p(i,6)==3
        c=0;
        for j=0:0.001:1
            c=c+1;
            fh_r(1,c,i)=trap(i,0.5,0.6,1,1);
            fh_r(2,c,i)=j;
        end
    end    
end


%out
fl_v=[];
fm_v=[];
fh_v=[];
fvh_v=[];
for i=1:length(p(:,4))
    if p(i,4)==1
        c=0;
        for j=0:0.1:120
            c=c+1;
            fl_v(1,c,i)=trap(j,0,0,10,40);
            fl_v(2,c,i)=j;
        end
    end
    if p(i,4)==2
        c=0;
        for j=0:0.1:120
            c=c+1;
            fm_v(1,c,i)=trap(j,30,40,50,60);
            fm_v(2,c,i)=j;
        end
    end
    if p(i,4)==3
        c=0;
        for j=0:0.1:120
            c=c+1;
            fh_v(1,c,i)=trap(j,50,60,70,80);
            fh_v(2,c,i)=j;
        end
    end
    if p(i,4)==4
        c=0;
        for j=0:0.1:120
            c=c+1;
            fvh_v(1,c,i)=trap(j,70,80,120,120);
            fvh_v(2,c,i)=j;
        end
    end
end
p
ayip1=[];
cp=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %1in1
% if isempty(fl_x)==0
%     i=1;
%     while i<=length(fl_x(1,1,:))
%         if sum(fl_x(1,:,i))==0
%             fl_x(:,:,i)=[];
%         else
%             i=i+1;
%         end
%     end
%     c=0;
%     for j=1:length(p(:,5))
%         if p(j,5)==1
%             c=c+1;
%             x=fl_x(2,:,c);
%             y1=fl_x(1,:,c);
%             y2=p(j,1)*ones(1,length(y1));
%             for i=1:length(y1)
%                 yf(i)=min(y1(i),y2(i));
%             end
%             cp=cp+1;
%             xtip1(cp,:)=x;
%             yip1(cp,:)=yf;
%             ytip1(cp,:)=y1;
%             ayip1=[ayip1;yf];
%         end
%     end
%     axip1=x;
% end
% %2in1
% if isempty(fm_x)==0
%     i=1;
%     while i<=length(fm_x(1,1,:))
%         if sum(fm_x(1,:,i))==0
%             fm_x(:,:,i)=[];
%         else
%             i=i+1;
%         end
%     end
%     c=0;
%     for j=1:length(p(:,5))
%         if p(j,5)==2
%             c=c+1;
%             x=fm_x(2,:,c);
%             y1=fm_x(1,:,c);
%             y2=p(j,1)*ones(1,length(y1));
%             for i=1:length(y1)
%                 yf(i)=min(y1(i),y2(i));
%             end
%             cp=cp+1;
%             xtip1(cp,:)=x;
%             yip1(cp,:)=yf;
%             ytip1(cp,:)=y1;
%             ayip1=[ayip1;yf];
%         end
%     end
%     axip1=x;
% end
% %3in1
% if isempty(fh_x)==0
%     i=1;
%     while i<=length(fh_x(1,1,:))
%         if sum(fh_x(1,:,i))==0
%             fh_x(:,:,i)=[];
%         else
%             i=i+1;
%         end
%     end
%     c=0;
%     for j=1:length(p(:,5))
%         if p(j,5)==3
%             c=c+1;
%             x=fh_x(2,:,c);
%             y1=fh_x(1,:,c);
%             y2=p(j,1)*ones(1,length(y1));
%             for i=1:length(y1)
%                 yf(i)=min(y1(i),y2(i));
%             end
%             cp=cp+1;
%             xtip1(cp,:)=x;
%             yip1(cp,:)=yf;
%             ytip1(cp,:)=y1;
%             ayip1=[ayip1;yf];
%         end
%     end
%     axip1=x;
% end
% %4in1
% if isempty(fvh_x)==0
%     i=1;
%     while i<=length(fvh_x(1,1,:))
%         if sum(fvh_x(1,:,i))==0
%             fvh_x(:,:,i)=[];
%         else
%             i=i+1;
%         end
%     end
%     c=0;
%     for j=1:length(p(:,5))
%         if p(j,5)==4
%             c=c+1;
%             x=fvh_x(2,:,c);
%             y1=fvh_x(1,:,c);
%             y2=p(j,1)*ones(1,length(y1));
%             for i=1:length(y1)
%                 yf(i)=min(y1(i),y2(i));
%             end
%             cp=cp+1;
%             xtip1(cp,:)=x;
%             yip1(cp,:)=yf;
%             ytip1(cp,:)=y1;
%             ayip1=[ayip1;yf];
%         end
%     end
%     axip1=x;
% end
% 
% 
% 
% 
% 
% 
% ayip2=[];
% cp=0;
% %1in2
% if isempty(fl_r)==0
%     i=1;
%     while i<=length(fl_r(1,1,:))
%         if sum(fl_r(1,:,i))==0
%             fl_r(:,:,i)=[];
%         else
%             i=i+1;
%         end
%     end
%     c=0;
%     for j=1:length(p(:,6))
%         if p(j,6)==1
%             c=c+1;
%             x=fl_r(2,:,c);
%             y1=fl_r(1,:,c);
%             y2=p(j,2)*ones(1,length(y1));
%             for i=1:length(y1)
%                 yf(i)=min(y1(i),y2(i));
%             end
%             cp=cp+1;
%             xtip2(cp,:)=x;
%             yip2(cp,:)=yf;
%             ytip2(cp,:)=y1;
%             ayip2=[ayip2;yf];
%         end
%     end
%     axip2=x;
% end
% %2in2
% if isempty(fm_r)==0
%     i=1;
%     while i<=length(fm_r(1,1,:))
%         if sum(fm_r(1,:,i))==0
%             fm_r(:,:,i)=[];
%         else
%             i=i+1;
%         end
%     end
%     c=0;
%     for j=1:length(p(:,6))
%         if p(j,6)==2
%             c=c+1;
%             x=fm_r(2,:,c);
%             y1=fm_r(1,:,c);
%             y2=p(j,2)*ones(1,length(y1));
%             for i=1:length(y1)
%                 yf(i)=min(y1(i),y2(i));
%             end
%             cp=cp+1;
%             xtip2(cp,:)=x;
%             yip12(cp,:)=yf;
%             ytip2(cp,:)=y1;
%             ayip2=[ayip2;yf];
%         end
%     end
%     axip2=x;
% end
% %3in2
% if isempty(fh_r)==0
%     i=1;
%     while i<=length(fh_r(1,1,:))
%         if sum(fh_r(1,:,i))==0
%             fh_r(:,:,i)=[];
%         else
%             i=i+1;
%         end
%     end
%     c=0;
%     for j=1:length(p(:,6))
%         if p(j,6)==3
%             c=c+1;
%             x=fh_r(2,:,c);
%             y1=fh_r(1,:,c);
%             y2=p(j,2)*ones(1,length(y1));
%             for i=1:length(y1)
%                 yf(i)=min(y1(i),y2(i));
%             end
%             cp=cp+1;
%             xtip2(cp,:)=x;
%             yip2(cp,:)=yf;
%             ytip2(cp,:)=y1;
%             ayip2=[ayip2;yf];
%         end
%     end
%     axip2=x;
% end











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ayip1=[];
cp=0;
%1ip2
if isempty(fl_x)==0
    i=1;
    while i<=length(fl_x(1,1,:))
        if sum(fl_x(1,:,i))==0
            fl_x(:,:,i)=[];
        else
            i=i+1;
        end
    end
    c=0;
    for j=1:length(p(:,6))
        if p(j,5)==1
            c=c+1;
            x=fl_x(2,:,c);
            y1=fl_x(1,:,c);
            y2=p(j,1)*ones(1,length(y1));
            for i=1:length(y1)
                yf(i)=min(y1(i),y2(i));
            end
            cp=cp+1;
            xtip1(cp,:)=x;
            yip1(cp,:)=yf;
            ytip1(cp,:)=y1;
            ayip1=[ayip1;yf];
        end
    end
    axip2=x;
end

%2ip2
if isempty(fm_x)==0
    i=1;
    while i<=length(fm_x(1,1,:))
        if sum(fm_x(1,:,i))==0
            fm_x(:,:,i)=[];
        else
            i=i+1;
        end
    end
    c=0;
    for j=1:length(p(:,6))
        if p(j,5)==2
            c=c+1;
            x=fm_x(2,:,c);
            y1=fm_x(1,:,c);
            y2=p(j,1)*ones(1,length(y1));
            for i=1:length(y1)
                yf(i)=min(y1(i),y2(i));
            end
            cp=cp+1;
            xtip2(cp,:)=x;
            yip2(cp,:)=yf;
            ytip2(cp,:)=y1;
            ayip1=[ayip1;yf];
        end
    end
    axip2=x;
end
%3ip2
if isempty(fh_x)==0
    i=1;
    while i<=length(fh_x(1,1,:))
        if sum(fh_x(1,:,i))==0
            fh_x(:,:,i)=[];
        else
            i=i+1;
        end
    end
    c=0;
    for j=1:length(p(:,6))
        if p(j,5)==3
            c=c+1;
            x=fh_x(2,:,c);
            y1=fh_x(1,:,c);
            y2=p(j,1)*ones(1,length(y1));
            for i=1:length(y1)
                yf(i)=min(y1(i),y2(i));
            end
            cp=cp+1;
            xtip2(cp,:)=x;
            yip2(cp,:)=yf;
            ytip2(cp,:)=y1;
            ayip1=[ayip1;yf];
        end
    end
    axip2=x;
end
%4ip2
if isempty(fvh_x)==0
    i=1;
    while i<=length(fvh_x(1,1,:))
        if sum(fvh_x(1,:,i))==0
            fvh_x(:,:,i)=[];
        else
            i=i+1;
        end
    end
    c=0;
    for j=1:length(p(:,6))
        if p(j,5)==4
            c=c+1;
            x=fvh_x(2,:,c);
            y1=fvh_x(1,:,c);
            y2=p(j,1)*ones(1,length(y1));
            for i=1:length(y1)
                yf(i)=min(y1(i),y2(i));
            end
            cp=cp+1;
            xtip2(cp,:)=x;
            yip2(cp,:)=yf;
            ytip2(cp,:)=y1;
            ayip1=[ayip1;yf];
        end
    end
    axip2=x;
end


ayip2=[];
cp=0;
%1ip2
if isempty(fl_r)==0
    i=1;
    while i<=length(fl_r(1,1,:))
        if sum(fl_r(1,:,i))==0
            fl_r(:,:,i)=[];
        else
            i=i+1;
        end
    end
    c=0;
    for j=1:length(p(:,6))
        if p(j,6)==1
            c=c+1;
            x=fl_r(2,:,c);
            y1=fl_r(1,:,c);
            y2=p(j,2)*ones(1,length(y1));
            for i=1:length(y1)
                yf(i)=min(y1(i),y2(i));
            end
            cp=cp+1;
            xtip2(cp,:)=x;
            yip2(cp,:)=yf;
            ytip2(cp,:)=y1;
            ayip2=[ayip2;yf];
        end
    end
    axip2=x;
end

%2ip2
if isempty(fm_r)==0
    i=1;
    while i<=length(fm_r(1,1,:))
        if sum(fm_r(1,:,i))==0
            fm_r(:,:,i)=[];
        else
            i=i+1;
        end
    end
    c=0;
    for j=1:length(p(:,6))
        if p(j,6)==2 && length(fm_r(1,1,:))~=0
            c=c+1;
            x=fm_r(2,:,c);
            y1=fm_r(1,:,c);
            y2=p(j,2)*ones(1,length(y1));
            for i=1:length(y1)
                yf(i)=min(y1(i),y2(i));
            end
            cp=cp+1;
            xtip2(cp,:)=x;
            yip2(cp,:)=yf;
            ytip2(cp,:)=y1;
            ayip2=[ayip2;yf];
        end
    end
    axip2=x;
end
%3ip2
if isempty(fh_r)==0
    i=1;
    while i<=length(fh_r(1,1,:))
        if sum(fh_r(1,:,i))==0
            fh_r(:,:,i)=[];
        else
            i=i+1;
        end
    end
    c=0;
    for j=1:length(p(:,6))
        if p(j,6)==3 && length(fh_r(1,1,:))~=0
            c=c+1;
            x=fh_r(2,:,c);
            y1=fh_r(1,:,c);
            y2=p(j,2)*ones(1,length(y1));
            for i=1:length(y1)
                yf(i)=min(y1(i),y2(i));
            end
            cp=cp+1;
            xtip2(cp,:)=x;
            yip2(cp,:)=yf;
            ytip2(cp,:)=y1;
            ayip2=[ayip2;yf];
        end
    end
    axip2=x;
end




%%% %%% %%% %%%%%%%%%%
ay=[];
cp=0;
%1out
if isempty(fl_v)==0
    i=1;
    while i<=length(fl_v(1,1,:))
        if sum(fl_v(1,:,i))==0
            fl_v(:,:,i)=[];
        else
            i=i+1;
        end
    end
    c=0;
    for j=1:length(p(:,4))
        if p(j,4)==1
            c=c+1;
            x=fl_v(2,:,c);
            y1=fl_v(1,:,c);
            y2=p(j,3)*ones(1,length(y1));
            for i=1:length(y1)
                yf(i)=min(y1(i),y2(i));
            end
            cp=cp+1;
            xt(cp,:)=x;
            y(cp,:)=yf;
            yt(cp,:)=y1;
            ay=[ay;yf];
        end
    end
    ax=x;
end
% 2out
if isempty(fm_v)==0
    i=1;
    while i<=length(fm_v(1,1,:))
        if sum(fm_v(1,:,i))==0
            fm_v(:,:,i)=[];
        else
            i=i+1;
        end
    end
    c=0;
    for j=1:length(p(:,4))
        if p(j,4)==2
            c=c+1;
            x=fm_v(2,:,c);
            y1=fm_v(1,:,c);
            y2=p(j,3)*ones(1,length(y1));
            for i=1:length(y1)
                yf(i)=min(y1(i),y2(i));
            end
            cp=cp+1;
            xt(cp,:)=x;
            y(cp,:)=yf;
            yt(cp,:)=y1;
            ay=[ay;yf];
        end
    end
    ax=x;
end
% 3out
if isempty(fh_v)==0
    i=1;
    while i<=length(fh_v(1,1,:))
        if sum(fh_v(1,:,i))==0
            fh_v(:,:,i)=[];
        else
            i=i+1;
        end
    end
    c=0;
    for j=1:length(p(:,4))
        if p(j,4)==3
            c=c+1;
            x=fh_v(2,:,c);
            y1=fh_v(1,:,c);
            y2=p(j,3)*ones(1,length(y1));
            for i=1:length(y1)
                yf(i)=min(y1(i),y2(i));
            end
            cp=cp+1;
            xt(cp,:)=x;
            y(cp,:)=yf;
            yt(cp,:)=y1;
            ay=[ay;yf];
        end
    end
    ax=x;
end
% 4out
if isempty(fvh_v)==0
    i=1;
    while i<=length(fvh_v(1,1,:))
        if sum(fvh_v(1,:,i))==0
            fvh_v(:,:,i)=[];
        else
            i=i+1;
        end
    end
    c=0;
    for j=1:length(p(:,4))
        if p(j,4)==4
            c=c+1;
            x=fvh_v(2,:,c);
            y1=fvh_v(1,:,c);
            y2=p(j,3)*ones(1,length(y1));
            for i=1:length(y1)
                yf(i)=min(y1(i),y2(i));
            end
            cp=cp+1;
            xt(cp,:)=x;
            y(cp,:)=yf;
            yt(cp,:)=y1;
            ay=[ay;yf];
        end
    end
    ax=x;
end

% Subplot
figure
for i=1:cp
%     subplot(cp,3,(1+3*(i-1)))
%     area(xtip1(i,:),yip1(i,:))
%     hold on
%     plot(xtip1(i,:),ytip1(i,:))
%     
%     subplot(cp,3,(2+3*(i-1)))
%     area(xtip2(i,:),yip2(i,:))
%     hold on
%     plot(xtip12(i,:),ytip2(i,:))
    
    subplot(cp,3,(3+3*(i-1)))
    area(xt(i,:),y(i,:))
    hold on
    plot(xt(i,:),yt(i,:))
end

% centroid
[ray,cay]=size(ay);
for i=1:cay
    ay(ray+1,i)=max(ay(1:ray,i));
end
figure
area(ax,ay(ray+1,:))
[xc,yc]=cntr(ax,ay(ray+1,:))
hold on
scatter(xc,yc,'filled','k')
msgbox(sprintf('Speed = %2.3g',xc))
