%%question 2-2
close all;
clear;
clc;

%part4

t = 0 : 0.01 : 4;

x = sin(2*pi.*t)+sin(100*pi.*t);

L = 300;
S = reflection(x,L);

for i = 1:1:30
m = i;
xf = moving_average_filter(S,m,L);
plot(t,xf);
hold on
plot(t,x);
pause(0.2);
hold off
title('m =30');
end

%%
%part 5
close all;
clear;
clc;

fs = 100;
Ts = 1/fs;
t = 0 : Ts : 4;
f0 = 1;
Xs = sin(2*pi*f0*t);

a = 0.2;
L = 401;
type1 = 'normal';
N1 = noise(a,L,type1);
Xn1 = Xs+N1;

L = 250;
S = reflection(Xn1,L);

for i = 1:1:5
m = i;
xf = moving_average_filter(S,m,L);
plot(t,xf);
hold on
plot(t,Xn1);
pause(0.2);
title('Xn1 = Xs+N1/m=5');
legend('without noise','with noise');
hold off
end

%%
%part 5
close all;
clear;
clc;

fs = 100;
Ts = 1/fs;
t = 0 : Ts : 4;
f0 = 1;
Xs = sin(2*pi*f0*t);

a = 0.2;
L = 401;
type2 = 'uniform';
N2 = noise(a,L,type2);
Xn2 = Xs+N2;

L = 250;
S = reflection(Xn2,L);

for i = 1:1:15
m = i;
xf = moving_average_filter(S,m,L);
plot(t,xf);
hold on
plot(t,Xn2);
pause(0.2);
title('Xn2 = Xs+N2/m=15');
legend('without noise','with noise');
hold off
end
%%
%part 5
close all;
clear;
clc;

fs = 100;
Ts = 1/fs;
t = 0 : Ts : 4;
f0 = 1;
Xs = sin(2*pi*f0*t);

a = 0.2;
L = 401;
type1 = 'normal';
N1 = noise(a,L,type1);
Xn1 = Xs.*N1;

L = 250;
S = reflection(Xn1,L);

for i = 1:1:15
m = i;
xf = moving_average_filter(S,m,L);
plot(t,xf);
hold on
plot(t,Xn1);
pause(0.2);
title('Xn1 = Xs.*N1/m=15');
legend('without noise','with noise');
hold off
end
%%
%part 5
close all;
clear;
clc;

fs = 100;
Ts = 1/fs;
t = 0 : Ts : 4;
f0 = 1;
Xs = sin(2*pi*f0*t);

a = 0.2;
L = 401;
type2 = 'uniform';
N2 = noise(a,L,type2);
Xn2 = Xs.*N2;

L = 250;
S = reflection(Xn2,L);

for i = 1:1:10
m = i;
xf = moving_average_filter(S,m,L);
plot(t,xf);
hold on
plot(t,Xn2);
pause(0.2);
title('Xn2 = Xs+N2/m=10');
legend('without noise','with noise');
hold off
end
%%
%part 5
close all;
clear;
clc;

fs = 100;
Ts = 1/fs;
t = 0 : Ts : 4;
f0 = 1;
Xs = sin(2*pi*f0*t);

a = 0.2;
L = 401;
type1 = 'normal';
N1 = noise(a,L,type1);
Xn1 = Xs.*(1+N1);

L = 250;
S = reflection(Xn1,L);

for i = 1:1:10
m = i;
xf = moving_average_filter(S,m,L);
plot(t,xf);
hold on
plot(t,Xn1);
pause(0.2);
title('Xs.*(1+N1)/m=10');
legend('without noise','with noise');
hold off
end
%%
%part 5
close all;
clear;
clc;

fs = 100;
Ts = 1/fs;
t = 0 : Ts : 4;
f0 = 1;
Xs = sin(2*pi*f0*t);

a = 0.2;
L = 401;
type2 = 'uniform';
N2 = noise(a,L,type2);
Xn2 = Xs.*(1+N2);

L = 250;
S = reflection(Xn2,L);

for i = 1:1:15
m = i;
xf = moving_average_filter(S,m,L);
plot(t,xf);
hold on
plot(t,Xn2);
pause(0.2);
title('Xn2 = Xs.*(1+N2)/m=15');
legend('without noise','with noise');
hold off
end
%%
%part 8
close all;
clear;
clc;

t = 0 : 0.01 : 4;

x = sin(2*pi.*t)+sin(100*pi.*t);

L = 300;
S = reflection(x,L);

for i = 1:1:250
m = i;
xf = median_filter(S,m,L);
plot(t,xf);
hold on
plot(t,x);
pause(0.2);
hold off
title('m =250');
end

%%
%part 9
close all;
clear;
clc;

fs = 100;
Ts = 1/fs;
t = 0 : Ts : 4;
f0 = 1;
Xs = sin(2*pi*f0*t);

a = 0.2;
L = 401;
type1 = 'normal';
N1 = noise(a,L,type1);
Xn1 = Xs+N1;

L = 100;
S = reflection(Xn1,L);

for i = 1:1:12
m = i;
xf = median_filter(S,m,L);
plot(t,xf);
hold on
plot(t,Xn1);
pause(0.2);
title('Xn1 = Xs+N1/m=12');
legend('without noise','with noise');
hold off
end
%%
%part 9
close all;
clear;
clc;

fs = 100;
Ts = 1/fs;
t = 0 : Ts : 4;
f0 = 1;
Xs = sin(2*pi*f0*t);

a = 0.2;
L = 401;
type2 = 'uniform';
N2 = noise(a,L,type2);
Xn2 = Xs+N2;

L = 250;
S = reflection(Xn2,L);

for i = 1:1:13
m = i;
xf = median_filter(S,m,L);
plot(t,xf);
hold on
plot(t,Xn2);
pause(0.2);
title('Xn2 = Xs+N2/m=13');
legend('without noise','with noise');
hold off
end
%%
%part 9
close all;
clear;
clc;

fs = 100;
Ts = 1/fs;
t = 0 : Ts : 4;
f0 = 1;
Xs = sin(2*pi*f0*t);

a = 0.2;
L = 401;
type1 = 'normal';
N1 = noise(a,L,type1);
Xn1 = Xs.*N1;

L = 250;
S = reflection(Xn1,L);

for i = 1:1:15
m = i;
xf = median_filter(S,m,L);
plot(t,xf);
hold on
plot(t,Xn1);
pause(0.2);
title('Xn1 = Xs.*N1/m=15');
legend('without noise','with noise');
hold off
end
%%
%part 9
close all;
clear;
clc;

fs = 100;
Ts = 1/fs;
t = 0 : Ts : 4;
f0 = 1;
Xs = sin(2*pi*f0*t);

a = 0.2;
L = 401;
type2 = 'uniform';
N2 = noise(a,L,type2);
Xn2 = Xs.*N2;

L = 250;
S = reflection(Xn2,L);

for i = 1:1:10
m = i;
xf = median_filter(S,m,L);
plot(t,xf);
hold on
plot(t,Xn2);
pause(0.2);
title('Xn2 = Xs+N2/m=10');
legend('without noise','with noise');
hold off
end
%%
%part 9
close all;
clear;
clc;

fs = 100;
Ts = 1/fs;
t = 0 : Ts : 4;
f0 = 1;
Xs = sin(2*pi*f0*t);

a = 0.2;
L = 401;
type1 = 'normal';
N1 = noise(a,L,type1);
Xn1 = Xs.*(1+N1);

L = 250;
S = reflection(Xn1,L);

for i = 1:1:10
m = i;
xf = median_filter(S,m,L);
plot(t,xf);
hold on
plot(t,Xn1);
pause(0.2);
title('Xs.*(1+N1)/m=10');
legend('without noise','with noise');
hold off
end
%%
%part 9
close all;
clear;
clc;

fs = 100;
Ts = 1/fs;
t = 0 : Ts : 4;
f0 = 1;
Xs = sin(2*pi*f0*t);

a = 0.2;
L = 401;
type2 = 'uniform';
N2 = noise(a,L,type2);
Xn2 = Xs.*(1+N2);

L = 250;
S = reflection(Xn2,L);

for i = 1:1:15
m = i;
xf = moving_average_filter(S,m,L);
plot(t,xf);
hold on
plot(t,Xn2);
pause(0.2);
title('Xn2 = Xs.*(1+N2)/m=15');
legend('without noise','with noise');
hold off
end
%%
%part 10
close all;
clear;
clc;

filename = '11.wav';
[y,Fs] = audioread(filename);
a = 1;
L = 947200;
type2 = 'uniform';
N2 = noise(a,L,type2);
Xn2 = y.*(1+N2);
%sound(Xn2,Fs);
audiowrite('D:\Zahra\new Xn21.a1.wav',Xn2,Fs);


a = 0.1;
L = 947200;
type1 = 'uniform';
N1 = noise(a,L,type1);
Xn1 = y+N1;
%sound(Xn1,Fs);
audiowrite('D:\Zahra\new Xn11.a0.1.wav',Xn1,Fs);

a = 0.1;
L = 947200;
type1 = 'normal';
N1 = noise(a,L,type1);
Xn1 = y+N1;
%sound(Xn1,Fs);
audiowrite('D:\Zahra\new Xn12.a0.1.wav',Xn1,Fs);

a = 1;
L = 947200;
type2 = 'normal';
N2 = noise(a,L,type2);
Xn2 = y.*(1+N2);
%sound(Xn2,Fs);
audiowrite('D:\Zahra\new Xn22.a1.wav',Xn2,Fs);


%%
function xf = moving_average_filter(x,m,L)
[n y] = size(x);
sum = 0;
xf = zeros(n,y-2*L);

for i =L+1:1:y-L
    for j = i-m:1:i+m
        sum = sum + x(1,j);
    end
    xf(1,i-L) = sum/(2*m+1);
    sum = 0;
end
end

function xf = moving_average_filtersound(x,m,L)
[n y] = size(x);
sum = 0;
xf = zeros(n,y-2*L);

for k = 1:1:n
for i =L+1:1:y-L
    for j = i-m:1:i+m
        sum = sum + x(k,j);
    end
    xf(1,i-L) = sum/(2*m+1);
    sum = 0;
end
end
end

function xf = median_filter(x,m,L)
[n y] = size(x);
za = zeros(1,2*m+1);
xf = zeros(n,y-2*L);
c = 1;

for i=L+1:1:y-L
    for j = i-m:1:i+m
        za(1,c) = x(1,j);
        c = c+1;
    end
    xf(1,i-L) = median(za);
    za = zeros(1,2*m+1);
    c = 1;
end
end

function S = reflection(I,L)
[x y]=size(I);
R=zeros(x,y+L);
R(1,1:L)=I(1,L:-1:1);
R(1,L+1:y+L)=I(1,1:y);

S = zeros(x,y+2*L);
S(1,1:y+L) = R(1,1:y+L);
S(1,y+L+1:y+2*L) = R(1,y+L:-1:y+1);
R=uint8(S);
end



function N = noise(a,L,type)
s1 = 'normal';
s2 = 'uniform';
if(strcmp(s1,type))
    N = normrnd(0,a^2,[L,2]);
else if(strcmp(s2,type))
    N =-abs(a) + (2*abs(a)).*rand(L,2);
end
end
end


