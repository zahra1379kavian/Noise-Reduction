close all;
clear;
clc;

%part 1
img = imread('D:\Zahra\EngineeringMath\project\image.png');
J = imresize(img, 0.5);

figure
imshow(img);
title('image origina');

%%
%part 2
close all;
clear;
clc;

img = imread('D:\Zahra\EngineeringMath\project\image.png');

k1 = 1;
k2 = 1;

for i = 1:2:573
    for j = 1:2:1028
        img2(k1,k2,[1 2 3]) = img(i,j,[1 2 3]);
        k2 = k2+1;
    end
    k1 = k1+1;
    k2 = 1;
end

img2=uint8(img2);

%low noise
figure
data = im2double(img2);
yl = awgn(data,20);
imshow(yl);
title('low noise');
%median noise
figure
ym = awgn(data,10);
imshow(ym);
title('median noise');
%high noise
figure
yh = awgn(data,1);
imshow(yh);
title('high noise');

%minimom SNR???
figure
y = awgn(data,0.0000000000000009);
imshow(y);

%%
%part 3
close all;
clear;
clc;

img = imread('D:\Zahra\EngineeringMath\project\image.png');
J = imnoise(img,'salt & pepper',0.2);
imshow(J)
im2double(img);

%%
%part 4
close all;
clear;
clc;

img = imread('D:\Zahra\EngineeringMath\project\image.png');

k1 = 1;
k2 = 1;

for i = 1:2:573
    for j = 1:2:1028
        img2(k1,k2,[1 2 3]) = img(i,j,[1 2 3]);
        k2 = k2+1;
    end
    k1 = k1+1;
    k2 = 1;
end

%data = im2double(img2);

R = reflectionR(img2,50);
LL = reflectionL(R,50);
U = reflectionU(LL,50);
D = reflectionD(U,50);
D = im2double(D);
L = 40;


figure
m1 = 30;
m2 = 30;
filtered_image1 = img_moving_average_filter(D,m1,m2,L);
subplot(2,1,1);
imshow(filtered_image1);
%title('lenght>width');
title('30.30');

m1 = 5;
m2 = 5;
filtered_image1 = img_moving_average_filter(D,m1,m2,L);
subplot(2,1,2);
imshow(filtered_image1);
%title('lenght<width');
title('5.5');

figure
m1 = 5;
m2 = 5;

filtered_image2 =  img_median_filter(D,m1,m2,L);
subplot(2,1,1);
imshow(filtered_image2);
%%title('lenght<width//10.10');
title('5.5');

m1 = 30;
m2 = 30;
filtered_image2 =  img_median_filter(D,m1,m2,L);
subplot(2,1,2);
imshow(filtered_image2);
%%title('lenght>width');
title('30.30');



%%
%part 4
close all;
clear;
clc;

img = imread('D:\Zahra\EngineeringMath\project\image.png');

k1 = 1;
k2 = 1;

for i = 1:2:573
    for j = 1:2:1028
        img2(k1,k2,[1 2 3]) = img(i,j,[1 2 3]);
        k2 = k2+1;
    end
    k1 = k1+1;
    k2 = 1;
end

img2=uint8(img2);

%low noise
data = im2double(img2);
ym = awgn(data,10);


R = reflectionR(ym,50);
LL = reflectionL(R,50);
U = reflectionU(LL,50);
D = reflectionD(U,50);
D = im2double(D);

L = 40;

figure
m1 = 10;
m2 = 10;
filtered_image1 = img_moving_average_filter(D,m1,m2,L);
subplot(2,1,1);
imshow(filtered_image1);
%title('lenght>width');
title('10.10');

m1 = 10;
m2 = 10;

filtered_image2 =  img_median_filter(D,m1,m2,L);
subplot(2,1,2);
imshow(filtered_image2);
%%title('lenght<width//10.10');
title('10.10');

%%
img = imread('coin.png');
J = imnoise(img,'salt & pepper',0.2);
R = reflectionR(J,50);
LL = reflectionL(R,50);
U = reflectionU(LL,50);
D = reflectionD(U,50);
D = im2double(D);

L = 40;

figure
m1 = 7;
m2 = 7;
filtered_image1 = img_moving_average_filter(D,m1,m2,L);
subplot(2,1,1);
imshow(filtered_image1);
%title('lenght>width');
title('7.7');

m1 = 10;
m2 = 10;

filtered_image2 =  img_median_filter(D,m1,m2,L);
subplot(2,1,2);
imshow(filtered_image2);
%%title('lenght<width//10.10');
title('10.10');



%%
%Program for one sided image reflection along a Line
function R = reflectionR(I,L)
[x y z]=size(I);
R=zeros(x,y+L,z);
R(:,1:y,:)=I(:,1:y,:);
R(:,y+1:y+L,:)=I(:,y:-1:y-L+1,:);
R=uint8(R);
end


function LL = reflectionL(I,L)
[x y z]=size(I);
R=zeros(x,y+L,z);
R(:,1:L,:)=I(:,L:-1:1,:);
R(:,L+1:y+L,:)=I(:,1:y,:);
LL = R;
LL=uint8(R);
end

function U = reflectionU(I,L)
[x y z]=size(I);
R=zeros(x+L,y,z);
R(1:L,:,:)=I(L:-1:1,:,:);
R(L+1:x+L,:,:)=I(1:x,:,:);
U = R;
U=uint8(R);
end

function D = reflectionD(I,L)
[x y z]=size(I);
R=zeros(x+L,y,z);
R(1:x,:,:)=I(1:x,:,:);
R(x+1:x+L,:,:)=I(x:-1:x-L+1,:,:);
D = R;
D=uint8(R);
end

%program for filtering image
function filtered_image = img_moving_average_filter(img,m1,m2,L)
X = img(:,:,1);
[n y] = size(X);
sum = 0;
X1 = zeros(n-2*L,y-2*L);

for k = L+1:1:n-L
    for g = L+1:1:y-L
for i =k-m1:1:k+m1
    for j = g-m2:1:g+m2
        sum = sum + X(i,j);
    end
    X1(k-L,g-L) = sum/(m1*m2);
    sum = 0;
end
    end
end

X = img(:,:,2);
[n y] = size(X);
sum = 0;
X2 = zeros(n-2*L,y-2*L);

for k = L+1:1:n-L
    for g = L+1:1:y-L
for i =k-m1:1:k+m1
    for j = g-m2:1:g+m2
        sum = sum + X(i,j);
    end
    X2(k-L,g-L) = sum/(m1*m2);
    sum = 0;
end
    end
end


X = img(:,:,3);
[n y] = size(X);
sum = 0;
X3 = zeros(n-2*L,y-2*L);

for k = L+1:1:n-L
    for g = L+1:1:y-L
for i =k-m1:1:k+m1
    for j = g-m2:1:g+m2
        sum = sum + X(i,j);
    end
    X3(k-L,g-L) = sum/(m1*m2);
    sum = 0;
end
    end
end


[x y z]=size(img);
R=zeros(x-2*L,y-2*L,z);

R(:,:,1) = X1;
R(:,:,2) = X2;
R(:,:,3) = X3;

filtered_image = R;
end

function R = img_median_filter(img,m1,m2,L)
X = img(:,:,1);
[n y] = size(X);

X1 = zeros(n-2*L,y-2*L);
za = zeros(2*m1+1,2*m2+1);
c1 = 1;
c2 = 1;

for k = L+1:1:n-L
    for g = L+1:1:y-L
for i =k-m1:1:k+m1
    for j = g-m2:1:g+m2
        za(c1,c2)= X(i,j);
        c2 = c2+1;
        c2 = 1;
    end
    c1 = c1 + 1;
end
c1 = 1;
X1(k-L,g-L) = median(za,'all');
    end
end

X = img(:,:,2);
[n y] = size(X);
sum = 0;
X2 = zeros(n-2*L,y-2*L);
za = zeros(2*m1+1,2*m2+1);
c1 = 1;
c2 = 1;

for k = L+1:1:n-L
    for g = L+1:1:y-L
for i =k-m1:1:k+m1
    for j = g-m2:1:g+m2
        za(c1,c2)= X(i,j);
        c2 = c2+1;
    end
    c1 = c1 + 1;
	c2 = 1;
end
c1 = 1;
X2(k-L,g-L) = median(za,'all');
    end
end

X = img(:,:,3);
[n y] = size(X);
sum = 0;
X3 = zeros(n-2*L,y-2*L);
za = zeros(2*m1+1,2*m2+1);
c1 = 1;
c2 = 1;

for k = L+1:1:n-L
    for g = L+1:1:y-L
for i =k-m1:1:k+m1
    for j = g-m2:1:g+m2
        za(c1,c2)= X(i,j);
        c2 = c2+1;
    end
    c1 = c1 + 1;
	c2 = 1;
end
c1 = 1;
X3(k-L,g-L) = median(za,'all');
    end
end


[x y z]=size(img);
R=zeros(x-2*L,y-2*L,z);

R(:,:,1) = X1;
R(:,:,2) = X2;
R(:,:,3) = X3;

filtered_image = R;

end

