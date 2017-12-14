clc;
clear all;
close all;
%% image data
% I(:,:,1) = imread('/home/hulio/Downloads/proj3_seq1/sphere2.pgm');
% I(:,:,2) = imread('/home/hulio/Downloads/proj3_seq1/sphere3.pgm');
% I(:,:,3) = imread('/home/hulio/Downloads/proj3_seq1/sphere4.pgm');
% I(:,:,4) = imread('/home/hulio/Downloads/proj3_seq1/sphere5.pgm');
% I(:,:,5) = imread('/home/hulio/Downloads/proj3_seq1/sphere6.pgm');

I(:,:,1) = imread('/home/hulio/Downloads/proj3_seq3/1.pgm');
I(:,:,2) = imread('/home/hulio/Downloads/proj3_seq3/2.pgm');
I(:,:,3) = imread('/home/hulio/Downloads/proj3_seq3/3.pgm');
I(:,:,4) = imread('/home/hulio/Downloads/proj3_seq3/4.pgm');
I(:,:,5) = imread('/home/hulio/Downloads/proj3_seq3/5.pgm');
%% symbolic representation of the equations
syms x y t a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20 F;

F = a1 + a2*x + a3*y + a4*t + a5*(x)^2 + a6*x*y + a7*(y)^2 + a8*y*t + a9*(t)^2 + a10*x*t + a11*(x)^3 + a12*((x)^2)*y + a13*x*(y)^2 + a14*(y)^3 + a15*((y)^2)*t + a16*y*(t)^2 + a17*(t)^3 + a18*((x)^2)*t + a19*x*(t)^2 + a20*x*y*t;
Ix = jacobian(F,x);
Iy = jacobian(F,y);
Ixx = jacobian(Ix,x);
Iyy = jacobian(Iy,y);
Ixy = jacobian(Ix,y);
Iyx = jacobian(Iy,x);
It = jacobian(F,t);
Itx = jacobian(It,x);
Ity = jacobian(It,y);
Iyt = jacobian(Iy,t);
Ixt = jacobian(Ix,t);
Itt = jacobian(It,t);
%% initialize parameters
%%%% all of the above are the complete equations which are solved per block
%%%% to get the coefficients
n = 3;  %% NxN = [3x3]
T = 5;  %% T = 5, i.e [NxNxT] = [3x3x5]
Y = n;
X = n;
[mm,nn,~] = size(I);

ii = 1;
%% precalculating all the coefficients per NxNxT voxel so that only the values can be used later as placeholders and it is not intensive computationally
for tt = -(round(T/2) - 1):(round(T/2) - 1)
    for yy = -(round(Y/2)-1):(round(Y/2) - 1)
        for xx = -(round(X/2)-1):(round(X/2)-1)
            
            
        A = [subs(subs(subs(Ix,x,xx),y,yy),t,tt) subs(subs(subs(Iy,x,xx),y,yy),t,tt); subs(subs(subs(Ixx,x,xx),y,yy),t,tt) subs(subs(subs(Ixy,x,xx),y,yy),t,tt); subs(subs(subs(Iyx,x,xx),y,yy),t,tt) subs(subs(subs(Iyy,x,xx),y,yy),t,tt); subs(subs(subs(Itx,x,xx),y,yy),t,tt) subs(subs(subs(Ity,x,xx),y,yy),t,tt)];
        b = -[subs(subs(subs(It,x,xx),y,yy),t,tt); subs(subs(subs(Ixt,x,xx),y,yy),t,tt); subs(subs(subs(Iyt,x,xx),y,yy),t,tt); subs(subs(subs(Itt,x,xx),y,yy),t,tt)];
        %%% A and b
        BC(ii,:) = {A , b, [xx yy tt]};
        ii = ii + 1;
        end
    end
end

%% Optical Flow calculated per NxNxT voxel and then segregated into U and V [X and Y] values
% opti_fxn = @(image_block) opti_flow_poly(image_block.data,BC,100);
% UV = blockproc(I,[n n],opti_fxn);

%%  Alternate way to get optical flow values
[M,N,O] = size(I);
ind = 1:O*2;
UV = zeros(M,N,O*2);

for i = round(n/2) : M - (n - round(n/2))
for j = round(n/2 ): N - (n - round(n/2))
    
    UUVV = opti_flow_poly(I(i-floor(n/2):i+floor(n/2),j-floor(n/2):j+floor(n/2),:),BC,100);
%     UV(i-floor(n/2):i+floor(n/2),j-floor(n/2):j+floor(n/2),ind) = UUVV(:,:,ind);
    UV(i,j,ind) = UUVV(round(n/2),round(n/2),ind);
    
    lin = sub2ind(I(:,:,1),i,j);
    
    if (mod(lin,10) == 0)
    (lin / numel(I(:,:,1)))*100
    end
    
end
end


U = UV(:,:,1:5);
V = UV(:,:,6:10);

quiv_plot(U(:,:,3),V(:,:,3),30);

