clc
clear all


% Logistic 初始参数
x0 = 0.54180; 
u = 3.678927;
N = 1000;


% Chen 初始参数
y = [0, 2, 9];

img = double(imread('lena256.bmp'));
[h, w, c] = size(img);
r = reshape(img(:,:,1)', 1, []);
g = reshape(img(:,:,2)', 1, []);
b = reshape(img(:,:,3)', 1, []);
m = [r, g, b];
% Logistic 置乱
[idx, cm] = LgsTs((h*w*c), u, x0, N);
m = m(:, idx);
r_lgs_cm = cm(1:h*w);
g_lgs_cm = cm(h*w+1:2*h*w);
b_lgs_cm = cm(2*h*w+1:end);

% 解码种类
r_ca_decode = floor(mod(sum(r_lgs_cm), 8));
g_ca_decode = floor(mod(sum(g_lgs_cm), 8));
b_ca_decode = floor(mod(sum(b_lgs_cm), 8));

r_lgs_cm = floor(mod(r_lgs_cm*10^8, 256));
g_lgs_cm = floor(mod(g_lgs_cm*10^8, 256));
b_lgs_cm = floor(mod(b_lgs_cm*10^8, 256));
m = bitxor(m, floor(mod(cm*10^8, 256)));

% Logistic 矩阵DNA编码
[r_lgs_flags, r_lgs_encode] = dna_encode(r_lgs_cm, 'dynamic');
[g_lgs_flags, g_lgs_encode] = dna_encode(g_lgs_cm, 'dynamic');
[b_lgs_flags, b_lgs_encode] = dna_encode(b_lgs_cm, 'dynamic');

r = reshape(m(1:h*w), h, w)'; 
g = reshape(m(h*w+1:2*h*w), h, w)'; 
b = reshape(m(2*h*w+1:3*h*w), h, w)';
% 原图动态DNA编码
[r_flags, r_encode] = dna_encode(r, 'dynamic');
[g_flags, g_encode] = dna_encode(g, 'dynamic');
[b_flags, b_encode] = dna_encode(b, 'dynamic');
% Lorenz 扩散矩阵
[t,y1] = ode45('Lorenz', [0,1200], y);
y1 = mod(round(y1*10^8), 256);
c_chaos = y1(:,1); g_chaos = y1(:,2); b_chaos = y1(:,3);
c_chaos=c_chaos(1:h*w)'; g_chaos=g_chaos(1:h*w)'; b_chaos=b_chaos(1:h*w)';

c_chaos = reshape(c_chaos, h, w);
g_chaos = reshape(g_chaos, h, w);
b_chaos = reshape(b_chaos, h, w);

% Lorenz 扩散矩阵动态DNA编码
[~, r_chaos_encode] = dna_encode(c_chaos, 'static', r_ca_decode);
[~, g_chaos_encode] = dna_encode(g_chaos, 'static', g_ca_decode);
[~, b_chaos_encode] = dna_encode(b_chaos, 'static', b_ca_decode);

% DNA 加法
r_dna_add = DnaOper(r_encode, r_lgs_encode, 'add');
g_dna_add = DnaOper(g_encode, g_lgs_encode, 'add');
b_dna_add = DnaOper(b_encode, b_lgs_encode, 'add');

% DNA 异或
r_dna_xor = DnaOper(r_dna_add, r_chaos_encode, 'xor');
g_dna_xor = DnaOper(g_dna_add, g_chaos_encode, 'xor');
b_dna_xor = DnaOper(b_dna_add, b_chaos_encode, 'xor');

% DNA固定解码

enImg(:,:,1) = dna_decode(r_dna_xor, r_ca_decode, h, w, 'static');
enImg(:,:,2) = dna_decode(r_dna_xor, g_ca_decode, h, w, 'static');
enImg(:,:,3) = dna_decode(r_dna_xor, b_ca_decode, h, w, 'static');
save('key.mat', 'r_ca_decode', 'g_ca_decode', 'b_ca_decode', 'x0',...
    'u', 'N', 'y', 'r_lgs_flags', 'g_lgs_flags', 'b_lgs_flags',...
    'r_flags', 'g_flags', 'b_flags');
imwrite(uint8(enImg), 'res.bmp');
