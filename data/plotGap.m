clear all

mat = 'TaAs';
ext = '_gap_100.dat';

filename = strcat(mat, ext);
data = load(filename);

ld = 100;
x = 1:ld;
y = x;
z = data(:,3);

min = min(z);

gap = zeros(ld, ld);

for i = 1:ld
    for j = 1:ld
       gap(i, j) = z(j + (i - 1) * ld);  
    end     
end    

[exe, why] = meshgrid(x, y);

surf(exe, why, gap);
shading interp
colorbar
title('TaAs band-gap (100 plane) (with spin-orbit coupling)')