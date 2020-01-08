%fid = fopen('BZpath.dat', 'w');

path = pwd;
addpath(path);
fname = 'selfjobs';
mkdir(fname);
cd(fname)

steps = 20;
counter = 0;

qxstart = 0;
qystart = 0;
qxend = 4*pi/3;
qyend = 0;
dqx = (qxend-qxstart)/steps;
dqy = (qyend-qystart)/steps;

figure

% the path from \Gamma -> K  
% (0, 0) -> (4*pi/3, 0)
for j = 1:steps
    counter = counter + 1;
    fjob = ['job', num2str(counter)];
    mkdir(fjob);
    cd(fjob);
    fid = fopen('input.txt', 'w');
    vec1 = qxstart + j*dqx;
    vec2 = qystart + j*dqy;
    [q1, q2] = qxytoq12(vec1, vec2);
    plot(vec1, vec2, 'r*');
    hold on;
    fprintf(fid, '%2.0f \n', counter-1);
    fprintf(fid, '%20.16f \n', q1);
    fprintf(fid, '%20.16f \n', q2);
    fclose(fid);
    cd ..
end

% the path from K -> M  
% (4*pi/3, 0) -> (2*pi/sqrt(3), 2*pi/3)
qxstart = 4*pi/3;
qystart = 0;
qxend = pi;
qyend = pi/sqrt(3);

steps = 10;
dqx = (qxend - qxstart)/steps;
dqy = (qyend - qystart)/steps;

for j = 1:steps
    counter = counter + 1;
    fjob = ['job', num2str(counter)];
    mkdir(fjob);
    cd(fjob);
    fid = fopen('input.txt', 'w');
    vec1 = qxstart + dqx*j;
    vec2 = qystart + dqy*j;
    [q1, q2] = qxytoq12(vec1, vec2);
    plot(vec1, vec2, 'g*');
    hold on;
    fprintf(fid, '%2.0f \n', counter-1);
    fprintf(fid, '%20.16f \n', q1);
    fprintf(fid, '%20.16f \n', q2);
    fclose(fid);
    cd ..
end


% the path from M -> \Gamma
qxstart = pi;
qystart = pi/sqrt(3);
qxend = 0;
qyend = 0;

steps = 17;
dqx = (qxend - qxstart)/steps;
dqy = (qyend - qystart)/steps;

for j = 1:steps
    counter = counter + 1;
    fjob = ['job', num2str(counter)];
    mkdir(fjob);
    cd(fjob);
    fid = fopen('input.txt', 'w');
    vec1 = qxstart + dqx*j;
    vec2 = qystart + dqy*j;
    [q1, q2] = qxytoq12(vec1, vec2);
    plot(vec1, vec2, 'b*');
    hold on;
    fprintf(fid, '%2.0f \n', counter-1);
    fprintf(fid, '%20.16f \n', q1);
    fprintf(fid, '%20.16f \n', q2);
    fclose(fid);
    cd ..
end

% the path from Y -> M

qxstart = 0;
qystart = pi/sqrt(3);
qxend = pi;
qyend = pi/sqrt(3);

steps = 15;
dqx = (qxend - qxstart)/steps;
dqy = (qyend - qystart)/steps;

for j = 1:steps
    counter = counter + 1;
    fjob = ['job', num2str(counter)];
    mkdir(fjob);
    cd(fjob);
    fid = fopen('input.txt', 'w');
    vec1 = qxstart + dqx*j;
    vec2 = qystart + dqy*j;
    [q1, q2] = qxytoq12(vec1, vec2);
    plot(vec1, vec2, 'k*');
    hold on;
    fprintf(fid, '%2.0f \n', counter-1);
    fprintf(fid, '%20.16f \n', q1);
    fprintf(fid, '%20.16f \n', q2);
    fclose(fid);
    cd ..
end

% the path from M -> X
qxstart = pi;
qystart = pi/sqrt(3);
qxend = pi;
qyend = 0;

steps = 43;
dqx = (qxend - qxstart)/steps;
dqy = (qyend - qystart)/steps;

for j = 1:steps
    counter = counter + 1;
    fjob = ['job', num2str(counter)];
    mkdir(fjob);
    cd(fjob);
    fid = fopen('input.txt', 'w');
    vec1 = qxstart + dqx*j;
    vec2 = qystart + dqy*j;
    [q1, q2] = qxytoq12(vec1, vec2);
    plot(vec1, vec2, 'y*');
    hold on;
    fprintf(fid, '%2.0f \n', counter-1);
    fprintf(fid, '%20.16f \n', q1);
    fprintf(fid, '%20.16f \n', q2);
    fclose(fid);
    cd ..
end

xlabel('k1')
ylabel('k2')
hold off



