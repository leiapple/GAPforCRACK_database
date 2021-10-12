% define some parameters
npoints=4000;  % total new points generated
ndim=9;        % dimention of new points
a_0=2.834;     % equilibrium lattice constant
I=eye(3);      % unit tensor
delta=0.2;     % overapped distance criterion 
%======================================================
% First, calculate the strain state from the database DB1.xyz.
% read data from DB1
delimiterIn   = ' '; 
headerlinesIn = 0;  
filename1='db1_shrinked.dat';
CC = importdata(filename1, delimiterIn, headerlinesIn);
% format of D: id_anumber, volume, energy, a1, a2, a3, b1, b2, b3, c1, c2, c3).

% undeformed primitive lattice constant
config=size(CC,1);
eps_all_db1=[0 0 0 0 0 0];

filename2 = fopen('db1_strain','w');
for i=1:config
    % deformed lattice vector
    a=CC(i,4:6);
    b=CC(i,7:9);
    c=CC(i,10:12);
    % calculate deformation tensor: F
    %{
    F=[a(1)/a_0,b(1)/a_0,(2*c(1)-a(1)-b(1))/a_0;
       a(2)/a_0,b(2)/a_0,(2*c(2)-a(2)-b(2))/a_0;
       a(3)/a_0,b(3)/a_0,(2*c(3)-a(3)-b(3))/a_0];
    %}
    F=[(a(1)-a(3))/a_0,(a(2)-a(3))/a_0,(2*a(3))/a_0;
       (b(1)-b(3))/a_0,(b(2)-b(3))/a_0,(2*b(3))/a_0;
       (c(1)-c(3))/a_0,(c(2)-c(3))/a_0,(2*c(3))/a_0] ;  
    epsilon = ((transpose(F)*F)-I)/2;
    eps_all_db1(end+1,:)=([epsilon(1,1),epsilon(1,2),epsilon(1,3),epsilon(2,2),epsilon(2,3),epsilon(3,3)]);
%    fprintf(filename2,'%.8f %.8f %.8f %.8f %.8f %.8f\n',epsilon(1,1),epsilon(1,2),epsilon(1,3),epsilon(2,2),epsilon(2,3),epsilon(3,3));
end

% Addtional operation can be added here to elimite the too closed data in
% original DB1.xyz

%======================================================

%======================================================
% Load the sample data set. 
p = sobolset(ndim,'Skip',1e4,'Leap',1e3);
% Apply a random linear scramble combined with a random digital 
% shift by using scramble.
p = scramble(p,'MatousekAffineOwen');
% Generate 1000 points by using net.
X0 = net(p,npoints);
% linear rescale of the random array 
DD(:,1:8) = X0(:,1:8)*0.6-0.07;
DD(:,9)= X0(:,9)*0.6-0.35;
% Calculate the strain
eps_all=[0 0 0 0 0 0];
% open file to write all the lattice input needed for DFT calculation
filename2 = fopen('Lattice_data.txt','w');
for i=1:npoints
    % deformed lattice vector
    a=DD(i,1:3);
    b=DD(i,4:6);
    c=DD(i,7:9);
    a(1)=a(1)+a_0;
    b(2)=b(2)+a_0;
    c(:)=c(:)+a_0/2 ;
    % calculate deformation tensor: F
    %{
    F=[a(1)/a_0,b(1)/a_0,(2*c(1)-a(1)-b(1))/a_0;
       a(2)/a_0,b(2)/a_0,(2*c(2)-a(2)-b(2))/a_0;
       a(3)/a_0,b(3)/a_0,(2*c(3)-a(3)-b(3))/a_0];
    %}
    F=[(a(1)-a(3))/a_0,(a(2)-a(3))/a_0,(2*a(3))/a_0;
       (b(1)-b(3))/a_0,(b(2)-b(3))/a_0,(2*b(3))/a_0;
       (c(1)-c(3))/a_0,(c(2)-c(3))/a_0,(2*c(3))/a_0] ;  
    epsilon = ((transpose(F)*F)-I)/2;
    %  delete the points that are too close to each other: the criterion is
    %  set to be "norm < 0.2".
    for j=1:config
        disp(j)=norm(epsilon-eps_all_db1(j));
    end
    if min(disp) > delta
        eps_all(end+1,:)=([epsilon(1,1),epsilon(1,2),epsilon(1,3),epsilon(2,2),epsilon(2,3),epsilon(3,3)]);
        fprintf(filename2,'%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n',a(:),b(:),c(:));
    end
end

[coeff1,score1,latent1,tsquared1] = pca(eps_all,'NumComponents',2);
[coeff2,score2,latent2,tsquared2] = pca(eps_all_db1,'NumComponents',2);

figure(1); 
scatter(score1(:,1),score1(:,2),5,'blue','filled','o')
hold on;
scatter(score2(:,1),score2(:,2),5,'red','filled','s')
hold on;