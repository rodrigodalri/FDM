%-------------------------------------------------------------------------%
%  This simple program computes the Electric Fields due to 
%  Parallel plate Capacitors using the Finite difference method (FDM)  
%-------------------------------------------------------------------------%

clc
close all; clear all;

%-------------------------------------------------------------------------%
%                   SYMBOLS USED IN THIS CODE                             
%-------------------------------------------------------------------------%

% E = Total electric field matrix using Poisson's equation
% V = Potential matrix
% Nx = Number of grid points in X- direction
% Ny = Number of grid points in Y-Direction
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
%                         INITIALIZATION                                  
%          Here, all the grid, size, charges, etc. are defined
%-------------------------------------------------------------------------%

% Enter the dimensions
coefAjuste= 20;
Nx = 3*coefAjuste;     % Number of X-grids
Ny = 20*coefAjuste;     % Number of Y-grids
mpx = ceil(Nx/2); % Mid-point of x
mpy = ceil(Ny/2); % Mid point of y
   

Ni = 1000;  % Number of iterations for the Poisson solver

V = zeros(Nx,Ny);   % Potential (Voltage) matrix

T = 0;            % Top-wall potential
B = 0;            % Bottom-wall potential
L = 0;            % Left-wall potential
R = 0;            % Right-wall potential

%-------------------------------------------------------------------------%
% Initializing edges potentials
%-------------------------------------------------------------------------%

V(1,:) = L;
V(Nx,:) = R;
V(:,1) = B;
V(:,Ny) = T;

%-------------------------------------------------------------------------%
% Initializing Corner potentials
%-------------------------------------------------------------------------%

V(1,1) = 0.5*(V(1,2)+V(2,1));
V(Nx,1) = 0.5*(V(Nx-1,1)+V(Nx,2));
V(1,Ny) = 0.5*(V(1,Ny-1)+V(2,Ny));
V(Nx,Ny) = 0.5*(V(Nx,Ny-1)+V(Nx-1,Ny));

%-------------------------------------------------------------------------%


length_plate = 10*coefAjuste;  % Length of plate in terms of number of grids  
lp = floor(length_plate/2);

position_plate = 0.05*coefAjuste; % Position of plate on x axis
pp1 = mpx+position_plate;
pp2 = mpx-position_plate;

for z = 1:Ni    % Number of iterations
    clc;
    fprintf('%0.2f',(z*100/Ni));
    disp('%');
        for i=2:Nx-1            
            for j=2:Ny-1      

                % The next two lines are meant to force the matrix to hold the 
                % potential values for all iterations


                    V(i,j)=0.25*(V(i+1,j)+V(i-1,j)+V(i,j+1)+V(i,j-1));
                    V(pp1,mpy-lp:mpy+lp) = 10;
                    V(pp2,mpy-lp:mpy+lp) = 0;
            end
        end
        
end

% Take transpose for proper x-y orientation
V = V';

[Ex,Ey]=gradient(V);
Ex = -Ex;
Ey = -Ey;



% Electric field Magnitude
 E = (sqrt(Ex.^2+Ey.^2)*coefAjuste*1000); %%convert from v/mm ti v/m  
 a = E(mpy, mpx);

x = ((1:Nx)-mpx)/(coefAjuste);
y = ((1:Ny)-mpy)/(coefAjuste);

% Contour Display for electric potential
figure(1);
contour_range_V = -101:0.5:101;
                
contour(x,y,V,contour_range_V,'linewidth',0.05);
line([(pp1-mpx)/coefAjuste,(pp1-mpx)/coefAjuste],[(-lp)/coefAjuste,(+lp)/coefAjuste],'linewidth', 2, 'color', 'k');
line([(pp2-mpx)/coefAjuste,(pp2-mpx)/coefAjuste],[(-lp)/coefAjuste,(+lp)/coefAjuste],'linewidth', 2, 'color', 'k');
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14);
xlabel('X','fontsize',14);
ylabel('Y','fontsize',14);
title('Distribuiçao potencial, V(x,y) em Volts','fontsize',14);
h1=gca;
set(h1,'fontsize',14);
fh1 = figure(1); 
set(fh1, 'color', 'white')


% Contour Display for electric field
figure(2);
contour_range_E = -100001:25000:100001;
contour(x,y,E,contour_range_E,'linewidth',0.05, 'showtext', 'on');
line([(pp1-mpx)/coefAjuste,(pp1-mpx)/coefAjuste],[(-lp)/coefAjuste,(+lp)/coefAjuste],'linewidth', 2, 'color', 'k');
line([(pp2-mpx)/coefAjuste,(pp2-mpx)/coefAjuste],[(-lp)/coefAjuste,(+lp)/coefAjuste],'linewidth', 2, 'color', 'k');
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14);
xlabel('X','fontsize',14);
ylabel('Y','fontsize',14);
title('Ditribuição do campo elétrico, E (x,y) V/m','fontsize',14);
h2=gca;
set(h2,'fontsize',14);
fh2 = figure(2); 
set(fh2, 'color', 'white')


for i=2:Ny-1
    if(mod(i,7)~=0)
        Ex(i,:)=0;
        Ey(i,:)=0;        
    end
end

% Quiver Display for electric field Lines
figure(3);
%contour(x,y,E,'linewidth',0.05);
hq = quiver(x,y,Ex,Ey,2);
title('Linhas do campo Elétrico, E (x,y) V/m','fontsize',14);
line([(pp1-mpx)/coefAjuste,(pp1-mpx)/coefAjuste],[(-lp)/coefAjuste,(+lp)/coefAjuste],'linewidth', 2, 'color', 'k');
line([(pp2-mpx)/coefAjuste,(pp2-mpx)/coefAjuste],[(-lp)/coefAjuste,(+lp)/coefAjuste],'linewidth', 2, 'color', 'k');
hold on;
txt1 = strcat('\leftarrow E[0,0] = ', num2str(E(mpy, mpx)));
text(0,0,txt1);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14);
xlabel('X','fontsize',14);
ylabel('Y','fontsize',14);
h3=gca;
set(h3,'fontsize',14);
fh3 = figure(3); 
set(fh3, 'color', 'white');