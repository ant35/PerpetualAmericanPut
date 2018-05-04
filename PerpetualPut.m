clc
format long g
clear
%% Parameters
tol = 0.01;
r = 0.05; %risk free rate
K =100;%strike
N = 40; % size
intrinsic_value = @(K,s) max(K-s,0);
order = N*N;% setting
theta = 0.5;
kappa = 0.5;
ksi = 0.5;
rho = 0;
%% Define coordinate vector
coordinates = [
    1,0; %1
    -1,0; %2
    0,-1; %3
    0,1; %4
    0,0; %5
    1,1; %6
    1,-1; %7
    -1,1; %8
    -1,-1; %9
    ];
c = zeros(length(coordinates),1);
%% Initialize Spot vs Volatility Grid
LargeA = zeros(N);
A = LargeA(2:N-1,2:N-1); %Take off the boundaries
dX = 3;
dY = .1;
%% Set mesh for spot on top row and vol on left column
SpotMesh = 0:dX:(N-1)*dX;
VolMesh = 0:dY:(N-1)*dY;
LargeA(1,1:N) = SpotMesh; %set spot coordinates
LargeA(1:N,1) = VolMesh; %set vol coordinates

%% Set Boundary minus corners for A using I
    for j = 1:N %top row, Left side
        LargeA(1,j)=GBM_Price(K,SpotMesh(j),VolMesh(1),r,GBM_Lstar(K,VolMesh(1),r));%I(mapping(j,1,N)); 
    end
    %Left side
    for j = 2:N
        LargeA(j,1)=GBM_Price(K,SpotMesh(1),VolMesh(j),r,GBM_Lstar(K,VolMesh(j),r));%I(mapping(1,j,N));
    end
    %bottom row
    for j = 2:N
        LargeA(N,j)=GBM_Price(K,SpotMesh(j),VolMesh(N),r,GBM_Lstar(K,VolMesh(N),r));%I(mapping(j,N,N));
    end
    for j = 2:N-1 %bottom row, right side
        LargeA(j,N)=GBM_Price(K,SpotMesh(N),VolMesh(j),r,GBM_Lstar(K,VolMesh(j),r));%I(mapping(N,j,N));
    end
       
    %% Setting f
    dt=1;
    dS = dX; dNu = dY;
  
    c1 = @(s,nu) s/(2*dS)*(nu*s/dS +r);
    c2 = @(s,nu) s/(2*dS)*(nu*s/dS - r);
    c3 = @(s,nu)(0.5/dNu)*(ksi^2*nu/dNu - kappa*(theta - nu));
    c4 = @(s,nu)(0.5/dNu)*(ksi^2*nu/dNu + kappa*(theta -nu));
    c5 = @(s,nu) -nu*(s^2/(dS^2)+ksi^2/(dNu^2))-r;
    c6 = @(s,nu) rho*ksi*nu*s/(4*dS*dNu);
    c7 = @(s,nu) -rho*ksi*nu*s/(4*dS*dNu);
    c8 = @(s,nu) -rho*ksi*nu*s/(4*dS*dNu);
    c9 = @(s,nu)  rho*ksi*nu*s/(4*dS*dNu);    
    
    f = zeros((N-2)^2,1);
    n = N-2;
    %Top left corner
    spot = SpotMesh(2); vol = VolMesh(2);
    f(1) = c8(spot,vol)*LargeA(1,1)+c4(spot,vol)*LargeA(1,2) ...
        +c6(spot,vol)*LargeA(1,3)+c2(spot,vol)*LargeA(2,1) ...
        +c9(spot,vol)*LargeA(3,1);
    
    %Top inner row
    for i =2:n-1
        spot = SpotMesh(i); vol = VolMesh(2);
        f(i)=c8(spot,vol)*LargeA(1,i-1)+c4(spot,vol)*LargeA(1,i)...
            +c6(spot,vol)*LargeA(1,i+1);
    end
    
    %Right top corner
    spot = SpotMesh(N-1); vol=VolMesh(2);
    f(N-2)= c6(spot,vol)*LargeA(1,N)...
        +c4(spot,vol)*LargeA(1,N-1)+c8(spot,vol)*LargeA(1,N-2)...
        +c1(spot,vol)*LargeA(2,N)+c7(spot,vol)*LargeA(3,N);
    
    %Left inner column
    for i = 2:n-1
        spot = SpotMesh(2); vol = VolMesh(i);
        f((i-1)*n+1) = c8(spot,vol)*LargeA(i-1,1)+c2(spot,vol)*LargeA(i,1)...
            +c9(spot,vol)*LargeA(i+1,1);
    end
    
    %Left bottom corner
    spot = SpotMesh(2); vol = VolMesh(N-1);
    f(n*(n-1)+1)= c9(spot,vol)*LargeA(N,1)+...
        c3(spot,vol)*LargeA(N,2)+c7(spot,vol)*LargeA(N,3)...
        +c2(spot,vol)*LargeA(N-1,1)+c8(spot,vol)*LargeA(N-2,1);
    
    %Right inner column
    for i = 2:n-1
        spot=SpotMesh(N-1);vol=VolMesh(i);
        f(i*n) = c6(spot,vol)*LargeA(i-1,N)+c1(spot,vol)*LargeA(i,N)...
            +c7(spot,vol)*LargeA(i+1,N);
    end
    
    %Bottom right corner
    spot = SpotMesh(N-1); vol = VolMesh(N-1);
    f(n^2) = c7(spot,vol)*LargeA(N,N)+c3(spot,vol)*LargeA(N,N-1)...
        +c9(spot,vol)*LargeA(N,N-2)+c1(spot,vol)*LargeA(N-1,N)+...
        c6(spot,vol)*LargeA(N-2,N);
    
    %Bottom inner row
    for i=2:n-1
        spot = SpotMesh(i); vol = VolMesh(N-1);
        f((N-3)*(N-2)+i-1) = c9(spot,vol)*LargeA(N,i-1)+c3(spot,vol)*LargeA(N,i)...
            +c7(spot,vol)*LargeA(N,i+1);
    end
    
    els = 1:(N-2)^2;
    els = els';
    
    %[f,els]
    %% Setting U matrix
    U = zeros(n^2);
    vhat = zeros((N-2)^2,1); %Initialize guess for put value
    f = -f;
    iv = zeros(n,1);
    row = 1;
    
    
    for k=1:n
        nu = VolMesh(1,k+1);
        for i = 1:n
            s = SpotMesh(1,i+1);
            iv(row) = intrinsic_value(K,s);
            
            %Reset stencil
            c(1) = s/(2*dS)*(nu*s/dS +r);
            c(2) = s/(2*dS)*(nu*s/dS - r);
            c(3) = (0.5/dNu)*(ksi^2*nu/dNu - kappa*(theta - nu));
            c(4) = (0.5/dNu)*(ksi^2*nu/dNu + kappa*(theta -nu));
            c(5) = -nu*(s^2/(dS^2)+ksi^2/(dNu^2))-r;
            c(6) = rho*ksi*nu*s/(4*dS*dNu);
            c(7) = -rho*ksi*nu*s/(4*dS*dNu);
            c(8) = -rho*ksi*nu*s/(4*dS*dNu);
            c(9) =  rho*ksi*nu*s/(4*dS*dNu);
            
            %Map stencil to V plane
            if k == 1 && i ==1 %Top Left corner
                U(row,1)=c(5);
                U(row,2)=c(1);
                U(row,n+1)=c(3);
                U(row,n+2)=c(7);
            elseif k==1 && i > 1 && i < N-2 %Top inner row
                U(row,i-1)=c(2);
                U(row,i)=c(5);
                U(row,i+1)=c(1);
                U(row,n+i-1)=c(9);
                U(row,n+i)=c(3);
                U(row,n+i+1)=c(7);

            elseif k==1 && i == N-2 %Top right corner
                U(row,i-1)=c(2);
                U(row,i) = c(5);
                U(row,2*i-1)=c(9);
                U(row,2*i)=c(3);

            elseif k > 1 && k < N-2 && i == 1 %inner Left column
                U(row,(k-2)*n+1)=c(4);
                U(row,(k-2)*n+2)=c(6);
                U(row,(k-1)*n+1)=c(5);
                U(row,(k-1)*n+2)=c(1);
                U(row,k*n+1)=c(3);
                U(row,k*n+2)=c(7);

            elseif k == N-2 && i == 1 %Bottom Left corner
                U(row,(k-2)*n + 1)=c(4);
                U(row,(k-2)*n + 2)=c(6);
                U(row,(k-1)*n+1)=c(5);
                U(row,(k-1)*n+2)=c(1);

            elseif k > 1 && k < N-2 && i == N-2 %Inner Right column
                U(row,(k-1)*n-1) = c(8);
                U(row,(k-1)*n) = c(4);
                U(row,(k)*n-1) = c(2);
                U(row,(k)*n) = c(5);
                U(row,(k+1)*n-1) = c(9);
                U(row,(k+1)*n) = c(3);
            
            elseif k == N-2 && i > 1 && i < N-2 %Bottom inner row
                U(row,n*(n-2)+i-1)=c(8);
                U(row,n*(n-2)+i)=c(4);
                U(row,n*(n-2)+i+1)=c(6);
                U(row,n*(n-1)+i-1)=c(2);
                U(row,n*(n-1)+i)=c(5);
                U(row,n*(n-1)+i+1)=c(1);
            
            elseif k == N-2 && i == N-2 %Bottom right corner
                U(row,n*(n-1)-1)=c(8);
                U(row,n*(n-1))=c(4);
                U(row,n*n-1) = c(2);
                U(row,n*n)=c(5);
               
            elseif k > 1 && k < n && i > 1 && i < n %interior point
               U(row,(k-2)*n + i-1)=c(8);
               U(row,(k-2)*n + i)=c(4);
               U(row,(k-2)*n + i+1)=c(6);
               U(row,(k-1)*n + i-1)=c(2);
               U(row,(k-1)*n + i)= c(5);
               U(row,(k-1)*n + i+1)=c(1);
               U(row,k*n + i-1)=c(9);
               U(row,k*n + i)=c(3);
               U(row,k*n + i+1)=c(7);
            else
                error(["Unknown point in f: k=" k " i=" i]);
            end
            row = row + 1;
        end
    end

%% Out put results
z=U\f;
z = max(z,iv); %optimal boundary
row = 1;

for i = 1:N-2
    for j=1:N-2
        A(i,j) = z(row);%intrinsic_value(K,K-O(row));
        row = row+1;
    end
end

LargeA(2:N-1,2:N-1) = A; %Implant A in LargeA

copy(2:N-1,2:N-1) = A; %Implant A in LargeA
surf(SpotMesh(1,:),VolMesh(1,:),LargeA)
xlabel("Spot Price")
ylabel("\nu")
zlabel("Perpetual American Put Value")
title(["Explicit Method, Strike " K])
figure

surf(SpotMesh(1,:),VolMesh(1,:),K-LargeA)
xlabel("Spot Price")
ylabel("\nu")
zlabel("Optimal L value")
title(["Explicit Method, Strike " K])


%% Helper Functions
function mesh = nonUniformMeshY(N,V,d)
    dEta = (1/N)*asinh(V/d);
    Eta = zeros(N,2);
    for j = 1:N
        Eta(j,1) = j*dEta;
    end
    mesh = zeros(N,2);
    mesh(1,2) = 1;
    for j = 2:N
        mesh(j,1)=d*sinh(Eta(j));
        mesh(j,2) = mesh(j,1) - mesh(j-1,1);
    end
end
function mesh = nonUniformMeshX(N,K,c,S_upper)
    ksi = zeros(N,1);
    dKsi = (1/N)*(asinh((S_upper-K)/c) - asinh(-K/c));
    for m = 1:N-1
        ksi(m+1,1) = asinh(-K/c) + m*dKsi;
    end

    mesh = zeros(N,2);
    mesh(1,2) = 1;
    for i = 2:N
        mesh(i,1)= K + c*sinh(ksi(i));
        mesh(i,2) = mesh(i,1) - mesh(i-1,1);
    end
end
