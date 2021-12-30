#initialisation paramètres financiers
S = 50; #sous jacent
K = 50; #strike
r = 0.05;
sigma = 0.2;
T = 3;
#initialisation paramètres numériques
M = 150;
N = 200;
Szero = 0;
Smax = 150;

#initialisation maillage et matrice du système linéaire
solution_mesh = zeros(N+1,M+1);
Smesh = 0:(Smax/M):Smax;
Tmesh = T:-(T/N):0;
dt = T/N;
solution_mesh(1,:) = max(K-Smesh,0);
solution_mesh(:,1) = K*exp(-r*(T-Tmesh));
solution_mesh(:,M+1) = 0;
A = @(i) 0.5*dt*(r*i-sigma^2*i^2);
B = @(i) 1 + (sigma^2*i^2 +r)*dt;
C = @(i) - 0.5*dt*(sigma^2*i^2+r*i);

#Matrice tridiagonale
Acoeffs = zeros(M+1,1);
Bcoeffs = zeros(M+1,1);
Ccoeffs = zeros(M+1,1);
for i=1:M+1 Acoeffs(i)=A(i-1);
  Bcoeffs(i)=B(i-1);
  Ccoeffs(i)=C(i-1);
end
Tri=diag(Acoeffs(2:end),-1) + diag(Bcoeffs) + diag(Ccoeffs(1:end:-1),+1);
Tri_Inv = inv(Tri);

#boucle temps et schéma implicite
for j=1:N temp=zeros(M+1,1);
  temp(1)=A(0)*solution_mesh(j+1,1);
  temp(end)=C(M)*solution_mesh(j+1,M+1);
  RHS=solution_mesh(j,:)'-temp;
  temp=Tri_Inv*RHS;
  solution_mesh(j+1,(2:end-1))=temp(2:end-1);
end

#sortie graphique
figure(1)
mesh(Smesh,Tmesh,solution_mesh)
xlabel('S');ylabel('t'); zlabel('V(t,S)','rot',0)

figure(2)
plot(Tmesh,solution_mesh(:,1),Tmesh,solution_mesh(:,10),Tmesh,solution_mesh(:,50),Tmesh)
xlabel('Temps');ylabel('Sous-jacent (S)'); legend('0','50','100','150')