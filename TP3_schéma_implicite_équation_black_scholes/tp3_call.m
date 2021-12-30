clear;

%Paramètres financiers
S = 50; %Valeur du sous-jacent
K = 50; %Strike
r = 0.05; %Taux d'interet sans risque
sigma = 0.2; %Volatilité
T = 3; %Maturité

%Paramètres numériques
M = 150; %Nombre de points de maillage d'actifs
N = 200; %Nombre de points de maillage en temps
Szero = 0; %Valeur extreme
Smax = 150; %Valeur extreme

%Initialisation du maillage et de la matrice du système linéaire
solution_mesh = zeros(N+1, M+1); %Tableau de la solution approchée
Smesh = 0 : (Smax/M) : Smax; %Maillage en espace S
Tmesh = T : -(T/N) : 0; %Maillage en temps T
dt = T/N; %Pas de temps, dt
solution_mesh(1,:) = max(Smesh-K, 0); %Payoff (Put)
solution_mesh(:,1) = 0; %Condition à la limite en S = 0
solution_mesh(:,M+1) = 150; %Condition à la limite en S = M
A = @(i) 0.5 * dt * (r * i - sigma^2 * i^2);
B = @(i) 1 + (sigma^2 * i^2 + r) * dt;
C = @(i) -0.5 * dt * (sigma^2 * i^2 + r * i);

%Construction de la matrice tri-diagonale et son inverse
Acoeffs = zeros(M+1, 1);
Bcoeffs = zeros(M+1, 1);
Ccoeffs = zeros(M+1, 1);

for i = 1 : M+1
  Acoeffs(i) = A(i-1);
  Bcoeffs(i) = B(i-1);
  Ccoeffs(i) = C(i-1);
endfor

Tri = diag(Acoeffs(2 : end), -1) + diag(Bcoeffs) + diag(Ccoeffs(1 : end-1), +1);
Tri_Inv = inv(Tri); %Calcul de la matrice inverse

%Boucle en temps et schéma implicite
for j = 1 : N
  temp = zeros(M+1, 1);
  temp(1) = A(0) * solution_mesh(j+1, 1);
  temp(end) = C(M) * solution_mesh(j+1, M+1); %Termes frontiers
  RHS = solution_mesh(j,:)'-temp;
  temp = Tri_Inv * RHS;
  solution_mesh(j+1, (2 : end-1)) = temp(2 : end-1);
endfor

%Sortie graphique
figure(1)
mesh(Smesh, Tmesh, solution_mesh)
xlabel('S'); 

figure(2)
plot(Tmesh, solution_mesh(:,1), Tmesh, solution_mesh(:,10), Tmesh, solution_mesh(:,50),Tmesh, solution_mesh(:,100),Tmesh, solution_mesh(:,150))
xlabel('Temps'); ylabel('Sous-jacent (S)'); legend('0', '50', '100', '150')