clear;

S = 50; % Valeur du sous-jacent
K = 50; % Strike
r = 0.05; % Taux d'intérêt sans risque
sigma = 0.25; % Volatilié
T = 3; % Maturité

M = 100; % Nombre de points de maillage d'actifs
N = 400; % Nombre de points de maillage en temps
Szero = 0; % Valeur extreme
Smax = 150; % Valeur extreme
omega=1.2;tol=0.001; % Paramètres SOR


solution_mesh = zeros(N+1,M+1); % Tableau de la solution approchée
Smesh=0 : (Smax/M) : Smax; % maillage en S
Tmesh=T:-(T/N):0; % maillage en temps t
dt=T/N; % pas de temps, dt
dS=Smax/M; % pas d'actifs, ds
solution_mesh(1,:)=max(K-Smesh,0); % Payoff (Put)
solution_mesh(:,1) = K*exp(-r*(T-Tmesh)); % Condition à la limite en S=0
solution_mesh(:,M+1)=0; % Condition à la limite en S=M
a = @(i) 0.5*dt*(r*i-sigma^2*i^2);
b = @(i) 1 + (sigma^2*i^2 + r)*dt; % Definition des fonctions A, B & C
c = @ (i) -0.5*dt*(sigma^2*i^2+r*i);


Acoeffs = zeros (M+1,1);Bcoeffs = zeros(M+1,1);
Ccoeffs=zeros(M+1,1);
for i=1: M+1
  Acoeffs(i) = a(i-1);
  Bcoeffs(i) = b(i-1);
  Ccoeffs(i) = c(i-1);
end

Tri=diag(Acoeffs(3:M),-1)+diag(Bcoeffs(2:M))+diag(Ccoeffs(2:M-1),+1);
Tri_Inv=inv(Tri); % Calcul de la matrice inverse


% Boucle en temps et Schéma implicite + Itération SOR

for p=1:N
  temp=zeros(M-1,1);
  temp(1)=a(0)*solution_mesh(p+1,1) ;
  temp(end)=c(M)*solution_mesh(p+1,M+1); % Terme à la frontière du domaine
  RHS=solution_mesh(p,2:M)'-temp;
  
% Méthode SOR & Gauss-Siedel :Résolution itèrative de Ax=b
% Noter que : A est (M-1) x (M-1), b est (M-1)

  A=Tri;b=RHS; % Définir la matrice A et le second membre RHS=b
  x=solution_mesh(p,2:M)'; % Définir x
  xold=10000*x; % Initialiser xold pour commencer la boucle d'itération
  n=length(x);
  
  while norm(xold-x)>tol
    xold=x;  % redéfinir xold
    for i=1:n % for: parcourir les lignes de la matrice
      if i==1
        z=(b(i)-A(i,i+1)*x(i+1))/A(i,i);
        x(i) = max(omega*z + (1-omega)*xold(i),K-i*dS);
      elseif i==n
        z=(b(i)-A(i,i-1)*x(i-1))/A(i,i);
        x(i) = max(omega*z + (1-omega)*xold(i),K-i*dS);
      else
        z = (b(i)-A(i,i-1)*x(i-1)-A(i,i+1)*x(i+1))/A(i,i);
        x(i) = max(omega*z + (1-omega)*xold(i),K-i*dS);
      end
    end
  end
  solution_mesh(p+1,(2:end-1))=x;
end

% Sortie graphique
figure (1)
mesh(Smesh,Tmesh,solution_mesh)
xlabel('S');
ylabel('t');
zlabel('V(t,S)','rotation',0)

figure (2)
plot(Tmesh, solution_mesh(:,1),Tmesh,solution_mesh(:,11),Tmesh,solution_mesh(:,41),
Tmesh,solution_mesh(:,101))
xlabel('Temps'); ylabel('Sous-jacent (S)');
legend('0','10','40','100')
