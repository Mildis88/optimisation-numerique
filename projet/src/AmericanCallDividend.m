clear;

% Parametres financiers
S = 50; % Valeur du sous-jacent
K = 50; % Strike
r = 0.05; % Taux d'interet sans risque
sigma = 0.2; % Volatilite
T = 3; % Maturite
delta = 0.02; % Taux de Dividende

% Parametres numeriques
M = 80; % Nombre de points de maillage d'actifs
N = 200; % Nombre de points de maillage en temps
Szero = 0; % Valeur extreme
Smax = 150; % Valeur extreme
omega=1.2;tol=0.001; % Parametres SOR

% Initialisation du maillage et de la matrice du systeme lineaire
solution_mesh = zeros(N+1,M+1); % Tableau de la solution approchee
Smesh=0 : (Smax/M) : Smax; % maillage en S
Tmesh=T:-(T/N):0; % maillage en temps t
dt=T/N; % pas de temps, dt
dS=Smax/M; % pas d'actifs, ds
solution_mesh(1,:)=max(Smesh - K,0); % Payoff (Call)
solution_mesh(:,1) = 0; % Condition a la limite en S=0
solution_mesh(:,M+1)= Smax; % Condition a la limite en S=M
a = @(i) 0.5*dt*((r-delta)*i-0.5*sigma^2*i^2);
b = @(i) 1 + (0.5*sigma^2*i^2 + r)*dt; % Definition des fonctions A, B & C
c = @ (i) -0.5*dt*((r-delta)*i + 0.5*sigma^2*i^2);


% Construction de la matrice tri-diagonale et son inverse
Acoeffs = zeros (M+1,1);Bcoeffs = zeros(M+1,1);
Ccoeffs=zeros(M+1,1);
for i=1: M+1
  Acoeffs(i) = a(i-1);
  Bcoeffs(i) = b(i-1);
  Ccoeffs(i) = c(i-1);
end

Tri=diag(Acoeffs(3:M),-1)+diag(Bcoeffs(2:M))+diag(Ccoeffs(2:M-1),+1);
Tri_Inv=inv(Tri); % Calcul de la matrice inverse

% Boucle en temps et Schema implicite + Iteration SOR
for p=1:N
  temp=zeros(M-1,1);
  temp(1)=a(0)*solution_mesh(p+1,1) ;
  temp(end)=c(M)*solution_mesh(p+1,M+1); % Terme a la frontiere du domaine
  RHS=solution_mesh(p,2:M)'-temp;
  
  % Methode SOR & Gauss-Siedel :Resolution iterative de Ax=b
  % Noter que : A est (M-1) x (M-1), b est (M-1)

  A=Tri;b=RHS; % Definir la matrice A et le second membre RHS=b
  x=solution_mesh(p,2:M)'; % Definir x
  xold=10000*x; % Initialiser xold pour commencer la boucle d'iteration
  n=length(x);
  
  while norm(xold-x)>tol
    xold=x;  % redefinir xold
    for i=1:n % for: parcourir les lignes de la matrice
      if i==1
        z=(b(i)-A(i,i+1)*x(i+1))/A(i,i);
        x(i) = max(omega*z + (1-omega)*xold(i),i-K*dS);
      elseif i==n
        z=(b(i)-A(i,i-1)*x(i-1))/A(i,i);
        x(i) = max(omega*z + (1-omega)*xold(i),i-K*dS);
      else
        z = (b(i)-A(i,i-1)*x(i-1)-A(i,i+1)*x(i+1))/A(i,i);
        x(i) = max(omega*z + (1-omega)*xold(i),i-K*dS);
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
plot(Tmesh, solution_mesh(:,1),Tmesh,solution_mesh(:,11),Tmesh,
solution_mesh(:,41),Tmesh,solution_mesh(:,81))
xlabel('Temps'); ylabel('Sous-jacent (S)');
legend('0','10','40','100')

