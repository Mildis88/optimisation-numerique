clear

# Paramètres physiques :

K = 2*10^-6; % Coefficient de diffusion K est 10e-6 m^2/s

# Paramètres numériques
dz = .25; % Chaque pas de la profondeur est 1/4 mètre
Nz = 400; % Choisir le pas de profondeur (devrait aller à au moins 100 m)
Nt = 5000; % Choisir le pas de temps
dt = (365*24*60*60)/Nt; % L'unité de pas de temps est en secondes
                        % (~ 6.3*10^3 secondes, or ~105 minutes)
                        
                        
ones(Nz+1, Nt+1);

% Créer une matrice de température avec Nz + 1 lignes,
% et Nt + 1 colonnes
% Donnée initiale est que u est de 15 partout.
# Initialisation du maillage et de la donnée initiale

temps = [0:12/Nt:12];
u = 15*ones(Nz+1,Nt+1);

# Condition à la surface:

u(1,:) = 15-10*sin(2*pi*temps/12); %Imposer la température de surface

# Boucle en temps :

tic

maxiter = 500;
for iter = 1:maxiter
  uold = u; % Sauver la dernière température calculée
  u(:,1) = uold (:,end); % Initialiser la température à t=0 au dernier
                           % temps évalué
  for i=2:Nt+1,
    profondeur = (u(1:end-2,i-1)-2*u(2:end-1,i-1)+u(3:end,i-1))/dz^2;
    temps_1D = K*profondeur;
    u(2:end-1,i) = temps_1D*dt + u(2:end-1,i-1) ;
    u(end,i) = u(end-1,i); % Condition à la limite au fond
  end
  err(iter) = max(abs(u(:)-uold(:))); % Trouver le maximum en valeur
                                      % absolue entre deux solutions
  
  if err(iter)<1E-4
    break; % Stop si les solutions sont semblables,
           % on a convergence
  end
end

# Sortie graphique

figure(1)
plot(log(err)), title ('Graphe de Convergence')

figure(2)
imagesc([0 12], [0 100],u); title ('Variation de Temperature (imagesc)')
colorbar

figure(3)
profondeur = [0:dz:Nz*dz];
contourf(temps,-profondeur,u) ; title('Variation de Temperature (contourf)')
colorbar

figure (4)
plot(temps,u(1,:),temps,u(21,:),temps,u(41,:),temps,u(61,:),temps,u(81,:))
xlabel ('Temps (mois) '); ylabel ('Temperature (C)');
legend ('Om','5m','10m','15m','20m')

toc