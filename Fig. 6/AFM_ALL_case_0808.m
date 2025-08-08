%% Figure 6
clc; clear all; close all

%Definition of parameters
k = linspace(-4e8,4e8,2e4+1);

d1 = 10e-9;
d2 = 10e-9;
Ms1 = 1.4e5;
Ms2 = 7.4e5;
Aex1 = 3.7e-12;
Aex2 = 8.7e-12;
u0 = 4*pi*1e-7;
gamma = 28e9*2*pi;
JSb = 1e6; %1e6 ~ 35GHz 
JSc = 1e6; %1e6 ~ 35GHz
N=3;

% The overall effective external field
Hbc = 5e5;
H1 = -2.75e6; %Fig. 6b
% H1 = -1e6; %Fig. 6c



N1 = (1-exp(-abs(k)*d1))./(abs(k)*d1);
N2 = (1-exp(-abs(k)*d2))./(abs(k)*d2);


a = 2/1.5e8;%

w1P = (2*Aex1*gamma*k.^2/Ms1 + gamma*u0*H1) + gamma*u0*Ms1/2 +0*k;% the last term accounts for intralayer dipolar interaction (self)
w1N = -w1P;
w2P = gamma*u0*(Hbc + JSc) + gamma*u0*Ms2/2 +0*k;%0*k is for for loop convenience

w2N = -w2P;
w3P = gamma*u0*(-Hbc + JSb) + gamma*u0*Ms2/2 +0*k;

w3N = -w3P;

B2P =gamma*u0*JSc*cos(k*a/2);
% B2P = 0*k;
B3P =gamma*u0*JSb*cos(k*a/2);
% B3P = 0*k;
B2N =-B2P;
B3N =-B3P;


f2 = Ms1*d1*N1.*N2/2;
f1 = Ms2*d2*N1.*N2/2;


theta_1 = pi/2;
theta_2 = pi/2;
theta_3 = -pi/2;

D1P = 1/2*gamma*u0.*f1.*(k*(sin(theta_1)-sin(theta_2))+abs(k)*(1-sin(theta_1)*sin(theta_2))); %intra b to a 

D2P = 1/2*gamma*u0.*f1.*(k*(sin(theta_1)-sin(theta_3))+abs(k)*(1-sin(theta_1)*sin(theta_3))); %intra c to a

D3P = 1/2*gamma*u0.*f2.*(k*(sin(theta_1)-sin(theta_2))+abs(k)*(1-sin(theta_1)*sin(theta_2))); %intra a to b

D4P = 1/2*gamma*u0.*f2.*(k*(sin(theta_1)-sin(theta_3))+abs(k)*(1-sin(theta_1)*sin(theta_3))); %intra a to c



D5P = 1/2*gamma*u0.*f1.*(k*(sin(theta_2)-sin(theta_1))+abs(k)*(1-sin(theta_2)*sin(theta_1))); %inter b to a

D6P = 1/2*gamma*u0.*f1.*(k*(sin(theta_3)-sin(theta_1))+abs(k)*(1-sin(theta_3)*sin(theta_1))); %inter c to a

D7P = 1/2*gamma*u0.*f2.*(k*(sin(theta_2)-sin(theta_1))+abs(k)*(1-sin(theta_2)*sin(theta_1))); %inter a to b

D8P = 1/2*gamma*u0.*f2.*(k*(sin(theta_3)-sin(theta_1))+abs(k)*(1-sin(theta_3)*sin(theta_1))); %inter a to c




D1N = 1/2*gamma*u0.*f1.*(k*(sin(theta_1)-sin(theta_2))-abs(k)*(1-sin(theta_1)*sin(theta_2))); %intra b to a

D2N = 1/2*gamma*u0.*f1.*(k*(sin(theta_1)-sin(theta_3))-abs(k)*(1-sin(theta_1)*sin(theta_3))); %intra c to a

D3N = 1/2*gamma*u0.*f2.*(k*(sin(theta_1)-sin(theta_2))-abs(k)*(1-sin(theta_1)*sin(theta_2))); %intra a to b

D4N = 1/2*gamma*u0.*f2.*(k*(sin(theta_1)-sin(theta_3))-abs(k)*(1-sin(theta_1)*sin(theta_3))); %intra a to c



D5N = 1/2*gamma*u0.*f1.*(k*(sin(theta_2)-sin(theta_1))-abs(k)*(1-sin(theta_2)*sin(theta_1))); %inter b to a

D6N = 1/2*gamma*u0.*f1.*(k*(sin(theta_3)-sin(theta_1))-abs(k)*(1-sin(theta_3)*sin(theta_1))); %inter c to a

D7N = 1/2*gamma*u0.*f2.*(k*(sin(theta_2)-sin(theta_1))-abs(k)*(1-sin(theta_2)*sin(theta_1))); %inter a to b

D8N = 1/2*gamma*u0.*f2.*(k*(sin(theta_3)-sin(theta_1))-abs(k)*(1-sin(theta_3)*sin(theta_1))); %inter a to c


%Self dipolar interation
Dself1P = gamma*u0*Ms1*(1-2*N1)/2; 
Dself2P = gamma*u0*Ms2*(1-2*N2)/2;
Dself3P = gamma*u0*Ms2*(1-2*N2)/2;

Dself1N = -Dself1P;
Dself2N = -Dself2P;
Dself3N = -Dself3P;



DbcP1 = gamma*u0*Ms2/2*(-1+N2+N2);% b+ to c+ or c+ to b+

DbcN1 = gamma*u0*Ms2/2*(-1)+0*k; %b- to c+ or c- to b+ 0*k to convert dimension

% simply add the minus of above two equation to get b- to c-, c- to b- and
% b+ to c-, c+ to b-






%intra 
delta1P = 1/2*gamma*u0*f1.*( k*(sin(theta_1)+sin(theta_2)) - abs(k)*(1+sin(theta_1)*sin(theta_2))); % b to a 异号

delta2P = 1/2*gamma*u0*f1.*( k*(sin(theta_1)+sin(theta_3)) - abs(k)*(1+sin(theta_1)*sin(theta_3))); % c to a ==0 异号

delta3P = 1/2*gamma*u0*f2.*( - k*(sin(theta_1)+sin(theta_2)) - abs(k)*(1+sin(theta_1)*sin(theta_2))); % a to b 同号

delta4P = 1/2*gamma*u0*f2.*( - k*(sin(theta_1)+sin(theta_3)) - abs(k)*(1+sin(theta_1)*sin(theta_3))); % a to c ==0 同号

%Inter
delta5P = 1/2*gamma*u0*f1.*( - k*(sin(theta_1)+sin(theta_2)) - abs(k)*(1+sin(theta_1)*sin(theta_2))); % b to a 同号

delta6P = 1/2*gamma*u0*f1.*( - k*(sin(theta_1)+sin(theta_3)) - abs(k)*(1+sin(theta_1)*sin(theta_3))); % c to a ==0 同号

delta7P = 1/2*gamma*u0*f2.*( k*(sin(theta_1)+sin(theta_2)) - abs(k)*(1+sin(theta_1)*sin(theta_2))); % a to b 异号

delta8P = 1/2*gamma*u0*f2.*( k*(sin(theta_1)+sin(theta_3)) - abs(k)*(1+sin(theta_1)*sin(theta_3))); % a to c ==0 异号




%intra
delta1N = -1/2*gamma*u0*f1.*( - k*(sin(theta_1)+sin(theta_2)) - abs(k)*(1+sin(theta_1)*sin(theta_2))); % b to a %同号

delta2N = -1/2*gamma*u0*f1.*( - k*(sin(theta_1)+sin(theta_3)) - abs(k)*(1+sin(theta_1)*sin(theta_3))); % c to a ==0 %同号

delta3N = -1/2*gamma*u0*f2.*( k*(sin(theta_1)+sin(theta_2)) - abs(k)*(1+sin(theta_1)*sin(theta_2))); % a to b 异号

delta4N = -1/2*gamma*u0*f2.*( k*(sin(theta_1)+sin(theta_3)) - abs(k)*(1+sin(theta_1)*sin(theta_3))); % a to c ==0 异号

%inter
delta5N = -1/2*gamma*u0*f1.*( k*(sin(theta_1)+sin(theta_2)) - abs(k)*(1+sin(theta_1)*sin(theta_2))); % b to a 异号

delta6N = -1/2*gamma*u0*f1.*( k*(sin(theta_1)+sin(theta_3)) - abs(k)*(1+sin(theta_1)*sin(theta_3))); % c to a ==0 异号

delta7N = -1/2*gamma*u0*f2.*( - k*(sin(theta_1)+sin(theta_2)) - abs(k)*(1+sin(theta_1)*sin(theta_2))); % a to b%同号

delta8N = -1/2*gamma*u0*f2.*( - k*(sin(theta_1)+sin(theta_3)) - abs(k)*(1+sin(theta_1)*sin(theta_3))); % a to c ==0 %同号


%number of cells



normal_f_obc=zeros(N*6,length(k));
H_obc = zeros(N*6,N*6);
H_all_obc = zeros(N*6,N*6,length(k));
eigenvector = zeros(N*6+1,N*6,length(k));

%Matrix
for m=1:length(k)
    H_obc = zeros(N*6,N*6);

    wP = [w1P(m);w2P(m);w3P(m)];
    wN = [w1N(m);w2N(m);w3N(m)];
    DelP_intra = [delta1P(m);delta2P(m);delta3P(m);delta4P(m)];
    DelP_inter = [delta5P(m);delta6P(m);delta7P(m);delta8P(m)];
    DelN_intra = [delta1N(m);delta2N(m);delta3N(m);delta4N(m)];
    DelN_inter = [delta5N(m);delta6N(m);delta7N(m);delta8N(m)];

    DP_intra = [D1P(m);D2P(m);D3P(m);D4P(m)];
    DP_INTER = [D5P(m);D6P(m);D7P(m);D8P(m)];
    DN_intra = [D1N(m);D2N(m);D3N(m);D4N(m)];
    DN_INTER = [D5N(m);D6N(m);D7N(m);D8N(m)];
    BP = [B2P(m);B3P(m)];
    BN = [B2N(m);B3N(m)];

    DselfP = [Dself1P(m);Dself2P(m);Dself3P(m)];
    DselfN = [Dself1N(m);Dself2N(m);Dself3N(m)];
    DbcP = [DbcP1(m);DbcN1(m)];
    DbcN = [-DbcP1(m);-DbcN1(m)];

    ki=k(m);
    %constructing H @ k=ki
    for i =1:6*N
        if i==1 %1
            H_obc(i,i) = wP(1); % No1
            H_obc(i,i+1) = DP_intra(1); % No1
            H_obc(i,i+2) = DP_intra(mod(i+2,3)+2); % No2
            % H_obc(i,3*N-1) = DP_INTER(mod(i+2,3)+1); pbc % No1
            % H_obc(i,3*N) = DP_INTER(mod(i+2,3)+2); pbc % No2

            H_obc(1,3*N+i) = DselfP(1); % self

            H_obc(1,3*N+i+1) = DelP_intra(mod(i+2,3)+1); % No1
            H_obc(1,3*N+i+2) = DelP_intra(mod(i+2,3)+2); % No2
            % H_obc(1,6*N-1) = DelP_inter(mod(i+2,3)+1); pbc % No1
            % H_obc(1,6*N) = DelP_inter(mod(i+2,3)+2); pbc % No2
        end

        if (i<3*N-1) && (i>1) && (mod(i+2,3)==1) %2
            H_obc(i,i) = wP(2); %N02
            H_obc(i,i-1) = DP_intra(mod(i+2,3)+2); %N03
            H_obc(i,i+2) = DP_INTER(mod(i+2,3)+2); %N03

            H_obc(i,i+1) = DbcP(1); %antiintralayer dipolar
            H_obc(i,3*N+i) = DselfP(2); %self

            H_obc(i,3*N+i-1) = DelP_intra(mod(i+2,3)+2); %N03
            H_obc(i,3*N+i+1) = BP(mod(i+2,3))+DbcP(2); % No.1
            H_obc(i,3*N+i+2) = DelP_inter(mod(i+2,3)+2); %N03
        end

        if (i<3*N-1) && (i>1) && (mod(i+2,3)==2) %3
            H_obc(i,i) = wP(3);  %N03
            H_obc(i,i-2) = DP_intra(mod(i+2,3)+2); %N04
            H_obc(i,i+1) = DP_INTER(mod(i+2,3)+2); %N04

            H_obc(i,i-1) = DbcP(1); %antiintralayer dipolar
            H_obc(i,3*N+i) = DselfP(3); %self

            H_obc(i,3*N+i-2) = DelP_intra(mod(i+2,3)+2); %N04
            H_obc(i,3*N+i-1) = BP(mod(i+2,3))+DbcP(2);% No.2
            H_obc(i,3*N+i+1) = DelP_inter(mod(i+2,3)+2); %N04
        end

        if (i<3*N-1) && (i>1) && (mod(i+2,3)==0) %1
            H_obc(i,i-2) = DP_INTER(mod(i+2,3)+1); %N01
            H_obc(i,i-1) = DP_INTER(mod(i+2,3)+2); %N02
            H_obc(i,i) = wP(1); %N01
            H_obc(i,i+1) = DP_intra(mod(i+2,3)+1); %N01
            H_obc(i,i+2) = DP_intra(mod(i+2,3)+2); %N02

            H_obc(i,3*N+i) = DselfP(1); % self

            H_obc(i,3*N+i-2) = DelP_inter(mod(i+2,3)+1); %N01
            H_obc(i,3*N+i-1) = DelP_inter(mod(i+2,3)+2); %N02
            H_obc(i,3*N+i+1) = DelP_intra(mod(i+2,3)+1); %N01 
            H_obc(i,3*N+i+2) = DelP_intra(mod(i+2,3)+2); %N02
        end

      
        if i==3*N-1 %2 mod(i+2,3)==1
            H_obc(i,i) = wP(2); %N02
            H_obc(i,i-1) = DP_intra(mod(i+2,3)+2); %N03
            % H_obc(i,1) = DP_INTER(mod(i+2,3)+2); pbc %N03

            H_obc(i,i+1) = DbcP(1); %antiintralayer dipolar
            H_obc(i,3*N+i) = DselfP(2); %self

            H_obc(i,3*N+i-1) = DelP_intra(mod(i+2,3)+2); %N03
            H_obc(i,3*N+i+1) = BP(mod(i+2,3))+DbcP(2);% No.1
            % H_obc(i,3*N+1) = DelP_inter(mod(i+2,3)+2); pbc %N03
        end

        if i==3*N %3 mod(i+2,3)==2
            H_obc(i,i) = wP(3); %N03
            H_obc(i,i-2) = DP_intra(mod(i+2,3)+2); %N04
            % H_obc(i,1) = DP_INTER(mod(i+2,3)+2); pbc %N04

            H_obc(i,i-1) = DbcP(1); %antiintralayer dipolar
            H_obc(i,3*N+i) = DselfP(3); %self

            H_obc(i,3*N+i-2) = DelP_intra(mod(i+2,3)+2); %N04
            H_obc(i,3*N+i-1) = BP(mod(i+2,3))+DbcP(2); % No.2 
            % H_obc(i,3*N+1) = DelP_inter(mod(i+2,3)+2); pbc %N04
        end

        if i==3*N+1 %1
            H_obc(i,i) = wN(mod(i+2,3)+1);
            H_obc(i,i+1) = DN_intra(mod(i+2,3)+1);
            H_obc(i,i+2) = DN_intra(mod(i+2,3)+2);
            % H_obc(i,6*N-1) = DN_INTER(mod(i+2,3)+1); pbc
            % H_obc(i,6*N) = DN_INTER(mod(i+2,3)+2); pbc

            H_obc(i,-3*N+i) = DselfN(1); % self

            H_obc(i,-3*N+i+1) = DelN_intra(mod(i+2,3)+1);
            H_obc(i,-3*N+i+2) = DelN_intra(mod(i+2,3)+2);
            % H_obc(i,3*N-1) = DelN_inter(mod(i+2,3)+1); pbc
            % H_obc(i,3*N) = DelN_inter(mod(i+2,3)+2); pbc
        end

        if (i<6*N-1) && (i>3*N+1) && (mod(i+2,3)==1) %2
            H_obc(i,i) = wN(mod(i+2,3)+1);
            H_obc(i,i-1) = DN_intra(mod(i+2,3)+2);
            H_obc(i,i+2) = DN_INTER(mod(i+2,3)+2); 

            H_obc(i,i+1) = DbcN(1); %antiintralayer dipolar
            H_obc(i,-3*N+i) = DselfN(2); %self
           
            H_obc(i,-3*N+i-1) = DelN_intra(mod(i+2,3)+2);
            H_obc(i,-3*N+i+1) = BN(mod(i+2,3))+DbcN(2);
            H_obc(i,-3*N+i+2) = DelN_inter(mod(i+2,3)+2);
        end

        if (i<6*N-1) && (i>3*N+1) && (mod(i+2,3)==2) %3
            H_obc(i,i) = wN(mod(i+2,3)+1);
            H_obc(i,i-2) = DN_intra(mod(i+2,3)+2);
            H_obc(i,i+1) = DN_INTER(mod(i+2,3)+2); 

            H_obc(i,i-1) = DbcN(1); %antiintralayer dipolar
            H_obc(i,-3*N+i) = DselfN(3); %self

            H_obc(i,-3*N+i-2) = DelN_intra(mod(i+2,3)+2);
            H_obc(i,-3*N+i-1) = BN(mod(i+2,3))+DbcN(2);
            H_obc(i,-3*N+i+1) = DelN_inter(mod(i+2,3)+2);
        end

        if (i<6*N-1) && (i>3*N+1) && (mod(i+2,3)==0) %1
            H_obc(i,i-2) = DN_INTER(mod(i+2,3)+1);
            H_obc(i,i-1) = DN_INTER(mod(i+2,3)+2);
            H_obc(i,i) = wN(mod(i+2,3)+1);
            H_obc(i,i+1) = DN_intra(mod(i+2,3)+1);
            H_obc(i,i+2) = DN_intra(mod(i+2,3)+2);

            H_obc(i,-3*N+i) = DselfN(1); % self

            H_obc(i,-3*N+i-2) = DelN_inter(mod(i+2,3)+1);
            H_obc(i,-3*N+i-1) = DelN_inter(mod(i+2,3)+2);
            H_obc(i,-3*N+i+1) = DelN_intra(mod(i+2,3)+1); 
            H_obc(i,-3*N+i+2) = DelN_intra(mod(i+2,3)+2); 
        end

        if i==6*N-1 %2
            H_obc(i,i) = wN(mod(i+2,3)+1);
            H_obc(i,i-1) = DN_intra(mod(i+2,3)+2);
            % H_obc(i,1+3*N) = DN_INTER(mod(i+2,3)+2); pbc

            H_obc(i,i+1) = DbcN(1); %antiintralayer dipolar
            H_obc(i,-3*N+i) = DselfN(2); %self

            H_obc(i,-3*N+i-1) = DelN_intra(mod(i+2,3)+2);
            H_obc(i,-3*N+i+1) = BN(mod(i+2,3))+DbcN(2);
            % H_obc(i,1) = DelN_inter(mod(i+2,3)+2); pbc
        end

        if i==6*N %3
            H_obc(i,i) = wN(mod(i+2,3)+1);
            H_obc(i,i-2) = DN_intra(mod(i+2,3)+2);
            % H_obc(i,3*N+1) = DN_INTER(mod(i+2,3)+2); pbc

            H_obc(i,i-1) = DbcN(1); %antiintralayer dipolar
            H_obc(i,-3*N+i) = DselfN(3); %self

            H_obc(i,-3*N+i-2) = DelN_intra(mod(i+2,3)+2);
            H_obc(i,-3*N+i-1) = BN(mod(i+2,3))+DbcN(2);
            % H_obc(i,1) = DelN_inter(mod(i+2,3)+2); pbc
        end

        
    end 
    H_all_obc(:,:,m) = H_obc;
    % Compute eigenvalues and eigenvectors
    [V, D] = eig(H_obc); % Compute eigenvalues and eigenvectors
    eigvals = diag(D);   % Extract eigenvalues into a vector
    eigvecs = V;         % Extract eigenvectors into columns of V
    
    % Store eigenvalues and eigenvectors
    ALL_eigenvector(1,:,m) = eigvals / (2 * pi * 1e9);
    ALL_eigenvector(2:6*N+1,:,m) = eigvecs;
    
    % Tolerance to check if eigenvalue is real
    tol = 0.1e9;
    
    % Preallocate storage for 6 groups
    grouped_eigenvalues = NaN(6, size(H_obc, 1)); % 6 rows for the 6 groups
    grouped_eigenvectors = NaN(6, size(H_obc, 1), size(H_obc, 2)); % Store eigenvectors in 6 groups
    
    % Process each eigenvalue and its corresponding eigenvector
    for eig_idx = 1:length(eigvals)
        current_eigval = eigvals(eig_idx);        % Current eigenvalue
        current_eigvec = eigvecs(:, eig_idx);     % Corresponding eigenvector
    
        % Check if the eigenvalue is real within the tolerance
        if abs(imag(current_eigval)) < tol
            % Find the index of the largest component in the eigenvector
            [~, max_idx] = max(abs(current_eigvec));
            
            % Determine the group based on max_idx mod 6
            % Determine the group based on the updated rules
        if max_idx <= 3*N
            % Index is in the range [1, 2N]
            if mod(max_idx, 3) == 1
                group_id = 1; % mode=1 indices -> Group 1
            elseif mod(max_idx, 3) == 2
                group_id = 2; % mod=2 indices -> Group 2
            else
                group_id = 3; % mod=0 indices -> Group 3
            end
        else
            % Index is in the range [2N+1, 4N]
            if mod(max_idx, 3) == 1
                group_id = 4; % mode=1 indices -> Group 4
            elseif mod(max_idx, 3) == 2
                group_id = 5; % mod=2 indices -> Group 5
            else
                group_id = 6; % mode=0 indices -> Group 6
            end
        end
            

            % Assign the eigenvalue to the group
            grouped_eigenvalues(group_id, eig_idx) = (current_eigval) / (2 * pi * 1e9);
            
            % Assign the eigenvector to the group
            grouped_eigenvectors(group_id, :, eig_idx) = current_eigvec;
            
        end
    end
    All_grouped_eigenvalues(:, :,m) = grouped_eigenvalues;
    All_grouped_eigenvectors(:, :, :,m) = grouped_eigenvectors;
    
end





% Define colors for each band
dark_blue = [68, 113, 196]/255;
dark_blue_anti = [161, 184, 225]/255;
light_blue = [0, 175, 239]/255;
light_blue_anti = [178, 230, 250]/255;
orange = [236, 124, 48]/255;
orange_anti = [249, 222, 189]/255;


% Define marker types for magnons and antimagnons
magnon_marker = 'o';          % Solid circle for magnons
antimagnon_marker = 'o';      % Hollow circle for antimagnons
antimagnon_face = 'none';     % No fill for antimagnons

% Define alpha values for transparency
antimagnon_alpha = 1;       % 90% transparent
magnon_alpha = 1;           % Fully opaque for magnons

% Plotting the band structure using scatter plots
figure;
ax = axes('Units', 'inches', 'Position', [2 2 4.5 3.5]); %left=1, bottom=1, width=4 inches, height=3 inches
ax.ActivePositionProperty = 'position';
hold on;
% Plot a-magnon
scatter(k, squeeze(All_grouped_eigenvalues(1,:,:)), 15, 'MarkerEdgeColor', dark_blue, ...
    'MarkerFaceColor', dark_blue, 'Marker', magnon_marker, ...
    'MarkerFaceAlpha', magnon_alpha, 'MarkerEdgeAlpha', magnon_alpha, ...
    'DisplayName', 'a-magnon');

% Plot b-magnon
scatter(k, squeeze(All_grouped_eigenvalues(2,:,:)), 15, 'MarkerEdgeColor', light_blue, ...
    'MarkerFaceColor', light_blue, 'Marker', magnon_marker, ...
    'MarkerFaceAlpha', magnon_alpha, 'MarkerEdgeAlpha', magnon_alpha, ...
    'DisplayName', 'c-magnon');

% Plot c-magnon
scatter(k, squeeze(All_grouped_eigenvalues(3,:,:)), 15, 'MarkerEdgeColor', orange, ...
    'MarkerFaceColor', orange, 'Marker', magnon_marker, ...
    'MarkerFaceAlpha', magnon_alpha, 'MarkerEdgeAlpha', magnon_alpha, ...
    'DisplayName', 'c-magnon');

% Plot a-antimagnon with 90% transparency
scatter(k, squeeze(All_grouped_eigenvalues(4,:,:)), 15, 'MarkerEdgeColor', dark_blue_anti, ...
    'MarkerFaceColor', antimagnon_face, 'Marker', antimagnon_marker, ...
    'MarkerFaceAlpha', antimagnon_alpha, 'MarkerEdgeAlpha', antimagnon_alpha, ...
    'DisplayName', 'a-antimagnon');

% Plot b-antimagnon
scatter(k, squeeze(All_grouped_eigenvalues(5,:,:)), 15, 'MarkerEdgeColor', light_blue_anti, ...
    'MarkerFaceColor', light_blue_anti, 'Marker', magnon_marker, ...
    'MarkerFaceAlpha', magnon_alpha, 'MarkerEdgeAlpha', magnon_alpha, ...
    'DisplayName', 'c-magnon');


% Plot c-antimagnon with 90% transparency
scatter(k, squeeze(All_grouped_eigenvalues(6,:,:)), 15, 'MarkerEdgeColor', orange_anti, ...
    'MarkerFaceColor', antimagnon_face, 'Marker', antimagnon_marker, ...
    'MarkerFaceAlpha', antimagnon_alpha, 'MarkerEdgeAlpha', antimagnon_alpha, ...
    'DisplayName', 'c-antimagnon');

hold off;


% === 5. Customize Axis Labels ===
xlabel('$k_x \ (m^{-1})$', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial', ...
    'FontWeight', 'bold', ...
    'FontSize', 28);

ylabel('$\omega / 2\pi \ (GHz$)', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial', ...
    'FontWeight', 'bold', ...
    'FontSize', 28);
% title('Band Structure', 'FontSize', 16);

% Create Dummy Plot Handles for Legend with Solid Markers
hold on; % Ensure these dummy plots don’t disrupt existing plots

antimagnon_alpha1 = 1;       % 90% transparent
magnon_alpha1 = 1;           % Fully opaque for magnons
dark_blue = [68, 113, 196]/255;
orange = [236, 124, 48]/255;

% Dummy plots for magnons (solid markers)
h_dummy_a_magnon = scatter(NaN, NaN, 36, dark_blue, magnon_marker, 'filled', ...
    'MarkerFaceAlpha', magnon_alpha1, 'MarkerEdgeAlpha', magnon_alpha, 'MarkerFaceColor', dark_blue, 'MarkerEdgeColor', dark_blue);

h_dummy_b_magnon = scatter(NaN, NaN, 36, light_blue, magnon_marker, 'filled', ...
    'MarkerFaceAlpha', magnon_alpha1, 'MarkerEdgeAlpha', magnon_alpha, 'MarkerFaceColor', light_blue, 'MarkerEdgeColor', light_blue);

h_dummy_c_magnon = scatter(NaN, NaN, 36, orange, magnon_marker, 'filled', ...
    'MarkerFaceAlpha', magnon_alpha1, 'MarkerEdgeAlpha', magnon_alpha, 'MarkerFaceColor', orange, 'MarkerEdgeColor', orange);

% Dummy plots for antimagnons (transparent markers)


h_dummy_a_antimagnon = scatter(NaN, NaN, 36, dark_blue_anti, antimagnon_marker, 'filled', ...
    'MarkerFaceAlpha', antimagnon_alpha1, 'MarkerEdgeAlpha', antimagnon_alpha1, 'MarkerFaceColor', dark_blue_anti, 'MarkerEdgeColor', dark_blue_anti);

h_dummy_b_antimagnon = scatter(NaN, NaN, 36, light_blue_anti, antimagnon_marker, 'filled', ...
    'MarkerFaceAlpha', antimagnon_alpha1, 'MarkerEdgeAlpha', antimagnon_alpha1, 'MarkerFaceColor', light_blue_anti, 'MarkerEdgeColor', light_blue_anti);


h_dummy_c_antimagnon = scatter(NaN, NaN, 36, orange_anti, antimagnon_marker, 'filled', ...
    'MarkerFaceAlpha', antimagnon_alpha1, 'MarkerEdgeAlpha', antimagnon_alpha1, 'MarkerFaceColor', orange_anti, 'MarkerEdgeColor', orange_anti);



% Turn off visibility of dummy plots
set([h_dummy_a_magnon, h_dummy_b_magnon, h_dummy_c_magnon, ...
     h_dummy_a_antimagnon, h_dummy_b_antimagnon,  h_dummy_c_antimagnon], 'Visible', 'on');

% Customize the legend using dummy handles
legend_handles = [h_dummy_a_magnon,  h_dummy_b_magnon, h_dummy_c_magnon, ...
                  h_dummy_a_antimagnon, h_dummy_b_antimagnon,  h_dummy_c_antimagnon];
legend_labels = {'a-magnon', 'b-magnon','c-magnon', ...
                'a-antimagnon', 'b-antimagnon','c-antimagnon'};

legend(legend_handles, legend_labels, ...
    'Location', 'best', ...
    'Box', 'off', ...                      % Remove the black frame
    'FontSize', 28, ...                    % Set desired font size
    'FontName', 'Arial', ...               % Set desired font type
    'FontWeight', 'bold', ...              % Optional: set font weight
    'TextColor', 'black', ...                  % Set text color to black
    'Interpreter', 'latex');               % Optional: use LaTeX interpreter



ylim([0 150])
xlim([-4e8 4e8])

hold off;
%%% End of plotting simulated data

set(gca, 'LineWidth', 2.5,...
         'FontName', 'Arial', ...
         'FontSize', 30, ...                    % Increased font size for ticks
         'FontWeight', 'bold', ...              % Make ticks bold
         'TickLabelInterpreter', 'latex', ...
         'XAxisLocation', 'bottom', ...
         'YAxisLocation', 'left');      % Use LaTeX interpreter for tick labels
box on;