%%Figure S3
clc; clear all; close all

%  Definition of parameters
kx = linspace(-1e8,1e8,1e3);
ky = linspace(0,0,1);
index=0;
d1 = 25e-9;
d2 = 25e-9;
Ms1 = 1.4e5;
Ms2 = Ms1; 
Aex1 = 3.7e-12*0;
Aex2 = 8.7e-12*0;
u0 = 4*pi*1e-7;
gamma = 28e9*2*pi;



% overall effective external field
H1 = 1e5;
H2 = 1e5 ;


%number of cells
N=20;

angluar_f_obc=zeros(length(kx),length(ky),N*4);


%Matrix
for m=1:length(ky)
    for n = 1:length(kx)
        k = sqrt((ky(m))^2 + (kx(n))^2);

        N1 = (1-exp(-abs(k)*d1))./(abs(k)*d1);
        N2 = (1-exp(-abs(k)*d2))./(abs(k)*d2);

        f2 = Ms1*d1*N1.*N2/2;
        f1 = Ms2*d2*N1.*N2/2;

        w1P = (2*Aex1*gamma*k^2/Ms1 + gamma*u0*H1) + gamma*u0*Ms1/2*((1-N1)*(kx(n))^2/k^2+N1);
        w1N = -(2*Aex1*gamma*k^2/Ms1 + gamma*u0*H1) - gamma*u0*Ms1/2*((1-N1)*(kx(n))^2/k^2+N1);
        w2P = (2*Aex2*gamma*k^2/Ms2 + gamma*u0*H2) + gamma*u0*Ms2/2*((1-N2)*(kx(n))^2/k^2+N2);
        w2N = -(2*Aex2*gamma*k^2/Ms2 + gamma*u0*H2) - gamma*u0*Ms2/2*((1-N2)*(kx(n))^2/k^2+N2);


        D1P = 0/2*gamma*u0.*f1.* ((k+(kx(n))^2/k)+2*kx(n));
        D2P = 0/2*gamma*u0.*f2.* ((k+(kx(n))^2/k)+2*kx(n));

        D3P = 0/2*gamma*u0.*f2.* ((k+(kx(n))^2/k)-2*kx(n));
        D4P = 0/2*gamma*u0.*f1.* ((k+(kx(n))^2/k)-2*kx(n));

        D1N = 0/2*gamma*u0.*f1.* (-(k+(kx(n))^2/k)+2*kx(n));
        D2N = 0/2*gamma*u0.*f2.* (-(k+(kx(n))^2/k)+2*kx(n));
        D3N = 0/2*gamma*u0.*f2.* (-(k+(kx(n))^2/k)-2*kx(n));
        D4N = 0/2*gamma*u0.*f1.* (-(k+(kx(n))^2/k)-2*kx(n));

        delta1P_ = 0/2*gamma*u0*(f1*(k-(kx(n))^2/k))+1/2*gamma*u0.*f1.* (-(k+(kx(n))^2/k)+2*kx(n)); %interlayer intracell %
        delta2P_ = 0/2*gamma*u0*(f2*(k-(kx(n))^2/k))+1/2*gamma*u0.*f2.* (-(k+(kx(n))^2/k)-2*kx(n));%
        delta3P_ = 0/2*gamma*u0*(f1*(k-(kx(n))^2/k))+1/2*gamma*u0.*f1.* (-(k+(kx(n))^2/k)-2*kx(n));%intralayer 
        delta4P_ = 0/2*gamma*u0*(f2*(k-(kx(n))^2/k))+1/2*gamma*u0.*f2.* (-(k+(kx(n))^2/k)+2*kx(n));%




        delta1N_ = -0/2*gamma*u0*(f1*(k-(kx(n))^2/k))-1/2*gamma*u0.*f1.* (-(k+(kx(n))^2/k)-2*kx(n));%
        delta2N_ = -0/2*gamma*u0*(f2*(k-(kx(n))^2/k))-1/2*gamma*u0.*f2.* (-(k+(kx(n))^2/k)+2*kx(n));%
        delta3N_ = -0/2*gamma*u0*(f1*(k-(kx(n))^2/k))-1/2*gamma*u0.*f1.* (-(k+(kx(n))^2/k)+2*kx(n));%
        delta4N_ = -0/2*gamma*u0*(f2*(k-(kx(n))^2/k))-1/2*gamma*u0.*f2.* (-(k+(kx(n))^2/k)-2*kx(n));%

        DELTA1P = gamma*u0*Ms1/2*((1-N1)*(kx(n)/k)^2-N1); %intralayer
        DELTA2P = gamma*u0*Ms2/2*((1-N2)*(kx(n)/k)^2-N2);
        DELTA1N = -DELTA1P;
        DELTA2N = -DELTA2P;
    
        wP = [w1P;w2P];
        wN = [w1N;w2N];

        DelP_intra_ = [delta1P_;delta2P_];
        DelP_inter_ = [delta3P_;delta4P_];
        DelN_intra_ = [delta1N_;delta2N_];
        DelN_inter_ = [delta3N_;delta4N_];

        DELP = [DELTA1P;DELTA2P];
        DELN = [DELTA1N;DELTA2N];

        DP_intra = [D1P;D2P];
        DP_INTER = [D3P;D4P];
        DN_intra = [D1N;D2N];
        DN_INTER = [D3N;D4N];

        H_obc = zeros(N*4,N*4);

        for i =1:4*N
            num = mod(i,2)+1; %i=even=1 i=odd=2
            nrest = 1.5+(-1)^(num+1)/2; %i=even=2 i=odd=1
            if i==1
                H_obc(i,1) = wP(nrest);
                H_obc(i,2) = DP_intra(nrest);
                % H_obc(i,2*N) = DP_INTER(num);%pbc i is odd so num=2 D4P
                H_obc(i,2*N+1) = DELP(nrest); 
                H_obc(i,2*N+2) = DelP_intra_(nrest); %index same as DP_intra
                %H_obc(i,4*N) = DelP_inter_(nrest);%pbc
            end

            if (i<2*N) && (i>1) && (mod(i,2)==0) %if i=even
                H_obc(i,i) = wP(nrest);
                H_obc(i,i-1) = DP_intra(nrest);
                H_obc(i,i+1) = DP_INTER(num); % row 2 is D3P at i+1 col
                H_obc(i,2*N+i) = DELP(nrest);
                H_obc(i,2*N+i-1) = DelP_intra_(nrest);
                H_obc(i,2*N+i+1) = DelP_inter_(nrest);           
            end

            if (i<2*N) && (i>1) && (mod(i,2)==1) %if i=odd
                H_obc(i,i) = wP(nrest); 
                H_obc(i,i-1) = DP_INTER(num); % row 3 is D4P at i-1 col
                H_obc(i,i+1) = DP_intra(nrest);
                H_obc(i,2*N+i) = DELP(nrest);
                H_obc(i,2*N+i-1) = DelP_inter_(nrest);
                H_obc(i,2*N+i+1) = DelP_intra_(nrest);  
            end

            if i==2*N
                H_obc(i,i) = wP(nrest);
                % H_obc(i,1)=DP_INTER(num);%pbc i is even so num =1 D3P
                H_obc(i,i-1) = DP_intra(nrest);
                H_obc(i,2*N+i) = DELP(nrest);
                H_obc(i,2*N+i-1) = DelP_intra_(nrest); 
                % H_obc(i,2*N+1) = DelP_inter_(nrest);%pbc
            end
    
            if i==2*N+1
                H_obc(i,i) = wN(nrest);
                H_obc(i,i+1) = DN_intra(nrest);
                % H_obc(i,4*N) = DN_INTER(num);%pbc
                H_obc(i,i-2*N) = DELN(nrest);
                H_obc(i,i-2*N+1) = DelN_intra_(nrest);
                % H_obc(i,2*N) = DelN_inter_(nrest);%pbc
            end
    
            if (i>2*N+1) && (i<4*N) && (mod(i,2)==0) %if i=even
                H_obc(i,i) = wN(nrest);
                H_obc(i,i-1) = DN_intra(nrest);
                H_obc(i,i+1) = DN_INTER(num);
                H_obc(i,i-2*N) = DELN(nrest);
                H_obc(i,i-2*N-1) = DelN_intra_(nrest);
                H_obc(i,i-2*N+1) = DelN_inter_(nrest);
            end
    
            if (i>2*N+1) && (i<4*N) && (mod(i,2)==1) %if i=odd
                H_obc(i,i) = wN(nrest);
                H_obc(i,i-1) = DN_INTER(num);
                H_obc(i,i+1) = DN_intra(nrest);
                H_obc(i,i-2*N) = DELN(nrest);
                H_obc(i,i-2*N-1) = DelN_inter_(nrest);
                H_obc(i,i-2*N+1) = DelN_intra_(nrest);
            end
    
            if i==4*N
                H_obc(i,i) = wN(nrest);
                H_obc(i,i-1) = DN_intra(nrest);
                % H_obc(i,2*N+1) = DN_INTER(num); %pbc
                H_obc(i,i-2*N) = DELN(nrest); 
                H_obc(i,i-2*N-1) = DelN_intra_(nrest);
                % H_obc(i,1) = DelN_inter_(nrest); %pbc
            end
            egv = eig(H_obc);
            angluar_f_obc(n,m,:) = egv/2/pi;
        end 
        H_all_obc(:,:,n) = H_obc;


        % Compute eigenvalues and eigenvectors
        [V, D] = eig(H_obc); % Compute eigenvalues and eigenvectors
        eigvals = diag(D);   % Extract eigenvalues into a vector
        eigvecs = V;         % Extract eigenvectors into columns of V
        
        % Store eigenvalues and eigenvectors
        ALL_eigenvector(1,:,n) = eigvals / (2 * pi * 1e9);
        ALL_eigenvector(2:4*N+1,:,n) = eigvecs;
        
        % Tolerance to check if eigenvalue is real
        tol = 1e-10;
        
        % Preallocate storage for 4 groups (mod 4 = 1, 2, 3, 0)
        grouped_eigenvalues = NaN(4, size(H_obc, 1)); % 4 rows for the 4 groups
        grouped_eigenvectors = NaN(4, size(H_obc, 1), size(H_obc, 2)); % Store eigenvectors in 4 groups


        %0425 second version: find max
        new_eigenvalues = NaN(2, size(H_obc, 1)); % 2 rows for the 2 groups
        new_eigenvectors = NaN(2, size(H_obc, 1), size(H_obc, 2)); % Store eigenvectors in 2 groups
        % Process each eigenvalue and its corresponding eigenvector
        [max_eigenvalue, max_value_idx] = max(eigvals);
        ALL_max_eigenvector(1,n) = max_eigenvalue / (2 * pi * 1e9);
        ALL_max_eigenvector(2:4*N+1,n) = eigvecs(:,max_value_idx);


        for eig_idx = 1:length(eigvals)
            current_eigval = eigvals(eig_idx);        % Current eigenvalue
            current_eigvec = eigvecs(:, eig_idx);     % Corresponding eigenvector
        
            % Check if the eigenvalue is real within the tolerance
            if abs(imag(current_eigval)) < tol
                % Find the index of the largest component in the eigenvector
                [max_value, max_idx] = max(abs(current_eigvec));
                
                % Determine the group based on max_idx mod 4
                % Determine the group based on the updated rules

            %0425 plot edge state    
            if max_idx == 1 & max_value>=0.48 
                new_group_id = 1;
                % Assign the eigenvalue to the group
                new_eigenvalues(new_group_id, eig_idx) = (current_eigval) / (2 * pi * 1e9);
                % Assign the eigenvector to the group
                new_eigenvectors(new_group_id, :, eig_idx) = current_eigvec;
            end
            if max_idx == 2*N & max_value>=0.48 
                new_group_id = 2;
                % Assign the eigenvalue to the group
                new_eigenvalues(new_group_id, eig_idx) = (current_eigval) / (2 * pi * 1e9);
                % Assign the eigenvector to the group
                new_eigenvectors(new_group_id, :, eig_idx) = current_eigvec;
            end

            

            if max_idx <= 2*N
                % Index is in the range [1, 2N]
                if mod(max_idx, 2) == 1
                    group_id = 1; % Odd indices -> Group 1
                else
                    group_id = 2; % Even indices -> Group 2
                end
            else
                % Index is in the range [2N+1, 4N]
                if mod(max_idx, 2) == 1
                    group_id = 3; % Odd indices -> Group 3
                else
                    group_id = 4; % Even indices -> Group 4
                end
            end
                

                % Assign the eigenvalue to the group
                grouped_eigenvalues(group_id, eig_idx) = (current_eigval) / (2 * pi * 1e9);
                
                % Assign the eigenvector to the group
                grouped_eigenvectors(group_id, :, eig_idx) = current_eigvec;
                
            end
        end
        All_grouped_eigenvalues(:, :,n) = grouped_eigenvalues;
        All_grouped_eigenvectors(:, :, :,n) = grouped_eigenvectors;

        NEW_All_eigenvalues(:, :,n) = new_eigenvalues;
        NEW_All_eigenvectors(:, :, :,n) = new_eigenvectors;
        
    end
end

%%%%%%%Finish Band Structure Calculation

% Define colors for each band
dark_blue = [68, 113, 196]/255;
dark_blue_anti = [161, 184, 225]/255;
light_blue = [0, 175, 239]/255;
light_blue_anti = [178, 230, 250]/255;
orange = [236, 124, 48]/255;
orange_anti = [249, 222, 189]/255;
green = [0, 255, 0]/255;
blue = [0, 0, 255]/255;
gray = [160 160 160]/255;

% Define marker types for magnons and antimagnons
magnon_marker = 'o';          % Solid circle for magnons
antimagnon_marker = 'o';      % Hollow circle for antimagnons
antimagnon_face = 'none';     % No fill for antimagnons

% Define alpha values for transparency
antimagnon_alpha = 1;       % 90% transparent
magnon_alpha = 1;           % Fully opaque for magnons




%%%%%
figure;
ax = axes('Units', 'inches', 'Position', [2 2 4.5 3.5]); %left=1, bottom=1, width=4 inches, height=3 inches
ax.ActivePositionProperty = 'position';
idx_negative = (kx < 0);
idx_positive = (kx > 0);
scatter(kx, squeeze(ALL_eigenvector(1,:,:)), 15, 'MarkerEdgeColor', gray, ...
    'MarkerFaceColor', gray, 'Marker', antimagnon_marker, ...
    'MarkerFaceAlpha', antimagnon_alpha, 'MarkerEdgeAlpha', antimagnon_alpha, ...
    'DisplayName', 'c-magnon');
hold on;
scatter(kx(idx_negative), squeeze(ALL_max_eigenvector(1, idx_negative)), 15, ...
    'MarkerEdgeColor', green, 'MarkerFaceColor', green, ...
    'Marker', magnon_marker, ...
    'MarkerFaceAlpha', magnon_alpha, 'MarkerEdgeAlpha', magnon_alpha, ...
    'DisplayName', 'a-magnon (kx<0)');
hold on;

% Plot kx > 0 part
scatter(kx(idx_positive), squeeze(ALL_max_eigenvector(1, idx_positive)), 15, ...
    'MarkerEdgeColor', blue, 'MarkerFaceColor', blue, ...
    'Marker', magnon_marker, ...
    'MarkerFaceAlpha', magnon_alpha, 'MarkerEdgeAlpha', magnon_alpha, ...
    'DisplayName', 'c-magnon (kx>0)');

k = linspace(0,1e8,1e3+1);
d = (d1+d2)*2;
wH = gamma*u0*H1;
wM = gamma*u0*Ms1;
MSSW_theory = sqrt((wH+wM/2)^2-(wM/2)^2*exp(-2*k*d))/1e9/2/pi;
plot(k, MSSW_theory, 'r--', 'LineWidth', 2.5);   

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
set(gca, 'LineWidth', 2.5,...
         'FontName', 'Arial', ...
         'FontSize', 28, ...                    % Increased font size for ticks 32 for small， 28 big
         'FontWeight', 'bold', ...              % Make ticks bold
         'TickLabelInterpreter', 'latex', ...
         'XAxisLocation', 'bottom', ...
         'YAxisLocation', 'left');      % Use LaTeX interpreter for tick labels
box on;
ylim([5 6])
%%%%%%
figure;
ax = axes('Units', 'inches', 'Position', [2 2 4.5 3.5]); %left=1, bottom=1, width=4 inches, height=3 inches
ax.ActivePositionProperty = 'position'; 
% data_flipped = flipud(ALL_max_eigenvector(2:2*N+1,:)); 
im = imagesc(kx,1:2*N,abs(ALL_max_eigenvector(2:2*N+1,:)).^2);
c = [linspace(0.5, 1, 256)', linspace(0, 1, 256)', linspace(0.5, 1, 256)'];
c = flipud(c);
colormap(c);
% im.AlphaData =0.25;
caxis([0 1]);
cb = colorbar;           % Adds a colorbar to the current axes
% cb.Ticks = [0 0.5 1];    % Sets tick positions at 0, 0.5, and 1
% cb.TickLabels = {'0', '0.5', '1'};  % (Optional) Sets custom labels for the ticks
cb.FontName = 'Arial';
cb.TickLabelInterpreter = 'latex';
cb.FontWeight = 'bold';
cb.FontSize = 28;
cb.Ticks = [0 0.5 1];
cb.TickLabels = {'0', '0.5', '1'};
% Correct the y-axis direction
set(gca, 'YDir', 'normal'); % Make y increase from bottom to top
xlabel('$k_x \ (m^{-1})$', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial', ...
    'FontWeight', 'bold', ...
    'FontSize', 28);

ylabel('Layer  index', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial', ...
    'FontWeight', 'bold', ...
    'FontSize', 28);
set(gca, 'LineWidth', 2.5,...
         'FontName', 'Arial', ...
         'FontSize', 28, ...                    % Increased font size for ticks 32 for small， 28 big
         'FontWeight', 'bold', ...              % Make ticks bold
         'TickLabelInterpreter', 'latex', ...
         'XAxisLocation', 'bottom', ...
         'YAxisLocation', 'left');      % Use LaTeX interpreter for tick labels
ylim([0 40.8])
box on;
%%%%%%%%



