%% Figure 3
clc; clear all; close all

%  Definition of parameters
kx = linspace(-2e8,2e8,3e4+1);
ky = linspace(0,0,1);
index=0;

d1 = 10e-9;
d2 = 10e-9;
Ms1 = 1.4e5;
Ms2 = 7.4e5;
Aex1 = 3.7e-12;
Aex2 = 8.7e-12;
u0 = 4*pi*1e-7;
gamma = 28e9*2*pi;



H1 = -7e5;
H2 = -8.5e5 ;

N=3;

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

        delta1P_ = 1/2*gamma*u0*(f1*(k-(kx(n))^2/k))+1/2*gamma*u0.*f1.* (-(k+(kx(n))^2/k)+2*kx(n)); %interlayer intracell 

        delta2P_ = 1/2*gamma*u0*(f2*(k-(kx(n))^2/k))+1/2*gamma*u0.*f2.* (-(k+(kx(n))^2/k)-2*kx(n));%
        delta3P_ = 1/2*gamma*u0*(f1*(k-(kx(n))^2/k))+1/2*gamma*u0.*f1.* (-(k+(kx(n))^2/k)-2*kx(n));%intralayer 
        delta4P_ = 1/2*gamma*u0*(f2*(k-(kx(n))^2/k))+1/2*gamma*u0.*f2.* (-(k+(kx(n))^2/k)+2*kx(n));%


        delta1N_ = -1/2*gamma*u0*(f1*(k-(kx(n))^2/k))-1/2*gamma*u0.*f1.* (-(k+(kx(n))^2/k)-2*kx(n));%
        delta2N_ = -1/2*gamma*u0*(f2*(k-(kx(n))^2/k))-1/2*gamma*u0.*f2.* (-(k+(kx(n))^2/k)+2*kx(n));%
        delta3N_ = -1/2*gamma*u0*(f1*(k-(kx(n))^2/k))-1/2*gamma*u0.*f1.* (-(k+(kx(n))^2/k)+2*kx(n));%
        delta4N_ = -1/2*gamma*u0*(f2*(k-(kx(n))^2/k))-1/2*gamma*u0.*f2.* (-(k+(kx(n))^2/k)-2*kx(n));%

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
        
        % Process each eigenvalue and its corresponding eigenvector
        for eig_idx = 1:length(eigvals)
            current_eigval = eigvals(eig_idx);        % Current eigenvalue
            current_eigvec = eigvecs(:, eig_idx);     % Corresponding eigenvector
        
            % Check if the eigenvalue is real within the tolerance
            if abs(imag(current_eigval)) < tol
                % Find the index of the largest component in the eigenvector
                [~, max_idx] = max(abs(current_eigvec));
                
                % Determine the group based on max_idx mod 4
                % Determine the group based on the updated rules
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
        
    end
end

%%%%%%%Finish Band Structure Calculation

%%%%%%%%Start to extract and fft COMSOL's data
%coordinate setting
ksCOM2 = 1/(3e-9); %sampling frequency for x
xCOM2 = linspace(-600*4,600*4,1601); %
% fCOM2 = 0.25:0.25:60;
fCOM2 = 39:-0.15:0;
kxCOM2 = 2*pi*((0:length(xCOM2)-1)*ksCOM2/length(xCOM2)-ksCOM2/2);


%Py top data reading
table2 = readtable('Real_YIG.xlsx');
table_im2 = readtable('Im_YIG.xlsx');

table1 = readtable('Real_Py.xlsx');
table_im1 = readtable('Im_Py.xlsx');

table4 = readtable('Real_YIG - next_cell.xlsx');
table_im4 = readtable('Im_YIG - next_cell.xlsx');

table3 = readtable('Real_Py - next_cell.xlsx');
table_im3 = readtable('Im_Py - next_cell.xlsx');



A1COM2 = table2array(table2);
A2COM2 = table2array(table_im2);
A1COM1 = table2array(table1);
A2COM1 = table2array(table_im1);
A1COM4 = table2array(table4);
A2COM4 = table2array(table_im4);
A1COM3 = table2array(table3);
A2COM3 = table2array(table_im3);

% Replace NaN values with 0 explicitly
A1COM2(isnan(A1COM2)) = 0;
A2COM2(isnan(A2COM2)) = 0;
A1COM1(isnan(A1COM1)) = 0;
A2COM1(isnan(A2COM1)) = 0;
A1COM4(isnan(A1COM4)) = 0;
A2COM4(isnan(A2COM4)) = 0;
A1COM3(isnan(A1COM3)) = 0;
A2COM3(isnan(A2COM3)) = 0;

% Combine real and imaginary parts after NaN removal
A1COM2 = A1COM2 + 1j * A2COM2;
A1COM1 = A1COM1 + 1j * A2COM1;
A1COM4 = A1COM4 + 1j * A2COM4;
A1COM3 = A1COM3 + 1j * A2COM3;


%FFT for A1
for i = 1:length(fCOM2)
    
    DP1COM2(:,i) = fftshift(fft(A1COM2(:,i)));
    DP1COM1(:,i) = fftshift(fft(A1COM1(:,i)));
    DP1COM4(:,i) = fftshift(fft(A1COM4(:,i)));
    DP1COM3(:,i) = fftshift(fft(A1COM3(:,i)));
end
z1COM3 = abs(flipud(transpose(DP1COM3)));%Top Py
z1COM4 = abs(flipud(transpose(DP1COM4)));%Top YIG
z1COM1 = abs(flipud(transpose(DP1COM1)));%Bottom Py
z1COM2 = abs(flipud(transpose(DP1COM2)));%Bottom YIG
%%%%%%%%End of extract and fft COMSOL's data




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
%BIG
ax = axes('Units', 'inches', 'Position', [2 2 4.5 3.5]); %left=1, bottom=1, width=4 inches, height=3 inches
%small
%ax = axes('Units', 'inches', 'Position', [2 2 4.5/2 3.5/2]); 
ax.ActivePositionProperty = 'position';

hold on;
% Plot a-magnon
scatter(kx, squeeze(All_grouped_eigenvalues(1,:,:)), 15, 'MarkerEdgeColor', dark_blue, ...
    'MarkerFaceColor', dark_blue, 'Marker', magnon_marker, ...
    'MarkerFaceAlpha', magnon_alpha, 'MarkerEdgeAlpha', magnon_alpha, ...
    'DisplayName', 'a-magnon');

% Plot c-magnon
scatter(kx, squeeze(All_grouped_eigenvalues(2,:,:)), 15, 'MarkerEdgeColor', orange, ...
    'MarkerFaceColor', orange, 'Marker', magnon_marker, ...
    'MarkerFaceAlpha', magnon_alpha, 'MarkerEdgeAlpha', magnon_alpha, ...
    'DisplayName', 'c-magnon');

% Plot a-antimagnon with 90% transparency
scatter(kx, squeeze(All_grouped_eigenvalues(3,:,:)), 15, 'MarkerEdgeColor', dark_blue_anti, ...
    'MarkerFaceColor', antimagnon_face, 'Marker', antimagnon_marker, ...
    'MarkerFaceAlpha', antimagnon_alpha, 'MarkerEdgeAlpha', antimagnon_alpha, ...
    'DisplayName', 'a-antimagnon');


% Plot c-antimagnon with 90% transparency
scatter(kx, squeeze(All_grouped_eigenvalues(4,:,:)), 15, 'MarkerEdgeColor', orange_anti, ...
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

h_dummy_c_magnon = scatter(NaN, NaN, 36, orange, magnon_marker, 'filled', ...
    'MarkerFaceAlpha', magnon_alpha1, 'MarkerEdgeAlpha', magnon_alpha, 'MarkerFaceColor', orange, 'MarkerEdgeColor', orange);

% Dummy plots for antimagnons (transparent markers)


h_dummy_a_antimagnon = scatter(NaN, NaN, 36, dark_blue, antimagnon_marker, 'filled', ...
    'MarkerFaceAlpha', antimagnon_alpha1, 'MarkerEdgeAlpha', antimagnon_alpha1, 'MarkerFaceColor', dark_blue_anti, 'MarkerEdgeColor', dark_blue);

h_dummy_c_antimagnon = scatter(NaN, NaN, 36, orange, antimagnon_marker, 'filled', ...
    'MarkerFaceAlpha', antimagnon_alpha1, 'MarkerEdgeAlpha', antimagnon_alpha1, 'MarkerFaceColor', orange_anti, 'MarkerEdgeColor', orange);

ylim([0 60]);
xlim([-2e8 2e8]);
% Turn off visibility of dummy plots
set([h_dummy_a_magnon,  h_dummy_c_magnon, ...
     h_dummy_a_antimagnon,  h_dummy_c_antimagnon], 'Visible', 'on');

% Customize the legend using dummy handles
legend_handles = [h_dummy_a_magnon,  h_dummy_c_magnon, ...
                  h_dummy_a_antimagnon,  h_dummy_c_antimagnon];
legend_labels = {'a-magnon', 'c-magnon', ...
                'a-antimagnon', 'c-antimagnon'};

legend(legend_handles, legend_labels, ...
    'Location', 'best', ...
    'Box', 'off', ...                      % Remove the black frame
    'FontSize', 28, ...                    % Set desired font size
    'FontName', 'Arial', ...               % Set desired font type
    'FontWeight', 'bold', ...              % Optional: set font weight
    'TextColor', 'black', ...                  % Set text color to black
    'Interpreter', 'latex');               % Optional: use LaTeX interpreter

set(gca, 'LineWidth', 2.5,...
         'FontName', 'Arial', ...
         'FontSize', 28, ...                    % Increased font size for ticks
         'FontWeight', 'bold', ...              % Make ticks bold
         'TickLabelInterpreter', 'latex', ...
         'XAxisLocation', 'bottom', ...
         'YAxisLocation', 'left');      % Use LaTeX interpreter for tick labels
box on;




figure;
%BIG
ax = axes('Units', 'inches', 'Position', [2 2 4.5 3.5]); %left=1, bottom=1, width=4 inches, height=3 inches
%small
%ax = axes('Units', 'inches', 'Position', [2 2 4.5/2 3.5/2]); 
ax.ActivePositionProperty = 'position';
% %%%% Start plotting simulated data
im = imagesc(kxCOM2,fCOM2,(z1COM2+z1COM1+z1COM4+z1COM3)/10);
c = gray; %[linspace(0.5, 1, 256)', linspace(0, 1, 256)', linspace(0.5, 1, 256)'];
c = flipud(c);
colormap(c);
im.AlphaData =1;
caxis([0 1]);
cb = colorbar;           % Adds a colorbar to the current axes
cb.Ticks = [0 0.5 1];    % Sets tick positions at 0, 0.5, and 1
cb.TickLabels = {'0', '0.5', '1'};  % (Optional) Sets custom labels for the ticks
cb.FontName = 'Arial';
cb.TickLabelInterpreter = 'latex';
% cb.FontWeight = 'bold';
cb.FontSize = 28;
cb.Ticks = [0 0.5 1];
cb.TickLabels = {'0', '0.5', '1'};
ylim([0 60]);
xlim([-2e8 2e8]);
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

% Reverse the Y-axis direction
set(gca, 'YDir', 'normal');

hold off;
%%%% End of plotting simulated data

set(gca, 'LineWidth', 2.5,...
         'FontName', 'Arial', ...
         'FontSize', 28, ...                    % Increased font size for ticks
         'FontWeight', 'bold', ...              % Make ticks bold
         'TickLabelInterpreter', 'latex', ...
         'XAxisLocation', 'bottom', ...
         'YAxisLocation', 'left');      % Use LaTeX interpreter for tick labels
box on;