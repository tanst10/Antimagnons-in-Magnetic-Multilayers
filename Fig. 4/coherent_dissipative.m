%% Figure 4

clc; clear all; close all

%  Definition of parameters
kx = 0.68e8;
ky = 0;

d1 = 10e-9;
d2 = 10e-9;
Ms1 = 1.4e5;
Ms2 = 7.4e5;
Aex1 = 3.7e-12;
Aex2 = 8.7e-12;
u0 = 4*pi*1e-7;
gamma = 28e9*2*pi;



H1 = linspace(-1e6,0.5e6,25);%
H2 = 0e5 ;% or -8.5e5 for Fig. 4b





%number of cells
N=1;

angluar_f_obc=zeros(N*4,length(H1),length(H2));
eigenvector = zeros(N*4+1,N*4,length(H1),length(H2));


% Initialize storage for frequencies
a_magnon = NaN(4*N, length(H1),length(H2));
c_magnon = NaN(4*N, length(H1),length(H2));
a_antimagnon = NaN(4*N, length(H1),length(H2));
c_antimagnon = NaN(4*N, length(H1),length(H2));

eigvals_pos = NaN(3, length(H1),length(H2));

% Define colors for each band
dark_blue = [68, 113, 196]/255;
dark_blue_anti = [161, 184, 225]/255;
light_blue = [0, 175, 239]/255;
light_blue_anti = [178, 230, 250]/255;
orange = [236, 124, 48]/255;
orange_anti = [249, 222, 189]/255;


fig1 = figure('Visible','off');       
ax1 = axes('Parent', fig1, 'Units', 'inches', 'Position', [2 2 4.5 3.5]);
ax1.ActivePositionProperty = 'position';
hold(ax1,'on');


fig2 = figure('Visible','off');       
ax2 = axes('Parent', fig2, 'Units', 'inches', 'Position', [2 2 4.5 3.5]);  
ax2.ActivePositionProperty = 'position';
hold(ax2,'on');


axes_list = [ax1, ax2];
for ax = axes_list
    set(ax, ...
        'LineWidth', 2.5, ...
        'FontName', 'Arial', ...
        'FontSize', 26, ...
        'FontWeight', 'bold', ...
        'TickLabelInterpreter', 'latex', ...
        'XAxisLocation', 'bottom', ...
        'YAxisLocation', 'left');
    box(ax, 'on');
end



axis(ax1,'equal');            
ylabel(ax1, '$f/\mathrm{GHz}$', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial', ...
    'FontWeight', 'bold', ...
    'FontSize', 26);

xlabel(ax1, '$$B_a^y/\mu_0/5\times10^4 \mathrm{A/m}$$', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial', ...
    'FontWeight', 'bold', ...
    'FontSize', 26);



axis(ax2,'equal');
xlabel(ax2, '$H_1$', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial', ...
    'FontWeight', 'bold', ...
    'FontSize', 26);

ylabel(ax2, '$H_2$', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial', ...
    'FontWeight', 'bold', ...
    'FontSize', 26);

title(ax2, 'Unnormalized Polarization Ellipse', ...
    'Interpreter', 'latex', ...
    'FontName', 'Arial', ...
    'FontWeight', 'bold', ...
    'FontSize', 26);

ylim(ax1, [0 30]);   
xlim(ax1, [min(H1)/5e4 max(H1)/5e4]);  
xlim(ax2, 'tight');   
ylim(ax2, [min(H1)/5e4 max(H1)/5e4]);   



%Matrix
for m = 1:length(H1)
    for n = 1:length(H2)

        N1 = (1-exp(-abs(kx)*d1))./(abs(kx)*d1);
        N2 = (1-exp(-abs(kx)*d2))./(abs(kx)*d2);

        f2 = Ms1*d1*N1.*N2/2;
        f1 = Ms2*d2*N1.*N2/2;

        w1P = (2*Aex1*gamma*abs(kx)^2/Ms1 + gamma*u0*H1(m)) + gamma*u0*Ms1/2*((1-N1)*(kx)^2/abs(kx)^2+N1);
        w1N = -(2*Aex1*gamma*abs(kx)^2/Ms1 + gamma*u0*H1(m)) - gamma*u0*Ms1/2*((1-N1)*(kx)^2/abs(kx)^2+N1);
        w2P = (2*Aex2*gamma*abs(kx)^2/Ms2 + gamma*u0*H2(n)) + gamma*u0*Ms2/2*((1-N2)*(kx)^2/abs(kx)^2+N2);
        w2N = -(2*Aex2*gamma*abs(kx)^2/Ms2 + gamma*u0*H2(n)) - gamma*u0*Ms2/2*((1-N2)*(kx)^2/abs(kx)^2+N2);


        D1P = 1/2*gamma*u0.*f1.* ((abs(kx)+(kx)^2/abs(kx))+2*kx);
        D2P = 1/2*gamma*u0.*f2.* ((abs(kx)+(kx)^2/abs(kx))+2*kx); 


        D3P = 1/2*gamma*u0.*f2.* ((abs(kx)+(kx)^2/abs(kx))-2*kx);
        D4P = 1/2*gamma*u0.*f1.* ((abs(kx)+(kx)^2/abs(kx))-2*kx);

        D1N = 1/2*gamma*u0.*f1.* (-(abs(kx)+(kx)^2/abs(kx))+2*kx);
        D2N = 1/2*gamma*u0.*f2.* (-(abs(kx)+(kx)^2/abs(kx))+2*kx);

        D3N = 1/2*gamma*u0.*f2.* (-(abs(kx)+(kx)^2/abs(kx))-2*kx);
        D4N = 1/2*gamma*u0.*f1.* (-(abs(kx)+(kx)^2/abs(kx))-2*kx);

        delta1P_ = 0;%1/2*gamma*u0*(f1*(abs(kx)-(kx)^2/abs(kx))+Hex1); %interlayer
        delta2P_ = 0;%1/2*gamma*u0*(f2*(abs(kx)-(kx)^2/abs(kx))+Hex2);
        delta1N_ = -delta1P_;
        delta2N_ = -delta2P_;

        DELTA1P = gamma*u0*Ms1/2*((1-N1)*(kx/abs(kx))^2-N1); %intralayer
        DELTA2P = gamma*u0*Ms2/2*((1-N2)*(kx/abs(kx))^2-N2);%
        DELTA1N = -DELTA1P;
        DELTA2N = -DELTA2P;
    
        wP = [w1P;w2P];
        wN = [w1N;w2N];

        DelP_ = [delta1P_;delta2P_];
        DelN_ = [delta1N_;delta2N_];

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
                H_obc(i,2*N+2) = DelP_(nrest);
                %H_obc(i,4*N) = DelP_(nrest);%pbc
            end

            if (i<2*N) && (i>1) && (mod(i,2)==0) %if i=even
                H_obc(i,i) = wP(nrest);
                H_obc(i,i-1) = DP_intra(nrest);
                H_obc(i,i+1) = DP_INTER(num); % row 2 is D3P at i+1 col
                H_obc(i,2*N+i) = DELP(nrest);
                H_obc(i,2*N+i-1) = DelP_(nrest);
                H_obc(i,2*N+i+1) = DelP_(nrest);           
            end

            if (i<2*N) && (i>1) && (mod(i,2)==1) %if i=odd
                H_obc(i,i) = wP(nrest); 
                H_obc(i,i-1) = DP_INTER(num); % row 3 is D4P at i-1 col
                H_obc(i,i+1) = DP_intra(nrest);
                H_obc(i,2*N+i) = DELP(nrest);
                H_obc(i,2*N+i-1) = DelP_(nrest);
                H_obc(i,2*N+i+1) = DelP_(nrest);  
            end

            if i==2*N
                H_obc(i,i) = wP(nrest);
                % H_obc(i,1)=DP_INTER(num);%pbc i is even so num =1 D3P
                H_obc(i,i-1) = DP_intra(nrest);
                H_obc(i,2*N+i) = DELP(nrest);
                H_obc(i,2*N+i-1) = DelP_(nrest); 
                % H_obc(i,2*N+1) = DelP_(nrest);%pbc
            end
    
            if i==2*N+1
                H_obc(i,i) = wN(nrest);
                H_obc(i,i+1) = DN_intra(nrest);
                % H_obc(i,4*N) = DN_INTER(num);%pbc
                H_obc(i,i-2*N) = DELN(nrest);
                H_obc(i,i-2*N+1) = DelN_(nrest);
                % H_obc(i,2*N) = DelN_(nrest);%pbc
            end
    
            if (i>2*N+1) && (i<4*N) && (mod(i,2)==0) %if i=even
                H_obc(i,i) = wN(nrest);
                H_obc(i,i-1) = DN_intra(nrest);
                H_obc(i,i+1) = DN_INTER(num);
                H_obc(i,i-2*N) = DELN(nrest);
                H_obc(i,i-2*N-1) = DelN_(nrest);
                H_obc(i,i-2*N+1) = DelN_(nrest);
            end
    
            if (i>2*N+1) && (i<4*N) && (mod(i,2)==1) %if i=odd
                H_obc(i,i) = wN(nrest);
                H_obc(i,i-1) = DN_INTER(num);
                H_obc(i,i+1) = DN_intra(nrest);
                H_obc(i,i-2*N) = DELN(nrest);
                H_obc(i,i-2*N-1) = DelN_(nrest);
                H_obc(i,i-2*N+1) = DelN_(nrest);
            end
    
            if i==4*N
                H_obc(i,i) = wN(nrest);
                H_obc(i,i-1) = DN_intra(nrest);
                % H_obc(i,2*N+1) = DN_INTER(num); %pbc
                H_obc(i,i-2*N) = DELN(nrest); 
                H_obc(i,i-2*N-1) = DelN_(nrest);
                % H_obc(i,1) = DelN_(nrest); %pbc
            end            
        end
        H_all_obc(:,:,m,n) = H_obc;

        % Compute eigenvalues and eigenvectors
        [V, D] = eig(H_obc); % Compute eigenvalues and eigenvectors
        eigvals = diag(D);   % Extract eigenvalues into a vector
        tol = 1e-1;
        for x = 1:4*N
            if abs(imag(eigvals(x)/2/pi/1e9)) > tol
                eigvals(x)= NaN;
            end
        end
        angluar_f_obc(:,m,n) = eigvals/2/pi/1e9;
        eigenvector (1,:,m,n) = eigvals/2/pi/1e9;
        eigenvector(2:4*N+1,:,m,n) = V;

        keep         = eigvals >= 0;       
        posvals = eigvals(keep)/2/pi/1e9;
        len_posvals = numel(posvals);   
        eigvals_pos(1:len_posvals,m,n)  = posvals;       
        V_pos        = V(:, keep);         
        m_plus = V_pos(1:2*N,:);
        m_minus = V_pos(2*N+1:4*N,:);
        mx = (m_plus + m_minus)/(sqrt(2));
        my = (m_plus - m_minus)/(sqrt(2)*1i);

        % Define a time array over one full period (0 to 2*pi) with sufficient points
        ntime = 500;
        t = linspace(0, 2*pi, ntime);    % time vector (radians)
        scale = reshape( exp(-1i * t), 1, 1, [] );
        
        % Compute the real time-domain field components using Euler's relation
        mx_t = real(mx .* scale);
        my_t = real(my .* scale);

        mean_mx = mean(mx_t, 3);   %1Darray size: ntime
        mean_my = mean(my_t, 3);

        mx_t = mx_t - mean_mx;
        my_t = my_t - mean_my;
        
        for num = 1:len_posvals
             crossval_a = mx_t(1,num,1)*my_t(1,num,2) - my_t(1,num,1)*mx_t(1,num,2);
             crossval_b = mx_t(2,num,1)*my_t(2,num,2) - my_t(2,num,1)*mx_t(2,num,2);

             if crossval_a < 0 %  m_plus(1,1) > m_minus(1,1)
                color_a = dark_blue;   % Dark blue for right-handed (CW rotation)
             else
                color_a = dark_blue_anti; % Light blue for left-handed (CCW rotation)
             end

             if crossval_b < 0 %  m_plus(2,1) > m_minus(2,1)
                color_b = orange;   % Dark blue for right-handed (CW rotation)
             else
                color_b = orange_anti; % Light blue for left-handed (CCW rotation)
             end

            amp = sqrt(mx_t.^2 + my_t.^2);
            norm = max(amp); % Normalize the ellipse so that its maximum amplitude is 1
            norm = 1; %don't normalize
            freq_offset  = reshape(posvals, 1, [], 1);   % 1 × 3 × 1
            mx_t0 = mx_t / norm + H1(m)/5e4;
            my_t0 = my_t / norm + freq_offset;



            plot(ax1, squeeze(mx_t0(1,num,:)), squeeze(my_t0(1,num,:)), ...
                   '-',  'Color', color_a, 'LineWidth',3);

        
            plot(ax1, squeeze(mx_t0(2,num,:)), squeeze(my_t0(2,num,:)), ...
                   '-', 'Color', color_b, 'LineWidth',3);


        end


        % Tolerance to check if eigenvalue is real
        tol = 1e-3*(2 * pi * 1e9); % unit is GHz 
        
        % Preallocate storage for 4 groups (mod 4 = 1, 2, 3, 0)
        grouped_eigenvalues = NaN(4, size(H_obc, 1)); % 4 rows for the 4 groups
        grouped_eigenvectors = NaN(4, size(H_obc, 1), size(H_obc, 2)); % Store eigenvectors in 4 groups
        
        % Process each eigenvalue and its corresponding eigenvector
        for eig_idx = 1:length(eigvals)
            current_eigval = eigvals(eig_idx);        % Current eigenvalue
            current_eigvec = V(:, eig_idx);     % Corresponding eigenvector
        
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
        All_grouped_eigenvalues(:, :, m, n) = grouped_eigenvalues;
        All_grouped_eigenvectors(:, :, :, m, n) = grouped_eigenvectors;
        
    end
end 

%%%%%%%%Finish Band Structure Calculation

set([fig1],'Visible','on');   
drawnow;                          

hold off;
%%% End of plotting simulated data


