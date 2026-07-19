%% TO DO
%{

1. decompose the code into functions for ali to make GUI
2. integrate omar yehia's code for part 1

%}

%% time axis properties initialization < bit period, sampling frequency, number of bits >

n_bits = 1e5;
Tb = 1; % time period of each individual bit (Tb = 1s)
fs = 100; % sampling frequency (should be >= 10/Tb, increase with caution!!)
t = 1/fs:1/fs:Tb*n_bits; % time axis definition [start:step:end]

%% symbols initialization < insert symbols here >

% nrz symbols (s1(t) = 1; s2(t) = 0) [0 <= t <= Tb]
s1_nrz = ones(1, Tb*fs);
s2_nrz = zeros(1, Tb*fs);

% manchester symbols (s1(t) = [1 -1]; s2(t) = [-1 1])
% [0<=t<=Tb/2 Tb/2<=t<=Tb]
s1_manchester = [ones(1, Tb*fs/2) -1*ones(1, Tb*fs/2)];
s2_manchester = [-1*ones(1, Tb*fs/2) ones(1, Tb*fs/2)];

% symbols containers
s_nrz = [s1_nrz; s2_nrz];
s_man = [s1_manchester; s2_manchester];

% phis containers
phi_nrz = zeros(2, fs*Tb);
s_ij_nrz = zeros(2, 2);
phi_man = zeros(2, fs*Tb);
s_ij_man = zeros(2, 2);

% pre-initialized phis and coefficients for both line codes
[phi_nrz, s_ij_nrz] = GSA(s_nrz, s_ij_nrz, phi_nrz, fs, 2);
[phi_man, s_ij_man] = GSA(s_man, s_ij_man, phi_man, fs, 2);

%% bit stream definition

bits=randi([0 1],1, n_bits); % 01010101011
bstream = repelem(bits,fs); % adding samples to each bit in time domain
manchester_bits = zeros(1, n_bits); % container for manchester bits (-1 instead of zeros)

for i=1:n_bits
    if bits(i) == 0
        manchester_bits(i) = -1;
    else
        manchester_bits(i) = 1;
    end
end

mstream = repelem(manchester_bits, fs); % manchester bits
s_nrz_t = bstream .* repelem(phi_nrz(1, :), n_bits); % nrz line coding
s_man_t = zeros(1, length(mstream)); % manchester line coding

for i=1:fs:fs*n_bits
    s_man_t(i:i+fs-1) = mstream(i:i+fs-1).*phi_man(1, :);
end

% ------------------- PICK THE LINE CODE HERE -------------------
line_codes = ["Unipolar NRZ", "Manchester"];
line_code = line_codes(1); % line code selection index

if line_code == line_codes(1)
    s_t = s_nrz_t;
    phi = phi_nrz;
    s_ij_t = s_ij_nrz;
else
    s_t = s_man_t;
    phi = phi_man;
    s_ij_t = s_ij_man;
end

%{ 
------------------- BLOCK DIAGRAM OF TRANSCEIVER -------------------

----------     -------------    s_t  /-\     s_t * φ(t)     /-\  channel   /-\                 
'  0101  '---->'NRZ or MAN '-------->'X'------------------->'+'----------> 'x' --------> [ ∫ ] -------> [ Decision Device (λ) ] -----> 0101                                 
'--------'     '-----------'         \-/                    \-/            \-/                     
                                      ^                      ^              ^   
                                      |                      |              |  
                                      φ(t)                   w(t)           φ(t)
%}

%% propagation in channel
Eb = 1; % bit energy for both line codes
SNR = -10:5:10; % SNR = Eb/No (in dB) axis

errs = zeros(1, length(SNR)); % container (ignore it)
ber = zeros(1, length(SNR)); % theoritical bit error rate container
ber_noise = zeros(1, length(SNR)); % practical bit error rate container

for k = 1:length(SNR) % sweeping over each SNR 

    N_o = (Eb/(10^(SNR(k)/10)));
    AWGN = normrnd(0,sqrt(N_o*fs/2),1,length(t)); % variance = N_o*fs/2
    x = s_t + AWGN; % signal propagating in channel
    figure;
    plot(t, x)
    
    %% slicing x(t) into symbols
    
    s_x = zeros(n_bits, fs*Tb); % x(t) symbols container
    s_x(1, :) = x(1:fs*Tb);
    
    for i=2:n_bits
        s_x(i, :) = x((i-1)*fs+1:(i-1+Tb)*fs); % decomposing x(t) into symbols
    end
    
    %% find phis at channel output
    
    M = n_bits; % number of symbols 
    N = 2; % number of basis functions (2 for 2d plot, 3 for 3d plot)
    
    % ignore this ( phis and coefficients init )
    s_ij = zeros(M, N);
    phi_j = zeros(N, fs*Tb);
    phi_j(1, :) = phi(1, :);
    s_ij(1, 1) = sum(s_x(1, :) .* phi_j(1, :)) * (1/fs); 
    
    % finding quadrature phase components of x(t)
    for i=2:N
        phi_s = [zeros(1, length(phi_j(1, :)))];
        for j=1:i-1
            s_ij(i, j) = sum(s_x(i, :) .* phi_j(j, :)) * (1/fs);
            phi_s = phi_s + s_ij(i, j)*phi_j(j, :);
        end
        g_i = s_x(i, :) - phi_s;
        if s_x(i, :) == phi_s
            phi_j(i, :) = 0;
        else
            Eg = sum(g_i.^2)*(1/fs);
            s_ij(i, i) = sum(s_x(i, :).*phi_j(i, :))*(1/fs);
            disp(Eg)
            phi_j(i, :) = (g_i / sqrt(Eg));
        end
    end
    
    for i=1:M 
        for j=1:N
            if (i == j)
                break
            else
            s_ij(i,j) = sum(s_x(i, :) .* phi_j(j, :)) * (1/fs);
            end
        end
    end
    
    ber(k) = qfunc((s_ij_t(1, 1) - s_ij_t(2, 1)) / sqrt(2 * N_o)); % theoritical BER calculation using Q-function

    %% receiver model

    Lamda=(s_ij_t(1, 1) + s_ij_t(2, 1))/2; % decision device

    out_bits=zeros(1,n_bits); % container for output data

    for i = 1:n_bits
        out_bits(i)=Lamda<=s_ij(i,1); % applying decision device to decoded data
    end 

    Error_in_bits=sum(abs(out_bits-bits)); % number of wrong bits at output
    errs(k) = Error_in_bits;

    ber_noise(k)=Error_in_bits/n_bits; % practical BER calculation (number of wrong bits / total number of bits)

    %% plotting constellation diagram
    s_ij_T = transpose(s_ij); % transposing so each row is an axis for each phi
    part = 2; % ignore this
    phi_1_ax = s_ij_T(1, :); % phi_1 axis
    phi_2_ax = s_ij_T(2, :); % phi_2 axis
    
    if ((length(s_ij_T(1, :)) < 3 || part == 2))     % 2d constellation diagram
        figure;
        scatter(phi_1_ax, phi_2_ax);
        xlabel('φ1')
        ylabel('φ2')

        if Lamda == 0.5
            xticks([-3 -2 -1 0 0.5 1 2 3])
            xticklabels({'-3' ,'-2', '-1', '0', sprintf('λ=0.5'), '1', '2', '3'})
            xline(0.5, 'LineWidth', 1, 'Color', 'black')
        else 
            xticks([-3 -2 -1 0 1 2 3])
            xticklabels({'-3' ,'-2', '-1', sprintf('λ=0'), '1', '2', '3'})
            xline(0, 'LineWidth', 1, 'Color','black')
        end
        xlim ([-3 3]);
        ylim ([-3 3]);
        hold on;
        plot(s_ij_t(1, 1), 0, '+', 'LineWidth', 2, 'Color', 'red', "MarkerSize", 10);
        plot(s_ij_t(2, 1), 0, '+', 'LineWidth', 2, 'Color', 'red', "MarkerSize", 10);
        grid on;
        title(sprintf('Channel Output @ SNR = %d for %s', SNR(k), line_code));
        
    else 
        if (part == 1)                          % 3d constellation diagram
            phi_3_ax = s_ij_T(3, :);
            figure();
            scatter3(phi_1_ax, phi_2_ax, phi_3_ax);
            hold on
            plot3(1, 0, 0, 'ro', 'MarkerSize', 10);
            plot3(0, 0, 0, 'ro', 'MarkerSize', 10);
            xlabel('φ1')
            ylabel('φ2')
            zlabel('φ3')
            title(sprintf('Channel Output @ SNR = %d', SNR))
            
            xlim([-1.1*max(max(s_ij)) 1.1*max(max(s_ij))])
            ylim([-1.1*max(max(s_ij)) 1.1*max(max(s_ij))])
            zlim([-1.1*max(max(s_ij)) 1.1*max(max(s_ij))])
            grid on
        
            %{
            for i=1:length(phi_1_ax)
                text(phi_1_ax(i), phi_2_ax(i), phi_3_ax(i), sprintf(' s%d', i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);
            end
            %}
        end
    end
end

%% BER plot
figure;
semilogy(SNR, ber); % theoritical
hold on;
semilogy(SNR, ber_noise); % practical
xlabel('Eb/No (dB)')
ylabel('BER')
legend("Theoritical", "Practical")
title(sprintf('Bit Error Rate of %s', line_code))