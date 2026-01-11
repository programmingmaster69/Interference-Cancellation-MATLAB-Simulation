clc
clear all

%This simulation will simulate communication without defense for 1/3rd part of
%trials, communication with IC for next 1/3rd part, and communication with
%IC and SSE for last part. So make sure num_trials is a multiple of 3
%Reactive jammer is used in all three cases.

%==== Defining parameters =====%
N = 10^4; % Total symbols per transmission
Num_scs = 64;
Occupied_scs = 48;
Num_pckts = 500;
N_sym_pckt = 14;
N_fft = Occupied_scs;
cp_len = 16;
num_trails = 30;
Mod_ord = 2;
N_bits_per_symb = Occupied_scs; % for BPSK
SNR_dB = 16;

%Pilot Info
pilot_pos = [1,2,4,6,8,10,12];
pilot_number = 0;
pilot_matrix = ones(Occupied_scs,1);
pilot_matrix(Occupied_scs/2 - 6: Occupied_scs/2 + 6) = -1;

%==== MIMO parameters ====%

n_tx = 2;
n_rx = 2;

%Sender signal Enchancement parameter
rotation_vector = [1;0];
Apply_SSE = false; % Set True if SSE is required

%Threshold of when to apply IC (Applying IC at low JSR might
%actually worsen the signal recovery, since ratio estimation becomes
%inaccurate)

threshold = 0.4; 
Apply_IC = false; % Set True if IC is required

% == Jammer Info == %
jammer_type = 'Reactive';
jamming_power = zeros(1,11);
for i = 1: length(jamming_power)
    jamming_power(i) = i-1; % Jamming power 0 - 10 W.
end
% PDR Results
pdr_results = zeros(num_trails,length(jamming_power));

%Array Initializations

equivalent_hs_freq_est = zeros(n_rx,1);

for iter = 1:num_trails
    if iter == 11
        Apply_IC = true;
    elseif iter == 21
        Apply_SSE = true;
    end
    for i = 1: length(jamming_power)
    
        ip_data = rand(1,N_bits_per_symb*N) > 0.5;
        ip_data = double(ip_data);
        ip_symbols = pskmod(ip_data,Mod_ord);
        parallel_data = reshape(ip_symbols, N_bits_per_symb,[]);
        
        Packet_success = 0;
        for p = 1:Num_pckts

        % == Defining channel per each coherence interval (each packet) ==%

            h_s = (randn(n_rx,n_tx) + 1j *randn(n_rx,n_tx))/sqrt(2);
            h_j =  (randn(n_rx,1) + 1j *randn(n_rx,1))/sqrt(2);
            equivalent_hs = h_s*rotation_vector;
            equivalent_hs = equivalent_hs / norm(equivalent_hs);
            if strcmp(jammer_type, 'Reactive') 
                    jammer_signal_after_detection = sqrt(jamming_power(i)/2)*(randn(1,Num_scs) + 1i*randn(1,Num_scs));
            end 
            if p ~= 1
                h_j = h_s(:,1) + 0.4*(randn(n_rx,1) + 1i*randn(n_rx,1));
                h_j = h_j / sqrt(mean(abs(h_j).^2));
            end
            Hj = zeros(n_rx,N_fft);
            for j = 1:n_rx
                Hj(j,:) = fft(h_j(j), N_fft);
            end
            symbol_success = 0;
            pilot_number = 0;
            for s = 1:N_sym_pckt  
            % === OFDM transmission starts here === %
                % Reactive jammer is in sensing mode for the first OFDM
                % symbol, after detection, it starts jamming.
                if s == 1 || s == 2
                    jammer_signal = zeros(1, Num_scs);
                else
                    jammer_signal = jammer_signal_after_detection;
                end

                if ismember(s,pilot_pos)
                    parallel_data(:,(p-1)*(N_sym_pckt)+ s) = pilot_matrix;
                end

                time_data = sqrt(N_fft)*ifft(parallel_data(:,(p-1)*(N_sym_pckt)+ s), N_fft);
                serial_data = time_data.';
                serial_tx_data = [serial_data(end-cp_len+1:end) serial_data]; % after cp addition

            % === Computing the received signal at rx antennas === %

                n = 10^(-SNR_dB/20)*sqrt(1/2)*(randn(n_rx,Num_scs) + 1j *randn(n_rx,Num_scs));
                rx_raw_data = equivalent_hs * serial_tx_data + h_j * jammer_signal + n;
                rx_serial_data = zeros(n_rx, Occupied_scs);
                rx_fft = rx_serial_data.';
                for j = 1:n_rx
                     rx_serial_data(j,:) = rx_raw_data(j,cp_len+1:end); %Remove CP
                     rx_fft(:,j) = sqrt(1/N_fft)*fft(rx_serial_data(j,:).',N_fft);
                end
                y = rx_fft;
    
            % === Alpha estimation === %
               equivalent_hs_freq = repmat(equivalent_hs,1,Occupied_scs);
                if ismember(s, pilot_pos)
                    alpha_m = zeros(N_bits_per_symb,1);
                    y1_pilot = y(:,1);  % Rx antenna 1
                    y2_pilot = y(:,2);
                    
                    if pilot_number == 0
                        equivalent_hs_est(1) = mean(y1_pilot./pilot_matrix);
                        equivalent_hs_est(2) = mean(y2_pilot./pilot_matrix);
                        alpha_avg = 0;
                        JSR=0;
                    elseif mod(pilot_number,2) == 1
                        expected_1 = equivalent_hs_freq_est(1,:).'.*pilot_matrix;
                        expected_2 = equivalent_hs_freq_est(2,:).'.*pilot_matrix;

                        residual_1 = y1_pilot - expected_1; %Residual = Jammer + Noise
                        residual_2 = y2_pilot - expected_2;
                        % Estimate total residual power = jammer + noise
                        P_residual = (mean(abs(residual_1).^2) + mean(abs(residual_2).^2))/2;

                        % Estimate useful signal power
                        P_signal = (mean(abs(expected_1).^2) + mean(abs(expected_2).^2)) / 2;

                        JSR = P_residual / P_signal;
                        if JSR > threshold && Apply_IC == true
                            
                            alpha = residual_1./residual_2;
                            alpha_avg = mean(alpha);
                            hj_orthogonal = [1,-alpha_avg].';
                            if Apply_SSE == true
                                rotation_vector = (h_s\hj_orthogonal)/norm(h_s\hj_orthogonal);
                            end
                        else 
                            alpha = zeros(N_bits_per_symb,1);
                            alpha_avg = 0;
                            rotation_vector = [1;0];
                        end
                    elseif mod(pilot_number,2) == 0
                        if mod(pilot_number,4) ~= 0 
                            if alpha_avg ~= 0
                                equivalent_hs_est(2) = equivalent_hs_est(1)/alpha_avg - mean((y1_pilot/alpha_avg - y2_pilot)./pilot_matrix);
                            else
                                equivalent_hs_est(2) = 0;
                            end
                        else
                            equivalent_hs_est(1) = alpha_avg*equivalent_hs_est(2) + mean((y1_pilot - alpha_avg*y2_pilot)./pilot_matrix);
                        end
                    end
                    pilot_number = pilot_number + 1;
                    true_alpha = mean(Hj(1,:) ./ Hj(2,:)); % Compare with true alpha (if known)
                    equivalent_hs_freq_est = [fft(equivalent_hs_est(1),N_fft);fft(equivalent_hs_est(2),N_fft)];
                end
                % disp([equivalent_hs_freq(1,1) equivalent_hs_freq_est(1,1)]);
                % disp([equivalent_hs_freq(2,1) equivalent_hs_freq_est(2,1)]);
                % disp([alpha_avg,true_alpha,JSR]);
                    
            % === Interference Cancellation === %

                rx_dec_data = zeros(N_bits_per_symb,1);
                for sub = 1 : Occupied_scs
                    rx_dec_data(sub) = (y(sub,1) - alpha_avg*y(sub,2))/(equivalent_hs_freq_est(1,sub) - alpha_avg*equivalent_hs_freq_est(2,sub));
                end

            % === Signal Demodulation and comparison === %

                dec_bits = pskdemod(rx_dec_data,Mod_ord);
                parallel_bits = pskdemod(parallel_data(:,(p-1)*(N_sym_pckt)+ s),Mod_ord);
                if all(dec_bits == parallel_bits)
                    symbol_success = symbol_success + 1; % Count successful symbols
                end
            end
            if symbol_success == N_sym_pckt
                Packet_success = Packet_success + 1;
            end
        end
       
        pdr_results(iter,i) = Packet_success/Num_pckts*100;
        fprintf("Trial %2d → PDR: %.2f%%\n", iter, pdr_results(iter,i));
    end
    iter %#ok<*NOPTS>
   
end

final_pdr_without_IC = mean(pdr_results(1:num_trails/3,:),1);
final_pdr_with_IC = mean(pdr_results(num_trails/3+1: 2*num_trails/3,:),1);
final_pdr_with_IC_and_SSE = mean(pdr_results(2*num_trails/3 : end,:),1);

% === Plotting PDR vs. Jamming Power ===

figure;
plot(jamming_power, final_pdr_without_IC, '-o', 'LineWidth', 2);
hold on
plot(jamming_power, final_pdr_with_IC, '-o', 'LineWidth', 2);
plot(jamming_power, final_pdr_with_IC_and_SSE, '-o', 'LineWidth', 2);
grid on;
xlabel('Jamming power');
ylabel('Packet Delivery Rate (PDR) [%]');
title('PDR vs. Jamming Power');
ylim([0 110]);
legend('Defenseless','With IC' , 'With IC and SSE')
hold off