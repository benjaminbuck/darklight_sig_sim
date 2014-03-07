%Script to simulate events in proton detector for darklight

close all;

%first define some constants:
tq = 100e-12; %time quanta in seconds

%event rate information (seconds)
rate_sep_mean = 1/5e6;
rate_sep_sd = rate_sep_mean/3;

%event magnitude information
mag_mean = 1;
mag_sd = 0.4;

%Sampeling information (seconds)
t_sample = 1/150e6;
sample_bits = 10;

%simulation time (in seconds)
t_simulation = 10000e-9;

%gaussian shaping properties
t_shape = 20e-9;%shaping time in seconds

%plots
gen_all_plots = 0;
gen_sig_gen_plots = 1;
gen_sig_shaping_plots = 1;
gen_dig_plots = 1;

%% First: generate simulation space
num_tq = ceil(t_simulation/tq);
if num_tq < 100
    error('sig_sim:time_scale_err', 'Check time scale')
end
%pre-allocate all simulation spaces (if additional simulation spaces are
%added, pre-allocate them here.
num_simulation_steps = 3;
sim_space = zeros(num_simulation_steps, num_tq);

tq_sample = ceil(t_sample / tq);
if tq_sample < 10
    error('sig_sim:sample_rate','Sample rate too high for time scale')
end
num_samples = ceil(num_tq / tq_sample);
num_digital_steps = 5;
dig_space = zeros(num_digital_steps, num_samples);

%% Step 1: generate event time locations
%We do this by saying that time quanta 0 is an event and then
%generating a gaussian amount of time with an average of
%rate_sep_mean and a standard deviation of rate_sep_sd.
rate_sep_mean_q = rate_sep_mean / tq;
rate_sep_sd_q = rate_sep_sd / tq;

if rate_sep_mean_q < 100
    error('sig_sim:rate_sep_err', 'Rate separation too small')
end

%this is ugly... is there a more creative way to do this?
index = 1;
while index < num_tq
    sim_space(1,index) = 1;
    index = index + floor(normrnd(rate_sep_mean_q, rate_sep_sd_q));
end

if (gen_all_plots == 1) || (gen_sig_gen_plots == 1)
    figure();
    plot(sim_space(1,:),'-o')
    hold on;
    title('First Sim Space: impulses at event locations');
end

%% Step 2: generate event time magnitudes
sim_space(2,:) = normrnd(mag_mean, mag_sd, [1, num_tq]) .* sim_space(1,:);

if (gen_all_plots == 1) || (gen_sig_gen_plots == 1)
    figure();
    plot(sim_space(2,:),'-o')
    hold on;
    title('Second Sim Space: impulses at event locations w gaussian magnitude');
end

%% Step 3: generate shaping
%first pass at shaping alg.

%reminder of gaussian : g(x) = (1/(sd*sqrt(2*pi)))*exp((-1/2)*((x-mean)/sd)^2)
%we say that the shaping time is equal to 6 sigma.  So:
shape_sd = t_shape / tq / 6;
%we put the mean 3 sigma in to give a 3 sigma leadup and a 6 sigma tail
shape_mean = shape_sd*3;
%now we generate the shape coeff
shape_coeff = (1/(shape_sd*sqrt(2*pi)))*exp((-1/2)*(((1:9*shape_sd)-shape_mean)/shape_sd).^2);
%now we normalize the shape_coeff
shape_coeff = shape_coeff / (1/(shape_sd*sqrt(2*pi)));

if (gen_all_plots == 1) || (gen_sig_gen_plots == 1)
    figure();
    plot(shape_coeff);
    hold on;
    title('Signal Shape');
end

for index = (1:length(shape_coeff))
    if index == 1
        sim_space(3,:) = sim_space(2,:)*shape_coeff(1);
    else
        sim_space(3,:) = sim_space(3,:) + horzcat(zeros(1,index-1),sim_space(2,(1:num_tq-(index-1))))*shape_coeff(index);
    end
end

if (gen_all_plots == 1) || (gen_sig_gen_plots == 1)
    figure();
    plot(sim_space(3,:),'o-')
    hold on;
    title('Third Sim Space: shaped events')
end
%TODO add gaussian shape_coeff depending on the shaping time, use 3 sigma
%shaping time on the early edge, and 6 sigma on the falling edge

%% Step 4: generate sampling times

dig_space(2,:) = sim_space(3,(1:tq_sample:num_tq));
dig_space(1,:) = (1:tq_sample:num_tq);
if (gen_all_plots == 1) || (gen_dig_plots == 1)
    figure();
    plot(dig_space(2,:),'-o')
    hold on;
    title('Dig Space 2: Sampled Signal');
    figure();
    plot(sim_space(3,:),'-o')
    hold on;
    plot(dig_space(1,:),dig_space(2,:), 'r*', 'MarkerSize', 10)
    title('Samples overlayed on shaped signal');
end

%% Step 5: quantize digital samples
%max value calculation
max_val = 2^sample_bits-1;
max_analog_val = mag_sd*3;

%insert clipping function here.
dig_space(3,:) = dig_space(2,:); %placeholder for clipped values

%perform quantization, using floor
dig_space(4,:) = floor(dig_space(3,:)/max_analog_val*max_val);

%quantization error:
dig_space(5,:) = (dig_space(4,:)-dig_space(2,:)/max_analog_val*max_val);

if (gen_all_plots == 1) || (gen_dig_plots == 1)
    figure();
    plot(dig_space(4,:), '-o');
    hold on;
    title('Dig space 4: Quantized sampled data');
    
    figure();
    plot(dig_space(2,:), '-o');
    hold on;
    plot(dig_space(4,:)/max_val*max_analog_val,'-r*', 'MarkerSize', 7);
    title('Quantized sampled data, normalized by max value, on the non quantized data');
    
    figure();
    plot(dig_space(5,:), '-o');
    hold on;
    title('Dig space 5: Quantization error in LSBs');
end