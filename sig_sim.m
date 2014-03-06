%Script to simulate events in proton detector for darklight

%first define some constants:
tq = 1e-10; %time quanta, 100pS

%event rate information (seconds)
rate_sep_mean = 1/5e6;
rate_sep_sd = rate_sep_mean/3;

%event magnitude information
mag_mean = 1;
mag_sd = 0.4;

%Sampeling information (seconds)
t_sample = 1/350e6;

%simulation time (in seconds)
t_simulation = 10000e-9;

%% First: generate simulation space
num_tq = t_simulation/tq;
if num_tq < 100
    error('sig_sim:time_scale_err', 'Check time scale')
end
%pre-allocate all simulation spaces (if additional simulation spaces are
%added, pre-allocate them here.
num_simulation_steps = 3;
sim_space = zeros(num_simulation_steps, num_tq);

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

%% Step 2: generate event time magnitudes
sim_space(2,:) = normrnd(mag_mean, mag_sd, [1, num_tq]) .* sim_space(1,:);

%% Step 3: generate shaping
%first pass at shaping alg.
shape_coeff = [0, .1, .4, .5, .6, .7, .8, .9, 1, .6, .5, .4, .3, .15, .05];
for index = (1:length(shape_coeff))
    if index == 1
        sim_space(3,:) = sim_space(2,:)*shape_coeff(1);
    else
        sim_space(3,:) = sim_space(3,:) + horzcat(zeros(1,index-1),sim_space(2,(1:num_tq-(index-1))))*shape_coeff(index);
    end
end
plot(sim_space(3,:),'o-')
%TODO add gaussian shape_coeff depending on the shaping time, use 3 sigma
%shaping time on the early edge, and 6 sigma on the falling edge

%% Step 4: generate sampling times
tq_sample = floor(t_sample / tq);
if tq_sample < 10
    error('sig_sim:sample_rate','Sample rate too high for time scale')
end
num_samples = num_tq / tq_sample;
sampled_signal = sim_space(3,(1:tq_sample:num_tq));
sample_time = (1:tq_sample:num_tq);
plot(sim_space(3,:),'-')
hold on;
plot(sample_time,sampled_signal, 'r*')
hold off;