%% General Script Parameters
%Script to simulate events in proton detector for darklight
clear;
close all;
tic;

%first define some constants:
tq = 100e-12; %time quanta in seconds

%event rate information (seconds)
rate_sep_mean = 1/100e3;

%event magnitude information
mag_mean = 1;
mag_sd = 0.4;
%RMS = sqrt(mean^2+sd^2)
rms_signal = sqrt(mag_mean^2+mag_sd^2);
%!!!!TODO CHECK THIS DIST...

%Sampeling information (seconds)
t_sample = 1/100e6;
sample_bits = 10;

sample_by_quanta = t_sample / tq;
if sample_by_quanta < 10
    error('gen_prop:quanta_size_err', 'Check quanta size')
end

%Noise information
noise_det_rms = rms_signal / 10000; %~80dB noise floor

noise_shaper_rms = rms_signal / 10000;

%simulation time (in seconds)
t_simulation = 1000e-6;

%gaussian shaping properties
t_shape = 40e-9;%shaping time in seconds

%plots
gen_all_plots = 0;
gen_sig_gen_plots = 1;
gen_sig_shaping_plots = 1;
gen_dig_plots = 1;
gen_hit_find_plots = 1;

fprintf('gen parameters duration: %d\n',toc);

%% First: generate simulation space
tic;
num_tq = ceil(t_simulation/tq);
if num_tq < 100
    error('sig_sim:time_scale_err', 'Check time scale')
end
%pre-allocate all simulation spaces (if additional simulation spaces are
%added, pre-allocate them here.
num_simulation_steps = 4;
sim_space = zeros(num_simulation_steps, num_tq);

tq_sample = ceil(t_sample / tq);
if tq_sample < 10
    error('sig_sim:sample_rate','Sample rate too high for time scale')
end
num_samples = ceil(num_tq / tq_sample);
num_digital_steps = 6;
dig_space = zeros(num_digital_steps, num_samples);

fprintf('sample space generation duration: %d\n',toc);

%% Step 1: generate event time locations
%We do this by saying that time quanta 0 is an event and then
%generating a gaussian amount of time with an average of
%rate_sep_mean and a standard deviation of rate_sep_sd.
tic;
rate_sep_mean_q = rate_sep_mean / tq;

if rate_sep_mean_q < 100
    error('sig_sim:rate_sep_err', 'Rate separation too small')
end

%this is ugly... is there a more creative way to do this?
index = 1;
event_space_tracker = [];
event_space = 0;
while index < num_tq
    sim_space(1,index) = 1;
    event_space = floor(exprnd(rate_sep_mean_q));
    event_space_tracker = horzcat(event_space_tracker, event_space);
    index = index + event_space;
end

if (gen_all_plots == 1) || (gen_sig_gen_plots == 1)
    figure();
    plot(sim_space(1,:),'-o')
    hold on;
    title('First Sim Space: impulses at event locations');
    figure();
    hist(event_space_tracker,30);
    hold on;
    title('Histogram of time between event locations');
end

num_input_pulses = sum(sim_space(1,:));

fprintf('Step 1 duration: %d\n',toc);

%% Step 2: generate event time magnitudes
tic;
sim_space(2,:) = abs(normrnd(mag_mean, mag_sd, [1, num_tq])) .* sim_space(1,:);

if (gen_all_plots == 1) || (gen_sig_gen_plots == 1)
    figure();
    plot(sim_space(2,:),'-o')
    hold on;
    title('Second Sim Space: impulses at event locations w gaussian magnitude');
end
fprintf('Step 2 duration: %d\n',toc);

%% Step 2a: generate event / time pairs
tic;
input_hits = [];
for index = 1:num_tq
    if (sim_space(2,index) ~= 0)
        input_hits = [input_hits, [index ; sim_space(2,index)]];
    end
end
fprintf('Step 2a duration: %d\n',toc);


%% Step 3: generate noise events
tic;
sim_space(3,:) = sim_space(2,:)+normrnd(0, noise_det_rms, [1, num_tq]);
if (gen_all_plots == 1) || (gen_sig_gen_plots == 1)
    figure();
    plot(sim_space(3,:),'-o')
    hold on;
    title('Third Sim Space: noise before shaper');
end
fprintf('Step 3 duration: %d\n',toc);
%%%TODO: investigate if this is the distro of noise we'd actually see.  May
%%%not make a difference in principle, but worth looking into. (pink noise)
%% Step 4: generate shaping
%first pass at shaping alg.
tic;

%reminder of gaussian : g(x) = (1/(sd*sqrt(2*pi)))*exp((-1/2)*((x-mean)/sd)^2)
%we say that the shaping time is equal to 6 sigma.  So:
shape_sd = t_shape / tq / 6;
%we put the mean 3 sigma in to give a 3 sigma leadup and a 6 sigma tail
shape_mean = shape_sd*3;
%now we generate the shape coeff
shape_coeff = (1/(shape_sd*sqrt(2*pi)))*exp((-1/2)*(((1:9*shape_sd)-shape_mean)/shape_sd).^2);
%now we normalize the shape_coeff
%shape_coeff = shape_coeff / (1/(shape_sd*sqrt(2*pi)));
%above was old.  Now we integrate.
shape_coeff = shape_coeff/sum(shape_coeff);
hit_scale_factor = max(shape_coeff);

if (gen_all_plots == 1) || (gen_sig_gen_plots == 1)
    figure();
    plot(shape_coeff);
    hold on;
    title('Signal Shape');
end

for index = (1:length(shape_coeff))
    if index == 1
        sim_space(4,:) = sim_space(3,:)*shape_coeff(1);
    else
        sim_space(4,:) = sim_space(4,:) + horzcat(zeros(1,index-1),sim_space(3,(1:num_tq-(index-1))))*shape_coeff(index);
    end
end

if (gen_all_plots == 1) || (gen_sig_shaping_plots == 1)
    figure();
    plot(sim_space(4,:),'o-')
    hold on;
    title('Fourth Sim Space: shaped events')
end
fprintf('Step 4 duration: %d\n',toc);

%% Step 5: add electronic noise to output
tic;
sim_space(5,:) = sim_space(4,:)+normrnd(0, noise_shaper_rms, [1, num_tq]);
if (gen_all_plots == 1) || (gen_sig_shaping_plots == 1)
    figure();
    plot(sim_space(5,:),'-o')
    hold on;
    title('Fifth Sim Space: noise after shaper');
end
fprintf('Step 5 duration: %d\n',toc);

%%%TODO: investigate if this is the distro of noise we'd actually see.  May
%%%not make a difference in principle, but worth looking into. (pink noise)
%% Step 6: generate sampling times
tic;
dig_space(2,:) = sim_space(5,(1:tq_sample:num_tq));
dig_space(1,:) = (1:tq_sample:num_tq);
if (gen_all_plots == 1) || (gen_dig_plots == 1)
    figure();
    plot(dig_space(2,:),'-o')
    hold on;
    title('Dig Space 2: Sampled Signal');
    figure();
    plot(sim_space(5,:),'-o')
    hold on;
    plot(dig_space(1,:),dig_space(2,:), 'r*', 'MarkerSize', 10)
    title('Samples overlayed on shaped signal');
end
fprintf('Step 6 duration: %d\n',toc);

%% Step 7: quantize digital samples
%max value calculation
tic;
max_val = 2^sample_bits-1;
max_analog_val = mag_sd*5*hit_scale_factor;
offset_val = mag_sd*0.8*hit_scale_factor;

%TODO insert clipping function here.
dig_space(3,:) = dig_space(2,:); %placeholder for clipped values

%perform quantization, using floor
dig_space(4,:) = floor((dig_space(3,:)+offset_val)/max_analog_val*max_val);

%quantization error:
dig_space(5,:) = (dig_space(4,:)-(dig_space(2,:)+offset_val)/max_analog_val*max_val);

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

fprintf('Quantization Error - Mean: %d counts, Sigma: %d counts\n',mean(dig_space(5,:)), std(dig_space(5,:)));

fprintf('Step 7 duration: %d\n',toc);

%% Step 8: non-linear cliping operation
tic;
dig_space(6,:) = dig_space(4,:);
for index = (1:num_samples)
    if dig_space(6,index) > max_val
        dig_space(6,index) = max_val;
    else if dig_space(6,index) < 0
            dig_space(6,index) = 0;
        end
    end
end

if (gen_all_plots == 1) || (gen_dig_plots == 1)
    figure();
    plot(dig_space(6,:), '-og');
    hold on;
    title('Dig space 6: clipped samples, approximation of what we would see from an ADC');
    
    figure();
    plot(dig_space(4,:), '-ob');
    hold on;
    plot(dig_space(6,:), '-*g', 'MarkerSize', 10);
    title('Dig space 4 and 6 (before and after non-linear clipping)');
end
fprintf('Step 8 duration: %d\n',toc);

%% Step 9: Pedestal Subtraction
tic;
dig_space(7,:) = dig_space(6,:) - offset_val/max_analog_val*max_val;

if (gen_all_plots == 1) || (gen_hit_find_plots == 1)
    figure();
    plot(dig_space(7,:), '-og');
    hold on;
    title('Dig space 7: Pedestal Subtraction');
end
fprintf('Step 9 duration: %d\n',toc);

%% Step 10: hit-finding
tic;
%
%Here we introduce the concept of the hit locker
%hit_locker{1,:} contains a vector.  The first element in the vector
%contains the sample number from dig_space(7,:).  The remaining elements in
%the vector contain values from dig_space(7,:) starting at the sample
%number.
%hit_locker{2,:} contains the fit data from each of the hits located above
%hit_locker{3,:} contains the rsquared value from the fit
%hit_locker{4,:} contains the equivilant analog value of the hit
%hit_locker{5,:} contains the equivilant time of the hit
%hit_locker{6,:} contains the index of the input_hit vector most closely
%hit_locker{7,:} contains the time difference between the hit and it's
%match
%hit_locker{8,:} contains the amplitude difference between the hit and it's
%match
%associated with the hit
%Set threshold based on the DC offset and the noise floor:
threshold = ((noise_det_rms+noise_shaper_rms)*3)/max_analog_val*max_val+10;
%threshold is the offset value, plus 3 sigma of the noise plus 10 counts.
%TODO eval if this threshold is actualy what it claims to be

%number of samples on the edge of a hit to grab:
edge_grab = 2;

if (gen_all_plots == 1) || (gen_hit_find_plots == 1)
    figure();
    plot(dig_space(7,:), '-o');
    hold on;
    line([1,num_samples],[threshold,threshold]);
    title('Samples with threshold drawn');
end

%Now we find hits
index = 1;  %index holds where in the sample space we are
hit_id = 1; %hit_id holds the ID of the hit to be written out
%hit_index = 1; %hit_index holds the index inside the hit of the sample
%end_index = 1; %makes coding a bit easier at the end, can probably be
%simplified away...
%the hit locker is not pre-allocated.  this is inefficient.  Cannot be
%pre-allocated because we do not pre-determine the size of the hit locker
%each hit locker array has a format like this:
%index 1 = index of first sample in the locker (hit location)
%all subsequent values are the values in the locker

while index <= num_samples
    %entering the loop we assume that we are not in the middle of a hit
    if dig_space(7,index) >= threshold;
        %if this is true, then we have just come to the first above
        %threshold sample of a hit.
        if index <= edge_grab %make sure we don't try to access samples which don't exist
            hit_locker{hit_id}(1) = 1;%assign 1 to the hit location
            hit_index = 2;
            while (hit_index-1 <= index-1)
                hit_locker{hit_id}(hit_index) = dig_space(7,hit_index-1);
                hit_index = hit_index + 1;
            end
        else
            hit_locker{hit_id}(1) = index-edge_grab;
            hit_index = 2;
            while (hit_index-1 <= edge_grab)
                hit_locker{hit_id}(hit_index) = dig_space(7,index-edge_grab+hit_index-2);    %hit_index-2 is zero for the first iteration
                hit_index = hit_index + 1;
            end
        end
        %ok, now we have gathered all of the previous hits, time to start
        %gathering the normal hits
        while (dig_space(7,index) >= threshold)
            hit_locker{hit_id}(hit_index) = dig_space(7,index);
            hit_index = hit_index + 1;
            index = index + 1;
        end
        
        %Now we have gathered all of the pulses above threshold, time for
        %the end grab at the other end
        end_index = 1;
        
        if (index + edge_grab) > num_samples
            while index + end_index <= num_samples
                hit_locker{hit_id}(hit_index) = dig_space(7,index+end_index-1);
                end_index = end_index + 1;
                hit_index = hit_index + 1;
            end
        else
            while end_index <= edge_grab
                hit_locker{hit_id}(hit_index) = dig_space(7,index+end_index-1);
                end_index = end_index + 1;
                hit_index = hit_index + 1;
            end
        end
        %note that the end grab doesn't touch the index.  We should be
        %ending on an index where we are below threshold, so theoretically
        %if there is another hit that starts during the edge grab this will
        %be picked up as a different event.
        hit_id = hit_id + 1;
    else
        index = index + 1;
    end
end

num_hits = length(hit_locker);
%now plot all of the hits...
%first generate some xy pairs
x_hit=[];
y_hit=[];
for index = (1:num_hits)
    for sub_index = (2:length(hit_locker{index}))
        y_hit = horzcat(y_hit,hit_locker{index}(sub_index));
        x_hit = horzcat(x_hit,hit_locker{index}(1)+sub_index-2);
    end
end

if (gen_all_plots == 1) || (gen_hit_find_plots == 1)
    figure();
    plot(x_hit, y_hit, '-xr');
    hold on;
    line([1,num_samples],[threshold,threshold]);
    title('Hits only with threshold line');
    figure();
    plot(dig_space(7,:), '-o');
    hold on;
    plot(x_hit, y_hit, 'xr','MarkerSize',15);
    line([1,num_samples],[threshold,threshold]);
    title('Hits only with threshold line');
end
fprintf('Step 10 duration: %d\n',toc);

%% Step 11: Hit re-construction
tic;
for index = (1:num_hits)
    [temp_curve, temp_gof] = fit((1:length(hit_locker{1,index})-1)',hit_locker{1,index}(2:length(hit_locker{1,index}))','gauss1', 'Upper', [max_val, inf, inf], 'Lower', [-1*max_val,0,0]);
    hit_locker{2,index} = temp_curve;
    hit_locker{3,index} = temp_gof.rsquare;
    %figure();
    %plot(hit_locker{2,index},(1:length(hit_locker{1,index})-1),hit_locker{1,index}(2:length(hit_locker{1,index})));
    %waitforbuttonpress;
end
r_squareds = [hit_locker{3,:}];
fprintf('Average R^2 = %d Average sigma = %d \n',mean(r_squareds),std(r_squareds));
fprintf('Step 11 duration: %d\n',toc);

%% Step 12: Hit re-construction 2
tic;
% hit_locker{4,:} = zeros(num_hits,1);%implement pre-allocation in the
% future.
% hit_locker{5,:} = zeros(num_hits,1);
for index = (1:num_hits)
    hit_locker{4,index} = hit_locker{2,index}.a1*max_analog_val/max_val; %y values
    hit_locker{5,index} = (hit_locker{2,index}.b1+hit_locker{1,index}(1))*tq_sample-shape_mean; %x values
end
if (gen_all_plots == 1) || (gen_hit_find_plots == 1)
    
    figure();
    plot([hit_locker{5,:}],[hit_locker{4,:}], '*r','MarkerSize', 10);
    hold on;
    plot(input_hits(1,:),input_hits(2,:)*hit_scale_factor, 'ob');
    title('Hit Locations overlayed on Orig. Hits');
    
    figure();
    plot([hit_locker{5,:}],[hit_locker{4,:}], '*r','MarkerSize', 10);
    hold on;
    plot(sim_space(4,:), '-o');
    title('Hit Locations overlayed on Orig. Shaped Hits');
end
fprintf('Step 12 duration: %d\n',toc);

%% Step 13: Hit reconstruction error estimation
tic;
%In this first pass we find the orig hit which is closest to the
%reconstructed hit
%we will locate in time only for the time being.

for index = 1:num_hits
    [~,hit_locker{6,index}] = min(abs(input_hits(1,:)-hit_locker{5,index})); %here we locate the smallest time difference
    hit_locker{7,index} = hit_locker{5,index}-input_hits(1,hit_locker{6,index}); %publish the time difference
    hit_locker{8,index} = hit_locker{4,index}-input_hits(2,hit_locker{6,index}); %publish amplitude difference
end

%now detect if any cells have the same index as their previous
error_locker=[];
hit_duplicates=[];
[~,hit_duplicates] = find([hit_locker{6,:},0]-[0,hit_locker{6,:}]==0);
%the above line finds duplicates, but only duplicates which are next to
%eachother.
num_errors = length(hit_duplicates);

while (~isempty(hit_duplicates))    %now go over each of the duplicates and remove the error hit.
    for index = 1:length(hit_duplicates)
        %we access the hit locker in the following way
        %hit_locker{7,hit_duplicates(index)-(index-1)}
        %The seven specifies the correct line in the hit locker
        %the hit_duplicates vector contains the list of index in the hit
        %locker of all hits which are duplicated.  As we cycle through the
        %for loop, we remove a hit from the hit locker with every
        %iteration, so we add the -(index-1) to account for this
        if (abs(hit_locker{7,hit_duplicates(index)-(index-1)}) < abs(hit_locker{7,hit_duplicates(index)-1-(index-1)}))
            error_locker = [error_locker,hit_locker(:,hit_duplicates(index)-1-(index-1))]; %move the error to the error locker
            hit_locker(:,hit_duplicates(index)-1-(index-1)) = []; %remove the data from the hit locker
        else
            error_locker = [error_locker,hit_locker(:,hit_duplicates(index)-(index-1))]; %move the error to the error locker
            hit_locker(:,hit_duplicates(index)-(index-1)) = []; %remove the data from the hit locker
        end
        num_hits = num_hits - 1;
    end
[~,hit_duplicates] = find([hit_locker{6,:},0]-[0,hit_locker{6,:}]==0);
end
fprintf('Step 13 duration: %d\n',toc);


%% Step 14: analysis
tic;
[~,num_good_hits]=size(hit_locker);
fprintf('Number of hits input: %d \n',num_input_pulses);
fprintf('Number of hits recovered: %d \n',num_good_hits);
fprintf('Number of errors detected: %d \n',num_errors);

timing_error_mean = mean([hit_locker{7,:}])*tq;
timing_error_std = std([hit_locker{7,:}])*tq;
amplitude_error_mean = mean([hit_locker{8,:}]);
amplitude_error_std = std([hit_locker{8,:}]);

fprintf('Timing error: %d, Sigma: %d\n',timing_error_mean,timing_error_std);
fprintf('Amplitude error: %d, Sigma: %d\n',amplitude_error_mean,amplitude_error_std);


if (gen_all_plots == 1) || (gen_hit_find_plots == 1)
    
    figure();
    hist([hit_locker{7,:}]*tq,30);
    hold on;
    title('Error in timing reconstruction');
    
    figure();
    hist([hit_locker{8,:}],30);
    hold on;
    title('Error in amplitude reconstruction');
end
fprintf('Step 14 duration: %d\n',toc);
