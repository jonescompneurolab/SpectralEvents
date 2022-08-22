% Run this after example.m in order to write the detected event times to disk

ts_data = [];

num_subj = numel(timeseries);
num_trials = 200;

event_times = cell(num_subj * num_trials, 1);

for subj_idx = 1:num_subj
    ts_data = [ts_data, timeseries{subj_idx}];
    trial_idxs = specEvents(subj_idx).Events.Events.trialind + ((subj_idx - 1) * num_trials);
    %trial_idxs = specEvents(subj_idx).Events.Events.trialind;



    for trial_idx=1:length(trial_idxs)
        event_time = specEvents(subj_idx).Events.Events.maximatiming(trial_idx);
        event_times{trial_idxs(trial_idx)} = [event_times{trial_idxs(trial_idx)}, event_time];
    end
end

save("beta_events_shin_2017.mat", "ts_data", "event_times", "Fs")