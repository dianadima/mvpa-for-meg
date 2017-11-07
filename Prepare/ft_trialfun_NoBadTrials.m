function [trl, event] = ft_trialfun_NoBadTrials(cfg);

hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

% search for "trigger" events
sample  = [event(strcmp(cfg.trialdef.eventtype, {event.type})).sample]';
pretrig  = -round(cfg.trialdef.prestim  * hdr.Fs);
posttrig =  round(cfg.trialdef.poststim * hdr.Fs);

trl = [];
for j = 1 : length(sample)
   trlbegin = sample(j) + pretrig;       
    trlend   = sample(j) + posttrig;       
    offset   = pretrig;
    newtrl   = [trlbegin trlend offset];
    trl      = [trl; newtrl];
end
fprintf('Found %d trials\n', length(sample));
cls = readClassFile([cfg.dataset '/ClassFile.cls']);

%This sectino of code rereads the markerfile to work out which trials need to be deleted 
%It covers the case when only a subset of trials are read in
Markers = readmarkerfile(cfg.dataset);
for i = 1 : length(Markers.marker_names)
   if strcmp(cfg.trialdef.eventtype, Markers.marker_names{i}) == 1
TrialIndexes = Markers.trial_times{i}(:,1);  
   end
end
Remove = [];
for i = 1 : length(TrialIndexes)
   if ismember(TrialIndexes(i), cls.trial)
      Remove = [Remove i];
   end
end
fprintf('Removing %d bad trials\n', length(Remove));

%%%
%%%%%%%%%%%%%%
trl(Remove, :) = [];


