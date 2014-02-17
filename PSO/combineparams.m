% combine continuous and discrete parameters to one vector
function params = combineparams(continuous, discrete)
continuous = squeeze(continuous); % squeeze N-D matrices
discrete = squeeze(discrete);
params = [discrete(1) continuous(1:7)' discrete(2) continuous(8:12)'];
