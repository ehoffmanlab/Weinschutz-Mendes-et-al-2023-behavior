% Function to start parallel pool with number of workers as currently
% accessible CPUs
%
% Jeff Eilbott, 2016, jeilbott@surveybott.com

function pool = start_parpool()
numCores = feature('numCores');
pool = gcp('nocreate');
if ~isempty(pool) && pool.NumWorkers < numCores
    delete(pool);
    pool = gcp('nocreate');
end
if isempty(pool)
    pool = parpool(numCores);
end
end