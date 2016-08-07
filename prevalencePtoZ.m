function prevalencePtoZ(prefix)

% convert p-value images to z-value images in results of prevalence
%
% prevalencePtoZ(prefix)
%
% See also prevalence
%
%
% Copyright (C) 2016 Carsten Allefeld
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.

if (nargin == 0) || isempty(prefix)
    prefix = 'prevalence_';
end

fprintf('prevalencePtoZ\n')
pattern = [prefix 'p*.nii'];
d = dir(pattern);
if numel(d) == 0
    fprintf('no files "%s" found\n', pattern)
end

[pathstr, name, ~] = fileparts(prefix);
for ni = 1 : numel(d)
    fprintf('  %s\n', d(ni).name)
    ifn = fullfile(pathstr, d(ni).name);
    ofn = [prefix 'z' d(ni).name(numel(name) + 2 : end)];
    Vi = spm_vol(ifn);
    p = spm_read_vols(Vi);
    z = sqrt(2) * erfinv(1 - 2 * p);
    saveMRImage(z, ofn, Vi.mat, ['z' d(ni).name(numel(name) + 2 : end - 4)])
end
