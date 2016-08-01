function V = saveMRImage(data, filename, mat, descrip)

% saves a 3d-matrix to an MR image file, using SPM routines
%
% V = saveMRImage(data, filename, mat = eye(4), descrip = 'saveMRImage')
%
% data:         image data to be written, 3d matrix
% filename:     name of the file to write to
% mat:          transformation from voxel indices to mm coordinates
% descrip:      string describing contents
% V:            spm volume struct
%
% The file format is determined by the filename extension.
% The data type is is determined by the class of data.
%
%
% This file is part of v1.0.0 of prevalence-permutation, see
% https://github.com/allefeld/prevalence-permutation/releases
%
% Copyright (C) 2013–2016 Carsten Allefeld
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.

if nargin < 3
    mat = eye(4);
end

if nargin < 4
    descrip = 'saveMRImage';
end

if size(data, 4) > 1
    warning('only the first volume is written')
end

if isfloat(data) && ~isreal(data)
    warning('only real part of complex data is written')
    data = real(data);
end

% determine data type
dt = class(data);
if isa(data, 'single')
    dt = 'float32';
elseif isa(data, 'double')
    dt = 'float64';
end
dt = spm_type(dt);
if isnan(dt)
    error('data type is not supported')
end

V = struct;
V.fname = filename;
V.dim = size(data);
V.dt = [dt 0];
V.mat = mat;
V.descrip = descrip;

V = spm_write_vol(V, data);

if nargout == 0
    clear V
end
