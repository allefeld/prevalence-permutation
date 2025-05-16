function N = spmWriteImage(Y, fname, mat, varargin)

% writes a data array to a NIfTI image, using SPM routines
%
% N = spmWriteImage(Y, fname, mat = eye(4), ...)
%
% Y:              data to be written, up to 7D array; alternatively for Y =
%                 {dim, init} an image of size dim is initialized to init
% fname:          name of the image file to write to (extension img or nii).
% mat:            transformation from voxel indices to mm coordinates (sform)
% ...:            parameters in the form 'Name', Value; see below
% N:              SPM nifti object, may be ignored
%
% The file format is determined by the filename extension (img or nii).
%
% The data type is determined by the class of Y (or of init). Supported are
% double, single, int8, uint8, int16, uint16, int32, and uint32. Data of
% type logical are written as uint8.
%
% The data in the NIfTI file can be read back and further modified by
% accessing the file_array object N.dat. The original dimensions of the
% data are preserved if N.dat is accessed using the syntax
%   N.dat(:, :, :, :, :, :, :)
% except that the data are nominally always at least 2D.
%
% There is no spmReadImage; an existing NIfTI file can opened for reading
% and writing by
%   N = nifti(fname);
%
%
% It should not normally be necessary to set the following parameters:
%   mat_intent:   intention of mat (sform_code), default 'Aligned'
%   mat0:         alternative transformations (qform), default mat
%   mat0_intent:  intention of mat0 (qform_code), default mat_intent
%   descrip:      description of data, default ''
%   scl_slope     linear transformation slope, default 1
%   scl_inter     linear transformation intercept, default 0
% For background information, see
%   http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/
%
% It is however possible to store floating-point data that only attain
% discrete and linearly-spaced values (e.g. ratios) using an integer data
% type using the linear transformation feature of the NIfTI format,
%   spmWriteImage(int, fname, mat, 'scl_slope', slope, 'scl_inter', inter)
% where int is of class uint8, uint16, or uint32. Upon reading the file,
% the integer data are automatically transformed according to
%   Y = double(int) * slope + inter;
% This way the image file may be kept substantially smaller. Note however
% that 'scl_slope' and 'scl_inter' are stored only with single precision.
%
% See also spmCoords, nifti/Contents


% allow Y to be {dim, init}, initializing data
if ~iscell(Y)
    dim = size(Y);
else
    dim = Y{1};
    Y = Y{2};
end
% default argument values
if nargin < 3 || isempty(mat)
    mat = eye(4);
end
% name-value parameters & defaults
p = inputParser;
p.PartialMatching = false;
p.addParameter('mat_intent', 'Aligned');
p.addParameter('mat0', mat);
p.addParameter('mat0_intent', []);
p.addParameter('descrip', '');
p.addParameter('scl_slope', 1);
p.addParameter('scl_inter', 0);
p.parse(varargin{:})
p = p.Results;
if isempty(p.mat0_intent)
    p.mat0_intent = p.mat_intent;
end

% check dimensionality
assert(numel(dim) <= 7, 'data can only be up to 7D')

% set datatype according to class(Y)
if islogical(Y)
    Y = uint8(Y);
end
assert(isnumeric(Y), 'data must be numeric')
assert(isreal(Y), 'complex data are not supported')
dtype = class(Y);
switch dtype
    case 'single'
        dtype = 'FLOAT32';
    case 'double'
        dtype = 'FLOAT64';
    case {'int64', 'uint64'}
        error('64-bit integers are not supported')
    otherwise
        dtype = upper(dtype);
end
dtype = [dtype '-LE'];

% delete file of this name if it exists
spm_unlink(fname)

% create nifti object containing file_array object
N = nifti;
N.dat = file_array(fname, dim, dtype, 0, p.scl_slope, p.scl_inter);
N.mat = mat;
N.mat0 = p.mat0;
N.mat_intent  = p.mat_intent;
N.mat0_intent = p.mat0_intent;
N.descrip = p.descrip;

% write header
create(N)

% write data
if (p.scl_slope == 1) && (p.scl_inter == 0)
    N.dat(:) = Y(:);
else
    fprintf('b\n')
    N.dat(:) = double(Y(:)) * p.scl_slope + p.scl_inter;
    % It would be better to write to N.raw, but that seems not to be
    % possible. This way we transform from raw values to intended values,
    % only to have the nifti object transform them back to raw. :(
end
if nargout == 0
    clear N
end


% Copyright © 2016–2023 Carsten Allefeld
% SPDX-License-Identifier: GPL-3.0-or-later

