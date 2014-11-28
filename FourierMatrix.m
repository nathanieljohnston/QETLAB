%%  FOURIERMATRIX  Generates the unitary matrix that implements the quantum Fourier transform
%   This function has one required argument:
%     DIM: the size of the Fourier matrix
%
%   F = FourierMatrix(DIM) is the DIM-by-DIM unitary matrix that implements
%   the quantum Fourier transform.
%
%   URL: http://www.qetlab.com/FourierMatrix

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 30, 2012

function F = FourierMatrix(dim)

w = exp(2i*pi/dim); % primitive root of unity
F = (w.^((0:dim-1).'*(0:dim-1)))/sqrt(dim);