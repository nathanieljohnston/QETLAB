%%  IsLocallyPSD    Determines whether or not a matrix is k-locally positive semidefinite
%   This function has two required arguments:
%     X: The matrix to check k-local PSDness of
%     K: The parameter K that determines which level of K-locally PSD
%        is being tested
%
%   Returns either 1 or 0, indicating "yes" or "no".
%
%   See paper "Absolutely k-Incoherent Quantum States and Spectral
%   Inequalities for Factor Width of a Matrix" by Johnston et al for
%   details.

%   requires: IsPSD.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: August 29, 2024

function [ilp,wit] = IsLocallyPSD(X,k)

sX = length(X);
wit = 0;

sub_ind = nchoosek(1:sX(1),k);
sub_ind_len = size(sub_ind,1);
for j = 1:sub_ind_len
    if(~IsPSD(X(sub_ind(j,:),sub_ind(j,:))))
        ilp = 0;
        wit = [sub_ind(j,:);sub_ind(j,:)];
        return;
    end
end

ilp = 1;