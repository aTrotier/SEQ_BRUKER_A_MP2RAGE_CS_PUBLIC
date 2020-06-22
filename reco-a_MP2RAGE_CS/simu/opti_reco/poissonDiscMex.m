%function [mask] = poissonDiscMex(fovy,fovz,sky,skz,ry,rz,ncal,cutcorners, pp)
%
% function computes a variable density poisson-disc sampling mask
%
% Inputs:
%	 fovx, fovy - FOV of the actual image in cm
% 	 sky, skz, - size of the kspace
%	 ry, rz - acceleration in x and y dimensions. 
%	 ncal -  number of autocalibration lines
%        cutcorners - flag to cut corners. 
%	 pp      - polynomial order of variable density (0 = uniform)
%
%
% Outputs:
%	mask - a binary mask. 
%
%
%   Original function written by Marcus Alley based on a Modification of a python script from 
%   http://devmag.org.za/2009/05/03/poisson-disk-sampling/ . The Mex interface was written by 
%   Michael Lustig 
%
% (c) Michael Lustig 2013 & Aurelien J Trotier 2020


