function [ JD ] = JD(Y, M, D, h, m, s)
% JD    Julian date calculator
%
%   JD(Y, M, D, h, m, s) will compute the Julian date given a date in put
%   formatted by year Y, month M, day D, hour h, minute m and second s. All
%   inputs except for Y are optional. If inputs are omitted, they will be
%   set to the earliest possible time defined by the inputs given. For 
%   example, JD(2018) will return the Julian date at January 1, 0:00. 
%
%   This function does not account for leap seconds.
%
%   Class support for inputs Y, M, D, h, m, s
%     float: double, single

if nargin < 6
    s = 0;
    if nargin < 5
        m = 0;
        if nargin < 4
            h = 0;
            if nargin < 3
                D = 1;
                if nargin < 2
                    M = 1;
                    if nargin < 1
                        error(message('MATLAB:narginchk:notEnoughInputs'));
                    end
                end
            end
        end
    end
end
JD = 1721013.5+367*Y - fix(7/4*(Y+fix((M+9)/12)))+fix(275*M/9)+D+...
    (60*h+m+s/60)/1440;
end