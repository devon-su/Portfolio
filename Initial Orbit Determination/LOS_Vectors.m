function LOS = LOS_Vectors(meas)
% LOS_VECTORS Calculation of Line-Of-Site Vectors
% 
%    [LOS] = LOS_Vectors(meas) Returns the line-of-site unit vectors 
%    given topocentric right ascension and declination observations.
%
% INPUTS:
% meas - nx2 matrix containing n observations with the first column holding
%        right ascention values and the second column holding declination
%        values [RA Dec] (degrees)
%
% OUTPUTS:
% LOS - 3xn matrix of LOS vectors with each column, n, corresponding to
%       each observation given in the input

LOS = [];
for i = 1:length(meas)
    RA = meas(i, 1);
    dec = meas(i, 2);
    LOS = [LOS; [cosd(dec)*cosd(RA) cosd(dec)*sind(RA) sind(dec)]];
end

LOS = LOS.';

end

