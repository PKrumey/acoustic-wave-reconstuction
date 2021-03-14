function RCnorm = norm2unp(RC)
% This function normalizes the RC matrix to the maximum of the unpumped RC,
% which is supposed to be in the first column.
% mini = min(RC(:));
% RCnorm = RC - mini;

maxi = max(RC);
maxOFF = maxi(1);

RCnorm = RC / maxOFF;
end

