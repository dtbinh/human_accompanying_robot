function ell_set = ellipseBound(rec_vert)
% find ellipses for bounding rectangles
ell_set = zeros(4,size(rec_vert,1)); %[c;a;b]
for ii = 1:size(rec_vert,1)
    cor = rec_vert{ii};
    c = sum(cor,2)/2; % center of an ellipse
    w = abs(cor(1,2)-cor(1,1)); % width
    h = abs(cor(2,2)-cor(2,1)); % height
    a = w/sqrt(2); % half-length of long-axis
    b = h/sqrt(2); % half-length of short-axis
    ell_set(:,ii) = [c;a;b];
end    