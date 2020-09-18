function [Ir, Ig, Ib, Irgb] = readJpg(filepathIn, mouseNum, sectionNum)
% readJpg
filepathTotal = [filepathIn, mouseNum, '\', sectionNum, '.jpg'];
Irgb = imread(filepathTotal);
Ir = Irgb(:,:,1);
Ig = Irgb(:,:,2);
Ib = Irgb(:,:,3);

% uncomment below to combine RGB image to test:
% Irgb = zeros(size(Ib,1),size(Ib,2),3);
% Irgb(:,:,1) = double(Ir)./ double(max(max(Ir))); 
% Irgb(:,:,2) = double(Ig)./ double(max(max(Ig)));
% Irgb(:,:,3) = double(Ib)./ double(max(max(Ib)));
% imshow(Irgb);
% keyboard
end %end function
