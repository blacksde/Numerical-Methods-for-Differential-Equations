% load the image from celine1.jpg
fimage = imread('iron_man_3.jpg');
PQ = size(fimage);

% fast fourier transform to get the initial condition of the equation
Finitial = zeros(PQ);
for i = 1:3
    Finitial(:,:,i) = fft2(fimage(:,:,i));
end




% construct matrix with entries kx^2+ky^2
Kx = 2*pi*[0:PQ(2)/2,-PQ(2)/2+1:-1]/PQ(2);
Ky = 2*pi*transpose([0:PQ(1)/2,-PQ(1)/2+1:-1])/PQ(1);
% M2x = ones(PQ(1),1)*Kx.^2;    % matrix with entries kx^2
% M2y = Ky.^2*ones(1,PQ(2));    % matrix with entries ky^2
[M2x,M2y] = meshgrid(Kx.^2,Ky.^2);
M2 = M2x+M2y;

% create a new avi file
vidObj = VideoWriter('ironman.avi');
vidObj.FrameRate = 100;
open(vidObj); 

% solve the heat equation after fourier transform
% equation is dF(u)/dt = -(kx^2+ky^2)F(u)
Fvalue = zeros(PQ);
imagevalue = zeros(PQ);

N = 300;

t = linspace(0,1000,N);

for i = 1:N
    for j = 1:3
        Fvalue(:,:,j) = exp(-M2*t(i)).*Finitial(:,:,j);
        % use the explict formula to get the value at time t(i)
        imagevalue(:,:,j) = round(ifft2(Fvalue(:,:,j)));
    end
    
    % write the image to video 
    currFrame = uint8(imagevalue);
    writeVideo(vidObj,currFrame);
    disp(num2str(i));
    
end

% get the last frame
imwrite(uint8(imagevalue), 'final.jpg','jpg');
close(vidObj); 

