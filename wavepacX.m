% Time evolution of a Gaussian with Split-Operator propagation method
%
clear
close all
% physical parameters
hbar = 1;
m = 1;
omega=1; % potential parameter
% Time
Nt = 400;       % time steps
T  = 20;        % Propagation time
time = linspace(0,T,Nt);
dt = time(2)-time(1);
nout=5; % 
% Coordinate space
xmax = 6;       % Maximum value in coordinate space
N = 2^9;        % Number of points
x = linspace(-xmax,xmax,N+1)'; x(length(x))= [];
dx = x(2)-x(1);
% Momentum space
k = 2*pi*linspace(-0.5,0.5,N+1)'/dx;k(length(k))=[];
%kmax=max(k)
% Exponential of the kinetic energy
expT = exp(-i/hbar*(hbar*k).^2/(2*m)*dt);

% Exponential of the potential energy
V = 0.5*m*omega^2*x.^2; % potential
expV = exp(-i/hbar*V*dt/2);

% Gaussian wavepacket
sigma =1;  % Initial width
x0 = 0;         % initial position
p0 = .5;       % initial momentum
psi = exp(-(x-x0).^2/4/sigma^2+i/hbar*p0*x)*(2*pi*sigma^2)^(-1/4);

% Time evolution
ncount=1;
Psi(ncount,:)=psi;
for n = 2:length(time)
  psi = expV.*fftshift(ifft(fftshift(expT.*fftshift(fft(fftshift((expV.*psi)))))));
  if mod(n,nout)==0
       ncount=ncount+1;
       Psi(ncount,:)=psi;
    end
end
%
imagesc(abs(Psi))





