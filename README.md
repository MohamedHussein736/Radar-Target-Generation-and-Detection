# Radar-Target-Generation-and-Detection
https://confirm.udacity.com/9R4N6A4H

![SensorFusion_Certificate](sensorFusion_certificate.png)

#### 1. FMCW Waveform Design
Using the given system requirements, design a FMCW waveform. Find its Bandwidth (B), chirp time (Tchirp) and slope of the chirp.

```Matlab
%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
range_max = 200
% Range Resolution = 1 m
delta_r = 1
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 3*10^8
pos = 80
vel = 40
B = c/(2*delta_r)
Tchirp = 5.5*(range_max*2/c)
slope = B/Tchirp
disp(slope)
```

#### 2. Simulation Loop
Simulate Target movement and calculate the beat or mixed signal for every timestamp.

```
%% Signal generation and Moving Target simulation
% Running the radar scenario over the time.

for i=1:length(t)
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity.
    r_t(i) = Range_of_target + (Velocity_of_target*t(i));
    td(i) = (2*r_t(i))/speed_of_light;
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal.
    Tx(i) = cos(2*pi*(fc*t(i) + 0.5*Slope*t(i)^2));
    Rx(i)  = cos(2*pi*(fc*(t(i)-td(i)) + 0.5*Slope*(t(i)-td(i))^2));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i);
end
```

#### 3. Range FFT (1st FFT)

Implement the Range FFT on the Beat or Mixed Signal and plot the result.

```Matlab
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
sig_fft1 = fft(Mix, Nr);

% sig_fft1 = sig_fft1./Nr;
% *%TODO* :
% Take the absolute value of FFT output
sig_fft1 = abs(sig_fft1);
sig_fft1 = sig_fft1./max(sig_fft1); % Normalize

% *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
sig_fft1 = sig_fft1(1:Nr/2-1);
```

#### 4. 2D CFAR
Implement the 2D CFAR process on the output of 2D FFT operation, i.e the Range Doppler Map.

Determine the number of Training cells for each dimension. Similarly, pick the number of guard cells.

```Matlab
Tr = 10;
Td = 8;

%Select the number of Guard Cells in both dimensions around the Cell under
%test (CUT) for accurate estimation
Gr = 4;
Gd = 4;
% offset the threshold by SNR value in dB
offset = 1.4;
```

Slide the cell under test across the complete matrix. Make sure the CUT has margin for Training and Guard cells from the edges.

```Matlab
RDM = RDM/max(max(RDM)); % Normalizing

for i = Tr+Gr+1:(Nr/2)-(Tr+Gr)
    for j = Td+Gd+1:(Nd)-(Td+Gd)
        %Create a vector to store noise_level for each iteration on training cells
        noise_level = zeros(1,1);
        
        for p = i-(Tr+Gr) : i+(Tr+Gr)
            for q = j-(Td+Gd) : j+(Td+Gd)
                if (abs(i-p) > Gr || abs(j-q) > Gd)
                    noise_level = noise_level + db2pow(RDM(p,q));
                end
            end
        end
```

For every iteration sum the signal level within all the training cells. To sum convert the value from logarithmic to linear using db2pow function.

Average thesummed values for all of the training cells used. After averaging convert it back to logarithmic using pow2db.
Further add the offset to it to determine the threshold.

Next, compare the signal under CUT against this threshold.
If the CUT level > threshold assign it a value of 1, else equate it to 0.


To keep the map size same as it was before CFAR, equate all the non-thresholded cells to 0.

```Matlab
% Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
% CFAR

RDM = RDM/max(max(RDM)); % Normalizing

for i = Tr+Gr+1:(Nr/2)-(Tr+Gr)
    for j = Td+Gd+1:(Nd)-(Td+Gd)
        %Create a vector to store noise_level for each iteration on training cells
        noise_level = zeros(1,1);
        
        for p = i-(Tr+Gr) : i+(Tr+Gr)
            for q = j-(Td+Gd) : j+(Td+Gd)
                if (abs(i-p) > Gr || abs(j-q) > Gd)
                    noise_level = noise_level + db2pow(RDM(p,q));
                end
            end
        end
        
        % Calculate threshould from noise average then add the offset
        threshold = pow2db(noise_level/(2*(Td+Gd+1)*2*(Tr+Gr+1)-(Gr*Gd)-1));
        threshold = threshold + offset;
        CUT = RDM(i,j);
        
        if (CUT < threshold)
            RDM(i,j) = 0;
        else
            RDM(i,j) = 1;
        end
        
    end
end

```
Result
![result](https://github.com/MohamedHussein736/Radar-Target-Generation-and-Detection/blob/master/img/3.jpg)
![result](https://github.com/MohamedHussein736/Radar-Target-Generation-and-Detection/blob/master/img/4.jpg)
![result](https://github.com/MohamedHussein736/Radar-Target-Generation-and-Detection/blob/master/img/1.jpg)
![result](https://github.com/MohamedHussein736/Radar-Target-Generation-and-Detection/blob/master/img/2.jpg)
![result](https://github.com/MohamedHussein736/Radar-Target-Generation-and-Detection/blob/master/img/5.jpg)
![result](https://github.com/MohamedHussein736/Radar-Target-Generation-and-Detection/blob/master/img/6.jpg)

