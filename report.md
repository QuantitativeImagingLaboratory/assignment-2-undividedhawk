1)
Forward_transform(): Since the input is an n*n matrix I make a complex matrix initialized to 0.
using 4 nested loops i iterated through n to calculate e^-angle according to the 2 dimensional formula
for the dft. then store the value of tempJ + matrix[i,j]*e^-angle into a complex variable called tempJ.
exiting the innermost loop store tempI + tempJ into a new complex variable called tempI. Then reset tempJ to 0.
This is the summation from the formula. exiting your second innermost loop, store temp I into finalmatrix[u,v].
and set temp I to 0. exiting all loops and finalMatrix is the DFT of the n*n matrix.


Inverse_transform(): Same exact process as before but use +angle instead of -angle and divide tempI.real by n^2 then round the result
store the result into your finalMatrix and finalMatrix is the IDFT of the input matrix.


Discrete_cosine_transform(): I set up 4 nested loops set to iterate through n. when your on the edge of the matrix set constants a, b, to there respective values
from the formula of the discrete cosine transform. then set angle1 = cos(pi/n*(i+.5)*u)
and angle2 = cos(pi/n*(j+.5)*v). and set e = angle1*angle2. use a temp variable tempJ = tempJ + matrix[i,j]*e*a*b. exit innermost loop and set a new temp variable
tempI = tempI + tempJ. set tempJ = 0. and exit the next innermost loop. set finalMatrix[u,v] = tempI and reset tempI = 0. exit all loops and finalMatrix is
the cosine transform of the input matrix.


Magnitude(): computes the DFT of the input matrix, and initializes a finalMatrix to all 0 entries. have a nested loop iterate through n
set a = real number of the dftMatrix[u,v]^2 and b = imaginary number of the dftMatrix[u,v]^2.
and set the finalMatrix[u,v] = squareroot(a+b). exit both loops and the result of finalMatrix is the magnitude
of the dftInputMatrix.


2)
get_ideal_low_pass_filter(): initialize a zero matrix of equal size of the image stored as mask. create a nested loop iterate through the range of the image
set variable d equal the distance from the current point in iteration to the center coordinate in the matrix. if the distance is less than or equal to the cutoff 
point specified set the mask[u,v] = 1. else set mask[u,v] = 0. exit both loops and return the mask which is the ideal low pass filter.


get_ideal_high_pass_filter(): initialize a matrix b = low pass filter of the image. and initialize mask to be a zero matrix of the same size. take the negative of 
b stored as mask, return mask which is the ideal_high_pass_filter 


get_butterworth_low_pass_filter(): initialize zero matrix called mask, nested for loop iterating through image dimensions. store distance between current interation values
x and y and the center of the image as d. set mask[x,y] = 1/(1+(d/cutoff_value)^(2*order)) exit both loops. return mask which is the butterworth high pass filter.


get_butterworth_high_pass_filter(): uses same steps as butterworth_low_pass_filter but switches the values of cutoff_value and d from the formula used to find the mask.


get_gaussian_low_pass_filter(): initialize zero matrix called mask, nested for loops iterate through dimensions of image. store the distance from the current iteration values to the 
center coordinate as the value d. set mask[x,y] = e^(-d**2/(2*(cutoff_value^2))). exit both loops and return value of mask, which is the gaussian low pass mask.


get_gaussian_high_pass_filter():  initialize b = the low pass gaussian filter of the input image. initialize a zero matrix of same dimensions as mask. nested for loops iterate through 
dimensions of the input image set mask[x,y] = 1-b[x,y]. exit both loops and return mask which is the gaussian high pass filter.

filtering(): set shape to be the coordinates of the image, cutoff to be the cutoff argument, and the order to be the order argument. store the DFT of the image as fftimage, shift the low frequencies 
of the DFT to be in the center of the image and store as fftimageshift. store the logarithmic compression of fftimageshift as magnitude_spectrum. set mask using the specified get_filter function. initialize
a zero matrix of the same size as the input image stored as maskedimage. nested for loops iterate through coordinates of the input image set maskedimage[x,y] = mask[x,y]*fftimageshift[x,y]. exit both loops
reverse the shift done previously using the ifftshift function and store that matrix as fftmaskshift. Take the inverse fourier transform of fftmaskshift and store that matrix as fftinverse. set the absolute value of
fftinverse as img_back. store the logarithmic compression of maskedimage and store it as disp_maskedimage. return an array of img_back, magnitude_spectrum, disp_maskedimage. img_back is the filtered image, magnitude_spectrum is the 
viewable dft of the input image, and disp_maskedimage is the magnitude of the DFT after filtering.





