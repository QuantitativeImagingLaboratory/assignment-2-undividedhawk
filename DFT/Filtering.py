# For this part of the assignment, You can use inbuilt functions to compute the fourier transform
# You are welcome to use fft that are available in numpy and opencv
import cv2
import math as c
import numpy as np

class Filtering:
    image = None
    filter = None
    cutoff = None
    order = None

    def __init__(self, image, filter_name, cutoff, order = 0):
        """initializes the variables frequency filtering on an input image
        takes as input:
        image: the input image
        filter_name: the name of the mask to use
        cutoff: the cutoff frequency of the filter
        order: the order of the filter (only for butterworth
        returns"""
        self.image = image
        if filter_name == 'ideal_l':
            self.filter = self.get_ideal_low_pass_filter
        elif filter_name == 'ideal_h':
            self.filter = self.get_ideal_high_pass_filter
        elif filter_name == 'butterworth_l':
            self.filter = self.get_butterworth_low_pass_filter
        elif filter_name == 'butterworth_h':
            self.filter = self.get_butterworth_high_pass_filter
        elif filter_name == 'gaussian_l':
            self.filter = self.get_gaussian_low_pass_filter
        elif filter_name == 'gaussian_h':
            self.filter = self.get_gaussian_high_pass_filter

        self.cutoff = cutoff
        self.order = order

    def get_ideal_low_pass_filter(self, shape, cutoff):
        """Computes a Ideal low pass mask
        takes as input:
        shape: the shape of the mask to be generated
        cutoff: the cutoff frequency of the ideal filter
        returns a ideal low pass mask"""
        rows = shape
        columns = shape

        for x in range(rows):
            for y in range(columns):
                if np.sqrt((x ** 2 + y ** 2)) <= cutoff:
                    self.image[x, y] = 1
                else:
                    self.image[x, y] = 0

        return self.image

    def get_ideal_high_pass_filter(self, shape, cutoff):
        """Computes a Ideal high pass mask
        takes as input:
        shape: the shape of the mask to be generated
        cutoff: the cutoff frequency of the ideal filter
        returns a ideal high pass mask"""

        # Hint: May be one can use the low pass filter function to get a high pass mask
        b = self.get_ideal_high_pass_filter(shape, cutoff)
        for x in shape[0]:
            for y in shape[1]:
                if b[x, y] == 0:
                    self.image[x, y] == 1
                else:
                    self.image[x, y] == 0

        return self.image

    def get_butterworth_low_pass_filter(self, shape, cutoff, order):
        """Computes a butterworth low pass mask
        takes as input:
        shape: the shape of the mask to be generated
        cutoff: the cutoff frequency of the butterworth filter
        order: the order of the butterworth filter
        returns a butterworth low pass mask"""
        for x in shape[0]:
            for y in shape[1]:
                return 1 / (1 + (self[x, y] / cutoff) ** (2 * order))
        return 0

    def get_butterworth_high_pass_filter(self, shape, cutoff, order):
        """Computes a butterworth high pass mask
        takes as input:
        shape: the shape of the mask to be generated
        cutoff: the cutoff frequency of the butterworth filter
        order: the order of the butterworth filter
        returns a butterworth high pass mask"""

        # Hint: May be one can use the low pass filter function to get a high pass mask

        return self.get_butterworth_low_pass_filter(cutoff, shape, order)

        return 0

    def get_gaussian_low_pass_filter(self, shape, cutoff):
        """Computes a gaussian low pass mask
        takes as input:
        shape: the shape of the mask to be generated
        cutoff: the cutoff frequency of the gaussian filter (sigma)
        returns a gaussian low pass mask"""

        for x in shape[0]:
            for y in shape[1]:
                return c.exp(self[x, y] ** 2 / (2 * (cutoff ** 2)))

        return 0

    def get_gaussian_high_pass_filter(self, shape, cutoff):
        """Computes a gaussian high pass mask
        takes as input:
        shape: the shape of the mask to be generated
        cutoff: the cutoff frequency of the gaussian filter (sigma)
        returns a gaussian high pass mask"""

        # Hint: May be one can use the low pass filter function to get a high pass mask
        self.get_gaussian_low_pass_filter(cutoff, shape)

        return 0

    def post_process_image(self, image):
        """Post process the image to create a full contrast stretch of the image
        takes as input:
        image: the image obtained from the inverse fourier transform
        return an image with full contrast stretch
        -----------------------------------------------------
        1. Full contrast stretch (fsimage)
        2. take negative (255 - fsimage)
        """


        return image


    def filtering(self):
        """Performs frequency filtering on an input image
        returns a filtered image, magnitude of DFT, magnitude of filtered DFT        
        ----------------------------------------------------------
        You are allowed to used inbuilt functions to compute fft
        There are packages available in numpy as well as in opencv
        Steps:
        1. Compute the fft of the image
        2. shift the fft to center the low frequencies
        3. get the mask (write your code in functions provided above) the functions can be called by self.filter(shape, cutoff, order)
        4. filter the image frequency based on the mask (Convolution theorem)
        5. compute the inverse shift
        6. compute the inverse fourier transform
        7. compute the magnitude

        8. You will need to do a full contrast stretch on the magnitude and depending on the algorithm you may also need to
        take negative of the image to be able to view it (use post_process_image to write this code)
        Note: You do not have to do zero padding as discussed in class, the inbuilt functions takes care of that
        filtered image, magnitude of DFT, magnitude of filtered DFT: Make sure all images being returned have grey scale full contrast stretch and dtype=uint8 
        """

        shape = self.image.shape[0]
        fftimage = np.fft.fft2(self.image)
        fftimageshift = np.fft.fftshift(fftimage)
        print("DFT:")
        print(fftimage)
        print("DFT Shifted:")
        print(fftimageshift)

        mask = self.filter(shape, 50)
        maskedimage = np.zeros((shape,shape), np.complex64)
        for x in range(shape):
            for y in range(shape):
                maskedimage[x,y] = mask[x,y] * fftimageshift[x,y]
        fftmaskshift = np.fft.ifftshift(maskedimage)
        fftinverse = np.fft.ifft2(fftmaskshift)
        print("mask")
        print(mask)
        print("masked image")
        print(maskedimage)
        print("fft inverse")
        print(fftinverse)
        return 0

        #return [fftinverse, dftmagnitude,idftmagnitude]
