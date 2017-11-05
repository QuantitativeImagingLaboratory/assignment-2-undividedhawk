# For this part of the assignment, please implement your own code for all computations,
# Do not use inbuilt functions like fft from either numpy, opencv or other libraries
import cmath as c
import math
import numpy as np



class DFT:

    def forward_transform(self, matrix):
        """Computes the forward Fourier transform of the input matrix
        takes as input:
        matrix: a 2d matrix
        returns a complex matrix representing fourier transform"""
        n = len(matrix[0])
        finalMatrix = np.zeros((n, n), np.complex64)
        tempI = complex(0)
        tempJ = complex(0)

        for u in range(n):
            for v in range(n):
                for i in range(n):
                    for j in range(n):
                        angle = 2j * c.pi / n * ((u * i) + (v * j))
                        e = c.exp(-angle)
                        tempJ += matrix[i, j] * e
                    tempI += tempJ
                    tempJ = 0

                finalMatrix[u, v] = tempI
                tempI = 0
        return finalMatrix






    def inverse_transform(self, matrix):
        """Computes the inverse Fourier transform of the input matrix
        matrix: a 2d matrix (DFT) usually complex
        takes as input:
        returns a complex matrix representing the inverse fourier transform"""
        n = len(matrix[0])
        finalMatrix = np.zeros((n, n), np.complex64)
        tempI = complex(0)
        tempJ = complex(0)

        for u in range(n):
            for v in range(n):
                for i in range(n):
                    for j in range(n):
                        angle = 2j * c.pi / n * ((u * i) + (v * j))
                        e = c.exp(angle)
                        tempJ += matrix[i, j] * e
                    tempI += tempJ
                    tempJ = 0

                finalMatrix[u, v] = round((tempI.real / (n ** 2)))
                tempI = 0
        return finalMatrix





    def discrete_cosine_tranform(self, matrix):

        """Computes the discrete cosine transform of the input matrix
        takes as input:
        matrix: a 2d matrix
        returns a matrix representing discrete cosine transform"""
        matrixDimensions = matrix.shape
        finalMatrix = np.zeros((matrixDimensions[0], matrixDimensions[1]), np.float64)
        n = matrixDimensions[0]
        tempI = float(0)
        tempJ = float(0)

        for u in range(n):
            for v in range(n):
                for i in range(n):
                    for j in range(n):
                        if u == 0:
                            a = 1 / float(np.sqrt(n))
                        else:
                            a = float(np.sqrt(2 / n))
                        if v == 0:
                            b = 1 / float(np.sqrt(n))
                        else:
                            b = float(np.sqrt(2 / n))
                        angle1 = math.cos(c.pi / float(n) * ((float(i + .5))) * float(u))
                        angle2 = math.cos(c.pi / float(n) * ((float(j + .5))) * float(v))
                        e = angle1 * angle2
                        tempJ += matrix[i, j] * e * a * b
                    tempI += tempJ
                    tempJ = 0
                finalMatrix[u, v] = tempI
                tempI = 0

        return finalMatrix


    def magnitude(self, matrix):
        """Computes the magnitude of the DFT
        takes as input:
        matrix: a 2d matrix
        returns a matrix representing magnitude of the dft"""

        matrixDimensions = matrix.shape
        n =matrixDimensions[0]
        dftMatrix = self.forward_transform(matrix)
        finalMatrix = np.zeros((matrixDimensions[0], matrixDimensions[1]), np.complex64)

        for u in range(n):
            for v in range(n):
                a = dftMatrix[u, v].real
                b = dftMatrix[u, v].imag
                cm = ((a * a) + (b * b))
                sqrtcm = float(np.sqrt(cm))
                finalMatrix[u, v] = sqrtcm
        return finalMatrix



