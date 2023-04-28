###########################################
# A tool for chosen ciphertexts using error
# correcting codes for PQ-algorithms Saber
# and Kyber.
#
# Linus Backlund - lbackl@kth.se - KTH
# 22-08-2022
###########################################


import numpy as np

from .CCT import CCT


class UncorrectableError(Exception):
    def __init__(self, message):            
        super().__init__(f"Uncorrectable error! ({message})")


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


class ECC_CCT_TOOL():
    def __init__(self, algorithm, code_distance):
        self.algorithm = algorithm
        self.cd = code_distance
        self.stats = {}
        self.DEBUG = False
        
        # Perform sanity-check on parameters
        if self.algorithm == "kyber768":
            if self.cd not in [2,3,4,5,6,7,8]:
                raise ValueError(f"Code distance {self.cd} not available for {self.algorithm}.")
        elif self.algorithm == "saber":
            if self.cd not in [3,4,5,6]:
                raise ValueError(f"Code distance {self.cd} not available for {self.algorithm}.")
        else:
            raise ValueError(f"Algorithm {self.algorithm} unknown.")
    
        
        if self.algorithm == "saber":
            self.data_length = 4
            if self.cd == 3:
                self.ct_table = [
                    [373, 15, 12],
                    [793, 10, 3],
                    [913, 15, 1],
                    [613, 4, 2],
                    [68, 4, 14],
                    [240, 10, 12],
                    [12, 3, 0],
                ]
                self.H = [
                    [1, 0, 1, 1, 1, 0, 0],
                    [1, 1, 1, 0, 0, 1, 0],
                    [1, 1, 0, 1, 0, 0, 1],
                ]
                self.coefficient_map = [15, 2, 11, 12, 5, 9, 0, 7, 14]
            elif self.cd == 4:
                self.ct_table = [
                    [613, 4, 2],
                    [68, 4, 14],
                    [373, 15, 12],
                    [793, 10, 3],
                    [456, 15, 5],
                    [913, 15, 1],
                    [240, 10, 12],
                    [12, 3, 0],
                ]
                self.H = [
                    [0, 1, 1, 1, 1, 0, 0, 0],
                    [1, 1, 1, 0, 0, 1, 0, 0],
                    [1, 1, 0, 1, 0, 0, 1, 0],
                    [1, 0, 1, 1, 0, 0, 0, 1],
                ]
                self.coefficient_map = [15, 4, 14, 7, 13, 10, 0, 9, 3]
            elif self.cd == 5:
                self.ct_table = [
                    [240, 10, 12],
                    [373, 15, 12],
                    [793, 10, 3],
                    [68, 4, 14],
                    [377, 10, 12],
                    [806, 10, 3],
                    [613, 4, 2],
                    [456, 15, 5],
                    [917, 10, 1],
                    [913, 15, 1],
                    [12, 3, 0],
                ]
                self.H = [
                    [0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0],
                    [1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0],
                    [1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0],
                    [0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0],
                    [1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                    [1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0],
                    [1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1],
                ]
                self.coefficient_map = [15, 9, 5, 7, 11, 12, 0, 2, 14]
            elif self.cd == 6:
                self.ct_table = [
                    [240, 10, 12],
                    [377, 10, 12],
                    [613, 4, 2],
                    [373, 15, 12],
                    [913, 15, 1],
                    [12, 3, 0],
                    [793, 10, 3],
                    [755, 4, 0],
                    [917, 10, 1],
                    [806, 10, 3],
                    [456, 15, 5],
                    [68, 4, 14],
                ]
                self.H = [
                    [1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                    [0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                    [0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0],
                    [1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                    [1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0],
                    [0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0],
                    [1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0],
                    [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1],
                ]
                self.coefficient_map = [11, 8, 7, 1, 14, 15, 0, 6, 9]
        elif self.algorithm == "kyber768":
            self.data_length = 3
            if self.cd == 2:
                self.ct_table = [
                    [335, 15, 2],
                    [77, 5, 0],
                    [432, 8, 3],
                    [606, 3, 3],
                ]
                self.H = [
                    [1, 1, 1, 1],
                ]
                self.coefficient_map = [7, 6, 3, 4, 5]
            elif self.cd == 3:
                self.ct_table = [
                    [153, 8, 1],
                    [335, 15, 2],
                    [606, 3, 3],
                    [432, 8, 3],
                    [77, 5, 0],
                    [432, 3, 3],
                ]
                self.H = [
                    [1, 1, 0, 1, 0, 0],
                    [1, 0, 1, 0, 1, 0],
                    [0, 1, 1, 0, 0, 1],
                ]
                self.coefficient_map = [3, 6, 4, 7, 2]
                self.map_lower_cd =  [1, 4, 3, 2]
            elif self.cd == 4:
                self.ct_table = [
                    [432, 8, 3],
                    [432, 3, 3],
                    [915, 5, 0],
                    [606, 3, 3],
                    [335, 15, 2],
                    [153, 8, 1],
                    [77, 5, 0],
                ]
                self.H = [
                    [1, 0, 1, 1, 0, 0, 0],
                    [1, 1, 1, 0, 1, 0, 0],
                    [0, 1, 1, 0, 0, 1, 0],
                    [1, 1, 0, 0, 0, 0, 1],
                ]
                self.coefficient_map = [4, 2, 5, 1, 7]
                self.map_lower_cd =  [5, 4, 3, 0, 6, 1]
            elif self.cd == 5:
                self.ct_table = [
                    [77, 5, 0],
                    [432, 3, 3],
                    [606, 3, 3],
                    [432, 8, 3],
                    [915, 5, 0],
                    [335, 15, 2],
                    [321, 15, 2],
                    [632, 3, 3],
                    [153, 8, 1],
                    [864, 8, 1],
                ]
                self.H = [
                    [1, 1, 0, 1, 0, 0, 0, 0, 0, 0],
                    [1, 1, 1, 0, 1, 0, 0, 0, 0, 0],
                    [0, 1, 1, 0, 0, 1, 0, 0, 0, 0],
                    [0, 1, 1, 0, 0, 0, 1, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0, 0, 1, 0, 0],
                    [1, 0, 1, 0, 0, 0, 0, 0, 1, 0],
                    [1, 0, 1, 0, 0, 0, 0, 0, 0, 1],
                ]
                self.coefficient_map = [5, 6, 4, 1, 2]
                self.map_lower_cd =  [3, 1, 4, 2, 5, 8, 0]
            elif self.cd == 6:
                self.ct_table = [
                    [335, 15, 2],
                    [153, 8, 1],
                    [606, 3, 3],
                    [915, 5, 0],
                    [898, 5, 0],
                    [864, 8, 1],
                    [632, 3, 3],
                    [77, 5, 0],
                    [432, 8, 3],
                    [321, 15, 2],
                    [432, 3, 3],
                ]
                self.H = [
                    [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
                    [1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0],
                    [0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0],
                    [1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                    [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                    [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1],
                ]
                self.coefficient_map = [5, 6, 2, 7, 4]
                self.map_lower_cd =  [7, 10, 2, 8, 3, 0, 9, 6, 1, 5]
            elif self.cd == 7:
                self.ct_table = [
                    [153, 8, 1],
                    [606, 3, 3],
                    [915, 5, 0],
                    [898, 5, 0],
                    [335, 15, 2],
                    [321, 15, 2],
                    [432, 3, 3],
                    [864, 8, 1],
                    [632, 3, 3],
                    [77, 5, 0],
                    [432, 8, 3],
                    [105, 5, 0],
                    [606, 8, 3],
                ]
                self.H = [
                    [0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                    [1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                    [1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                    [1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                    [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                    [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                    [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                ]
                self.coefficient_map = [2, 4, 5, 7, 1]
                self.map_lower_cd =  [4, 0, 1, 2, 3, 7, 8, 9, 10, 5, 6]
            elif self.cd == 8:
                self.ct_table = [
                    [153, 8, 1],
                    [77, 5, 0],
                    [335, 15, 2],
                    [432, 3, 3],
                    [606, 3, 3],
                    [864, 8, 1],
                    [915, 5, 0],
                    [898, 5, 0],
                    [432, 8, 3],
                    [105, 5, 0],
                    [632, 3, 3],
                    [386, 3, 3],
                    [606, 8, 3],
                    [321, 15, 2],
                ]
                self.H = [
                    [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                    [0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                    [1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                    [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                    [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                    [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                ]
                self.coefficient_map = [3, 7, 6, 5, 1]
                self.map_lower_cd = [0, 4, 6, 7, 2, 13, 3, 5, 10, 1, 8, 9, 12]

        self.__nttZetasInv = [
            1701, 1807, 1460, 2371, 2338, 2333, 308, 108, 2851, 870, 854, 1510, 2535,
            1278, 1530, 1185, 1659, 1187, 3109, 874, 1335, 2111, 136, 1215, 2945, 1465,
            1285, 2007, 2719, 2726, 2232, 2512, 75, 156, 3000, 2911, 2980, 872, 2685,
            1590, 2210, 602, 1846, 777, 147, 2170, 2551, 246, 1676, 1755, 460, 291, 235,
            3152, 2742, 2907, 3224, 1779, 2458, 1251, 2486, 2774, 2899, 1103, 1275, 2652,
            1065, 2881, 725, 1508, 2368, 398, 951, 247, 1421, 3222, 2499, 271, 90, 853,
            1860, 3203, 1162, 1618, 666, 320, 8, 2813, 1544, 282, 1838, 1293, 2314, 552,
            2677, 2106, 1571, 205, 2918, 1542, 2721, 2597, 2312, 681, 130, 1602, 1871,
            829, 2946, 3065, 1325, 2756, 1861, 1474, 1202, 2367, 3147, 1752, 2707, 171,
            3127, 3042, 1907, 1836, 1517, 359, 758, 1441
            ]

        self.secret_coefficients = np.zeros((768,), dtype=int)

       
    def CCT(self, CT_set, part, rotation=0, k2_indexes=None):
        return CCT(CT_set, part, rotation, k2_indexes)


    def __barrett(self, a):
        v = ((1<<24) + 3329//2) // 3329
        t = v * a >> 24
        t = t * 3329

        return a - t


    def __int16(self, n):
        end = -32768
        start = 32767

        if n >= end and n <= start:
            return n
        elif n < end:
            return ((n + 32769) % 65536) + start
        else:
            return ((n - 32768) % 65536) + end


    def __int32(self, n):
        end = -2147483648
        start = 2147483647

        if n >= end and n <= start:
            return n
        elif n < end:
            return ((n + 2147483649) % 4294967296) + start
        else:
            return ((n - 2147483648) % 4294967296) + end


    def __byteopsMontgomeryReduce(self, a):
        u = self.__int16(self.__int32(a) * 62209)
        t = u * 3329
        t = a - t
        t = t >> 16

        return self.__int16(t)


    def __nttFqMul(self, a, b):
        return self.__byteopsMontgomeryReduce(a * b)


    def __inv_ntt(self, r):
        j = 0
        k = 0
        for l in [2,4,8,16,32,64,128]:
            start = 0
            while start < 256:
                zeta = self.__nttZetasInv[k]
                k += 1
                j = start
                while j < start + l:
                    t = r[j]
                    r[j] = self.__barrett(t + r[j + l])
                    r[j + l] = t - r[j + l]
                    r[j + l] = self.__nttFqMul(zeta, r[j + l])
                    j += 1
                start = j + l
        for j in range(256):
            r[j] = self.__nttFqMul(r[j], self.__nttZetasInv[127])

        return r
    

    def secret_bytes_to_coefficients(self, sk):

        if self.algorithm == "kyber768":
            coeff_len = 12
        elif self.algorithm == "saber":
            coeff_len = 13

        # Convert bytearray into bitstring
        bins = [int(bit) for bit in ''.join([f"{byte:08b}"[::-1] for byte in sk[:int(768*coeff_len/8)]])]

        # Cut bitstring into 12-bit segments
        binstrings = []
        for i in range(int(768)):
            binstrings.append("".join([str(bit) for bit in bins[i*coeff_len:(i+1)*coeff_len]]))

        # Convert bitstrings to integers
        coefficients = []
        if self.algorithm == "kyber768":
            for i in range(768):
                coefficients.append(int(binstrings[i][::-1], 2))

            # Perform inverse NTT
            inv1 = self.__inv_ntt(coefficients[:256])
            inv2 = self.__inv_ntt(coefficients[256:512])
            inv3 = self.__inv_ntt(coefficients[512:])
            all_inv = inv1+inv2+inv3

            # Remove noise and map to fixed values
            all_inv = np.array([((i + 3*3329) % 3329) // 17 for i in all_inv])

            # Map values to their corresponding coefficient
            all_inv = np.where(all_inv==134, 1, all_inv)
            all_inv = np.where(all_inv==61, -1, all_inv)
            all_inv = np.where(all_inv==73, 2, all_inv)
            all_inv = np.where(all_inv==122, -2, all_inv)

            return all_inv

        elif self.algorithm == "saber":
            for i in range(768):
                coefficients.append(int(binstrings[i][:3][::-1], 2))
                if int(binstrings[i][3]):
                    coefficients[i] -= 8

            return coefficients


    def print_coefficients(self, predicted_coefficients, secret_coefficients=None):
        if secret_coefficients is None:
            print("Predicted coefficients:")
        else:
            print("A posteriori:")
            print(bcolors.OKGREEN+"  GREEN"+bcolors.ENDC+" = SUCCESS")
            print(bcolors.FAIL+"  RED"+bcolors.ENDC+" = FAILURE")
        print(bcolors.WARNING+"  YELLOW"+bcolors.ENDC+" = ENUMERABLE")
        print()

        ascii_coefficients = []
        for i in range(768):
            coefficient = predicted_coefficients[i]
            true_coefficient = None
            if not secret_coefficients is None: true_coefficient = secret_coefficients[i]

            # If error detected, print '?'
            if coefficient == 100 or coefficient == -100:
                ascii_coefficients.append(bcolors.WARNING+bcolors.BOLD+" ?"+bcolors.ENDC)
            else:
                str_coeff = str(coefficient)
                # If true coefficients not provided, print in white
                if true_coefficient is None:
                    ascii_coefficients.append(" "*(2-len(str_coeff))+str_coeff)
                else:
                    # If coefficient is correct, print green, otherwise, print red
                    if coefficient == true_coefficient:
                        ascii_coefficients.append(bcolors.OKGREEN+" "*(2-len(str_coeff))+str_coeff+bcolors.ENDC)
                    else:
                        ascii_coefficients.append(bcolors.BOLD+bcolors.FAIL+" "*(2-len(str_coeff))+str_coeff+bcolors.ENDC)

        for i in range(0,768,32):
            print("\t"+" ".join(ascii_coefficients[i:i+32]))


    def generate_error_dict(self):
        # For each codeword that do not produce a syndrome
        for i in range(2**len(self.H[0])):
            codeword = np.array([int(bit) for bit in f"{i:016b}"[:len(self.H[0])]])  # int to binary array
            if self.DEBUG: print("codeword", codeword)
            syndrome = self.H@codeword%2
            if self.DEBUG: print("syndrome", syndrome)
            if not np.sum(syndrome):
                error_dict = {}
        
                # Build dict of correctable errors
                for error1 in range(len(self.H[0])):
                    error_codeword = np.array(codeword)
                    error_codeword[error1] = (codeword[error1] + 1)%2  # Flip bit #error in codeword
                    syndrome = self.H@error_codeword%2
                    syndrome_int = np.sum([bit*2**i for i, bit in enumerate(syndrome)])  # Binary array to int
                    error_dict[syndrome_int] = error1
                    if self.cd >= 5:
                        for error2 in range(len(self.H[0])):
                            if error1 == error2: continue
                            error_codeword = np.array(codeword)
                            error_codeword[error1] = (codeword[error1] + 1)%2  # Flip bit #error1 in codeword
                            error_codeword[error2] = (codeword[error2] + 1)%2  # Flip bit #error2 in codeword
                            syndrome = self.H@error_codeword%2
                            syndrome_int = np.sum([bit*2**i for i, bit in enumerate(syndrome)])  # Binary array to int
                            error_dict[syndrome_int] = [error1, error2]  # Save map from int representation of syndrome to erroneous bits
                            if self.cd >= 7:
                                for error3 in range(len(self.H[0])):
                                    if error1 == error3: continue
                                    if error2 == error3: continue
                                    error_codeword = np.array(codeword)
                                    error_codeword[error1] = (codeword[error1] + 1)%2  # Flip bit #error1 in codeword
                                    error_codeword[error2] = (codeword[error2] + 1)%2  # Flip bit #error2 in codeword
                                    error_codeword[error3] = (codeword[error3] + 1)%2  # Flip bit #error3 in codeword
                                    syndrome = self.H@error_codeword%2
                                    syndrome_int = np.sum([bit*2**i for i, bit in enumerate(syndrome)])  # Binary array to int
                                    error_dict[syndrome_int] = [error1, error2, error3]  # Save map from int representation of syndrome to erroneous bits

                    # If dict not complete, the codeword was not valid
                    if self.cd >= 7:
                        if len(error_dict.keys()) == (len(self.H[0]) + np.math.factorial(len(self.H[0]))/(np.math.factorial(len(self.H[0])-3)*6)):
                            break
                    elif self.cd >= 5:
                        if len(error_dict.keys()) == (len(self.H[0]) + np.math.factorial(len(self.H[0]))/(np.math.factorial(len(self.H[0])-2)*2)):
                            break
                    else:
                        if len(error_dict.keys()) == len(self.H[0]):
                            break

        return error_dict


    def predict_secret_key(self, recovered_messages):
        # Generate a dictionary to hold statistics
        self.stats = {
                "single_errors_corrected": 0,
                "double_errors_corrected": 0,
                "triple_errors_corrected": 0,
                "detected_failures": 0,
                "undetected_failures": 0
                }

        # Generate a dictionary with known corrections
        error_dict = self.generate_error_dict()

        # Recover one part of the secret key at a time
        for part, part_msgs in enumerate(recovered_messages):
            coefficients = np.zeros((256,))

            # Recover one coefficient at a time
            for coeff_i, pattern in enumerate(part_msgs.T):
                if self.DEBUG: print("pattern =", pattern)

                # Check if an error has occured
                syndrome = self.H@pattern%2
                if self.DEBUG: print("syndrome =", syndrome)
                error = np.sum([bit*2**i for i, bit in enumerate(syndrome)])
                
                try:
                    # If error detected, attempt to correct it
                    if error:
                        if not error in error_dict.keys(): raise UncorrectableError(f"No correction known for error: {error}")
                        error_index = error_dict[error]
                        if type(error_index) == int:
                            pattern[error_index] = (pattern[error_index] + 1)%2
                            if self.DEBUG: print("corrected bit =", len(pattern)-1-error_index)
                            self.stats["single_errors_corrected"] += 1
                        else:
                            for e_index in error_index:
                                pattern[e_index] = (pattern[e_index] + 1)%2
                                if self.DEBUG: print("corrected bit =", len(pattern)-1-e_index)
                            if len(error_index) == 2:
                                self.stats["double_errors_corrected"] += 1
                            elif len(error_index) == 3:
                                self.stats["triple_errors_corrected"] += 1

                    # Check if an error remains after correction
                    syndrome = self.H@pattern%2
                    error = np.sum([bit*2**i for i, bit in enumerate(syndrome)])
                    if error: raise UncorrectableError(f"Still error after correction: {error}")
                    
                    # Map received data to a coefficient
                    if self.DEBUG: print("data =", pattern[:self.data_length][::-1])
                    data = np.sum([bit*2**i for i, bit in enumerate(pattern[:self.data_length][::-1])])
                    if self.DEBUG: print("data =", data)
                    if data not in self.coefficient_map: raise UncorrectableError(f"No coefficient maping for data: {data}")
                    coefficient = self.coefficient_map.index(data)-((len(self.coefficient_map)-1)//2)
                    
                    if self.DEBUG: print("coefficient =", coefficient)

                # Handle the occurence of an uncorrectable error
                except UncorrectableError as e:
                    if self.DEBUG: print(e)
                    coefficient = 100 # This represent an uncorrectable error
                    self.stats["detected_failures"] += 1
                
                coefficients[coeff_i] = coefficient
                if self.DEBUG: print()

            self.secret_coefficients[part*256:(part+1)*256] = coefficients[::-1]


    def compare_against_true_secret_key(self, sk):
        if type(sk) == bytearray:
            true_coefficients = self.secret_bytes_to_coefficients(sk)
        elif type(sk) == np.ndarray:
            true_coefficients = sk
        else:
            raise TypeError('Incorrect type for true secret key')

        # Print the predicted key compared to the true secret key
        self.print_coefficients(self.secret_coefficients, true_coefficients)
        print()

        # Count the undetected failures
        correct = 0
        for b, pb in zip(true_coefficients, self.secret_coefficients):
            if b==pb: correct += 1
        if correct == 768: print(bcolors.OKGREEN, end='')
        else: print(bcolors.FAIL, end='')
        print(bcolors.UNDERLINE+f"{correct}/768{bcolors.ENDC} coefficients correctly predicted")
        self.stats["undetected_failures"] = 768 - correct - self.stats["detected_failures"]

        print()
        print("Corrected single errors    :", self.stats["single_errors_corrected"])
        print("Corrected double errors    :", self.stats["double_errors_corrected"])
        print("Corrected triple errors    :", self.stats["triple_errors_corrected"])
        print("Detected failures          :", self.stats["detected_failures"])
        print("Undetected failures        :", self.stats["undetected_failures"])
