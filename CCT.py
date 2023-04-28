###########################################
# Functions to simplify chosen ciphertext
# packing for Kyber768 and Saber. Partly
# based on work by Kalle Ngo - ngo@kth.se - KTH
#
# Linus Backlund - lbackl@kth.se - KTH
# 22-08-2022
###########################################


from struct import pack

def __pack_b_part(b_part):
    packed = []
    for i in range(0, 256, 4):
        packed.append(b_part[i] & 0xFF)
        packed.append((b_part[i] >> 8) + ((b_part[i+1] & 0x3F) << 2))
        packed.append((b_part[i+1] >> 6) + ((b_part[i+2] & 0x0F) << 4))
        packed.append((b_part[i+2] >> 4) + ((b_part[i+3] & 0x03) << 6))
        packed.append(b_part[i+3] >> 2)

    return pack('<320B', *packed)


def __pack_cm(cm):
    packed = []
    for i in range(0, 256, 2):
        packed.append(cm[i] + (cm[i+1] << 4))
    
    return pack('<128B', *packed)


def __pack_b_and_cm(b, cm):
    ct = bytearray()
    for b_part in b:
        ct.extend(__pack_b_part(b_part))
    ct.extend(__pack_cm(cm))
    
    return ct


# Pack a ciphertext using a CT_set on the form (k1, k0) or (k1, k0, k2)
# If using a k2, k2 indexes must also be provided
def CCT(CT_set, part, rotation=0, k2_indexes=None):
    # Unpack CT_set
    if len(CT_set) == 2:
        k1, k0 = CT_set
        k2 = None
    elif len(CT_set) == 3:
        k1, k0, k2 = CT_set

    # Construct ciphertext polynomials
    b = [[0]*256, [0]*256, [0]*256]
    b[part][rotation] = k1
    cm = [k0]*256
    
    # If k2 provided, insert it into its assigned indexes
    if not k2 is None and not k2_indexes is None:
        for k2_index in k2_indexes:
            cm[k2_index] = k2

    # Pack polynomials into ciphertext
    return __pack_b_and_cm(b, cm)
    

if __name__ == '__main__':
    ct = CCT((1, 0, 1), part=0, rotation=1, k2_indexes=[253,254,255])
    print(ct)
