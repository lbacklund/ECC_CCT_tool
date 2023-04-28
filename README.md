# ECC_CCT_tool
Tools for working with chosen ciphertexts (CCT) based on error-correcting codes. (ECC)

## Credit
The ciphertext packing is a direct python-port of the one imlemented by Kalle Ngo.
The rotation technique is the work of Joel Gärtner.

## Behind the scenes
This tool takes care of a lot of gnarly stuff. Here is a list of the most notable ones:
- CCT-sets for linear codes of several code distances
- CCT-packing
- CCT rotation 
- Message hamming weight minimization
- Deriving secret key from recovered messages
- Error correction handling
- Conversion from Saber and Kyber secret key formats to coefficients
- Recovery statistics collection
- Incremental codes, where each code is a subset of the ones with higher code distance

## Supported algorithms
- saber
- kyber768

## Supported code distances
- 2 (only kyber768)
- 3
- 4
- 5
- 6
- 7 (only kyber768)
- 8 (only kyber768)

## Basic usage
Clone repository into you working directory.
```console
git clone https://gits-15.sys.kth.se/hardware-security-group/ECC_CCT_tool.git
```

Import module.
```python
from ECC_CCT_tool.ECC_CCT_tool import ECC_CCT_TOOL
```

Initialize an ECC_CCT_tool object.
```python
ECC_tool = ECC_CCT_tool("kyber768", code_distance=4)
```

Generate ciphertexts
```python
ct = ECC_tool.CCT(CT_set, part, rotation, k2_indexes)
```
where `CT_set` is either a tuple (k1, k0) or, if message hamming-weight minimization is desired, (k1, k0, k2). `part` is an integer 0 to 2. `rotation` (optional) is an integer 0 to 255. If message hamming weight minimization is desired, `k2_indexes` (optional) can be provided as a list of the messeage bit indexes to minimize.

Recover the messages for each part and CT as a `numpy` array and store them
```python
recovered_messages = np.zeros((3, len(ECC_tool.ct_table), 256), dtype=int)

for part in range(3):
  for CT_num, CT_set in enumerate(ECC_tool.ct_table):
    # Your code to recover and store the message
    # into recovered_messages[part, CT_num]
```
containing all the messages for each CT-set in ECC_tool.ct_table

Predict the secret key from recovered messages.
```python
ECC_tool.predict_secret_key(recovered_messages)
```

Load the true secret key from a file
```python
sk = bytearray()
sk_filename = "GNDTruth/SecKey.bin"
with open(sk_filename, "rb") as f:
    byte = f.read(1)
    while byte != b"":
        sk.extend(byte)
        temp = int.from_bytes(byte,byteorder='little')
        byte = f.read(1)
```
or, if you are using the pq_kem tool, just use the secret key generated.

Compare predicted key against true secret key.
```python
ECC_tool.compare_against_true_secret_key(sk)
```
This will provide a nice printout of the recovered coefficients and statistics.

If you want the statistics dictionary for saving to a file it can be accessed through
```python
ECC_tool.stats
```

Should you want to map your recovered traces fit the code of a lower code distance. You can use the provided index-array.
```python
recovered_messages = recovered_messages[:, ECC_tool.map_lower_cd]
ECC_tool = ECC_CCT_TOOL(ECC_tool.algorithm, ECC_tool.cd-1)
```
