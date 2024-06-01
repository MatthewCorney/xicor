# Introduction


# Dependencies
The function requires the following Python packages:

- numpy
- scipy

# Usage
The function takes three parameters:
xvec= x vector or list of values
yvec = secondary vector of values
method = either "permutation" or a "asymptotic", default "asymptotic"
nperm = if using "permutation" then the number of permutations

# Example
```
pip install git+https://github.com/MatthewCorney/xicor.git
```

```
from xicor import xicor

# From two lists/vectors
xvec=[10, 8, 13, 9, 11, 14, 6, 4, 12, 7, 5],
yvec=[8.04, 6.95, 7.58, 8.81, 8.33, 9.96, 7.24, 4.26, 10.84, 4.82, 5.68]
# Calcuate with the permutation method
result = xicor(xvec=xvec, yvec=yvec, method='permutation', nperm=1000)


# Print the results
print("Test Statistic:", result["xi"])
print("P Value:", result["pval"])
print("Deviation:", result["sd"])
```

```
Test Statistic: 1.0
P Value: 0.0
Deviation: 0.16075038880201817
```
```
# Or with the asymptotic method
result = xicor(xvec=xvec, yvec=yvec, method='asymptotic')


# Print the results
print("Test Statistic:", result["xi"])
print("P Value:", result["pval"])
print("Deviation:", result["sd"])
```