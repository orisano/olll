# oLLL
A Python3 Implementation of LLL.

## Installation
```bash
python3 -m pip install olll
```
or
```bash
curl -O https://raw.githubusercontent.com/orisano/olll/master/olll.py
```

## Example
https://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm#Example

```python
import olll

reduced_basis = olll.reduction([
  [1, 1, 1],
  [-1, 0, 2],
  [3, 5, 6],
], 0.75)

print(reduced_basis)
# [[0, 1, 0], [1, 0, 1], [-1, 0, 2]]
```

## Reference
https://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm
