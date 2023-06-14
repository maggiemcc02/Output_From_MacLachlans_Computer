# import needed stuff

from needed_operations import *

# parameters

x = 10
y = 2
z = 13

a = add_numbers(x, y)
b = sub_numbers(y, z)
c = divide_numbers(x, y)
d = multiply_numbers(y, z)
e = multiply_numbers(add_numbers(10, 2), sub_numbers(z, x))

print(a, b, c, d, e)
print(12, -11, 5, 26, 36)
