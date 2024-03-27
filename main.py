import numpy as np

def cubic_coefficients(xi, fi):
    n = len(xi)
    h = np.diff(xi)
    delta = np.diff(fi) / h
    A = np.zeros((n, n))
    B = np.zeros(n)

    A[0, 0] = 1
    A[-1, -1] = 1
    for i in range(1, n - 1):
        A[i, i - 1] = h[i - 1]
        A[i, i] = 2 * (h[i - 1] + h[i])
        A[i, i + 1] = h[i]
        B[i] = 3 * (delta[i] - delta[i - 1])

    c = np.linalg.solve(A, B)

    d = np.zeros(n - 1)
    b = np.zeros(n - 1)
    for i in range(n - 1):
        d[i] = (c[i + 1] - c[i]) / (3 * h[i])
        b[i] = delta[i] - h[i] * (2 * c[i] + c[i + 1]) / 3

    return b, c, d

def cubic_spline(x, xi, fi, b, c, d):
    i = np.searchsorted(xi, x) - 1
    if i == len(xi) - 1:
        i -= 1
    h = x - xi[i]
    return fi[i] + b[i] * h + c[i] * h**2 + d[i] * h**3 if i < len(d) else fi[i] + b[i] * h + c[i] * h**2

xi = np.array([-0.5, -0.2, 0.1, 0.4, 0.7])
fi = np.array([-0.81152, -0.20017, 0.40136, 1.0236, 1.7273])

b, c, d = cubic_coefficients(xi, fi)

X_star = 1
result = cubic_spline(X_star, xi, fi, b, c, d)
print("Значення функції в точці X* = 1:", result)