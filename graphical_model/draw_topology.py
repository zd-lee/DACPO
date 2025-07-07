import numpy as np
'''
    Parameters:
    - mesh: A mesh represented as a tuple (vertices, faces), where vertices is a list of 3D points and faces is a list of vertex indices for each face
    - topology: A mesh represented as a tuple (vertices, faces) to add to the mesh
    Returns:
    - None
'''
def add_topology(mesh,topology):
    start = len(mesh[0])
    for i in topology[0]:
        mesh[0].append(i)
    for i in topology[1]:
        t = []
        for j in i:
            t.append(j + start)
        mesh[1].append(t)
    return 0

def normalize(v):
    norm = np.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    v[0] /= norm
    v[1] /= norm
    v[2] /= norm

'''
    Parameters:
    - start: The starting point of the arrow
    - end: The ending point of the arrow
    - radius: The radius of the arrow
    - n: the longtitude precision of the arrow
    Returns:
    - A mesh representing the arrow
'''
def get_arrow(start, end, radius = 0.01, n = 10):
    rate = 2
    dir = end - start
    z = dir.copy()
    normalize(z)
    x = np.array([1, 0, 0])
    if np.linalg.norm(x - z) < 0.01:
        x = np.array([0, 1, 0])
    y = np.cross(z, x)
    normalize(y)
    x = np.cross(y, z)
    normalize(x)

    end = start + dir * 0.95
    start = start + dir * 0.05
    dir = end - start

    cylinder_end = start + dir * 0.9

    res = ([], [])
    for i in range(n):
        theta = 2 * np.pi * i / n
        p = start + radius / rate * (np.cos(theta) * x + np.sin(theta) * y)
        res[0].append(p)

    for i in range(n):
        theta = 2 * np.pi * i / n
        p = cylinder_end + radius / rate * (np.cos(theta) * x + np.sin(theta) * y)
        res[0].append(p)

    for i in range(n):
        res[1].append([n + i, i, (i + 1) % n])
        res[1].append([n + i, (i + 1) % n, (i + 1) % n + n])
    cend = len(res[0])

    for i in range(n):
        theta = 2 * np.pi * i / n
        p = cylinder_end + radius * (np.cos(theta) * x + np.sin(theta) * y)
        res[0].append(p)

    res[0].append(end)
    top = len(res[0]) - 1

    for i in range(n):
        res[1].append([top, i + cend, (i + 1) % n + cend])
    return res

'''
    Parameters:
    - center: The center of the sphere
    - radius: The radius of the sphere
    - n: The latitude precision of the sphere
    - m: The longitude precision of the sphere 
    Returns:
    - A mesh representing the sphere
'''
def get_sphere(center, radius = 0.08, n = 10, m = 10):
    res = ([], [])
    for i in range(n):
        for j in range(m):
            theta = 2 * np.pi * i / n
            phi = np.pi * j / m
            p = center + np.array([radius * np.sin(phi) * np.cos(theta), radius * np.sin(phi) * np.sin(theta), radius * np.cos(phi)])
            res[0].append(p)
    for i in range(n):
        for j in range(m):
            a = i * m + j
            b = i * m + (j + 1) % m
            c = ((i + 1) % n) * m + j
            d = ((i + 1) % n) * m + (j + 1) % m
            res[1].append([a, b, c])
            res[1].append([b, d, c])
    return res