from math import acos, pi, sqrt
from typing import overload, Union, List
from collections import namedtuple


class Vector(object):
    def __init__(self, vec):
        super(Vector, self).__init__()
        self.__vec__: List[float] = vec

    @classmethod
    def undefinedVector(cls):
        return cls([float("nan"), float("nan"), float("nan")])

    def __mul__(self, other: Union['Vector', int, float]) -> Union['Vector', float]:
        if isinstance(other, (int, float)) and (other == other):
            return Vector([other * n for n in self.__vec__])
        elif isinstance(other, Vector):
            return sum([self[i] * other[i] for i in range(len(self))])
        return float('nan')

    def __rmul__(self, other: 'Vector') -> Union['Vector', float]:
        return self.__mul__(other)

    # cross product operator
    def __pow__(self, other: 'Vector') -> 'Vector':
        if isinstance(other, Vector) and len(self) == 3:
            return Vector([1, 1, 1]).crossProduct(self, other)
        return Vector.undefinedVector()

    # cross product function which allows to set up ijk vector
    def crossProduct(self, other1: 'Vector', other2: 'Vector') -> 'Vector':
        if isinstance(other1, Vector) and isinstance(other2, Vector) and len(self) == 3:
            return Vector([self[0] * (other1[1] * other2[2] - other1[2] * other2[1]),
                           -self[1] * (other1[0] * other2[2] - other1[2] * other2[0]),
                           self[2] * (other1[0] * other2[1] - other1[1] * other2[0])])
        return Vector.undefinedVector()

    def __truediv__(self, other: 'Vector') -> Union['Vector', float]:
        if isinstance(other, (int, float)) and (other == other):
            return Vector([n / other for n in self.__vec__])
        elif isinstance(other, Vector):
            return sum([self[i] / other[i] for i in range(len(self))])
        return float('nan')

    def __add__(self, other: 'Vector') -> 'Vector':
        if isinstance(other, Vector):
            return Vector([self[i] + other[i] for i in range(len(self))])
        return Vector.undefinedVector()

    def __radd__(self, other: 'Vector') -> 'Vector':
        return self.__add__(other)

    def __sub__(self, other: 'Vector') -> 'Vector':
        if isinstance(other, Vector):
            return self + (-other)
        return Vector.undefinedVector()

    def __rsub__(self, other: 'Vector') -> 'Vector':
        return self.__sub__(other)

    def __neg__(self) -> 'Vector':
        return self * (-1)

    def __repr__(self):
        return f"V({self.__vec__})"

    def __str__(self):
        return str(self.__vec__)

    def __len__(self):
        return len(self.__vec__)

    def __getitem__(self, key):
        return self.__vec__[key]

    def __setitem__(self, key, value):
        self.__vec__[key] = value
        return value

    def append(self, key):
        self.__vec__.append(key)
        return self

    def angleWith(self, other: 'Vector') -> float:
        if isinstance(other, Vector):
            if self * other == 0:
                return pi / 2
            a = acos((self * other) / (self.length() * other.length()))
            return a if a > 0 else pi - a

    def length(self) -> float:
        return sqrt(self * self)

    def projection(self, other: 'Vector') -> 'Vector':
        return self * other * other / (other * other)

    @overload
    def distance(self, plane: 'Vector') -> float:
        normal = Vector(plane[:-1])
        return (self * normal + plane[-1]) / normal.length()

    def distance(self, point: 'Vector', vector: 'Vector') -> float:
        return ((point - self) ** vector).length() / vector.length()

    def isDefined(self) -> bool:
        return all([self[i] == self[i] for i in range(len(self))])

    def isInRectangle(self, *args) -> bool:
        l1 = args[1] - args[0]
        l2 = args[1] - args[2]
        d1 = max([self.distance(args[0], l1), self.distance(args[2], l1)]) <= l2.length()
        d2 = max([self.distance(args[0], l2), self.distance(args[2], l2)]) <= l1.length()
        return d1 and d2

    def reflect(self, plane: 'Vector') -> 'Vector':
        normal = Vector(plane[:-1])
        return self - self.projection(normal) * 2


def createPlane(a: Vector, b: Vector, c: Vector) -> Vector:
    return ((b - a) ** (c - a)).append(sum(-a.crossProduct(b - a, c - a)))


def getIntersection(M: Vector, directionVector: Vector, P: Vector) -> Vector:
    normal = Vector(P[:-1])
    if normal * directionVector == 0:
        return Vector([float("nan"), float("nan"), float("nan")])
    return M - directionVector * ((normal * M + P[-1]) / (normal * directionVector))


Collision = namedtuple('Collision', ['plane', 'intersection', 'rectangle'])

with open("input.txt", "r") as file:
    A = Vector([float(i) for i in file.readline().split()])
    B = Vector([float(i) for i in file.readline().split()])
    C = Vector([float(i) for i in file.readline().split()])
    D = Vector([float(i) for i in file.readline().split()])
    E = B + (A - B) + (C - B)
    F = C + (D - C) + (B - C)
    G = C + (D - C) + (E - C)

    cubePlanes: List[Vector] = [
        createPlane(A, B, C),
        createPlane(B, C, D),
        createPlane(C, D, E),
        createPlane(A, E, G),
        createPlane(G, D, F),
        createPlane(A, B, F)
    ]

    V = Vector([float(i) for i in file.readline().split()])
    V0 = Vector([float(i) for i in file.readline().split()])

    energy = int(file.readline())
    n = int(file.readline())

    planes: List[Vector] = []
    mirrors: List[List[Vector]] = []

    for i in range(n):
        x1 = Vector([float(i) for i in file.readline().split()])
        x2 = Vector([float(i) for i in file.readline().split()])
        x3 = Vector([float(i) for i in file.readline().split()])
        x4 = x2 + (x1 - x2) + (x3 - x2)
        planes.append(createPlane(x1, x2, x3))
        mirrors.append([x1, x2, x3, x4])

    lastMirror = Vector.undefinedVector()

    while energy > 0:
        mirrorCollisions: List[Collision] = [
            Collision(
                createPlane(*mirror[0:3]),
                getIntersection(V0, V, createPlane(*mirror[0:3])),
                mirror
            ) for mirror in mirrors
        ]

        cubeCollisions: List[Collision] = [
            Collision(
                cubePlane,
                getIntersection(V0, V, cubePlane),
                Vector.undefinedVector()
            ) for cubePlane in cubePlanes
        ]

        predicate = lambda x: x.intersection.isDefined() and \
                              (V * (x.intersection - V0) > 0) and \
                              (x.intersection.isInRectangle(*x.rectangle))


        cubeCollisions = [*filter(predicate, cubeCollisions)]
        mirrorCollisions = [*filter(predicate, mirrorCollisions)]

        if mirrorCollisions:
            closestMirror = min([(i.intersection - V0).length() for i in mirrorCollisions])
            closestEdge = min([(i.intersection - V0).length() for i in cubeCollisions])
            closestPlane = min(closestEdge, closestMirror)

            if (closestEdge == closestPlane):
                V0 = cubeCollisions[0].intersection
                break

            collision = [i for i in mirrorCollisions if (i.intersection - V0).length() == closestMirror][0]
            V0 = collision.intersection
            V = V.reflect(collision.plane)
            energy -= 1
        else:
            cubeCollisions: List[Collision] = [
                Collision(
                    cubePlane,
                    getIntersection(V0, V, cubePlane),
                    Vector.undefinedVector()
                ) for cubePlane in cubePlanes
            ]

            cubeCollisions = [*filter(lambda x: (x.intersection - V0) * V > 0, cubeCollisions)]
            closestPlane = min([(i.intersection - V0).length() for i in cubeCollisions])
            cubeCollisions = [i for i in cubeCollisions if (i.intersection - V0).length() == closestPlane]
            V0 = cubeCollisions[0].intersection
            break

with open("output.txt", "w") as result:
    if energy == 0:
        result.write("0\n")
        result.write(' '.join([str(i) for i in V0]) + "\n")
    else:
        result.write("1\n")
        result.write(str(energy) + "\n")
        result.write(' '.join([str(i) for i in V0]) + "\n")
        result.write(' '.join([str(i) for i in V]) + "\n")
