import numpy as np

def rotator_to_quaternion(mat):
    trace = np.trace(mat)
    if trace > 0.0:
        s = np.sqrt(trace + 1.0) * 2.0
        qw = 0.25 * s
        qx = (mat[1, 2] - mat[2, 1]) / s
        qy = (mat[2, 0] - mat[0, 2]) / s
        qz = (mat[0, 1] - mat[1, 0]) / s
    elif mat[0, 0] > mat[1, 1] and mat[0, 0] > mat[2, 2]:
        s = np.sqrt(1.0 + mat[0, 0] - mat[1, 1] - mat[2, 2]) * 2.0
        qw = (mat[1, 2] - mat[2, 1]) / s
        qx = 0.25 * s
        qy = (mat[0, 1] + mat[1, 0]) / s
        qz = (mat[2, 0] + mat[0, 2]) / s
    elif mat[1, 1] > mat[2, 2]:
        s = np.sqrt(1.0 + mat[1, 1] - mat[0, 0] - mat[2, 2]) * 2.0
        qw = (mat[2, 0] - mat[0, 2]) / s
        qx = (mat[0, 1] + mat[1, 0]) / s
        qy = 0.25 * s
        qz = (mat[1, 2] + mat[2, 1]) / s
    else:
        s = np.sqrt(1.0 + mat[2, 2] - mat[0, 0] - mat[1, 1]) * 2.0
        qw = (mat[0, 1] - mat[1, 0]) / s
        qx = (mat[2, 0] + mat[0, 2]) / s
        qy = (mat[1, 2] + mat[2, 1]) / s
        qz = 0.25 * s
    return (qw, qx, qy, qz)

def quaternion_multiply(q1, q2):
    w = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3]
    x = q1[0] * q2[1] + q1[1] * q2[0] - q1[2] * q2[3] + q1[3] * q2[2]
    y = q1[0] * q2[2] + q1[1] * q2[3] + q1[2] * q2[0] - q1[3] * q2[1]
    z = q1[0] * q1[3] - q1[1] * q2[2] + q1[2] * q2[1] + q1[3] * q2[0]
    return (w, x, y, z)

if __name__ == "__main__":
    A = np.array((0., 0., 0.))
    B = np.array((0., 3., 0.))
    position = (A + B) / 2.0
    vec1 = A - B
    norm = np.linalg.norm(vec1)
    vec1 /= norm
    vec2 = np.array([0, 0, 1])
    if (vec1 == vec2).all():
        rotator = np.identity(3)
    elif (vec1 == -vec2).all():
        rotator = np.identity(3) * -1.0
        rotator[0, 0] += 2.0
    else:
        vec3 = np.cross(vec1, vec2)
        dot = np.dot(vec1, vec2)
        skew = np.array([[0, -vec3[2], vec3[1]],
                        [vec3[2], 0, -vec3[0]],
                        [-vec3[1], vec3[0], 0]])
        rotator = np.identity(3) + skew + (np.matmul(skew, skew) * (1.0 / (1.0 + dot)))

    print(rotator)
    
    quat = rotator_to_quaternion(rotator)
    print(quat)

    meta_mod = 1 / (np.sqrt(2))
    quat = quaternion_multiply(quat, (meta_mod, 0., meta_mod, 0.))
    print(quat)