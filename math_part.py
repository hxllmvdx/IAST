from math import sqrt, sin, cos
from itertools import count


def get_curr_mass(t):
    if 0 <= t <= 115:
        return 389052 - 1673.348 * t
    elif 115 < t <= 262:
        return 162819 - 756.0544 * (t - 115)
    elif 262 < t <= 467:
        return 46246 - 129.439 * (t - 262)
    elif 467 < t <= 911:
        return 16258 - 30.691 * (t - 467)


def calc_ax(x, y, z, vx, vy, vz, m, k):
    return (-2 * m * vz * 0.879 + k * vx / sqrt(vx**2 + vy**2 + vz**2) - 6.67 * 10**(-11) * (m * 5.97 * 10**24 * x) /
          sqrt(x**2 + y**2 + z**2)**3 - 0.001 * vx) / m


def calc_ay(x, y, z, vx, vy, vz, m, k):
    return (2 * m * vz * 7.29 * 10 ** -5 * sin(28.5) + k * vy / sqrt(vx**2 + vy**2 + vz**2) - 6.67 * 10**(-11) *
          (m * 5.97 * 10 ** 24 * y) / sqrt(x ** 2 + y ** 2 + z ** 2)**3 - 0.001 * vy) / m


def calc_az(x, y, z, vx, vy, vz, m, k):
    return (2 * m * 7.29 * 10**-5 * vx * 0.879 - 2 * m * vy * 7.29**-5 * 0.477 + k * vz / sqrt(vx**2 + vy**2 + vz**2) -
            6.67 * 10**(-11) * (m * 5.97 * 10**24 * z) / sqrt(x**2 + y**2 + z**2)**3 - 0.001 * vz) / m

def runge_kutta(x_0, y_0, z_0, vx_0, vy_0, vz_0, ax_0, ay_0, az_0, m, k, dt):
    data = {}

    x = x_0 + dt / 2 * vx_0 + dt ** 2 / 8 * ax_0
    vx = vx_0 + dt / 2 * ax_0
    y = y_0 + dt / 2 * vy_0 + dt ** 2 / 8 * ay_0
    vy = vy_0 + dt / 2 * ay_0
    z = z_0 + dt / 2 * vz_0 + dt ** 2 / 8 * az_0
    vz = vz_0 + dt / 2 * az_0

    ax = calc_ax(x, y, z, vx, vy, vz, m, k)
    ay = calc_ay(x, y, z, vx, vy, vz, m, k)
    az = calc_az(x, y, z, vx, vy, vz, m, k)

    vx_05 = vx_0 + dt / 2 * ax
    vy_05 = vy_0 + dt / 2 * ay
    vz_05 = vz_0 + dt / 2 * az
    x_05 = x_0 + dt / 2 * vx_05 + dt ** 2 / 8 * ax
    y_05 = y_0 + dt / 2 * vy_05 + dt ** 2 / 8 * ay
    z_05 = z_0 + dt / 2 * vz_05 + dt ** 2 / 8 * az

    ax_05 = calc_ax(x_05, y_05, z_05, vx_05, vy_05, vx_05, m, k)
    ay_05 = calc_ay(x_05, y_05, z_05, vx_05, vy_05, vx_05, m, k)
    az_05 = calc_az(x_05, y_05, z_05, vx_05, vy_05, vx_05, m, k)

    vx_ = vx_0 + dt * ax_05
    vy_ = vy_0 + dt * ay_05
    vz_ = vz_0 + dt * az_05
    x_ = x_0 + dt * vx_ + dt**2 / 2 * ax_05
    y_ = y_0 + dt * vy_ + dt**2 / 2 * ay_05
    z_ = z_0 + dt * vz_ + dt**2 / 2 * az_05

    ax_ = calc_ax(x_, y_, z_, vx_, vy_, vz_, m, k)
    ay_ = calc_ay(x_, y_, z_, vx_, vy_, vz_, m, k)
    az_ = calc_az(x_, y_, z_, vx_, vy_, vz_, m, k)

    vx_1 = vx_0 + dt / 6 * (ax_0 + 2 * ax + 2 * ax_05 + ax_)
    vy_1 = vy_0 + dt / 6 * (ay_0 + 2 * ay + 2 * ay_05 + ay_)
    vz_1 = vz_0 + dt / 6 * (az_0 + 2 * az + 2 * az_05 + az_)
    x_1 = x_0 + dt * vx_1 + dt**2 / 12 * (ax_0 + 2 * ax + 2 * ax_05 + ax_)
    y_1 = y_0 + dt * vy_1 + dt**2 / 12 * (ay_0 + 2 * ay + 2 * ay_05 + ay_)
    z_1 = z_0 + dt * vz_1 + dt**2 / 12 * (az_0 + 2 * az + 2 * az_05 + az_)

    ax_1 = calc_ax(x_1, y_1, z_1, vx_1, vy_1, vz_1, m, k)
    ay_1 = calc_ay(x_1, y_1, z_1, vx_1, vy_1, vz_1, m, k)
    az_1 = calc_az(x_1, y_1, z_1, vx_1, vy_1, vz_1, m, k)

    speed = sqrt(vx_1**2 + vy_1**2 + vz_1**2)

    data['x_0'] = x_1
    data['y_0'] = y_1
    data['z_0'] = z_1
    data['vx_0'] = vx_1
    data['vy_0'] = vy_1
    data['vz_0'] = vz_1
    data['ax_0'] = ax_1
    data['ay_0'] = ay_1
    data['az_0'] = az_1
    data['speed'] = speed

    return data


def get_speeds_from_stage_0_to_stage_3(dt):
    speeds = [0]
    times = [0]
    distance = [0]
    acceleration = []

    x_0 = 6378000
    y_0 = 0
    z_0 = 0
    vx_0 = 0
    vy_0 = 0
    vz_0 = 0
    m = get_curr_mass(0)

    ax_0 = (-2 * m * vz_0 * 0.879 + 5849411 - 6.67 * 10**(-11) * (m * 5.97 * 10**24 * x_0) / sqrt(x_0**2 + y_0**2 + z_0**2)**3 -
            0.001 * vx_0) / m
    ay_0 = (2 * m * vz_0 * 7.29 * 10**-5 * sin(28.5) - 6.67 * 10**(-11) * (m * 5.97 * 10**24 * y_0) /
            sqrt(x_0**2 + y_0**2 + z_0**2)**3 - 0.001 * vy_0) / m
    az_0 = (2 * m * 7.29 * 10**-5 * vx_0 * 0.879 - 2 * m * vy_0 * 7.29**-5 * 0.477 - 6.67 * 10**(-11) * (m * 5.97 * 10**24 * z_0)
            / sqrt(x_0**2 + y_0**2 + z_0**2)**3 - 0.001 * vz_0) / m

    acceleration.append(sqrt(ax_0**2 + ay_0**2 + az_0**2))

    for t in count(dt, dt):
        if t > 115:
            break

        m = get_curr_mass(round(t, 1))
        data = runge_kutta(x_0, y_0, z_0, vx_0, vy_0, vz_0, ax_0, ay_0, az_0, m, 5849411, dt)
        x_0, y_0, z_0 = data['x_0'], data['y_0'], data['z_0']
        vx_0, vy_0, vz_0 = data['vx_0'], data['vy_0'], data['vz_0']
        ax_0, ay_0, az_0 = data['ax_0'], data['ay_0'], data['az_0']

        times.append(round(t, 1))
        speeds.append(data['speed'])
        acceleration.append(sqrt(ax_0 ** 2 + ay_0 ** 2 + az_0 ** 2))
        dist = sqrt((x_0 - 6378000)**2 + y_0**2 + z_0**2)
        distance.append(dist / 1000)

    for t in count(115 + dt, dt):
        if t > 262:
            break

        m = get_curr_mass(round(t, 1))
        data = runge_kutta(x_0, y_0, z_0, vx_0, vy_0, vz_0, ax_0, ay_0, az_0, m, 2339760, dt)
        x_0, y_0, z_0 = data['x_0'], data['y_0'], data['z_0']
        vx_0, vy_0, vz_0 = data['vx_0'], data['vy_0'], data['vz_0']
        ax_0, ay_0, az_0 = data['ax_0'], data['ay_0'], data['az_0']

        times.append(round(t, 1))
        speeds.append(data['speed'])
        acceleration.append(sqrt(ax_0 ** 2 + ay_0 ** 2 + az_0 ** 2))
        dist = sqrt((x_0 - 6378000) **2 + y_0 ** 2 + z_0 ** 2)
        distance.append(dist / 1000)

    for t in count(262 + dt, dt):
        if t > 467:
            break

        m = get_curr_mass(round(t, 1))
        data = runge_kutta(x_0, y_0, z_0, vx_0, vy_0, vz_0, ax_0, ay_0, az_0, m, 453714, dt)
        x_0, y_0, z_0 = data['x_0'], data['y_0'], data['z_0']
        vx_0, vy_0, vz_0 = data['vx_0'], data['vy_0'], data['vz_0']
        ax_0, ay_0, az_0 = data['ax_0'], data['ay_0'], data['az_0']

        times.append(round(t, 1))
        speeds.append(data['speed'])
        acceleration.append(sqrt(ax_0 ** 2 + ay_0 ** 2 + az_0 ** 2))
        dist = sqrt((x_0 - 6378000)**2 + y_0**2 + z_0**2)
        distance.append(dist / 1000)

    for t in count(467 + dt, dt):
        if t > 911:
            break

        m = get_curr_mass(round(t, 1))
        data = runge_kutta(x_0, y_0, z_0, vx_0, vy_0, vz_0, ax_0, ay_0, az_0, m, 131222, dt)
        x_0, y_0, z_0 = data['x_0'], data['y_0'], data['z_0']
        vx_0, vy_0, vz_0 = data['vx_0'], data['vy_0'], data['vz_0']
        ax_0, ay_0, az_0 = data['ax_0'], data['ay_0'], data['az_0']

        times.append(round(t, 1))
        speeds.append(data['speed'])
        acceleration.append(sqrt(ax_0 ** 2 + ay_0 ** 2 + az_0 ** 2))
        dist = sqrt((x_0 - 6378000)**2 + y_0 ** 2 + z_0**2)
        distance.append(dist / 1000)

    return times, speeds, acceleration, distance


def get_v4(v3):
    v4 = sqrt((86 * (v3)**2 / 2 - 6.67 * 10**-11 * 86 * 5.97 * 10**24 / 6378000 + 6.67 * 10**-11 * 86 * 5.97 * 10**24 /
          928 / 10**6) * 2 / 86)
    return v4

def get_v5(v4):
    v5 = sqrt((-6.67 * 10**-11 * 86 * 1.99 * 10**30 / 152.1 / 10**9 + 86 * v4**2 / 2 + 6.67 * 10**-11 * 86 * 1.99 * 10**30 /
          778.5 / 10**9) * 2 / 86)
    return v5


def get_v6(v5):
    v6 = 2 * (sqrt(13000**2 - cos(67.17)**2 + 13000**2 - v5**2) - 13000 * cos(67.17)) * sin((67.17 + 21.83) / 2) + v5
    return v6


def get_v7(v6):
    v7 = sqrt((86 * v6**2 / 2 - 6.67 * 10**-11 * 86 * 1.99 * 10**30 / 778.5 / 10**9 + 6.67 * 10**-11 * 86 * 1.99 * 10**30 /
          1430.39 / 10**9) * 2 / 86)
    return v7


def get_v8(v7):
    v8 = (2 * (sqrt(9690**2 * cos(64.5)**2 * cos(107.8)**2 + 9690**2 * cos(64.5)**2 - v7**2) - 9690 * cos(64.5) * cos(107.8)) *
          sin((107.8 + 18.7) / 2) + v7)
    return v8

