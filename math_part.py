import json
from math import sqrt, sin, cos
from itertools import count


def get_curr_mass(t):
    if 0 <= t <= 80:
        return 182400 - 1632.5 * t
    elif 80 < t <= 136:
        return 30600 - 132.142 * (t - 80)
    elif 136 < t <= 185:
        return 17300 - 36 * (t - 136)
    elif 185 < t <= 451:
        return 8586 - 14.184 * (t - 185)


def calc_ax(x, y, z, vx, vy, vz, m, k):
    return (-2 * m * vz * cos(1) + k * vx / sqrt(vx**2 + vy**2 + vz**2) - 6.67 * 10**(-11) * (m * 5.291 * 10**22 * x) /
          sqrt(x**2 + y**2 + z**2)**3 - 0.001 * vx) / m


def calc_ay(x, y, z, vx, vy, vz, m, k):
    return (2 * m * vz * 2.91 * 10**-4 * sin(1) + k * vy / sqrt(vx**2 + vy**2 + vz**2) - 6.67 * 10**(-11) *
          (m * 5.291 * 10 ** 22 * y) / sqrt(x ** 2 + y ** 2 + z ** 2)**3 - 0.001 * vy) / m


def calc_az(x, y, z, vx, vy, vz, m, k):
    return (2 * m * 2.91 * 10**-4 * vx * cos(1) - 2 * m * vy * 2.91 * 10**-4 * sin(1) + k * vz / sqrt(vx**2 + vy**2 + vz**2) -
            6.67 * 10**(-11) * (m * 5.291 * 10 ** 22 * z) / sqrt(x**2 + y**2 + z**2)**3 - 0.001 * vz) / m

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

    x_0 = 600000
    y_0 = 0
    z_0 = 0
    vx_0 = 0
    vy_0 = 0
    vz_0 = 0
    m = get_curr_mass(0)

    ax_0 = (-2 * m * vz_0 * cos(1) + 3800000 - 6.67 * 10**(-11) * (m * 5.291 * 10 ** 22 * x_0) / sqrt(x_0**2 + y_0**2 + z_0**2)**3 -
            0.001 * vx_0) / m
    ay_0 = (2 * m * vz_0 * 2.91 * 10**-4 * sin(1) - 6.67 * 10**(-11) * (m * 5.291 * 10 ** 22 * y_0) /
            sqrt(x_0**2 + y_0**2 + z_0**2)**3 - 0.001 * vy_0) / m
    az_0 = (2 * m * 2.91 * 10**-4 * vx_0 * cos(1) - 2 * m * vy_0 * 2.91 * 10**-4 * sin(1) - 6.67 * 10**(-11) * (m * 5.291 * 10 ** 22 * z_0)
            / sqrt(x_0**2 + y_0**2 + z_0**2)**3 - 0.001 * vz_0) / m

    acceleration.append(sqrt(ax_0**2 + ay_0**2 + az_0**2))

    for t in count(dt, dt):
        if t > 80:
            break

        m = get_curr_mass(round(t, 1))
        data = runge_kutta(x_0, y_0, z_0, vx_0, vy_0, vz_0, ax_0, ay_0, az_0, m, 3800000, dt)
        x_0, y_0, z_0 = data['x_0'], data['y_0'], data['z_0']
        vx_0, vy_0, vz_0 = data['vx_0'], data['vy_0'], data['vz_0']
        ax_0, ay_0, az_0 = data['ax_0'], data['ay_0'], data['az_0']

        times.append(round(t, 1))
        speeds.append(data['speed'])
        acceleration.append(sqrt(ax_0 ** 2 + ay_0 ** 2 + az_0 ** 2))
        dist = sqrt((x_0 - 600000)**2 + y_0**2 + z_0**2)
        distance.append(dist)

    for t in count(80 + dt, dt):
        if t > 136:
            break

        m = get_curr_mass(round(t, 1))
        data = runge_kutta(x_0, y_0, z_0, vx_0, vy_0, vz_0, ax_0, ay_0, az_0, m, 400000, dt)
        x_0, y_0, z_0 = data['x_0'], data['y_0'], data['z_0']
        vx_0, vy_0, vz_0 = data['vx_0'], data['vy_0'], data['vz_0']
        ax_0, ay_0, az_0 = data['ax_0'], data['ay_0'], data['az_0']

        times.append(round(t, 1))
        speeds.append(data['speed'])
        acceleration.append(sqrt(ax_0 ** 2 + ay_0 ** 2 + az_0 ** 2))
        dist = sqrt((x_0 - 600000) **2 + y_0 ** 2 + z_0 ** 2)
        distance.append(dist)

    for t in count(136 + dt, dt):
        if t > 185:
            break

        m = get_curr_mass(round(t, 1))
        data = runge_kutta(x_0, y_0, z_0, vx_0, vy_0, vz_0, ax_0, ay_0, az_0, m, 125000, dt)
        x_0, y_0, z_0 = data['x_0'], data['y_0'], data['z_0']
        vx_0, vy_0, vz_0 = data['vx_0'], data['vy_0'], data['vz_0']
        ax_0, ay_0, az_0 = data['ax_0'], data['ay_0'], data['az_0']

        times.append(round(t, 1))
        speeds.append(data['speed'])
        acceleration.append(sqrt(ax_0 ** 2 + ay_0 ** 2 + az_0 ** 2))
        dist = sqrt((x_0 - 600000)**2 + y_0**2 + z_0**2)
        distance.append(dist)

    for t in count(185 + dt, dt):
        if t > 451:
            break

        m = get_curr_mass(round(t, 1))
        data = runge_kutta(x_0, y_0, z_0, vx_0, vy_0, vz_0, ax_0, ay_0, az_0, m, 48750, dt)
        x_0, y_0, z_0 = data['x_0'], data['y_0'], data['z_0']
        vx_0, vy_0, vz_0 = data['vx_0'], data['vy_0'], data['vz_0']
        ax_0, ay_0, az_0 = data['ax_0'], data['ay_0'], data['az_0']

        times.append(round(t, 1))
        speeds.append(data['speed'])
        acceleration.append(sqrt(ax_0 ** 2 + ay_0 ** 2 + az_0 ** 2))
        dist = sqrt((x_0 - 600000)**2 + y_0 ** 2 + z_0**2)
        distance.append(dist)

    return times, speeds, acceleration, distance


def get_v4(v3):
    v4 = sqrt((960 * v3**2 / 2 + 6.67 * 10**-11 * 960 * 5.291 * 10**22 / 600000 - 6.67 * 10**-11 * 960 * 5.291 * 10**22 /
          84159286) * 2 / 86)
    return v4

def get_v5(v4):
    v5 = sqrt((6.67 * 10**-11 * 960 * 1.756 * 10**28 / 13599840256 + 960 * v4**2 / 2 - 6.67 * 10**-11 * 960 * 1.756 * 10**28 /
          68773560320) * 2 / 86)
    return v5


def get_v6(v5):
    v_in = sqrt(v5**2 + 4134**2 - 2 * v5 * 4134 * cos(60.1))
    v6 = sqrt(4134**2 + v_in**2 - 2 * v_in * 4134 * cos(110))
    return v6


def compare_speeds():
    with open('flight_data.json') as f:
        data = json.load(f)
    res = get_speeds_from_stage_0_to_stage_3(0.1)
    speed = res[1][-1]
    v4 = get_v4(speed + 9285)
    v5 = get_v5(v4 + 230000)
    v6 = get_v6(v5)
    print(f'{"По модели"}     {"По KSP"}')
    print(f'{v4:.3f}     {data["4"]:.3f}')
    print(f'{v5:.3f}    {data["5"]:.3f}')
    print(f'{v6:.3f}    {data["6"]:.3f}')

