import matplotlib.pyplot as plt
from json import load
import math_part as mp


def get_data_from_ksp():
    with open('flight_data.json') as f:
        data = load(f)
    return data


def draw_speed_plots(x1, y1, x2, y2):
    plt.plot(x1, y1, x2, y2)
    plt.xlabel('Time in seconds')
    plt.ylabel('Speed in m/s')
    plt.title('Speed graphs')
    plt.show()


def draw_altitude_plots(x1, y1, x2, y2):
    plt.plot(x1, y1, x2, y2)
    plt.xlabel('Time in seconds')
    plt.ylabel('Altitude in meters')
    plt.title('Altitude graphs')
    plt.show()


def draw_acceleration_plot(x, y):
    plt.plot(x, y)
    plt.xlabel('Time in seconds')
    plt.ylabel('Acceleration in m/s^2')
    plt.title('Acceleration graphs')
    plt.show()


def main():
    data = get_data_from_ksp()
    times_ksp = data['stage_1to3']['time']
    speeds_ksp = data['stage_1to3']['speed']
    altitudes_ksp = data['stage_1to3']['alt']
    res = mp.get_speeds_from_stage_0_to_stage_3(0.1)
    times_math = res[0]
    speeds_math = res[1]
    accelerations_math = res[2]
    altitudes_math = res[3]
    draw_altitude_plots(times_ksp, altitudes_ksp, times_math, altitudes_math)
    draw_speed_plots(times_ksp, speeds_ksp, times_math, speeds_math)
    draw_acceleration_plot(times_math, accelerations_math)

