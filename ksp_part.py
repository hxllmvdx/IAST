import krpc
from time import sleep
from json import dump
import threading


def guide_vessel_to_maneuver(conn, vessel):
    vessel.control.activate_next_stage()
    vessel.control.throttle = 1.0
    vessel.control.sas = True
    sleep(1.5)
    pitch = conn.get_call(getattr, vessel.flight(), 'pitch')
    vessel.control.activate_next_stage()

    mean_altitude = conn.get_call(getattr, vessel.flight(), 'mean_altitude')
    expr = conn.krpc.Expression.greater_than(
        conn.krpc.Expression.call(mean_altitude),
        conn.krpc.Expression.constant_double(17700))
    event = conn.krpc.add_event(expr)
    with event.condition:
        event.wait()
    vessel.auto_pilot.target_pitch_and_heading(0, 90)
    sleep(1)
    vessel.auto_pilot.engage()
    expr = conn.krpc.Expression.less_than(
        conn.krpc.Expression.call(pitch),
        conn.krpc.Expression.constant_float(3))
    event = conn.krpc.add_event(expr)
    with event.condition:
        event.wait()
    vessel.auto_pilot.disengage()
    vessel.control.sas = True

    solid_fuel_amount = conn.get_call(vessel.resources.amount, 'SolidFuel')
    expr = conn.krpc.Expression.less_than(
        conn.krpc.Expression.call(solid_fuel_amount),
        conn.krpc.Expression.constant_float(0.1))
    event = conn.krpc.add_event(expr)
    with event.condition:
        event.wait()
    vessel.auto_pilot.target_pitch_and_heading(-25, 90)
    sleep(0.5)
    vessel.control.activate_next_stage()
    vessel.auto_pilot.engage()
    expr = conn.krpc.Expression.less_than(
        conn.krpc.Expression.call(pitch),
        conn.krpc.Expression.constant_float(-24))
    event = conn.krpc.add_event(expr)
    with event.condition:
        event.wait()
    vessel.auto_pilot.disengage()
    vessel.control.sas = True

    stage_6_resources = vessel.resources_in_decouple_stage(stage=6, cumulative=False)
    liquid_fuel_amount = conn.get_call(stage_6_resources.amount, 'LiquidFuel')
    expr = conn.krpc.Expression.less_than(
        conn.krpc.Expression.call(liquid_fuel_amount),
        conn.krpc.Expression.constant_float(0.1))
    event = conn.krpc.add_event(expr)
    with event.condition:
        event.wait()
    vessel.auto_pilot.target_pitch_and_heading(-45, 90)
    vessel.control.activate_next_stage()
    sleep(0.5)
    vessel.auto_pilot.engage()
    expr = conn.krpc.Expression.less_than(
        conn.krpc.Expression.call(pitch),
        conn.krpc.Expression.constant_float(-42))
    event = conn.krpc.add_event(expr)
    with event.condition:
        event.wait()
    vessel.auto_pilot.disengage()
    vessel.control.sas = True

    apoapsis_altitude = conn.get_call(getattr, vessel.orbit, 'apoapsis_altitude')
    expr = conn.krpc.Expression.greater_than(
        conn.krpc.Expression.call(apoapsis_altitude),
        conn.krpc.Expression.constant_double(250000))
    event = conn.krpc.add_event(expr)
    with event.condition:
        event.wait()
    vessel.control.throttle = 0
    sleep(1)
    vessel.control.activate_next_stage()
    vessel.control.activate_next_stage()
    vessel.control.sas_mode = vessel.control.sas_mode.prograde

    vessel.control.add_node(ut=94576455.18315232, prograde=1919.6012901481956, normal=364.49585151672363, radial=-102.89438092231738)
    time_to_node = conn.get_call(getattr, vessel.control.nodes[0], 'time_to')
    expr = conn.krpc.Expression.less_than(
        conn.krpc.Expression.call(time_to_node),
        conn.krpc.Expression.constant_double(182))
    event = conn.krpc.add_event(expr)
    with event.condition:
        event.wait()
    vessel.control.sas_mode = vessel.control.sas_mode.maneuver
    sleep(0.5)
    vessel.control.throttle = 1

    remaining_delta = conn.get_call(getattr, vessel.control.nodes[0], 'remaining_delta_v')
    expr = conn.krpc.Expression.less_than(
        conn.krpc.Expression.call(remaining_delta),
        conn.krpc.Expression.constant_double(10))
    event = conn.krpc.add_event(expr)
    with event.condition:
        event.wait()
    vessel.control.throttle = 0
    sleep(0.5)
    vessel.control.activate_next_stage()
    vessel.control.nodes[0].remove()


def get_statistics(conn, vessel):
    data = []
    sleep(1)
    start_time = conn.space_center.ut
    reached_jool = False
    stage = 0
    last_mass = vessel.dry_mass
    srf_frame = vessel.orbit.body.reference_frame
    sun_frame = conn.space_center.bodies['Sun'].non_rotating_reference_frame

    while conn:
        time_from_launch = conn.space_center.ut - start_time
        if last_mass - vessel.dry_mass > 10:
            last_mass = vessel.dry_mass
            stage += 1
        if vessel.orbit.body.name == 'Jool':
            reached_jool = True
        if reached_jool and vessel.orbit.body.name == 'Sun':
            conn = False

        if stage in [0, 1, 2, 3]:
            print(time_from_launch, vessel.flight(srf_frame).speed, vessel.flight().pitch, stage,
                  vessel.orbit.body.name)
            data.append([time_from_launch, vessel.flight(srf_frame).speed, vessel.flight().pitch, stage,
                         vessel.orbit.body.name])
            sleep(0.5)
        else:
            print(time_from_launch, vessel.flight(sun_frame).speed, stage, vessel.orbit.body.name)
            data.append([time_from_launch, vessel.flight(sun_frame).speed, stage, vessel.orbit.body.name])
            sleep(0.5)

    with open('flight_data.json', 'w') as f:
        dump(data, f)


def main():
    conn = krpc.connect(name='Voyger-1')
    vessel = conn.space_center.active_vessel
    t1 = threading.Thread(target=guide_vessel_to_maneuver, args=(conn, vessel), daemon=True)
    t2 = threading.Thread(target=get_statistics, args=(conn, vessel), daemon=True)
    t1.start()
    t2.start()
    t1.join()
    t2.join()
