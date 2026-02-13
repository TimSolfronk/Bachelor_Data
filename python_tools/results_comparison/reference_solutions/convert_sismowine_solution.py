scenario_name = "loh1"

tracers = {
    "loh1":
    [
        [0.000, 0., 0.693], [0.000, 0., 5.543], [0.000, 0., 10.392],
        [0.490, 0., 0.490], [3.919, 0., 3.919], [7.348, 0.,  7.348],
        [0.577, 0., 0.384], [4.612, 0., 3.075], [8.647, 0.,  5.764]
    ]
}

'''
This is a script that takes in SISMOWINE reference data located in the folder scenario_name.upper()
It then converts the data into a csv format that is very similar to ExaHyPE's CSV format
'''


def read_xyz_data(prefix):
    all_lines = []

    # read all files for v_x, v_y and v_z simultaneously to put fuse them into one array with all velocity values at a certain point and time
    with open(prefix+"/x.dat", "r", encoding="utf-8") as x, open(prefix+"/y.dat", "r", encoding="utf-8") as y, open(prefix+"/z.dat", "r", encoding="utf-8") as z:
        lines_per_tracer = 0
        tracer_lines_read = 0
        cur_tracer_id = 0
        cur_time = 0.0

        for x_line,y_line,z_line in zip(x.read().splitlines(),y.read().splitlines(),z.read().splitlines()):
            if x_line.startswith("#"):
                lines_per_tracer = int(x_line.split(" ")[0][1:])
                delta_time = float(x_line.split(" ")[1])
                continue
            if x_line == "":
                continue

            all_lines.append([str(cur_tracer_id),"0",str(cur_time),
                              str(tracers[scenario_name][cur_tracer_id][0]),
                              str(tracers[scenario_name][cur_tracer_id][1]),
                              str(tracers[scenario_name][cur_tracer_id][2]),
                              x_line,z_line,y_line])
            
            cur_time += delta_time
            tracer_lines_read += 1
            if tracer_lines_read == lines_per_tracer:
                cur_time = 0.0
                tracer_lines_read = 0
                cur_tracer_id += 1

    return all_lines


# read data of all tracers and put them into an array of ExaHyPE-like CSV lines 
lines_array = read_xyz_data(scenario_name.upper())

# write the lines into an output file
with open(scenario_name.upper() + "_ref_sismowine.csv", "w") as output_file:
    output_file.write("# data_on_fault = -\n")
    output_file.write("# data_off_fault = x_vel, y_vel, z_vel\n")
    output_file.write("number(0), number(1), t, x(0), x(1), x(2), data \n")
    for line in lines_array:
        output_file.write(", ".join(line) + "\n")