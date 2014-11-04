f = open('ftle.submit', 'w')

f.write('universe = vanilla\n')
f.write('executable = ftle\n')
f.write('should_transfer_files = YES\n')
f.write('when_to_transfer_output = ON_EXIT\n\n')

for t in range(6, 4946, 4):
    for direction in ["pos", "neg"]:
        t_str = str(t).zfill(4)
        f.write('output = ftle_' + t_str + '_' + direction + '.out\n')
        f.write('error = ftle_' + t_str + '_' + direction + '.error\n')
        f.write('log = ftle_' + t_str + '_' + direction + '.log\n')
        trans = ''
        for i in range(t-5, t+6):
            trans += 'data/u_' + str(i).zfill(4) + '.ascii,' + \
                     'data/v_' + str(i).zfill(4) + '.ascii,'
        trans = trans[:-1]
        f.write('transfer_input_files = ' + trans + '\n')
        if (direction == "pos"):
            f.write('arguments = ' + str(t) + ' 1 3 .5 5 2 \n')
        if (direction == "neg"):
            f.write('arguments = ' + str(t) + ' -1 3 .5 5 2 \n')

        f.write('queue\n\n')

f.close()
