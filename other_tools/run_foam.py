# -- coding: utf-8 --
from other_tools.setup_foam import Config, load_patternList
from other_tools.generate_runOF import gen_fsource, compile
from grid_generator.naca_4digit_test import Naca_4_digit, Naca_5_digit
import os

def main(env):
    shapeHeader = 'NACA'
    target_name = 'wing_5degrees.obj'
    naca4_list = [True, False]
    id = 0
    for naca4 in naca4_list:
        if naca4:
            naca_wing = Naca_4_digit
        else:
            naca_wing = Naca_5_digit

        reynoldsList, machList, angleList, i1_list, i2_list, i3_list, i2_len = load_patternList(naca4)

        for i1 in i1_list:
            for i2 in i2_list:
                if ((i1 == 0) and (i2 == 0)) or (i1 != 0):
                    for i34 in i3_list:
                        int_4 = str(i1) + str(i2).zfill(i2_len) + str(i34)
                        shape = shapeHeader + int_4
                        for angle in angleList:
                            naca = naca_wing(int_4, attack_angle_deg = angle, resolution = 100, quasi_equidistant=False)

                            for reynolds in reynoldsList:
                                for mach in machList:
                                    offset = id
                                    split = 1
                                    con = Config(reynolds, mach, shape, angle, env=env)
                                    fname = os.path.join(con.dir_tri, target_name)
                                    if not os.path.exists(fname):
                                        naca.generate_obj(fname)
                                        print('created')
                                    con.gen_all()

                                    gen_if_block = lambda offset, split: 'directory = "' + con.casename + '"\n'
                                    gen_fsource(offset=offset, split=split, use_offset=False, path='', mach=1.4, gen_if_block_all=gen_if_block)
                                    compile(split=split, offset=offset, run=True)
                                    id += 1


if __name__ == '__main__':
    main(env='local')