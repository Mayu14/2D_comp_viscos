# coding: utf-8
import math
# from pathlib import Path
import os
import shutil


class FoamDict(object):
    def __init__(self, fields, values, depth=0):
        self.fields = fields
        self.values = values
        self.depth = depth

    def write(self):
        text = ''
        for i, field in enumerate(self.fields):
            text += self.depth * '\t' + field + "\t"
            if not isinstance(self.values[i], FoamDict):
                text += self.__get_internalField_value(self.values[i]) + " ;\n"
            else:
                indent = self.depth * '\t'
                text += '\n' + indent + '{\n' + self.values[i].write() + indent + '}\n'
        return text

    def __get_internalField_value(self, values):
        if (type(values) == type(1.0)) or (type(values) == type(1)):
            internalField = str(values)
        elif type(values) == type("1"):
            internalField = values
        else:
            internalField = "("
            for value in values:
                internalField += str(value) + " "
            internalField += ")"
        return internalField


class Config(object):
    def __init__(self, reynolds, mach, shape, angle, env="server"):
        self.mkdir = True
        self.output = True
        self.reynolds = reynolds
        self.mach = mach
        self.shape = shape
        self.angle_deg = angle
        self.angle = angle / 180.0 * math.pi
        self.env = env
        # self.rotation = gen_rotation_matrix(self.angle)
        # self.velocity = self.mach * np.array([np.cos(self.angle), np.sin(self.angle), 0.0])
        self.velocity = [self.mach * math.cos(self.angle), self.mach * math.sin(self.angle), 0.0]
        self.gamma = 1.4
        self.density = self.gamma * (mach ** 2)
        if math.isinf(reynolds):
            self.static_viscosity = (self.density * self.mach) / self.reynolds
            self.kinematic_viscosity = self.static_viscosity / self.density
        else:
            self.static_viscosity = 0.0
            self.kinematic_viscosity = 0.0
        self.load_default()

    def gen_casename(self):
        return self.shape + "_aoa" + str(self.angle_deg) + "_ma" + str(self.mach) + "_re" + str(self.reynolds)

    def __makedirs(self):
        """
        self.dir_0.mkdir(parents=True, exist_ok=True)
        self.dir_sys.mkdir(parents=True, exist_ok=True)
        self.dir_tri.mkdir(parents=True, exist_ok=True)
        """
        """
            os.makedirs(self.dir_0, exist_ok=True)
            os.makedirs(self.dir_sys, exist_ok=True)
            os.makedirs(self.dir_tri, exist_ok=True)
        """

        def check_and_make(dir):
            if not os.path.exists(dir):
                os.mkdir(dir)

        check_and_make(self.remote_dir)
        check_and_make(self.case_dir)
        check_and_make(self.dir_0)
        check_and_make(self.dir_sys)
        check_and_make(self.dir_const)
        check_and_make(self.dir_tri)
        return

    def load_default(self):
        def set_vector(x, y, z):
            # return np.array([x, y, z])
            return [x, y, z]

        self.casename = self.gen_casename()
        """
        self.source_dir = Path('/work/A/FMa/FMa037/source')
        self.remote_dir = Path('/work/A/FMa/FMa037/openFoam')
        self.dir_0 = (self.remote_dir / self.casename / '0')
        self.dir_sys = (self.remote_dir / self.casename / 'system')
        self.dir_const = (self.remote_dir / self.casename / 'constant')
        self.dir_tri = (self.dir_const / 'triSurface')
        """
        self.source_dir = '/work/A/FMa/FMa037/source'
        self.remote_dir = '/work/A/FMa/FMa037/openFoam_2'
        if self.env == "local":
            self.source_dir = 'D:\\Toyota\\source'
            self.remote_dir = 'D:\\Toyota\\work3'
        elif self.env == 'wsl':
            self.source_dir = '/mnt/d/Toyota/source'
            self.remote_dir = '/mnt/d/Toyota/work3'
        self.case_dir = os.path.join(self.remote_dir, self.casename)
        self.dir_0 = os.path.join(self.case_dir, '0')
        self.dir_sys = os.path.join(self.case_dir, 'system')
        self.dir_const = os.path.join(self.case_dir, 'constant')
        self.dir_tri = os.path.join(self.dir_const, 'triSurface')
        if self.mkdir:
            self.__makedirs()
        self.startTime = 0
        self.endTime = int(300 * 0.8 / self.mach)   # non-dimensionalize
        self.deltaT = 0.000001
        self.maxDeltaT = 0.1
        self.maxCourant = 0.9
        self.writeInterval = 5

        self.writeIntervalForces = 0.1
        self.CofR = set_vector(0, 0.5, 0)
        self.liftDir = set_vector(0, 1.0, 0)
        self.dragDir = set_vector(1.0, 0.0, 0.0)
        self.pitchAxis = set_vector(0, 0, 1)
        self.lRef = 1.0
        self.Aref = 0.05

    def __save_text(self, generator, obj_name, dir):
        text = generator(obj_name)
        if self.output:
            with open((os.path.join(dir, obj_name)), 'w') as f:
                f.write(text)
        return text

    def __get_internalField_value(self, values):
        if type(values) == type(1.0):
            internalField = str(values)
        elif type(values) == type("1"):
            internalField = values
        else:
            internalField = "("
            for value in values:
                internalField += str(value) + " "
            internalField += ")"
        return internalField

    def __write_from_fields_and_values(self, fields, values):
        dict = FoamDict(fields, values)
        return dict.write()

    def __gen_0(self, object="U"):
        def gen_0_body(inletVariable, inletValue, dimension, freestream_field, freestream_value, wall_field,
                       wall_value):
            boundary_field = ['freestream', 'wall', '#includeEtc']
            boundary_value = [FoamDict(freestream_field, freestream_value, depth=2),
                              FoamDict(wall_field, wall_value, depth=2), '"caseDicts/setConstraintTypes"']
            field_0 = [inletVariable, 'dimensions', 'internalField', 'boundaryField']
            value_0 = [inletValue, dimension, 'uniform $' + inletVariable,
                       FoamDict(boundary_field, boundary_value, depth=1)]
            return self.__write_from_fields_and_values(field_0, value_0)

        if object == 'U':
            class_type = 'volVectorField'
            inletVariable = 'Uinlet'
            inletValue = self.__get_internalField_value(self.velocity)
            dimensions = '[0 1 -1 0 0 0 0]'
            freestream_field = ['type', 'freestreamValue', 'value']
            freestream_value = ['freestreamVelocity', 'uniform $' + inletVariable, 'uniform $' + inletVariable]
            wall_field = ['type']
            wall_value = ['noSlip']
        elif object == 'nuTilda':
            class_type = 'volScalarField'
            inletVariable = 'nuTildaIn'
            inletValue = self.__get_internalField_value(3.0 * self.kinematic_viscosity)
            dimensions = '[0 2 -1 0 0 0 0]'
            freestream_field = ['type', 'fixedValue']
            freestream_value = ['nutkWallFunction', 'uniform $' + inletVariable]
            wall_field = ['type', 'value']
            wall_value = ['nutkWallFunction', 'uniform 0']
        else:
            print("please check gen_0's input")
            exit()
        zero = self.__gen_header(class_type, object)
        zero += gen_0_body(inletVariable, inletValue, dimensions, freestream_field, freestream_value, wall_field,
                           wall_value)
        zero += self.__gen_footer()
        return zero

    def __gen_constant(self, object='thermophysicalProperties'):
        class_type = 'dictionary'
        locaton = 'constant'
        const = self.__gen_header(class_type, object, locaton)
        if object == 'thermophysicalProperties':
            specie_value = FoamDict(['molWeight'], ['11640.3'], depth=2)
            thermo_value = FoamDict(['Cp', 'Hf'], ['2.5', '0.0'], depth=2)
            transport_value = FoamDict(['mu', 'Pr'], [str(self.static_viscosity), '1.0'], depth=2)
            mixture_field = ['specie', 'thermodynamics', 'transport']
            thermoType_field = ['type', 'mixture', 'transport', 'thermo', 'equationOfState', 'specie', 'energy']
            thermoType_value = ['hePsiThermo', 'pureMixture', 'const', 'hConst', 'perfectGas', 'specie',
                                'sensibleInternalEnergy']
            field_0 = ['thermoType', 'mixture']
            value_0 = [FoamDict(thermoType_field, thermoType_value, depth=1),
                       FoamDict(mixture_field, [specie_value, thermo_value, transport_value], depth=1)]
        else:
            print('please check __get_constant')
            exit()
        const += self.__write_from_fields_and_values(field_0, value_0)
        const += self.__gen_footer()
        return const

    def __gen_system(self, object='controlDict'):
        class_type = 'dictionary'
        system = self.__gen_header(class_type, object)
        if object == 'controlDict':
            force_field = ['type', 'libs', 'writeContarol', 'writeInterval',
                           'patches', 'log', 'rhoInf', 'CofR',
                           'liftDir', 'dragDir', 'pitchAxis', 'magUInf', 'lRef', 'Aref']
            force_value = ['forceCoeffs', '("libforces.so")', 'adjustableRunTime', str(self.writeIntervalForces),
                           '\n\t\t(\n\t\t\twall\n\t\t)', 'true', str(self.density), self.CofR,
                           self.liftDir, self.dragDir, self.pitchAxis, self.mach, self.lRef, self.Aref]

            func_field = ['#includeFunc', '#includeFunc', 'forces']
            func_value = ['MachNo', 'residuals', FoamDict(force_field, force_value, depth=2)]
            field = ['application', 'startFrom', 'startTime', 'stopAt', 'endTime',
                     'deltaT', 'maxDeltaT', 'writeControl', 'writeInterval',
                     'maxCo', 'adjustTimeStep', 'functions']
            value = ['rhoCentralFoam', 'latestTime', str(self.startTime), 'endTime', str(self.endTime),
                     str(self.deltaT), str(self.maxDeltaT), 'adjustableRunTime', str(self.writeInterval),
                     str(self.maxCourant), 'yes', FoamDict(func_field, func_value, depth=1)]
        elif object == 'extrudeMeshDict':
            field = ['constructFrom', 'sourceCase', 'sourcePatches',
                     'exposedPatchName', 'flipNormals', 'extrudeModel',
                     'nLayers', 'expansionRatio', 'linearNormalCoeffs',
                     'mergeFaces', 'mergeTol']
            value = ['patch', '"../' + self.casename + '"', '(symFront)',
                     'symBack', 'false', 'linearNormal',
                     1, 1.0, FoamDict(['thickness'], [self.Aref], depth=1),
                     'false', 0]
        else:
            print('please check __gen_system')
            exit()
        system += self.__write_from_fields_and_values(field, value)
        system += self.__gen_footer()
        return system

    def gen_U(self):
        return self.__save_text(generator=self.__gen_0, obj_name='U', dir=self.dir_0)

    def gen_nuTilda(self):
        return self.__save_text(generator=self.__gen_0, obj_name='nuTilda', dir=self.dir_0)

    def gen_thremoPhysicalPropaties(self):
        return self.__save_text(generator=self.__gen_constant, obj_name='thermophysicalProperties', dir=self.dir_const)

    def gen_controlDict(self):
        return self.__save_text(generator=self.__gen_system, obj_name='controlDict', dir=self.dir_sys)

    def gen_extrudeMeshDict(self):
        return self.__save_text(generator=self.__gen_system, obj_name='extrudeMeshDict', dir=self.dir_sys)

    def __gen_header(self, class_type="volVectorField", object_name="U", location=None):
        header = "/*--------------------------------*- C++ -*----------------------------------*\n"
        header += "  =========                 |\n"
        header += "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n"
        header += "   \\    /   O peration     | Website:  https://openfoam.org\n"
        header += "    \\  /    A nd           | Version:  6\n"
        header += "     \\/     M anipulation  |\n"
        header += "\*---------------------------------------------------------------------------*/\n"
        header += "FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       "
        header += class_type + ";\n"
        if type(location) != type(None):
            header += '    location   "' + location + '";\n'
        header += "    object      " + object_name + ";\n}\n"
        header += self.__gen_footer()
        return header

    def __gen_footer(self):
        return "//*************************************//\n"

    def gen_all(self):
        self.gen_U()
        self.gen_nuTilda()
        self.gen_thremoPhysicalPropaties()
        self.gen_controlDict()
        self.gen_extrudeMeshDict()
        self.copy_all()

    def copy_all(self):
        def get_src_and_dst(dir, name):
            src = os.path.join(self.source_dir, dir, name)
            dst = os.path.join(self.case_dir, dir, name)
            return src, dst

        zero = ['nut','p', 'T']
        const = ['turbulenceProperties']
        sys = ['blockMeshDict', 'createPatchDict', 'fvOptions',
               'fvSchemes', 'fvSolution', 'residuals', 'snappyHexMeshDict']
        nameList = [zero, const, sys]
        dirs = ['0', 'constant', 'system']
        for i, dir in enumerate(dirs):
            for name in nameList[i]:
                src, dst = get_src_and_dst(dir, name)
                shutil.copy(src, dst)

def main():
    reynoldsList = [1000, 10000]
    machList = [0.5, 0.8, 1.1, 1.4]
    shapeHeader = 'NACA'
    angleList = [-5.0, 0.0, 5.0, 10.0, 15.0, 20.0]

    #from naca_4digit_test import Naca_4_digit
    for i1 in range(1, 10, 3):
        for i2 in range(1, 10, 3):
            for i34 in range(12, 100, 20):
                int_4 = str(i1) + str(i2) + str(i34)
                shape = shapeHeader + int_4
                for angle in angleList:
                    """
                    naca = Naca_4_digit(int_4, attack_angle_deg = angle, resolution = 100, quasi_equidistant=False)
                    dir = 'D:\\Toyota\\work2\\obj\\'
                    fname = shapeHeader + int_4 + "_" + str(angle)
                    naca.generate_obj(dir + fname)
                    """
                    for reynolds in reynoldsList:
                        for mach in machList:
                            con = Config(reynolds, mach, shape, angle)
                            #print(con.casename)
                            #exit()
                            con.gen_all()

def load_patternList(naca4=True):
    reynoldsList = [1000, 10000]
    machList = [0.5, 0.8, 1.1, 1.4]
    angleList = [-5.0, 0.0, 5.0, 10.0, 15.0, 20.0]
    i3_list = [12, 22, 32, 42, 52, 62, 72, 82]
    # from naca_4digit_test import Naca_4_digit, Naca_5_digit
    if naca4:
        i2_len = 1
        i1_list = [0, 2, 5, 7, 9]
        i2_list = [0, 1, 3, 4, 6, 7]
        # naca_wing = Naca_4_digit
    else:
        i2_len = 2
        i1_list = [1, 2, 3, 4, 5]
        i2_list = [10, 21, 30, 41, 50]
        # naca_wing = Naca_5_digit
    return reynoldsList, machList, angleList, i1_list, i2_list, i3_list, i2_len

def main_remake(naca4=True):
    shapeHeader = 'NACA'
    #from naca_4digit_test import Naca_4_digit, Naca_5_digit
    if naca4:
        i2_len = 1
        #naca_wing = Naca_4_digit
    else:
        i2_len = 2
        #naca_wing = Naca_5_digit
    reynoldsList, machList, angleList, i1_list, i2_list, i3_list, i2_len = load_patternList(naca4)
    for i1 in i1_list:
        for i2 in i2_list:
            if ((i1 == 0) and (i2 == 0)) or (i1 != 0):
                for i34 in i3_list:
                    int_4 = str(i1) + str(i2).zfill(i2_len) + str(i34)
                    shape = shapeHeader + int_4
                    for angle in angleList:
                        #naca = naca_wing(int_4, attack_angle_deg = angle, resolution = 100, quasi_equidistant=False)
                        #dir = 'D:\\Toyota\\work2\\obj\\'
                        #fname = shapeHeader + int_4 + "_" + str(angle)
                        #naca.generate_obj(dir + fname)

                        for reynolds in reynoldsList:
                            for mach in machList:
                                con = Config(reynolds, mach, shape, angle)
                                #print(con.casename)
                                #exit()
                                con.gen_all()


def move_objs():
    obj_dir = "obj"
    dst_internal_dir = 'constant/triSurface/'
    parent = 'openFoam_2'
    dirList = os.listdir()
    target_name = 'wing_5degrees.obj'
    for dir in dirList:
        param = dir.split('_')
        if 'NACA' in param[0]:
            name = param[0] + "_" + param[1][3:] + ".obj"
            src = os.path.join(obj_dir, name)
            dst = os.path.join(parent, dir, dst_internal_dir, target_name)
            shutil.copy(src, dst)

if __name__ == '__main__':
    #main()
    #move_objs()
    main_remake(naca4=True)