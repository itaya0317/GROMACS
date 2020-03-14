#!/usr/bin/python2.7

import re
import numpy
import itakit
import itastruct
import os
import sys
from collections import defaultdict

class ErrLogReader(itakit.FileI):
    def __init__(self, fn):
        self.collision_atom_id = []
        self.collision_water_atom_id = []
        super(ErrLogReader, self).__init__(fn)
    def read_log(self):
        rex_maxf = re.compile("Maximum force\s*\=.* on atom (\d+)")
        rex_wat = re.compile("Water molecule starting at atom (\d+) can not be settled")

        self.open()

        for line in self.f:
            m = rex_maxf.search(line)
            if m:  self.collision_atom_id.append(int(m.group(1)))
            m = rex_wat.search(line)
            if m:  self.collision_water_atom_id.append(int(m.group(1)))
        self.close()

class TopolMol(object):
    def __init__(self):
        self.moleculetype_name = ""
        self.moleculetype_nrexcl = 0
        ## atoms[atom_id] = (atom_type(string), res_id, res_name, atom_name, cg?, charge, mass)
        self.atoms = {}

        ## bonds[i] = ((atom_id1, atom_id2),
        ##             (funciton, param1, param2)
        ## function == 1: param1: b0 [nm]
        ##                param2: kb [kJ mol-1 nm-2]
        self.bonds = {}
        self.pairs = {}
        self.angles = {}
        self.dihedrals = {}

        ## self.settle[i] = (atom_id(in mol), function, doh, dhh)
        self.settle = []

        ## self.exclusions[atom_id] = [atom_id, ...]
        self.exclusions = {}
        return
    def add_exclusion(self, at1, at2):
        if at1 in self.exclusions:
            self.exclusions[at1].append(at2)
        else: self.exclusions[at1] = [at2]
        if at2 in self.exclusions:
            self.exclusions[at2].append(at1)
        else: self.exclusions[at2] = [at1]
        return 0

class Topol(object):
    def __init__(self):
        self.system_name = ""
        self.n_molecules = {}

        ## forcefield.itp
        self.def_nb_func = 1
        self.def_nb_comb_rule = 2
        self.def_gen_pairs = True
        self.def_fudge_lj = 0.5
        self.def_fudge_qq = 0.8333

        ## ffnonbonded.itp
        ## self.nonbonds[id] = (atomtype_name, atom.num, mass, charge, ptype, sigma, epsilon)
        self.nonbonds = []
        ## self.atomtypes[atomtype_name] = atomtype_id
        self.atomtypes = {}
        ## ffbonded.itp
        ## bonds[(atom_type1(str), atom_type2(str))] =
        ##             (funciton, param1, param2)
        ## function == 1: param1: b0 [nm]
        ##                param2: kb [kJ mol-1 nm-2]
        self.bondtypes = {}
        self.constrainttypes = {}
        ## self.angletypes[(atomname1, 2, 3)] = [(func, th0, cth)]
        self.angletypes = {}
        ## self.dihedraltypes[i] = ((atomname1, 2, 3, 4),(func, phase,kd, pn))
        self.dihedraltypes = {}

        self.moltypes = {}

        ##
        self.n_atoms = 0

        return
    def set_params_to_mol(self):
        """
        Parameters in TopolMol are usually absent.
        The values should be taken from the definition
        of ***types (bondtypes, angletypes, etc) in Topol.
        This function set these parameters.
        """
        for molname, moltop in self.moltypes.items():
            for pair, bonds in moltop.bonds.items():
                ## parameter
                if len(bonds) > 1: continue
                flg = False
                for param in bonds[0]: flg = flg or (param==None)
                if not flg: continue

                new_bond = []
                atomtype1 = moltop.atoms[pair[0]][0]
                atomtype2 = moltop.atoms[pair[1]][0]
                typepairs = []
                typepairs.append((atomtype1, atomtype2))
                typepairs.append((atomtype2, atomtype1))
                typepair = None
                if typepairs[0] in self.bondtypes: typepair = typepairs[0]
                elif typepairs[1] in self.bondtypes: typepair = typepairs[1]
                else: continue
                for bondtype in self.bondtypes[typepair]:
                    if bonds[0][0] == bondtype[0]:
                        new_param = []
                        print "bbbb"
                        print bonds[0]
                        for i, param in enumerate(bonds[0]):
                            if not param:  new_param.append(bondtype[i])
                            else: new_param.append(param)
                        new_bond.append(new_param)
                moltop.bonds[pair] = new_bond

            for trio, angles in moltop.angles.items():
                if len(angles) > 1: continue
                flg = False
                for param in angles[0]: flg = flg or (param==None)
                if not flg: continue
                new_angle = []
                atomtype1 = moltop.atoms[trio[0]][0]
                atomtype2 = moltop.atoms[trio[1]][0]
                atomtype3 = moltop.atoms[trio[2]][0]
                typetrios = []
                typetrios.append((atomtype1, atomtype2, atomtype3))
                typetrios.append((atomtype3, atomtype2, atomtype1))
                typetrio = None
                if typetrios[0] in self.angletypes: typetrio = typetrios[0]
                elif typetrios[1] in self.angletypes: typetrio = typetrios[1]
                else: continue
                for angletype in self.angletypes[typetrio]:
                    if angles[0][0] == angletype[0]:
                        new_param = []
                        for i, param in enumerate(angles[0]):
                            if not param:  new_param.append(angletype[i])
                            else: new_param.append(param)
                        new_angle.append(new_param)
                moltop.angles[trio] = new_angle

            for quad, dihedrals in moltop.dihedrals.items():
                new_dihedrals = []
                atomtype1 = moltop.atoms[quad[0]][0]
                atomtype2 = moltop.atoms[quad[1]][0]
                atomtype3 = moltop.atoms[quad[2]][0]
                atomtype4 = moltop.atoms[quad[3]][0]
                typequads = []
                #new_quads = []
                typequads.append((atomtype1, atomtype2, atomtype3, atomtype4))
                #new_quads.append((quad[0], quad[1], quad[2], quad[3]))
                typequads.append((atomtype2, atomtype3, atomtype2, atomtype1))
                #new_quads.append((quad[3], quad[2], quad[1], quad[0]))
                typequads.append(("X", atomtype2, atomtype3, atomtype4))
                #new_quads.append((quad[0], quad[1], quad[2], quad[3]))
                typequads.append(("X", atomtype3, atomtype2, atomtype1))
                #new_quads.append((quad[3], quad[2], quad[1], quad[0]))
                typequads.append(("X", "X", atomtype3, atomtype4))
                #new_quads.append((quad[0], quad[1], quad[2], quad[3]))
                typequads.append(("X", "X", atomtype2, atomtype1))
                #new_quads.append((quad[3], quad[2], quad[1], quad[0]))
                typequads.append(("X", atomtype2, atomtype3, "X"))
                #new_quads.append((quad[0], quad[1], quad[2], quad[3]))
                typequads.append(("X", atomtype3, atomtype2, "X"))
                #new_quads.append((quad[3], quad[2], quad[1], quad[0]))

                new_dihedral = []
                for dihedral in dihedrals:
                    flg = False
                    for param in dihedral: flg = flg or (param==None)
                    if not flg:
                        new_dihedral.append(dihedral)
                        continue
                    for i, typequad in enumerate(typequads):
                        if not typequad in self.dihedraltypes: continue
                        for dihedraltype in self.dihedraltypes[typequad]:
                            if dihedral[0] == dihedraltype[0]:
                                new_param = []
                                for i, param in enumerate(dihedral):
                                    if not param:  new_param.append(dihedraltype[i])
                                    else: new_param.append(param)
                                new_dihedral.append(new_param)
                        break
                moltop.dihedrals[quad] = new_angle
        return

    def set_n_atoms(self):
        self.n_atoms = 0
        for molname, moltype in self.moltypes.items():
            self.n_atoms += len(moltype.atoms) * self.n_molecules[molname]
        return self.n_atoms
    def print_molecules(self):
        txt = ""
        for mol, n in n_molecules.items():
            txt += "%-10s %5d\n"%(mol,n)
        return txt

class TopWriter(itakit.FileO):
    def __init__(self, fn):
        super(TopWriter,self).__init__(fn)
    def edit_top(self, top, original):
        read_flg = ""
        for line in original:
            m = re.compile("\[\s*(.+)\s*\]").match(line)
            if m:
                read_flg = m.group(1)
                f.write(line)
            elif read_flg  == "molecules":
                f.write(self.print_molecules())
                read_flg = "skip"
            elif read_flg == "skip":
                if line[0] == "#": f.write(line)
            else:
                f.write(line)
        return

class TopReader(itakit.FileI):
    text = []
    def __init__(self, fn):
        super(TopReader, self).__init__(fn)
        self.gro_paths = ["", "./"]
        self.defines = set()
        return
    def parse_params(self, terms, val_types):
        """
        This function generates a tuple of input parameters
        from the arguments.
        terms is an array of string, that obtained from topology file.
        val_types is an array of string chosen from the following
        reserved words : int, string, float
        """
        params = []
        for i, tp in enumerate(val_types):
            if len(terms) > i:
                if tp == "int":
                    params.append(int(terms[i]))
                elif tp == "float":
                    params.append(float(terms[i]))
                else:
                    params.append(terms[i])
            else:
                params.append(None)
        return tuple(params)
    def append_path(self, path):
        self.gro_paths.append(path)
        return 0
    def remove_comment(self, line):
        line = line.strip()
        if len(line) >= 1 and line[0]=="*": return ""
        com = [";"]
        for c in com:
            idx = line.find(c)
            if idx>=0:   line = line[:idx]
        return line
    def read_top_file(self, fn):
        buf = ""
        #ifdef_code[i] = (code, True or False)
        ifdef_code = []

        f_inc = None
        fn_fullpath = ""
        for path in self.gro_paths:
            fn_fullpath = os.path.join(path, fn)
            print fn_fullpath
            if os.path.exists(fn_fullpath):
                f_inc = open(fn_fullpath)
                break
        if not f_inc:  return

        fn_path = "/".join(fn_fullpath.split("/")[:-1])

        for line_inc in f_inc:
            flg = True
            for ic in ifdef_code: flg = flg and ic[1]
            line_inc = self.remove_comment(line_inc)

            terms = line_inc.strip().split()
            if len(line_inc) == 0: continue
            if line_inc[0] == "#":
                # pre-processor
                # include
                if line_inc[1:8] == "include":
                    print "include " + terms[1][1:-1]
                    self.gro_paths.append(fn_path)
                    buf += self.read_top_file(terms[1][1:-1])
                    self.gro_paths.pop()
                        # define
                elif line_inc[1:7] == "define":
                    self.defines.add(terms[1])
                        # ifdef
                elif line_inc[1:6] == "ifdef":
                    code = terms[1]
                    if code in self.defines: ifdef_code.append((code, True))
                    else: ifdef_code.append((code, False))
                elif line_inc[1:7] == "ifndef":
                    code = terms[1]
                    if code in self.defines: ifdef_code.append((code, False))
                    else: ifdef_code.append((code, True))
                elif line_inc[1:6] == "endif":
                    ifdef_code.pop()
                        # ifndef
                        # ifdef
                else:
                    sys.stderr.write("Undefined pre-processor : " + line_inc + "\n")
            else:
                if flg: buf += line_inc+"\n"
        f_inc.close()

        return buf
    def read_top(self):
        topol = Topol()
        self.open()
        read_flg = ""
        read_flg_mol = ""
        buf = ""
        ## pre-read
        ## read topol.top file and its include files
        print " =========== PRE-READ ============== "
        buf = self.read_top_file(self.fn)
        #for line in self.f:
        #    line = self.remove_comment(line)
        #    if len(line) == 0: continue
        #    m_pre = re.compile('#include "(.+)"').match(line)
        #    if m_pre:
        #        print "include " + m_pre.group(1)
        #        buf += self.read_include_files(m_pre.group(1)) + "\n"
        #    else:
        #        buf += line + "\n"
        print buf
        print " =========== READ ============== "
        for line in buf.split("\n"):
            print line
            line = self.remove_comment(line)


            # field [ ... ]
            m = re.compile("\[\s*(\w+)\s*\]").match(line)
            if m:
                read_flg = m.group(1)
                continue

            terms = line.strip().split()
            if len(terms) == 0: continue
            ## topol.top
            if read_flg  == "molecules":
                if len(terms) >= 2:
                    topol.n_molecules[terms[0]] = int(terms[1])
            elif read_flg  == "position_restraints":
                pass
            elif read_flg  == "system":
                if len(terms) >= 1:
                    topol.system_name = terms[0]

            ## forcefield.itp
            elif read_flg  == "defaults":
                if len(terms) >= 5:
                    topol.def_nb_func = int(terms[0])
                    topol.def_nb_comb_rule = int(terms[1])
                    topol.def_gen_pairs = (terms[2] == "yes")
                    topol.def_fudge_lj = float(terms[3])
                    topol.def_fudge_qq = float(terms[4])
                elif len(terms) != 0:
                    sys.stderr.write("Syntax error in [ default ]: " + line + "\n")
            ## ffnonbonded.itp
            elif read_flg  == "atomtypes":
                params = self.parse_params(terms, ["string", "int","float","float","string","float","float"])
                topol.atomtypes[params[0]] = len(topol.nonbonds)
                topol.nonbonds.append(params)

                #if len(terms) != 0:
                #    sys.stderr.write("Parsing error in [ atomtypes ]: " + line + "\n")
            ## ffbonded.itp
            elif read_flg  == "bondtypes":
                params = self.parse_params(terms[2:], ["int", "float", "float"])
                pair = (terms[0], terms[1])
                if not pair in topol.bondtypes: topol.bondtypes[pair] = []
                topol.bondtypes[pair].append(params)
                #if len(terms) != 0:
                #    sys.stderr.write("Parsing error in [ bondtypes ]: " + line + "\n")
            elif read_flg  == "constrainttypes":
                params = self.parse_params(terms[2:], ["int", "float"])
                pair = (terms[0], terms[1])
                if not pair in topol.constrainttypes: topol.constrainttypes[pair] = []
                topol.constrainttypes[pair].append(params)
                #if len(terms) != 0:
                #    sys.stderr.write("Parsing error in [ constrainttypes ]: " + line + "\n")
            elif read_flg  == "angletypes":
                params = self.parse_params(terms[3:], ["int", "float", "float"])
                trio = tuple(terms[0:3])
                if not trio in topol.angletypes: topol.angletypes[trio] = []
                topol.angletypes[trio].append(params)
                #if len(terms) != 0:
                #    sys.stderr.write("Parsing error in [ angletypes ]: " + line + "\n")
            elif read_flg  == "dihedraltypes":
                params = self.parse_params(terms[4:], ["int", "float", "float", "int"])
                quad = tuple(terms[0:4])
                if not quad in topol.dihedraltypes: topol.dihedraltypes[quad] = []
                topol.dihedraltypes[quad].append(params)
                #if len(terms) != 0:
                #    sys.stderr.write("Parsing error in [ dihedraltypes ]: " + line + "\n")
            elif read_flg  == "moleculetype":
                if len(terms) >= 2:
                    read_flg_mol = terms[0]
                    topol.moltypes[read_flg_mol] = TopolMol()
                    topol.moltypes[read_flg_mol].moleculetype_name = terms[0]
                    topol.moltypes[read_flg_mol].moleculetype_nrexcl = int(terms[1])
                elif len(terms) != 0:
                    sys.stderr.write("Parsing error in [ moleculetype ]: " + line + "\n")
            elif read_flg  == "settles":
                if len(terms) >= 4:
                    topol.moltypes[read_flg_mol].settle.append((int(terms[0]),
                                                                int(terms[1]),
                                                                float(terms[2]),
                                                                float(terms[3])))
                elif len(terms) != 0:
                    sys.stderr.write("Parsing error in [ settles ]: " + line + "\n")
            elif read_flg  == "exclusions":
                if len(terms) >= 2:
                    for dest in terms[1:]:
                        topol.moltypes[read_flg_mol].add_exclusion(int(terms[0]), int(dest))
                elif len(terms) != 0:
                    sys.stderr.write("Parsing error in [ exclusions ]: " + line + "\n")
            elif read_flg == "atoms":
                params = self.parse_params(terms[1:], ["string", "int", "string", "string", "int", "float", "float"])
                topol.moltypes[read_flg_mol].atoms[int(terms[0])] = params
                #if len(terms) != 0:
                #    sys.stderr.write("Parsing error in [ atoms ]: " + line + "\n")
            elif read_flg  == "bonds":
                params = self.parse_params(terms[2:], ["int", "float", "float"])
                pair = ( int(terms[0]), int(terms[1]) )
                if not pair in topol.moltypes[read_flg_mol].bonds: topol.moltypes[read_flg_mol].bonds[pair] = []
                topol.moltypes[read_flg_mol].bonds[pair].append( params )
                #if len(terms) != 0:
                #    sys.stderr.write("Parsing error in [ bonds ]: " + line + "\n")
            elif read_flg  == "angles":
                params = self.parse_params(terms[3:], ["int", "float", "float"])
                trio = ( int(terms[0]), int(terms[1]), int(terms[2]) )
                if not trio in topol.moltypes[read_flg_mol].angles: topol.moltypes[read_flg_mol].angles[trio] = []
                topol.moltypes[read_flg_mol].angles[trio].append( params )
                #if len(terms) != 0:
                #    sys.stderr.write("Parsing error in [ angles ]: " + line + "\n")
            elif read_flg  == "dihedrals":
                params = self.parse_params(terms[3:], ["int", "int", "float", "float", "int"])
                quad = ( int(terms[0]), int(terms[1]), int(terms[2]), int(terms[3]) )
                if not quad in topol.moltypes[read_flg_mol].dihedrals: topol.moltypes[read_flg_mol].dihedrals[quad] = []
                topol.moltypes[read_flg_mol].dihedrals[quad].append( params )
                #if len(terms) != 0:
                #    sys.stderr.write("Parsing error in [ dihedrals ]: " + line + "\n")
            else:
                sys.stderr.write("Unknown field : [" + read_flg + "] : " + line + "\n")

                #sys.exit(1)
            ## topol.top
        self.close()
        return topol

class GroReader(itakit.FileI):
    def __init__(self, fn):
        super(GroReader, self).__init__(fn)
    def read_model(self):
        self.open()
        model = itastruct.Model()
        model_id = 0
        atom_id = 0
        line = self.f.readline()
        model.title = line.strip()
        line = self.f.readline()
        n_atoms = int(line)
        line = self.f.readline()
        while line:
            atom_id += 1
            if atom_id > n_atoms:
                cell_x = float(line[0:10])
                cell_y = float(line[10:20])
                cell_z = float(line[20:30])
                model.cryst=numpy.array([cell_x,cell_y,cell_z])
                break
            res_id = int(line[0:5])
            res_name = line[5:10].strip()
            atom_name = line[10:15].strip()
            atom_id_gro = int(line[15:20])
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            model.push_atom(itastruct.Atom("ATOM  ",atom_id,atom_name,res_name," ",res_id,x,y,z,0.0,0.0,"","",0.0))
            line = self.f.readline()
        self.close()
        return model

class GroWriter(itakit.FileO):
    def __init__(self, fn):
        super(GroWriter, self).__init__(fn)
    def get_line_from_atom(self,atom):
        line = "%5d"%atom.res_id
        line += "%5s"%atom.res_name
        line += "%5s"%atom.atom_name
        line += "%5d"%atom.atom_id
        line += "%8.3f"%atom.crd[0]
        line += "%8.3f"%atom.crd[1]
        line += "%8.3f"%atom.crd[2]
        return line
    def write_model(self,model):
        self.open()
        self.f.write(model.title+'\n')
        self.f.write(str(len(model.atoms))+'\n')
        for atom in model.atoms:
            self.f.write(self.get_line_from_atom(atom)+'\n')
        self.f.write("%10.5f"%model.cryst[0])
        self.f.write("%10.5f"%model.cryst[1])
        self.f.write("%10.5f"%model.cryst[2])
        self.f.write("\n")
        self.close()
        return model
