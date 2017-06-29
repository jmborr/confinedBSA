import re

buffer = ''
iatom = 1
for line in open('SiO2carved_ovl1.5.pdb').readlines():
    if 'ATOM' in line:
        line = re.sub('ATOM\s+(\d+)\s+', 'ATOM{:7d} '.format(iatom), line)
        iatom += 1
        line = re.sub('SIL  1202', 'SIL     1', line)
    buffer += line
open('confiningSIL.pdb', 'w').write(buffer)

