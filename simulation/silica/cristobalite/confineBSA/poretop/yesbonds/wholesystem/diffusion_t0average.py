import argparse
from amber.amber16.cpptraj.command import exec_cpptraj
from utilities.path import which
from tempfile import mkstemp
import numpy
from pdb import set_trace as tr

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Obtain MSD(t) averaged over different chunks of the trajectory.
    Example: diffusion_t0average.py pdb equil_rms2first.dcd 5000 '(:1-32)&(@H=)' diffusion_sty.dat""")
    parser.add_argument('topfile', help='input topology/PDB file')
    parser.add_argument('trajfile', help='DCD trajectory')
    parser.add_argument('nframe', type=int, help='number of frames in trajectory')
    parser.add_argument('dt', type=float, help='time in between frames in the input trajectory')
    parser.add_argument('span', type=int, help='number of frames for the MSD(t) chunk')
    parser.add_argument('nt0', type=int, help='number of chunks')
    parser.add_argument('mask', type=str, help='ptraj mask defining the atoms for which to calculate the MSD')
    parser.add_argument('outf', type=str, help='output diffusion file')
    parser.add_argument('--rms2t0', type=str, default='yes', help='for each chunk, do rms to its first snapshot. Default is yes')
    args=parser.parse_args()

    # Check appropriate cpptraj exists
    if not which('cpptraj'):
        raise Exception('cpptraj not found')
    msd = numpy.zeros(args.span)  #we substract 1 because cpptraj diffusion command does not calculate MSD(0), which is zero by definition
    template='''trajin _trajfile_  _t0_  _t0+t_
_rms_
diffusion time _dt_ _mask_ separateout xvg
'''
    handle,tmp_suffix = mkstemp(suffix='junk.xvg', dir='/tmp')
    jump = max(1, int((args.nframe - args.span - 1) / args.nt0))  # consecutive chunks start at positions that differ by 'jump'
    for it0 in range(args.nt0):
        t0 = 1 + it0*jump       # first frame of the chunk
        te = t0 + args.span -1  # last frame of the chunk
        script = template.replace('_trajfile_', args.trajfile)
        script = script.replace('_t0_', '{0}'.format(t0))
        script = script.replace('_t0+t_', '{0}'.format(te))
        script = script.replace('_mask_', args.mask)
        script = script.replace('_dt_', '{0}'.format(args.dt))
        if args.rms2t0 == 'y':
            script = script.replace('_rms_', 'rms first')
        else:
            script = script.replace('_rms_', '#')
        exec_cpptraj(args.topfile, script)
        file2mem = numpy.transpose(numpy.loadtxt('r_xvg')) #read MSD(t) for the chunk and store
        msd += file2mem[1]
    msd /= args.nt0 # average for all the chunks
    buf = '# time  MSD\n'
    for it in range(1,args.span):
        buf += '{0:5f} {1:8f}\n'.format( it*args.dt, msd[it-1] )
    open(args.outf, 'w').write(buf)

