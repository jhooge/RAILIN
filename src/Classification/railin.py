#!/usr/bin/env python
'''
Created on 04.12.2012

@author: jhooge
'''
import subprocess as p
import argparse
from Pasta import Pasta
import sys, os
homedir = os.path.expanduser('~')
sys.path.append(homedir+'/Thesis/raip/src/Optimization')
sys.path.append(homedir+'/Thesis/raip/src/Visualization')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    ## positional arguments
    parser.add_argument("resfile", help="Path to list of pseudo-residues in ASCII format.",
                    type=str)
    parser.add_argument("seqfile", help="Path to sequence file in FASTA format.",
                    type=str)
    
    ## optional arguments
    parser.add_argument("--reffile",
                        type=str,   
                        default='shiftPresets/bmrb.shift',
                        help="Path to BMRB reference shifts in ASCII format.")
    parser.add_argument("-a", "--assignments", 
                        type=int,
                        default=100,
                        help="Number of assignments to be computed (default: 100).")
    parser.add_argument("-t", "--tolerance", 
                        type=float,
                        default=0.6,
                        help="Linking Tolerance in ppm (default: 0.6).")
    parser.add_argument("-s", "--strategy", type=str, 
                        choices={'Single','JointNew','Joint','ILP','CplexILP'},
                        default='CplexILP',
                        help="Assignment strategy (default: CplexILP).")
    parser.add_argument("-v", "--verbose", help="Turn on verbosity option.",
                        action="store_true")
    
    args = parser.parse_args()
    pasta = Pasta(args.resfile, args.reffile, args.seqfile)
    pasta.run(no_assignments=args.assignments,
              tolerance=args.tolerance,
              verbose=args.verbose, 
              map_strat=args.strategy)