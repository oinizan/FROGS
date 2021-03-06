#!/usr/bin/env python2.7
#
# Copyright (C) 2017 INRA
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
__author__ = 'Ta Thi Ngan & Maria Bernard INRA - SIGENAE'
__copyright__ = 'Copyright (C) 2017 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'frogs@inra.fr'
__status__ = 'prod'

import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']
# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from frogsUtils import *
##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################

class Rscript(Cmd):
    """
    @summary: Launch Rmarkdown script to present the alpha diversity of data with phyloseq.
    @see: http://rmarkdown.rstudio.com/
          https://joey711.github.io/phyloseq/
    @return: html file containing the plots
             alpha divesity table
    """
    def __init__(self, data, html, varExp, measures, alphaOut):
        """
        @param data: [str] path to phyloseq object in RData file, the result of FROGS Phyloseq Import Data.
        @param html: [str] Path to store resulting html file.
        @param varExp: [str] Experiment variable used to aggregate sample diversities.
        @param measures: [str] The indexes of alpha diversity, in list (Observed, Chao1, Shannon, Simpson, InvSimpson, ACE, Fisher).
        @param alphaOut: [str] Path to store resulting data file containing alpha diversity table.
        """ 
        rmd = os.path.join(CURRENT_DIR, "r_alpha_diversity.Rmd")
        Cmd.__init__( self,
                      'Rscript',
                      'Run 1 code Rmarkdown',
                       '-e "rmarkdown::render('+"'"+rmd+"',output_file='"+html+"', params=list(data='"+data+"', measures='"+measures+"', varExp='"+varExp+"',fileAlpha='"+alphaOut+"'), intermediates_dir='"+os.path.dirname(html)+"')"+'" 2> /dev/null',
                       "-e '(sessionInfo()[[1]][13])[[1]][1]; paste(\"Rmarkdown version: \",packageVersion(\"rmarkdown\")) ; library(phyloseq); paste(\"Phyloseq version: \",packageVersion(\"phyloseq\"))'")
    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stdout')

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
    
if __name__ == "__main__":
   
    # Manage parameters
    parser = argparse.ArgumentParser( description='To compute and present the data alpha diversity with plot_richness of Phyloseq.' )
    parser.add_argument('-v', '--varExp', type=str, required=True, default=None, help='The experiment variable used to aggregate sample diversities.' )
    parser.add_argument('-m', '--alpha-measures', type=str, nargs="*", default=['Observed','Chao1','Shannon','InvSimpson'], help='The indices of alpha diversity. Available indices : Observed, Chao1, Shannon, InvSimpson, Simpson, ACE, Fisher. [Default: %(default)s]')
  
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('-r','--rdata', required=True, default=None, help="The path of RData file containing a phyloseq object-the result of FROGS Phyloseq Import Data" )

    # output
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('-o','--html', default='alpha_diversity.html', help="The path to store resulting html file. [Default: %(default)s]" )
    group_output.add_argument('-a','--alpha-out', default='alpha_diversity.tsv', help="The path to store resulting data file containing alpha diversity table. [Default: %(default)s]" )    
    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help="This output file will contain several information on executed commands.")    
    args = parser.parse_args()
    prevent_shell_injections(args)   
    # Process 
    Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
    if len(args.alpha_measures) > 1 :
        for m in args.alpha_measures:
            if m not in ["Observed", "Chao1", "Shannon", "InvSimpson", "Simpson", "ACE", "Fisher"] :
                raise Exception("Measure, " + m + " is not a valid alpha diversity indice\n")
    else :
        for m in args.alpha_measures[0].split(","):
            if m not in ["Observed", "Chao1", "Shannon", "InvSimpson", "Simpson", "ACE", "Fisher"] :
                raise Exception("Measure, " + m + " is not a valid alpha diversity indice\n")

    data=os.path.abspath(args.rdata)
    html=os.path.abspath(args.html)
    alphaOut=os.path.abspath(args.alpha_out)
    measures=",".join(args.alpha_measures)
    Rscript(data, html, args.varExp, measures, alphaOut).submit( args.log_file )
