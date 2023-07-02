
import sys

from . import simu, utils

args = utils._parse_args(sys.argv[1:])

simu.simu_treated(args)
