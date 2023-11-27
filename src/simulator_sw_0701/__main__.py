
import sys
from . import simu, utils

args = utils._parse_args(sys.argv[1:])

if args.ART == "treated":
    simu.simu_treated(args)
elif args.ART == "untreated":
    simu.simu_untreated(args)
else:
    print("Invalid value for ART. Please specify 'treated' or 'untreated'.")
