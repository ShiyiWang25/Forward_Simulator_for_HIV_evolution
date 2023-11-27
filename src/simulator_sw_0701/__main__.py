import sys
from . import treated, untreated

def main():
    if len(sys.argv) < 2:
        print("Usage: python -m simulator_sw_0701 <treated|untreated> [args]")
        return

    script = sys.argv[1]
    if script == "treated":
        treated.simu(utils._parse_args(sys.argv[2:]))
    elif script == "untreated":
        untreated.simu(utils._parse_args(sys.argv[2:]))
    else:
        print("Invalid script name. Please specify 'treated' or 'untreated'")
        return

if __name__ == "__main__":
    main()
