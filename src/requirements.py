"""Command line requirements verification script."""

import sys
from subprocess import DEVNULL, check_call


def test_cmd(cmd):
    """Test if `cmd` is available to run as a subprocess."""
    try:
        check_call(cmd, stdout=DEVNULL, stderr=DEVNULL)
        return None
    except Exception as e:
        return e


def check_requirements():
    """Check mds_clustering external tools requirements."""
    print("Checking requirements...")
    errors = []
    for cmd in [
        ["git", "--version"],
        ["java", "-version"],
        ["javac", "-version"],
        ["R", "--version"],
        ["gcc", "--version"],
        ["Rscript", "--version"],
    ]:
        print(f"checking if {cmd[0]} is available... ", end="")
        ret = test_cmd(cmd)
        if ret is not None:
            errors.append(ret)
            print("no")
        else:
            print("yes")
    if errors:
        print("Runtime errors:", file=sys.stderr)
        for e in errors:
            print(e, file=sys.stderr)
        sys.exit(
            "Some requirements are not met. Please check your installation"
            "or follow the manual installation procedure."
        )


if __name__ == "__main__":
    check_requirements()
