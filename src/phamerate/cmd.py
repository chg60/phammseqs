from subprocess import Popen, PIPE, DEVNULL
import shlex


def run_command(command, verbose=False):
    """
    Run the indicated command as a subprocess.

    :param command: the command to run
    :type command: str
    :param verbose: want stdout/stderr?
    :type verbose: bool
    :return: out, err
    """
    command = shlex.split(command)
    if verbose:
        with Popen(command, stdout=PIPE, stderr=PIPE) as sp:
            out = sp.stdout.read().decode("utf-8")
            err = sp.stderr.read().decode("utf-8")
            return out, err
    else:
        with Popen(command, stdout=DEVNULL, stderr=DEVNULL) as sp:
            sp.wait()
